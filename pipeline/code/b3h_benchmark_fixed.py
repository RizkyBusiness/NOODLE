"""The mouse model error floor, with the crystal preparation fixed.

B3e reused benchmark_threshold.py's crystal machinery unchanged and produced a p95 of
10.2 A - not an error floor but a broken measurement. The diagnostic (B3e_diagnostic.csv)
showed framework fit and site RMSD correlating at 0.992, i.e. the large values are
superposition failures. Two defects in that machinery:

1. REGISTER. imgt_renumber() walks the OBSERVED residues of a chain in lockstep with an
   ANARCI numbering computed on the ENTITY (SEQRES) sequence, starting at query_start.
   Crystals routinely have unmodelled residues, so any gap in or before the V domain
   throws the two out of register; the residue-identity guard then skips nearly every
   residue, leaving chains with 2-10 mapped residues (7Z50, 5M02, 6BGA, 9NV6, 3ARB...).
   FIX: number each candidate chain on its OWN observed sequence, so the numbering and
   the residues being mapped are the same string and cannot slip.

2. PAIRING. Chains are taken as the first auth_asym_id of the alpha entity and the first
   of the beta entity, independently. In an entry with several TCR copies in the
   asymmetric unit that can pair an alpha from one copy with a beta from another, and the
   lattice separation is then measured as model error (3E3Q 62.5 A with a 59.8 A framework
   fit; 2E7L 27.0 A with 21.9 A).
   FIX: enumerate every TCR chain in the coordinates, then choose the alpha-beta pair
   whose V-domain centroids are closest. A genuinely paired alpha-beta TCR has its two V
   domains in contact, so proximity identifies same-copy partners without reference to
   the model being tested.

The model is still folded from the COMPLETE entity-derived V-domain sequence, as the
packaged benchmark does and as the repertoire models are, so disorder in the crystal
reduces the positions compared rather than the sequence folded. Site RMSD is measured in
the framework frame shared by that crystal-model pair (min_frac=1.0), unchanged.

Acceptance gates, fixed in advance and reported with their counts:
  - both chains mapped: >= 90 IMGT-numbered residues each (a V domain is ~110)
  - both chains' site loops present: >= 24 site positions each
  - V-domain centroid separation < 40 A (same-copy sanity check)
  - framework fit < 2.0 A (a fit worse than this means the superposition, not the
    model, is what is being measured)

Writes: tables/B3h_benchmark_structures.csv  per entry, with all gate quantities
        tables/threshold.json                mouse-derived, what clustering reads
        tables/B3h_threshold_provenance.json selection, gates and exclusion counts
"""
import hashlib
import json
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

import benchmark_threshold as BT  # noqa: E402
import imgt  # noqa: E402
import site_geometry as G  # noqa: E402

MIN_MAPPED = 90
MIN_SITE_PER_CHAIN = 24
MAX_CENTROID_SEP = 40.0
MAX_FR_FIT = 2.0


def renumber_observed(residues, numbering, query_start, label):
    """Rewrite ATOM lines using a numbering computed on this chain's OWN sequence."""
    lines, i = [], query_start
    for (num, ins), aa in numbering:
        if aa == "-":
            continue
        if i >= len(residues):
            break
        _key, obs_aa, atoms = residues[i]
        i += 1
        if obs_aa != aa or num > 128:
            continue
        for L in atoms:
            lines.append(L[:21] + label + "%4d%s" % (num, ins if ins != " " else " ")
                         + L[27:])
    return lines


def ca_centroid(lines):
    P = [(float(L[30:38]), float(L[38:46]), float(L[46:54]))
         for L in lines if L[12:16].strip() == "CA"]
    return np.array(P).mean(0) if P else None


def prepare(eid, entities, res, cfg, chains, types):
    """-> {chain: (lines, vseq, n_mapped)}, centroid separation. Raises on failure."""
    # entity-level numbering supplies the COMPLETE V-domain sequence to fold
    ent_num = BT.number_chains([("%s|e%d" % (eid, k), s)
                                for k, (_ch, s) in enumerate(entities)], cfg)
    # candidate chains: every auth chain of every entity that numbered as a TCR chain
    cands = {c: [] for c in chains}
    for k, (auth_ids, _s) in enumerate(entities):
        key = "%s|e%d" % (eid, k)
        if key not in ent_num:
            continue
        ct, ent_numbering, _qs = ent_num[key]
        lab = next((c for c in chains if ct in types[c]), None)
        if lab is None:
            continue
        vseq = "".join(a for _, a in ent_numbering if a != "-")
        for auth in auth_ids:
            if auth not in res:
                continue
            obs = "".join(aa for _k, aa, _at in res[auth])
            if len(obs) < 60:
                continue
            onum = BT.number_chains([("%s|%s" % (eid, auth), obs)], cfg)
            got = onum.get("%s|%s" % (eid, auth))
            if not got:
                continue
            _ct2, numbering, qs = got
            lines = renumber_observed(res[auth], numbering, qs, lab)
            cen = ca_centroid(lines)
            if cen is None:
                continue
            n_mapped = len({L[22:27] for L in lines})
            cands[lab].append(dict(auth=auth, lines=lines, vseq=vseq,
                                   n_mapped=n_mapped, centroid=cen))
    for c in chains:
        if not cands[c]:
            raise KeyError("no usable %s chain in coordinates" % c)
    # pair the alpha and beta candidates whose V domains are closest = same copy
    best, sep = None, None
    for a in cands[chains[0]]:
        for b in cands[chains[1]]:
            d = float(np.linalg.norm(a["centroid"] - b["centroid"]))
            if sep is None or d < sep:
                best, sep = (a, b), d
    return {chains[0]: best[0], chains[1]: best[1]}, sep


def main(config_path="config_main.json"):
    cfg = imgt.load(config_path)
    out, b = cfg["out_dir"], cfg["benchmark"]
    chains = imgt.chain_labels(cfg)
    types = {c: imgt.SHARED_CHAIN_TYPES.get(c, (c,)) for c in chains}
    work = os.path.join(out, "benchmark")
    os.makedirs(work, exist_ok=True)

    T = pd.read_csv(os.path.join(out, "B3b_benchmark_structure_survey.csv"))
    sel = T[(T.tcr_species == "mouse") & (~T.redundant_vseq)].copy()
    print("mouse non-redundant entries: %d (res %.2f-%.2f, median %.2f)"
          % (len(sel), sel.resolution.min(), sel.resolution.max(),
             sel.resolution.median()), flush=True)
    resol = dict(zip(sel.pdb_id, sel.resolution))

    import ib_patches
    ib_patches.apply()
    from ImmuneBuilder import TCRBuilder2
    pred = TCRBuilder2(weights_dir=cfg["folding"]["weights_dir"])

    ents = BT.entity_sequences(list(sel.pdb_id))
    rows = []
    for eid in sel.pdb_id:
        rec = dict(entry=eid, resolution=resol.get(eid))
        try:
            raw = os.path.join(work, "%s.pdb" % eid)
            if not os.path.exists(raw):
                urllib.request.urlretrieve(BT.PDB_FILE % eid, raw)
            res = BT.read_pdb_chain_residues(raw)
            pick, sep = prepare(eid, ents.get(eid) or [], res, cfg, chains, types)
            rec.update(centroid_sep=round(sep, 2),
                       **{"auth_%s" % c: pick[c]["auth"] for c in chains},
                       **{"mapped_%s" % c: pick[c]["n_mapped"] for c in chains})
            vseqs = {c: pick[c]["vseq"] for c in chains}
            cryst = os.path.join(work, "fixed_%s_crystal.pdb" % eid)
            with open(cryst, "w") as fh:
                fh.writelines(sum((pick[c]["lines"] for c in chains), []) + ["END\n"])
            tag = hashlib.sha1("|".join(vseqs[c] for c in chains).encode()).hexdigest()[:10]
            model = os.path.join(work, "seq_%s_model.pdb" % tag)
            if not os.path.exists(model):
                pred.predict(vseqs).save(model, n_threads=cfg["folding"]["threads"])
            ca_c, ca_m = G.read_ca(cryst, chains), G.read_ca(model, chains)
            fr = G.framework_positions([ca_c, ca_m], chains, min_frac=1.0)
            Fc, Fm = G.framework_matrix(ca_c, fr), G.framework_matrix(ca_m, fr)
            R, cB, cA = G.kabsch(Fm, Fc)
            fit = float(np.sqrt((((Fm - cB) @ R + cA - Fc) ** 2).sum(1).mean()))
            rmsd, n_site = G.site_rmsd_between(ca_c, ca_m, fr, chains)
            per = {c: len([k for k in ca_c[c]
                           if any(k[0] in imgt.SITE_RANGES[nm]
                                  for _c, nm in imgt.site_loops(cfg) if _c == c)])
                   for c in chains}
            rec.update(site_rmsd=rmsd, n_site_positions=n_site, n_framework=len(fr),
                       fr_fit=round(fit, 3),
                       **{"site_%s" % c: per[c] for c in chains},
                       **{"len_%s" % c: len(vseqs[c]) for c in chains})
            print("  %s  site RMSD %s A | fit %.2f | site %d (%d/%d) | sep %.1f"
                  % (eid, ("%.2f" % rmsd) if rmsd else "nan", fit, n_site or 0,
                     per[chains[0]], per[chains[1]], sep), flush=True)
        except Exception as e:
            rec.update(error="%s: %s" % (type(e).__name__, str(e)[:140]))
            print("  %s FAILED %s: %s" % (eid, type(e).__name__, str(e)[:90]), flush=True)
        rows.append(rec)

    B = pd.DataFrame(rows)
    for c in chains:
        for p in ("mapped_", "site_"):
            if p + c not in B:
                B[p + c] = np.nan
    gate = (B.site_rmsd.notna()
            & (B["mapped_%s" % chains[0]] >= MIN_MAPPED)
            & (B["mapped_%s" % chains[1]] >= MIN_MAPPED)
            & (B["site_%s" % chains[0]] >= MIN_SITE_PER_CHAIN)
            & (B["site_%s" % chains[1]] >= MIN_SITE_PER_CHAIN)
            & (B.centroid_sep < MAX_CENTROID_SEP)
            & (B.fr_fit < MAX_FR_FIT))
    B["accepted"] = gate
    B.to_csv(os.path.join(out, "B3h_benchmark_structures.csv"), index=False)
    ok = B[gate]

    excl = dict(
        no_measurement=int(B.site_rmsd.isna().sum()),
        one_chain_poorly_mapped=int((B.site_rmsd.notna()
                                     & ((B["mapped_%s" % chains[0]] < MIN_MAPPED)
                                        | (B["mapped_%s" % chains[1]] < MIN_MAPPED))).sum()),
        site_loops_incomplete=int((B.site_rmsd.notna()
                                   & ((B["site_%s" % chains[0]] < MIN_SITE_PER_CHAIN)
                                      | (B["site_%s" % chains[1]] < MIN_SITE_PER_CHAIN))).sum()),
        chains_not_same_copy=int((B.site_rmsd.notna()
                                  & (B.centroid_sep >= MAX_CENTROID_SEP)).sum()),
        framework_fit_too_poor=int((B.site_rmsd.notna() & (B.fr_fit >= MAX_FR_FIT)).sum()))

    if len(ok) < 5:
        raise SystemExit("only %d structures passed the gates" % len(ok))
    p = b["percentile"]
    val = float(np.percentile(ok.site_rmsd, p))
    thr = max(round(val * 20) / 20.0, b["min_threshold"])
    rec = dict(threshold=thr,
               estimator="p%d of model-vs-crystal combining-site RMSD" % p,
               percentile=p, raw_percentile_value=round(val, 3),
               median=round(float(ok.site_rmsd.median()), 3),
               twice_median=round(2 * float(ok.site_rmsd.median()), 3),
               n_structures=int(len(ok)), n_attempted=int(len(B)),
               min_threshold_floor=b["min_threshold"],
               max_resolution=b["max_resolution"], species=cfg["species"],
               locus_pair=cfg["locus_pair"], builder=cfg["folding"]["builder"],
               entries=sorted(ok.entry))
    json.dump(rec, open(os.path.join(out, "threshold.json"), "w"), indent=1)
    prov = dict(
        selection="all non-redundant mouse alpha-beta entries at <= %.1f A" % b["max_resolution"],
        non_redundancy="deduplicated on the exact (V-domain alpha, V-domain beta) sequence pair",
        species_assignment="source organism of the TCR chains themselves, not the entry",
        crystal_preparation=("numbering computed per chain on its own observed sequence "
                             "(fixes the SEQRES/observed register slip); alpha-beta pair "
                             "chosen by V-domain centroid proximity (fixes cross-copy "
                             "pairing in multi-copy asymmetric units)"),
        gates=dict(min_mapped_residues_per_chain=MIN_MAPPED,
                   min_site_positions_per_chain=MIN_SITE_PER_CHAIN,
                   max_centroid_separation_A=MAX_CENTROID_SEP,
                   max_framework_fit_A=MAX_FR_FIT),
        entries_attempted=int(len(B)), entries_accepted=int(len(ok)),
        exclusions=excl,
        resolution_of_accepted=dict(
            min=round(float(ok.resolution.min()), 2),
            median=round(float(ok.resolution.median()), 2),
            max=round(float(ok.resolution.max()), 2)),
        human_structures_measured=0,
        substitution_flag=("no substitution: mouse crystals only. No human structure "
                           "enters the floor."))
    json.dump(prov, open(os.path.join(out, "B3h_threshold_provenance.json"), "w"), indent=1)
    print("\naccepted %d of %d | exclusions %s" % (len(ok), len(B), excl))
    print("site RMSD  median %.3f  p75 %.3f  p90 %.3f  p%d %.3f  max %.3f"
          % (ok.site_rmsd.median(), *np.percentile(ok.site_rmsd, [75, 90]), p, val,
             ok.site_rmsd.max()))
    print("THRESHOLD %.2f A" % thr)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config_main.json")
