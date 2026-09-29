"""The model error floor in the metric the clustering actually thresholds.

The packaged pipeline has a units mismatch. cluster_sites.py clusters on SITE-DESCRIPTOR
RMSD: each receptor's site loops are resampled to n_resample points per loop (80x3 here)
after superposing that receptor on ONE global reference framework (ref_fr, the framework
of the first model in folding-set order, over the 166 anchor positions). But
benchmark_threshold.py measures RMSD over TRUE IMGT site positions in a frame shared by
just that crystal-model pair, and the pipeline then feeds that number in as the
descriptor cut. Those are different quantities: B3h's position-space p95 came out at
2.50 A while the repertoire's within-class DESCRIPTOR background has a median of 2.285 A,
so taking the position-space p95 as the descriptor cut would merge most of the
repertoire.

This measures the same crystal-model pairs in descriptor space:
  - fr_keys and ref_fr are rebuilt exactly as extract_sites.py does (min_frac=1.0 over
    all models, reference = first model in folding-set order), so the crystal, its
    model and all repertoire receptors sit in one frame;
  - the model descriptor uses all 166 anchors, as the repertoire descriptors do;
  - the crystal descriptor uses the anchors it actually has (crystals have disorder),
    superposed onto the corresponding rows of ref_fr; the count used is recorded;
  - the distance is the descriptor RMSD, the same quantity cluster_sites.py thresholds.

Crystal preparation (per-chain observed-sequence numbering, proximity-based alpha-beta
pairing) is imported from b3h_benchmark_fixed rather than duplicated.

Writes: tables/B3i_descriptor_floor.csv   per entry
        tables/threshold.json             descriptor-space, what clustering reads
        tables/B3i_threshold_provenance.json
"""
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import benchmark_threshold as BT  # noqa: E402
import imgt  # noqa: E402
import site_geometry as G  # noqa: E402
import b3h_benchmark_fixed as B3H  # noqa: E402

MIN_LOOP_RESIDUES = 4      # a resampled loop needs a real trace
MAX_LOOP_DEFICIT = 2       # crystal loop may be at most this many residues short


def descriptor_on_ref(ca, fr_keys, ref_fr, loops, n_res):
    """Descriptor in the global reference frame, using whichever anchors are present.

    Returns (descriptor, n_anchors_used, framework RMSD, {loop: n_residues}).
    """
    idx = [j for j, (ch, k) in enumerate(fr_keys) if k in ca[ch]]
    if len(idx) < 40:
        return None, len(idx), None, None
    F = np.array([ca[fr_keys[j][0]][fr_keys[j][1]] for j in idx])
    Rref = ref_fr[idx]
    R, cF, cR = G.kabsch(F, Rref)
    fit = (F - cF) @ R + cR
    fr_rmsd = float(np.sqrt(((fit - Rref) ** 2).sum(1).mean()))
    blocks, lens = [], {}
    for ch, nm in loops:
        keys = sorted([k for k in ca[ch] if k[0] in imgt.SITE_RANGES[nm]])
        if len(keys) < MIN_LOOP_RESIDUES:
            return None, len(idx), fr_rmsd, None
        pts = np.array([ca[ch][k] for k in keys])
        blocks.append(G.resample_loop((pts - cF) @ R + cR, n_res))
        lens["%s_%s" % (ch, nm)] = len(keys)
    return np.concatenate(blocks), len(idx), fr_rmsd, lens


def main(config_path="config_main.json"):
    cfg = imgt.load(config_path)
    out, b = cfg["out_dir"], cfg["benchmark"]
    chains = imgt.chain_labels(cfg)
    types = {c: imgt.SHARED_CHAIN_TYPES.get(c, (c,)) for c in chains}
    loops, n_res = imgt.site_loops(cfg), cfg["clustering"]["n_resample"]
    mdir, work = cfg["folding"]["model_dir"], os.path.join(out, "benchmark")

    # --- rebuild the repertoire frame exactly as extract_sites.py does
    fs = pd.read_csv(os.path.join(out, "folding_set.csv.gz"), low_memory=False)
    ids = [c for c in fs.clone_id if os.path.exists(os.path.join(mdir, "%s.pdb" % c))]
    print("reading %d repertoire models to rebuild the frame" % len(ids), flush=True)
    ca_rep = {i: G.read_ca(os.path.join(mdir, "%s.pdb" % i), chains) for i in ids}
    fr_keys = G.framework_positions(list(ca_rep.values()), chains, min_frac=1.0)
    ref_fr = G.framework_matrix(ca_rep[ids[0]], fr_keys)
    print("anchors %d | reference model %s" % (len(fr_keys), ids[0]), flush=True)

    # sanity: reproduce a repertoire descriptor and check it matches the stored one
    z = np.load(os.path.join(out, "combining_sites.npz"), allow_pickle=True)
    stored = dict(zip([str(x) for x in z["clone_id"]], z["descriptor"]))
    d0, _, fr0, _ = descriptor_on_ref(ca_rep[ids[0]], fr_keys, ref_fr, loops, n_res)
    dev = float(np.abs(d0 - stored[str(ids[0])]).max())
    print("frame check: max deviation from stored descriptor %.2e A (fr_rmsd %.3f)"
          % (dev, fr0), flush=True)
    assert dev < 1e-3, "rebuilt frame does not reproduce the stored descriptors"

    acc = pd.read_csv(os.path.join(out, "B3h_benchmark_structures.csv"))
    acc = acc[acc.accepted.astype(bool)]
    ents = BT.entity_sequences(list(acc.entry))
    rows = []
    for eid in acc.entry:
        rec = dict(entry=eid)
        try:
            res = BT.read_pdb_chain_residues(os.path.join(work, "%s.pdb" % eid))
            pick, sep = B3H.prepare(eid, ents.get(eid) or [], res, cfg, chains, types)
            vseqs = {c: pick[c]["vseq"] for c in chains}
            tag = hashlib.sha1("|".join(vseqs[c] for c in chains).encode()).hexdigest()[:10]
            cryst = os.path.join(work, "fixed_%s_crystal.pdb" % eid)
            model = os.path.join(work, "seq_%s_model.pdb" % tag)
            ca_c, ca_m = G.read_ca(cryst, chains), G.read_ca(model, chains)
            dc, na_c, frc, lc = descriptor_on_ref(ca_c, fr_keys, ref_fr, loops, n_res)
            dm, na_m, frm, lm = descriptor_on_ref(ca_m, fr_keys, ref_fr, loops, n_res)
            if dc is None or dm is None:
                raise ValueError("descriptor unavailable (anchors %d/%d, loops %s)"
                                 % (na_c, len(fr_keys), lc is not None))
            deficit = max(lm[k] - lc[k] for k in lm)
            rec.update(desc_rmsd=float(np.sqrt(((dc - dm) ** 2).sum(1).mean())),
                       anchors_crystal=na_c, anchors_model=na_m,
                       fr_rmsd_crystal=round(frc, 3), fr_rmsd_model=round(frm, 3),
                       max_loop_deficit=int(deficit), centroid_sep=sep,
                       **{"loop_%s" % k: lc[k] for k in sorted(lc)})
            print("  %s  descriptor RMSD %.2f A | anchors %d/%d | loop deficit %d"
                  % (eid, rec["desc_rmsd"], na_c, len(fr_keys), deficit), flush=True)
        except Exception as e:
            rec.update(error="%s: %s" % (type(e).__name__, str(e)[:140]))
            print("  %s FAILED %s: %s" % (eid, type(e).__name__, str(e)[:90]), flush=True)
        rows.append(rec)

    B = pd.DataFrame(rows)
    gate = B.desc_rmsd.notna() & (B.max_loop_deficit <= MAX_LOOP_DEFICIT)
    B["accepted"] = gate
    B.to_csv(os.path.join(out, "B3i_descriptor_floor.csv"), index=False)
    ok = B[gate]
    if len(ok) < 5:
        raise SystemExit("only %d structures usable in descriptor space" % len(ok))

    p = b["percentile"]
    val = float(np.percentile(ok.desc_rmsd, p))
    thr = max(round(val * 20) / 20.0, b["min_threshold"])
    pos = pd.read_csv(os.path.join(out, "B3c_extraction_summary.csv"))
    getq = lambda q: float(pos.loc[pos.quantity == q, "value"].iloc[0])  # noqa: E731
    rec = dict(threshold=thr,
               estimator="p%d of model-vs-crystal SITE-DESCRIPTOR RMSD" % p,
               metric="site-descriptor RMSD, the quantity cluster_sites.py thresholds",
               percentile=p, raw_percentile_value=round(val, 3),
               median=round(float(ok.desc_rmsd.median()), 3),
               twice_median=round(2 * float(ok.desc_rmsd.median()), 3),
               n_structures=int(len(ok)), n_attempted=int(len(B)),
               min_threshold_floor=b["min_threshold"],
               max_resolution=b["max_resolution"], species=cfg["species"],
               locus_pair=cfg["locus_pair"], builder=cfg["folding"]["builder"],
               entries=sorted(ok.entry))
    json.dump(rec, open(os.path.join(out, "threshold.json"), "w"), indent=1)
    prov = dict(
        metric_note=("cluster_sites.py thresholds site-descriptor RMSD; "
                     "benchmark_threshold.py measures RMSD at true IMGT site positions "
                     "in a pair-local frame. Those are different quantities and the "
                     "packaged pipeline feeds the second into the first. This threshold "
                     "is measured in descriptor space, so cut and distances share units."),
        position_space_p95_for_reference=2.50,
        selection="all non-redundant mouse alpha-beta entries at <= %.1f A" % b["max_resolution"],
        crystal_preparation=("per-chain observed-sequence numbering; alpha-beta pair by "
                             "V-domain centroid proximity (see b3h_benchmark_fixed.py)"),
        gates=dict(inherited_from_B3h=True, max_loop_deficit_residues=MAX_LOOP_DEFICIT,
                   min_loop_residues=MIN_LOOP_RESIDUES),
        entries_attempted=int(len(B)), entries_accepted=int(len(ok)),
        crystal_anchors_used=dict(min=int(ok.anchors_crystal.min()),
                                  median=float(ok.anchors_crystal.median()),
                                  of=len(fr_keys)),
        repertoire_reference=dict(anchors=len(fr_keys), reference_model=str(ids[0])),
        comparison=dict(
            near_identical_pair_median=getq("near-identical pair site RMSD, median (A)"),
            near_identical_pair_p90=getq("near-identical pair site RMSD, 90th pct (A)"),
            background_within_class_median=getq("background within-class RMSD, median (A)"),
            background_within_class_p5=getq("background within-class RMSD, 5th pct (A)")),
        human_structures_measured=0,
        substitution_flag="no substitution: mouse crystals only.")
    json.dump(prov, open(os.path.join(out, "B3i_threshold_provenance.json"), "w"), indent=1)

    print("\nDESCRIPTOR-SPACE FLOOR  n=%d  median %.3f  p75 %.3f  p90 %.3f  p%d %.3f  max %.3f"
          % (len(ok), ok.desc_rmsd.median(), *np.percentile(ok.desc_rmsd, [75, 90]), p,
             val, ok.desc_rmsd.max()))
    print("THRESHOLD %.2f A" % thr)
    print("for context: near-identical pairs median %.3f p90 %.3f | background median "
          "%.3f p5 %.3f"
          % (getq("near-identical pair site RMSD, median (A)"),
             getq("near-identical pair site RMSD, 90th pct (A)"),
             getq("background within-class RMSD, median (A)"),
             getq("background within-class RMSD, 5th pct (A)")))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config_main.json")
