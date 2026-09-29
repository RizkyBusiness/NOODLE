"""VXV0 (the voxel pipeline, A12.1): audit for the descriptor validity panel.

The accepted pairs (vxlib.load_lm3) -> their receptors (distinct models); where a receptor has several crystals, the
highest-resolution one (B3h resolution; ties: first PDB id), fixed here before any distance. Per chosen crystal:
shared atoms with its model, positions missing in the crystal, complete-CDR3 check (every CDR3 Ca of the model present
in the crystal, as the vector arcs need), species of the TCR chains and release date (RCSB, public ids only; A12.9.1),
bound / unbound (a contacting peptide of 5-30 residues within 5 A of the TCR, local raw PDB, as H15). Disk.
Writes out/VXV0_receptors.csv, out/VXV0_audit.json, checks/VXV0_checks.csv. No state data.
usage: python pipeline/code/vxv0_audit.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, shutil, urllib.request
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NBP = expected("n_benchmark_pairs")
EXP_NBM = expected("n_benchmark_models")
from vxv_common import BEN, missing_positions, shared_keys, residues

if __name__ == "__main__":
    C = Checks("VXV0", stop_on_fail=True)
    C.info("status", "A12 VXV0 audit; no state data read")
    lm3 = load_lm3(); B = lm3["B"].copy()
    C.add("lm3 accepted pairs", len(B), len(B) == EXP_NBP, "== %s" % EXP_NBP)
    B3h = pd.read_csv(os.path.join(STB, "B3h_benchmark_structures.csv")).set_index("entry")
    B["resolution"] = B3h.reindex(B.entry).resolution.values
    B["auth_A"] = B3h.reindex(B.entry).auth_A.values; B["auth_B"] = B3h.reindex(B.entry).auth_B.values
    C.add("resolution known for all %s crystals" % EXP_NBP, int(B.resolution.notna().sum()), B.resolution.notna().all(), "%s" % EXP_NBP)
    ch = (B.sort_values(["model", "resolution", "entry"]).groupby("model", as_index=False).first())
    ch["n_crystals"] = ch.model.map(B.model.value_counts())
    ch["all_entries"] = ch.model.map(B.groupby("model").entry.apply(lambda s: " ".join(sorted(s))))
    C.add("receptors (distinct models)", len(ch), len(ch) == EXP_NBM, "== %s" % EXP_NBM)

    rows = []
    for _, r in ch.iterrows():
        pc, pm = os.path.join(BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(BEN, r.model)
        mc, mm = parse_heavy(pc)[1], parse_heavy(pm)[1]
        miss = missing_positions(pc, pm)
        cdr3_m = {(a[0], a[1], a[2]) for a in mm if a[4] == "CA" and 105 <= a[1] <= 117}
        cdr3_c = {(a[0], a[1], a[2]) for a in mc if a[4] == "CA" and 105 <= a[1] <= 117}
        rows.append(dict(model=r.model, entry=r.entry, resolution=r.resolution, identity=r.identity,
                         n_crystals=int(r.n_crystals), all_entries=r.all_entries, auth_A=r.auth_A, auth_B=r.auth_B,
                         model_atoms=len(mm), crystal_atoms=len(mc), shared_atoms=len(shared_keys(pc, pm)),
                         missing_positions=len(miss), missing_list=" ".join("%s%d%s" % m for m in miss),
                         cdr3_complete=cdr3_m <= cdr3_c))
    R = pd.DataFrame(rows)
    C.info("receptors with several crystals", "%d (crystal = highest resolution)" % int((R.n_crystals > 1).sum()))
    C.info("crystals with missing residues / total missing positions",
           "%d / %d" % (int((R.missing_positions > 0).sum()), int(R.missing_positions.sum())))
    C.info("complete CDR3 (vector arcs usable)", "%d / %s; failing: %s" % (int(R.cdr3_complete.sum()), EXP_NBM,
           " ".join(R.entry[~R.cdr3_complete]) or "none"))
    V3 = np.load(VXP("reference/out/V3_crystal_coords.npz"), allow_pickle=True)
    C.add("every chosen crystal has V3 arcs", "", set(R.entry) <= set(V3["entry"]), "subset")

    # bound / unbound from the local raw PDB (H15 rule: a 5-30 residue chain within 5 A of the TCR)
    bound = []
    for _, r in R.iterrows():
        raw = os.path.join(BEN, "%s.pdb" % r.entry)
        at = [(L[21], L[17:20], L[12:16].strip(), float(L[30:38]), float(L[38:46]), float(L[46:54]))
              for L in open(raw) if L.startswith("ATOM") and L[16] in " A"]
        df = pd.DataFrame(at, columns=["ch", "res", "atom", "x", "y", "z"])
        tcr = df[df.ch.isin([r.auth_A, r.auth_B])][["x", "y", "z"]].values
        pep = False
        for c_, g in df[~df.ch.isin([r.auth_A, r.auth_B])].groupby("ch"):
            nres = int((g.atom == "CA").sum())
            if 5 <= nres <= 30:
                xyz = g[["x", "y", "z"]].values
                dmin = min(np.sqrt(((tcr[i:i + 400, None] - xyz[None]) ** 2).sum(2)).min() for i in range(0, len(tcr), 400))
                pep |= dmin <= 5.0
        bound.append(pep)
    R["bound_pMHC"] = bound
    C.info("bound (contacting peptide) / unbound", "%d / %d" % (int(R.bound_pMHC.sum()), int((~R.bound_pMHC).sum())))

    # species and release date from RCSB (public PDB ids only)
    q = """{ entries(entry_ids: [%s]) { rcsb_id rcsb_accession_info { initial_release_date }
             polymer_entities { rcsb_polymer_entity_container_identifiers { auth_asym_ids }
                                rcsb_entity_source_organism { scientific_name } } } }""" % ",".join('"%s"' % e for e in R.entry)
    species, rel, net = {}, {}, "ok"
    try:
        req = urllib.request.Request("https://data.rcsb.org/graphql", data=json.dumps({"query": q}).encode(),
                                     headers={"Content-Type": "application/json"})
        js = json.load(urllib.request.urlopen(req, timeout=60))
        for e in js["data"]["entries"]:
            rel[e["rcsb_id"]] = e["rcsb_accession_info"]["initial_release_date"][:10]
            r = R[R.entry == e["rcsb_id"]].iloc[0]
            for pe in e["polymer_entities"]:
                ids_ = pe["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] or []
                if r.auth_A in ids_ or r.auth_B in ids_:
                    sp = {o["scientific_name"] for o in (pe["rcsb_entity_source_organism"] or [])}
                    species.setdefault(e["rcsb_id"], set()).update(sp)
    except Exception as ex:
        net = "failed: %s" % ex
    C.info("RCSB metadata query", net)
    R["species"] = [" + ".join(sorted(species.get(e, {"unknown"}))) for e in R.entry]
    R["release_date"] = [rel.get(e, "") for e in R.entry]
    C.info("species of TCR chains (chosen crystals)", "; ".join("%s: %d" % kv for kv in R.species.value_counts().items()))
    C.add("species known for every chosen crystal", int((R.species != "unknown").sum()), (R.species != "unknown").all(),
          "%s (strata are report-only)" % EXP_NBM)
    R.to_csv(os.path.join(OUT, "VXV0_receptors.csv"), index=False)

    du = shutil.disk_usage(TMP)
    C.info("free disk (project volume)", "%.1f GB" % (du.free / 1e9))
    save_json(dict(n_receptors=len(R), crystal_rule="highest resolution; ties first PDB id",
                   cdr3_incomplete=list(R.entry[~R.cdr3_complete]), rcsb=net, free_disk_GB=round(du.free / 1e9, 1)),
              os.path.join(OUT, "VXV0_audit.json"))
    C.write()
