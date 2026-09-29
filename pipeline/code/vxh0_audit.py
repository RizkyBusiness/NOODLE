"""VXH0 (the voxel pipeline, A12.3): audit for the hinge module.

Inputs (the landmark file anchors, D3, <reference arm> thresholds, controls; vector features are rebuilt from the landmark file as in VXV1, which
reproduced the V6 w = 0.5 row, because the reference method/out/V1_features.npz is not present); disk (Arm B ~39 GB, A12.10.1;
nothing in tmp/ is deleted); the separate STCRpy env (A12.11) and a TRangle test; anchors per chain for every benchmark
structure (fits need >= 50); leakage: release date, species of the TCR chains (RCSB) and bound / unbound (local raw
PDB, as VXV0) for all benchmark crystals. TCRBuilder2+ training cutoff: not stated in Quast et al. 2025 (training = TCRBuilder2's
704 STCRDab structures + 204 Immunocore; test = 45 newer non-identical structures) and TCRBuilder2's own data are in the
ImmuneBuilder supplement (Abanades et al. 2023, not accessible); the antibody model there used SAbDab to 31 Jul 2021.
Per A12.3, leakage is therefore reported by release year only.
Writes out/VXH0_crystals.csv, out/VXH0_audit.json; checks/VXH0_checks.csv. No state data.
usage: python pipeline/code/vxh0_audit.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, shutil, subprocess, urllib.request
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NANC = expected("n_anchors")
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column

import vx5b_frame_diagnostic as FB

STPY = __import__("vxpaths").env_python("stcrpy")   # config "envs"

if __name__ == "__main__":
    C = Checks("VXH0", stop_on_fail=True)
    C.info("status", "A12 VXH0 audit; no state data read")
    for lg in ("landmarks/out/LM1_internal_coords.npz", "descriptors/out/D3_property_descriptor.npz", _REF2.get("n10_thresholds", ""),
               "tables/B3c_near_identical_pairs.csv.gz", _N10F,
               "reference/out/V6_orientation_weight_grid.csv", "voxel_out/out/VXV1_descriptors.npz"):
        f = VXP(lg)                                   # check name: the logical path (config "paths")
        C.add("input exists: %s" % lg, os.path.exists(f), os.path.exists(f), "exists")
    C.info("reference/out/V1_features.npz", "absent; vector features rebuilt from LM1 (VXV1 reproduced V6 w=0.5 to 5e-5)")
    L = np.load(LM1, allow_pickle=True)
    C.add("LM1 anch shape", str(L["anch"].shape), L["anch"].shape == (EXP_NS, EXP_NANC, 3), "(%s,%s,3)" % (EXP_NS, EXP_NANC))
    du = shutil.disk_usage(TMP)
    C.add("free disk for Arm B (~39 GB)", "%.1f GB" % (du.free / 1e9), du.free / 1e9 > 60, "> 60 GB (39 + margin)")
    C.info("tmp/ contents kept (nothing deleted)", "%.1f GB" % (sum(os.path.getsize(os.path.join(TMP, f)) for f in os.listdir(TMP)) / 1e9))

    # STCRpy in its own env
    code = ("import warnings; warnings.filterwarnings('ignore'); import stcrpy, numpy; "
            "t=stcrpy.load_TCR(%r); t=t[0] if isinstance(t,list) else t; "
            "a=t.get_TCR_angles(); print(numpy.__version__, *['%%s=%%.3f'%%kv for kv in a.items()])" % VXP("structures/%s.pdb" % REF_ID))
    r = subprocess.run([STPY, "-c", code], capture_output=True, text=True)
    ok = r.returncode == 0 and "BA=" in r.stdout
    C.add("STCRpy env computes TRangle on the reference", r.stdout.strip().splitlines()[-1] if ok else r.stderr[-300:], ok, "runs")

    # benchmark: anchors per chain, bound, species, release
    lm3 = load_lm3(); B = lm3["B"].sort_values("entry").reset_index(drop=True)
    B3h = pd.read_csv(os.path.join(STB, "B3h_benchmark_structures.csv")).set_index("entry")
    FB.PA.lm3_LMPOS = lm3["LMPOS"]
    ia = np.array([k[0] == "A" for k in FB.AKEYS])
    rows = []
    for _, r_ in B.iterrows():
        pc, pm = os.path.join(FB.PA.BEN, "fixed_%s_crystal.pdb" % r_.entry), os.path.join(FB.PA.BEN, r_.model)
        okc, okm = FB.bench_points(pc, lm3["parse"])[2], FB.bench_points(pm, lm3["parse"])[2]
        sh = okc & okm
        aA, aB = B3h.loc[r_.entry, "auth_A"], B3h.loc[r_.entry, "auth_B"]
        at = [(L_[21], L_[12:16].strip(), float(L_[30:38]), float(L_[38:46]), float(L_[46:54]))
              for L_ in open(os.path.join(FB.PA.BEN, "%s.pdb" % r_.entry)) if L_.startswith("ATOM") and L_[16] in " A"]
        df = pd.DataFrame(at, columns=["ch", "atom", "x", "y", "z"])
        tcr = df[df.ch.isin([aA, aB])][["x", "y", "z"]].values
        pep = False
        for _, g in df[~df.ch.isin([aA, aB])].groupby("ch"):
            if 5 <= int((g.atom == "CA").sum()) <= 30:
                xyz = g[["x", "y", "z"]].values
                pep |= min(np.sqrt(((tcr[i:i + 400, None] - xyz[None]) ** 2).sum(2)).min() for i in range(0, len(tcr), 400)) <= 5.0
        rows.append(dict(entry=r_.entry, model=r_.model, resolution=B3h.loc[r_.entry, "resolution"], auth_A=aA, auth_B=aB,
                         anchors_crystal_A=int((okc & ia).sum()), anchors_crystal_B=int((okc & ~ia).sum()),
                         anchors_shared_A=int((sh & ia).sum()), anchors_shared_B=int((sh & ~ia).sum()), bound_pMHC=pep))
    X = pd.DataFrame(rows)
    mn = int(X[["anchors_crystal_A", "anchors_crystal_B", "anchors_shared_A", "anchors_shared_B"]].min().min())
    C.add("anchors per chain per fit, all %s pairs: min" % EXP_NBP, mn, mn >= 50, ">= 50")
    q = """{ entries(entry_ids: [%s]) { rcsb_id rcsb_accession_info { initial_release_date }
             polymer_entities { rcsb_polymer_entity_container_identifiers { auth_asym_ids }
                                rcsb_entity_source_organism { scientific_name } } } }""" % ",".join('"%s"' % e for e in X.entry)
    rel, spc, net = {}, {}, "ok"
    try:
        req = urllib.request.Request("https://data.rcsb.org/graphql", data=json.dumps({"query": q}).encode(),
                                     headers={"Content-Type": "application/json"})
        js = json.load(urllib.request.urlopen(req, timeout=60))
        for e in js["data"]["entries"]:
            rel[e["rcsb_id"]] = e["rcsb_accession_info"]["initial_release_date"][:10]
            r_ = X[X.entry == e["rcsb_id"]].iloc[0]
            for pe in e["polymer_entities"]:
                ids_ = pe["rcsb_polymer_entity_container_identifiers"]["auth_asym_ids"] or []
                if r_.auth_A in ids_ or r_.auth_B in ids_:
                    spc.setdefault(e["rcsb_id"], set()).update(o["scientific_name"] for o in (pe["rcsb_entity_source_organism"] or []))
    except Exception as ex:
        net = "failed: %s" % ex
    C.add("RCSB metadata for the %s crystals" % EXP_NBP, net, net == "ok" and len(rel) == EXP_NBP, "%s entries" % EXP_NBP)
    X["release_date"] = X.entry.map(rel); X["species"] = [" + ".join(sorted(spc.get(e, {"unknown"}))) for e in X.entry]
    X["release_year"] = X.release_date.str[:4].astype(int)
    X["year_bin"] = pd.cut(X.release_year, [0, 2015, 2019, 2100], labels=["<=2015", "2016-2019", ">=2020"]).astype(str)
    C.info("release years (bins)", "; ".join("%s: %d" % kv for kv in X.year_bin.value_counts().sort_index().items()))
    C.info("release year range", "%d-%d" % (X.release_year.min(), X.release_year.max()))
    C.info("bound / unbound (%s)" % EXP_NBP, "%d / %d" % (int(X.bound_pMHC.sum()), int((~X.bound_pMHC).sum())))
    C.info("species (%s)" % EXP_NBP, "; ".join("%s: %d" % kv for kv in X.species.value_counts().items()))
    C.info("TCRBuilder2+ training cutoff", "not established (see docstring); leakage reported by release year only")
    X.to_csv(os.path.join(OUT, "VXH0_crystals.csv"), index=False)
    save_json(dict(free_disk_GB=round(du.free / 1e9, 1), stcrpy_env=STPY, rcsb=net,
                   training_cutoff="not established; ABodyBuilder2 used SAbDab to 2021-07-31 (context only)",
                   release_years=[int(X.release_year.min()), int(X.release_year.max())]),
              os.path.join(OUT, "VXH0_audit.json"))
    C.write()
