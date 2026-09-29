"""VXH7 (the voxel pipeline, A12.3): confounds for the surviving arm (Arm B, grid alone). No state data.

VX6 items 1-2 on the reference cluster-test procedure's 60,000 background pairs: Spearman of distance with |delta CDR3 length| (alpha, beta, total);
median distance for length-matched vs mismatched pairs; median distance for same- vs different-V-pair; share of clusters
that are single-V-pair. ARI of the Arm B partition against <reference arm> (only clone_id and label columns read), the
whole-molecule grid (VX4) and the F0 loops grid (VX6L). Writes out/VXH7_diagnostics.csv, out/VXH7_summary.json;
checks/VXH7_checks.csv.
usage: python pipeline/code/vxh7_confounds.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column

from vx4_cluster import ari

if __name__ == "__main__":
    C = Checks("VXH7", stop_on_fail=True)
    C.info("status", "A12 VXH7 confounds; no state data read")
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; n = len(ids); idx = {c: i for i, c in enumerate(ids)}
    DB = np.load(os.path.join(TMP, "VXH5_D_receptors.f32.npy"), mmap_mode="r")
    DW = np.load(os.path.join(TMP, "VX3_D_receptors.f32.npy"), mmap_mode="r")
    DL = np.load(os.path.join(TMP, "VX6L_D_receptors.f32.npy"), mmap_mode="r")
    rng = np.random.default_rng(0); ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
    ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz")).set_index("clone_id").reindex(ids)
    vpair = (ID.v_A_prot.astype(str) + "|" + ID.v_B_prot.astype(str)).values
    la, lb = R.cdr3len_A.values, R.cdr3len_B.values
    dA, dB = np.abs(la[ba] - la[bb]), np.abs(lb[ba] - lb[bb]); sameV = vpair[ba] == vpair[bb]
    LB = pd.read_csv(os.path.join(OUT, "VXH6_molecule_labels.csv.gz"))
    L4 = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"), usecols=["clone_id", "prop_voxel"])
    L6 = pd.read_csv(os.path.join(OUT, "VX6L_molecule_labels.csv.gz"), usecols=["clone_id", "prop_voxel_loops"])
    LE = pd.read_csv(VXP(_N10F), usecols=["clone_id", _N10C])
    M = LB.merge(L4, on="clone_id").merge(L6, on="clone_id").merge(LE, on="clone_id", how="left")
    C.add("partitions aligned on the {:,} molecules".format(EXP_NMOL), len(M), len(M) == EXP_NMOL, "== %s" % EXP_NMOL)
    M[_N10C] = M[_N10C].fillna(-1).astype(int)
    rep = np.array([idx[c] for c in M.clone_id])
    rows = []
    for name, Dx, lab in (("Arm B: F1 loops, 7 ch", DB, M.prop_voxel_armB.values), ("whole molecule (VX4)", DW, M.prop_voxel.values),
                          ("F0 loops (VX6L)", DL, M.prop_voxel_loops.values)):
        d = np.asarray(Dx[ba, bb], dtype=np.float64)
        o = dict(arm=name)
        for k, x in (("alpha", dA), ("beta", dB), ("total", dA + dB)):
            o["spearman_dlen_" + k] = float(spearmanr(d, x).correlation)
        m0 = (dA + dB) == 0
        o["median_len_matched"], o["median_len_mismatched"] = float(np.median(d[m0])), float(np.median(d[~m0]))
        o["median_same_Vpair"], o["median_diff_Vpair"] = float(np.median(d[sameV])), float(np.median(d[~sameV]))
        vp = pd.Series(vpair[rep]); cl = pd.Series(lab); g = vp[cl >= 0].groupby(cl[cl >= 0]).nunique()
        o["clusters"] = int(len(g)); o["share_single_Vpair_clusters"] = float((g == 1).mean())
        o["ARI_vs_reference"] = ari(lab, M[_N10C].values)
        rows.append(o)
    le = M[_N10C].values; ge = pd.Series(vpair[rep])[le >= 0].groupby(pd.Series(le)[le >= 0]).nunique()
    rows.append(dict(arm="reference method (context)", clusters=int(len(ge)), share_single_Vpair_clusters=float((ge == 1).mean())))
    DG = pd.DataFrame(rows); DG.to_csv(os.path.join(OUT, "VXH7_diagnostics.csv"), index=False)
    for _, o in DG.iterrows():
        C.info("diagnostics: %s" % o.arm, "; ".join("%s=%.4g" % (k, v) for k, v in o.items() if k != "arm" and pd.notna(v)))
    C.info("same-V-pair background pairs", "%d of %d" % (int(sameV.sum()), len(sameV)))
    C.info("ARI Arm B vs whole-molecule grid / vs F0 loops grid", "%.4f / %.4f"
           % (ari(M.prop_voxel_armB.values, M.prop_voxel.values), ari(M.prop_voxel_armB.values, M.prop_voxel_loops.values)))
    save_json(dict(rows=DG.to_dict("records")), os.path.join(OUT, "VXH7_summary.json"))
    C.write()
