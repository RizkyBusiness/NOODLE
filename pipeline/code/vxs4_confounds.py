"""VXS4 (the voxel pipeline, A14.2.4): confounds for the sigma 1.5 A arms D and E beside Arms B and C and <reference arm>. No state data.

As VXH7 / VXC3, on the reference cluster-test procedure's 60,000 background pairs (default_rng(0)): Spearman of distance with |delta CDR3 length| (alpha,
beta, total); length-matched vs mismatched medians; same- vs different-V-pair medians; share of single-V-pair clusters;
ARI with <reference arm> (only clone_id and label columns read); ARI D vs B, E vs C, D vs E.
Check: the Arm B and Arm C rows reproduce VXC3_diagnostics.csv (1e-12). Writes out/VXS4_diagnostics.csv,
out/VXS4_summary.json; checks/VXS4_checks.csv.
usage: python pipeline/code/vxs4_confounds.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys
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
    C = Checks("VXS4", stop_on_fail=True)
    C.info("status", "A14 VXS4 confounds; no state data read")
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; n = len(ids); idx = {c: i for i, c in enumerate(ids)}
    rng = np.random.default_rng(0); ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
    ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz")).set_index("clone_id").reindex(ids)
    vpair = (ID.v_A_prot.astype(str) + "|" + ID.v_B_prot.astype(str)).values
    la, lb = R.cdr3len_A.values, R.cdr3len_B.values
    L = {k: pd.read_csv(os.path.join(OUT, f), usecols=["clone_id", c]) for k, f, c in (
        ("D", "VXS2_D_molecule_labels.csv.gz", "prop_voxel_armD"), ("E", "VXS2_E_molecule_labels.csv.gz", "prop_voxel_armE"),
        ("B", "VXH6_molecule_labels.csv.gz", "prop_voxel_armB"), ("C", "VXC2_molecule_labels.csv.gz", "prop_voxel_armC"))}
    M = L["D"]
    for k in "EBC":
        M = M.merge(L[k], on="clone_id")
    LE = pd.read_csv(VXP(_N10F), usecols=["clone_id", _N10C])
    M = M.merge(LE, on="clone_id", how="left"); M[_N10C] = M[_N10C].fillna(-1).astype(int)
    C.add("partitions aligned on the {:,} molecules".format(EXP_NMOL), len(M), len(M) == EXP_NMOL, "== %s" % EXP_NMOL)
    rep = np.array([idx[c] for c in M.clone_id])

    def diag(name, Dx, lab):
        d = np.asarray(Dx[ba, bb], dtype=np.float64); dA, dB = np.abs(la[ba] - la[bb]), np.abs(lb[ba] - lb[bb]); sv = vpair[ba] == vpair[bb]
        o = dict(arm=name, pairs=len(d))
        for k, x in (("alpha", dA), ("beta", dB), ("total", dA + dB)):
            o["spearman_dlen_" + k] = float(spearmanr(d, x).correlation)
        m0 = (dA + dB) == 0
        o["median_len_matched"], o["median_len_mismatched"] = float(np.median(d[m0])), float(np.median(d[~m0]))
        o["median_same_Vpair"], o["median_diff_Vpair"] = float(np.median(d[sv])), float(np.median(d[~sv]))
        vp = pd.Series(vpair[rep]); cl = pd.Series(lab); g = vp[cl >= 0].groupby(cl[cl >= 0]).nunique()
        o["clusters"] = int(len(g)); o["share_single_Vpair_clusters"] = float((g == 1).mean())
        o["ARI_vs_reference"] = ari(lab, M[_N10C].values)
        return o

    rows = []
    for name, f, col in (("Arm D: F1 loops, 7 ch, sigma 1.5", "VXS1_D_D_receptors.f32.npy", "prop_voxel_armD"),
                         ("Arm E: F1 CDR3, 7 ch, sigma 1.5", "VXS1_E_D_receptors.f32.npy", "prop_voxel_armE"),
                         ("Arm B: F1 loops, 7 ch", "VXH5_D_receptors.f32.npy", "prop_voxel_armB"),
                         ("Arm C: F1 CDR3, 7 ch", "VXC1_D_receptors.f32.npy", "prop_voxel_armC")):
        rows.append(diag(name, np.load(os.path.join(TMP, f), mmap_mode="r"), M[col].values))
    le = M[_N10C].values; ge = pd.Series(vpair[rep])[le >= 0].groupby(pd.Series(le)[le >= 0]).nunique()
    rows.append(dict(arm="reference method (context)", clusters=int(len(ge)), share_single_Vpair_clusters=float((ge == 1).mean())))
    DG = pd.DataFrame(rows)
    ref = pd.read_csv(os.path.join(OUT, "VXC3_diagnostics.csv")).set_index("arm")
    num = [c for c in DG.columns if c not in ("arm",)]
    e = max(abs(float(DG.set_index("arm").loc[a, c]) - float(ref.loc[a, c])) for a in ("Arm B: F1 loops, 7 ch", "Arm C: F1 CDR3, 7 ch") for c in num)
    C.add("Arm B and Arm C rows reproduce VXC3_diagnostics.csv", "max |diff| %.1e" % e, e < 1e-12, "< 1e-12")
    aris = {}
    for a_, b_ in (("prop_voxel_armD", "prop_voxel_armB"), ("prop_voxel_armE", "prop_voxel_armC"), ("prop_voxel_armD", "prop_voxel_armE")):
        aris["%s_vs_%s" % (a_[-1], b_[-1])] = ari(M[a_].values, M[b_].values); C.info("ARI %s vs %s" % (a_, b_), "%.4f" % aris["%s_vs_%s" % (a_[-1], b_[-1])])
    DG.to_csv(os.path.join(OUT, "VXS4_diagnostics.csv"), index=False)
    for _, o in DG.iterrows():
        C.info("diagnostics: %s" % o.arm, "; ".join("%s=%.4g" % (k, v) for k, v in o.items() if k != "arm" and pd.notna(v)))
    save_json(dict(rows=DG.to_dict("records"), ari=aris), os.path.join(OUT, "VXS4_summary.json"))
    C.write()
