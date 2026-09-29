"""VXH9 (the voxel pipeline, A12.3): the comparison table against the existing arms, assembled from the step outputs only (every
number is read from a file in <voxel out>/out/ or the reference method/), and the scoring of the recorded predictions. No new analysis.
Writes out/VXH9_comparison.csv, out/VXH9_predictions.csv; checks/VXH9_checks.csv.
usage: python pipeline/code/vxh9_report.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
_NBP = expected("n_benchmark_pairs")   # number of benchmark crystal-model pairs
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column
_N10NAME = _REF2.get("n10_name", "reference method (secondary)")   # its row label in the comparison table

from vxpaths import CFG as _CFG
_REF = _CFG.get("reference_method", {})
_refp = lambda k, d="": _REF.get(k, d)
_RD = _refp("dir", "reference/out/")
_RNAME = _refp("name", "reference method")


J = lambda f: json.load(open(os.path.join(OUT, f)))
if __name__ == "__main__":
    C = Checks("VXH9", stop_on_fail=True)
    E = pd.read_csv(os.path.join(OUT, "VXV2_endpoints.csv")).set_index("descriptor")
    V3 = pd.read_csv(os.path.join(OUT, "VXV3_cdr3.csv")).set_index("descriptor")
    DG = pd.read_csv(os.path.join(OUT, "VXC3_diagnostics.csv")).set_index("arm")
    S8 = pd.read_csv(os.path.join(OUT, "VXH8_arm_comparison.csv")).set_index("arm")
    V4s, V5s, H6, C2 = J("VX4_summary.json"), J("VX5_summary.json"), J("VXH6_summary.json"), J("VXC2_summary.json")
    T0 = pd.read_csv(VXP(_REF2.get("n10_thresholds", ""))).set_index("metric").loc[_N10C]
    V6 = pd.read_csv(VXP("reference/out/V6_orientation_weight_grid.csv")).set_index("w").loc[0.5]
    AC = pd.read_csv(VXP(_REF2.get("n10_arm_comparison", ""))).iloc[0]
    KT = pd.read_csv(VXP(_REF2.get("n10_cluster_tests", "")))
    LE = pd.read_csv(VXP(_N10F), usecols=[_N10C])[_N10C]
    A5 = pd.read_csv(VXP(_RD + _refp("survivor_sets", {}).get("summary", "A5_summary.csv"))).set_index("arm").loc[_refp("label_column")]
    F5 = pd.read_csv(os.path.join(OUT, "VX5_crystal_floor.csv"))
    rows = [
        dict(arm="voxel whole molecule (VX1-VX5)", scope="whole V-domain pair, F0 frame, 7 ch", cut=V4s["cut"],
             control_within_cut=V4s["control_within_cut"], crystal_err_over_cut=V5s["ratio_shared"],
             crystal_within_cut=int((F5.over_cut_shared <= 1).sum()), clusters=V4s["clusters"], clustered=V4s["clustered"],
             single_Vpair_clusters=DG.loc["whole molecule (VX4)", "share_single_Vpair_clusters"],
             ARI_vs_vc_ori_w050=DG.loc["whole molecule (VX4)", "ARI_vs_reference"],
             VXV_E1="%.3f [%.3f, %.3f] (occupancy)" % tuple(E.loc["grid_F0W_occ", ["E1_partial", "E1_lo", "E1_hi"]]),
             VXV_E4="%.3f" % E.loc["grid_F0W_occ", "E4_rel"], state="not tested (VX5 hard stop)"),
        dict(arm="voxel Arm B", scope="all four loops, F1 frame, 7 ch", cut=H6["cut"], control_within_cut=H6["control_within_cut"],
             crystal_err_over_cut=H6["crystal_error_over_cut"], crystal_within_cut=np.nan, clusters=H6["clusters"], clustered=H6["clustered"],
             single_Vpair_clusters=DG.loc["Arm B: F1 loops, 7 ch", "share_single_Vpair_clusters"],
             ARI_vs_vc_ori_w050=DG.loc["Arm B: F1 loops, 7 ch", "ARI_vs_reference"],
             VXV_E1="%.3f [%.3f, %.3f]" % tuple(E.loc["grid_F1L_7ch", ["E1_partial", "E1_lo", "E1_hi"]]),
             VXV_E4="%.3f" % E.loc["grid_F1L_7ch", "E4_rel"], state="tested"),
        dict(arm="voxel Arm C", scope="CDR3 only, F1 frame, 7 ch", cut=C2["cut"], control_within_cut=C2["control_within_cut"],
             crystal_err_over_cut=C2["crystal_error_over_cut"], crystal_within_cut=np.nan, clusters=C2["clusters"], clustered=C2["clustered"],
             single_Vpair_clusters=DG.loc["Arm C: F1 CDR3, 7 ch", "share_single_Vpair_clusters"],
             ARI_vs_vc_ori_w050=DG.loc["Arm C: F1 CDR3, 7 ch", "ARI_vs_reference"],
             VXV_E1="%.3f [%.3f, %.3f] (VXV3)" % tuple(V3.loc["grid_F1C3_7ch", ["E1_partial", "E1_lo", "E1_hi"]]),
             VXV_E4="%.3f" % V3.loc["grid_F1C3_7ch", "E4_rel"], state="tested"),
        dict(arm=_N10NAME, scope="CDR3 arcs, per-chain frames", cut=float(T0.bg_p1), control_within_cut=float(T0.sens_at_bg_p1),
             crystal_err_over_cut=float(V6.F_crystal_over_cut), crystal_within_cut=(int(round(float(V6.C_crystal_within_cut) * _NBP)) if isinstance(_NBP, (int, float)) else None),
             clusters=int(LE[LE >= 0].nunique()), clustered=int((LE >= 0).sum()),
             single_Vpair_clusters=DG.loc["reference method (context)", "share_single_Vpair_clusters"], ARI_vs_vc_ori_w050=1.0,
             VXV_E1="%.3f [%.3f, %.3f] (geometry); full %.3f" % (*E.loc["vec_vcori_geo", ["E1_partial", "E1_lo", "E1_hi"]], E.loc["vec_vcori_full", "E1_partial"]),
             VXV_E4="%.3f (geometry)" % E.loc["vec_vcori_geo", "E4_rel"], state="tested (%s)" % _REF2.get("n10_name", "reference")),
        dict(arm=_RNAME, scope="CDR3 arcs, per-chain frames", cut=float(A5.cut), control_within_cut=float(A5.sens_at_cut),
             crystal_err_over_cut=float(A5.crystal_err_over_cut_chem0), crystal_within_cut=int(A5.crystal_within_cut_chem0),
             clusters=np.nan, clustered=np.nan, single_Vpair_clusters=np.nan, ARI_vs_vc_ori_w050=np.nan, VXV_E1="not in the panel (N = 10 arcs used)",
             VXV_E4="", state="tested (%s)" % _RNAME),
    ]
    T = pd.DataFrame(rows)
    st = {"voxel Arm B": "prop_voxel_armB", "voxel Arm C": "prop_voxel_armC",
          _RNAME: "%s (%s, reproduced; 3 rungs)" % (_refp("label_column"), _RNAME)}
    for lv in ("unstratified", "mouse", "mouse_V", "mouse_V_len"):
        T["excess_" + lv] = [S8.loc[st[a], "delta_" + lv] if a in st and "delta_" + lv in S8.columns and pd.notna(S8.loc[st[a], "delta_" + lv]) else np.nan for a in T.arm]
        T["z_" + lv] = [S8.loc[st[a], "z_" + lv] if a in st and "z_" + lv in S8.columns and pd.notna(S8.loc[st[a], "z_" + lv]) else np.nan for a in T.arm]
    T.loc[T.arm == _N10NAME, ["excess_unstratified", "z_unstratified", "excess_mouse", "z_mouse", "excess_mouse_V", "z_mouse_V"]] = \
        [AC.delta_unstratified, AC.z_unstratified, AC.delta_mouse, AC.z_mouse, AC.delta_mouse_V, AC.z_mouse_V]
    T["survivors_q015"] = [S8.loc[st[a], "survivors_q015"] if a in st else (int((KT.q_bh <= 0.15).sum()) if a == _N10NAME else np.nan) for a in T.arm]
    T.to_csv(os.path.join(OUT, "VXH9_comparison.csv"), index=False)
    C.add("comparison table assembled from outputs", "%d rows" % len(T), len(T) == 5 and T.cut.notna().all(), "5 rows, cuts present")
    C.add("reference state figures in the table equal the reference method's own file", "",
          abs(T.loc[4, "z_mouse_V"] - float(A5.z_mouse_V)) < 1e-9 and abs(T.loc[4, "excess_mouse_V"] - float(A5.excess_mouse_V)) < 1e-9, "equal")

    # A12.5 / A13.5 recorded predictions were scored here against this run's outputs. The scored table is specific to the
    # data set it was run on and is not shipped; re-create it for your own run by listing your pre-registered predictions
    # as rows dict(source=..., prediction=..., result=..., held="yes"/"no"/"partly") and filling `result` from out/.
    P = pd.DataFrame(columns=["source", "prediction", "result", "held"])
    if not len(P):
        C.info("prediction scoring", "no predictions registered in this copy (see the comment above)")
    P.to_csv(os.path.join(OUT, "VXH9_predictions.csv"), index=False)
    if len(P):
        C.info("predictions scored", "%d; held %d, not held %d, other %d" % (len(P), int((P.held == "yes").sum()),
               int(P.held.str.startswith("no").sum()), int(len(P) - (P.held == "yes").sum() - P.held.str.startswith("no").sum())))
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    print(T.to_string(index=False))
    C.write()
