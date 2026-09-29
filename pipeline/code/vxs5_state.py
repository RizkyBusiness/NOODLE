"""VXS5 (the voxel pipeline, A14.2.5): state tests for the sigma 1.5 A arms D and E, run ONCE each - the VXH8 functions unchanged.

Validation first (as VXH8): the same function on the the reference method primary's labels reproduces its V2_arm_comparison exactly, and
re-running it on Arm B's labels reproduces VXH8's Arm B row (all fields). Then per arm: purity over clusters >= 3; the
four-rung 200-permutation null ladder (default_rng(0) per arm); per-cluster one-sided binomial; BH q <= 0.15.
Overlap with the reference method's survivor sets and with the survivors of the sigma 2.0 arm on the same
atoms (B for D, C for E), both directions, by the reference method's match() rule. Report only (A14.3): an additional look, BH within arm,
no arm called best. Writes out/VXS5_*; checks/VXS5_checks.csv.
usage: python pipeline/code/vxs5_state.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column

from vxh8_state import state_tests, match, A30, REF_LABELS, REF_COL, REF_CMP, REF_TESTS, REF_NAME, _refp

if __name__ == "__main__":
    C = Checks("VXS5", stop_on_fail=True)
    C.info("status", "A14 VXS5 state tests, once per arm (D, E); report only")
    ST = pd.read_csv(os.path.join(STB, "_receptor_states_all.csv.gz")).drop_duplicates("clone_key")[["clone_key", "state", "donor"]]
    AL = pd.read_csv(A30 + REF_LABELS)
    L4 = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"), usecols=["clone_id", "length_class"])
    AL = AL.merge(L4, on="clone_id", how="left")
    ra, Ka, _ = state_tests(AL, REF_COL, fourth=False)
    ref = pd.read_csv(A30 + REF_CMP).iloc[0]
    diff = max(abs(float(ra[k]) - float(ref[k])) for k in ref.index if k != "arm")
    C.add("mirror reproduces %s primary V2_arm_comparison (all fields)" % REF_NAME, "max |diff| %.2e" % diff, diff < 1e-9, "exact")
    S8 = pd.read_csv(os.path.join(OUT, "VXH8_arm_comparison.csv")).set_index("arm")
    MB = pd.read_csv(os.path.join(OUT, "VXH6_molecule_labels.csv.gz")).merge(ST, on="clone_key", how="left")
    rb, Kb, _ = state_tests(MB, "prop_voxel_armB")
    e = max(abs(float(rb[k]) - float(S8.loc["prop_voxel_armB", k])) for k in rb if k != "arm")
    C.add("re-run on Arm B labels reproduces VXH8 Arm B row (all fields)", "max |diff| %.2e" % e, e < 1e-9, "exact")
    _SS = _refp("survivor_sets", {}); A5 = pd.read_csv(A30 + _SS["survivors"]); R8 = pd.read_csv(A30 + _SS["secondary"])
    labA = AL.set_index("clone_id")[REF_COL]; tstA = Ka.set_index("cluster")
    _N10F, _N10C = _refp("n10_labels", ""), _refp("n10_column", "")
    N10 = pd.read_csv(VXP(_N10F), usecols=["clone_id", _N10C]).set_index("clone_id")[_N10C]
    sets = {"%s survivors (%d)" % (REF_NAME, len(A5)): [list(labA.index[labA == c]) for c in A5.cluster],
            "%s robust (%d)" % (REF_NAME, int(A5.robust_n30.sum())): [list(labA.index[labA == c]) for c in A5.cluster[A5.robust_n30]],
            "%s secondary robust (%d)" % (REF_NAME, len(R8)): [list(N10.index[N10 == c]) for c in R8.n10_cluster]}
    rows, overl = [], []
    for aname, f, rname, rf in (("prop_voxel_armD", "VXS2_D_molecule_labels.csv.gz", "prop_voxel_armB", "VXH6_molecule_labels.csv.gz"),
                                ("prop_voxel_armE", "VXS2_E_molecule_labels.csv.gz", "prop_voxel_armC", "VXC2_molecule_labels.csv.gz")):
        M = pd.read_csv(os.path.join(OUT, f)).merge(ST, on="clone_key", how="left")
        C.add("%s: {:,} molecules, order as VX4".format(EXP_NMOL) % aname, len(M), len(M) == EXP_NMOL and
              list(M.clone_id) == list(pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"), usecols=["clone_id"]).clone_id), "identical")
        r, K, nulls = state_tests(M, aname)
        K.to_csv(os.path.join(OUT, "VXS5_cluster_tests_%s.csv" % aname), index=False)
        r["clusters_tested"] = len(K); r["survivors_q015"] = int((K.q_bh <= 0.15).sum()); rows.append(r)
        C.info("%s: clusters >= 3 / molecules / observed purity" % aname, "%d / %d / %.4f" % (r["clusters"], r["molecules"], r["observed"]))
        for lv in ("unstratified", "mouse", "mouse_V", "mouse_V_len"):
            C.info("%s: rung %s - null mean (sd) / excess / z" % (aname, lv), "%.4f (%.4f) / %+.4f / %.1f"
                   % (nulls[lv]["mean"], nulls[lv]["sd"], r["delta_" + lv], r["z_" + lv]))
        C.info("%s: clusters tested (>= 2) / survivors q <= 0.15" % aname, "%d / %d" % (r["clusters_tested"], r["survivors_q015"]))
        lab = M.set_index("clone_id")[aname]; tst = K.set_index("cluster")
        for sname, groups in sets.items():
            res = [match(g, lab, tst) for g in groups]
            overl.append(dict(arm=aname, reference=sname, clusters=len(groups), recovered=sum(x[3] for x in res),
                              median_frac_kept=float(np.median([x[1] for x in res]))))
            C.info("%s: %s recovered (kept >= 75 %% and q <= 0.15)" % (aname, sname), "%d of %d" % (sum(x[3] for x in res), len(groups)))
        sv = K[K.q_bh <= 0.15]
        res = [match(list(lab.index[lab == c]), labA, tstA) for c in sv.cluster]
        overl.append(dict(arm=aname, reference="own survivors recovered in %s primary" % REF_NAME, clusters=len(sv),
                          recovered=sum(x[3] for x in res), median_frac_kept=float(np.median([x[1] for x in res])) if res else np.nan))
        C.info("%s: own survivors recovered in the reference method" % aname, "%d of %d" % (sum(x[3] for x in res), len(sv)))
        # the sigma 2.0 arm on the same atoms, both directions
        RL = pd.read_csv(os.path.join(OUT, rf)).set_index("clone_id")[rname]
        RK = pd.read_csv(os.path.join(OUT, "VXH8_cluster_tests_%s.csv" % rname)).set_index("cluster")
        rs = RK[RK.q_bh <= 0.15]
        res1 = [match(list(RL.index[RL == c]), lab, tst) for c in rs.index]
        res2 = [match(list(lab.index[lab == c]), RL, RK) for c in sv.cluster]
        overl.append(dict(arm=aname, reference="%s survivors recovered here" % rname, clusters=len(rs), recovered=sum(x[3] for x in res1),
                          median_frac_kept=float(np.median([x[1] for x in res1])) if res1 else np.nan))
        overl.append(dict(arm=aname, reference="own survivors recovered in %s" % rname, clusters=len(sv), recovered=sum(x[3] for x in res2),
                          median_frac_kept=float(np.median([x[1] for x in res2])) if res2 else np.nan))
        C.info("%s: %s survivors recovered here / own survivors recovered in %s" % (aname, rname, rname),
               "%d of %d / %d of %d" % (sum(x[3] for x in res1), len(rs), sum(x[3] for x in res2), len(sv)))
        if len(sv):
            C.info("%s: survivors - members / top state / frac / mice / V pairs / q" % aname,
                   " | ".join("%d:%d %s %.2f m%d v%d q%.3f" % (x.cluster, x.members, x.top_state, x.frac, x.mice, x.v_pairs, x.q_bh)
                              for x in sv.itertuples()))
    ref_rows = [dict(S8.loc[a].to_dict(), arm=a) for a in ("prop_voxel_armB", "prop_voxel_armC")]
    pd.DataFrame(rows + ref_rows).to_csv(os.path.join(OUT, "VXS5_arm_comparison.csv"), index=False)
    pd.DataFrame(overl).to_csv(os.path.join(OUT, "VXS5_overlap.csv"), index=False)
    rows_j = [{k: (v.item() if hasattr(v, "item") else v) for k, v in r_.items()} for r_ in rows]
    save_json(dict(seed=0, n_perm=200, rungs=["unstratified", "mouse", "mouse_V", "mouse_V_len"], bh_q=0.15, arms=rows_j,
                   status="report only (A14.3): an additional look; BH within arm"), os.path.join(OUT, "VXS5_summary.json"))
    C.write()
