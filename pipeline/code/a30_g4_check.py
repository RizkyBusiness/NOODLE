"""arc30 gate G4 check: compare the a30_3 N = 10 run (checks/G4/) with the existing vec2 outputs of vc_ori_w050.
Writes reference/arc30/checks/G4_vs_V2.csv.
All reads go through paths.src(<rel>), all writes through paths.dst(<rel>).
usage: python a30_g4_check.py"""
import numpy as np, pandas as pd
import paths
from math import comb
# data-set expected values of gate G4 (voxel config "expected", paths.expected); an unset one fails its row, named in the row
E_CUT, E_LAM, E_LAMO = (paths.expected(k) for k in ("reference_cut_vc_ori_w050", "reference_lam_vc_ori_w050", "reference_lam_ori_vc_ori_w050"))
E_OBS, E_DMV, E_ZMV = (paths.expected(k) for k in ("reference_observed_purity_vc_ori_w050", "reference_excess_mouse_V_vc_ori_w050",
                                                   "reference_z_mouse_V_vc_ori_w050"))
E_NSURV = paths.expected("reference_n_survivors_vc_ori_w050")
VO, G = "reference/out/", "reference/arc30/checks/G4/"
def ari(a, b):
    ct = pd.crosstab(a, b).values; s = sum(comb(int(x), 2) for x in ct.ravel()); sa = sum(comb(int(x), 2) for x in ct.sum(1)); sb = sum(comb(int(x), 2) for x in ct.sum(0))
    n = comb(int(ct.sum()), 2); e = sa * sb / n; return (s - e) / (0.5 * (sa + sb) - e)
def sing(x):
    x = np.asarray(x).copy(); m = x < 0; x[m] = -(np.arange(m.sum()) + 1); return x
a = "vc_ori_w050"; c = "prop_" + a; rows = []
T0, T1 = pd.read_csv(paths.src(VO + "V2_thresholds_%s.csv" % a)), pd.read_csv(paths.src(G + "V2_thresholds_%s.csv" % a))
rows.append(dict(check="thresholds table identical", value=bool(T0.equals(T1)), ok=bool(T0.equals(T1))))
r = T1.set_index("metric").loc[c]
ok = abs(r.bg_p1 - E_CUT) < 1e-4 and abs(r.lam - E_LAM) < 1e-4 and abs(r.lam_ori - E_LAMO) < 1e-4
rows.append(dict(check="cut / lam / lam_o = %s / %s / %s" % (E_CUT, E_LAM, E_LAMO), value="%.4f / %.4f / %.4f" % (r.bg_p1, r.lam, r.lam_ori), ok=bool(ok)))
L0, L1 = pd.read_csv(paths.src(VO + "V2_molecule_labels_%s.csv.gz" % a)), pd.read_csv(paths.src(G + "V2_molecule_labels_%s.csv.gz" % a))
same_ids = bool((L0.clone_id.values == L1.clone_id.values).all())
A = ari(sing(L0[c].values), sing(L1[c].values))
sz0 = sorted(pd.Series(L0[c][L0[c] >= 0]).value_counts().tolist()); sz1 = sorted(pd.Series(L1[c][L1[c] >= 0]).value_counts().tolist())
rows.append(dict(check="molecule order identical", value=same_ids, ok=same_ids))
rows.append(dict(check="partition ARI (singleton-aware)", value=round(A, 6), ok=bool(round(A, 3) == 1.0)))
rows.append(dict(check="cluster size multiset identical", value="%d clusters" % len(sz1), ok=bool(sz0 == sz1)))
rows.append(dict(check="labels identical element-wise", value=bool((L0[c].values == L1[c].values).all()), ok=True))
C0, C1 = pd.read_csv(paths.src(VO + "V2_arm_comparison_%s.csv" % a)), pd.read_csv(paths.src(G + "V2_arm_comparison_%s.csv" % a))
rc = C1.iloc[0]
rows.append(dict(check="arm comparison identical", value=bool(C0.equals(C1)), ok=bool(C0.equals(C1))))
rows.append(dict(check="observed / delta mouse+V / z = %s / %s / %s" % (E_OBS, E_DMV, E_ZMV), value="%.4f / %.4f / %.1f" % (rc.observed, rc.delta_mouse_V, rc.z_mouse_V),
                 ok=bool(rc.observed == E_OBS and rc.delta_mouse_V == E_DMV and rc.z_mouse_V == E_ZMV)))
K0, K1 = pd.read_csv(paths.src(VO + "V2_cluster_tests_%s.csv" % a)), pd.read_csv(paths.src(G + "V2_cluster_tests_%s.csv" % a))
n1 = int((K1.q_bh <= 0.15).sum())
rows.append(dict(check="clusters at q <= 0.15 = %s" % E_NSURV, value=n1, ok=n1 == E_NSURV))
# clusters are matched by their member sets, not by id: ids depend on the order of exact ties in the linkage, which
# last-bit float differences between machines can change without changing any cluster (approved 29 Sep 2026)
def by_members(L, K):
    m = L[L[c] >= 0].groupby(c).clone_id.apply(frozenset)
    K = K.assign(member_set=K.cluster.map(m))
    return K.set_index("member_set").drop(columns="cluster")
kk = ["members", "top_state", "top_n"]
M0, M1 = by_members(L0, K0), by_members(L1, K1)
same_set = bool(M0.index.notna().all() and M1.index.notna().all() and set(M0.index) == set(M1.index) and len(M0) == len(M1))
if same_set:
    M1 = M1.loc[M0.index]
    same_k = bool(M0[kk].reset_index(drop=True).equals(M1[kk].reset_index(drop=True)))
    qd = float(np.abs(M0.q_bh.values - M1.q_bh.values).max())
else:
    same_k, qd = False, float("nan")
rows.append(dict(check="cluster tests identical, clusters matched by member set (members, top state, max |dq|)",
                 value="%s, max|dq| %.1e" % (same_set and same_k, qd), ok=bool(same_set and same_k and qd < 1e-9)))
R = pd.DataFrame(rows); R.to_csv(paths.dst("reference/arc30/checks/G4_vs_V2.csv"), index=False)
pd.set_option("display.width", 200); print(R.to_string(index=False)); print("GATE G4:", "PASS" if R.ok.all() else "FAIL")
