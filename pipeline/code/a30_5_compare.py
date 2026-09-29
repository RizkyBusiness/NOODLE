"""arc30 step 5 (27 Sep 2026): compare the single N = 30 state-test run (5 arms) with the N = 10 method (vec4/vec7 style).
Reads only arc30 outputs (A1 crystal features + chemistry sd for the crystal error, A2 grid, the Stage 4 tables and
labels) and the existing N = 10 files (reference/out/V2_*, D3 chem_sd). Nothing is re-clustered or re-tested.
Writes reference/arc30/out/
  A5_summary.csv                per N = 30 arm + the N = 10 primary: lam, lam_o, cut, control sensitivity, crystal
                                error / cut (chemistry = 0 and with chemistry), crystal pairs within cut, clusters >= 3,
                                excess and z at the three null levels, survivors (q <= 0.15)
  A5_ARI_matrix.csv             singleton-aware ARI (vec7) between the 5 N = 30 arms and the 5 N = 10 vector arms
                                (vc_ori_w050 = the N = 10 primary, vc_ori, vc, vg, vs)
  A5_survivors.csv              each N = 30 primary survivor, its main cluster / fraction kept / q in each N = 30 variant
                                and in the N = 10 primary, and robust_n30 (rule 3.4: kept >= 75 % and q <= 0.15 in all
                                four N = 30 variants)
  A5_n10_robust8_recovery.csv   the robust N = 10 clusters: best-matching N = 30 primary cluster, fraction kept, q
All reads go through paths.src(<rel>), all writes through paths.dst(<rel>).
usage: python a30_5_compare.py
"""
import numpy as np, pandas as pd, json
from math import comb
import paths
AO, VO = "reference/arc30/out/", "reference/out/"
W = json.load(open(paths.src(AO + "A2_chosen_weight.json"))); PRIM = W["primary_arm"]
A30 = [PRIM, "a30_vc_ori", "a30_vc", "a30_vg", "a30_vs"]; VAR = A30[1:]
A10 = ["vc_ori_w050", "vc_ori", "vc", "vg", "vs"]; P10 = "vc_ori_w050"
# the robust N = 10 clusters of the development record, each named by one member molecule (voxel config "expected"
# reference_robust_members_N10) and mapped to its cluster in this run's N = 10 labels: cluster ids change when exact ties
# in the linkage reorder. Context only: unset -> the recovery table is written without rows
ROBM = paths.expected("reference_robust_members_N10")


def ari(a, b):
    ct = pd.crosstab(a, b).values; s = sum(comb(int(x), 2) for x in ct.ravel()); sa = sum(comb(int(x), 2) for x in ct.sum(1)); sb = sum(comb(int(x), 2) for x in ct.sum(0))
    n = comb(int(ct.sum()), 2); e = sa * sb / n; return (s - e) / (0.5 * (sa + sb) - e)


def sing(x):
    x = np.asarray(x).copy(); m = x < 0; x[m] = -(np.arange(m.sum()) + 1); return x


def arm_files(a):
    d = VO if a in A10 and not a.startswith("a30_") else AO
    return dict(thr=pd.read_csv(paths.src(d + "V2_thresholds_%s.csv" % a)).set_index("metric"), cmp=pd.read_csv(paths.src(d + "V2_arm_comparison_%s.csv" % a)).iloc[0],
                tst=pd.read_csv(paths.src(d + "V2_cluster_tests_%s.csv" % a)), lab=pd.read_csv(paths.src(d + "V2_molecule_labels_%s.csv.gz" % a)))


def key(a, n):
    return ("N%d:" % n) + a


F = {key(a, 30): arm_files(a) for a in A30}; F.update({key(a, 10): arm_files(a) for a in A10})
ref_ids = F[key(PRIM, 30)]["lab"].clone_id.values
for k_, f in F.items(): assert (f["lab"].clone_id.values == ref_ids).all(), k_
LAB = {k_: f["lab"]["prop_" + k_.split(":")[1]].values for k_, f in F.items()}
TST = {k_: f["tst"].set_index("cluster") for k_, f in F.items()}

# ---------------------------------------------------------------- crystal error per arm
CR = {30: np.load(paths.src(AO + "A1_crystal_features_N30.npz"), allow_pickle=True), 10: np.load(paths.src(AO + "A1_crystal_features_N10.npz"), allow_pickle=True)}
SD = {30: np.load(paths.src(AO + "A1_property_descriptor_N30.npz"))["chem_sd"],
      10: np.load(paths.src("descriptors/out/D3_property_descriptor.npz"), allow_pickle=True)["chem_sd"]}   # the N = 10 primary used the D3 chemz
NV = {"vg": 40, "vc": 20, "vs": 16}


def crystal(a, N, lam, lamo):
    C = CR[N]; g = (a[4:] if a.startswith("a30_") else a).split("_")[0]
    dg = np.sqrt(((C[g + "_c"].astype(float) - C[g + "_m"].astype(float)) ** 2).sum(1) / (NV[g] * N))
    do = np.sqrt(((C["u%s_c" % g[1]].astype(float) - C["u%s_m" % g[1]].astype(float)) ** 2).sum(1) / (2 * N)) if "_ori" in a else np.zeros(len(dg))
    dc = np.sqrt((((C["PR_c"].astype(float) - C["PR_m"].astype(float)) / SD[N]) ** 2).sum(-1).mean(-1))
    return np.sqrt(dg ** 2 + (lamo * do) ** 2), np.sqrt(dg ** 2 + (lam * dc) ** 2 + (lamo * do) ** 2)


rows = []
for N, arms in ((30, A30), (10, [P10])):
    for a in arms:
        f = F[key(a, N)]; t = f["thr"].loc["prop_" + a]; c = f["cmp"]; lam, lamo, cut = float(t.lam), float(t.lam_ori), float(t.bg_p1)
        c0, c1 = crystal(a, N, lam, lamo)
        r = dict(arm=("" if N == 30 else "N10 primary: ") + "prop_" + a, N=N, lam=lam, lam_o=lamo, cut=cut, sens_at_cut=float(t.sens_at_bg_p1),
                 pos_median_over_cut=round(float(t.pos_median) / cut, 3), crystal_err_over_cut_chem0=round(float(np.median(c0) / cut), 4),
                 crystal_err_over_cut_chem=round(float(np.median(c1) / cut), 4), crystal_within_cut_chem0=int((c0 <= cut).sum()),
                 crystal_within_cut_chem=int((c1 <= cut).sum()), crystal_pairs=len(c0), clusters_ge3=int(c.clusters), molecules_in_ge3=int(c.molecules),
                 observed_purity=float(c.observed))
        for lv in ("unstratified", "mouse", "mouse_V"): r["excess_" + lv] = float(c["delta_" + lv]); r["z_" + lv] = float(c["z_" + lv])
        r["clusters_tested"] = len(f["tst"]); r["survivors_q015"] = int((f["tst"].q_bh <= 0.15).sum())
        rows.append(r)
S = pd.DataFrame(rows)
# consistency with Stage 2 (primary, chemistry = 0) and with the V6 grid (N = 10 primary)
G2 = pd.read_csv(paths.src(AO + "A2_weight_grid.csv")).set_index("w").loc[W["chosen_w"]]
assert abs(S.iloc[0].crystal_err_over_cut_chem0 - G2.F_crystal_over_cut) < 1e-3 and abs(S.iloc[0].cut - G2.cut) < 1e-3, "primary crystal error differs from A2"
V6 = pd.read_csv(paths.src(VO + "V6_orientation_weight_grid.csv")).set_index("w").loc[paths.expected_or_stop("reference_orientation_weight_N10")]
assert abs(S.iloc[-1].crystal_err_over_cut_chem0 - V6.F_crystal_over_cut) < 1e-3, "N = 10 crystal error differs from V6"

# ---------------------------------------------------------------- ARI matrix
names = [key(a, 30) for a in A30] + [key(a, 10) for a in A10]
M = pd.DataFrame([[round(ari(sing(LAB[x]), sing(LAB[y])), 3) for y in names] for x in names], index=names, columns=names)


# ---------------------------------------------------------------- survivors and robustness (rule 3.4)
def match(mem, k_):
    """main cluster of the members `mem` in arm k_, fraction of mem in it, its q, kept (>= 75 %) and significant."""
    b = pd.Series(LAB[k_][mem]); b = b[b >= 0]
    if not len(b): return -1, 0.0, np.nan, False
    main = int(b.value_counts().index[0]); fr = float(b.value_counts().iloc[0] / len(mem))
    q = float(TST[k_].loc[main, "q_bh"]) if main in TST[k_].index else np.nan
    return main, fr, q, bool(fr >= 0.75 and q == q and q <= 0.15)


x = LAB[key(PRIM, 30)]; KP = F[key(PRIM, 30)]["tst"]; pr = []
for _, r in KP[KP.q_bh <= 0.15].sort_values("q_bh").iterrows():
    mem = np.where(x == r.cluster)[0]
    row = dict(cluster=int(r.cluster), members_tested=int(r.members), members_all=len(mem), top_state=r.top_state, frac=r.frac, q=round(float(r.q_bh), 4),
               mice=int(r.mice), v_pairs=int(r.v_pairs))
    ok = []
    for k_ in [key(a, 30) for a in VAR] + [key(P10, 10)]:
        main, fr, q, ks = match(mem, k_); tag = k_.replace(":", "_")
        row.update({tag + "_main": main, tag + "_frac_kept": round(fr, 3), tag + "_q": round(q, 4) if q == q else np.nan, tag + "_kept_and_significant": ks})
        if k_.startswith("N30"): ok.append(ks)
    row["n30_variants_kept_and_significant"] = int(sum(ok)); row["robust_n30"] = bool(all(ok)); pr.append(row)
PR = pd.DataFrame(pr)
base = LAB[key(P10, 10)]; T10 = TST[key(P10, 10)]; rec = []
if paths.missing(ROBM):
    print("NOTE: %s -- no robust N = 10 clusters to follow; A5_n10_robust8_recovery.csv is written without rows" % ROBM); ROB8 = []
else:
    pos = {c: i for i, c in enumerate(ref_ids)}; bad = [m for m in ROBM if m not in pos or base[pos[m]] < 0]
    if bad:
        raise SystemExit("reference_robust_members_N10: not in a cluster of the N = 10 labels: %s" % ", ".join(map(str, bad)))
    ROB8 = [int(base[pos[m]]) for m in ROBM]
for c in ROB8:
    mem = np.where(base == c)[0]; main, fr, q, ks = match(mem, key(PRIM, 30))
    rec.append(dict(n10_cluster=c, n10_members_all=len(mem), n10_top_state=T10.loc[c, "top_state"], n10_q=round(float(T10.loc[c, "q_bh"]), 4),
                    n30_best_cluster=main, frac_kept=round(fr, 3), n30_q=round(q, 4) if q == q else np.nan,
                    n30_top_state=TST[key(PRIM, 30)].loc[main, "top_state"] if main in TST[key(PRIM, 30)].index else "",
                    n30_best_is_robust_n30=bool(main in set(PR[PR.robust_n30].cluster)), recovered=ks))
RC = pd.DataFrame(rec, columns=["n10_cluster", "n10_members_all", "n10_top_state", "n10_q", "n30_best_cluster", "frac_kept", "n30_q",
                               "n30_top_state", "n30_best_is_robust_n30", "recovered"])
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 60)
for nm, T_ in (("A5_summary", S), ("A5_ARI_matrix", M.reset_index().rename(columns={"index": "arm"})), ("A5_survivors", PR), ("A5_n10_robust8_recovery", RC)):
    T_.to_csv(paths.dst(AO + nm + ".csv"), index=False); print("\n== %s ==\n%s" % (nm, T_.to_string(index=False)))
print("\nrobust clusters at N = 30: %d of %d survivors | N = 10 robust-8 recovered: %d of %d" % (int(PR.robust_n30.sum()), len(PR), int(RC.recovered.sum()), len(ROB8)))
