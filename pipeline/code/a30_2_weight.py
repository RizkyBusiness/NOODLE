"""arc30 step 2 (27 Sep 2026): choose the side-chain orientation weight w at N arc points WITHOUT state data.
vec6 logic; crystal features from a30_1; lam, lam_o and the cut recomputed for every w (runbook section 2):
  prop(w) = sqrt( d_vec^2 + (lam d_chem)^2 + (w lam_o d_ori)^2 ),  lam = med_bg(d_vec)/med_bg(d_chem),
  lam_o = med_bg(d_vec)/med_bg(d_ori), d_vec = ||dG||/sqrt(20N), d_ori = ||dU||/sqrt(2N), cut = bg 1st percentile
  S(w)      share of the control pairs within the cut
  F(w)      median crystal-model distance / cut over the benchmark pairs, chemistry = 0 (as vec6)
  F_chem(w) the same with the chemistry term computed
  C(w)      share of crystal-model pairs within the cut
Rule 3.2 (frozen): the largest w with S(w) >= S(0) - 0.01 and F(w) <= F(0) + 0.03 (on the grid values rounded to
4 decimals, as vec6). Read: A1 features, chemz, crystal features, B3c control pairs, _slim_receptors (not used
beyond the vec6 stream). No state, mouse, cluster or purity information is read.
N = 10 -> gate G2 vs V6_orientation_weight_grid.csv (checks/); N = 30 -> A2_weight_grid.csv, A2_chosen_weight.json,
and gate G3 vs ARCN3_mode_metrics.csv (checks/).
All reads go through paths.src(<rel>), all writes through paths.dst(<rel>).
usage: python a30_2_weight.py N
"""
import numpy as np, pandas as pd, json, os, sys
import paths
# data-set expected values (voxel config "expected", paths.expected); an unset one fails the assert / gate that uses it
EXP_NCP, EXP_NBG, EXP_NBP = paths.expected("n_control_pairs"), paths.expected("n_background_pairs"), paths.expected("n_benchmark_pairs")
EXP_W10 = paths.expected("reference_orientation_weight_N10")
N = int(sys.argv[1]); AO = "reference/arc30/out/"; CHK = "reference/arc30/checks/"; VO = "reference/out/"; STB = "tables"
V1 = np.load(paths.src(AO + "A1_features_N%d.npz" % N), allow_pickle=True); ids = list(V1["clone_id"]); n = len(ids); idx = {c: i for i, c in enumerate(ids)}
P = np.load(paths.src(AO + "A1_property_descriptor_N%d.npz" % N), allow_pickle=True); assert list(P["clone_id"]) == ids
NV = dict(zip([str(x) for x in V1["n_vec_keys"]], [int(x) for x in V1["n_vec"]])); k, ku = NV["vc"], 2 * N
CZ = P["chemz"].astype(np.float64); G = V1["vc"].astype(np.float64); U = V1["uc"].astype(np.float64); sd = P["chem_sd"]
# pair sets exactly as vec2 / vec6 (same rng stream)
R = pd.read_csv(paths.src(os.path.join(STB, "_slim_receptors.csv.gz"))).set_index("clone_id").reindex(ids)
rng = np.random.default_rng(0)
NP = pd.read_csv(paths.src(os.path.join(STB, "B3c_near_identical_pairs.csv.gz")))
pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
assert len(pa) == EXP_NCP and len(ba) == EXP_NBG, (len(pa), len(ba), str(EXP_NCP), str(EXP_NBG))
dv = lambda a, b: np.sqrt(((G[a] - G[b]) ** 2).sum(-1) / k)
dc = lambda a, b: np.sqrt(((CZ[a] - CZ[b]) ** 2).sum(-1).mean(-1))
do = lambda a, b: np.sqrt(((U[a] - U[b]) ** 2).sum(-1) / ku)
bg = (dv(ba, bb), dc(ba, bb), do(ba, bb)); po = (dv(pa, pb), dc(pa, pb), do(pa, pb))
lam = float(np.median(bg[0]) / np.median(bg[1])); lamo = float(np.median(bg[0]) / np.median(bg[2]))
C = np.load(paths.src(AO + "A1_crystal_features_N%d.npz" % N), allow_pickle=True); assert len(C["entry"]) == EXP_NBP, (len(C["entry"]), str(EXP_NBP))
cg = np.sqrt(((C["vc_c"].astype(float) - C["vc_m"].astype(float)) ** 2).sum(1) / k)
co = np.sqrt(((C["uc_c"].astype(float) - C["uc_m"].astype(float)) ** 2).sum(1) / ku)
cc = np.sqrt((((C["PR_c"].astype(float) - C["PR_m"].astype(float)) / sd) ** 2).sum(-1).mean(-1))
rows = []
for w in np.round(np.arange(0, 1.01, 0.1), 2):
    comb = lambda t: np.sqrt(t[0] ** 2 + (lam * t[1]) ** 2 + (w * lamo * t[2]) ** 2)
    b, p = comb(bg), comb(po); cut = float(np.percentile(b, 1))
    cr = np.sqrt(cg ** 2 + (w * lamo * co) ** 2); crc = np.sqrt(cg ** 2 + (lam * cc) ** 2 + (w * lamo * co) ** 2)
    share = float(np.median((w * lamo * bg[2]) ** 2 / b ** 2))
    rows.append(dict(w=w, lam_o_effective=round(w * lamo, 4), cut=round(cut, 4), S_sensitivity=round(float((p <= cut).mean()), 4),
                     F_crystal_over_cut=round(float(np.median(cr) / cut), 4), C_crystal_within_cut=round(float((cr <= cut).mean()), 4),
                     orientation_share_of_distance=round(share, 3), F_chem=round(float(np.median(crc) / cut), 4),
                     C_chem_within_cut=round(float((crc <= cut).mean()), 4), lam=round(lam, 4), cut_full=cut))
Gd = pd.DataFrame(rows)
S0, F0 = Gd.S_sensitivity.iloc[0], Gd.F_crystal_over_cut.iloc[0]
okw = Gd[(Gd.S_sensitivity >= S0 - 0.01) & (Gd.F_crystal_over_cut <= F0 + 0.03)]
wstar = float(okw.w.max()) if len(okw) else np.nan
Gd["admissible"] = Gd.w.isin(okw.w)
pd.set_option("display.width", 250); print(Gd.to_string(index=False)); print("N=%d lam %.4f lam_o(w=1) %.4f -> chosen w = %s" % (N, lam, lamo, wstar))
if N == 10:
    Gd.to_csv(paths.dst(CHK + "G2_weight_grid_N10.csv"), index=False)
    V6 = pd.read_csv(paths.src(VO + "V6_orientation_weight_grid.csv"))
    rr = []
    for (_, a), (_, b) in zip(Gd.iterrows(), V6.iterrows()):
        rr.append(dict(w=a.w, cut=a.cut, cut_V6=b.cut, S=a.S_sensitivity, S_V6=b.S_sensitivity, F=a.F_crystal_over_cut, F_V6=b.F_crystal_over_cut,
                       ok=bool(a.w == b.w and abs(a.cut_full - b.cut) <= 1e-3 and a.S_sensitivity == b.S_sensitivity and abs(a.F_crystal_over_cut - b.F_crystal_over_cut) <= 1e-3)))
    G2 = pd.DataFrame(rr); G2["chosen_w"] = wstar; G2.to_csv(paths.dst(CHK + "G2_vs_V6.csv"), index=False)
    print(G2.to_string(index=False)); print("GATE G2:", "PASS" if (G2.ok.all() and len(G2) == 11 and wstar == EXP_W10) else "FAIL", *([str(EXP_W10)] if paths.missing(EXP_W10) else []))
else:
    Gd.to_csv(paths.dst(AO + "A2_weight_grid.csv"), index=False)
    if N == 30:
        M3 = pd.read_csv(paths.src(VO + "ARCN3_mode_metrics.csv")).query("N == 30 and mode == 'lin'").set_index("w")
        rr = []
        for w in (0.0, 0.5, 1.0):
            a = Gd.set_index("w").loc[w]; b = M3.loc[w]
            rr.append(dict(w=w, cut=a.cut_full, cut_ARCN3=b.cut, S=a.S_sensitivity, S_ARCN3=b.S, F=a.F_crystal_over_cut, F_ARCN3=b.F,
                           ok=bool(abs(a.cut_full - b.cut) <= 1e-3 and abs(a.S_sensitivity - b.S) <= 1e-3 and abs(a.F_crystal_over_cut - b.F) <= 1e-3)))
        G3 = pd.DataFrame(rr); G3.to_csv(paths.dst(CHK + "G3_vs_ARCN3.csv"), index=False)
        print(G3.to_string(index=False)); g3 = bool(G3.ok.all()); print("GATE G3:", "PASS" if g3 else "FAIL")
        if not g3 or wstar != wstar:
            print("STOP: G3 failed or no admissible w -- A2_chosen_weight.json NOT written"); sys.exit(2)
        arm = "a30_vc" if wstar == 0 else "a30_vc_ori_w%03d" % round(100 * wstar)
        json.dump(dict(rule="largest w with S(w) >= S(0)-0.01 and F(w) <= F(0)+0.03", N=N, S0=float(S0), F0=float(F0), chosen_w=wstar,
                       admissible_w=[float(x) for x in okw.w], lam=round(lam, 4), lam_o_default=round(lamo, 4), lam_o_effective=round(wstar * lamo, 4),
                       primary_arm=arm, state_data_read=False), open(paths.dst(AO + "A2_chosen_weight.json"), "w"), indent=1)
        print("wrote A2_chosen_weight.json: primary arm", arm)
