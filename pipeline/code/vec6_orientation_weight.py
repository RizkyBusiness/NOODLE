"""reference step 6 -- choose the side-chain orientation weight for prop_vc_ori WITHOUT using transcriptional states.

prop(w) = sqrt( d_vec^2 + (lam d_chem)^2 + (w * lam_o * d_ori)^2 ),  w in {0, 0.1, ..., 1.0}
(w = 1 is the equal-median default used in the first run; w = 0 is prop_vc). lam, lam_o as in vec2 (fixed).
For each w, using ONLY the pair sets and the crystal benchmark:
  S(w) = share of the one-amino-acid positive-control pairs within the cut (cut = 1st percentile of the
         same 60,000 background pairs, re-derived for every w)
  F(w) = median crystal-model distance / cut over the benchmark pairs (chemistry = 0: same sequence)
  C(w) = share of crystal-model pairs within the cut
Selection rule, fixed before computing: the LARGEST w with
  S(w) >= S(0) - 0.01   (a loss of at most 0.01 of the control pairs relative to no orientation)   and
  F(w) <= F(0) + 0.03   (crystal-model error at most 0.03 of the cut above no orientation).
No state, mouse, cluster or purity information is read by this script.
Writes reference/out/V6_orientation_weight_grid.csv and V6_chosen_weight.json
"""
import numpy as np, pandas as pd, json, os, sys
import paths  # package: config (reference receptor, data-set settings)
D3F, V1F, V3F, LM1F, STB, VO = sys.argv[1:7]
P = np.load(D3F, allow_pickle=True); ids = list(P["clone_id"]); n = len(ids); idx = {c: i for i, c in enumerate(ids)}
V1 = np.load(V1F, allow_pickle=True); assert list(V1["clone_id"]) == ids
CZ = P["chemz"].astype(np.float64); G = V1["vc"].astype(np.float64); U = V1["uc"].astype(np.float64)
th = pd.read_csv(os.path.join(VO, "V2_thresholds_vc_ori.csv")); lam, lamo = float(th.lam.iloc[0]), float(th.lam_ori.iloc[0])
# pair sets exactly as vec2 / lm2 / d5 (same rng stream)
R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
rng = np.random.default_rng(0)
NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
dv = lambda a, b: np.sqrt(((G[a] - G[b]) ** 2).sum(-1) / 200)
dc = lambda a, b: np.sqrt(((CZ[a] - CZ[b]) ** 2).sum(-1).mean(-1))
do = lambda a, b: np.sqrt(((U[a] - U[b]) ** 2).sum(-1) / 20)
bg = (dv(ba, bb), dc(ba, bb), do(ba, bb)); po = (dv(pa, pb), dc(pa, pb), do(pa, pb))
# check: the w = 1 background cut equals the first run's cut
# crystal features (same construction as vec1 / vec4)
L1 = np.load(LM1F, allow_pickle=True); T = L1["lm"][list(L1["clone_id"]).index(paths.reference_receptor())].astype(float)
def kab(X, Y):
    cx, cy = X.mean(0), Y.mean(0); Uu, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy)); d = np.sign(np.linalg.det(Uu @ Vt))
    return Uu @ np.diag([1, 1, d]) @ Vt
PA = list(range(10)) + list(range(20, 30)); PB = list(range(10, 20)) + list(range(30, 40))
def orient(arc, Rm, pos):
    v = arc[[20 + p for p in pos]] - arc[pos]; ln = np.linalg.norm(v, axis=1, keepdims=True)
    return np.where(ln >= 0.5, v / np.maximum(ln, 1e-9), 0.0) @ Rm
def feats(arc, lm):
    Ra = kab(lm[:5], T[:5]); Rb = kab(lm[5:], T[5:])
    v = np.concatenate([((arc[PA][:, None] - lm[:5][None]) @ Ra).reshape(-1), ((arc[PB][:, None] - lm[5:][None]) @ Rb).reshape(-1)])
    return v, np.concatenate([orient(arc, Ra, list(range(10))), orient(arc, Rb, list(range(10, 20)))]).reshape(-1)
for i in [r for r in (0, 7000) if r < n]:   # spot-check rows present in this model set
    v, u = feats(L1["arc"][i].astype(float), L1["lm"][i].astype(float)); assert np.abs(v - G[i]).max() < 1e-3 and np.abs(u - U[i]).max() < 1e-4
C = np.load(V3F, allow_pickle=True); cg, co = [], []
for j in range(len(C["entry"])):
    vc_, uc_ = feats(C["arc_c"][j].astype(float), C["lm_c"][j].astype(float)); vm_, um_ = feats(C["arc_m"][j].astype(float), C["lm_m"][j].astype(float))
    cg.append(np.sqrt(((vc_ - vm_) ** 2).sum() / 200)); co.append(np.sqrt(((uc_ - um_) ** 2).sum() / 20))
cg, co = np.array(cg), np.array(co)
rows = []
for w in np.round(np.arange(0, 1.01, 0.1), 2):
    comb = lambda t: np.sqrt(t[0] ** 2 + (lam * t[1]) ** 2 + (w * lamo * t[2]) ** 2)
    b, p = comb(bg), comb(po); cut = float(np.percentile(b, 1))
    cr = np.sqrt(cg ** 2 + (w * lamo * co) ** 2)
    # share of the background-pair distance (squared) carried by the orientation term: how much it contributes
    share = float(np.median((w * lamo * bg[2]) ** 2 / b ** 2))
    rows.append(dict(w=w, lam_o_effective=round(w * lamo, 4), cut=round(cut, 4), S_sensitivity=round(float((p <= cut).mean()), 4),
                     F_crystal_over_cut=round(float(np.median(cr) / cut), 4), C_crystal_within_cut=round(float((cr <= cut).mean()), 4),
                     orientation_share_of_distance=round(share, 3)))
Gd = pd.DataFrame(rows); Gd.to_csv(os.path.join(VO, "V6_orientation_weight_grid.csv"), index=False)
assert abs(Gd.cut.iloc[-1] - float(th.loc[th.metric == "prop_vc_ori", "bg_p1"].iloc[0])) < 1e-3, "w=1 must reproduce the first run's cut"
assert abs(Gd.cut.iloc[0] - float(pd.read_csv(os.path.join(VO, "V2_thresholds_vc.csv")).query("metric=='prop_vc'").bg_p1.iloc[0])) < 1e-3, "w=0 must reproduce prop_vc"
S0, F0 = Gd.S_sensitivity.iloc[0], Gd.F_crystal_over_cut.iloc[0]
okw = Gd[(Gd.S_sensitivity >= S0 - 0.01) & (Gd.F_crystal_over_cut <= F0 + 0.03)]
wstar = float(okw.w.max())
json.dump(dict(rule="largest w with S(w) >= S(0)-0.01 and F(w) <= F(0)+0.03", S0=S0, F0=F0, chosen_w=wstar,
               lam_o_effective=round(wstar * lamo, 4), lam=lam, lam_o_default=lamo), open(os.path.join(VO, "V6_chosen_weight.json"), "w"), indent=1)
print(Gd.to_string(index=False)); print("chosen w = %.1f (lam_o = %.4f)" % (wstar, wstar * lamo))
