"""arc30 step 6 self-review (27 Sep 2026): re-derive three Stage 4 numbers by routes independent of a30_3.
  1 primary cut: background distances as row norms of one pre-scaled, concatenated feature matrix (float64; lam and
    lam_o recomputed here from the medians), 1st percentile via np.quantile(method='linear')
  2 primary control sensitivity: count of the control pairs within that cut (integer count / number of control pairs)
  3 primary survivor count: per-cluster one-sided binomial p recomputed from the molecule label file (labels + state
    columns, base rates over all state-annotated molecules) and BH q from scipy.stats.false_discovery_control
Writes reference/arc30/checks/SR_independent.csv.
All reads go through paths.src(<rel>), all writes through paths.dst(<rel>).
usage: python a30_6_selfreview.py"""
import numpy as np, pandas as pd, json
from scipy import stats
import paths
AO = "reference/arc30/out/"; W = json.load(open(paths.src(AO + "A2_chosen_weight.json"))); arm = W["primary_arm"]; w = W["chosen_w"]
F = np.load(paths.src(AO + "A1_features_N30.npz"), allow_pickle=True); P = np.load(paths.src(AO + "A1_property_descriptor_N30.npz"), allow_pickle=True)
ids = list(F["clone_id"]); n = len(ids); pos = {c: i for i, c in enumerate(ids)}; N = int(F["N"])
X = np.concatenate([F["vc"].astype(float) / np.sqrt(20 * N)], 1); U = F["uc"].astype(float) / np.sqrt(2 * N)
C = P["chemz"].astype(float).reshape(n, -1) / np.sqrt(2 * N)            # d_chem^2 = sum over props, mean over points
rng = np.random.default_rng(0)
NP = pd.read_csv(paths.src("tables/B3c_near_identical_pairs.csv.gz")); pa = np.array([pos[c] for c in NP.clone_a]); pb = np.array([pos[c] for c in NP.clone_b])
ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); keep = ba != bb; ba, bb = ba[keep], bb[keep]
nrm = lambda M, a, b: np.linalg.norm(M[a] - M[b], axis=1)
lam = np.median(nrm(X, ba, bb)) / np.median(nrm(C, ba, bb)); lamo = w * np.median(nrm(X, ba, bb)) / np.median(nrm(U, ba, bb))
A = np.concatenate([X, lam * C, lamo * U], 1)
cut = float(np.quantile(nrm(A, ba, bb), 0.01, method="linear")); k_in = int((nrm(A, pa, pb) <= cut).sum())
T = pd.read_csv(paths.src(AO + "V2_thresholds_%s.csv" % arm)).set_index("metric").loc["prop_" + arm]
L = pd.read_csv(paths.src(AO + "V2_molecule_labels_%s.csv.gz" % arm)); lab = "prop_" + arm
D = L[L.state.notna()]; base = D.state.value_counts(normalize=True); d = D[D[lab] >= 0]; sz = d[lab].value_counts(); d = d[d[lab].isin(sz[sz >= 2].index)]
p = []
for c, g in d.groupby(lab):
    vc = g.state.value_counts(); top = vc.index[0]; p.append(stats.binom.sf(vc.iloc[0] - 1, len(g), base[top]))
q = stats.false_discovery_control(np.array(p), method="bh"); nsurv = int((q <= 0.15).sum())
K = pd.read_csv(paths.src(AO + "V2_cluster_tests_%s.csv" % arm))
R = pd.DataFrame([dict(quantity="primary cut", independent=round(cut, 4), stage4=float(T.bg_p1), ok=abs(cut - T.bg_p1) < 1e-3),
                  dict(quantity="lam / lam_o", independent="%.4f / %.4f" % (lam, lamo), stage4="%.4f / %.4f" % (T.lam, T.lam_ori), ok=abs(lam - T.lam) < 1e-3 and abs(lamo - T.lam_ori) < 1e-3),
                  dict(quantity="control sensitivity", independent="%d/%d = %.4f" % (k_in, len(pa), k_in / len(pa)), stage4=float(T.sens_at_bg_p1), ok=abs(k_in / len(pa) - T.sens_at_bg_p1) < 1e-4),
                  dict(quantity="clusters tested", independent=len(p), stage4=len(K), ok=len(p) == len(K)),
                  dict(quantity="survivors q <= 0.15", independent=nsurv, stage4=int((K.q_bh <= 0.15).sum()), ok=nsurv == int((K.q_bh <= 0.15).sum()))])
R.to_csv(paths.dst("reference/arc30/checks/SR_independent.csv"), index=False); print(R.to_string(index=False)); print("SELF-REVIEW numbers:", "PASS" if R.ok.all() else "FAIL")
