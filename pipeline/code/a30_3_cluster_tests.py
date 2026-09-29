"""arc30 step 3 (27 Sep 2026): a copy of vec2_cluster_tests.py. The ONLY changes (runbook Stage 3):
  (a) k read from the features' n_vec (vg 40N, vc 20N, vs 16N) and the orientation count 2N from the features' N,
      instead of the hard-coded NVEC = {vg: 400, vc: 200, vs: 160} and 20
  (b) input/output paths into reference/arc30/ (given on the command line, as before)
  (c) arm names may carry the prefix a30_ (stripped when choosing the geometry block; kept in output names)
Path arguments: used as given when absolute (the pipeline runner resolves its src:/dst:/dstdir: tokens); a relative
argument is a project-relative path, resolved per file with paths.src(<rel>) for reads and paths.dst_dir(<rel>) for OUT_DIR.
usage: python a30_3_cluster_tests.py CHEMZ.npz A1_features_N.npz TABLES_DIR OUT_DIR ARM [ARM ...]
------------------------------------------------------------------------------------------------------------------
PRIMARY METHOD (from 25 Sep 2026): arm vc_ori_w050 -- per-chain frames + orientation at w = 0.5 + chemistry.
reference-method step 2 (vec2) -- thresholds, clustering and state tests for the vector arms.
Identical machinery to lm2_cluster_tests.py (validated to reproduce v2 prop_both exactly);
only the geometry block changes, plus an optional third (orientation) term:
  prop = sqrt( d_vec^2 + (lam_c d_chem)^2 [+ (lam_o d_ori)^2] ),  each lam = median_bg(d_vec) / median_bg(term)
Arms: vg, vc, vs (geometry + chemistry); vc_ori (vc + side-chain orientation, equal-median weight);
vc_ori_wNNN (orientation weight NNN/100 x default; the weight is chosen by vec6 without state data).
usage: python vec2_cluster_tests.py D3.npz V1_features.npz TABLES_DIR OUT_DIR ARM [ARM ...]
"""
import numpy as np, pandas as pd, os, sys, time
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist
from scipy import stats
import paths

D3F, LM1F, STB, OUT = sys.argv[1:5]
ARMS = sys.argv[5:]
NPERM = 200
_src = lambda p: p if os.path.isabs(p) else paths.src(p)                                # port: path handling only
_in = lambda d, f: os.path.join(d, f) if os.path.isabs(d) else paths.src(os.path.join(d, f))
if os.path.isabs(OUT):
    os.makedirs(OUT, exist_ok=True)
else:
    OUT = paths.dst_dir(OUT)

P = np.load(_src(D3F), allow_pickle=True)
ids = list(P["clone_id"])
n = len(ids)
idx = {c: i for i, c in enumerate(ids)}
CZ = P["chemz"].astype(np.float32)
V1 = np.load(_src(LM1F), allow_pickle=True)
assert list(V1["clone_id"]) == ids, "V1 order differs from D3"
NVEC = dict(zip([str(x) for x in V1["n_vec_keys"]], [int(x) for x in V1["n_vec"]]))   # (a)
NORI = 2 * int(V1["N"])                                                                   # (a)


def geometry(arm):
    """flat G, k with d_vec = ||G_a - G_b|| / sqrt(k); orientation block U, ku (or None)."""
    g = arm[4:].split("_")[0] if arm.startswith("a30_") else arm.split("_")[0]           # (c)
    G = V1[g].astype(np.float64)
    U = V1["u" + g[1]].astype(np.float64) if "_ori" in arm else None
    return G, NVEC[g], U, NORI                                                            # (a)


# ---------------------------------------------------------------- pair sets (as d5)
R = pd.read_csv(_in(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
lclass = (R.cdr3_A.str.len().astype(str) + "-" + R.cdr3_B.str.len().astype(str)).values
rng = np.random.default_rng(0)
NP = pd.read_csv(_in(STB, "B3c_near_identical_pairs.csv.gz"))
pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
NB = 60000
ba, bb = rng.integers(0, n, NB), rng.integers(0, n, NB)
ok = ba != bb; ba, bb = ba[ok], bb[ok]
lc_of = pd.Series(lclass)
groups = {k: v.values for k, v in lc_of.groupby(lc_of).groups.items()}
big = [k for k, v in groups.items() if len(v) >= 2]
w = np.array([len(groups[k]) for k in big], float); w /= w.sum()
wa, wb = [], []
for k in rng.choice(len(big), 40000, p=w):
    v = groups[big[k]]
    i, j = rng.integers(0, len(v), 2)
    if i != j:
        wa.append(v[i]); wb.append(v[j])
wa, wb = np.array(wa), np.array(wb)
CHK_A, CHK_B = rng.integers(0, n, 200), rng.integers(0, n, 200)

# ---------------------------------------------------------------- molecules (as d6)
ID = pd.read_csv(_in(STB, "C1w_receptor_identity.csv.gz"))
SL = pd.read_csv(_in(STB, "_slim_receptors.csv.gz"))
ST = pd.read_csv(_in(STB, "_receptor_states_all.csv.gz")).drop_duplicates("clone_key")
M0 = (pd.DataFrame({"clone_id": ids})
      .merge(SL[["clone_id", "clone_key"]], on="clone_id", how="left")
      .merge(ID.drop_duplicates("clone_key")[["clone_key", "prot_key", "v_A_prot", "v_B_prot"]],
             on="clone_key", how="left")
      .merge(ST[["clone_key", "state", "donor"]], on="clone_key", how="left"))
keep = M0.prot_key.notna().values & ~M0.prot_key.duplicated().values
rep = np.where(keep)[0]
M0 = M0.iloc[rep].reset_index(drop=True)
nm = len(M0)
print("pairs: positive %d | background %d | within %d | molecules %d"
      % (len(pa), len(ba), len(wa), nm), flush=True)


def summ(metric, p, b, w_):
    return dict(metric=metric, pos_median=round(float(np.median(p)), 4),
                pos_p90=round(float(np.percentile(p, 90)), 4),
                bg_p1=round(float(np.percentile(b, 1)), 4),
                bg_median=round(float(np.median(b)), 4),
                within_p1=round(float(np.percentile(w_, 1)), 4),
                within_median=round(float(np.median(w_)), 4),
                sens_at_bg_p1=round(float((p <= np.percentile(b, 1)).mean()), 4),
                sens_at_within_p1=round(float((p <= np.percentile(w_, 1)).mean()), 4))


def cluster(A, cut):
    D = pdist(A.astype(np.float32), metric="euclidean").astype(np.float32)
    f = fcluster(linkage(D, method="complete"), t=cut, criterion="distance")
    del D
    cnt = pd.Series(f).value_counts()
    out = np.full(len(A), -1)
    for j, c in enumerate(cnt[cnt >= 2].index):
        out[f == c] = j
    return out


def purity(codes, states, ncl):
    tot = np.bincount(codes, minlength=ncl).astype(float)
    best = np.zeros(ncl)
    for s in range(states.max() + 1):
        best = np.maximum(best, np.bincount(codes[states == s], minlength=ncl))
    return float((best / tot).mean())


def dchem(a, b):
    return np.sqrt(((CZ[a] - CZ[b]) ** 2).sum(-1).mean(-1))


for arm in ARMS:
    t0 = time.time()
    G, k, U, ku = geometry(arm)
    gname = arm + "_geom"
    aname = "prop_" + arm

    def dg(a, b):
        return np.sqrt(((G[a] - G[b]) ** 2).sum(-1) / k)

    def do(a, b):
        return np.sqrt(((U[a] - U[b]) ** 2).sum(-1) / ku) if U is not None else np.zeros(len(a))

    rows = [summ(gname, dg(pa, pb), dg(ba, bb), dg(wa, wb)),
            summ("prop_chem", dchem(pa, pb), dchem(ba, bb), dchem(wa, wb))]
    if U is not None:
        rows.append(summ("orientation", do(pa, pb), do(ba, bb), do(wa, wb)))
    lam = float(np.median(dg(ba, bb)) / np.median(dchem(ba, bb)))
    lamo = float(np.median(dg(ba, bb)) / np.median(do(ba, bb))) if U is not None else 0.0
    if "_ori_w" in arm:                     # e.g. vc_ori_w050: orientation weight 0.50 x the equal-median default
        lamo *= int(arm.split("_ori_w")[1]) / 100.0   # chosen by vec6 from controls + crystals only
    lamo *= float(os.environ.get("ORI_MULT", "1"))    # orientation multiplier m (vec6; default 1 = equal weight)

    def dp(a, b):
        return np.sqrt(dg(a, b) ** 2 + (lam * dchem(a, b)) ** 2 + (lamo * do(a, b)) ** 2)

    rows.append(summ(aname, dp(pa, pb), dp(ba, bb), dp(wa, wb)))
    T = pd.DataFrame(rows); T["lam"] = round(lam, 4); T["lam_ori"] = round(lamo, 4)
    T.to_csv(os.path.join(OUT, "V2_thresholds_%s.csv" % arm), index=False)
    print("\n== %s (k=%d, lam=%.4f) ==\n%s" % (arm, k, lam, T.to_string(index=False)), flush=True)
    cut = float(T.loc[T.metric == aname, "bg_p1"].iloc[0])

    blocks = [G / np.sqrt(k), lam * CZ.reshape(n, -1) / np.sqrt(CZ.shape[1])] + ([lamo * U / np.sqrt(ku)] if U is not None else [])
    A = np.concatenate(blocks, 1)
    got = np.linalg.norm(A[CHK_A] - A[CHK_B], axis=1)
    assert np.allclose(dp(CHK_A, CHK_B), got, rtol=1e-4, atol=1e-4)

    # receptor-level labels (as d5)
    lab_r = cluster(A, cut)
    pd.DataFrame({"clone_id": ids, "length_class": lclass, aname: lab_r}) \
      .to_csv(os.path.join(OUT, "V2_receptor_labels_%s.csv.gz" % arm), index=False)
    # molecule-level labels (as d6: lambda from the ROUNDED medians written to the threshold
    # table, features in float64 and cast to float32 only for pdist)
    Am = A[rep]
    M = M0.copy()
    M[aname] = cluster(Am, cut)
    del A, Am
    sz = pd.Series(M[aname][M[aname] >= 0]).value_counts()
    print("%s cut %.4f -> %d clusters, %d molecules, largest %d, >=3 members %d"
          % (aname, cut, len(sz), int((M[aname] >= 0).sum()), int(sz.max()), int((sz >= 3).sum())),
          flush=True)
    M.to_csv(os.path.join(OUT, "V2_molecule_labels_%s.csv.gz" % arm), index=False)

    # null ladder (rng re-seeded per arm)
    rng2 = np.random.default_rng(0)
    D = M[M.state.notna()].reset_index(drop=True)
    d = D[D[aname] >= 0]
    s_ = d[aname].value_counts(); d = d[d[aname].isin(s_[s_ >= 3].index)]
    cl, _ = pd.factorize(d[aname]); st, _ = pd.factorize(d.state)
    ncl = cl.max() + 1
    obs = purity(cl, st, ncl)
    r = dict(arm=aname, clusters=ncl, molecules=len(d), observed=round(obs, 4))
    for lab, strata in (("unstratified", np.zeros(len(d), int)),
                        ("mouse", pd.factorize(d.donor)[0]),
                        ("mouse_V", pd.factorize(d.donor.astype(str) + "|" + d.v_A_prot.astype(str)
                                                 + "|" + d.v_B_prot.astype(str))[0])):
        order = np.argsort(strata, kind="stable")
        gs, base = strata[order], st[order]
        nl = np.empty(NPERM)
        for kk in range(NPERM):
            i2 = np.lexsort((rng2.random(len(gs)), gs))
            s = np.empty_like(st); s[order] = base[i2]
            nl[kk] = purity(cl, s, ncl)
        r["null_" + lab] = round(float(nl.mean()), 4)
        r["delta_" + lab] = round(float(obs - nl.mean()), 4)
        r["z_" + lab] = round(float((obs - nl.mean()) / max(nl.std(), 1e-9)), 1)
    pd.DataFrame([r]).to_csv(os.path.join(OUT, "V2_arm_comparison_%s.csv" % arm), index=False)
    print(pd.DataFrame([r]).to_string(index=False), flush=True)

    # per cluster
    d = D[D[aname] >= 0]
    s_ = d[aname].value_counts(); d = d[d[aname].isin(s_[s_ >= 2].index)]
    basef = D.state.value_counts(normalize=True)
    crows = []
    for cid, g in d.groupby(aname):
        top = g.state.value_counts().index[0]
        kk = int((g.state == top).sum())
        crows.append(dict(arm=aname, cluster=int(cid), members=len(g), top_state=top, top_n=kk,
                          frac=round(kk / len(g), 3), mice=g.donor.nunique(),
                          v_pairs=g.v_A_prot.astype(str).add("|").add(g.v_B_prot.astype(str)).nunique(),
                          p=stats.binomtest(kk, len(g), float(basef[top]), alternative="greater").pvalue))
    K = pd.DataFrame(crows)
    pv = K.p.values; o = np.argsort(pv)
    q = np.empty(len(pv))
    q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    K["q_bh"] = q
    K.sort_values("q_bh").to_csv(os.path.join(OUT, "V2_cluster_tests_%s.csv" % arm), index=False)
    print("%s clusters tested %d | q <= 0.15 %d | %.0f s"
          % (aname, len(K), int((K.q_bh <= 0.15).sum()), time.time() - t0), flush=True)
