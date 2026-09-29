"""VX4 (the voxel pipeline): threshold, clustering, controls - procedure mirrors the reference method/code/vec2_cluster_tests.py.

Pairs drawn exactly as the reference cluster-test procedure (default_rng(0), same order): the control pairs; 60,000 random background pairs over the
N receptors (self-pairs dropped); 40,000 within-length-class draws; then the reference cluster-test procedure's 200 check pairs (not used).
Cut = 1st percentile of the background grid distances, rounded to 4 decimals as the reference cluster-test procedure writes and re-reads it.
Molecules: the M the reference cluster-test procedure-rule representatives; complete linkage; fcluster at the cut; singletons dropped;
labels in the reference cluster-test procedure's order (cluster size, descending). Control statistics: share of the control pairs within the cut (the reference cluster-test procedure's
sens_at_bg_p1) and the median percentile of the control distances in the background (d_chem review).
Stability: ARI (Hubert & Arabie 1985) of the partitions at cut x0.98 and x1.02 against the reported one, singletons
counted as their own clusters. No state data are used here.
usage: python pipeline/code/vx4_cluster.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.special import comb
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
EXP_NBG = expected("n_background_pairs")

ARM = "prop_voxel"


def summ(metric, p, b, w_):
    return dict(metric=metric, pos_median=round(float(np.median(p)), 4), pos_p90=round(float(np.percentile(p, 90)), 4),
                bg_p1=round(float(np.percentile(b, 1)), 4), bg_median=round(float(np.median(b)), 4),
                within_p1=round(float(np.percentile(w_, 1)), 4), within_median=round(float(np.median(w_)), 4),
                sens_at_bg_p1=round(float((p <= np.percentile(b, 1)).mean()), 4),
                sens_at_within_p1=round(float((p <= np.percentile(w_, 1)).mean()), 4))


def labels_at(Z, cut, n):
    f = fcluster(Z, t=cut, criterion="distance")
    cnt = pd.Series(f).value_counts()
    out = np.full(n, -1)
    for j, c in enumerate(cnt[cnt >= 2].index):
        out[f == c] = j
    return out


def ari(a, b):
    """adjusted Rand index; label -1 = singleton, each its own cluster."""
    a = np.where(a < 0, -1 - np.arange(len(a)), a); b = np.where(b < 0, -1 - np.arange(len(b)), b)
    _, ai = np.unique(a, return_inverse=True); _, bi = np.unique(b, return_inverse=True)
    pairs = pd.Series(ai.astype(np.int64) * (bi.max() + 1) + bi).value_counts().values
    sij = comb(pairs, 2).sum(); sa = comb(np.bincount(ai), 2).sum(); sb = comb(np.bincount(bi), 2).sum()
    exp = sa * sb / comb(len(a), 2)
    return float((sij - exp) / (0.5 * (sa + sb) - exp))


if __name__ == "__main__":
    C = Checks("VX4")
    J = json.load(open(os.path.join(OUT, "VX3_distances.json")))
    D = np.load(os.path.join(ROOT, J["receptor_matrix"]), mmap_mode="r")
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]
    n = len(ids); idx = {c: i for i, c in enumerate(ids)}

    # ---------------------------------------------------------------- pair sets, exactly as the reference cluster-test procedure
    R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
    lclass = (R.cdr3_A.str.len().astype(str) + "-" + R.cdr3_B.str.len().astype(str)).values
    rng = np.random.default_rng(0)
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
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
    C.info("pairs: control / background / within-length-class", "%d / %d / %d" % (len(pa), len(ba), len(wa)))
    dp_, db_, dw_ = (np.asarray(D[x, y], dtype=np.float64) for x, y in ((pa, pb), (ba, bb), (wa, wb)))
    T = pd.DataFrame([summ(ARM, dp_, db_, dw_)])
    T.to_csv(os.path.join(OUT, "VX4_thresholds.csv"), index=False)
    cut = float(T.bg_p1.iloc[0])
    C.info("cut (1st pct of background, 4 dp as vec2)", "%.4f" % cut)
    C.info("background median / within-class p1 / within-class median", "%.4f / %.4f / %.4f"
           % (T.bg_median.iloc[0], T.within_p1.iloc[0], T.within_median.iloc[0]))
    bgs = np.sort(db_)
    pct = 100 * np.searchsorted(bgs, dp_, side="right") / len(bgs)
    C.info("control pairs: share within cut (comparable to existing arms)", "%.4f" % T.sens_at_bg_p1.iloc[0])
    C.info("control pairs: median percentile in background (d_chem statistic)", "%.4f %% (p90 %.3f %%)"
           % (np.median(pct), np.percentile(pct, 90)))
    C.info("control pairs: median distance / cut", "%.4f" % (np.median(dp_) / cut))
    pd.DataFrame({"clone_a": NP.clone_a, "clone_b": NP.clone_b, "d": dp_, "over_cut": dp_ / cut, "bg_pct": pct}) \
      .to_csv(os.path.join(OUT, "VX4_control_pairs.csv"), index=False)
    C.add("background pairs after dropping self-pairs (vec2 draw, seed 0)", len(ba), len(ba) == EXP_NBG, "== %s" % EXP_NBG)

    # ---------------------------------------------------------------- molecules and clustering
    ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz"))
    SL = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz"))
    ST = pd.read_csv(os.path.join(STB, "_receptor_states_all.csv.gz")).drop_duplicates("clone_key")
    M0 = (pd.DataFrame({"clone_id": ids}).merge(SL[["clone_id", "clone_key"]], on="clone_id", how="left")
          .merge(ID.drop_duplicates("clone_key")[["clone_key", "prot_key", "v_A_prot", "v_B_prot"]], on="clone_key", how="left")
          .merge(ST[["clone_key", "state", "donor"]], on="clone_key", how="left"))
    keep = M0.prot_key.notna().values & ~M0.prot_key.duplicated().values
    rep = np.where(keep)[0]
    C.add("molecule rows == VX3 molecule rows", len(rep), list(rep) == J["molecule_rows"], "identical, %s" % EXP_NMOL)
    M = M0.iloc[rep].reset_index(drop=True)
    Dm = np.load(os.path.join(ROOT, J["molecule_condensed"]))
    Z = linkage(Dm, method="complete")
    lab = labels_at(Z, cut, len(M))
    M[ARM] = lab
    M["length_class"] = lclass[rep]
    M.drop(columns=["state", "donor"]).to_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"), index=False)

    Sq = squareform(Dm)
    worst = 0.0
    for c in range(lab.max() + 1):
        mm = np.where(lab == c)[0]
        worst = max(worst, float(Sq[np.ix_(mm, mm)].max()))
    C.add("every within-cluster pair within the cut", "max intra-cluster d %.4f vs cut %.4f" % (worst, cut),
          worst <= cut + 1e-9, "<= cut")
    sz = pd.Series(lab[lab >= 0]).value_counts()
    C.info("clusters (>= 2 members) / molecules clustered / largest / clusters >= 3",
           "%d / %d / %d / %d" % (len(sz), int((lab >= 0).sum()), int(sz.max()), int((sz >= 3).sum())))
    C.info("cluster size distribution (size: count)", ", ".join("%d: %d" % (k, v) for k, v in
           sz.value_counts().sort_index().items()))
    for f_ in (0.98, 1.02):
        l2 = labels_at(Z, cut * f_, len(M))
        s2 = pd.Series(l2[l2 >= 0]).value_counts()
        C.info("stability: cut x%.2f -> clusters / clustered / ARI vs reported" % f_,
               "%d / %d / %.4f" % (len(s2), int((l2 >= 0).sum()), ari(lab, l2)))
    # context from the existing arms (their V2 files), for the record only
    E = []
    for a in ("vc_ori_w050", "vc_ori", "vc"):
        t = pd.read_csv(VXP("reference/out/V2_thresholds_%s.csv") % a).set_index("metric").loc["prop_" + a]
        l = pd.read_csv(VXP("reference/out/V2_molecule_labels_%s.csv.gz") % a)["prop_" + a]
        s = l[l >= 0].value_counts()
        E.append("%s: sens %.3f, clusters %d, clustered %d" % (a, t.sens_at_bg_p1, len(s), int((l >= 0).sum())))
    C.info("existing arms (context)", " | ".join(E))
    save_json(dict(arm=ARM, cut=cut, clusters=int(len(sz)), clustered=int((lab >= 0).sum()),
                   largest=int(sz.max()), clusters_ge3=int((sz >= 3).sum()),
                   control_within_cut=float(T.sens_at_bg_p1.iloc[0]), control_median_bg_pct=float(np.median(pct))),
              os.path.join(OUT, "VX4_summary.json"))
    C.write()
