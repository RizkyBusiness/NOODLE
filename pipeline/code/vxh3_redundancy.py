"""VXH3 (the voxel pipeline, A12.3 / A12.4 rule H2): is the hinge redundant with germline and CDR3 length? No state data.

- 5-fold CV R^2 of h (123 x 3, flattened) on one-hot V pair (C1w v_A_prot|v_B_prot) + CDR3a and CDR3b lengths, ridge
  alpha 1.0 (closed form, unpenalised intercept), folds from default_rng(0), over the M molecules. R^2 pooled over all
  outputs: 1 - sum SSE / sum SST, SST about the training-fold mean. Also V pair only and lengths only (report).
  H2: redundant if R^2 >= 0.90.
- Spearman(d_hinge, D_vc_ori_w050) on the reference cluster-test procedure's 60,000 background pairs (the existing arm's full distance).
- Spearman(hinge angle to the reference, CDR3a length) over the molecules (Dunbar et al. 2014 predict positive).
- d_hinge within the existing <reference arm> clusters (>= 2 molecules; only clone_id and label columns read) against
  size-matched random groups (200 draws per cluster, seed 0).
Writes out/VXH3_summary.json, out/VXH3_cluster_hinge.csv; checks/VXH3_checks.csv.
usage: python pipeline/code/vxh3_redundancy.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys
import numpy as np, pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
EXP_NBG = expected("n_background_pairs")
EXP_CUT_REF = expected("reference_cut_vc_ori_w050")
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column

import vxv_common as VC


def ridge_cv(X, Y, alpha=1.0, k=5, seed=0):
    idx = np.random.default_rng(seed).permutation(len(X)); folds = np.array_split(idx, k)
    sse = sst = 0.0
    for f in folds:
        tr = np.setdiff1d(idx, f)
        xm, ym = X[tr].mean(0), Y[tr].mean(0)
        Xt, Yt = X[tr] - xm, Y[tr] - ym
        Wt = np.linalg.solve(Xt.T @ Xt + alpha * np.eye(X.shape[1]), Xt.T @ Yt)
        pred = (X[f] - xm) @ Wt + ym
        sse += ((Y[f] - pred) ** 2).sum(); sst += ((Y[f] - ym) ** 2).sum()
    return float(1 - sse / sst)


if __name__ == "__main__":
    C = Checks("VXH3", stop_on_fail=True)
    C.info("status", "A12 VXH3 redundancy; no state data read")
    P = np.load(os.path.join(OUT, "VXH1_pose.npz"), allow_pickle=True)
    nr = int(P["n_repertoire"]); ids = [str(x) for x in P["names"][:nr]]
    FR = np.load(os.path.join(OUT, "VX1_frames.npz"), allow_pickle=True); mol = np.where(FR["is_molecule"])[0]
    C.add("molecules", len(mol), len(mol) == EXP_NMOL, "== %s" % EXP_NMOL)
    ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz"), usecols=["clone_id", "v_A_prot", "v_B_prot"]).set_index("clone_id").reindex(ids)
    SL = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz"), usecols=["clone_id", "cdr3len_A", "cdr3len_B"]).set_index("clone_id").reindex(ids)
    vp = (ID.v_A_prot.astype(str) + "|" + ID.v_B_prot.astype(str)).values[mol]
    lens = SL[["cdr3len_A", "cdr3len_B"]].values[mol].astype(float)
    OH = pd.get_dummies(pd.Series(vp)).values.astype(float)
    Y = P["h"][mol].reshape(len(mol), -1)
    r2 = ridge_cv(np.hstack([OH, lens]), Y); r2v = ridge_cv(OH, Y); r2l = ridge_cv(lens, Y)
    C.info("5-fold CV R^2 of h: V pair + CDR3 lengths / V pair only / lengths only",
           "%.4f / %.4f / %.4f (%d V pairs)" % (r2, r2v, r2l, OH.shape[1]))
    H2 = r2 >= 0.90
    C.info("RULE H2: R^2 >= 0.90 (hinge essentially germline)", "%.4f -> %s" % (r2, "REDUNDANT" if H2 else "not redundant"))

    # vs the existing arm's full distance on the 60,000 background pairs
    L1 = VC.L1; lam, lo = VC.vec_params()
    allf = [VC.vec_feats(L1["arc"][i].astype(float), L1["lm"][i].astype(float)) for i in range(nr)]
    Vc = np.array([f[0] for f in allf]); Uc = np.array([f[1] for f in allf]); CZ = VC.D3["chemz"].astype(np.float64)
    rng0 = np.random.default_rng(0); ba, bb = rng0.integers(0, nr, 60000), rng0.integers(0, nr, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    Dv = np.sqrt(((Vc[ba] - Vc[bb]) ** 2).sum(-1) / 200 + lam ** 2 * ((CZ[ba] - CZ[bb]) ** 2).sum(-1).mean(-1)
                 + lo ** 2 * ((Uc[ba] - Uc[bb]) ** 2).sum(-1) / 20)
    C.add("vc_ori_w050 background cut reproduced", "%.4f" % np.percentile(Dv, 1), abs(np.percentile(Dv, 1) - EXP_CUT_REF) < 5e-5, "%s +- 5e-5" % EXP_CUT_REF)
    Hh = P["h"]; dhb = np.sqrt(((Hh[ba] - Hh[bb]) ** 2).sum(-1).mean(-1))
    rho = spearmanr(dhb, Dv).correlation
    C.info("Spearman(d_hinge, D_vc_ori_w050), {:,} pairs".format(EXP_NBG), "%.4f" % rho)
    ang = P["hinge_angle_to_ref_deg"][mol]
    ra = spearmanr(ang, lens[:, 0]).correlation; rb = spearmanr(ang, lens[:, 1]).correlation
    C.info("Spearman(hinge angle to reference, CDR3a length) / CDR3b", "%.4f / %.4f" % (ra, rb))

    # within existing clusters vs size-matched random groups
    Lb = pd.read_csv(VXP(_N10F), usecols=["clone_id", _N10C])
    idx = {c: k for k, c in enumerate(ids)}
    Lb = Lb[Lb.prop_vc_ori_w050 >= 0]
    molset = np.array(mol); rng = np.random.default_rng(0); rows = []
    for cl, g in Lb.groupby(_N10C):
        mem = np.array([idx[c] for c in g.clone_id]); k = len(mem)
        if k < 2:
            continue
        iu = np.triu_indices(k, 1)
        obs = float(np.sqrt(((Hh[mem][:, None] - Hh[mem][None]) ** 2).sum(-1).mean(-1))[iu].mean())
        rnd = []
        for _ in range(200):
            r = rng.choice(molset, k, replace=False)
            rnd.append(np.sqrt(((Hh[r][:, None] - Hh[r][None]) ** 2).sum(-1).mean(-1))[iu].mean())
        rows.append(dict(cluster=int(cl), size=k, mean_dhinge=obs, random_mean=float(np.mean(rnd)),
                         ratio=obs / float(np.mean(rnd)), pct_below=float((np.array(rnd) <= obs).mean())))
    CL = pd.DataFrame(rows); CL.to_csv(os.path.join(OUT, "VXH3_cluster_hinge.csv"), index=False)
    C.info("d_hinge within vc_ori_w050 clusters vs size-matched random: median ratio (IQR); clusters",
           "%.3f (%.3f-%.3f); %d" % (CL.ratio.median(), CL.ratio.quantile(.25), CL.ratio.quantile(.75), len(CL)))
    C.info("share of clusters with within-cluster d_hinge below the random median", "%.3f" % (CL.pct_below < 0.5).mean())
    save_json(dict(cv_r2=r2, cv_r2_vpair=r2v, cv_r2_lengths=r2l, H2_redundant=bool(H2), spearman_dhinge_vs_vcori=rho,
                   spearman_angle_cdr3a=ra, spearman_angle_cdr3b=rb, cluster_ratio_median=float(CL.ratio.median()),
                   n_clusters=len(CL)), os.path.join(OUT, "VXH3_summary.json"))
    C.write()
