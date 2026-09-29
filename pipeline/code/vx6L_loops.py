"""VX6-L (the voxel pipeline): EXPLORATORY loops-only view after the VX5 hard stop (amendment A8; not pre-registered).

Grids rebuilt from the heavy atoms of IMGT CDR1 27-38, CDR2 56-65, HV4 81-86 and CDR3 105-117 (both chains, insertions
included); frame, box, typing, sigma 2.0 A, voxel 1.0 A and channels 1-7 exactly as production. Then the VX3 distance
route (float64 blocks), the VX4 procedure (the reference cluster-test procedure's pairs, 1 % cut, complete linkage), the crystal floor on loop atoms,
and the VX6 confound diagnostics for this mask and for the whole molecule. No state data are used.
usage: python pipeline/code/vx6L_loops.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NBP = expected("n_benchmark_pairs")
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column

import vxgrid as vg
import vx2a_pilot as PA
from vx4_cluster import summ, labels_at, ari

LOOPS = {"CDR1": (27, 38), "CDR2": (56, 65), "HV4": (81, 86), "CDR3": (105, 117)}
SIGMA, H = 2.0, 1.0
GRID = os.path.join(TMP, "VX6L_h1.0_s2.00.f16.npy")
NF = vg.NPRIMARY * int(np.prod(PA.SHAPE))


def loop_of(num):
    for k, (a, b) in LOOPS.items():
        if a <= num <= b:
            return k
    return None


def loop_structure(path, Rb, tb, keep=None):
    xyz, meta = parse_heavy(path)
    W, bad = vg.type_atoms(meta)                       # typed on the whole chain, then selected
    m = np.array([loop_of(a[1]) is not None and (keep is None or tuple(a) in keep) for a in meta])
    return xyz[m] @ Rb + tb, W[m], [a for a, k in zip(meta, m) if k], len(bad)


def build_one(i):
    X, W, meta, nb = loop_structure(VXP("structures/%s.pdb") % PA.IDS[i], PA.FR["R_box"][i], PA.FR["t_box"][i])
    seg = {(a[0], loop_of(a[1])) for a in meta}
    return vg.build(X, W, SIGMA, H, PA.LO, PA.SHAPE).astype(np.float16), len(X), len(seg), nb


def gram_distances(G):
    n = G.shape[0]; flat = G.reshape(n, -1); Gm = np.zeros((n, n))
    for s in range(0, NF, 20000):
        X = np.asarray(flat[:, s:min(s + 20000, NF)], dtype=np.float64)
        Gm += X @ X.T
    nr = np.diag(Gm).copy()
    D2 = (nr[:, None] + nr[None] - 2 * Gm) / H ** 3
    return D2, nr


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VX6L")
    C.info("status", "EXPLORATORY, not pre-registered (A8); VX5 hard stop stands; no state data used")
    ids = PA.IDS; n = len(ids); idx = {c: i for i, c in enumerate(ids)}
    ismol = PA.FR["is_molecule"]

    # ------------------------------------------------------------ A. loops-only grids
    G = np.lib.format.open_memmap(GRID, mode="w+", dtype=np.float16, shape=(n, vg.NCH) + PA.SHAPE)
    natoms, nseg, nbad = np.zeros(n, int), np.zeros(n, int), 0
    with Pool(8) as pool:
        for i, (g, na, ns, nb) in enumerate(pool.imap(build_one, range(n), chunksize=8)):
            G[i] = g; natoms[i] = na; nseg[i] = ns; nbad += nb
    G.flush()
    C.info("build time", "%.1f min" % ((time.time() - t0) / 60))
    full = pd.read_csv(os.path.join(OUT, "VX2c_manifest.csv")).heavy_atoms.values
    C.info("loop heavy atoms per receptor: median (min-max) / share of all heavy atoms",
           "%d (%d-%d) / %.1f %%" % (np.median(natoms), natoms.min(), natoms.max(), 100 * natoms.sum() / full.sum()))
    C.add("typing complete (whole chains)", nbad, nbad == 0, "== 0")
    C.add("all 8 loop segments (4 loops x 2 chains) present in every receptor", "min %d" % nseg.min(),
          nseg.min() == 8, "== 8")
    rng = np.random.default_rng(61)
    mols = np.where(ismol)[0]
    rows = rng.choice(mols, 200, replace=False)
    rel = max(abs(float(G[r, 0].astype(np.float64).sum()) / natoms[r] - 1) for r in rows)
    C.add("mass conservation, 200 molecules (loop atoms)", "%.4f %%" % (100 * rel), rel < 5e-3, "< 0.5 %")
    errs = []
    for a, b in rng.choice(mols, (12, 2), replace=False):
        XA, WA, _, _ = loop_structure(VXP("structures/%s.pdb") % ids[a], PA.FR["R_box"][a], PA.FR["t_box"][a])
        XB, WB, _, _ = loop_structure(VXP("structures/%s.pdb") % ids[b], PA.FR["R_box"][b], PA.FR["t_box"][b])
        errs.append(vg.grid_d2(G[a], G[b], H) / vg.analytic_d2(XA, WA, XB, WB, SIGMA) - 1)
    C.add("analytic grid-free distance, 12 pairs", "%.4f (%+.4f..%+.4f)" % (np.abs(errs).max(), min(errs), max(errs)),
          np.abs(errs).max() < 0.02, "< 2 %")

    # ------------------------------------------------------------ B. distances (VX3 route)
    D2, nr = gram_distances(G)
    off = D2[~np.eye(n, dtype=bool)]; mn = float(off.min()); del off
    C.add("no negative squared distances", "min off-diagonal D^2 %.3e" % mn, mn > -1e-6, "> -1e-6")
    np.fill_diagonal(D2, 0.0)
    DL = np.sqrt(np.clip(D2, 0, None)); del D2
    DL = 0.5 * (DL + DL.T)
    sub = np.sort(rng.choice(mols, 200, replace=False))
    ref = pdist(np.asarray(G.reshape(n, -1)[sub, :NF], dtype=np.float64)) / H ** 1.5
    r_ = float(np.max(np.abs(squareform(DL[np.ix_(sub, sub)], checks=False) - ref) / ref))
    C.add("blocked Gram == scipy pdist, 200-molecule subset", "%.2e" % r_, r_ < 1e-4, "< 1e-4")
    np.save(os.path.join(TMP, "VX6L_D_receptors.f32.npy"), DL.astype(np.float32))

    # ------------------------------------------------------------ C. VX4 procedure
    R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
    rng0 = np.random.default_rng(0)
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
    ba, bb = rng0.integers(0, n, 60000), rng0.integers(0, n, 60000)
    ok = ba != bb; ba, bb = ba[ok], bb[ok]
    DW = np.load(os.path.join(TMP, "VX3_D_receptors.f32.npy"), mmap_mode="r")
    cutW = json.load(open(os.path.join(OUT, "VX4_summary.json")))["cut"]
    dbW = np.asarray(DW[ba, bb], dtype=np.float64)
    C.add("pair draw reproduces the whole-molecule VX4 cut", "%.4f vs %.4f" % (round(float(np.percentile(dbW, 1)), 4), cutW),
          round(float(np.percentile(dbW, 1)), 4) == cutW, "identical")
    dp, db = DL[pa, pb], DL[ba, bb]
    cut = round(float(np.percentile(db, 1)), 4)
    bgs = np.sort(db); pct = 100 * np.searchsorted(bgs, dp, side="right") / len(bgs)
    C.info("loops cut / background median", "%.4f / %.4f" % (cut, float(np.median(db))))
    C.info("control pairs: share within cut / median bg percentile / median d over cut",
           "%.4f / %.4f %% / %.4f" % (float((dp <= cut).mean()), float(np.median(pct)), float(np.median(dp) / cut)))

    ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz")).set_index("clone_id").reindex(ids)
    vpair = (ID.v_A_prot.astype(str) + "|" + ID.v_B_prot.astype(str)).values
    Mw = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"))
    rep = np.array([idx[c] for c in Mw.clone_id])
    assert list(rep) == list(mols)
    Dm = squareform(DL[np.ix_(rep, rep)], checks=False).astype(np.float64)
    Z = linkage(Dm, method="complete")
    lab = labels_at(Z, cut, len(rep))
    Sq = squareform(Dm)
    worst = max(float(Sq[np.ix_(np.where(lab == c)[0], np.where(lab == c)[0])].max()) for c in range(lab.max() + 1))
    C.add("every within-cluster pair within the cut", "max %.4f vs cut %.4f" % (worst, cut), worst <= cut + 1e-9, "<= cut")
    del Sq
    sz = pd.Series(lab[lab >= 0]).value_counts()
    C.info("clusters / molecules clustered / largest / clusters >= 3",
           "%d / %d / %d / %d" % (len(sz), int((lab >= 0).sum()), int(sz.max()), int((sz >= 3).sum())))
    for f_ in (0.98, 1.02):
        l2 = labels_at(Z, cut * f_, len(rep))
        C.info("stability: cut x%.2f -> clusters / ARI vs reported" % f_,
               "%d / %.4f" % (len(pd.Series(l2[l2 >= 0]).unique()), ari(lab, l2)))
    Mw["prop_voxel_loops"] = lab
    Mw.to_csv(os.path.join(OUT, "VX6L_molecule_labels.csv.gz"), index=False)

    # ------------------------------------------------------------ D. crystal floor on loop atoms
    lm3 = load_lm3(); B = lm3["B"].sort_values("entry")
    cf = []
    for _, r in B.iterrows():
        pc, pm = os.path.join(PA.BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(PA.BEN, r.model)
        keep = {tuple(a) for a in parse_heavy(pc)[1]} & {tuple(a) for a in parse_heavy(pm)[1]}
        g = []
        for p in (pc, pm):
            Rb, tb = PA.bench_frame(p, lm3["parse"], lm3["LMPOS"])
            X, W, _, _ = loop_structure(p, Rb, tb, keep)
            g.append(vg.build(X, W, SIGMA, H, PA.LO, PA.SHAPE).astype(np.float16))
        cf.append(dict(entry=r.entry, model=r.model, d_loops=np.sqrt(vg.grid_d2(g[0], g[1], H))))
    CF = pd.DataFrame(cf); CF["over_cut"] = CF.d_loops / cut
    CF.to_csv(os.path.join(OUT, "VX6L_crystal_floor.csv"), index=False)
    rL = float(CF.d_loops.median() / cut)
    C.info("crystal-model error over cut, loops (shared loop atoms)", "%.4f (whole-molecule value and the pre-registered limit as context)" % rL)
    C.info("crystal pairs within the loops cut", "%d / %s" % (int((CF.over_cut <= 1).sum()), EXP_NBP))

    # ------------------------------------------------------------ E. VX6 confound diagnostics, both masks
    la, lb = R.cdr3len_A.values, R.cdr3len_B.values
    dA, dB = np.abs(la[ba] - la[bb]), np.abs(lb[ba] - lb[bb])
    sameV = vpair[ba] == vpair[bb]
    Mex = pd.read_csv(VXP(_N10F))[["clone_id", _N10C]]
    rows_ = []
    for name, d, labs in (("whole molecule", dbW, Mw.prop_voxel.values), ("loops only", db, lab)):
        o = dict(mask=name, cut=cutW if name == "whole molecule" else cut, bg_median=float(np.median(d)))
        for k, x in (("alpha", dA), ("beta", dB), ("total", dA + dB)):
            o["spearman_dlen_" + k] = float(spearmanr(d, x).correlation)
        m0 = (dA + dB) == 0
        o["median_len_matched"], o["median_len_mismatched"] = float(np.median(d[m0])), float(np.median(d[~m0]))
        o["median_same_Vpair"], o["median_diff_Vpair"] = float(np.median(d[sameV])), float(np.median(d[~sameV]))
        vp = pd.Series(vpair[rep]); cl = pd.Series(labs)
        g = vp[cl >= 0].groupby(cl[cl >= 0]).nunique()
        o["clusters"] = int(len(g)); o["share_single_Vpair_clusters"] = float((g == 1).mean())
        rows_.append(o)
    # context: the existing primary arm
    E = Mw[["clone_id"]].merge(Mex, on="clone_id", how="left")
    el = E[_N10C].fillna(-1).astype(int).values
    ge = pd.Series(vpair[rep])[el >= 0].groupby(pd.Series(el)[el >= 0]).nunique()
    rows_.append(dict(mask="reference method (context)", clusters=int(len(ge)),
                      share_single_Vpair_clusters=float((ge == 1).mean())))
    DG = pd.DataFrame(rows_)
    DG.to_csv(os.path.join(OUT, "VX6L_diagnostics.csv"), index=False)
    for _, o in DG.iterrows():
        C.info("diagnostics: %s" % o["mask"], "; ".join("%s=%.4g" % (k, v) for k, v in o.items()
                                                        if k != "mask" and pd.notna(v)))
    C.info("same-V-pair background pairs", "%d of %d" % (int(sameV.sum()), len(sameV)))
    C.info("ARI loops vs whole molecule", "%.4f" % ari(lab, Mw.prop_voxel.values))
    C.info("ARI loops vs existing vc_ori_w050", "%.4f" % ari(lab, el))
    C.info("ARI whole molecule vs existing vc_ori_w050", "%.4f" % ari(Mw.prop_voxel.values, el))
    save_json(dict(status="exploratory (A8)", cut=cut, clusters=int(len(sz)), clustered=int((lab >= 0).sum()),
                   control_within_cut=float((dp <= cut).mean()), control_median_bg_pct=float(np.median(pct)),
                   crystal_ratio=rL, minutes=round((time.time() - t0) / 60, 1)), os.path.join(OUT, "VX6L_summary.json"))
    C.write()
