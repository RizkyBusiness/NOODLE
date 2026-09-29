"""VXV1 (the voxel pipeline, A12.1): build every descriptor for the validity panel.

Receptors: VXV0's 56 minus those failing the complete-CDR3 check (dropped from every descriptor).
Grids (5 configurations, vxv_common.CONFIGS): per receptor a crystal and a model base grid on shared atoms (A3.1), and,
for every pair (i, j), removal of the IMGT positions missing in the other receptor's crystal (A12.9.2), done as
base minus a sparse removal grid. C (crystal-crystal) and M (model-model) 55 x 55 matrices; own crystal-model
distances. Pilot (VX2a 2,000 molecules; background = their pairs minus the control pairs among them, A12.8.1): 1 % cut and
mean squared distance, per configuration.
Vectors: primary <reference arm> geometry-only, full <reference arm>, vc, vg (vxv_common), crystal and model from
the reference method/out/V3_crystal_coords.npz, pilot from the landmark file + D3. Sequence baseline S (per-loop BLOSUM62).
Writes <voxel out>/out/VXV1_descriptors.npz; checks/VXV1_checks.csv. No state data.
usage: python pipeline/code/vxv1_build.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial.distance import pdist, squareform
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected value of the checks (config "expected"; vxpaths.expected)
EXP_REFI = expected("reference_receptor_index")
import vxgrid as vg
from vxv_common import *

TMPG = os.path.join(TMP, "VXV1_pilot.f16.npy")


def sparse(v):
    nz = np.nonzero(v)[0]
    return nz, v[nz]


def build_job(args):
    """removal grids (invert=True) are returned sparse from the worker: dense they would need tens of GB"""
    v, na, nb, used = build_grid(args)
    return (sparse(v) if args[7] else v), na, nb, used


def pilot_background():
    MP = pd.read_csv(os.path.join(OUT, "VX2a_manifest.csv"))
    pil = list(MP.clone_id[MP.kind == "rep"])
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    pos = {c: k for k, c in enumerate(pil)}
    iu = np.triu_indices(len(pil), 1)
    keep = np.ones(len(iu[0]), bool)
    for a, b in zip(NP.clone_a, NP.clone_b):
        if a in pos and b in pos:
            x, y = sorted((pos[a], pos[b]))
            keep[np.where((iu[0] == x) & (iu[1] == y))[0]] = False
    return pil, iu, keep


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VXV1", stop_on_fail=True)
    C.info("status", "A12 VXV1 build; no state data read")
    lm3 = load_lm3()
    R0 = pd.read_csv(os.path.join(OUT, "VXV0_receptors.csv"))
    R = R0[R0.cdr3_complete].reset_index(drop=True)
    n = len(R)
    C.info("receptors in the panel (complete CDR3)", "%d (dropped: %s)" % (n, " ".join(R0.entry[~R0.cdr3_complete])))
    PC = [os.path.join(BEN, "fixed_%s_crystal.pdb" % e) for e in R.entry]
    PM = [os.path.join(BEN, m) for m in R.model]
    SH = [shared_keys(a, b) for a, b in zip(PC, PM)]
    MISS = [frozenset(missing_positions(a, b)) for a, b in zip(PC, PM)]
    iu = np.triu_indices(n, 1)
    pil, piu, pkeep = pilot_background()
    ids = [str(c) for c in L1["clone_id"]]; idx = {c: k for k, c in enumerate(ids)}
    pil_i = np.array([idx[c] for c in pil])
    C.info("pilot background pairs", "%d (removed %d control pairs)" % (int(pkeep.sum()), int((~pkeep).sum())))
    OUTD = {}

    # ------------------------------------------------------------ grids
    for name, (frame, atomset, nch) in CONFIGS.items():
        ts = time.time()
        boxes = boxes_for(frame)
        nvox = [int(np.prod(s)) for _, s in boxes.values()]
        occ_idx = np.concatenate([np.arange(sum(nvox[:k]) * nch, sum(nvox[:k]) * nch + nvox[k]) for k in range(len(nvox))])
        plc = [bench_place(p, frame, lm3) for p in PC]; plm = [bench_place(p, frame, lm3) for p in PM]
        jobs = [(PC[i], plc[i], boxes, atomset, nch, SH[i], None, False) for i in range(n)] + \
               [(PM[i], plm[i], boxes, atomset, nch, SH[i], None, False) for i in range(n)]
        rem_keys = []                                     # (receptor, P) for every non-empty removal
        for i in range(n):
            for P in sorted({MISS[j] for j in range(n) if j != i and MISS[j]}, key=lambda s: sorted(s)):
                rem_keys.append((i, P))
                jobs.append((PC[i], plc[i], boxes, atomset, nch, SH[i], P, True))
                jobs.append((PM[i], plm[i], boxes, atomset, nch, SH[i], P, True))
        with Pool(8) as pool:
            res = pool.map(build_job, jobs, chunksize=8)
        nbad = sum(r[2] for r in res)
        C.add("%s typing complete" % name, nbad, nbad == 0, "== 0")
        Bc = np.array([res[i][0] for i in range(n)]); Bm = np.array([res[n + i][0] for i in range(n)])
        for i in range(n):
            assert res[i][3] == res[n + i][3] == frozenset(k for k in SH[i] if atomset == "whole" or loop_of(k[1]))
        mass = max(abs(float(res[k][0][occ_idx].sum()) / res[k][1] - 1) for k in range(2 * n) if res[k][1])
        C.add("%s mass conservation, %d benchmark base grids" % (name, 2 * n), "%.4f %%" % (100 * mass), mass < 5e-3, "< 0.5 %")
        remc, remm, nmask = {}, {}, 0
        for k, (i, P) in enumerate(rem_keys):
            rc, rm = res[2 * n + 2 * k], res[2 * n + 2 * k + 1]
            if rc[3] != rm[3]:
                raise SystemExit("mask differs between crystal and model side for receptor %d" % i)
            if rc[1]:
                remc[(i, P)], remm[(i, P)] = rc[0], rm[0]
        del res

        def masked_d(Bx, rem, i, j):
            d = Bx[i] - Bx[j]
            if (i, MISS[j]) in rem:
                ix, v = rem[(i, MISS[j])]; d[ix] -= v
            if (j, MISS[i]) in rem:
                ix, v = rem[(j, MISS[i])]; d[ix] += v
            return float(np.sqrt((d * d).sum() / H ** 3))
        Cm, Mm = np.zeros((n, n)), np.zeros((n, n))
        for i, j in zip(*iu):
            Cm[i, j] = Cm[j, i] = masked_d(Bc, remc, i, j)
            Mm[i, j] = Mm[j, i] = masked_d(Bm, remm, i, j)
            nmask += int((i, MISS[j]) in remc or (j, MISS[i]) in remc)
        own = np.sqrt(((Bc - Bm) ** 2).sum(1) / H ** 3)
        C.info("%s pairs with a non-empty mask" % name, "%d / %d" % (nmask, len(iu[0])))
        C.add("%s matrices symmetric, zero diagonal" % name, "", np.allclose(Cm, Cm.T) and np.allclose(Mm, Mm.T)
              and not np.diag(Cm).any() and not np.diag(Mm).any(), "exact")
        if name == "grid_F1L_7ch":
            V5 = pd.read_csv(os.path.join(OUT, "VX5b_crystal_pairs.csv"))
            V5 = V5[(V5.frame == "F1") & (V5.atoms == "loops") & (V5.sigma == 2.0)].set_index("entry").loc[R.entry]
            e = float(np.max(np.abs(own - V5.d.values) / V5.d.values))
            C.add("F1 loops 7ch own-pair distances reproduce VX5b_crystal_pairs.csv", "%.2e rel" % e, e < 1e-3, "< 1e-3")
        del Bc, Bm, remc, remm

        # pilot: cut and mean squared background distance
        pj = [(VXP("structures/%s.pdb") % c, rep_place(idx[c], frame), boxes, atomset, nch, None, None, False) for c in pil]
        F_ = nch * sum(nvox)
        G = np.lib.format.open_memmap(TMPG, mode="w+", dtype=np.float16, shape=(len(pil), F_))
        na = np.zeros(len(pil), int)
        with Pool(8) as pool:
            for k, (v, a, b, _) in enumerate(pool.imap(build_grid, pj, chunksize=8)):
                G[k] = v; na[k] = a
        G.flush()
        rng = np.random.default_rng(81)
        mr = rng.choice(len(pil), 200, replace=False)
        pm_ = max(abs(float(G[k][occ_idx].astype(np.float64).sum()) / na[k] - 1) for k in mr)
        C.add("%s pilot mass conservation, 200 molecules" % name, "%.4f %%" % (100 * pm_), pm_ < 5e-3, "< 0.5 %")
        errs = []
        for a, b in rng.choice(len(pil), (12, 2), replace=False):
            g = float(((G[a].astype(np.float64) - G[b].astype(np.float64)) ** 2).sum()) / H ** 3
            A_ = structure_atoms(VXP("structures/%s.pdb") % pil[a], pj[a][1], atomset, nch)
            B_ = structure_atoms(VXP("structures/%s.pdb") % pil[b], pj[b][1], atomset, nch)
            ex = sum(vg.analytic_d2(xa, wa, xb, wb, SIGMA, channels=range(nch)) for (xa, wa), (xb, wb) in zip(A_, B_))
            errs.append(g / ex - 1)
        C.add("%s pilot analytic grid-free distance, 12 pairs" % name, "%.4f" % np.abs(errs).max(),
              np.abs(errs).max() < 0.02, "< 2 %")
        Gm_ = np.zeros((len(pil), len(pil)))
        for s in range(0, F_, 20000):
            X = np.asarray(G[:, s:min(s + 20000, F_)], dtype=np.float64); Gm_ += X @ X.T
        nr = np.diag(Gm_).copy()
        D2 = np.clip((nr[:, None] + nr[None] - 2 * Gm_) / H ** 3, 0, None)
        bg2 = D2[piu][pkeep]
        OUTD[name] = dict(C=Cm, M=Mm, own=own, cut=float(np.percentile(np.sqrt(bg2), 1)), bg_mean_d2=float(bg2.mean()))
        del G, Gm_, D2
        os.remove(TMPG)
        C.info("%s done" % name, "cut %.4f, own median %.4f (%.1f min)" % (OUTD[name]["cut"], np.median(own), (time.time() - ts) / 60))

    # ------------------------------------------------------------ vectors
    lam, lo = vec_params()
    _ri = [] if MISSING_TAG in str(EXP_REFI) else [int(EXP_REFI)]    # the reference receptor's row (config expected)
    chk = [i for i in [0, 4000] + _ri + [7000, 123] if i < len(ids)]   # rows beyond the landmark file are dropped
    ez = max(float(np.abs(chemz_of(VXP("structures/%s.pdb") % ids[i]) - D3["chemz"][i]).max()) for i in chk)
    C.add("chemistry recipe reproduces D3 chemz (%d models)%s" % (len(chk), "" if _ri else "; %s" % EXP_REFI), "%.2e" % ez, ez < 1e-4, "< 1e-4")
    allf = [vec_feats(L1["arc"][i].astype(float), L1["lm"][i].astype(float)) for i in range(len(ids))]
    VC = np.array([f[0] for f in allf]); UC = np.array([f[1] for f in allf]); VG = np.array([f[2] for f in allf])
    CZ = D3["chemz"].astype(np.float64)
    assert list(D3["clone_id"]) == ids
    # reproduction of the V6 grid row w = 0.5 (A12.8.2a, tolerance 5e-5)
    rng0 = np.random.default_rng(0)
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
    ba, bb = rng0.integers(0, len(ids), 60000), rng0.integers(0, len(ids), 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    def comb(a, b):
        return np.sqrt(((VC[a] - VC[b]) ** 2).sum(-1) / 200 + (lam ** 2) * ((CZ[a] - CZ[b]) ** 2).sum(-1).mean(-1)
                       + (lo ** 2) * ((UC[a] - UC[b]) ** 2).sum(-1) / 20)
    cut6 = float(np.percentile(comb(ba, bb), 1)); S6 = float((comb(pa, pb) <= cut6).mean())
    V3 = np.load(VXP("reference/out/V3_crystal_coords.npz"), allow_pickle=True)
    cr = []
    for j in range(len(V3["entry"])):
        vc_, uc_, _ = vec_feats(V3["arc_c"][j].astype(float), V3["lm_c"][j].astype(float))
        vm_, um_, _ = vec_feats(V3["arc_m"][j].astype(float), V3["lm_m"][j].astype(float))
        cr.append(np.sqrt(((vc_ - vm_) ** 2).sum() / 200 + lo ** 2 * ((uc_ - um_) ** 2).sum() / 20))
    cr = np.array(cr); F6, C6 = float(np.median(cr) / cut6), float((cr <= cut6).mean())
    V6 = pd.read_csv(VXP("reference/out/V6_orientation_weight_grid.csv")).set_index("w").loc[0.5]
    for nm_, got, tgt in (("cut", cut6, V6.cut), ("S (controls within cut)", S6, V6.S_sensitivity),
                          ("F (crystal error over cut)", F6, V6.F_crystal_over_cut), ("C (crystal within cut)", C6, V6.C_crystal_within_cut)):
        C.add("vector reproduces V6 w=0.5 %s" % nm_, "%.5f vs %.4f" % (got, tgt), abs(got - tgt) < 5e-5, "< 5e-5")

    vrow = {e: k for k, e in enumerate(V3["entry"])}
    fc = [vec_feats(V3["arc_c"][vrow[e]].astype(float), V3["lm_c"][vrow[e]].astype(float)) for e in R.entry]
    fm = [vec_feats(V3["arc_m"][vrow[e]].astype(float), V3["lm_m"][vrow[e]].astype(float)) for e in R.entry]
    zc = [chemz_of(p) for p in PC]; zm = [chemz_of(p) for p in PM]

    def vdist(kind, f1, z1, f2, z2):
        dv = ((f1[0] - f2[0]) ** 2).sum() / 200; do = ((f1[1] - f2[1]) ** 2).sum() / 20
        if kind == "vec_vcori_geo":
            return np.sqrt(dv + lo ** 2 * do)
        if kind == "vec_vcori_full":
            return np.sqrt(dv + lam ** 2 * ((z1 - z2) ** 2).sum(-1).mean() + lo ** 2 * do)
        if kind == "vec_vc":
            return np.sqrt(dv)
        return np.sqrt(((f1[2] - f2[2]) ** 2).sum() / 400)
    pf = [(VC[i], UC[i], VG[i]) for i in pil_i]; pz = [CZ[i] for i in pil_i]
    for kind in VECTORS:
        Cm, Mm = np.zeros((n, n)), np.zeros((n, n))
        for i, j in zip(*iu):
            Cm[i, j] = Cm[j, i] = vdist(kind, fc[i], zc[i], fc[j], zc[j])
            Mm[i, j] = Mm[j, i] = vdist(kind, fm[i], zm[i], fm[j], zm[j])
        own = np.array([vdist(kind, fc[i], zc[i], fm[i], zm[i]) for i in range(n)])
        # pilot background, vectorised
        A1 = np.array([p[0] for p in pf]) / np.sqrt(200); A2 = np.array([p[1] for p in pf]) * lo / np.sqrt(20)
        A3 = np.array([p[2] for p in pf]) / np.sqrt(400); A4 = np.array(pz).reshape(len(pz), -1) * lam / np.sqrt(20)
        X = {"vec_vcori_geo": np.hstack([A1, A2]), "vec_vcori_full": np.hstack([A1, A4, A2]), "vec_vc": A1,
             "vec_vg": A3}[kind]
        bgd = squareform(pdist(X))[piu][pkeep]
        chk_ = squareform(pdist(np.vstack([np.hstack([fc[0][0] / np.sqrt(200), fc[0][1] * lo / np.sqrt(20)]),
                                          np.hstack([fc[1][0] / np.sqrt(200), fc[1][1] * lo / np.sqrt(20)])])))[0, 1]
        if kind == "vec_vcori_geo":
            C.add("vector pilot feature layout == pairwise definition", "%.2e" % abs(chk_ - Cm[0, 1]),
                  abs(chk_ - Cm[0, 1]) < 1e-9, "< 1e-9")
        OUTD[kind] = dict(C=Cm, M=Mm, own=own, cut=float(np.percentile(bgd, 1)), bg_mean_d2=float((bgd ** 2).mean()))
        C.info("%s done" % kind, "cut %.4f, own median %.4f" % (OUTD[kind]["cut"], np.median(own)))

    # ------------------------------------------------------------ sequence baseline and strata
    al = aligner()
    seqs = [loop_seqs(p) for p in PM]
    def sc(a, b):
        return sum(al.score(a[k], b[k]) for k in a if a[k] and b[k])
    self_ = [sc(s, s) for s in seqs]
    S = np.zeros((n, n))
    for i, j in zip(*iu):
        S[i, j] = S[j, i] = 1 - sc(seqs[i], seqs[j]) / np.sqrt(self_[i] * self_[j])
    C.add("sequence baseline: self-scores positive", "min %.1f" % min(self_), min(self_) > 0, "> 0")
    VID = np.zeros((n, n, 2))
    for i, j in zip(*iu):
        VID[i, j] = VID[j, i] = vregion_identity(PM[i], PM[j])
    np.savez(os.path.join(OUT, "VXV1_descriptors.npz"), entry=np.array(R.entry), model=np.array(R.model),
             names=np.array(list(OUTD)), S=S, vregion_identity=VID,
             species=np.array(R.species), bound=np.array(R.bound_pMHC),
             **{"%s__%s" % (k, f): v[f] for k, v in OUTD.items() for f in ("C", "M", "own", "cut", "bg_mean_d2")})
    C.info("total time", "%.1f min" % ((time.time() - t0) / 60))
    C.write()
