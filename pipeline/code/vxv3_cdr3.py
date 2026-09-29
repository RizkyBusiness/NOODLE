"""VXV3 (the voxel pipeline, A12.10.2): report-only CDR3-only grid check. It cannot change V2 or any other rule.

Same construction as the VXV primary (vxv1_build.py) except the atom set: heavy atoms of IMGT CDR3 105-117 on both chains;
F1 frame, VX5b F1 boxes, sigma 2.0 A, 1.0 A voxel, shared atoms, A12.9.2 pair masks; occupancy and channels 1-7.
Compared with the VXV vector primary (<reference arm> geometry, from VXV1) using the VXV2 statistics.
Writes out/VXV3_cdr3.csv, out/VXV3_summary.json; checks/VXV3_checks.csv. No state data.
usage: python pipeline/code/vxv3_cdr3.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
import vxgrid as vg
from vxv_common import *
from vxv1_build import build_job, pilot_background
from vxv2_panel import pspear

CFG3 = {"grid_F1C3_occ": ("F1", "cdr3", 1), "grid_F1C3_7ch": ("F1", "cdr3", 7)}
TMPG = os.path.join(TMP, "VXV3_pilot.f16.npy")
NB, SEED = 2000, 0

if __name__ == "__main__":
    C = Checks("VXV3", stop_on_fail=True)
    C.info("status", "A12.10.2 report-only CDR3-only check; cannot change any rule; no state data read")
    lm3 = load_lm3()
    Z = np.load(os.path.join(OUT, "VXV1_descriptors.npz"), allow_pickle=True)
    R = pd.read_csv(os.path.join(OUT, "VXV0_receptors.csv")); R = R[R.cdr3_complete].reset_index(drop=True)
    assert list(R.entry) == [str(e) for e in Z["entry"]]
    n = len(R); iu = np.triu_indices(n, 1)
    PC = [os.path.join(BEN, "fixed_%s_crystal.pdb" % e) for e in R.entry]; PM = [os.path.join(BEN, m) for m in R.model]
    SH = [shared_keys(a, b) for a, b in zip(PC, PM)]; MISS = [frozenset(missing_positions(a, b)) for a, b in zip(PC, PM)]
    pil, piu, pkeep = pilot_background()
    ids = [str(c) for c in L1["clone_id"]]; idx = {c: k for k, c in enumerate(ids)}
    OUTD = {}
    for name, (frame, atomset, nch) in CFG3.items():
        boxes = boxes_for(frame); nvox = [int(np.prod(s)) for _, s in boxes.values()]
        occ_idx = np.concatenate([np.arange(sum(nvox[:k]) * nch, sum(nvox[:k]) * nch + nvox[k]) for k in range(len(nvox))])
        plc = [bench_place(p, frame, lm3) for p in PC]; plm = [bench_place(p, frame, lm3) for p in PM]
        jobs = [(PC[i], plc[i], boxes, atomset, nch, SH[i], None, False) for i in range(n)] + \
               [(PM[i], plm[i], boxes, atomset, nch, SH[i], None, False) for i in range(n)]
        rem_keys = []
        for i in range(n):
            for P in sorted({MISS[j] for j in range(n) if j != i and MISS[j]}, key=lambda s: sorted(s)):
                rem_keys.append((i, P))
                jobs += [(PC[i], plc[i], boxes, atomset, nch, SH[i], P, True), (PM[i], plm[i], boxes, atomset, nch, SH[i], P, True)]
        with Pool(8) as pool:
            res = pool.map(build_job, jobs, chunksize=8)
        C.add("%s typing complete" % name, sum(r[2] for r in res), sum(r[2] for r in res) == 0, "== 0")
        Bc = np.array([res[i][0] for i in range(n)]); Bm = np.array([res[n + i][0] for i in range(n)])
        for i in range(n):
            assert res[i][3] == res[n + i][3] == frozenset(k for k in SH[i] if 105 <= k[1] <= 117)
        mass = max(abs(float(res[k][0][occ_idx].sum()) / res[k][1] - 1) for k in range(2 * n))
        C.add("%s mass conservation, %d base grids" % (name, 2 * n), "%.4f %%" % (100 * mass), mass < 5e-3, "< 0.5 %")
        remc, remm = {}, {}
        for k, (i, P) in enumerate(rem_keys):
            rc, rm = res[2 * n + 2 * k], res[2 * n + 2 * k + 1]
            if rc[3] != rm[3]:
                raise SystemExit("mask differs between crystal and model side")
            if rc[1]:
                remc[(i, P)], remm[(i, P)] = rc[0], rm[0]
        del res

        def md(Bx, rem, i, j):
            d = Bx[i] - Bx[j]
            if (i, MISS[j]) in rem:
                ix, v = rem[(i, MISS[j])]; d[ix] -= v
            if (j, MISS[i]) in rem:
                ix, v = rem[(j, MISS[i])]; d[ix] += v
            return float(np.sqrt((d * d).sum() / H ** 3))
        Cm, Mm = np.zeros((n, n)), np.zeros((n, n))
        for i, j in zip(*iu):
            Cm[i, j] = Cm[j, i] = md(Bc, remc, i, j); Mm[i, j] = Mm[j, i] = md(Bm, remm, i, j)
        own = np.sqrt(((Bc - Bm) ** 2).sum(1) / H ** 3)
        C.info("%s pairs with a non-empty mask" % name, "%d / %d" % (sum(1 for i, j in zip(*iu) if (i, MISS[j]) in remc or (j, MISS[i]) in remc), len(iu[0])))
        pj = [(VXP("structures/%s.pdb") % c, rep_place(idx[c], frame), boxes, atomset, nch, None, None, False) for c in pil]
        F_ = nch * sum(nvox)
        G = np.lib.format.open_memmap(TMPG, mode="w+", dtype=np.float16, shape=(len(pil), F_)); na = np.zeros(len(pil), int)
        with Pool(8) as pool:
            for k, (v, a, b, _) in enumerate(pool.imap(build_grid, pj, chunksize=8)):
                G[k] = v; na[k] = a
        G.flush()
        rng = np.random.default_rng(83)
        pm_ = max(abs(float(G[k][occ_idx].astype(np.float64).sum()) / na[k] - 1) for k in rng.choice(len(pil), 200, replace=False))
        C.add("%s pilot mass conservation, 200 molecules" % name, "%.4f %%" % (100 * pm_), pm_ < 5e-3, "< 0.5 %")
        errs = []
        for a, b in rng.choice(len(pil), (12, 2), replace=False):
            g = float(((G[a].astype(np.float64) - G[b].astype(np.float64)) ** 2).sum()) / H ** 3
            A_ = structure_atoms(VXP("structures/%s.pdb") % pil[a], pj[a][1], atomset, nch)
            B_ = structure_atoms(VXP("structures/%s.pdb") % pil[b], pj[b][1], atomset, nch)
            errs.append(g / sum(vg.analytic_d2(xa, wa, xb, wb, SIGMA, channels=range(nch)) for (xa, wa), (xb, wb) in zip(A_, B_)) - 1)
        C.add("%s pilot analytic grid-free distance, 12 pairs" % name, "%.4f" % np.abs(errs).max(), np.abs(errs).max() < 0.02, "< 2 %")
        Gm_ = np.zeros((len(pil), len(pil)))
        for s in range(0, F_, 20000):
            X = np.asarray(G[:, s:min(s + 20000, F_)], dtype=np.float64); Gm_ += X @ X.T
        nr = np.diag(Gm_).copy(); D2 = np.clip((nr[:, None] + nr[None] - 2 * Gm_) / H ** 3, 0, None)
        bg2 = D2[piu][pkeep]
        OUTD[name] = dict(C=Cm, M=Mm, own=own, cut=float(np.percentile(np.sqrt(bg2), 1)), bg=float(bg2.mean()))
        del G, Gm_, D2; os.remove(TMPG)

    # endpoints against the vector primary (VXV2 statistics)
    S = Z["S"]; Cv, Mv, ownv = Z["vec_vcori_geo__C"], Z["vec_vcori_geo__M"], Z["vec_vcori_geo__own"]
    up = lambda A: A[iu]
    rows = []
    for name, D in OUTD.items():
        e1 = pspear(up(D["C"]), up(D["M"]), (up(S),)); e1v = pspear(up(Cv), up(Mv), (up(S),))
        e2 = pspear(up(D["M"]), up(D["C"]), (up(Mv), up(S))); e2m = pspear(up(Mv), up(Cv), (up(D["M"]), up(S)))
        rng = np.random.default_rng(SEED); b1, bd, b2, b2m, b4 = [], [], [], [], []
        for _ in range(NB):
            r = rng.integers(0, n, n); a_, b_ = iu; k = r[a_] != r[b_]; ia, ib = r[a_][k], r[b_][k]
            g1 = pspear(D["C"][ia, ib], D["M"][ia, ib], (S[ia, ib],)); v1 = pspear(Cv[ia, ib], Mv[ia, ib], (S[ia, ib],))
            b1.append(g1); bd.append(g1 - v1)
            b2.append(pspear(D["M"][ia, ib], D["C"][ia, ib], (Mv[ia, ib], S[ia, ib])))
            b2m.append(pspear(Mv[ia, ib], Cv[ia, ib], (D["M"][ia, ib], S[ia, ib])))
            b4.append(1 - 2 * float((D["own"][r] ** 2).mean()) / D["bg"])
        ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
        o = dict(descriptor=name, E1_partial=e1, E1_lo=ci(b1)[0], E1_hi=ci(b1)[1], E1_plain=pspear(up(D["C"]), up(D["M"])),
                 delta_vs_vector=e1 - e1v, delta_lo=ci(bd)[0], delta_hi=ci(bd)[1], E2_grid=e2, E2_lo=ci(b2)[0], E2_hi=ci(b2)[1],
                 E2_vector_mirror=e2m, E2m_lo=ci(b2m)[0], E2m_hi=ci(b2m)[1], E4_rel=1 - 2 * float((D["own"] ** 2).mean()) / D["bg"],
                 E4_lo=ci(b4)[0], E4_hi=ci(b4)[1], error_over_cut=float(np.median(D["own"]) / D["cut"]),
                 within_cut=int((D["own"] <= D["cut"]).sum()), pilot_cut=D["cut"])
        rows.append(o)
        C.info("%s E1 partial [CI] / Delta vs vector [CI]" % name, "%.3f [%.3f, %.3f] / %.3f [%.3f, %.3f]"
               % (e1, o["E1_lo"], o["E1_hi"], o["delta_vs_vector"], o["delta_lo"], o["delta_hi"]))
        C.info("%s E2 grid beyond vector [CI] / mirror [CI]" % name, "%.3f [%.3f, %.3f] / %.3f [%.3f, %.3f]"
               % (e2, o["E2_lo"], o["E2_hi"], e2m, o["E2m_lo"], o["E2m_hi"]))
        C.info("%s E4 [CI]; error over cut (within)" % name, "%.3f [%.3f, %.3f]; %.3f (%d/%d)"
               % (o["E4_rel"], o["E4_lo"], o["E4_hi"], o["error_over_cut"], o["within_cut"], n))
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "VXV3_cdr3.csv"), index=False)
    save_json(dict(status="report only (A12.10.2)", rows=rows, seed=SEED, n_boot=NB), os.path.join(OUT, "VXV3_summary.json"))
    C.write()
