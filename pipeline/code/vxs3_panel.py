"""VXS3 (the voxel pipeline, A14.2.3): validity panel at sigma 1.5 A beside sigma 2.0 A - report only, cannot change any rule.

Construction exactly as VXV3 (vxv3_cdr3.py, itself the VXV1 construction): the 55 crystal-benchmarked receptors, F1 frame,
VX5b F1 boxes, 1.0 A voxels, shared atoms, A12.9.2 pair masks; the pilot cut and mean squared background distance from the
2,000-molecule pilot. Atom sets loops (CDR1, CDR2, HV4, CDR3) and CDR3 (105-117), occupancy and channels 1-7, each at
sigma 2.0 and 1.5 (8 grids). Checks (implementation; stop on failure): the sigma 2.0 loops grids reproduce VXV1's stored
crystal, model and own-distance matrices and pilot cut (1e-9 relative); the sigma 2.0 CDR3 grids reproduce VXV3's point
estimates (1e-9). Endpoints (VXV2 statistics): E1 partial Spearman beyond sequence (S), E2 beyond the vector primary,
E4 reliability, error over cut; paired differences sigma 1.5 - sigma 2.0 of E1 and E4, and Delta vs the vector primary,
with 2,000 bootstrap resamples over receptors (seed 0; the same resamples for every descriptor). Nothing here is a gate
on the arms (A14.3). Writes out/VXS3_panel.csv, out/VXS3_sigma_diff.csv, out/VXS3_summary.json; checks/VXS3_checks.csv.
usage: python pipeline/code/vxs3_panel.py
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

CFG = [(a, nch, sig) for a in ("loops", "cdr3") for nch in (1, 7) for sig in (2.0, 1.5)]
NAME = lambda a, nch, sig: "grid_F1%s_%s_s%s" % ("L" if a == "loops" else "C3", "occ" if nch == 1 else "7ch", str(sig).replace(".", ""))
TMPG = os.path.join(TMP, "VXS3_pilot.f16.npy")
NB, SEED = 2000, 0

if __name__ == "__main__":
    C = Checks("VXS3", stop_on_fail=True)
    C.info("status", "A14 VXS3 validity panel, sigma 1.5 beside 2.0; report only; no state data read")
    lm3 = load_lm3()
    Z = np.load(os.path.join(OUT, "VXV1_descriptors.npz"), allow_pickle=True)
    R = pd.read_csv(os.path.join(OUT, "VXV0_receptors.csv")); R = R[R.cdr3_complete].reset_index(drop=True)
    assert list(R.entry) == [str(e) for e in Z["entry"]]
    n = len(R); iu = np.triu_indices(n, 1)
    PC = [os.path.join(BEN, "fixed_%s_crystal.pdb" % e) for e in R.entry]; PM = [os.path.join(BEN, m) for m in R.model]
    SH = [shared_keys(a, b) for a, b in zip(PC, PM)]; MISS = [frozenset(missing_positions(a, b)) for a, b in zip(PC, PM)]
    pil, piu, pkeep = pilot_background()
    ids = [str(c) for c in L1["clone_id"]]; idx = {c: k for k, c in enumerate(ids)}
    boxes = boxes_for("F1"); OUTD = {}
    plc = [bench_place(p, "F1", lm3) for p in PC]; plm = [bench_place(p, "F1", lm3) for p in PM]
    for atomset, nch, sig in CFG:
        name = NAME(atomset, nch, sig)
        nvox = [int(np.prod(s)) for _, s in boxes.values()]
        occ_idx = np.concatenate([np.arange(sum(nvox[:k]) * nch, sum(nvox[:k]) * nch + nvox[k]) for k in range(len(nvox))])
        jobs = [(PC[i], plc[i], boxes, atomset, nch, SH[i], None, False, sig) for i in range(n)] + \
               [(PM[i], plm[i], boxes, atomset, nch, SH[i], None, False, sig) for i in range(n)]
        rem_keys = []
        for i in range(n):
            for P in sorted({MISS[j] for j in range(n) if j != i and MISS[j]}, key=lambda s: sorted(s)):
                rem_keys.append((i, P))
                jobs += [(PC[i], plc[i], boxes, atomset, nch, SH[i], P, True, sig), (PM[i], plm[i], boxes, atomset, nch, SH[i], P, True, sig)]
        with Pool(8) as pool:
            res = pool.map(build_job, jobs, chunksize=8)
        C.add("%s typing complete" % name, sum(r[2] for r in res), sum(r[2] for r in res) == 0, "== 0")
        Bc = np.array([res[i][0] for i in range(n)]); Bm = np.array([res[n + i][0] for i in range(n)])
        for i in range(n):
            assert res[i][3] == res[n + i][3]
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
        own = np.sqrt(((Bc - Bm) ** 2).sum(1) / H ** 3); del Bc, Bm, remc, remm
        pj = [(VXP("structures/%s.pdb") % c, rep_place(idx[c], "F1"), boxes, atomset, nch, None, None, False, sig) for c in pil]
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
            errs.append(g / sum(vg.analytic_d2(xa, wa, xb, wb, sig, channels=range(nch)) for (xa, wa), (xb, wb) in zip(A_, B_)) - 1)
        C.add("%s pilot analytic grid-free distance, 12 pairs" % name, "%.4f" % np.abs(errs).max(), np.abs(errs).max() < 0.02, "< 2 %")
        Gm_ = np.zeros((len(pil), len(pil)))
        for s in range(0, F_, 20000):
            X = np.asarray(G[:, s:min(s + 20000, F_)], dtype=np.float64); Gm_ += X @ X.T
        nr = np.diag(Gm_).copy(); D2 = np.clip((nr[:, None] + nr[None] - 2 * Gm_) / H ** 3, 0, None)
        bg2 = D2[piu][pkeep]
        OUTD[name] = dict(C=Cm, M=Mm, own=own, cut=float(np.percentile(np.sqrt(bg2), 1)), bg=float(bg2.mean()), atoms=atomset, nch=nch, sigma=sig)
        del G, Gm_, D2; os.remove(TMPG)
        print("%s done" % name, flush=True)

    # reproduction checks at sigma 2.0
    for nch, k in ((1, "grid_F1L_occ"), (7, "grid_F1L_7ch")):
        mine = OUTD[NAME("loops", nch, 2.0)]
        e = max(float(np.max(np.abs(mine[f][iu] - Z[k + "__" + g][iu]) / Z[k + "__" + g][iu])) for f, g in (("C", "C"), ("M", "M")))
        e = max(e, float(np.max(np.abs(mine["own"] - Z[k + "__own"]) / Z[k + "__own"])), abs(mine["cut"] / float(Z[k + "__cut"]) - 1))
        C.add("sigma 2.0 loops %s reproduces VXV1 (C, M, own, pilot cut)" % ("occ" if nch == 1 else "7ch"), "%.2e rel" % e, e < 1e-9, "< 1e-9")
    S = Z["S"]; Cv, Mv, ownv = Z["vec_vcori_geo__C"], Z["vec_vcori_geo__M"], Z["vec_vcori_geo__own"]
    up = lambda A: A[iu]
    V3 = pd.read_csv(os.path.join(OUT, "VXV3_cdr3.csv")).set_index("descriptor")
    for nch, k in ((1, "grid_F1C3_occ"), (7, "grid_F1C3_7ch")):
        mine = OUTD[NAME("cdr3", nch, 2.0)]
        pt = dict(E1_partial=pspear(up(mine["C"]), up(mine["M"]), (up(S),)), E1_plain=pspear(up(mine["C"]), up(mine["M"])),
                  E4_rel=1 - 2 * float((mine["own"] ** 2).mean()) / mine["bg"], error_over_cut=float(np.median(mine["own"]) / mine["cut"]),
                  pilot_cut=mine["cut"])
        e = max(abs(pt[f] / float(V3.loc[k, f]) - 1) for f in pt)
        C.add("sigma 2.0 CDR3 %s reproduces VXV3 point estimates" % ("occ" if nch == 1 else "7ch"), "%.2e rel" % e, e < 1e-9, "< 1e-9")

    # endpoints and paired differences, same bootstrap resamples for every descriptor
    names = list(OUTD); rows = []
    rng = np.random.default_rng(SEED); draws = [rng.integers(0, n, n) for _ in range(NB)]
    def stats_(Dd, ia, ib, r=None):
        own = Dd["own"] if r is None else Dd["own"][r]
        return (pspear(Dd["C"][ia, ib], Dd["M"][ia, ib], (S[ia, ib],)), 1 - 2 * float((own ** 2).mean()) / Dd["bg"])
    boot = {k: [] for k in names}; bvec = []
    for r in draws:
        a_, b_ = iu; kk = r[a_] != r[b_]; ia, ib = r[a_][kk], r[b_][kk]
        for k in names:
            boot[k].append(stats_(OUTD[k], ia, ib, r))
        bvec.append(pspear(Cv[ia, ib], Mv[ia, ib], (S[ia, ib],)))
    boot = {k: np.array(v) for k, v in boot.items()}; bvec = np.array(bvec)
    ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    e1v = pspear(up(Cv), up(Mv), (up(S),))
    for k in names:
        Dd = OUTD[k]; e1, e4 = stats_(Dd, iu[0], iu[1])
        e2 = pspear(up(Dd["M"]), up(Dd["C"]), (up(Mv), up(S)))
        rows.append(dict(descriptor=k, atoms=Dd["atoms"], channels=Dd["nch"], sigma=Dd["sigma"], E1_partial=e1, E1_lo=ci(boot[k][:, 0])[0], E1_hi=ci(boot[k][:, 0])[1],
                         E1_plain=pspear(up(Dd["C"]), up(Dd["M"])), delta_vs_vector=e1 - e1v, delta_lo=ci(boot[k][:, 0] - bvec)[0], delta_hi=ci(boot[k][:, 0] - bvec)[1],
                         E2_grid=e2, E4_rel=e4, E4_lo=ci(boot[k][:, 1])[0], E4_hi=ci(boot[k][:, 1])[1],
                         error_over_cut=float(np.median(Dd["own"]) / Dd["cut"]), within_cut=int((Dd["own"] <= Dd["cut"]).sum()), pilot_cut=Dd["cut"]))
        C.info("%s: E1 [CI] / Delta vs vector [CI] / E4 [CI] / error over cut" % k, "%.3f [%.3f, %.3f] / %+.3f [%.3f, %.3f] / %.3f [%.3f, %.3f] / %.3f"
               % (e1, rows[-1]["E1_lo"], rows[-1]["E1_hi"], rows[-1]["delta_vs_vector"], rows[-1]["delta_lo"], rows[-1]["delta_hi"],
                  e4, rows[-1]["E4_lo"], rows[-1]["E4_hi"], rows[-1]["error_over_cut"]))
    P = pd.DataFrame(rows).set_index("descriptor"); diff = []
    for atomset in ("loops", "cdr3"):
        for nch in (1, 7):
            a, b = NAME(atomset, nch, 1.5), NAME(atomset, nch, 2.0)
            dE1, dE4 = boot[a][:, 0] - boot[b][:, 0], boot[a][:, 1] - boot[b][:, 1]
            o = dict(atoms=atomset, channels=nch, dE1=P.loc[a, "E1_partial"] - P.loc[b, "E1_partial"], dE1_lo=ci(dE1)[0], dE1_hi=ci(dE1)[1],
                     dE4=P.loc[a, "E4_rel"] - P.loc[b, "E4_rel"], dE4_lo=ci(dE4)[0], dE4_hi=ci(dE4)[1],
                     error_over_cut_s15=P.loc[a, "error_over_cut"], error_over_cut_s20=P.loc[b, "error_over_cut"])
            diff.append(o)
            C.info("sigma 1.5 - 2.0, %s %s: dE1 [CI] / dE4 [CI] / error over cut 1.5 vs 2.0" % (atomset, "occ" if nch == 1 else "7ch"),
                   "%+.3f [%.3f, %.3f] / %+.3f [%.3f, %.3f] / %.3f vs %.3f" % (o["dE1"], o["dE1_lo"], o["dE1_hi"], o["dE4"], o["dE4_lo"], o["dE4_hi"],
                                                                         o["error_over_cut_s15"], o["error_over_cut_s20"]))
    P.reset_index().to_csv(os.path.join(OUT, "VXS3_panel.csv"), index=False)
    pd.DataFrame(diff).to_csv(os.path.join(OUT, "VXS3_sigma_diff.csv"), index=False)
    save_json(dict(status="report only (A14.3)", n_boot=NB, seed=SEED, vector_E1=e1v, panel=rows, sigma_diff=diff), os.path.join(OUT, "VXS3_summary.json"))
    C.write()
