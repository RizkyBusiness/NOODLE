"""VXH5 (the voxel pipeline, A12.3; A12.10.1, A12.12): Arm B grid build and distances - the grid alone (no hinge).

Configuration: F1 per-chain frame (five landmarks per chain onto the reference receptor (config), then the H15 axes), loops only (IMGT CDR1
27-38, CDR2 56-65, HV4 81-86, CDR3 105-117), sigma 2.0 A, 1.0 A voxels, channels 1-7, A2/A3 typing, the A10 per-chain
boxes frozen in out/VX5b_boxes.json (7 A margin), every receptor in landmark-file order, float16 storage.
Distances: blocked Gram in float64 (A7), as VX3. Checks as VX2c and VX3, and the pilot subset reproduces the VX5b F1-loops
pilot cut (1e-3 relative). Writes tmp/VXH5_armB.f16.npy, tmp/VXH5_D_receptors.f32.npy, tmp/VXH5_D_molecules_condensed.f64.npy,
out/VXH5_manifest.csv, out/VXH5_build.json; checks/VXH5_checks.csv. No state data.
usage: python pipeline/code/vxh5_armB.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial.distance import pdist, squareform
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
EXP_NCP = expected("n_control_pairs")
import vxgrid as vg
from vxv_common import boxes_for, rep_place, build_grid, structure_atoms, SIGMA, H
from vxv1_build import pilot_background

FRAME, ATOMS, NCH = "F1", "loops", 7
GRID = os.path.join(TMP, "VXH5_armB.f16.npy")
BOXES = boxes_for(FRAME)
NVOX = [int(np.prod(s)) for _, s in BOXES.values()]
F = NCH * sum(NVOX)
OCC = np.concatenate([np.arange(sum(NVOX[:k]) * NCH, sum(NVOX[:k]) * NCH + NVOX[k]) for k in range(len(NVOX))])


def job(i, ids):
    return (VXP("structures/%s.pdb") % ids[i], rep_place(i, FRAME), BOXES, ATOMS, NCH, None, None, False)


def build16(args):
    v, na, nb, _ = build_grid(args)
    return v.astype(np.float16), na, nb, float(v[OCC].sum())


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VXH5", stop_on_fail=True)
    C.info("status", "A12 VXH5 Arm B (grid alone, A12.12); no state data read")
    L = np.load(LM1, allow_pickle=True); ids = [str(c) for c in L["clone_id"]]; n = len(ids)
    FR = np.load(os.path.join(OUT, "VX1_frames.npz"), allow_pickle=True); ismol = FR["is_molecule"]
    C.info("configuration", "F1 per chain, loops, sigma %.1f, voxel %.1f, 7 channels, boxes %s, %d features"
           % (SIGMA, H, {k: s for k, (_, s) in BOXES.items()}, F))
    G = np.lib.format.open_memmap(GRID, mode="w+", dtype=np.float16, shape=(n, F))
    na, nb, mass = np.zeros(n, int), 0, np.zeros(n)
    with Pool(8) as pool:
        for i, (v, a, b, m) in enumerate(pool.imap(build16, [job(i, ids) for i in range(n)], chunksize=8)):
            G[i] = v; na[i] = a; nb += b; mass[i] = m
            if i % 2000 == 0:
                print("built %d / %d (%.1f min)" % (i, n, (time.time() - t0) / 60), flush=True)
    G.flush()
    C.info("build", "%.1f min, %.1f GB" % ((time.time() - t0) / 60, os.path.getsize(GRID) / 1e9))
    C.add("typing complete ({:,} receptors)".format(EXP_NS), nb, nb == 0, "== 0")
    rel_all = np.abs(mass / na - 1)
    C.info("mass loss over all receptors (float64): max / n > 0.5 %", "%.4f %% / %d" % (100 * rel_all.max(), int((rel_all > 5e-3).sum())))
    rng = np.random.default_rng(101); mols = np.where(ismol)[0]
    rows = rng.choice(mols, 200, replace=False)
    rel = max(abs(float(G[r][OCC].astype(np.float64).sum()) / na[r] - 1) for r in rows)
    C.add("mass conservation, 200 molecules (stored)", "%.4f %%" % (100 * rel), rel < 5e-3, "< 0.5 %")
    errs = []
    for a, b in rng.choice(mols, (12, 2), replace=False):
        g = float(((G[a].astype(np.float64) - G[b].astype(np.float64)) ** 2).sum()) / H ** 3
        A_ = structure_atoms(job(a, ids)[0], job(a, ids)[1], ATOMS, NCH); B_ = structure_atoms(job(b, ids)[0], job(b, ids)[1], ATOMS, NCH)
        errs.append(g / sum(vg.analytic_d2(xa, wa, xb, wb, SIGMA, channels=range(NCH)) for (xa, wa), (xb, wb) in zip(A_, B_)) - 1)
    C.add("analytic grid-free distance, 12 pairs", "%.4f (%+.4f..%+.4f)" % (np.abs(errs).max(), min(errs), max(errs)),
          np.abs(errs).max() < 0.02, "< 2 %")
    worst = 0.0
    for r in rng.choice(mols, 3, replace=False):
        for (X, W), (key, (lo, shp)) in zip(structure_atoms(job(r, ids)[0], job(r, ids)[1], ATOMS, NCH), BOXES.items()):
            worst = max(worst, float(np.abs(vg.build(X, W, SIGMA, H, lo, shp) - vg.build_slow(X, W, SIGMA, H, lo, shp)).max()))
    C.add("fast grid == slow loop reference, 3 molecules, both chain boxes", "%.2e" % worst, worst < 1e-5, "< 1e-5")
    (X, W), = structure_atoms(job(int(mols[0]), ids)[0], {"A": job(int(mols[0]), ids)[1]["A"]}, ATOMS, NCH)
    Wp = W[np.random.default_rng(102).permutation(len(W))]
    lo, shp = BOXES["A"]
    G0, Gp = vg.build(X, W, SIGMA, H, lo, shp), vg.build(X, Wp, SIGMA, H, lo, shp)
    dch = [float(np.abs(G0[c] - Gp[c]).max()) for c in range(1, 7)]
    C.add("permuted atom identities: channel 1 unchanged", "%.1e" % float(np.abs(G0[0] - Gp[0]).max()),
          float(np.abs(G0[0] - Gp[0]).max()) < 1e-12, "< 1e-12")
    C.add("permuted atom identities: channels 2-7 all change", " ".join("%.3f" % v for v in dch), min(dch) > 1e-4, "each > 1e-4")
    pd.DataFrame({"row": range(n), "clone_id": ids, "is_molecule": ismol, "loop_atoms": na, "mass_occupancy": mass}) \
      .to_csv(os.path.join(OUT, "VXH5_manifest.csv"), index=False)

    # ------------------------------------------------------------ distances (VX3 route, float64 blocks)
    t1 = time.time()
    Gm = np.zeros((n, n))
    for s in range(0, F, 20000):
        Xb = np.asarray(G[:, s:min(s + 20000, F)], dtype=np.float64); Gm += Xb @ Xb.T
    C.info("Gram", "%.1f min" % ((time.time() - t1) / 60))
    nr = np.diag(Gm).copy()
    D2 = (nr[:, None] + nr[None] - 2 * Gm) / H ** 3; del Gm
    C.add("self-distances zero", "%.2e" % float(np.abs(np.diag(D2)).max()), float(np.abs(np.diag(D2)).max()) < 1e-6, "< 1e-6")
    C.add("symmetric", "%.2e" % float(np.abs(D2 - D2.T).max()), float(np.abs(D2 - D2.T).max()) < 1e-6, "< 1e-6")
    off = D2[~np.eye(n, dtype=bool)]; mn = float(off.min()); del off
    C.add("no negative squared distances", "min off-diagonal D^2 %.3e" % mn, mn > -1e-6, "> -1e-6")
    np.fill_diagonal(D2, 0.0); D = np.sqrt(np.clip(D2, 0, None)); del D2; D = 0.5 * (D + D.T)
    sub = np.sort(rng.choice(mols, 200, replace=False))
    ref = pdist(np.asarray(G[sub], dtype=np.float64)) / H ** 1.5
    r_ = float(np.max(np.abs(squareform(D[np.ix_(sub, sub)], checks=False) - ref) / ref))
    C.add("blocked Gram == scipy pdist, 200-molecule subset", "%.2e" % r_, r_ < 1e-4, "< 1e-4")
    idx = {c: k for k, c in enumerate(ids)}
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    ce = max(abs(D[idx[a], idx[b]] - np.sqrt(float(((G[idx[a]].astype(np.float64) - G[idx[b]].astype(np.float64)) ** 2).sum()) / H ** 3))
             / D[idx[a], idx[b]] for a, b in zip(NP.clone_a, NP.clone_b))
    C.add("control pairs: blocked Gram vs exact float64 (%s)" % EXP_NCP, "%.2e" % ce, ce < 1e-3, "< 1e-3")
    t = rng.integers(0, n, (10000, 3)); t = t[(t[:, 0] != t[:, 1]) & (t[:, 1] != t[:, 2]) & (t[:, 0] != t[:, 2])]
    viol = float((D[t[:, 0], t[:, 2]] - (D[t[:, 0], t[:, 1]] + D[t[:, 1], t[:, 2]])).max())
    C.add("triangle inequality, %d triples: max violation" % len(t), "%.2e" % viol, viol < 1e-6, "<= 1e-6")
    pil, piu, pkeep = pilot_background()
    pr = np.array([idx[c] for c in pil])
    pcut = float(np.percentile(D[np.ix_(pr, pr)][piu][pkeep], 1))
    V5 = pd.read_csv(os.path.join(OUT, "VX5b_configs.csv"))
    ref_cut = float(V5[(V5.frame == "F1") & (V5.atoms == "loops") & (V5.sigma == 2.0)].cut.iloc[0])
    C.add("pilot subset reproduces the VX5b F1-loops pilot cut", "%.5f vs %.5f" % (pcut, ref_cut),
          abs(pcut / ref_cut - 1) < 1e-3, "< 1e-3 rel")
    np.save(os.path.join(TMP, "VXH5_D_receptors.f32.npy"), D.astype(np.float32))
    np.save(os.path.join(TMP, "VXH5_D_molecules_condensed.f64.npy"), squareform(D[np.ix_(mols, mols)], checks=False))
    save_json(dict(frame=FRAME, atoms=ATOMS, channels=NCH, sigma=SIGMA, h=H, boxes="out/VX5b_boxes.json F1_A, F1_B",
                   features=F, grid=os.path.relpath(VXP("voxel_out/tmp/VXH5_armB.f16.npy"), ROOT), receptor_matrix=os.path.relpath(VXP("voxel_out/tmp/VXH5_D_receptors.f32.npy"), ROOT),
                   molecule_condensed=os.path.relpath(VXP("voxel_out/tmp/VXH5_D_molecules_condensed.f64.npy"), ROOT), molecule_rows=mols.tolist(),
                   minutes=round((time.time() - t0) / 60, 1)), os.path.join(OUT, "VXH5_build.json"))
    C.info("total", "%.1f min" % ((time.time() - t0) / 60))
    C.write()
