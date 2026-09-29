"""VX2c (the voxel pipeline): full build at the production settings (DESIGN VX2c + amendments A5.1, A6.1).

All N receptors in landmark-file order (row i = the landmark file index i; the M molecules are flagged), sigma 2.0 A on the
1.0 A voxel, frozen VX1 box, 8 channels (channels 1-7 primary), float16 memmap
<voxel out>/tmp/VX2c_h1.0_s2.00.f16.npy, shape (N, 8, 58, 50, 50). Manifest <voxel out>/out/VX2c_manifest.csv.
Same frame / typing / smearing code as the pilot (vx2a_pilot.structure, vxgrid.build).
usage: python pipeline/code/vx2c_build.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
import vxgrid as vg
import vx2a_pilot as PA

SIGMA, H = 2.0, 1.0
GRID = os.path.join(TMP, "VX2c_h1.0_s2.00.f16.npy")


def job(i):
    return ("rep", VXP("structures/%s.pdb") % PA.IDS[i], PA.FR["R_box"][i], PA.FR["t_box"][i], None)


def build_one(i):
    X, W, meta, bad = PA.structure(job(i))
    G = vg.build(X, W, SIGMA, H, PA.LO, PA.SHAPE)
    out = int((~((X >= PA.LO) & (X < PA.LO + np.array(PA.SHAPE) * H)).all(1)).sum())
    return G.astype(np.float16), len(X), G.reshape(G.shape[0], -1).sum(1), out, len(bad)


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VX2c")
    n = len(PA.IDS)
    ismol = PA.FR["is_molecule"]
    C.info("sigma / voxel / box / rows", "%.1f A / %.1f A / %s / %d (molecules %d)" % (SIGMA, H, PA.SHAPE, n, ismol.sum()))
    G = np.lib.format.open_memmap(GRID, mode="w+", dtype=np.float16, shape=(n, vg.NCH) + PA.SHAPE)
    man = []
    with Pool(8) as pool:
        for i, (g, na, ms, out, nb) in enumerate(pool.imap(build_one, range(n), chunksize=8)):
            G[i] = g
            man.append(dict(row=i, clone_id=PA.IDS[i], is_molecule=bool(ismol[i]), heavy_atoms=na, outside=out,
                            unmatched=nb, **{"mass_%s" % c: float(m) for c, m in zip(vg.CHANNELS, ms)}))
            if i % 1000 == 0:
                print("built %d / %d  (%.0f s)" % (i, n, time.time() - t0), flush=True)
    G.flush()
    M = pd.DataFrame(man)
    M.to_csv(os.path.join(OUT, "VX2c_manifest.csv"), index=False)
    C.info("build time", "%.1f min" % ((time.time() - t0) / 60))
    C.info("grid file", "%s (%.1f GB)" % (os.path.relpath(GRID, ROOT), os.path.getsize(GRID) / 1e9))
    C.add("atom typing complete: unmatched atoms, all {:,}".format(EXP_NS), int(M.unmatched.sum()), M.unmatched.sum() == 0, "== 0")
    for c, nm in enumerate(vg.CHANNELS):
        C.info("total mass channel %d %s" % (c + 1, nm), "%.1f" % M["mass_%s" % nm].sum())
    C.info("heavy atoms outside box (all rows)", "%d / %d" % (M.outside.sum(), M.heavy_atoms.sum()))

    # mass conservation on stored grids
    rng = np.random.default_rng(31)
    mols = np.where(ismol)[0]
    rows = rng.choice(mols, 200, replace=False)
    rel = np.array([abs(float(G[r, 0].astype(np.float64).sum()) / M.heavy_atoms[r] - 1) for r in rows])
    C.add("mass conservation, 200 random molecules: max", "%.4f %%" % (100 * rel.max()), rel.max() < 5e-3, "< 0.5 %")
    relall = np.abs(M.mass_occupancy / M.heavy_atoms - 1)
    C.info("mass loss over all {:,} (float64, before storage): max / n > 0.5 %".format(EXP_NS),
           "%.4f %% / %d" % (100 * relall.max(), int((relall > 5e-3).sum())))

    # analytic grid-free reference, 12 pairs
    pr = rng.choice(mols, (12, 2), replace=False)
    errs = []
    for a, b in pr:
        XA, WA, _, _ = PA.structure(job(a)); XB, WB, _, _ = PA.structure(job(b))
        errs.append(vg.grid_d2(G[a], G[b], H) / vg.analytic_d2(XA, WA, XB, WB, SIGMA) - 1)
    errs = np.array(errs)
    C.add("analytic grid-free distance, 12 pairs: max |grid/exact - 1|",
          "%.4f (%+.4f..%+.4f)" % (np.abs(errs).max(), errs.min(), errs.max()), np.abs(errs).max() < 0.02, "< 2 %")

    # independent slow reimplementation, 3 molecules
    worst = 0.0
    for r in rng.choice(mols, 3, replace=False):
        X, W, _, _ = PA.structure(job(r))
        worst = max(worst, float(np.abs(vg.build(X, W, SIGMA, H, PA.LO, PA.SHAPE)
                                        - vg.build_slow(X, W, SIGMA, H, PA.LO, PA.SHAPE)).max()))
    C.add("fast grid == slow loop reference, 3 molecules (max abs)", "%.2e" % worst, worst < 1e-5, "< 1e-5")

    # permuted-atom control
    X, W, _, _ = PA.structure(job(int(mols[0])))
    Wp = W[np.random.default_rng(32).permutation(len(W))]
    G0, Gp = vg.build(X, W, SIGMA, H, PA.LO, PA.SHAPE), vg.build(X, Wp, SIGMA, H, PA.LO, PA.SHAPE)
    d1 = float(np.abs(G0[0] - Gp[0]).max()); dch = [float(np.abs(G0[c] - Gp[c]).max()) for c in range(1, 7)]
    C.add("permuted atom identities: channel 1 unchanged", "%.1e" % d1, d1 < 1e-12, "< 1e-12")
    C.add("permuted atom identities: channels 2-7 all change", " ".join("%.3f" % v for v in dch), min(dch) > 1e-4,
          "each > 1e-4")

    # the production grid is the one the pilot measured: VX2b sigma-2.0 pilot rows are bitwise identical
    MP = pd.read_csv(os.path.join(OUT, "VX2a_manifest.csv"))
    Gs = np.load(os.path.join(TMP, "VX2b_h1.0_s2.00.f16.npy"), mmap_mode="r")
    prow = MP[MP.kind == "rep"].sample(5, random_state=33)
    same = all(np.array_equal(Gs[int(r.row)], G[int(r.lm1_index)]) for _, r in prow.iterrows())
    C.add("full-build rows == VX2b sigma 2.0 pilot rows (5 molecules, bitwise)", same, same, "identical")
    save_json(dict(sigma=SIGMA, h=H, shape=PA.SHAPE, lo=PA.LO.tolist(), truncation_sigma=3.0, rows=n,
                   row_order="LM1 clone_id order", grid=os.path.relpath(GRID, ROOT), channels=vg.CHANNELS,
                   primary_channels=vg.CHANNELS[:vg.NPRIMARY], value="mass per voxel; L2 inner product = sum/h^3"),
              os.path.join(OUT, "VX2c_build.json"))
    C.write()
