"""VX3 (the voxel pipeline): pairwise grid distances by the blocked-Gram route (DESIGN VX3).

D^2(a, b) = (|a|^2 + |b|^2 - 2 a.b) / h^3 over channels 1-7 of the VX2c grids, the Gram matrix accumulated over
feature blocks streamed from the memmap, each block converted to float64 for the matmul (A7; the design said float32).
All N receptors (A6.1). Writes <voxel out>/tmp/VX3_D_receptors.f32.npy (N x N, landmark-file order) and
<voxel out>/tmp/VX3_D_molecules_condensed.f64.npy (the M molecules, scipy condensed order, for linkage).
usage: python pipeline/code/vx3_distances.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, time, json
import numpy as np, pandas as pd
from scipy.spatial.distance import pdist, squareform
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NCP = expected("n_control_pairs")
import vxgrid as vg

B = json.load(open(os.path.join(OUT, "VX2c_build.json")))
H = B["h"]
G = np.load(os.path.join(ROOT, B["grid"]), mmap_mode="r")   # "grid" is stored relative to project_root
n = G.shape[0]
NV = int(np.prod(B["shape"])); NF = vg.NPRIMARY * NV
FLAT = G.reshape(n, -1)
BLOCK = 20000

if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VX3")
    Gm = np.zeros((n, n))
    for k, s in enumerate(range(0, NF, BLOCK)):
        X = np.asarray(FLAT[:, s:min(s + BLOCK, NF)], dtype=np.float64)
        Gm += X @ X.T
        if k % 10 == 0:
            print("block %d / %d (%.1f min)" % (k, -(-NF // BLOCK), (time.time() - t0) / 60), flush=True)
    C.info("Gram time", "%.1f min, %d feature blocks of %d" % ((time.time() - t0) / 60, -(-NF // BLOCK), BLOCK))
    nr = np.diag(Gm).copy()
    D2 = (nr[:, None] + nr[None] - 2 * Gm) / H ** 3
    del Gm
    self_ = float(np.abs(np.diag(D2)).max())
    C.add("self-distances zero: max |D^2_ii|", "%.2e (median |a|^2 %.1f)" % (self_, np.median(nr)), self_ < 1e-6 * np.median(nr) or self_ < 1e-6,
          "< 1e-6 (abs, or relative to median |a|^2)")
    asym = float(np.abs(D2 - D2.T).max())
    C.add("symmetric", "%.2e" % asym, asym < 1e-6 or asym < 1e-6 * np.median(nr), "< 1e-6")
    off = D2[~np.eye(n, dtype=bool)]
    mn = float(off.min())
    C.add("no negative squared distances beyond rounding", "min off-diagonal D^2 %.3e" % mn, mn > -1e-6, "> -1e-6")
    C.info("off-diagonal D^2 range", "%.3f .. %.3f" % (mn, float(off.max())))
    del off
    np.fill_diagonal(D2, 0.0)
    D = np.sqrt(np.clip(D2, 0, None)); del D2
    D = 0.5 * (D + D.T)

    # 200-molecule subset: blocked Gram vs scipy pdist on dense float64 grids
    rng = np.random.default_rng(41)
    ismol = pd.read_csv(os.path.join(OUT, "VX2c_manifest.csv")).is_molecule.values
    sub = np.sort(rng.choice(np.where(ismol)[0], 200, replace=False))
    dense = np.asarray(FLAT[sub, :NF], dtype=np.float64)
    ref = pdist(dense) / H ** 1.5
    del dense
    got = squareform(D[np.ix_(sub, sub)], checks=False)
    rel = float(np.max(np.abs(got - ref) / ref))
    C.add("blocked Gram == scipy pdist, 200-molecule subset (max rel)", "%.2e" % rel, rel < 1e-4, "< 1e-4")

    # extra: the near-identical control pairs, exact float64 recomputation (the smallest distances)
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; idx = {c: i for i, c in enumerate(ids)}
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    ce = []
    for a, b in zip(NP.clone_a, NP.clone_b):
        i, j = idx[a], idx[b]
        ex = np.sqrt(vg.grid_d2(G[i], G[j], H))
        ce.append(abs(D[i, j] - ex) / ex)
    C.add("control pairs: blocked Gram vs exact float64 (%s pairs, max rel)" % EXP_NCP, "%.2e" % max(ce), max(ce) < 1e-3,
          "< 1e-3 (extra check)")

    # triangle inequality on 10,000 random triples
    t = rng.integers(0, n, (10000, 3))
    t = t[(t[:, 0] != t[:, 1]) & (t[:, 1] != t[:, 2]) & (t[:, 0] != t[:, 2])]
    viol = D[t[:, 0], t[:, 2]] - (D[t[:, 0], t[:, 1]] + D[t[:, 1], t[:, 2]])
    C.add("triangle inequality, %d random triples: max violation" % len(t), "%.2e" % viol.max(), viol.max() < 1e-6,
          "<= 1e-6")

    np.save(os.path.join(TMP, "VX3_D_receptors.f32.npy"), D.astype(np.float32))
    mol = np.where(ismol)[0]
    Dm = squareform(D[np.ix_(mol, mol)], checks=False).astype(np.float64)
    np.save(os.path.join(TMP, "VX3_D_molecules_condensed.f64.npy"), Dm)
    C.info("molecule condensed distances", "%d pairs, median %.4f" % (len(Dm), float(np.median(Dm))))
    save_json(dict(receptor_matrix=os.path.relpath(VXP("voxel_out/tmp/VX3_D_receptors.f32.npy"), ROOT), molecule_condensed=
                   os.path.relpath(VXP("voxel_out/tmp/VX3_D_molecules_condensed.f64.npy"), ROOT), molecule_rows=mol.tolist(),
                   channels="1-7", minutes=round((time.time() - t0) / 60, 1)), os.path.join(OUT, "VX3_distances.json"))
    C.write()
