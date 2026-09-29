"""Arm C helpers (the voxel pipeline, A13): build a grid dataset for a set of receptors in given boxes, and blocked-Gram distances
(float64, A7). Reuses vxv_common.build_grid (F1 frame, typing A2/A3, sigma 2.0 A, 1.0 A voxels)."""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
from vxv_common import rep_place, build_grid, SIGMA, H


def load_boxes(path=None):
    J = json.load(open(path or os.path.join(OUT, "VXC0_boxes.json")))["boxes"]
    return {k: (np.array(J["F1_%s" % k]["lo"]), tuple(J["F1_%s" % k]["shape"])) for k in "AB"}


def occ_index(boxes, nch):
    nv = [int(np.prod(s)) for _, s in boxes.values()]
    return np.concatenate([np.arange(sum(nv[:k]) * nch, sum(nv[:k]) * nch + nv[k]) for k in range(len(nv))]), nch * sum(nv)


def _b16(args):
    """args: the build_grid arguments (8, or 9 with sigma); the occupancy index is recomputed here (cheap) rather than shipped per job"""
    occ = occ_index(args[2], args[4])[0]
    v, na, nb, _ = build_grid(args)
    return v.astype(np.float16), na, nb, float(v[occ].sum())


def build_dataset(rows, ids, boxes, atomset, nch, path, positions=None, sigma=SIGMA):
    """rows: LM1 indices. positions: optional {clone_id: set of (chain, num, ins)} - build only atoms at those positions.
    sigma: A14 (default 2.0)."""
    occ, F = occ_index(boxes, nch)
    jobs = []
    for i in rows:
        pos = None if positions is None else positions[ids[i]]
        jobs.append((VXP("structures/%s.pdb") % ids[i], rep_place(i, "F1"), boxes, atomset, nch, None, pos,
                     positions is not None, sigma))
    G = np.lib.format.open_memmap(path, mode="w+", dtype=np.float16, shape=(len(rows), F))
    na, nb, mass = np.zeros(len(rows), int), 0, np.zeros(len(rows))
    with Pool(8) as pool:
        for k, (v, a, b, m) in enumerate(pool.imap(_b16, jobs, chunksize=8)):
            G[k] = v; na[k] = a; nb += b; mass[k] = m
    G.flush()
    return G, na, nb, mass, occ, F


def gram_distances(G, F):
    n = G.shape[0]; Gm = np.zeros((n, n))
    for s in range(0, F, 20000):
        X = np.asarray(G[:, s:min(s + 20000, F)], dtype=np.float64); Gm += X @ X.T
    nr = np.diag(Gm).copy()
    D2 = (nr[:, None] + nr[None] - 2 * Gm) / H ** 3; del Gm
    stats = dict(self=float(np.abs(np.diag(D2)).max()), asym=float(np.abs(D2 - D2.T).max()))
    off = D2[~np.eye(n, dtype=bool)]; stats["min_off"] = float(off.min()); del off
    np.fill_diagonal(D2, 0.0); D = np.sqrt(np.clip(D2, 0, None)); del D2
    return 0.5 * (D + D.T), stats
