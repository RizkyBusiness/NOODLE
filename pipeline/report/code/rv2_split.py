"""the voxel pipeline report, stage 2: the exact chain x channel split of within-cluster distances, Arms B, C (sigma 2.0) and D, E (sigma 1.5, A14).

The voxel distance is D^2 = sum over chain boxes (alpha, beta) and channels (1-7) of |a - b|^2 / h^3, so it splits exactly
into 14 non-negative parts. For every cluster (>= 2 molecules) of each arm the members' grids are rebuilt exactly as the arm
(float64 build, stored as float16, as VXH5 / VXC1), every within-cluster pair is split, and for each member the mean over the
other members of each part is kept (their sum is the member's mean squared distance to the rest of its cluster).
Gate RG2: for every within-cluster pair the 14 parts sum to the stored D^2 (tmp/VXH5_D_receptors / VXC1_D_receptors) to
1e-5 relative. Writes work/split_{arm}.npz, checks/RG2_split.csv (rows of the arms run replaced). Paths come from config/voxel_config.json.:
python pipeline/report/code/rv2_split.py [arms, default all].
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", "..", "code"))
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, VXP("voxel_out/code"))
from vxv_common import boxes_for, rep_place, build_grid, H
from vxc_common import load_boxes

RP = VXP("voxel_out/report/")
ARMS = {"B": dict(lab="VXH6_molecule_labels.csv.gz", col="prop_voxel_armB", atoms="loops", boxes=boxes_for("F1"), sigma=2.0,
                  D=VXP("voxel_out/tmp/VXH5_D_receptors.f32.npy")),
        "C": dict(lab="VXC2_molecule_labels.csv.gz", col="prop_voxel_armC", atoms="cdr3", boxes=load_boxes(), sigma=2.0,
                  D=VXP("voxel_out/tmp/VXC1_D_receptors.f32.npy")),
        "D": dict(lab="VXS2_D_molecule_labels.csv.gz", col="prop_voxel_armD", atoms="loops", boxes=boxes_for("F1"), sigma=1.5,     # A14
                  D=VXP("voxel_out/tmp/VXS1_D_D_receptors.f32.npy")),
        "E": dict(lab="VXS2_E_molecule_labels.csv.gz", col="prop_voxel_armE", atoms="cdr3", boxes=load_boxes(), sigma=1.5,
                  D=VXP("voxel_out/tmp/VXS1_E_D_receptors.f32.npy"))}
L1 = np.load(VXP("landmarks/out/LM1_internal_coords.npz"), allow_pickle=True); IDS = [str(c) for c in L1["clone_id"]]
IDX = {c: i for i, c in enumerate(IDS)}


def blocks(boxes):
    """(start, stop) of the 14 parts in the flattened grid: box A channels 1-7, then box B channels 1-7"""
    out, s = [], 0
    for k in "AB":
        nv = int(np.prod(boxes[k][1]))
        for c in range(7):
            out.append((s + c * nv, s + (c + 1) * nv))
        s += 7 * nv
    return out


def cluster_split(args):
    arm, mem_rows = args
    A = ARMS[arm]; BL = blocks(A["boxes"])
    G = np.array([build_grid((VXP("structures/%s.pdb") % IDS[r], rep_place(r, "F1"), A["boxes"], A["atoms"], 7, None, None, False, A["sigma"]))[0]
                  .astype(np.float16).astype(np.float64) for r in mem_rows])
    k = len(mem_rows); P = np.zeros((k, k, 14))
    for b, (s, e) in enumerate(BL):
        X = G[:, s:e]; nr = (X * X).sum(1); Gm = X @ X.T
        P[:, :, b] = np.clip(nr[:, None] + nr[None] - 2 * Gm, 0, None) / H ** 3
    for b in range(14):
        np.fill_diagonal(P[:, :, b], 0.0)
    mean = P.sum(1) / (k - 1)                                      # per member, mean over the other members, per part
    iu = np.triu_indices(k, 1)
    return mem_rows, mean, P.sum(2)[iu], iu


if __name__ == "__main__":
    rows = []; todo = sys.argv[1:] or list(ARMS)                 # e.g. "D E"
    for arm in todo:
        A = ARMS[arm]
        M = pd.read_csv(os.path.join(VXP("voxel_out/out"), A["lab"]), usecols=["clone_id", A["col"]])
        lab = M[A["col"]].values; rep = np.array([IDX[c] for c in M.clone_id])
        cl = pd.Series(lab[lab >= 0]).value_counts()
        jobs = [(arm, list(rep[lab == c])) for c in cl.index]              # largest first
        Dm = np.load(A["D"], mmap_mode="r")
        split = np.full((len(M), 14), np.nan); worst, npairs = 0.0, 0
        pos = {r: j for j, r in enumerate(rep)}
        with Pool(8) as pool:
            for mem_rows, mean, pt, iu in pool.imap_unordered(cluster_split, jobs, chunksize=4):
                for j, r in enumerate(mem_rows):
                    split[pos[r]] = mean[j]
                ref = np.array([float(Dm[mem_rows[a], mem_rows[b]]) ** 2 for a, b in zip(*iu)])
                worst = max(worst, float(np.max(np.abs(pt - ref) / np.maximum(ref, 1e-12)))); npairs += len(ref)
        np.savez_compressed(RP + "work/split_%s.npz" % arm, clone_id=np.array(M.clone_id), split=split,
                            parts=np.array(["%s:%d" % (k, c + 1) for k in "AB" for c in range(7)]))
        ok = worst < 1e-5
        rows.append(dict(arm=arm, clusters=len(jobs), clustered=int((lab >= 0).sum()), pairs=npairs, max_rel_err=worst, passed=ok))
        print("Arm %s: %d clusters, %d pairs, parts sum == stored D^2: max rel err %.2e -> %s" % (arm, len(jobs), npairs, worst, "PASS" if ok else "FAIL"), flush=True)
    old = pd.read_csv(RP + "checks/RG2_split.csv") if os.path.exists(RP + "checks/RG2_split.csv") else pd.DataFrame(columns=["arm"])
    pd.concat([old[~old.arm.isin(todo)], pd.DataFrame(rows)]).sort_values("arm").to_csv(RP + "checks/RG2_split.csv", index=False)
    if not all(r["passed"] for r in rows):
        raise SystemExit("RG2 FAILED")
