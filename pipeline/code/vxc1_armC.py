"""VXC1 (the voxel pipeline, A13.2 / A13.3): Arm C build and distances - CDR3 only (IMGT 105-117, both chains), F1 per-chain frame,
channels 1-7, sigma 2.0 A, 1.0 A voxels, the frozen VXC0 boxes, every receptor, float16; blocked Gram in float64.
Checks as VXH5 plus the box-cropping check (A13.8.4): for the 55 VXV receptors, crystal-crystal and model-model matrices
in the Arm C boxes equal those in the larger VX5b F1 boxes (1e-6 relative), and the pilot subset reproduces the VXV3
grid_F1C3_7ch pilot cut (1e-3 relative). Writes tmp/VXC1_*; out/VXC1_build.json; checks/VXC1_checks.csv. No state data.
usage: python pipeline/code/vxc1_armC.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.spatial.distance import pdist, squareform
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NCP = expected("n_control_pairs")
import vxgrid as vg
from vxv_common import boxes_for, rep_place, bench_place, build_grid, structure_atoms, shared_keys, SIGMA, H, BEN
from vxv1_build import pilot_background
from vxc_common import load_boxes, build_dataset, gram_distances

ATOMS, NCH = "cdr3", 7
GRID = os.path.join(TMP, "VXC1_armC.f16.npy")

if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VXC1", stop_on_fail=True)
    C.info("status", "A13 VXC1 Arm C (CDR3 only, 7 channels); no state data read")
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; n = len(ids); idx = {c: k for k, c in enumerate(ids)}
    FR = np.load(os.path.join(OUT, "VX1_frames.npz"), allow_pickle=True); ismol = FR["is_molecule"]; mols = np.where(ismol)[0]
    BX = load_boxes()
    G, na, nb, mass, occ, F = build_dataset(range(n), ids, BX, ATOMS, NCH, GRID)
    C.info("build", "%.1f min, %.1f GB, %d features" % ((time.time() - t0) / 60, os.path.getsize(GRID) / 1e9, F))
    C.add("typing complete", nb, nb == 0, "== 0")
    rel_all = np.abs(mass / na - 1)
    C.info("mass loss over all receptors (float64): max / n > 0.5 %", "%.4f %% / %d" % (100 * rel_all.max(), int((rel_all > 5e-3).sum())))
    rng = np.random.default_rng(111)
    rows = rng.choice(mols, 200, replace=False)
    rel = max(abs(float(G[r][occ].astype(np.float64).sum()) / na[r] - 1) for r in rows)
    C.add("mass conservation, 200 molecules (stored)", "%.4f %%" % (100 * rel), rel < 5e-3, "< 0.5 %")
    job = lambda i: (VXP("structures/%s.pdb") % ids[i], rep_place(i, "F1"))
    errs = []
    for a, b in rng.choice(mols, (12, 2), replace=False):
        g = float(((G[a].astype(np.float64) - G[b].astype(np.float64)) ** 2).sum()) / H ** 3
        A_ = structure_atoms(*job(a), ATOMS, NCH); B_ = structure_atoms(*job(b), ATOMS, NCH)
        errs.append(g / sum(vg.analytic_d2(xa, wa, xb, wb, SIGMA, channels=range(NCH)) for (xa, wa), (xb, wb) in zip(A_, B_)) - 1)
    C.add("analytic grid-free distance, 12 pairs", "%.4f (%+.4f..%+.4f)" % (np.abs(errs).max(), min(errs), max(errs)),
          np.abs(errs).max() < 0.02, "< 2 %")
    worst = 0.0
    for r in rng.choice(mols, 3, replace=False):
        for (X, W), (key, (lo, shp)) in zip(structure_atoms(*job(r), ATOMS, NCH), BX.items()):
            worst = max(worst, float(np.abs(vg.build(X, W, SIGMA, H, lo, shp) - vg.build_slow(X, W, SIGMA, H, lo, shp)).max()))
    C.add("fast grid == slow loop reference, 3 molecules, both boxes", "%.2e" % worst, worst < 1e-5, "< 1e-5")
    # A13.9: the first molecule whose CDR3 atoms carry non-zero mass in every channel 2-7; both chains, both boxes
    for pm_ in mols:
        parts = structure_atoms(*job(int(pm_)), ATOMS, NCH)
        if (np.vstack([w for _, w in parts])[:, 1:7].sum(0) > 0).all():
            break
    C.info("permuted-atom control molecule (A13.9)", ids[int(pm_)])
    prng = np.random.default_rng(112)
    G0 = np.concatenate([vg.build(X, W, SIGMA, H, lo, shp).reshape(NCH, -1) for (X, W), (lo, shp) in zip(parts, BX.values())], 1)
    Gp = np.concatenate([vg.build(X, W[prng.permutation(len(W))], SIGMA, H, lo, shp).reshape(NCH, -1)
                         for (X, W), (lo, shp) in zip(parts, BX.values())], 1)
    dch = [float(np.abs(G0[c] - Gp[c]).max()) for c in range(1, 7)]
    C.add("permuted atom identities: channel 1 unchanged", "%.1e" % float(np.abs(G0[0] - Gp[0]).max()),
          float(np.abs(G0[0] - Gp[0]).max()) < 1e-12, "< 1e-12")
    C.add("permuted atom identities: channels 2-7 all change", " ".join("%.3f" % v for v in dch), min(dch) > 1e-4, "each > 1e-4")

    # box-cropping check (A13.8.4): 55 VXV receptors, Arm C boxes vs VX5b F1 boxes (float64)
    lm3 = load_lm3()
    R = pd.read_csv(os.path.join(OUT, "VXV0_receptors.csv")); R = R[R.cdr3_complete].reset_index(drop=True)
    PC = [os.path.join(BEN, "fixed_%s_crystal.pdb" % e) for e in R.entry]; PM = [os.path.join(BEN, m) for m in R.model]
    SH = [shared_keys(a, b) for a, b in zip(PC, PM)]
    # A13.11: reference = the Arm C boxes expanded by 10 A on every face (strictly larger, no truncation)
    BIG = {k: (lo - 10.0, tuple(int(x) + 20 for x in shp)) for k, (lo, shp) in BX.items()}; worst_c = 0.0
    from vxc_common import occ_index
    for boxes_pair in [(BX, BIG)]:
        mats = []
        for nm_, bx in zip(("Arm C boxes", "expanded reference"), boxes_pair):
            jobs = [(p, bench_place(p, "F1", lm3), bx, ATOMS, NCH, sh, None, False) for p, sh in zip(PC + PM, SH + SH)]
            with Pool(8) as pool:
                res = pool.map(build_grid, jobs)
            oc = occ_index(bx, NCH)[0]
            ml = max(1 - float(r[0][oc].sum()) / r[1] for r in res)
            C.add("box cropping: benchmark mass loss in %s (%d grids)" % (nm_, len(PC) + len(PM)), "%.1e" % ml, ml <= 1e-9, "<= 1e-9")
            V = np.array([r[0] for r in res]); nr_ = len(PC)
            Cm = squareform(pdist(V[:nr_])) / H ** 1.5; Mm = squareform(pdist(V[nr_:])) / H ** 1.5
            mats.append((Cm, Mm))
        iu = np.triu_indices(len(PC), 1)
        for a_, b_ in zip(mats[0], mats[1]):
            worst_c = max(worst_c, float(np.max(np.abs(a_[iu] - b_[iu]) / b_[iu])))
    C.add("box cropping: %d VXV receptors, Arm C boxes == expanded boxes (C and M)" % len(PC), "%.2e rel" % worst_c, worst_c < 1e-6, "< 1e-6")

    # distances
    t1 = time.time()
    D, st = gram_distances(G, F)
    C.info("Gram", "%.1f min" % ((time.time() - t1) / 60))
    C.add("self-distances zero", "%.2e" % st["self"], st["self"] < 1e-6, "< 1e-6")
    C.add("symmetric", "%.2e" % st["asym"], st["asym"] < 1e-6, "< 1e-6")
    C.add("no negative squared distances", "min off-diagonal D^2 %.3e" % st["min_off"], st["min_off"] > -1e-6, "> -1e-6")
    sub = np.sort(rng.choice(mols, 200, replace=False))
    ref = pdist(np.asarray(G[sub], dtype=np.float64)) / H ** 1.5
    r_ = float(np.max(np.abs(squareform(D[np.ix_(sub, sub)], checks=False) - ref) / ref))
    C.add("blocked Gram == scipy pdist, 200-molecule subset", "%.2e" % r_, r_ < 1e-4, "< 1e-4")
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    ce = max(abs(D[idx[a], idx[b]] - np.sqrt(float(((G[idx[a]].astype(np.float64) - G[idx[b]].astype(np.float64)) ** 2).sum()) / H ** 3))
             / D[idx[a], idx[b]] for a, b in zip(NP.clone_a, NP.clone_b))
    C.add("control pairs: blocked Gram vs exact float64 (%s)" % EXP_NCP, "%.2e" % ce, ce < 1e-3, "< 1e-3")
    t = rng.integers(0, n, (10000, 3)); t = t[(t[:, 0] != t[:, 1]) & (t[:, 1] != t[:, 2]) & (t[:, 0] != t[:, 2])]
    viol = float((D[t[:, 0], t[:, 2]] - (D[t[:, 0], t[:, 1]] + D[t[:, 1], t[:, 2]])).max())
    C.add("triangle inequality, %d triples: max violation" % len(t), "%.2e" % viol, viol < 1e-6, "<= 1e-6")
    pil, piu, pkeep = pilot_background(); pr = np.array([idx[c] for c in pil])
    pcut = float(np.percentile(D[np.ix_(pr, pr)][piu][pkeep], 1))
    V3 = pd.read_csv(os.path.join(OUT, "VXV3_cdr3.csv")).set_index("descriptor").loc["grid_F1C3_7ch", "pilot_cut"]
    C.info("pilot cut vs VXV3 grid_F1C3_7ch pilot cut (VXV3 used the truncating VX5b boxes; reported, A13.11.2)",
           "%.5f vs %.5f (%+.2e rel)" % (pcut, V3, pcut / V3 - 1))
    np.save(os.path.join(TMP, "VXC1_D_receptors.f32.npy"), D.astype(np.float32))
    np.save(os.path.join(TMP, "VXC1_D_molecules_condensed.f64.npy"), squareform(D[np.ix_(mols, mols)], checks=False))
    save_json(dict(atoms="CDR3 105-117", frame="F1", channels=NCH, sigma=SIGMA, h=H, boxes="out/VXC0_boxes.json", features=F,
                   grid=os.path.relpath(VXP("voxel_out/tmp/VXC1_armC.f16.npy"), ROOT), receptor_matrix=os.path.relpath(VXP("voxel_out/tmp/VXC1_D_receptors.f32.npy"), ROOT),
                   molecule_condensed=os.path.relpath(VXP("voxel_out/tmp/VXC1_D_molecules_condensed.f64.npy"), ROOT), molecule_rows=mols.tolist(),
                   minutes=round((time.time() - t0) / 60, 1)), os.path.join(OUT, "VXC1_build.json"))
    C.write()
