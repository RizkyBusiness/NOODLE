"""VXS1 (the voxel pipeline, A14.2.1): sigma 1.5 A arms - grid build and distances. Arm D = Arm B at sigma 1.5 (F1 loops, 7 channels,
VX5b F1 boxes); Arm E = Arm C at sigma 1.5 (F1 CDR3 105-117, 7 channels, VXC0 full-envelope boxes). All N receptors
in landmark-file order, float16 storage, blocked Gram in float64 (A7). Checks as VXH5 / VXC1 at sigma 1.5: typing, mass
conservation, analytic grid-free distance, fast == slow grid, permuted atom identities, Gram == pdist, control pairs exact,
triangle inequality; Arm E also the A13.11 box-cropping check. The pilot cut is reported (no sigma 1.5 reference exists).
The grid file is deleted at the end (A13.10); the distance matrices are kept. No state data.
Writes tmp/VXS1_{D,E}_D_receptors.f32.npy, tmp/VXS1_{D,E}_D_molecules_condensed.f64.npy, out/VXS1_{D,E}_build.json,
out/VXS1_{D,E}_manifest.csv; checks/VXS1_{D,E}_checks.csv.
usage: python pipeline/code/vxs1_build.py D|E
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
from vxv_common import boxes_for, rep_place, bench_place, build_grid, structure_atoms, shared_keys, H, BEN
from vxv1_build import pilot_background
from vxc_common import load_boxes, build_dataset, gram_distances, occ_index

SIG, NCH = 1.5, 7
ARMS = {"D": dict(atoms="loops", boxes="VX5b", like="Arm B (VXH5)"), "E": dict(atoms="cdr3", boxes="VXC0", like="Arm C (VXC1)")}

if __name__ == "__main__":
    arm = sys.argv[1]; A = ARMS[arm]; ATOMS = A["atoms"]; t0 = time.time()
    BX = boxes_for("F1") if A["boxes"] == "VX5b" else load_boxes()
    GRID = os.path.join(TMP, "VXS1_%s.f16.npy" % arm)
    C = Checks("VXS1_%s" % arm, stop_on_fail=True)
    C.info("status", "A14 VXS1 Arm %s = %s at sigma %.1f A; no state data read" % (arm, A["like"], SIG))
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; n = len(ids); idx = {c: k for k, c in enumerate(ids)}
    FR = np.load(os.path.join(OUT, "VX1_frames.npz"), allow_pickle=True); ismol = FR["is_molecule"]; mols = np.where(ismol)[0]
    C.info("configuration", "F1 per chain, %s, sigma %.1f, voxel %.1f, 7 channels, %s boxes %s" % (ATOMS, SIG, H, A["boxes"], {k: s for k, (_, s) in BX.items()}))
    G, na, nb, mass, occ, F = build_dataset(range(n), ids, BX, ATOMS, NCH, GRID, sigma=SIG)
    C.info("build", "%.1f min, %.1f GB, %d features" % ((time.time() - t0) / 60, os.path.getsize(GRID) / 1e9, F))
    C.add("typing complete", nb, nb == 0, "== 0")
    rel_all = np.abs(mass / na - 1)
    C.info("mass loss over all receptors (float64): max / n > 0.5 %", "%.4f %% / %d" % (100 * rel_all.max(), int((rel_all > 5e-3).sum())))
    rng = np.random.default_rng(141 if arm == "D" else 151)
    rows = rng.choice(mols, 200, replace=False)
    rel = max(abs(float(G[r][occ].astype(np.float64).sum()) / na[r] - 1) for r in rows)
    C.add("mass conservation, 200 molecules (stored)", "%.4f %%" % (100 * rel), rel < 5e-3, "< 0.5 %")
    job = lambda i: (VXP("structures/%s.pdb") % ids[i], rep_place(i, "F1"))
    errs = []
    for a, b in rng.choice(mols, (12, 2), replace=False):
        g = float(((G[a].astype(np.float64) - G[b].astype(np.float64)) ** 2).sum()) / H ** 3
        A_ = structure_atoms(*job(a), ATOMS, NCH); B_ = structure_atoms(*job(b), ATOMS, NCH)
        errs.append(g / sum(vg.analytic_d2(xa, wa, xb, wb, SIG, channels=range(NCH)) for (xa, wa), (xb, wb) in zip(A_, B_)) - 1)
    C.add("analytic grid-free distance, 12 pairs", "%.4f (%+.4f..%+.4f)" % (np.abs(errs).max(), min(errs), max(errs)),
          np.abs(errs).max() < 0.02, "< 2 %")
    worst = 0.0
    for r in rng.choice(mols, 3, replace=False):
        for (X, W), (key, (lo, shp)) in zip(structure_atoms(*job(r), ATOMS, NCH), BX.items()):
            worst = max(worst, float(np.abs(vg.build(X, W, SIG, H, lo, shp) - vg.build_slow(X, W, SIG, H, lo, shp)).max()))
    C.add("fast grid == slow loop reference, 3 molecules, both boxes", "%.2e" % worst, worst < 1e-5, "< 1e-5")
    # permuted atom identities (A13.9 rule: first molecule with mass in every channel 2-7, both boxes)
    for pm_ in mols:
        parts = structure_atoms(*job(int(pm_)), ATOMS, NCH)
        if (np.vstack([w for _, w in parts])[:, 1:7].sum(0) > 0).all():
            break
    C.info("permuted-atom control molecule", ids[int(pm_)])
    prng = np.random.default_rng(142)
    G0 = np.concatenate([vg.build(X, W, SIG, H, lo, shp).reshape(NCH, -1) for (X, W), (lo, shp) in zip(parts, BX.values())], 1)
    Gp = np.concatenate([vg.build(X, W[prng.permutation(len(W))], SIG, H, lo, shp).reshape(NCH, -1)
                         for (X, W), (lo, shp) in zip(parts, BX.values())], 1)
    dch = [float(np.abs(G0[c] - Gp[c]).max()) for c in range(1, 7)]
    C.add("permuted atom identities: channel 1 unchanged", "%.1e" % float(np.abs(G0[0] - Gp[0]).max()),
          float(np.abs(G0[0] - Gp[0]).max()) < 1e-12, "< 1e-12")
    C.add("permuted atom identities: channels 2-7 all change", " ".join("%.3f" % v for v in dch), min(dch) > 1e-4, "each > 1e-4")
    pd.DataFrame({"row": range(n), "clone_id": ids, "is_molecule": ismol, "atoms": na, "mass_occupancy": mass}) \
      .to_csv(os.path.join(OUT, "VXS1_%s_manifest.csv" % arm), index=False)
    if arm == "E":                                         # A13.11 box-cropping check, at sigma 1.5
        lm3 = load_lm3()
        R = pd.read_csv(os.path.join(OUT, "VXV0_receptors.csv")); R = R[R.cdr3_complete].reset_index(drop=True)
        PC = [os.path.join(BEN, "fixed_%s_crystal.pdb" % e) for e in R.entry]; PM = [os.path.join(BEN, m) for m in R.model]
        SH = [shared_keys(a, b) for a, b in zip(PC, PM)]
        BIG = {k: (lo - 10.0, tuple(int(x) + 20 for x in shp)) for k, (lo, shp) in BX.items()}; mats = []
        for nm_, bx in (("Arm E boxes", BX), ("expanded reference", BIG)):
            jobs = [(p, bench_place(p, "F1", lm3), bx, ATOMS, NCH, sh, None, False, SIG) for p, sh in zip(PC + PM, SH + SH)]
            with Pool(8) as pool:
                res = pool.map(build_grid, jobs)
            oc = occ_index(bx, NCH)[0]; ml = max(1 - float(r[0][oc].sum()) / r[1] for r in res)
            C.add("box cropping: benchmark mass loss in %s (%d grids)" % (nm_, len(PC) + len(PM)), "%.1e" % ml, ml <= 1e-9, "<= 1e-9")
            V = np.array([r[0] for r in res]); nr_ = len(PC)
            mats.append((squareform(pdist(V[:nr_])) / H ** 1.5, squareform(pdist(V[nr_:])) / H ** 1.5))
        iu = np.triu_indices(len(PC), 1)
        wc = max(float(np.max(np.abs(a_[iu] - b_[iu]) / b_[iu])) for a_, b_ in zip(mats[0], mats[1]))
        C.add("box cropping: %d VXV receptors, Arm E boxes == expanded boxes (C and M)" % len(PC), "%.2e rel" % wc, wc < 1e-6, "< 1e-6")
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
    C.info("pilot cut (sigma 1.5; reported, no reference)", "%.5f" % pcut)
    np.save(os.path.join(TMP, "VXS1_%s_D_receptors.f32.npy" % arm), D.astype(np.float32))
    np.save(os.path.join(TMP, "VXS1_%s_D_molecules_condensed.f64.npy" % arm), squareform(D[np.ix_(mols, mols)], checks=False))
    del G; os.remove(GRID)
    C.info("grid file deleted (A13.10)", os.path.basename(GRID))
    save_json(dict(arm=arm, atoms=ATOMS, frame="F1", channels=NCH, sigma=SIG, h=H, boxes=A["boxes"], features=F, pilot_cut=pcut,
                   receptor_matrix=os.path.relpath(VXP("voxel_out/tmp/VXS1_%s_D_receptors.f32.npy") % arm, ROOT),
                   molecule_condensed=os.path.relpath(VXP("voxel_out/tmp/VXS1_%s_D_molecules_condensed.f64.npy") % arm, ROOT), molecule_rows=mols.tolist(),
                   minutes=round((time.time() - t0) / 60, 1)), os.path.join(OUT, "VXS1_%s_build.json" % arm))
    C.write()
