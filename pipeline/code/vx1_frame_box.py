"""VX1 (the voxel pipeline): one whole-molecule frame per model, and the frozen box.

Fit: Kabsch on all ten framework landmarks (the landmark file.lm, Ca of IMGT 23 41 89 104 118 on A then B) onto the reference receptor (config)'s,
row convention R = U diag(1,1,sign det(U Vt)) Vt, applied as (X - c_mol) @ R + c_ref  (Kabsch 1976).
Then one fixed transform, the H15 intrinsic frame of the reference (local = (p - O) @ AX.T; z toward pMHC).
H15 is defined in the D2 frame, which is the reference receptor (config)'s own coordinates (checked at VX0), so it applies directly.
Box coordinates:  X_box = ((X - c_mol) @ R + c_ref - O) @ AX.T.  The H15 part is the same for every model and
cannot change any distance.

Box: heavy atoms of 500 random molecules (seed 0, drawn from the M the reference cluster-test procedure-rule representatives), 0.1-99.9
percentile envelope per axis, + 4 A margin, rounded outward to whole voxels at 1.0 A (the hard floor), frozen.
Outside-box fraction is then measured on all N models.

For the record (not used downstream): 123-anchor RMSD after the ten-landmark fit; the per-chain five-landmark
fits exactly as the reference method/code/vec1_features.py; and the Va/Vb relative rotation angle vs the reference
(angle of Ra^-1 Rb from the two per-chain fits; cf. ABangle/TRangle, Dunbar et al. 2014).

Writes <voxel out>/out/VX1_frames.npz, VX1_box.json, VX1_frame_table.csv.gz, VX1_heavy_atoms.csv.gz;
<voxel out>/checks/VX1_checks.csv.   usage: python pipeline/code/vx1_frame_box.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
EXP_NMOL = expected("n_molecules")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NBM = expected("n_benchmark_models")
EXP_NANC = expected("n_anchors")
EXP_FIT_A = expected("reference_perchain_fit_alpha_A")
EXP_FIT_B = expected("reference_perchain_fit_beta_A")

MARGIN, HMIN, NBOX, SEED = 4.0, 1.0, 500, 0
RES = (1.0, 1.5, 2.0, 2.5)
NCH = 7


def to_box(X, R, t, O, AX):
    return ((X @ R + t) - O) @ AX.T


def outside_worker(args):
    cid, R, t, O, AX, lo, hi = args
    xyz, _ = parse_heavy(VXP("structures/%s.pdb") % cid)
    B = to_box(xyz, R, t, O, AX)
    out = ~((B >= lo) & (B < hi)).all(1)
    return cid, len(xyz), int(out.sum())


def box_worker(args):
    cid, R, t, O, AX = args
    xyz, _ = parse_heavy(VXP("structures/%s.pdb") % cid)
    return to_box(xyz, R, t, O, AX).astype(np.float32)


def rot_angle(R):
    return float(np.degrees(np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))))


if __name__ == "__main__":
    C = Checks("VX1")
    L = np.load(LM1, allow_pickle=True)
    ids = [str(c) for c in L["clone_id"]]; n = len(ids)
    lm, an = L["lm"].astype(float), L["anch"].astype(float)
    ri = ids.index(REF_ID); T, TA = lm[ri], an[ri]
    H = np.load(H15, allow_pickle=True); O, AX = H["origin"].astype(float), H["axes"].astype(float)

    # molecules: the the reference cluster-test procedure collapse rule (first clone_id in landmark-file order per prot_key)
    ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz"))
    SL = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz"))
    M0 = (pd.DataFrame({"clone_id": ids}).merge(SL[["clone_id", "clone_key"]], on="clone_id", how="left")
          .merge(ID.drop_duplicates("clone_key")[["clone_key", "prot_key"]], on="clone_key", how="left"))
    is_mol = (M0.prot_key.notna().values & ~M0.prot_key.duplicated().values)
    assert is_mol.sum() == EXP_NMOL, "molecules %d; config: %s" % (is_mol.sum(), EXP_NMOL)

    # ------------------------------------------------------------ fits
    Rs, ts = np.zeros((n, 3, 3)), np.zeros((n, 3))
    f_lm, f_an, f_lmA, f_lmB, f_a, f_b, ang_ab = (np.zeros(n) for _ in range(7))
    for i in range(n):
        R, cx, cy = kabsch(lm[i], T)
        Rs[i], ts[i] = R, cy - cx @ R
        P = lm[i] @ R + ts[i]
        d2 = ((P - T) ** 2).sum(1)
        f_lm[i] = np.sqrt(d2.mean()); f_lmA[i] = np.sqrt(d2[:5].mean()); f_lmB[i] = np.sqrt(d2[5:].mean())
        f_an[i] = np.sqrt((((an[i] @ R + ts[i]) - TA) ** 2).sum(1).mean())
        Ra, ca, cya = kabsch(lm[i][:5], T[:5]); Rb, cb, cyb = kabsch(lm[i][5:], T[5:])
        f_a[i] = np.sqrt((((lm[i][:5] - ca) @ Ra - (T[:5] - cya)) ** 2).sum(1).mean())
        f_b[i] = np.sqrt((((lm[i][5:] - cb) @ Rb - (T[5:] - cyb)) ** 2).sum(1).mean())
        ang_ab[i] = rot_angle(Ra.T @ Rb)

    # reference maps to itself
    C.add("reference R == identity", "%.2e" % np.abs(Rs[ri] - np.eye(3)).max(),
          np.abs(Rs[ri] - np.eye(3)).max() < 1e-6, "< 1e-6")
    C.add("reference translation == 0", "%.2e A" % np.abs(ts[ri]).max(), np.abs(ts[ri]).max() < 1e-6, "< 1e-6 A")
    C.add("all R proper rotations", "max |det-1| %.1e" % np.abs(np.linalg.det(Rs) - 1).max(),
          np.abs(np.linalg.det(Rs) - 1).max() < 1e-6, "< 1e-6")

    def dist(name, v, sub):
        for lab, m in (("molecules_%s" % EXP_NMOL, is_mol), ("models_%s" % EXP_NS, np.ones(n, bool))):
            x = v[m]
            C.info("%s [%s] median / p95 / max" % (name, lab),
                   "%.3f / %.3f / %.3f A" % (np.median(x), np.percentile(x, 95), x.max()))
    dist("ten-landmark fit RMSD", f_lm, None)
    dist("  of which alpha landmarks", f_lmA, None)
    dist("  of which beta landmarks", f_lmB, None)
    dist("%s-anchor RMSD after ten-landmark fit" % EXP_NANC, f_an, None)
    dist("per-chain 5-landmark fit alpha (vec1 code)", f_a, None)
    dist("per-chain 5-landmark fit beta (vec1 code)", f_b, None)
    C.info("Va/Vb relative rotation vs reference [molecules] median / p95 / max",
           "%.2f / %.2f / %.2f deg" % (np.median(ang_ab[is_mol]), np.percentile(ang_ab[is_mol], 95),
                                       ang_ab[is_mol].max()))
    rho = pd.Series(f_lm[is_mol]).corr(pd.Series(ang_ab[is_mol]), method="spearman")
    C.info("Spearman(ten-landmark RMSD, Va/Vb rotation) [molecules]", "%.3f" % rho)
    ok_a, ok_b = abs(np.median(f_a) - EXP_FIT_A) < 0.01, abs(np.median(f_b) - EXP_FIT_B) < 0.01
    C.add("per-chain fits reproduce the stated %s / %s A" % (EXP_FIT_A, EXP_FIT_B), "%.3f / %.3f A" % (np.median(f_a), np.median(f_b)),
          ok_a and ok_b, "within 0.01 A")
    worst = np.argsort(-np.where(is_mol, f_an, -1))[:10]
    C.info("ten worst molecules by %s-anchor RMSD" % EXP_NANC,
           "; ".join("%s %.2f" % (ids[i], f_an[i]) for i in worst))

    # ------------------------------------------------------------ box from 500 random molecules
    rng = np.random.default_rng(SEED)
    mol_idx = np.where(is_mol)[0]
    samp = np.sort(rng.choice(mol_idx, NBOX, replace=False))
    with Pool(8) as pool:
        S = pool.map(box_worker, [(ids[i], Rs[i], ts[i], O, AX) for i in samp])
    P = np.concatenate(S)
    lo_p, hi_p = np.percentile(P, 0.1, axis=0), np.percentile(P, 99.9, axis=0)
    lo = np.floor((lo_p - MARGIN) / HMIN) * HMIN
    hi = np.ceil((hi_p + MARGIN) / HMIN) * HMIN
    ext = hi - lo
    C.info("box sample: molecules / heavy atoms", "%d / %d (seed %d)" % (NBOX, len(P), SEED))
    C.info("atom envelope 0.1-99.9 pct, lo", np.round(lo_p, 2).tolist())
    C.info("atom envelope 0.1-99.9 pct, hi", np.round(hi_p, 2).tolist())
    C.info("frozen box lo / hi (A)", "%s / %s" % (lo.tolist(), hi.tolist()))
    C.info("frozen box extent (A)", ext.tolist())
    sizes = {}
    for h in RES:
        shp = [int(np.ceil(e / h - 1e-9)) for e in ext]
        nv = int(np.prod(shp))
        sizes[str(h)] = dict(shape=shp, voxels=nv, features_7ch=nv * NCH,
                             dataset_GB_7ch_float16=round(EXP_NMOL * NCH * nv * 2 / 1e9, 2),
                             dataset_GB_8ch_float16=round(EXP_NMOL * 8 * nv * 2 / 1e9, 2),
                             exact_multiple=bool(all(abs(e / h - round(e / h)) < 1e-9 for e in ext)))
        C.info("at %.1f A: shape / voxels / GB (7ch, 8ch)" % h,
               "%s / %d / %.2f, %.2f" % (shp, nv, sizes[str(h)]["dataset_GB_7ch_float16"],
                                          sizes[str(h)]["dataset_GB_8ch_float16"]))
    nv1 = sizes["1.0"]["voxels"]
    for npil in (1500, 2000):
        C.info("pilot at 1.0 A, %d + %s benchmark structures, 7ch float16" % (npil, EXP_NBP + EXP_NBM),
               "{:.2f} GB".format((npil + EXP_NBP + EXP_NBM) * NCH * nv1 * 2 / 1e9))

    # ------------------------------------------------------------ outside-box fraction on all N
    with Pool(8) as pool:
        res = pool.map(outside_worker, [(ids[i], Rs[i], ts[i], O, AX, lo, hi) for i in range(n)], chunksize=20)
    HA = pd.DataFrame(res, columns=["clone_id", "heavy_atoms", "outside"])
    assert list(HA.clone_id) == ids
    tot, outn = int(HA.heavy_atoms.sum()), int(HA.outside.sum())
    insamp = np.isin(np.arange(n), samp)
    fo = outn / tot
    C.add("heavy atoms outside frozen box (all {:,} models)".format(EXP_NS), "%d / %d = %.4f %%" % (outn, tot, 100 * fo),
          fo < 1e-3, "< 0.1 %")
    C.info("outside fraction, models NOT in the box sample",
           "%.4f %%" % (100 * HA.outside[~insamp].sum() / HA.heavy_atoms[~insamp].sum()))
    C.info("models with any atom outside", int((HA.outside > 0).sum()))
    C.info("heavy atoms per model median / min / max",
           "%d / %d / %d" % (HA.heavy_atoms.median(), HA.heavy_atoms.min(), HA.heavy_atoms.max()))
    HA.to_csv(os.path.join(OUT, "VX1_heavy_atoms.csv.gz"), index=False)

    # ------------------------------------------------------------ invariance to input rigid motion
    rng2 = np.random.default_rng(1)
    worst_inv = 0.0
    for i in [ri] + list(rng2.choice(n, 4, replace=False)):
        xyz, _ = parse_heavy(VXP("structures/%s.pdb") % ids[i])
        q = rng2.normal(size=4); q /= np.linalg.norm(q)
        a, b, c, d = q
        Q = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
        sh = rng2.uniform(-50, 50, 3)
        lm2, xyz2 = lm[i] @ Q + sh, xyz @ Q + sh
        R2, cx2, cy2 = kabsch(lm2, T)
        B1 = to_box(xyz, Rs[i], ts[i], O, AX); B2 = to_box(xyz2, R2, cy2 - cx2 @ R2, O, AX)
        worst_inv = max(worst_inv, float(np.abs(B1 - B2).max()))
    C.add("fit invariant to input rotation+translation (5 models)", "%.2e A" % worst_inv, worst_inv < 1e-4,
          "< 1e-4 A")

    # ------------------------------------------------------------ outputs
    Rbox = Rs @ AX.T                       # combined: X_box = X @ Rbox + tbox
    tbox = (ts - O) @ AX.T
    # numpy 2.2 matmul on this arm64 build emits divide/overflow/invalid RuntimeWarnings on these small
    # products; confirm the numbers are nonetheless exact by recomputing without BLAS (einsum).
    e1 = max(float(np.abs(Rbox - np.einsum("nij,kj->nik", Rs, AX)).max()),
             float(np.abs(tbox - np.einsum("ij,kj->ik", ts - O, AX)).max()))
    e2 = 0.0
    for i in [ri] + list(rng2.choice(n, 4, replace=False)):
        xyz, _ = parse_heavy(VXP("structures/%s.pdb") % ids[i])
        ref_ = np.einsum("ij,kj->ik", np.einsum("ij,jk->ik", xyz, Rs[i]) + ts[i] - O, AX)
        e2 = max(e2, float(np.abs(to_box(xyz, Rs[i], ts[i], O, AX) - ref_).max()))
    fin = bool(np.isfinite(Rbox).all() and np.isfinite(tbox).all() and np.isfinite(f_lm).all()
               and np.isfinite(f_an).all() and np.isfinite(P).all())
    C.add("matmul warnings spurious: frames == einsum, all finite", "frames %.1e, atoms %.1e, finite %s"
          % (e1, e2, fin), fin and e1 < 1e-12 and e2 < 1e-9, "< 1e-12 / < 1e-9 A, finite")
    np.savez_compressed(os.path.join(OUT, "VX1_frames.npz"), clone_id=np.array(ids), is_molecule=is_mol,
                        R_fit=Rs, t_fit=ts, R_box=Rbox, t_box=tbox, h15_origin=O, h15_axes=AX,
                        rmsd_lm10=f_lm, rmsd_lm10_alpha=f_lmA, rmsd_lm10_beta=f_lmB, rmsd_anch123=f_an,
                        fit_alpha5=f_a, fit_beta5=f_b, vavb_rot_deg=ang_ab, box_sample=samp,
                        note=np.array("X_fit = X @ R_fit + t_fit ; X_box = X @ R_box + t_box = (X_fit - h15_origin) @ h15_axes.T"))
    pd.DataFrame({"clone_id": ids, "is_molecule": is_mol, "rmsd_lm10": f_lm, "rmsd_lm10_alpha": f_lmA,
                  "rmsd_lm10_beta": f_lmB, "rmsd_anch123": f_an, "fit_alpha5": f_a, "fit_beta5": f_b,
                  "vavb_rot_deg": ang_ab}).to_csv(os.path.join(OUT, "VX1_frame_table.csv.gz"), index=False)
    save_json(dict(frozen=True, frame="intrinsic axes of the reference receptor (z toward pMHC), origin at H15 origin",
                   lo=lo.tolist(), hi=hi.tolist(), extent_A=ext.tolist(), grid_step_rounded_to_A=HMIN,
                   margin_A=MARGIN, envelope_pct=[0.1, 99.9], envelope_lo=lo_p.tolist(),
                   envelope_hi=hi_p.tolist(), sample_n=NBOX, sample_seed=SEED,
                   sample_from="{:,} vec2-rule molecules".format(EXP_NMOL), voxel_convention="voxel k spans [lo+k*h, lo+(k+1)*h)",
                   outside_fraction_all_models=fo, sizes_by_resolution=sizes),
              os.path.join(OUT, "VX1_box.json"))
    C.write()
