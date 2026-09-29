"""VX2b (the voxel pipeline): blur sweep on the pilot, and the resolution decision (DESIGN VX2b + amendments A3.3, A4).

For sigma in {0.6, 0.8, 1.0, 1.4, 2.0} A, each built directly from the atoms (A4; sigma 0.6 reuses the VX2a grid),
on the 1.0 A pilot grid, channels 1-7:
  - 1 % cut = 1st percentile of the grid distances over all pairs of the 2,000 pilot repertoire molecules
  - crystal-model grid distance for the accepted pairs, on shared atoms (primary, A2.1/A3.1) and on all atoms
  - ratio = median crystal-model distance / cut;  control pairs inside the pilot: distance / cut and background pct
Rule (fixed in DESIGN 6.1, before the sweep): smallest sigma whose shared-atom ratio is within 0.02 of the sweep
minimum; frame floor 0.63 A (A3.3); production voxel = largest of {1.0,1.5,2.0,2.5} <= sigma, never below 1.0 (A1.4).
Distances: D^2 = |a|^2 + |b|^2 - 2 a.b accumulated over feature blocks (the VX3 blocked-Gram route), / h^3.

usage: python pipeline/code/vx2b_sweep.py            (the one real run)
                             python pipeline/code/vx2b_sweep.py --dryrun DIR (code-path test: 30 molecules,
                                                                          sigma 0.7/1.2, writes only to DIR)
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, time
import numpy as np, pandas as pd

CTX_EXISTING_RATIO = __import__("vxpaths").CFG.get("reference_method", {}).get("existing_method_ratio")   # optional: the reference method's crystal-error/cut, context only
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NBM = expected("n_benchmark_models")
import vxgrid as vg
import vx2a_pilot as PA

DRY = "--dryrun" in sys.argv
SWEEP = [0.7, 1.2] if DRY else [0.6, 0.8, 1.0, 1.4, 2.0]
FLOOR, TOL, HARD_STOP, H = 0.63, 0.02, 0.60, PA.H
VOXELS = (1.0, 1.5, 2.0, 2.5)
ODIR = sys.argv[sys.argv.index("--dryrun") + 1] if DRY else None
OUTD, CHKD, TMPD = (ODIR, ODIR, ODIR) if DRY else (OUT, CHK, TMP)
M = pd.read_csv(os.path.join(OUT, "VX2a_manifest.csv"))
NV = int(np.prod(PA.SHAPE)); NF = vg.NPRIMARY * NV


def make_jobs():
    """reconstruct the VX2a job list, row for row, from the manifest"""
    lm3 = load_lm3()
    B = lm3["B"]
    frames = {}
    def fr(p):
        if p not in frames:
            frames[p] = PA.bench_frame(p, lm3["parse"], lm3["LMPOS"])
        return frames[p]
    jobs = []
    for _, r in M.iterrows():
        if r.kind == "rep":
            i = int(r.lm1_index)
            jobs.append(("rep", VXP("structures/%s.pdb") % r.clone_id, PA.FR["R_box"][i], PA.FR["t_box"][i], None))
        elif r.kind in ("crystal_full", "model_full"):
            p = os.path.join(PA.BEN, "fixed_%s_crystal.pdb" % r.entry if r.kind == "crystal_full" else r.model)
            jobs.append(("full", p, *fr(p), None))
        else:
            pc, pm = os.path.join(PA.BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(PA.BEN, r.model)
            keep = {tuple(a) for a in parse_heavy(pc)[1]} & {tuple(a) for a in parse_heavy(pm)[1]}
            p = pc if r.kind == "crystal_shared" else pm
            jobs.append(("shared", p, *fr(p), keep))
    return jobs


def build_one(args):
    job, sigma = args
    X, W, meta, _ = PA.structure(job)
    return vg.build(X, W, sigma, H, PA.LO, PA.SHAPE).astype(np.float16), len(X)


def gram_d(G, rows, block=72500):
    """blocked Gram distances among `rows` of memmap G, channels 1-7; returns D (n,n) and min raw D^2."""
    n = len(rows)
    Gm = np.zeros((n, n))
    flat = G.reshape(G.shape[0], -1)
    for s in range(0, NF, block):
        Xc = np.asarray(flat[rows, s:min(s + block, NF)], dtype=np.float32)
        Gm += (Xc @ Xc.T).astype(np.float64)
    nr = np.diag(Gm).copy()
    D2 = (nr[:, None] + nr[None] - 2 * Gm) / H ** 3
    np.fill_diagonal(D2, 0.0)
    mn = float(D2.min())
    return np.sqrt(np.clip(D2, 0, None)), mn, nr


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VX2b_dryrun" if DRY else "VX2b")
    if DRY:
        CHK_SAVE = CHK
        import vxlib; vxlib.CHK = CHKD
    jobs = make_jobs()
    rep = np.where(M.kind == "rep")[0]
    if DRY:
        rep = rep[:30]
    use = np.concatenate([rep, np.where(M.kind != "rep")[0]])
    # the manifest reconstruction must reproduce the stored pilot exactly (to float16 rounding)
    Gp = np.load(PA.GRID + ".npy", mmap_mode="r")
    w = 0.0
    for r in [int(rep[0]), int(np.where(M.kind == "crystal_shared")[0][0]), int(np.where(M.kind == "model_full")[0][0])]:
        g, _ = build_one((jobs[r], PA.SIGMA))
        w = max(w, float(np.abs(g.astype(np.float32) - Gp[r].astype(np.float32)).max()))
    C.add("job reconstruction reproduces stored VX2a pilot (3 rows)", "%.1e" % w, w == 0.0, "== 0 (bitwise)")

    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    pos = {c: j for j, c in enumerate(M.clone_id[rep])}
    ctrl = [(pos[a], pos[b], a, b) for a, b in zip(NP.clone_a, NP.clone_b) if a in pos and b in pos]
    Bsh = M[M.kind == "crystal_shared"].merge(M[M.kind == "model_shared"], on=["entry", "model"], suffixes=("_c", "_m"))
    Bfu = (M[M.kind == "crystal_full"][["entry", "row"]].merge(M[M.kind == "crystal_shared"][["entry", "model"]], on="entry")
           .merge(M[M.kind == "model_full"][["model", "row"]], on="model", suffixes=("_c", "_m")))
    assert len(Bsh) == EXP_NBP and len(Bfu) == EXP_NBP, "benchmark pairs %d / %d; config: %s" % (len(Bsh), len(Bfu), EXP_NBP)

    rng = np.random.default_rng(21)
    ana_pairs = rng.choice(rep, (12, 2), replace=len(rep) < 24)
    mass_rows = rng.choice(rep, min(200, len(rep)), replace=False)
    sweep, pairs_out, ctrl_out = [], [], []

    def run_sigma(sigma, label):
        ts = time.time()
        if abs(sigma - PA.SIGMA) < 1e-12 and not DRY:
            G = Gp
        else:
            path = os.path.join(TMPD, "VX2b_h1.0_s%.2f.f16.npy" % sigma)
            G = np.lib.format.open_memmap(path, mode="w+", dtype=np.float16, shape=(len(M), vg.NCH) + PA.SHAPE)
            with Pool(8) as pool:
                for j, (g, na) in zip(use, pool.imap(build_one, [(jobs[j], sigma) for j in use], chunksize=4)):
                    G[j] = g
            G.flush()
        # checks at this sigma (A4.2)
        rel = max(abs(float(G[r, 0].astype(np.float64).sum()) / M.heavy_atoms[r] - 1) for r in mass_rows)
        relb = max(abs(float(G[r, 0].astype(np.float64).sum()) / M.heavy_atoms[r] - 1) for r in np.where(M.kind != "rep")[0])
        C.add("%s mass conservation, %d molecules: max" % (label, len(mass_rows)), "%.4f %%" % (100 * rel), rel < 5e-3, "< 0.5 %")
        C.add("%s mass conservation, %s benchmark grids: max" % (label, 3 * EXP_NBP + EXP_NBM), "%.4f %%" % (100 * relb), relb < 5e-3, "< 0.5 %")
        errs = []
        for a, b in ana_pairs:
            XA, WA, _, _ = PA.structure(jobs[a]); XB, WB, _, _ = PA.structure(jobs[b])
            errs.append(vg.grid_d2(G[a], G[b], H) / vg.analytic_d2(XA, WA, XB, WB, sigma) - 1)
        errs = np.array(errs)
        C.add("%s analytic grid-free distance, 12 pairs: max |grid/exact - 1|" % label,
              "%.4f (%+.4f..%+.4f)" % (np.abs(errs).max(), errs.min(), errs.max()), np.abs(errs).max() < 0.02, "< 2 %")
        # background and cut
        D, mn, nr = gram_d(G, rep)
        C.add("%s blocked Gram: min raw D^2 relative to median |a|^2" % label, "%.1e" % (mn / np.median(nr)),
              mn / np.median(nr) > -1e-6, "> -1e-6")
        chk = [(0, 1), (2, 5), (len(rep) - 1, 3)]
        gerr = max(abs(D[i, j] - np.sqrt(vg.grid_d2(G[rep[i]], G[rep[j]], H))) / D[i, j] for i, j in chk)
        C.add("%s blocked Gram == direct grid distance (3 pairs)" % label, "%.1e" % gerr, gerr < 1e-4, "< 1e-4 rel")
        iu = np.triu_indices(len(rep), 1)
        bg = D[iu]
        cut = float(np.percentile(bg, 1))
        dsh = np.array([np.sqrt(vg.grid_d2(G[a], G[b], H)) for a, b in zip(Bsh.row_c, Bsh.row_m)])
        dfu = np.array([np.sqrt(vg.grid_d2(G[a], G[b], H)) for a, b in zip(Bfu.row_c, Bfu.row_m)])
        bgs = np.sort(bg)
        cd = np.array([D[i, j] for i, j, _, _ in ctrl]) if ctrl else np.array([])
        row = dict(sigma=sigma, n_background_pairs=len(bg), cut_p1=cut, bg_median=float(np.median(bg)),
                   crystal_median_shared=float(np.median(dsh)), crystal_median_full=float(np.median(dfu)),
                   ratio_shared=float(np.median(dsh)) / cut, ratio_full=float(np.median(dfu)) / cut,
                   n_control_in_pilot=len(ctrl),
                   control_within_cut=float((cd <= cut).mean()) if len(cd) else np.nan,
                   control_median_over_cut=float(np.median(cd) / cut) if len(cd) else np.nan,
                   control_median_bg_pct=float(np.median(100 * np.searchsorted(bgs, cd) / len(bgs))) if len(cd) else np.nan,
                   minutes=round((time.time() - ts) / 60, 2))
        for e, m_, a, b in zip(Bsh.entry, Bsh.model, dsh, dfu):
            pairs_out.append(dict(sigma=sigma, entry=e, model=m_, d_shared=a, d_full=b, over_cut_shared=a / cut))
        for (i, j, a, b), d in zip(ctrl, cd):
            ctrl_out.append(dict(sigma=sigma, clone_a=a, clone_b=b, d=d, over_cut=d / cut,
                                 bg_pct=100 * np.searchsorted(bgs, d) / len(bgs)))
        print("sigma %.2f done: cut %.4f ratio_shared %.4f ratio_full %.4f (%.1f min)"
              % (sigma, cut, row["ratio_shared"], row["ratio_full"], row["minutes"]), flush=True)
        return row

    for s in SWEEP:
        sweep.append(run_sigma(s, "sigma %.2f" % s))
    S = pd.DataFrame(sweep)

    # ------------------------------------------------------------ the rule, mechanically
    rmin = S.ratio_shared.min()
    pick = float(S.sigma[S.ratio_shared <= rmin + TOL].min())
    sigma_prod, floor_applied = pick, False
    if pick < FLOOR and not DRY:
        floor_applied = True
        sigma_prod = FLOOR
        S = pd.concat([S, pd.DataFrame([run_sigma(FLOOR, "floor sigma %.2f" % FLOOR)])], ignore_index=True)
    voxel = max([v for v in VOXELS if v <= sigma_prod + 1e-9], default=1.0)
    voxel = max(voxel, 1.0)
    rs = S.set_index("sigma").ratio_shared
    ratio_prod = float(rs[sigma_prod])
    seq = S[S.sigma.isin(SWEEP)].sort_values("sigma").ratio_shared.values
    dirs = np.sign(np.diff(seq)); dirs = dirs[dirs != 0]
    changes = int((np.diff(dirs) != 0).sum()) if len(dirs) > 1 else 0
    shape_ok = changes == 0 or (changes == 1 and dirs[0] < 0)
    C.add("ratio curve monotone or single minimum", " ".join("%.4f" % v for v in seq), shape_ok,
          "<= 1 direction change, down then up")
    C.info("rule: sweep minimum ratio (shared atoms)", "%.4f" % rmin)
    C.info("rule: smallest sigma within %.2f of the minimum" % TOL, "%.2f A" % pick)
    C.info("frame floor %.2f A applied" % FLOOR, floor_applied)
    C.info("production sigma / voxel", "%.2f A / %.1f A" % (sigma_prod, voxel))
    C.info("ratio at production sigma (shared / full)", "%.4f / %.4f"
           % (ratio_prod, float(S.set_index("sigma").ratio_full[sigma_prod])))
    shp = [int(np.ceil(e / voxel - 1e-9)) for e in PA.BX["extent_A"]]
    C.info("full-build dataset at production voxel, 8 ch float16", "{} -> {:.1f} GB"
           .format(shp, EXP_NMOL * 8 * np.prod(shp) * 2 / 1e9))
    C.add("provisional viability (DESIGN 6.5): ratio at production sigma", "%.4f" % ratio_prod,
          ratio_prod <= HARD_STOP, "<= 0.60")
    def downsample_check(sig, vox, label):
        """1.0 A grid at sig, block-summed to vox, vs a direct build at vox (3 molecules, rel L2)."""
        f = int(round(vox / H)); worst = 0.0
        for r in rep[:3]:
            X, W, _, _ = PA.structure(jobs[r])
            g1 = vg.build(X, W, sig, H, PA.LO, PA.SHAPE)[:7]
            s_ = [n // f for n in PA.SHAPE]
            ds = g1[:, :s_[0] * f, :s_[1] * f, :s_[2] * f].reshape(7, s_[0], f, s_[1], f, s_[2], f).sum((2, 4, 6))
            gd = vg.build(X, W, sig, vox, PA.LO, s_)[:7]
            worst = max(worst, float(np.linalg.norm(ds - gd) / np.linalg.norm(gd)))
        C.add("%s downsampled 1.0 A grid == direct build at %.1f A (3 molecules)" % (label, vox), "%.4f" % worst,
              worst < 0.01, "< 1 % rel L2")

    if DRY:
        downsample_check(2.0, 2.0, "dry-run code test, sigma 2.0:")
    elif voxel != 1.0:
        downsample_check(sigma_prod, voxel, "production")
    else:
        C.info("downsampling check", "not applicable: production voxel is 1.0 A (the pilot grid itself)")

    S.to_csv(os.path.join(OUTD, "VX2b_sweep.csv"), index=False)
    pd.DataFrame(pairs_out).to_csv(os.path.join(OUTD, "VX2b_crystal_pairs.csv"), index=False)
    pd.DataFrame(ctrl_out).to_csv(os.path.join(OUTD, "VX2b_control_pairs.csv"), index=False)
    save_json(dict(sweep=SWEEP, rule="smallest sigma with ratio_shared <= min + %.2f" % TOL, sweep_min=rmin,
                   sigma_rule=pick, frame_floor=FLOOR, floor_applied=floor_applied, sigma_production=sigma_prod,
                   voxel_production=voxel, ratio_production_shared=ratio_prod, hard_stop=HARD_STOP,
                   existing_method_ratio=CTX_EXISTING_RATIO, run_once=not DRY, minutes=round((time.time() - t0) / 60, 1)),
              os.path.join(OUTD, "VX2b_decision.json"))
    C.write()
