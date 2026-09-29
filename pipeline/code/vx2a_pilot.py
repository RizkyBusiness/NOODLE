"""VX2a (the voxel pipeline): pilot build at 1.0 A, sigma 0.6 A, on the VX1 frozen box.

Contents (DESIGN VX2a + amendments A1, A2):
  - 2,000 molecules drawn once at random from the M the reference cluster-test procedure-rule molecules (seed 0)
  - benchmark: the accepted crystals (lm3 selection, loaded via vxlib.load_lm3) and their distinct models,
    each on all its atoms ("full"), plus each crystal-model pair restricted to their shared atoms ("shared", A2.1)
Every structure goes through the same frame (ten-landmark Kabsch onto the reference receptor (config) + H15), box, typing and smearing.
Crystal/model landmarks come from lm3's parse(); heavy atoms from vxlib.parse_heavy (lm3's filters).

Storage: float16 memmap <voxel out>/tmp/VX2a_pilot_h1.0_s0.6.f16, shape (N, 8, 58, 50, 50); manifest
<voxel out>/out/VX2a_manifest.csv. Checks: <voxel out>/checks/VX2a_checks.csv.
usage: python pipeline/code/vx2a_pilot.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, time
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NBM = expected("n_benchmark_models")
import vxgrid as vg

H, SIGMA, NPILOT, SEED = 1.0, 0.6, 2000, 0
BEN = VXP("benchmark/")
BX = json.load(open(os.path.join(OUT, "VX1_box.json")))
LO = np.array(BX["lo"], float)
SHAPE = tuple(int(round(e / H)) for e in BX["extent_A"])
FR = np.load(os.path.join(OUT, "VX1_frames.npz"), allow_pickle=True)
IDS = [str(c) for c in FR["clone_id"]]
LMK = np.load(LM1, allow_pickle=True)
T = LMK["lm"][IDS.index(REF_ID)].astype(float)
O, AX = FR["h15_origin"], FR["h15_axes"]
GRID = os.path.join(TMP, "VX2a_pilot_h1.0_s0.6.f16")


def bench_frame(path, lm3parse, lmpos):
    fr = lm3parse(path)[0]
    L = np.array([fr[k] for k in lmpos], float)
    R, cx, cy = kabsch(L, T)
    t = cy - cx @ R
    return R @ AX.T, (t - O) @ AX.T


def structure(job):
    """job -> boxed coords, typing weights, meta. kinds: rep | full | shared_crystal | shared_model"""
    kind, path, Rb, tb, keep = job
    xyz, meta = parse_heavy(path)
    W, bad = vg.type_atoms(meta)
    if keep is not None:
        m = np.array([(a[0], a[1], a[2], a[3], a[4]) in keep for a in meta])
        xyz, W, meta = xyz[m], W[m], [a for a, k in zip(meta, m) if k]
    return xyz @ Rb + tb, W, meta, bad


def build_job(job):
    X, W, meta, bad = structure(job)
    G = vg.build(X, W, SIGMA, H, LO, SHAPE)
    out = int((~((X >= LO) & (X < LO + np.array(SHAPE) * H)).all(1)).sum())
    return G.astype(np.float16), len(X), G.reshape(G.shape[0], -1).sum(1), out, len(bad)


def type_only(path):
    _, meta = parse_heavy(path)
    W, bad = vg.type_atoms(meta)
    return W.sum(0), (W > 0).sum(0), bad


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VX2a")
    C.info("sigma / voxel / box shape", "%.2f A / %.1f A / %s" % (SIGMA, H, SHAPE))

    # ------------------------------------------------------------ typing completeness on all structures
    lm3 = load_lm3()
    B = lm3["B"]
    C.add("lm3 selection (via ast, no side effects): accepted pairs", len(B), len(B) == EXP_NBP, "== %s" % EXP_NBP)
    bench_paths = sorted({os.path.join(BEN, "fixed_%s_crystal.pdb" % e) for e in B.entry} |
                         {os.path.join(BEN, m) for m in B.model})
    rep_paths = [VXP("structures/%s.pdb") % c for c in IDS]
    with Pool(8) as pool:
        TY = pool.map(type_only, rep_paths + bench_paths, chunksize=20)
    bad = [b for r in TY for b in r[2]]
    C.add("atom typing complete: unmatched atoms ({:,} models + {} benchmark)".format(EXP_NS, EXP_NBP + EXP_NBM), len(bad), len(bad) == 0, "== 0")
    mass = np.sum([r[0] for r in TY[:len(rep_paths)]], 0); cnt = np.sum([r[1] for r in TY[:len(rep_paths)]], 0)
    for c, nm in enumerate(vg.CHANNELS):
        C.info("repertoire atoms / mass in channel %d %s" % (c + 1, nm), "%d / %.1f (%.1f %% of heavy atoms)"
               % (cnt[c], mass[c], 100 * cnt[c] / cnt[0]))

    # ------------------------------------------------------------ pilot sample and benchmark jobs
    mol_idx = np.where(FR["is_molecule"])[0]
    samp = np.sort(np.random.default_rng(SEED).choice(mol_idx, NPILOT, replace=False))
    jobs, man = [], []
    for i in samp:
        jobs.append(("rep", rep_paths[i], FR["R_box"][i], FR["t_box"][i], None))
        man.append(dict(kind="rep", clone_id=IDS[i], lm1_index=int(i), entry="", model=""))
    frames = {p: bench_frame(p, lm3["parse"], lm3["LMPOS"]) for p in bench_paths}
    for e in sorted(B.entry):
        p = os.path.join(BEN, "fixed_%s_crystal.pdb" % e)
        jobs.append(("full", p, *frames[p], None)); man.append(dict(kind="crystal_full", clone_id="", lm1_index=-1, entry=e, model=""))
    for m in sorted(B.model.unique()):
        p = os.path.join(BEN, m)
        jobs.append(("full", p, *frames[p], None)); man.append(dict(kind="model_full", clone_id="", lm1_index=-1, entry="", model=m))
    typing_mismatch, nshared = 0, []
    for _, r in B.sort_values("entry").iterrows():
        pc, pm = os.path.join(BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(BEN, r.model)
        _, mc = parse_heavy(pc); _, mm = parse_heavy(pm)
        keep = {tuple(a) for a in mc} & {tuple(a) for a in mm}          # chain, num, ins, resname, atom (A3.1)
        Wc, _ = vg.type_atoms(mc); Wm, _ = vg.type_atoms(mm)
        tc = {tuple(a): tuple(w) for a, w in zip(mc, Wc)}
        tm = {tuple(a): tuple(w) for a, w in zip(mm, Wm)}
        typing_mismatch += sum(tc[k] != tm[k] for k in keep)
        nshared.append(len(keep))
        for kind, p in (("crystal_shared", pc), ("model_shared", pm)):
            jobs.append(("shared", p, *frames[p], keep))
            man.append(dict(kind=kind, clone_id="", lm1_index=-1, entry=r.entry, model=r.model))
    C.info("shared atoms per pair median / min", "%d / %d" % (np.median(nshared), min(nshared)))
    C.add("crystal vs model typing identical on shared atoms (%s pairs)" % EXP_NBP, "%d atoms differ" % typing_mismatch,
          typing_mismatch == 0, "== 0")
    N = len(jobs)
    C.info("pilot rows: repertoire / crystal full / model full / shared pairs x2", "%d / %d / %d / %d (total %d)"
           % (NPILOT, len(B), B.model.nunique(), len(B), N))
    C.info("pilot grid file size", "%.2f GB" % (N * vg.NCH * np.prod(SHAPE) * 2 / 1e9))

    # ------------------------------------------------------------ build
    G = np.lib.format.open_memmap(GRID + ".npy", mode="w+", dtype=np.float16, shape=(N, vg.NCH) + SHAPE)
    with Pool(8) as pool:
        for j, (g, na, ms, out, nb) in enumerate(pool.imap(build_job, jobs, chunksize=4)):
            G[j] = g
            man[j].update(row=j, heavy_atoms=na, outside=out, **{"mass_%s" % c: float(m) for c, m in zip(vg.CHANNELS, ms)})
            if j % 250 == 0:
                print("built %d / %d  (%.0f s)" % (j, N, time.time() - t0), flush=True)
    G.flush()
    M = pd.DataFrame(man)
    M.to_csv(os.path.join(OUT, "VX2a_manifest.csv"), index=False)
    C.info("build time", "%.0f s" % (time.time() - t0))
    for k in ("crystal_full", "model_full"):
        s = M[M.kind == k]
        C.info("heavy atoms outside box: %s" % k, "%d / %d" % (s.outside.sum(), s.heavy_atoms.sum()))

    # ------------------------------------------------------------ mass conservation (stored float16 grids)
    rng = np.random.default_rng(11)
    rows = rng.choice(np.where(M.kind == "rep")[0], 200, replace=False)
    rel = np.array([abs(float(G[r, 0].astype(np.float64).sum()) / M.heavy_atoms[r] - 1) for r in rows])
    C.add("mass conservation, 200 random molecules: max |sum occ / atoms - 1|", "%.4f %%" % (100 * rel.max()),
          rel.max() < 5e-3, "< 0.5 %")
    relb = np.array([abs(float(G[r, 0].astype(np.float64).sum()) / M.heavy_atoms[r] - 1) for r in np.where(M.kind != "rep")[0]])
    C.add("mass conservation, all %s benchmark grids: max" % (3 * EXP_NBP + EXP_NBM), "%.4f %%" % (100 * relb.max()), relb.max() < 5e-3, "< 0.5 %")

    # ------------------------------------------------------------ analytic, grid-free reference (12 pairs)
    reprow = np.where(M.kind == "rep")[0]
    pr = rng.choice(reprow, (12, 2), replace=False)
    errs, errs64, errs4, f16 = [], [], [], []
    for a, b in pr:
        XA, WA, _, _ = structure(jobs[a]); XB, WB, _, _ = structure(jobs[b])
        exact = vg.analytic_d2(XA, WA, XB, WB, SIGMA)
        g16 = vg.grid_d2(G[a], G[b], H)
        g64 = vg.grid_d2(vg.build(XA, WA, SIGMA, H, LO, SHAPE), vg.build(XB, WB, SIGMA, H, LO, SHAPE), H)
        g4 = vg.grid_d2(vg.build(XA, WA, SIGMA, H, LO, SHAPE, 4.0), vg.build(XB, WB, SIGMA, H, LO, SHAPE, 4.0), H)
        errs.append(g16 / exact - 1); errs64.append(g64 / exact - 1); errs4.append(g4 / exact - 1); f16.append(g16 / g64 - 1)
    errs = np.array(errs)
    C.add("analytic grid-free distance, 12 pairs: max |grid/exact - 1| (stored grids)",
          "%.4f (range %+.4f..%+.4f)" % (np.abs(errs).max(), errs.min(), errs.max()), np.abs(errs).max() < 0.02, "< 2 %")
    C.info("  same, float64 grids before storage", "%+.4f..%+.4f" % (min(errs64), max(errs64)))
    C.info("  diagnostic: truncation at 4 sigma instead of 3 (float64)", "%+.4f..%+.4f" % (min(errs4), max(errs4)))
    C.info("  float16 storage effect on D^2", "%+.2e..%+.2e" % (min(f16), max(f16)))

    # ------------------------------------------------------------ independent slow reimplementation (3 molecules)
    worst = 0.0
    for r in rng.choice(reprow, 3, replace=False):
        X, W, _, _ = structure(jobs[r])
        Gf, Gs = vg.build(X, W, SIGMA, H, LO, SHAPE), vg.build_slow(X, W, SIGMA, H, LO, SHAPE)
        worst = max(worst, float(np.abs(Gf - Gs).max()))
    C.add("fast grid == slow loop reference, 3 molecules (max abs, mass units)", "%.2e" % worst, worst < 1e-5, "< 1e-5")

    # ------------------------------------------------------------ permuted-atom control
    X, W, _, _ = structure(jobs[reprow[0]])
    Wp = W[np.random.default_rng(12).permutation(len(W))]
    G0, Gp = vg.build(X, W, SIGMA, H, LO, SHAPE), vg.build(X, Wp, SIGMA, H, LO, SHAPE)
    d1 = float(np.abs(G0[0] - Gp[0]).max())
    dch = [float(np.abs(G0[c] - Gp[c]).max()) for c in range(1, 7)]
    C.add("permuted atom identities: channel 1 unchanged", "%.1e" % d1, d1 < 1e-12, "< 1e-12")
    C.add("permuted atom identities: channels 2-7 all change", " ".join("%.2f" % v for v in dch),
          min(dch) > 1e-3, "each > 1e-3")

    # ------------------------------------------------------------ record
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    inpil = set(M.clone_id[M.kind == "rep"])
    ncp = int((NP.clone_a.isin(inpil) & NP.clone_b.isin(inpil)).sum())
    C.info("control pairs with both members in the pilot", ncp)
    save_json(dict(sigma=SIGMA, h=H, shape=SHAPE, lo=LO.tolist(), truncation_sigma=3.0, seed=SEED, n_pilot=NPILOT,
                   grid=os.path.relpath(GRID + ".npy", ROOT), manifest=os.path.relpath(VXP("voxel_out/out/VX2a_manifest.csv"), ROOT), channels=vg.CHANNELS,
                   primary_channels=vg.CHANNELS[:vg.NPRIMARY], value="mass per voxel; L2 inner product = sum/h^3"),
              os.path.join(OUT, "VX2a_pilot.json"))
    C.write()
