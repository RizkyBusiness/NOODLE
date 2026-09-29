"""VX5b (the voxel pipeline): frame decomposition diagnostic (DESIGN amendment A9.1; decision rules A9.2). Diagnostic and exploratory:
no state data are read, VX7 is not run, and nothing here changes the VX5 verdict.

Frames (Kabsch 1976, row convention; then the fixed H15 axes as production):
  F0 whole molecule, 10 landmarks (production)      F1 per chain, 5 landmarks each (vc_ori frames)
  F2 per chain, that chain's framework anchors      F3 whole molecule, all 123 anchors
Per-chain frames put alpha and beta atoms into two separate boxes (concatenated for the distance). Each chain box is the
chain's envelope over the VX1 500-molecule sample (0.1-99.9 pct, + 7 A per A10, whole voxels), frozen
before any distance is computed. Atom sets: whole molecule (as VX2c) and loops only (as A8). sigma 2.0 A, voxel 1.0 A,
channels 1-7, typing A2/A3, float64 Gram (A7), shared-atom crystal comparison (A2.1/A3.1), the accepted pairs.
Structures: the 2,000 VX2a pilot molecules (background: all their pairs except those that are control pairs), the
control-pair receptors (those not in the pilot added; control statistics only), the crystal/model pairs on shared atoms.
Hinge: crystal onto its own model per chain on the shared framework anchors; angle of R_alpha^T R_beta.
Bracketing (report only): F0 whole molecule at sigma 2.8 and 3.5 A, run last.
usage: python pipeline/code/vx5b_frame_diagnostic.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, time
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NANC_A = expected("n_anchors_alpha")
EXP_NANC_B = expected("n_anchors_beta")
import vxgrid as vg
import vx2a_pilot as PA
from vx6L_loops import loop_of

H, SIG = 1.0, 2.0
NP7 = vg.NPRIMARY
L = np.load(LM1, allow_pickle=True)
IDS = [str(c) for c in L["clone_id"]]; IDX = {c: i for i, c in enumerate(IDS)}
LMK, ANC = L["lm"].astype(float), L["anch"].astype(float)
AKEYS = [str(k) for k in L["anch_keys"]]
AK_T = [(k[0], int(k[1:])) for k in AKEYS]
RI = IDS.index(REF_ID)
T_LM, T_AN = LMK[RI], ANC[RI]
O, AX = PA.O, PA.AX
VX1LO, VX1SH = PA.LO, np.array(PA.SHAPE)
FRAMES = ["F0", "F1", "F2", "F3"]
PERCHAIN = {"F1", "F2"}
TMPG = os.path.join(TMP, "VX5b_grid.f16.npy")


def to_box(R, t):
    """fit-space (X @ R + t) -> box-space (Rb, tb) with the fixed H15 axes"""
    return R @ AX.T, (t - O) @ AX.T


def fit(X, Y):
    R, cx, cy = kabsch(X, Y)
    return R, cy - cx @ R


def placements(frame, lm, an, an_ok):
    """lm (10,3); an (123,3) with an_ok mask of anchors present. Returns {'AB': (Rb,tb)} or {'A':..., 'B':...}."""
    ia = np.array([k[0] == "A" for k in AKEYS])
    if frame == "F0":
        return {"AB": to_box(*fit(lm, T_LM))}
    if frame == "F3":
        return {"AB": to_box(*fit(an[an_ok], T_AN[an_ok]))}
    if frame == "F1":
        return {"A": to_box(*fit(lm[:5], T_LM[:5])), "B": to_box(*fit(lm[5:], T_LM[5:]))}
    return {"A": to_box(*fit(an[an_ok & ia], T_AN[an_ok & ia])), "B": to_box(*fit(an[an_ok & ~ia], T_AN[an_ok & ~ia]))}


def bench_points(path, lm3parse):
    fr = lm3parse(path)[0]
    lm = np.array([fr[k] for k in PA.lm3_LMPOS], float)
    ok = np.array([k in fr for k in AK_T])
    an = np.array([fr[k] if k in fr else (0.0, 0.0, 0.0) for k in AK_T], float)
    return lm, an, ok


def structure(spec):
    """spec = (path, place, keep, atomset). Returns per-box list of (X_box, W) and bookkeeping."""
    path, place, keep, atomset = spec
    xyz, meta = parse_heavy(path)
    W, bad = vg.type_atoms(meta)
    m = np.ones(len(meta), bool)
    if atomset == "loops":
        m &= np.array([loop_of(a[1]) is not None for a in meta])
    if keep is not None:
        m &= np.array([tuple(a) in keep for a in meta])
    ch = np.array([a[0] for a in meta])
    out = []
    for key, (Rb, tb) in place.items():
        s = m & np.isin(ch, list(key))
        out.append((key, xyz[s] @ Rb + tb, W[s]))
    return out, int(m.sum()), len(bad)


BOXES = None      # set in main, frozen before any distance; passed to workers explicitly


def box_of(frame, key):
    return BOXES[frame][key]


def build(args):
    spec, boxes, sigma = args               # boxes passed explicitly: workers are spawned, not forked
    parts, na, nb = structure(spec)
    vecs, outside = [], 0
    for key, X, W in parts:
        lo, shp = boxes[key]
        g = vg.build(X, W, sigma, H, lo, shp)[:NP7]
        vecs.append(g.reshape(-1))
        outside += int((~((X >= lo) & (X < lo + np.array(shp) * H)).all(1)).sum())
    v = np.concatenate(vecs)
    return v.astype(np.float16), na, outside, nb


def occ_mass(v, frame):
    s, tot = 0, 0.0
    for key in (["AB"] if frame not in PERCHAIN else ["A", "B"]):
        nv = int(np.prod(box_of(frame, key)[1]))
        tot += float(v[s:s + nv].astype(np.float64).sum()); s += NP7 * nv
    return tot


def analytic(specA, specB, sigma):
    pa, _, _ = structure(specA); pb, _, _ = structure(specB)
    return sum(vg.analytic_d2(XA, WA, XB, WB, sigma) for (_, XA, WA), (_, XB, WB) in zip(pa, pb))


def spearman_ci(x, y, rng, nb=2000):
    r = spearmanr(x, y).correlation
    bs = []
    for _ in range(nb):
        i = rng.integers(0, len(x), len(x))
        bs.append(spearmanr(x[i], y[i]).correlation)
    bs = np.array(bs)
    return float(r), float(np.nanpercentile(bs, 2.5)), float(np.nanpercentile(bs, 97.5))


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VX5b", stop_on_fail=True)       # stop at the first failed check (A10)
    C.info("status", "diagnostic and exploratory (A9); no state data; VX7 not run; VX5 verdict unchanged")
    lm3 = load_lm3(); PA.lm3_LMPOS = lm3["LMPOS"]
    Bp = lm3["B"].sort_values("entry").reset_index(drop=True)
    FRB = np.load(os.path.join(OUT, "VX1_frames.npz"), allow_pickle=True)
    C.add("anchors per chain (A9.1: %s alpha / %s beta)" % (EXP_NANC_A, EXP_NANC_B), "%d / %d" % (sum(k[0] == "A" for k in AKEYS),
          sum(k[0] == "B" for k in AKEYS)), sum(k[0] == "A" for k in AKEYS) == EXP_NANC_A and sum(k[0] == "B" for k in AKEYS) == EXP_NANC_B,
          "%s / %s" % (EXP_NANC_A, EXP_NANC_B))

    # ------------------------------------------------------------ fit checks
    allok = np.ones(123, bool)
    for f in ("F1", "F2", "F3"):
        P = placements(f, T_LM, T_AN, allok)
        e = max(max(float(np.abs(Rb - AX.T).max()), float(np.abs(tb + O @ AX.T).max())) for Rb, tb in P.values())
        C.add("%s reference maps to itself (R = I, t = 0 before the H15 axes)" % f, "%.2e" % e, e < 1e-6, "< 1e-6")
    rng = np.random.default_rng(71)
    for f in ("F0", "F1", "F2", "F3"):
        worst = 0.0
        for i in rng.choice(len(IDS), 5, replace=False):
            xyz, meta = parse_heavy(VXP("structures/%s.pdb") % IDS[i]); ch = np.array([a[0] for a in meta])
            q = rng.normal(size=4); q /= np.linalg.norm(q); a, b, c, d = q
            Q = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                          [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
            sh = rng.uniform(-50, 50, 3)
            P1 = placements(f, LMK[i], ANC[i], allok); P2 = placements(f, LMK[i] @ Q + sh, ANC[i] @ Q + sh, allok)
            for key in P1:
                s = np.isin(ch, list(key))
                worst = max(worst, float(np.abs((xyz[s] @ P1[key][0] + P1[key][1]) - ((xyz[s] @ Q + sh) @ P2[key][0] + P2[key][1])).max()))
        C.add("%s fit invariant to input rigid transform (5 models)" % f, "%.2e A" % worst, worst < 1e-4, "< 1e-4 A")

    # ------------------------------------------------------------ per-chain boxes, frozen before any distance
    samp = FRB["box_sample"]
    MARGIN = 7.0                                # A10.1: 3 sigma + 1 A, every VX5b box, no cropping
    BOXES = {"F0_vx1box": {"AB": (VX1LO, tuple(VX1SH))}}       # A10.2: reproduction check only
    boxrec = {}
    parsed = {}
    for i in samp:
        xyz, meta = parse_heavy(VXP("structures/%s.pdb") % IDS[i]); parsed[i] = (xyz, np.array([a[0] for a in meta]))
    for f in FRAMES:
        BOXES[f] = {}
        keys = ["A", "B"] if f in PERCHAIN else ["AB"]
        pts = {k: [] for k in keys}
        for i in samp:
            xyz, ch = parsed[i]
            P = placements(f, LMK[i], ANC[i], allok)
            for key in keys:
                s_ = np.isin(ch, list(key))
                pts[key].append(xyz[s_] @ P[key][0] + P[key][1])
        for key in keys:
            Pk = np.concatenate(pts[key])
            lo = np.floor((np.percentile(Pk, 0.1, 0) - MARGIN) / H) * H
            hi = np.ceil((np.percentile(Pk, 99.9, 0) + MARGIN) / H) * H
            shp = tuple(int(round(x)) for x in (hi - lo) / H)
            BOXES[f][key] = (lo, shp)
            boxrec["%s_%s" % (f, key)] = dict(lo=lo.tolist(), hi=hi.tolist(), shape=shp, margin_A=MARGIN)
    save_json(dict(frozen=True, written=pd.Timestamp.now().isoformat(), sample="VX1 box_sample (500, seed 0)",
                   boxes=boxrec, reproduction_box="F0_vx1box = VX1 box (4 A margin), A10.2 only"),
              os.path.join(OUT, "VX5b_boxes.json"))
    t_boxes = time.time()
    for k, v in boxrec.items():
        C.info("box %s lo / shape" % k, "%s / %s" % (v["lo"], v["shape"]))

    # ------------------------------------------------------------ structures and placements for every frame
    MP = pd.read_csv(os.path.join(OUT, "VX2a_manifest.csv"))
    pil = list(MP.clone_id[MP.kind == "rep"])
    NPp = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    ctrl_ids = sorted(set(NPp.clone_a) | set(NPp.clone_b))
    extra = [c for c in ctrl_ids if c not in set(pil)]
    rep_ids = pil + extra
    C.info("rows: pilot molecules / added control receptors / crystal-model shared pairs",
           "%d / %d / %d" % (len(pil), len(extra), len(Bp)))
    place = {f: {} for f in FRAMES}
    for c in rep_ids:
        i = IDX[c]
        for f in FRAMES:
            place[f][c] = placements(f, LMK[i], ANC[i], allok) if f != "F0" else {"AB": (FRB["R_box"][i], FRB["t_box"][i])}
    bench = []
    for _, r in Bp.iterrows():
        pc, pm = os.path.join(PA.BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(PA.BEN, r.model)
        keep = {tuple(a) for a in parse_heavy(pc)[1]} & {tuple(a) for a in parse_heavy(pm)[1]}
        pts = {p: bench_points(p, lm3["parse"]) for p in (pc, pm)}
        pl = {f: {} for f in FRAMES}
        for p in (pc, pm):
            lm_, an_, ok_ = pts[p]
            for f in FRAMES:
                pl[f][p] = placements(f, lm_, an_, ok_) if f != "F0" else {"AB": PA.bench_frame(p, lm3["parse"], lm3["LMPOS"])}
        # hinge: crystal onto its own model per chain, shared framework anchors
        okb = pts[pc][2] & pts[pm][2]; ia = np.array([k[0] == "A" for k in AKEYS])
        Ra, _ = fit(pts[pc][1][okb & ia], pts[pm][1][okb & ia]); Rb_, _ = fit(pts[pc][1][okb & ~ia], pts[pm][1][okb & ~ia])
        hinge = float(np.degrees(np.arccos(np.clip((np.trace(Ra.T @ Rb_) - 1) / 2, -1, 1))))
        bench.append(dict(entry=r.entry, model=r.model, pc=pc, pm=pm, keep=keep, place=pl, hinge_deg=hinge,
                          anchors_shared=int(okb.sum())))
    nrep = len(rep_ids)
    pos = {c: j for j, c in enumerate(rep_ids)}
    ipil = np.arange(len(pil))
    cpa = np.array([pos[a] for a in NPp.clone_a]); cpb = np.array([pos[b] for b in NPp.clone_b])
    iu = np.triu_indices(len(pil), 1)
    ctrl_in_pil = {(min(a, b), max(a, b)) for a, b in zip(cpa, cpb) if a < len(pil) and b < len(pil)}
    bgmask = np.ones(len(iu[0]), bool)
    for a, b in ctrl_in_pil:
        bgmask[np.where((iu[0] == a) & (iu[1] == b))[0]] = False
    C.info("background pairs (pilot pairs minus control pairs)", "%d (removed %d)" % (int(bgmask.sum()), len(ctrl_in_pil)))

    # ------------------------------------------------------------ run one configuration
    rows, pairs_out = [], []
    NVOX = {f: sum(int(np.prod(s)) for _, s in BOXES[f].values()) for f in BOXES}

    def run(frame, atomset, sigma, label, report_only=False, boxkey=None):
        bk = boxkey or frame
        ts = time.time()
        specs = [(VXP("structures/%s.pdb") % c, place[frame][c], None, atomset) for c in rep_ids]
        for b in bench:
            specs += [(b["pc"], b["place"][frame][b["pc"]], b["keep"], atomset),
                      (b["pm"], b["place"][frame][b["pm"]], b["keep"], atomset)]
        N, F = len(specs), NP7 * NVOX[bk]
        G = np.lib.format.open_memmap(TMPG, mode="w+", dtype=np.float16, shape=(N, F))
        na, outside, nbad = np.zeros(N, int), 0, 0
        with Pool(8) as pool:
            for j, (v, n_, o_, b_) in enumerate(pool.imap(build, [(s, BOXES[bk], sigma) for s in specs], chunksize=4)):
                G[j] = v; na[j] = n_; outside += o_; nbad += b_
        G.flush()
        C.add("%s typing complete" % label, nbad, nbad == 0, "== 0")
        if frame in PERCHAIN and atomset == "whole" and not report_only:
            C.add("%s per-chain boxes: atoms outside" % label, "%d / %d = %.4f %%" % (outside, int(na.sum()), 100 * outside / na.sum()),
                  outside / na.sum() < 1e-4, "< 0.01 %")
        else:
            C.info("%s atoms outside their box" % label, "%d / %d" % (outside, int(na.sum())))
        r2 = np.random.default_rng(72)
        mrow = r2.choice(len(pil), 200, replace=False)
        rel = max(abs(occ_mass(G[j], bk) / na[j] - 1) for j in mrow)
        C.add("%s mass conservation, 200 molecules" % label, "%.4f %%" % (100 * rel), rel < 5e-3, "< 0.5 %")
        errs = []
        for a, b in r2.choice(len(pil), (12, 2), replace=False):
            g = float(((G[a].astype(np.float64) - G[b].astype(np.float64)) ** 2).sum()) / H ** 3
            errs.append(g / analytic(specs[a], specs[b], sigma) - 1)
        C.add("%s analytic grid-free distance, 12 pairs" % label, "%.4f (%+.4f..%+.4f)" % (np.abs(errs).max(), min(errs), max(errs)),
              np.abs(errs).max() < 0.02, "< 2 %")
        if not rows:
            bt = os.path.getmtime(os.path.join(OUT, "VX5b_boxes.json"))
            C.add("per-chain boxes frozen (written) before any distance is computed", "boxes %.0f s before first Gram" % (time.time() - bt),
                  bt < time.time(), "written first")
        Gm = np.zeros((N, N))
        for s in range(0, F, 20000):
            X = np.asarray(G[:, s:min(s + 20000, F)], dtype=np.float64)
            Gm += X @ X.T
        nr = np.diag(Gm).copy()
        D = np.sqrt(np.clip((nr[:, None] + nr[None] - 2 * Gm) / H ** 3, 0, None)); np.fill_diagonal(D, 0.0)
        del Gm
        bg = D[:len(pil), :len(pil)][iu][bgmask]
        cut = float(np.percentile(bg, 1)); bgs = np.sort(bg)
        cd = D[cpa, cpb]; cpct = 100 * np.searchsorted(bgs, cd, side="right") / len(bgs)
        dx = np.array([D[nrep + 2 * k, nrep + 2 * k + 1] for k in range(len(bench))])
        o = dict(frame=bk, atoms=atomset, sigma=sigma, report_only=report_only, cut=cut, bg_median=float(np.median(bg)),
                 crystal_ratio_median=float(np.median(dx) / cut), crystal_ratio_mean=float(np.mean(dx) / cut),
                 crystal_within_cut=int((dx <= cut).sum()), control_within_cut=float((cd <= cut).mean()),
                 control_median_bg_pct=float(np.median(cpct)), control_median_over_cut=float(np.median(cd) / cut),
                 minutes=round((time.time() - ts) / 60, 2))
        rows.append(o)
        for b, d_ in zip(bench, dx):
            pairs_out.append(dict(frame=bk, atoms=atomset, sigma=sigma, entry=b["entry"], model=b["model"],
                                  d=d_, over_cut=d_ / cut, hinge_deg=b["hinge_deg"]))
        print("%s: cut %.4f ratio %.4f (%.1f min)" % (label, cut, o["crystal_ratio_median"], o["minutes"]), flush=True)
        return dx

    dx = run("F0", "whole", SIG, "F0 whole, VX1 box (reproduction check only)", report_only=True, boxkey="F0_vx1box")
    V5 = pd.read_csv(os.path.join(OUT, "VX5_crystal_floor.csv")).set_index("entry").loc[[b["entry"] for b in bench]]
    e = float(np.max(np.abs(dx - V5.d_shared.values) / V5.d_shared.values))
    C.add("F0 whole in the VX1 box reproduces VX5_crystal_floor.csv distances (%s pairs)" % EXP_NBP, "%.2e rel" % e, e < 1e-6, "< 1e-6")
    for f in FRAMES:
        for a in ("whole", "loops"):
            run(f, a, SIG, "%s %s" % (f, a))
    S = pd.DataFrame(rows); PR = pd.DataFrame(pairs_out)
    S.to_csv(os.path.join(OUT, "VX5b_configs.csv"), index=False)
    PR.to_csv(os.path.join(OUT, "VX5b_crystal_pairs.csv"), index=False)

    # ------------------------------------------------------------ hinge
    HG = pd.DataFrame([dict(entry=b["entry"], model=b["model"], hinge_deg=b["hinge_deg"], anchors_shared=b["anchors_shared"])
                       for b in bench])
    HG.to_csv(os.path.join(OUT, "VX5b_hinge.csv"), index=False)
    rep_rot = FRB["vavb_rot_deg"][FRB["is_molecule"]]
    q = lambda x: "median %.2f, IQR %.2f-%.2f, p95 %.2f, max %.2f" % (np.median(x), np.percentile(x, 25), np.percentile(x, 75),
                                                                     np.percentile(x, 95), np.max(x))
    C.info("hinge: crystal vs own model (%s pairs), deg" % EXP_NBP, q(HG.hinge_deg.values))
    C.info("hinge: repertoire Va/Vb rotation vs reference (VX1, {:,}), deg".format(EXP_NMOL), q(rep_rot))
    C.info("hinge: shared framework anchors per pair, min", int(HG.anchors_shared.min()))
    hin = {}
    for f in ("F0", "F1", "F2"):
        p = PR[(PR.frame == f) & (PR.atoms == "whole") & (PR.sigma == SIG)].set_index("entry").loc[HG.entry]   # F0 = widened box
        hin[f] = spearman_ci(HG.hinge_deg.values, p.over_cut.values, np.random.default_rng(0))
        C.info("Spearman(hinge, %s whole error over cut) [95 %% bootstrap CI]" % f, "%.3f [%.3f, %.3f]" % hin[f])

    # ------------------------------------------------------------ A9.2 decision quantities (reported, not interpreted)
    S8 = S[~S.report_only]
    g = S8[(S8.atoms == "whole") & (S8.sigma == SIG)].set_index("frame").crystal_ratio_median
    best_pc = "F1" if g["F1"] <= g["F2"] else "F2"
    rule2 = bool(min(g["F1"], g["F2"]) <= g["F0"] - 0.05)
    choice = "F1" if abs(g["F1"] - g["F2"]) <= 0.02 else best_pc
    C.info("A9.2 rule 2: F0 / F1 / F2 / F3 whole-molecule error over cut", "%.4f / %.4f / %.4f / %.4f" % tuple(g[f] for f in FRAMES))
    C.info("A9.2 rule 2: better per-chain frame vs F0 - 0.05", "%s %.4f vs %.4f -> %s" % (best_pc, g[best_pc], g["F0"] - 0.05,
                                                                                    "HOLDS" if rule2 else "FAILS"))
    C.info("A9.2 rule 2 context: better per-chain frame below 0.60", bool(g[best_pc] < 0.60))
    C.info("A9.2 rule 3: F0 Spearman positive and CI excludes 0", bool(hin["F0"][0] > 0 and hin["F0"][1] > 0))
    C.info("A9.2 rule 3: weakens under F1 / F2", "%s / %s" % (abs(hin["F1"][0]) < abs(hin["F0"][0]), abs(hin["F2"][0]) < abs(hin["F0"][0])))
    C.info("A9.2 rule 4: per-chain frame by rule (if v2 is proposed)", "%s (|F1-F2| = %.4f)" % (choice, abs(g["F1"] - g["F2"])))
    save_json(dict(ratios_whole={f: float(g[f]) for f in FRAMES},
                   ratios_loops={f: float(v) for f, v in S8[S8.atoms == "loops"].set_index("frame").crystal_ratio_median.items()},
                   rule2_holds=rule2, better_per_chain=best_pc, v2_frame_by_rule=choice,
                   spearman_hinge={f: dict(rho=v[0], ci=[v[1], v[2]]) for f, v in hin.items()},
                   hinge_crystal_vs_model_median=float(HG.hinge_deg.median()), hinge_repertoire_median=float(np.median(rep_rot))),
              os.path.join(OUT, "VX5b_summary.json"))

    # ------------------------------------------------------------ bracketing, report only, last
    for s_ in (2.8, 3.5):
        run("F0", "whole", s_, "bracket F0 whole sigma %.1f (report only)" % s_, report_only=True)
    S = pd.DataFrame(rows); S.to_csv(os.path.join(OUT, "VX5b_configs.csv"), index=False)
    pd.DataFrame(pairs_out).to_csv(os.path.join(OUT, "VX5b_crystal_pairs.csv"), index=False)
    br = S[(S.frame == "F0") & (S.atoms == "whole")].sort_values("sigma")   # widened F0 box at 2.0 / 2.8 / 3.5
    C.info("bracketing (report only): F0 whole error over cut by sigma", "; ".join("%.1f: %.4f" % (a, b) for a, b in
                                                                                  zip(br.sigma, br.crystal_ratio_median)))
    os.remove(TMPG)
    C.info("total time", "%.1f min" % ((time.time() - t0) / 60))
    C.write()
