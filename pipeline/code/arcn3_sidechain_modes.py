"""ARC-N step 3 (27 Sep 2026): side-chain representation at the arc points -- steps P0-P2 of
the arc30 side-chain design note (27 Sep 2026). No state, mouse, cluster or purity information is read.

Modes (backbone Ca is ALWAYS linearly interpolated along the CDR3 Ca arc, IMGT 105-117, IMGT order):
  lin     side-chain centroid, orientation and properties linearly interpolated (the existing construction)
  step    each arc point takes the side chain of the residue NEAREST in arc length: that residue's centroid,
          its own orientation unit(centroid - Ca) (zero if < 0.5 A, e.g. Gly) and its 5 properties
  stepcb  as step, but orientation = unit(Cb - Ca) of that residue (Gly: virtual Cb from N, Ca, C with the
          ideal-geometry coefficients used by trRosetta / ProteinMPNN) -- a backbone-determined direction
Features exactly as vec1 (chain frames: Kabsch of each chain's 5 landmarks onto the reference receptor's (config reference_receptor); vectors from the
N Ca and N centroid points to that chain's 5 landmarks), chemistry z-scored over all points, as v2 d3.
Metrics as vec6, per mode, lam = med_bg(d_vec)/med_bg(d_chem), lam_o = w x med_bg(d_vec)/med_bg(d_ori), cut = 1st
percentile of the vec2 background pairs:
  S  share of the near-identical control pairs within the cut
  F  median crystal-model distance / cut over the benchmark pairs, chemistry = 0 (as vec6); F_chem with the
     chemistry term computed (in step modes a nearest-residue flip caused by model error shows up there)
  C  share of crystal-model pairs within the cut
Validation: N = 10 lin must reproduce V6_orientation_weight_grid.csv at w = 0, 0.5, 1.
Rule (fixed 27 Sep, before running): step replaces lin unless F(step) > F(lin) + 0.03 or S(step) < S(lin) - 0.01
(at w = 0.5 and w = 0); stepcb replaces step only if F improves by >= 0.03 without losing > 0.01 of S.
Writes reference/out/ARCN3_mode_metrics.csv, ARCN3_checks.csv; features cached to FEATDIR.
All reads go through paths.src(<rel>), all writes through paths.dst(<rel>). FEATDIR: used as given when absolute (the
pipeline runner resolves a dstdir: token), otherwise a project-relative folder created with paths.dst_dir(<rel>).
usage: python arcn3_sidechain_modes.py FEATDIR
"""
import numpy as np, pandas as pd, pickle, glob, os, sys
import paths
FD = sys.argv[1]; FD = (os.makedirs(FD, exist_ok=True) or FD) if os.path.isabs(FD) else paths.dst_dir(FD); OUTD = "reference/out/"
BB = {"N", "CA", "C", "O", "OXT"}
LMPOS = [("A", p) for p in (23, 41, 89, 104, 118)] + [("B", p) for p in (23, 41, 89, 104, 118)]
AAT = {r[0]: np.array(r[1:], float) for r in np.load(paths.src("descriptors/out/D2_cdr3_atoms.npz"), allow_pickle=True)["aa_table"]}
def ikey(n, i):
    o = 0 if i == "" else ord(i.upper()) - 64
    return (n, -o) if n in (33, 61, 112) else (n, o)
def vcb(n_, ca, c_):
    b = ca - n_; c = c_ - ca; a = np.cross(b, c); return -0.58273431 * a + 0.56802827 * b - 0.54067466 * c + ca
def residues(atoms):
    """atoms: iterable of (ch, num, ins, resname, atom, xyz) for CDR3 105-117 -> per chain arrays in IMGT order."""
    res = {}
    for ch, num, ins, rn, at, xyz in atoms:
        r = res.setdefault((ch, num, ins), {"rn": rn, "sc": []}); xyz = np.asarray(xyz, float)
        if at in ("N", "CA", "C", "CB"): r[at] = xyz
        if at not in BB: r["sc"].append(xyz)
    out = {}
    for ch in "AB":
        ks = sorted([k for k in res if k[0] == ch and "CA" in res[k]], key=lambda k: ikey(k[1], k[2]))
        R = [res[k] for k in ks]
        ca = np.array([r["CA"] for r in R]); sc = np.array([np.mean(r["sc"], 0) if r["sc"] else r["CA"] for r in R])
        cb = np.array([r["CB"] if "CB" in r else vcb(r["N"], r["CA"], r["C"]) for r in R])
        realcb = np.array(["CB" in r for r in R]); vcb_all = np.array([vcb(r["N"], r["CA"], r["C"]) for r in R])
        out[ch] = dict(ca=ca, sc=sc, cb=cb, pr=np.array([AAT[r["rn"]] for r in R]), realcb=realcb, vcb=vcb_all)
    return out
def unit(v, zero_short=True):
    L = np.linalg.norm(v, axis=1, keepdims=True)
    return np.where(L >= 0.5, v / np.maximum(L, 1e-9), 0.0) if zero_short else v / np.maximum(L, 1e-9)
def arcpts(r, N, mode):
    ca, sc, pr = r["ca"], r["sc"], r["pr"]
    if len(ca) == 1:
        s = np.zeros(1); t = np.zeros(N)
    else:
        s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(ca, axis=0), axis=1))]; t = np.linspace(0, s[-1], N)
    lin = lambda V: np.stack([np.interp(t, s, V[:, k]) for k in range(V.shape[1])], 1)
    CA = lin(ca)
    if mode == "lin":
        SC = lin(sc); PR = lin(pr); U = unit(SC - CA)
    else:
        j = np.abs(t[:, None] - s[None]).argmin(1)
        SC, PR = sc[j], pr[j]
        U = unit(sc[j] - ca[j]) if mode == "step" else unit(r["cb"][j] - ca[j], zero_short=False)
    return CA, SC, PR, U
def kab(X, Y):
    cx, cy = X.mean(0), Y.mean(0); Uu, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy)); d = np.sign(np.linalg.det(Uu @ Vt))
    return Uu @ np.diag([1, 1, d]) @ Vt
def feats(res, lm, T, N, mode):
    g, u, p = [], [], []
    for c, ch in enumerate("AB"):
        ix = slice(5 * c, 5 * c + 5); R = kab(lm[ix], T[ix]); CA, SC, PR, U = arcpts(res[ch], N, mode)
        g.append(((np.concatenate([CA, SC])[:, None] - lm[ix][None]) @ R).reshape(-1)); u.append((U @ R).reshape(-1)); p.append(PR)
    return np.concatenate(g), np.concatenate(u), np.concatenate(p)
# ---------------------------------------------------------------- inputs
MODELS = []
for f in paths.src_glob("atoms/D1_chunk_*.pkl"):
    d = pickle.load(open(f, "rb"))
    for cid, fr, cd in zip(d["ids"], d["fr"], d["cdr3"]):
        MODELS.append((cid, residues(cd), np.array([fr[k] for k in LMPOS], float)))
    del d
mid = [m[0] for m in MODELS]; T = MODELS[mid.index(paths.reference_receptor())][2]
BEN = "tables/benchmark"
B = pd.read_csv(paths.src("descriptors/partI/out/B3i_descriptor_floor_v2.csv")); B = B[B.accepted.astype(bool) & (B.identity >= 0.95)]
def parse(path):                                  # as lm3_crystal_floor.parse, keeping N/C/CB and residue names
    fr, at = {}, []
    for L in open(path):
        if not L.startswith("ATOM") or L[21] not in "AB" or L[16] not in " A": continue
        el = (L[76:78].strip() or L[12:16].strip()[:1])
        if el == "H": continue
        ch, num, ins, a = L[21], int(L[22:26]), L[26].strip(), L[12:16].strip()
        xyz = np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
        if a == "CA" and ins == "" and not (27 <= num <= 38 or 56 <= num <= 65 or 105 <= num <= 117): fr[(ch, num)] = xyz
        elif 105 <= num <= 117: at.append((ch, num, ins, L[17:20], a, xyz))
    return fr, residues(at)
CRY = []
for _, r in B.iterrows():
    fc, rc = parse(paths.src(os.path.join(BEN, "fixed_%s_crystal.pdb" % r.entry))); fm, rm = parse(paths.src(os.path.join(BEN, r.model)))
    assert all(k in fc and k in fm for k in LMPOS)
    CRY.append((r.entry, rc, np.array([fc[k] for k in LMPOS]), rm, np.array([fm[k] for k in LMPOS])))
checks = []
# check: virtual Cb vs real Cb on non-Gly residues of 500 models
dd = np.concatenate([np.linalg.norm(m[1][ch]["vcb"][m[1][ch]["realcb"]] - m[1][ch]["cb"][m[1][ch]["realcb"]], axis=1) for m in MODELS[:500] for ch in "AB"])
checks.append(dict(check="virtual_vs_real_Cb_A_median_p95", value="%.3f / %.3f" % (np.median(dd), np.percentile(dd, 95))))
# ---------------------------------------------------------------- pair sets (vec2 RNG stream)
P3 = np.load(paths.src("descriptors/out/D3_property_descriptor.npz"), allow_pickle=True); ids = list(P3["clone_id"]); n = len(ids); idx = {c: i for i, c in enumerate(ids)}
order = [mid.index(c) for c in ids]
rng = np.random.default_rng(0)
NP = pd.read_csv(paths.src("tables/B3c_near_identical_pairs.csv.gz")); pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
rows = []
for N, mode in [(10, "lin"), (30, "lin"), (30, "step"), (30, "stepcb")]:
    F = [feats(MODELS[o][1], MODELS[o][2], T, N, mode) for o in order]
    G = np.array([f[0] for f in F]); U = np.array([f[1] for f in F]); PR = np.array([f[2] for f in F])
    mu, sd = PR.reshape(-1, 5).mean(0), PR.reshape(-1, 5).std(0); CZ = (PR - mu) / sd
    np.savez(os.path.join(FD, "feat3_%s_N%d.npz" % (mode, N)), clone_id=np.array(ids), G=G.astype(np.float32), U=U.astype(np.float32), CZ=CZ.astype(np.float32), mu=mu, sd=sd)
    k, ku = 20 * N, 2 * N
    dv = lambda a, b: np.sqrt(((G[a] - G[b]) ** 2).sum(-1) / k); do = lambda a, b: np.sqrt(((U[a] - U[b]) ** 2).sum(-1) / ku)
    dc = lambda a, b: np.sqrt(((CZ[a] - CZ[b]) ** 2).sum(-1).mean(-1))
    bg = (dv(ba, bb), dc(ba, bb), do(ba, bb)); po = (dv(pa, pb), dc(pa, pb), do(pa, pb))
    lam = np.median(bg[0]) / np.median(bg[1]); lo1 = np.median(bg[0]) / np.median(bg[2])
    cg, co, cc = [], [], []
    for e, rc, lc, rm, lmm in CRY:
        gc, uc, pc = feats(rc, lc, T, N, mode); gm, um, pm = feats(rm, lmm, T, N, mode)
        cg.append(np.sqrt(((gc - gm) ** 2).sum() / k)); co.append(np.sqrt(((uc - um) ** 2).sum() / ku)); cc.append(np.sqrt((((pc - pm) / sd) ** 2).sum(-1).mean()))
    cg, co, cc = map(np.array, (cg, co, cc))
    if mode != "lin": checks.append(dict(check="%s_N%d_property_values_on_residue_lattice" % (mode, N), value=bool(np.isin(PR[:, :, 0].round(3), np.round([v[0] for v in AAT.values()], 3)).all())))
    for w in (0.0, 0.5, 1.0):
        lo = w * lo1; comb = lambda t: np.sqrt(t[0] ** 2 + (lam * t[1]) ** 2 + (lo * t[2]) ** 2)
        b, p = comb(bg), comb(po); cut = np.percentile(b, 1)
        cr = np.sqrt(cg ** 2 + (lo * co) ** 2); crc = np.sqrt(cg ** 2 + (lam * cc) ** 2 + (lo * co) ** 2)
        rows.append(dict(N=N, mode=mode, w=w, lam=lam, lam_o=lo, cut=cut, S=(p <= cut).mean(), F=np.median(cr) / cut, F_chem=np.median(crc) / cut,
                         C=(cr <= cut).mean(), crystal_chem_nonzero=int((cc > 1e-9).sum()), orient_share=np.median((lo * bg[2]) ** 2 / b ** 2)))
    print(N, mode, "done", flush=True)
Mt = pd.DataFrame(rows); Mt.to_csv(paths.dst(OUTD + "ARCN3_mode_metrics.csv"), index=False)
V6 = pd.read_csv(paths.src(OUTD + "V6_orientation_weight_grid.csv")).set_index("w")
for w in (0.0, 0.5, 1.0):
    r = Mt[(Mt.N == 10) & (Mt["mode"] == "lin") & (Mt.w == w)].iloc[0]
    checks.append(dict(check="N10_lin_reproduces_V6_w%.1f" % w, value="cut %.4f/%.4f S %.4f/%.4f F %.4f/%.4f" % (r.cut, V6.loc[w, "cut"], r.S, V6.loc[w, "S_sensitivity"], r.F, V6.loc[w, "F_crystal_over_cut"]),
                       ok=bool(abs(r.cut - V6.loc[w, "cut"]) < 2e-3 and abs(r.F - V6.loc[w, "F_crystal_over_cut"]) < 2e-3 and abs(r.S - V6.loc[w, "S_sensitivity"]) < 1e-4)))
Ck = pd.DataFrame(checks); Ck.to_csv(paths.dst(OUTD + "ARCN3_checks.csv"), index=False)
pd.set_option("display.width", 250); print(Ck.to_string(index=False)); print(Mt.round(4).to_string(index=False))
