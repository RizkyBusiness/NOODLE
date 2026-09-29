"""arc30 step 1 (27 Sep 2026): vector, orientation and chemistry features at N arc points per CDR3, built
straight from the D1 chunks -- lm1_extract.py + vec1_features.py + the `lin` branch of arcn3_sidechain_modes.py
in one script (functions copied, originals untouched). No state, mouse, cluster or purity information is read.

Per model, in its own coordinates:
  CDR3 (IMGT 105-117, IMGT order incl. 112B, 112A, 112) Ca arc -> N points by arc length (np.linspace(0, s_end, N));
  side-chain centroid (mean of side-chain heavy atoms, Gly -> Ca) and the 5 residue properties linearly
  interpolated on the same arc; orientation u = unit(centroid - Ca) at each point (0 if < 0.5 A).
  Point order (as lm1): Ca chain A (N), Ca chain B (N), centroids chain A (N), centroids chain B (N).
  Landmarks: Ca of IMGT A23 A41 A89 A104 A118 B23 ... B118.
Frames (as vec1; Kabsch onto the reference receptor's landmarks (config reference_receptor), row convention (X - c) @ R):
  vg  10-landmark frame, 4N points x 10 landmarks -> 40N vectors;  ug 2N unit vectors
  vc  per-chain 5-landmark frames, 2N points x 5 per chain -> 20N;  uc
  vs  stem landmarks A104 A118 B104 B118, 4N x 4 -> 16N;           us
Chemistry z-scored over all 2N points of all molecules (as v2 d3); order chain A (N) then chain B (N).
Crystal benchmark: the accepted pairs (identity >= 0.95), parsed as arcn3_sidechain_modes.parse.
Writes reference/arc30/out/A1_features_N{N}.npz, A1_property_descriptor_N{N}.npz, A1_crystal_features_N{N}.npz;
for N = 10 also the gate-G1 comparison reference/arc30/checks/G1_features_N10.csv.
All reads go through paths.src(<rel>), all writes through paths.dst(<rel>).
usage: python a30_1_features.py N
"""
import numpy as np, pandas as pd, pickle, glob, os, sys, time
import paths
EXP_NBP = paths.expected("n_benchmark_pairs")   # data-set expected value (voxel config "expected")
N = int(sys.argv[1]); OUT = "reference/arc30/out/"; CHK = "reference/arc30/checks/"
t0 = time.time()
BB = {"N", "CA", "C", "O", "OXT"}
LMPOS = [("A", p) for p in (23, 41, 89, 104, 118)] + [("B", p) for p in (23, 41, 89, 104, 118)]
AAT = {r[0]: np.array(r[1:], float) for r in np.load(paths.src("descriptors/out/D2_cdr3_atoms.npz"), allow_pickle=True)["aa_table"]}


def ikey(n, i):
    o = 0 if i == "" else ord(i.upper()) - 64
    return (n, -o) if n in (33, 61, 112) else (n, o)


def residues(atoms):
    """atoms: iterable of (ch, num, ins, resname, atom, xyz) for CDR3 105-117 -> per chain arrays in IMGT order."""
    res = {}
    for ch, num, ins, rn, at, xyz in atoms:
        r = res.setdefault((ch, num, ins), {"rn": rn, "sc": []}); xyz = np.asarray(xyz, float)
        if at == "CA": r["CA"] = xyz
        if at not in BB: r["sc"].append(xyz)
    out = {}
    for ch in "AB":
        ks = sorted([k for k in res if k[0] == ch and "CA" in res[k]], key=lambda k: ikey(k[1], k[2]))
        R = [res[k] for k in ks]
        ca = np.array([r["CA"] for r in R]); sc = np.array([np.mean(r["sc"], 0) if r["sc"] else r["CA"] for r in R])
        out[ch] = dict(ca=ca, sc=sc, pr=np.array([AAT[r["rn"].strip()] for r in R]))
    return out


def unit(v):
    L = np.linalg.norm(v, axis=1, keepdims=True)
    return np.where(L >= 0.5, v / np.maximum(L, 1e-9), 0.0)


def arcpts(r):                                    # mode `lin` of arcn3
    ca, sc, pr = r["ca"], r["sc"], r["pr"]
    if len(ca) == 1:
        s = np.zeros(1); t = np.zeros(N)
    else:
        s = np.r_[0, np.cumsum(np.linalg.norm(np.diff(ca, axis=0), axis=1))]; t = np.linspace(0, s[-1], N)
    lin = lambda V: np.stack([np.interp(t, s, V[:, k]) for k in range(V.shape[1])], 1)
    CA, SC, PR = lin(ca), lin(sc), lin(pr)
    return CA, SC, PR


def kab(X, Y):
    cx, cy = X.mean(0), Y.mean(0); Uu, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy)); d = np.sign(np.linalg.det(Uu @ Vt))
    return Uu @ np.diag([1, 1, d]) @ Vt


A_IDX, B_IDX, SI = list(range(5)), list(range(5, 10)), [3, 4, 8, 9]


def feats(res, lm, T):
    """-> dict vg vc vs ug uc us (flat) and PR (2N, 5)."""
    CAa, SCa, PRa = arcpts(res["A"]); CAb, SCb, PRb = arcpts(res["B"])
    arc = np.concatenate([CAa, CAb, SCa, SCb])                      # (4N, 3), lm1 order
    ori = unit(np.concatenate([SCa, SCb]) - np.concatenate([CAa, CAb]))   # (2N, 3), A then B
    o = {}
    Rg = kab(lm, T); o["vg"] = ((arc[:, None] - lm[None]) @ Rg).reshape(-1); o["ug"] = (ori @ Rg).reshape(-1)
    Ra = kab(lm[A_IDX], T[A_IDX]); Rb = kab(lm[B_IDX], T[B_IDX])
    PA = np.concatenate([CAa, SCa]); PB = np.concatenate([CAb, SCb])
    o["vc"] = np.concatenate([((PA[:, None] - lm[A_IDX][None]) @ Ra).reshape(-1), ((PB[:, None] - lm[B_IDX][None]) @ Rb).reshape(-1)])
    o["uc"] = np.concatenate([ori[:N] @ Ra, ori[N:] @ Rb]).reshape(-1)
    Rs = kab(lm[SI], T[SI]); o["vs"] = ((arc[:, None] - lm[SI][None]) @ Rs).reshape(-1); o["us"] = (ori @ Rs).reshape(-1)
    o["PR"] = np.concatenate([PRa, PRb])
    return o


KEYS = ("vg", "vc", "vs", "ug", "uc", "us")
# ---------------------------------------------------------------- repertoire models
MODELS = {}
for f in paths.src_glob("atoms/D1_chunk_*.pkl"):
    d = pickle.load(open(f, "rb"))
    for cid, fr, cd in zip(d["ids"], d["fr"], d["cdr3"]):
        MODELS[cid] = (residues(cd), np.array([fr[k] for k in LMPOS], float))
    del d
P3 = np.load(paths.src("descriptors/out/D3_property_descriptor.npz"), allow_pickle=True); ids = list(P3["clone_id"]); n = len(ids)
assert set(ids) == set(MODELS), "D3 order and D1 models differ"
T = MODELS[paths.reference_receptor()][1]
F = {k: [] for k in KEYS}; PRs = []
for c in ids:
    o = feats(*MODELS[c], T)
    for k in KEYS: F[k].append(o[k])
    PRs.append(o["PR"])
F = {k: np.array(v, np.float32) for k, v in F.items()}
PR = np.array(PRs, float)                                           # (n, 2N, 5)
mu, sd = PR.reshape(-1, 5).mean(0), PR.reshape(-1, 5).std(0)
CZ = ((PR - mu) / sd).astype(np.float32)
nvec = np.array([40 * N, 20 * N, 16 * N])
for k, m in zip(("vg", "vc", "vs"), nvec): assert F[k].shape == (n, 3 * m), (k, F[k].shape)
for k in ("ug", "uc", "us"): assert F[k].shape == (n, 3 * 2 * N)
np.savez_compressed(paths.dst(OUT + "A1_features_N%d.npz" % N), clone_id=np.array(ids), **F, n_vec=nvec, n_vec_keys=np.array(["vg", "vc", "vs"]), N=N)
np.savez_compressed(paths.dst(OUT + "A1_property_descriptor_N%d.npz" % N), clone_id=np.array(ids), chemz=CZ, chem_mean=mu, chem_sd=sd)
print("N=%d models %d | %s | %.0f s" % (N, n, {k: v.shape for k, v in F.items()}, time.time() - t0), flush=True)

# ---------------------------------------------------------------- crystal benchmark (as arcn3 / lm3 parse)
BEN = "tables/benchmark"
B = pd.read_csv(paths.src("descriptors/partI/out/B3i_descriptor_floor_v2.csv")); B = B[B.accepted.astype(bool) & (B.identity >= 0.95)]


def parse(path):
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


C = {"entry": [], "model": []}; C.update({"%s_%s" % (k, s): [] for k in KEYS + ("PR",) for s in "cm"})
missing = []
for _, r in B.iterrows():
    fc, rc = parse(paths.src(os.path.join(BEN, "fixed_%s_crystal.pdb" % r.entry))); fm, rm = parse(paths.src(os.path.join(BEN, r.model)))
    miss = [k for k in LMPOS if k not in fc or k not in fm]
    if miss: missing.append((r.entry, miss)); continue
    oc = feats(rc, np.array([fc[k] for k in LMPOS]), T); om = feats(rm, np.array([fm[k] for k in LMPOS]), T)
    C["entry"].append(r.entry); C["model"].append(r.model)
    for k in KEYS + ("PR",): C[k + "_c"].append(oc[k]); C[k + "_m"].append(om[k])
if missing or not len(C["entry"]) == EXP_NBP:
    print("STOP (runbook 7.2): benchmark pairs parsed %d (expected: %s), missing landmarks %s" % (len(C["entry"]), EXP_NBP, missing)); sys.exit(2)
np.savez_compressed(paths.dst(OUT + "A1_crystal_features_N%d.npz" % N), **{k: np.array(v) for k, v in C.items()}, N=N)
print("crystal pairs parsed %d (all 10 landmarks present in every crystal and model)" % len(C["entry"]), flush=True)

# ---------------------------------------------------------------- gate G1 (N = 10 only)
if N == 10:
    V1 = np.load(paths.src("reference/out/V1_features.npz")   , allow_pickle=True); assert list(V1["clone_id"]) == ids
    rows = []
    for k in KEYS:
        dif = np.abs(F[k].astype(float) - V1[k].astype(float))
        rows.append(dict(array=k, shape=str(F[k].shape), ref_shape=str(V1[k].shape), max_abs_diff=float(dif.max()), tol=1e-4, ok=bool(F[k].shape == V1[k].shape and dif.max() <= 1e-4)))
    ref = P3["chemz"].astype(float); dif = np.abs(CZ.astype(float) - ref)
    rows.append(dict(array="chemz", shape=str(CZ.shape), ref_shape=str(ref.shape), max_abs_diff=float(dif.max()), tol=2e-3, ok=bool(CZ.shape == ref.shape and dif.max() <= 2e-3)))
    G1 = pd.DataFrame(rows); G1.to_csv(paths.dst(CHK + "G1_features_N10.csv"), index=False)
    print(G1.to_string(index=False)); print("GATE G1:", "PASS" if G1.ok.all() else "FAIL")
