"""landmarks step 3 (25 Sep 2026): model-error floor of the landmark descriptor on the crystal benchmark.

For each accepted crystal matched to its TCRBuilder2+ model (identity >= 0.95; the v2_b3i_offline matching),
both structures are parsed exactly as d1 (heavy atoms, chains A/B, IMGT numbering; altloc ' ' or 'A' only),
the CDR3 arc points are built exactly as lm1, and three crystal-vs-model errors are computed:
  d_lm10   RMS difference of the 40x10 arc-to-landmark distances (no superposition)
  d_lm123  same with the 123 v2 anchors (only anchors present in both; crystal may lack some)
  d_geom   40-point RMSD after superposing the crystal on its model over the shared 123 anchors
           (the v2 prop_geom quantity, measured directly between crystal and model)
Each is then expressed as a fraction of its own clustering-scale background (1st percentile cut), so the
floors of the two descriptors are comparable. Writes landmarks/out/LM3_crystal_floor.csv
[package: definitions part only -- the parser, the landmark list and the accepted crystal-model pairs (B), which
vec3_crystal_coords.py executes; the error table LM3_crystal_floor.csv is not computed in this package]
"""
import numpy as np, pandas as pd, os
import paths
BEN = paths.src("tables/benchmark")
B = pd.read_csv(paths.src("descriptors/partI/out/B3i_descriptor_floor_v2.csv"))
B = B[B.accepted.astype(bool) & (B.identity >= 0.95)]
D2 = np.load(paths.src("descriptors/out/D2_cdr3_atoms.npz"), allow_pickle=True)
AK = [(str(k)[0], int(str(k)[1:])) for k in D2["anchor_keys"]]
LMPOS = [("A", p) for p in (23, 41, 89, 104, 118)] + [("B", p) for p in (23, 41, 89, 104, 118)]
BBA = {"N", "CA", "C", "O", "OXT"}
def ikey(n, i):
    o = 0 if i == "" else ord(i.upper()) - 64
    return (n, -o) if n in (33, 61, 112) else (n, o)
def arc(pts, vals, n=10):
    if len(pts) == 1: return np.repeat(vals[:1], n, 0)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]); t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, vals[:, k]) for k in range(vals.shape[1])], 1)
def parse(path):
    fr, res = {}, {}
    for L in open(path):
        if not L.startswith("ATOM") or L[21] not in "AB" or L[16] not in " A": continue
        el = (L[76:78].strip() or L[12:16].strip()[:1])
        if el == "H": continue
        ch, num, ins, at = L[21], int(L[22:26]), L[26].strip(), L[12:16].strip()
        xyz = np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
        if at == "CA" and ins == "" and not (27 <= num <= 38 or 56 <= num <= 65 or 105 <= num <= 117):
            fr[(ch, num)] = xyz
        elif 105 <= num <= 117:
            r = res.setdefault((ch, num, ins), {"ca": None, "sc": []})
            if at == "CA": r["ca"] = xyz
            elif at not in BBA: r["sc"].append(xyz)
    ca_b, sc_b = [], []
    for ch in "AB":
        ks = sorted([k for k in res if k[0] == ch and res[k]["ca"] is not None], key=lambda k: ikey(k[1], k[2]))
        ca = np.array([res[k]["ca"] for k in ks]); sc = np.array([np.mean(res[k]["sc"], 0) if res[k]["sc"] else res[k]["ca"] for k in ks])
        A = arc(ca, np.concatenate([ca, sc], 1)); ca_b.append(A[:, :3]); sc_b.append(A[:, 3:])
    return fr, np.concatenate(ca_b + sc_b), [len([k for k in res if k[0] == c]) for c in "AB"]
def kabsch_fit(X, Y):          # returns X moved onto Y (row vectors)
    cx, cy = X.mean(0), Y.mean(0); H = (X - cx).T @ (Y - cy); U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt)); R = U @ np.diag([1, 1, d]) @ Vt
    return lambda P: (P - cx) @ R + cy
def dvec(arcp, fr, keys):
    return np.linalg.norm(arcp[:, None] - np.array([fr[k] for k in keys])[None], axis=2).ravel()
