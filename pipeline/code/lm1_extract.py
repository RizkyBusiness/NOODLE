"""landmarks step 1 (25 Sep 2026): per-model internal coordinates for the framework-landmark descriptor.

For every model, in the model's OWN coordinates (no superposition of any kind):
  ARC   (40, 3)  CDR3a then CDR3b, 10 arc-length points each: first the 20 Ca points, then the 20 side-chain
                 centroid points -- exactly the v2 d3 construction (IMGT sequence order incl. 112B,112A,112;
                 centroid = mean of side-chain heavy atoms, glycine -> its Ca; linear interpolation along the Ca arc)
  LM    (10, 3)  Ca of the IMGT-conserved framework positions 23, 41, 89, 104, 118 on chain A then chain B
  ANCH  (123, 3) Ca of the v2 D2 123 rigid-core anchors (sensitivity variant)
Source: atoms/D1_chunk_*.pkl (framework Ca + all CDR3 heavy atoms, straight from the PDBs).
Writes landmarks/out/LM1_internal_coords.npz
"""
import numpy as np, pickle, glob, os
import paths
BB = {"N", "CA", "C", "O", "OXT"}
LMPOS = [("A", p) for p in (23, 41, 89, 104, 118)] + [("B", p) for p in (23, 41, 89, 104, 118)]
D2 = np.load(paths.src("descriptors/out/D2_cdr3_atoms.npz"), allow_pickle=True)
AK = [(str(k)[0], int(str(k)[1:])) for k in D2["anchor_keys"]]
def ikey(n, i):
    o = 0 if i == "" else ord(i.upper()) - 64
    return (n, -o) if n in (33, 61, 112) else (n, o)
def arc(pts, vals, n=10):
    if len(pts) == 1: return np.repeat(vals[:1], n, 0)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]); t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, vals[:, k]) for k in range(vals.shape[1])], 1)
ids, ARC, LM, AN = [], [], [], []
for f in paths.src_glob("atoms/D1_chunk_*.pkl"):
    d = pickle.load(open(f, "rb"))
    for cid, fr, cd in zip(d["ids"], d["fr"], d["cdr3"]):
        res = {}
        for ch, num, ins, rn, at, xyz in cd:
            r = res.setdefault((ch, num, ins), {"ca": None, "sc": []})
            if at == "CA": r["ca"] = xyz
            elif at not in BB: r["sc"].append(xyz)
        blocks_ca, blocks_sc = [], []
        for ch in "AB":
            ks = sorted([k for k in res if k[0] == ch], key=lambda k: ikey(k[1], k[2]))
            ca = np.array([res[k]["ca"] for k in ks], float)
            sc = np.array([np.mean(res[k]["sc"], 0) if res[k]["sc"] else res[k]["ca"] for k in ks], float)
            A = arc(ca, np.concatenate([ca, sc], 1)); blocks_ca.append(A[:, :3]); blocks_sc.append(A[:, 3:])
        ARC.append(np.concatenate(blocks_ca + blocks_sc)); LM.append([fr[k] for k in LMPOS]); AN.append([fr[k] for k in AK]); ids.append(cid)
np.savez_compressed(paths.dst("landmarks/out/LM1_internal_coords.npz"), clone_id=np.array(ids), arc=np.array(ARC, np.float32),
                    lm=np.array(LM, np.float32), anch=np.array(AN, np.float32),
                    lm_keys=np.array(["%s%d" % k for k in LMPOS]), anch_keys=D2["anchor_keys"])
print("models", len(ids), "arc", np.array(ARC).shape)
# check: internal geometry of the arc points equals the v2 D3 descriptor's (which is the same points in the D2 frame)
P = np.load(paths.src("descriptors/out/D3_property_descriptor.npz"), allow_pickle=True)
assert list(P["clone_id"]) == ids
G = np.concatenate([P["geom"][:, :, :3], P["geom"][:, :, 3:]], 1).astype(float)
rng = np.random.default_rng(0); worst = 0.0
for i in rng.choice(len(ids), 300, replace=False):
    a = np.array(ARC[i]); b = G[i]
    Da = np.linalg.norm(a[:, None] - a[None], axis=2); Db = np.linalg.norm(b[:, None] - b[None], axis=2)
    worst = max(worst, float(np.abs(Da - Db).max()))
print("check: max |internal distance difference| vs v2 D3 points over 300 models: %.2e A" % worst)
assert worst < 1e-3
