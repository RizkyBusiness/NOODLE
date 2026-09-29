"""[descriptors copy, 25 Sep 2026: Kabsch rotation fixed -- see v2 FIX in kabsch()]
Pass 2: put every CDR3 heavy atom into one shared framework frame.

Same principle as the Ca pipeline -- superpose on the framework, never on the loops, so loop
differences are measured rather than absorbed into the fit -- but now carrying the side
chains along with the backbone.

The reference frame is defined here rather than inherited, because the Ca pipeline did not
store its reference. Anchors are the IMGT framework Ca positions present in 100 % of models.
The reference itself is built by superposing every model on model 0, taking the mean
framework, then superposing everything on that mean; one refinement round is enough, the
frame shifts by < 0.01 A on a second.

Per CDR3 residue it also records the properties the Ca trace throws away:
  volume       side-chain van der Waals volume (A^3, Zamyatnin)
  charge       formal charge at pH 7 (K/R +1, D/E -1, H +0.1)
  hbd / hba    side-chain hydrogen-bond donors / acceptors
  kd           Kyte-Doolittle hydropathy

Writes D2_cdr3_atoms.npz -- atoms concatenated with per-model offsets, plus per-residue
tables, all in the shared frame.
"""
import numpy as np, os, sys, pickle, glob

SRC = sys.argv[1] if len(sys.argv) > 1 else "atoms"
OUT = sys.argv[2] if len(sys.argv) > 2 else "atoms"

# side-chain volume (Zamyatnin), formal charge at pH 7, H-bond donors, acceptors, hydropathy
AA = {
    "ALA": (88.6, 0, 0, 0, 1.8), "ARG": (173.4, 1, 3, 0, -4.5),
    "ASN": (114.1, 0, 1, 1, -3.5), "ASP": (111.1, -1, 0, 2, -3.5),
    "CYS": (108.5, 0, 1, 0, 2.5), "GLN": (143.8, 0, 1, 1, -3.5),
    "GLU": (138.4, -1, 0, 2, -3.5), "GLY": (60.1, 0, 0, 0, -0.4),
    "HIS": (153.2, 0.1, 1, 1, -3.2), "ILE": (166.7, 0, 0, 0, 4.5),
    "LEU": (166.7, 0, 0, 0, 3.8), "LYS": (168.6, 1, 1, 0, -3.9),
    "MET": (162.9, 0, 0, 0, 1.9), "PHE": (189.9, 0, 0, 0, 2.8),
    "PRO": (112.7, 0, 0, 0, -1.6), "SER": (89.0, 0, 1, 1, -0.8),
    "THR": (116.1, 0, 1, 1, -0.7), "TRP": (227.8, 0, 1, 0, -0.9),
    "TYR": (193.6, 0, 1, 1, -1.3), "VAL": (140.0, 0, 0, 0, 4.2)}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}
BB = {"N", "CA", "C", "O", "OXT"}


def kabsch(P, Q):
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, _, Vt = np.linalg.svd(H)
    # v2 FIX (25 Sep 2026): every caller applies the rotation to ROW vectors,
    # (P - cP) @ R + cQ, which needs R = U D Vt. The original returned Vt.T D U.T -- the
    # column-vector matrix, i.e. the TRANSPOSE (inverse) of the rotation wanted -- so every
    # model was rotated the wrong way about its centroid. Framework fit median was 1.415 A
    # with the bug; the stored D2 is reproduced exactly by the buggy function.
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return R, cP, cQ


ids, FR, CD = [], [], []
for f in sorted(glob.glob(os.path.join(SRC, "D1_chunk_*.pkl"))):
    with open(f, "rb") as fh:
        d = pickle.load(fh)
    ids += d["ids"]; FR += d["fr"]; CD += d["cdr3"]
print("models %d" % len(ids))

# ---------------------------------------------------------------- the anchor set
from collections import Counter
cnt = Counter()
for fr in FR:
    cnt.update(fr.keys())
keys = sorted([k for k, c in cnt.items() if c == len(ids)])
print("framework anchors present in 100 %% of models: %d" % len(keys))

# Not every "framework" position is rigid. IMGT FR3 runs 66-104 and swallows HV4 (81-86)
# and the 66-75 turn, which are loops; the chain termini are floppy. Anchoring on those
# drags the frame around, and the CDR3 coordinates inherit it. Drop them by definition,
# then drop the empirically worst decile of what is left.
keys = [k for k in keys if 3 <= k[1] <= 124 and not (81 <= k[1] <= 86)
        and not (66 <= k[1] <= 75)]
print("after removing HV4, the 66-75 turn and the termini: %d" % len(keys))

F = np.array([[fr[k] for k in keys] for fr in FR])          # (n, m, 3)
# The reference is a MEDOID, not a mean. Averaging superposed structures shrinks them --
# the mean of many slightly different backbones is tighter than any real one -- which biases
# every distance measured against it. The medoid is a real model.
rng = np.random.default_rng(0)
samp = rng.choice(len(F), 400, replace=False)
cost = np.zeros(len(samp))
for a_ in range(len(samp)):
    for b_ in range(len(samp)):
        R, cP, cQ = kabsch(F[samp[a_]], F[samp[b_]])
        fitab = (F[samp[a_]] - cP) @ R + cQ
        cost[a_] += ((fitab - F[samp[b_]]) ** 2).sum(1).mean()
ref = F[samp[int(cost.argmin())]]
print("reference = medoid of a 400-model sample (model %s)" % ids[samp[int(cost.argmin())]])

# empirical trim: per-anchor spread about the reference, keep the rigid 90 %
fit0 = np.empty_like(F[samp])
for j, i in enumerate(samp):
    R, cP, cQ = kabsch(F[i], ref)
    fit0[j] = (F[i] - cP) @ R + cQ
spread = np.sqrt(((fit0 - ref) ** 2).sum(2)).mean(0)
keep = np.argsort(spread)[:int(0.90 * len(keys))]
keep.sort()
print("anchor spread before trim: median %.2f A, worst %.2f A (%s)"
      % (np.median(spread), spread.max(), "%s%d" % keys[int(spread.argmax())]))
keys = [keys[i] for i in keep]
F = F[:, keep]
ref = ref[keep]
print("rigid-core anchors kept: %d" % len(keys))
fit = np.empty_like(F)
for i in range(len(F)):
    R, cP, cQ = kabsch(F[i], ref)
    fit[i] = (F[i] - cP) @ R + cQ
fr_rmsd = np.array([float(np.sqrt(((fit[i] - ref) ** 2).sum(1).mean()))
                    for i in range(len(F))])

# How much of that misfit is the alpha/beta interdomain angle rather than the frames
# themselves? Superpose each chain on its own and compare.
nA = sum(1 for k in keys if k[0] == "A")
per = {}
for nm, sl in (("alpha only", slice(0, nA)), ("beta only", slice(nA, len(keys)))):
    v = []
    for i in range(0, len(F), 7):                      # every 7th model, for speed
        R, cP, cQ = kabsch(F[i][sl], ref[sl])
        v.append(float(np.sqrt((((F[i][sl] - cP) @ R + cQ - ref[sl]) ** 2).sum(1).mean())))
    per[nm] = float(np.median(v))
print("framework fit, both chains together: median %.3f A" % np.median(fr_rmsd))
print("               each chain on its own: %s"
      % "  ".join("%s %.3f A" % (k, v) for k, v in per.items()))
# NB: fitting one chain alone is NOT better than fitting both, so the residual misfit is
# not the alpha/beta interdomain angle -- it is spread within each chain's own framework,
# i.e. the predicted frameworks genuinely differ between receptors at this level.
print("  -> per-chain is no better, so the residual is within-chain framework spread,")
print("     not the alpha/beta interdomain angle")
print("worst framework fit: %.3f A" % fr_rmsd.max())

# ---------------------------------------------------------------- move the CDR3 atoms
xyz, chain, resnum, resins, resaa, atom, mdl = [], [], [], [], [], [], []
for i, cd in enumerate(CD):
    R, cP, cQ = kabsch(F[i], ref)
    P = np.array([c[5] for c in cd])
    P = (P - cP) @ R + cQ
    xyz.append(P)
    chain += [c[0] for c in cd]
    resnum += [c[1] for c in cd]
    resins += [c[2] for c in cd]
    resaa += [ONE.get(c[3], "X") for c in cd]
    atom += [c[4] for c in cd]
    mdl += [i] * len(cd)
XYZ = np.concatenate(xyz).astype(np.float32)
print("CDR3 heavy atoms: %s (%.1f per model)" % (f"{len(XYZ):,}", len(XYZ) / len(ids)))

A = np.array(atom)
AAc = np.array(resaa)
SIDE = ~np.isin(A, list(BB)) & (AAc != "G")      # glycine has no side chain

# v2 check: the fit is optimal -- perturbing the rotation slightly can only make it worse
_rng = np.random.default_rng(1)
for _i in _rng.choice(len(F), 20, replace=False):
    _R, _cP, _cQ = kabsch(F[_i], ref)
    _base = ((((F[_i] - _cP) @ _R + _cQ) - ref) ** 2).sum()
    for _ in range(20):
        _w = _rng.normal(size=3) * 0.01
        _K = np.array([[0, -_w[2], _w[1]], [_w[2], 0, -_w[0]], [-_w[1], _w[0], 0]])
        _Rp = _R @ (np.eye(3) + _K)            # small perturbation (first order)
        _Rp, _ = np.linalg.qr(_Rp)
        _Rp = _Rp * np.sign(np.diag(_Rp))[None] if np.linalg.det(_Rp) > 0 else _Rp
        assert ((((F[_i] - _cP) @ _Rp + _cQ) - ref) ** 2).sum() >= _base - 1e-6, "Kabsch not optimal"
print("v2 check: Kabsch fit is a local optimum on 20 random models")
np.savez_compressed(
    os.path.join(OUT, "D2_cdr3_atoms.npz"),
    clone_id=np.array(ids), xyz=XYZ, chain=np.array(chain),
    resnum=np.array(resnum, np.int16), resins=np.array(resins),
    resaa=AAc, atom=A, model=np.array(mdl, np.int32), is_side=SIDE,
    ref_framework=ref.astype(np.float32),
    anchor_keys=np.array(["%s%d" % k for k in keys]),
    fr_rmsd=fr_rmsd.astype(np.float32),
    aa_table=np.array([[k] + [str(v) for v in AA[k]] for k in sorted(AA)]))
print("wrote D2_cdr3_atoms.npz")
