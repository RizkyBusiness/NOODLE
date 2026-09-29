"""[descriptors copy, 25 Sep 2026: IMGT insertion-order fix, see v2 FIX below]
Arm D: the CDR3 descriptor with side chains and their chemistry attached.

The Ca arms describe a CDR3 by where its backbone runs. This describes it by where the
backbone runs, where each side chain points, how much space that side chain fills, and what
it can do chemically. Per CDR3 residue:

    Ca              backbone position, in the shared frame            3
    side centroid   mean of the side-chain heavy atoms                3
    volume          side-chain van der Waals volume, A^3              1
    charge          formal charge at pH 7                             1
    hbd / hba       side-chain H-bond donors / acceptors              2
    kd              Kyte-Doolittle hydropathy                         1

Each chain's CDR3 is resampled to 10 positions along its Ca arc, so receptors of different
CDR3 lengths are directly comparable -- the length-class restriction the Ca arms need does
not apply here. Properties are interpolated along the same arc as the coordinates.

Distance between two receptors has two parts, kept separate so they can be read apart:

    d_geom   RMSD over the 20 Ca points and 20 side-centroid points (A)
    d_chem   RMS difference over the standardised per-position properties

and are combined as d = sqrt(d_geom^2 + (lam * d_chem)^2), with lam set so that the two
terms have the same median in the background distribution. That is a choice, and it is the
reason both terms are also reported on their own.

Writes D3_property_descriptor.npz
"""
import numpy as np, os, sys

SRC = sys.argv[1] if len(sys.argv) > 1 else "atoms"
OUT = sys.argv[2] if len(sys.argv) > 2 else "atoms"
NRES = 10

z = np.load(os.path.join(SRC, "D2_cdr3_atoms.npz"), allow_pickle=True)
ids = z["clone_id"]
xyz, chain, resnum, resins = z["xyz"], z["chain"], z["resnum"], z["resins"]
resaa, atom, model, is_side = z["resaa"], z["atom"], z["model"], z["is_side"]
AA = {r[0]: tuple(float(x) for x in r[1:]) for r in z["aa_table"]}
ONE = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN",
       "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
       "M": "MET", "F": "PHE", "P": "PRO", "S": "SER", "T": "THR", "W": "TRP",
       "Y": "TYR", "V": "VAL"}
PROPS = ["volume", "charge", "hbd", "hba", "kd"]
print("models %d | CDR3 heavy atoms %s" % (len(ids), f"{len(xyz):,}"))

# ---------------------------------------------------------------- group atoms by residue
# v2 FIX (25 Sep 2026): residues must follow IMGT *sequence* order, not a plain sort. IMGT
# insertions run upward after 111/32/60 (111, 111A, 111B) but DOWNWARD before 112/33/61
# (112B, 112A, 112). The original lexsort on (resnum, resins) put 112, 112A, 112B and made the
# arc zig-zag through the CDR3 apex in every model with a 112 insertion.
_ins = np.array([0 if s == "" else ord(s.upper()) - 64 for s in resins], np.int16)
ins_rank = np.where(np.isin(resnum, (33, 61, 112)), -_ins, _ins)
order = np.lexsort((atom, ins_rank, resnum, chain, model))
m, c, n, i_ = model[order], chain[order], resnum[order], resins[order]
newres = np.ones(len(order), bool)
newres[1:] = (m[1:] != m[:-1]) | (c[1:] != c[:-1]) | (n[1:] != n[:-1]) | (i_[1:] != i_[:-1])
rid = np.cumsum(newres) - 1
nres = rid[-1] + 1
X = xyz[order]
isCA = atom[order] == "CA"
side = is_side[order]
print("CDR3 residues: %s (%.1f per model)" % (f"{nres:,}", nres / len(ids)))

ca = np.zeros((nres, 3), np.float32)
np.add.at(ca, rid[isCA], X[isCA])
cen = np.zeros((nres, 3), np.float64)
cnt = np.zeros(nres)
np.add.at(cen, rid[side], X[side])
np.add.at(cnt, rid[side], 1.0)
has = cnt > 0
cen[has] /= cnt[has, None]
cen[~has] = ca[~has]                      # glycine: the side centroid is the Ca itself
r_model = np.zeros(nres, np.int32); r_model[rid] = m
r_chain = np.empty(nres, "<U1"); r_chain[rid] = c
r_aa = np.empty(nres, "<U1"); r_aa[rid] = resaa[order]
P = np.array([AA[ONE[a]] for a in r_aa], np.float32)      # (nres, 5)

# ---------------------------------------------------------------- resample per chain
def arc_resample(pts, vals, n=NRES):
    """resample values along the arc length of an ordered Ca trace"""
    if len(pts) == 1:
        return np.repeat(vals[:1], n, axis=0)
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    s = np.concatenate([[0], np.cumsum(seg)])
    if s[-1] == 0:
        return np.repeat(vals[:1], n, axis=0)
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, vals[:, d]) for d in range(vals.shape[1])], 1)


GEOM = np.zeros((len(ids), 2 * NRES, 6), np.float32)      # Ca xyz + centroid xyz
CHEM = np.zeros((len(ids), 2 * NRES, 5), np.float32)
lenA = np.zeros(len(ids), np.int16)
lenB = np.zeros(len(ids), np.int16)
for ci, ch in enumerate("AB"):
    sel = r_chain == ch
    for mi in range(len(ids)):
        k = np.where(sel & (r_model == mi))[0]
        if len(k) == 0:
            continue
        (lenA if ch == "A" else lenB)[mi] = len(k)
        blk = np.concatenate([ca[k], cen[k]], 1)          # (L, 6)
        GEOM[mi, ci * NRES:(ci + 1) * NRES] = arc_resample(ca[k], blk)
        CHEM[mi, ci * NRES:(ci + 1) * NRES] = arc_resample(ca[k], P[k])
# v2 check: the residue order used equals the order the residues appear in the model file
# (D1/D2 keep atoms in file order). Distance is not used as the test: some models have a genuinely
# stretched CDR3 Ca-Ca (5.3-5.5 A) between consecutively numbered residues.
_first = {}
for j in range(len(model)):
    _first.setdefault((int(model[j]), str(chain[j]), int(resnum[j]), str(resins[j])), j)
_r_key = [None] * nres
for j, r in zip(order, rid):
    if _r_key[r] is None:
        _r_key[r] = (int(model[j]), str(chain[j]), int(resnum[j]), str(resins[j]))
_bad = 0
for mi in range(len(ids)):
    for ch in "AB":
        k = np.where((r_chain == ch) & (r_model == mi))[0]
        fo = [_first[_r_key[x]] for x in k]
        _bad += int(any(fo[i] > fo[i + 1] for i in range(len(fo) - 1)))
print("v2 check: CDR3 chains whose residue order differs from file order: %d" % _bad)
assert _bad == 0, "residue order does not follow the model file"
print("CDR3 lengths: alpha %d-%d, beta %d-%d"
      % (lenA.min(), lenA.max(), lenB.min(), lenB.max()))

mu, sd = CHEM.reshape(-1, 5).mean(0), CHEM.reshape(-1, 5).std(0)
CHEMZ = (CHEM - mu) / sd
print("property means: " + "  ".join("%s %.2f" % (p, m_) for p, m_ in zip(PROPS, mu)))

np.savez_compressed(os.path.join(OUT, "D3_property_descriptor.npz"),
                    clone_id=ids, geom=GEOM, chem=CHEM, chemz=CHEMZ.astype(np.float32),
                    props=np.array(PROPS), len_A=lenA, len_B=lenB,
                    chem_mean=mu, chem_sd=sd, n_resample=np.array([NRES]))
print("wrote D3_property_descriptor.npz  geom %s  chem %s"
      % (GEOM.shape, CHEM.shape))
