"""the voxel pipeline grid library: per-atom channel typing (DESIGN VX2 + amendment A2) and Gaussian smearing onto the frozen box.

Channels (index: name):
  0 occupancy     every heavy atom, mass 1
  1 hydrophobic   carbon atoms not bonded to N or O
  2 aromatic      ring atoms of Phe, Tyr, Trp, His
  3 donor         backbone N except Pro (A3.2); Arg NE/NH1/NH2, Asn ND2, Gln NE2,
                  Lys NZ, Trp NE1, His ND1/NE2 (tautomer-agnostic); Ser OG, Thr OG1, Tyr OH
  4 acceptor      all O; His ND1/NE2
  5 positive      Lys NZ 1; Arg CZ/NE/NH1/NH2 0.25; His ND1/NE2 0.05        (per formal charge, A2.2)
  6 negative      Asp OD1/OD2 0.5; Glu OE1/OE2 0.5; OXT 1
  7 side chain    every non-backbone heavy atom (stored; not in the primary distance, A1.1)
Primary distance: channels 0-6.

Smearing (following the voxelised pharmacophore channels of Jimenez et al. 2017, DeepSite): each atom is an
isotropic Gaussian of width sigma, point-sampled at voxel centres, truncated per axis at |c - x| <= 3 sigma, and each
1-D factor renormalised to sum 1, so every atom deposits exactly its channel mass unless the box clips it. Stored
values are mass per voxel; the continuous L2 inner product is recovered as sum(v_a * v_b) / h^3.
"""
import math
import numpy as np

NCH, NPRIMARY = 8, 7
CHANNELS = ["occupancy", "hydrophobic", "aromatic", "donor", "acceptor", "positive", "negative", "sidechain"]
BB = {"N", "CA", "C", "O", "OXT"}
SIDE = {
    "ALA": "CB", "ARG": "CB CG CD NE CZ NH1 NH2", "ASN": "CB CG OD1 ND2", "ASP": "CB CG OD1 OD2", "CYS": "CB SG",
    "GLN": "CB CG CD OE1 NE2", "GLU": "CB CG CD OE1 OE2", "GLY": "", "HIS": "CB CG ND1 CD2 CE1 NE2",
    "ILE": "CB CG1 CG2 CD1", "LEU": "CB CG CD1 CD2", "LYS": "CB CG CD CE NZ", "MET": "CB CG SD CE",
    "PHE": "CB CG CD1 CD2 CE1 CE2 CZ", "PRO": "CB CG CD", "SER": "CB OG", "THR": "CB OG1 CG2",
    "TRP": "CB CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2", "TYR": "CB CG CD1 CD2 CE1 CE2 CZ OH", "VAL": "CB CG1 CG2"}
SIDE = {k: set(v.split()) for k, v in SIDE.items()}
# carbons bonded to N or O (everything else that is carbon is hydrophobic)
POLAR_C = {("*", "CA"), ("*", "C"), ("ARG", "CD"), ("ARG", "CZ"), ("ASN", "CG"), ("ASP", "CG"), ("GLN", "CD"),
           ("GLU", "CD"), ("HIS", "CG"), ("HIS", "CD2"), ("HIS", "CE1"), ("LYS", "CE"), ("PRO", "CD"),
           ("SER", "CB"), ("THR", "CB"), ("TRP", "CD1"), ("TRP", "CE2"), ("TYR", "CZ")}
AROM = {"PHE": "CG CD1 CD2 CE1 CE2 CZ", "TYR": "CG CD1 CD2 CE1 CE2 CZ",
        "TRP": "CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2", "HIS": "CG ND1 CD2 CE1 NE2"}
AROM = {k: set(v.split()) for k, v in AROM.items()}
DON_SC = {("ARG", "NE"), ("ARG", "NH1"), ("ARG", "NH2"), ("ASN", "ND2"), ("GLN", "NE2"), ("LYS", "NZ"),
          ("TRP", "NE1"), ("HIS", "ND1"), ("HIS", "NE2"), ("SER", "OG"), ("THR", "OG1"), ("TYR", "OH")}
POS = {("LYS", "NZ"): 1.0, ("ARG", "CZ"): 0.25, ("ARG", "NE"): 0.25, ("ARG", "NH1"): 0.25, ("ARG", "NH2"): 0.25,
       ("HIS", "ND1"): 0.05, ("HIS", "NE2"): 0.05}
NEG = {("ASP", "OD1"): 0.5, ("ASP", "OD2"): 0.5, ("GLU", "OE1"): 0.5, ("GLU", "OE2"): 0.5}


def type_atoms(meta):
    """meta: list of (chain, resnum, ins, resname, atom). Returns W (n, 8) channel masses and the list of atoms that
    match no rule (unknown residue or atom name)."""
    W = np.zeros((len(meta), NCH))
    bad = []
    for i, (ch, num, ins, rn, at) in enumerate(meta):
        if rn not in SIDE or not (at in BB or at in SIDE[rn]):
            bad.append((ch, num, ins, rn, at)); continue
        el = at[0]
        w = W[i]
        w[0] = 1.0
        if el == "C" and ("*", at) not in POLAR_C and (rn, at) not in POLAR_C:
            w[1] = 1.0
        if at in AROM.get(rn, ()):
            w[2] = 1.0
        if at == "N" and rn != "PRO":
            w[3] = 1.0
        if (rn, at) in DON_SC:
            w[3] = 1.0
        if el == "O" or (rn, at) in (("HIS", "ND1"), ("HIS", "NE2")):
            w[4] = 1.0
        w[5] = POS.get((rn, at), 0.0)
        w[6] = 1.0 if at == "OXT" else NEG.get((rn, at), 0.0)
        if at not in BB:
            w[7] = 1.0
    return W, bad


def _kernels(X, sigma, h, lo, shape, trunc):
    """per-axis sampled, truncated, renormalised 1-D Gaussian factors and their voxel indices."""
    k = int(math.ceil(2 * trunc * sigma / h)) + 2
    idx, ker = [], []
    for d in range(3):
        u = (X[:, d] - lo[d]) / h - 0.5                       # continuous voxel coordinate of the atom
        i0 = np.floor(u - trunc * sigma / h).astype(int)
        ii = i0[:, None] + np.arange(k)[None]
        off = (ii - u[:, None]) * h                            # voxel centre minus atom, A
        w = np.exp(-0.5 * (off / sigma) ** 2)
        w[np.abs(off) > trunc * sigma] = 0.0
        w /= w.sum(1, keepdims=True)
        out = (ii < 0) | (ii >= shape[d])
        w[out] = 0.0
        idx.append(np.clip(ii, 0, shape[d] - 1)); ker.append(w)
    return idx, ker


def build(X, W, sigma, h, lo, shape, trunc=3.0):
    """fast grid: returns float64 (8, nx, ny, nz) of mass per voxel."""
    shape = tuple(int(s) for s in shape)
    (ix, iy, iz), (kx, ky, kz) = _kernels(np.asarray(X, float), sigma, h, np.asarray(lo, float), shape, trunc)
    K = (kx[:, :, None, None] * ky[:, None, :, None] * kz[:, None, None, :]).reshape(len(X), -1)
    I = ((ix[:, :, None, None] * shape[1] + iy[:, None, :, None]) * shape[2] + iz[:, None, None, :]).reshape(len(X), -1)
    nv = int(np.prod(shape))
    G = np.empty((W.shape[1], nv))
    If = I.ravel()
    for c in range(W.shape[1]):
        G[c] = np.bincount(If, (K * W[:, c:c + 1]).ravel(), minlength=nv)
    return G.reshape((W.shape[1],) + shape)


def build_slow(X, W, sigma, h, lo, shape, trunc=3.0):
    """obviously-correct reference: explicit loops over atoms and the voxels within trunc*sigma on each axis."""
    G = np.zeros((W.shape[1],) + tuple(shape))
    for a in range(len(X)):
        fac = []
        for d in range(3):
            vals = {}
            for j in range(int(shape[d])):
                c = lo[d] + (j + 0.5) * h
                if abs(c - X[a][d]) <= trunc * sigma:
                    vals[j] = math.exp(-0.5 * ((c - X[a][d]) / sigma) ** 2)
            # normalise over the full (unclipped) window: include centres beyond the box edge
            tot, j = 0.0, math.floor((X[a][d] - trunc * sigma - lo[d]) / h - 0.5) - 1
            while lo[d] + (j + 0.5) * h <= X[a][d] + trunc * sigma + h:
                c = lo[d] + (j + 0.5) * h
                if abs(c - X[a][d]) <= trunc * sigma:
                    tot += math.exp(-0.5 * ((c - X[a][d]) / sigma) ** 2)
                j += 1
            fac.append({j: v / tot for j, v in vals.items()})
        for i, wi in fac[0].items():
            for j, wj in fac[1].items():
                for k, wk in fac[2].items():
                    for c in range(W.shape[1]):
                        if W[a, c]:
                            G[c, i, j, k] += W[a, c] * wi * wj * wk
    return G


def analytic_d2(XA, WA, XB, WB, sigma, channels=range(NPRIMARY)):
    """exact ||rho_a - rho_b||^2 of the untruncated Gaussian-smeared atom sets, summed over channels:
    sum_ij G(a_i - a_j) + sum_ij G(b_i - b_j) - 2 sum_ij G(a_i - b_j), G a normalised Gaussian of width sigma*sqrt(2)."""
    from scipy.spatial.distance import cdist
    c = (4 * math.pi * sigma * sigma) ** -1.5
    def S(P, wp, Q, wq):
        E = c * np.exp(-cdist(P, Q, "sqeuclidean") / (4 * sigma * sigma))
        return sum(float(wp[:, ch] @ E @ wq[:, ch]) for ch in channels)
    return S(XA, WA, XA, WA) + S(XB, WB, XB, WB) - 2 * S(XA, WA, XB, WB)


def grid_d2(GA, GB, h, channels=range(NPRIMARY)):
    ch = list(channels)
    d = GA[ch].astype(np.float64) - GB[ch].astype(np.float64)
    return float((d * d).sum()) / h ** 3
