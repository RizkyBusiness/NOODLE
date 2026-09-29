"""Combining-site geometry for IMGT-numbered paired receptor models.

Represents the antigen-binding surface of a two-chain V module as a
framework-anchored, length-independent point set:

  1. superpose every model onto a common reference using framework C-alpha
     atoms at conserved IMGT positions (CDRs excluded), so the two V domains
     sit in a shared frame and the loops keep their real spatial relationship;
  2. resample each CDR C-alpha trace by arc length to a fixed number of points,
     giving a descriptor of identical size for loops of any length.

Structural distance between two receptors is then the RMSD between their
resampled combining-site point sets in the common frame. Nothing here is
locus- or species-specific: chain labels come from the config, positions from
IMGT unique numbering.
"""
import numpy as np

import imgt

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V"}

# Kyte-Doolittle hydropathy, used as a property channel of the atom cloud
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}   # Bondi radii, heavy atoms
POSITIVE, NEGATIVE, POLAR = set("KRH"), set("DE"), set("STNQYHW")
CHANNELS = ("all", "hydrophobic", "positive", "negative", "polar")


def channel_mask(aas, channel):
    a = np.asarray(list(aas))
    if channel == "all":
        return np.ones(len(a), bool)
    if channel == "hydrophobic":
        return np.array([KD.get(x, 0) > 0 for x in a])
    if channel == "positive":
        return np.isin(a, list(POSITIVE))
    if channel == "negative":
        return np.isin(a, list(NEGATIVE))
    return np.isin(a, list(POLAR))


def usr_moments(pts):
    """USR shape moments: 3 statistics for each of 4 reference locations."""
    if len(pts) < 4:
        return np.zeros(12)
    ctd = pts.mean(0)
    d0 = np.linalg.norm(pts - ctd, axis=1)
    cst, fct = pts[d0.argmin()], pts[d0.argmax()]
    ftf = pts[np.linalg.norm(pts - fct, axis=1).argmax()]
    out = []
    for ref in (ctd, cst, fct, ftf):
        d = np.linalg.norm(pts - ref, axis=1)
        out += [d.mean(), d.std(), np.cbrt(((d - d.mean()) ** 3).mean())]
    return np.array(out)


def usrcat(pts, aas):
    """Rotation- and translation-invariant descriptor of a property-labelled cloud."""
    return np.concatenate([usr_moments(pts[channel_mask(aas, c)]) for c in CHANNELS])


def read_ca(path, chains):
    """C-alpha coordinates of an IMGT-numbered model -> {chain: {(num, ins): xyz}}."""
    out = {c: {} for c in chains}
    with open(path) as fh:
        for L in fh:
            if L.startswith("ATOM") and L[12:16].strip() == "CA":
                ch = L[21]
                if ch in out:
                    out[ch][(int(L[22:26]), L[26].strip())] = np.array(
                        [float(L[30:38]), float(L[38:46]), float(L[46:54])])
    return out


def kabsch(P, Q):
    """Rotation R and centroids mapping P onto Q (both (n,3))."""
    cP, cQ = P.mean(0), Q.mean(0)
    V, _, W = np.linalg.svd((P - cP).T @ (Q - cQ))
    R = V @ np.diag([1, 1, np.sign(np.linalg.det(V @ W))]) @ W
    return R, cP, cQ


def framework_positions(ca_dicts, chains, min_frac=0.995):
    """IMGT framework keys present in at least min_frac of the given models."""
    from collections import Counter
    cnt = Counter()
    for ca in ca_dicts:
        for ch in chains:
            for k in ca[ch]:
                if k[0] in imgt.FR_POS and k[1] == "":
                    cnt[(ch, k)] += 1
    n = len(ca_dicts)
    return sorted([k for k, c in cnt.items() if c >= min_frac * n])


def framework_matrix(ca, keys):
    """(n,3) framework coordinates in fixed key order; None if any key is missing."""
    try:
        return np.array([ca[ch][k] for ch, k in keys])
    except KeyError:
        return None


def resample_loop(pts, n):
    """Resample an ordered C-alpha trace to n points equally spaced along arc length."""
    pts = np.asarray(pts, float)
    if len(pts) == 1:
        return np.repeat(pts, n, axis=0)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    if s[-1] == 0:
        return np.repeat(pts[:1], n, axis=0)
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, pts[:, d]) for d in range(3)], 1)


def site_descriptor(ca, fr_keys, ref_fr, loops, n_resample=10):
    """Superpose on the reference framework, then return the resampled loop points.

    Returns (descriptor (len(loops)*n_resample, 3), framework RMSD, loop lengths).
    """
    F = framework_matrix(ca, fr_keys)
    if F is None:
        return None, None, None
    R, cF, cR = kabsch(F, ref_fr)
    fit = (F - cF) @ R + cR
    fr_rmsd = float(np.sqrt(((fit - ref_fr) ** 2).sum(1).mean()))
    blocks, lens = [], {}
    for ch, nm in loops:
        keys = sorted([k for k in ca[ch] if k[0] in imgt.SITE_RANGES[nm]])
        if not keys:
            return None, None, None
        pts = np.array([ca[ch][k] for k in keys])
        blocks.append(resample_loop((pts - cF) @ R + cR, n_resample))
        lens["%s_%s" % (ch, nm)] = len(keys)
    return np.concatenate(blocks), fr_rmsd, lens


def read_atoms(path):
    """All heavy atoms: (coords, chains, keys, atom names, one-letter residues)."""
    xyz, ch, keys, names, aas = [], [], [], [], []
    with open(path) as fh:
        for L in fh:
            if not L.startswith("ATOM"):
                continue
            el = L[76:78].strip() or L[12:16].strip()[:1]
            rn = L[17:20].strip()
            if el == "H" or rn not in THREE_TO_ONE:
                continue
            xyz.append((float(L[30:38]), float(L[38:46]), float(L[46:54])))
            ch.append(L[21])
            keys.append((int(L[22:26]), L[26].strip()))
            names.append(L[12:16].strip())
            aas.append(THREE_TO_ONE[rn])
    return np.array(xyz), np.array(ch), keys, np.array(names), np.array(aas)


def site_atom_mask(chains_of_atom, keys, chains, include_hv4=True):
    """Atoms belonging to the antigen-binding site on either chain."""
    rngs = imgt.SITE_RANGES if include_hv4 else imgt.CDR_RANGES
    ok = np.zeros(len(keys), bool)
    for i, (c, k) in enumerate(zip(chains_of_atom, keys)):
        if c in chains and any(k[0] in r for r in rngs.values()):
            ok[i] = True
    return ok


def exposed_site_cloud(path, chains, min_sasa=1.0, include_hv4=True):
    """Solvent-exposed heavy atoms of the combining site, with residue labels.

    SASA is computed on the whole two-chain model (Shrake-Rupley), so burial by
    the partner domain is accounted for; only site atoms are then returned.
    """
    import freesasa
    xyz, ch, keys, names, aas = read_atoms(path)
    rad = [VDW.get(n[0], 1.7) for n in names]
    st = freesasa.calcCoord(xyz.flatten().tolist(), rad)
    sasa = np.array([st.atomArea(i) for i in range(len(xyz))])
    m = site_atom_mask(ch, keys, chains, include_hv4) & (sasa >= min_sasa)
    return xyz[m], aas[m], ch[m], [keys[i] for i in np.where(m)[0]], sasa[m]


def site_rmsd_between(ca_a, ca_b, fr_keys, chains):
    """Combining-site CA RMSD between two models over shared site positions.

    Model b is superposed on model a using the framework, then the RMSD is taken
    over site positions (CDR1/2/3 + HV4, both chains) present in both. Used for
    the crystal benchmark, where the two structures are not length-matched.
    """
    Fa, Fb = framework_matrix(ca_a, fr_keys), framework_matrix(ca_b, fr_keys)
    if Fa is None or Fb is None:
        return None, 0
    R, cB, cA = kabsch(Fb, Fa)
    shared = [(c, k) for c in chains for k in ca_a[c]
              if k in ca_b[c] and any(k[0] in r for r in imgt.SITE_RANGES.values())]
    if len(shared) < 10:
        return None, len(shared)
    A = np.array([ca_a[c][k] for c, k in shared])
    B = np.array([(ca_b[c][k] - cB) @ R + cA for c, k in shared])
    return float(np.sqrt(((A - B) ** 2).sum(1).mean())), len(shared)
