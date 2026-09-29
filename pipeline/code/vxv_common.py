"""VXV (the voxel pipeline, amendment A12.1): shared helpers for the descriptor validity panel.

Grid configurations (A12.1): primary F1 per chain, loops only, occupancy; secondary F0 whole / F0 loops / F1 whole
(occupancy) and F1 loops with channels 1-7. Frames and boxes as VX5b (A9, A10: 7 A margins, `out/VX5b_boxes.json`).
Vector descriptors rebuilt from arcs and landmarks exactly as `the reference method/code/vec6_orientation_weight.py` / `vec1`:
  vc  : CDR3 points in each chain's five-landmark frame, 200 vectors;   vg : global ten-landmark frame, 400 vectors
  ori : side-chain unit directions in the chain frames, 20 positions
  primary vector  = sqrt(d_vc^2 + (0.5 lam_o d_ori)^2)            (geometry only, A12.1)
  full <reference arm> = sqrt(d_vc^2 + (lam d_chem)^2 + (0.5 lam_o d_ori)^2)
d_chem as D3 (v2_d3_property_descriptor.py): per CDR3 residue (volume, charge, hbd, hba, kd), resampled to 10 points
along each chain's Ca arc, standardised with the stored D3 mean / sd.
Sequence baseline (A12.9.1): per-loop global alignment, BLOSUM62 (Henikoff & Henikoff 1992), gap open -10, extend -1.
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column

import vxgrid as vg
import vx2a_pilot as PA
import vx5b_frame_diagnostic as FB
from vx6L_loops import loop_of

SIGMA, H = 2.0, 1.0
BEN = PA.BEN
CONFIGS = {                       # name: (frame, atom set, channels)
    "grid_F1L_occ": ("F1", "loops", 1),      # primary grid
    "grid_F0W_occ": ("F0", "whole", 1),
    "grid_F0L_occ": ("F0", "loops", 1),
    "grid_F1W_occ": ("F1", "whole", 1),
    "grid_F1L_7ch": ("F1", "loops", 7),
}
PRIMARY_GRID, PRIMARY_VEC = "grid_F1L_occ", "vec_vcori_geo"
VECTORS = ["vec_vcori_geo", "vec_vcori_full", "vec_vc", "vec_vg"]
LOOPS = {"CDR1": (27, 38), "CDR2": (56, 65), "CDR2.5": (81, 86), "CDR3": (105, 117)}


def boxes_for(frame):
    J = json.load(open(os.path.join(OUT, "VX5b_boxes.json")))["boxes"]
    keys = ["A", "B"] if frame in FB.PERCHAIN else ["AB"]
    return {k: (np.array(J["%s_%s" % (frame, k)]["lo"]), tuple(J["%s_%s" % (frame, k)]["shape"])) for k in keys}


def ikey(n, i):
    o = 0 if i == "" else ord(i.upper()) - 64
    return (n, -o) if n in (33, 61, 112) else (n, o)


def residues(meta):
    """ordered residue list [(chain, num, ins, resname)] in IMGT sequence order, per chain A then B"""
    seen = {}
    for a in meta:
        seen.setdefault((a[0], a[1], a[2]), a[3])
    out = []
    for ch in "AB":
        ks = sorted([k for k in seen if k[0] == ch], key=lambda k: ikey(k[1], k[2]))
        out += [(k[0], k[1], k[2], seen[k]) for k in ks]
    return out


def missing_positions(pc, pm):
    """IMGT positions whose residue is in the model but has no atom in the crystal"""
    mc, mm = parse_heavy(pc)[1], parse_heavy(pm)[1]
    posc = {(a[0], a[1], a[2]) for a in mc}
    return sorted({(a[0], a[1], a[2]) for a in mm} - posc)


def shared_keys(pc, pm):
    return {tuple(a) for a in parse_heavy(pc)[1]} & {tuple(a) for a in parse_heavy(pm)[1]}


def bench_place(path, frame, lm3):
    if frame == "F0":
        return {"AB": PA.bench_frame(path, lm3["parse"], lm3["LMPOS"])}
    FB.PA.lm3_LMPOS = lm3["LMPOS"]
    lm_, an_, ok_ = FB.bench_points(path, lm3["parse"])
    return FB.placements(frame, lm_, an_, ok_)


def rep_place(i, frame):
    if frame == "F0":
        return {"AB": (PA.FR["R_box"][i], PA.FR["t_box"][i])}
    return FB.placements(frame, FB.LMK[i], FB.ANC[i], np.ones(123, bool))


def select(meta, atomset, keys=None, drop_pos=None):
    m = np.ones(len(meta), bool)
    if atomset == "loops":
        m &= np.array([loop_of(a[1]) is not None for a in meta])
    elif atomset == "cdr3":                      # VXV3 (A12.10.2): IMGT CDR3 105-117 only
        m &= np.array([105 <= a[1] <= 117 for a in meta])
    if keys is not None:
        m &= np.array([tuple(a) in keys for a in meta])
    if drop_pos:
        m &= np.array([(a[0], a[1], a[2]) not in drop_pos for a in meta])
    return m


def build_grid(args):
    """(path, place, boxes, atomset, nch, keys, drop_pos, invert[, sigma]) -> flattened float64 grid over the boxes, n atoms,
    atom keys used. invert=True builds ONLY the atoms that drop_pos would remove (the removal grid). sigma: optional 9th
    element (A14; default SIGMA = 2.0, so every earlier call is unchanged)."""
    path, place, boxes, atomset, nch, keys, drop_pos, invert = args[:8]
    sigma = args[8] if len(args) > 8 else SIGMA
    xyz, meta = parse_heavy(path)
    W, bad = vg.type_atoms(meta)
    m = select(meta, atomset, keys)
    if drop_pos is not None:
        inpos = np.array([(a[0], a[1], a[2]) in drop_pos for a in meta])
        m &= inpos if invert else ~inpos
    ch = np.array([a[0] for a in meta])
    vecs = []
    for key, (Rb, tb) in place.items():
        s = m & np.isin(ch, list(key))
        lo, shp = boxes[key]
        if not s.any():                  # e.g. a removal grid with no atoms in this chain's box
            vecs.append(np.zeros(nch * int(np.prod(shp)))); continue
        vecs.append(vg.build(xyz[s] @ Rb + tb, W[s][:, :nch], sigma, H, lo, shp).reshape(-1))
    used = frozenset(tuple(a) for a, k in zip(meta, m) if k)
    return np.concatenate(vecs), int(m.sum()), len(bad), used


def structure_atoms(path, place, atomset, nch, keys=None):
    """per-box atom sets in the box frame, for the analytic check"""
    xyz, meta = parse_heavy(path)
    W, _ = vg.type_atoms(meta)
    m = select(meta, atomset, keys); ch = np.array([a[0] for a in meta])
    return [(xyz[m & np.isin(ch, list(k))] @ Rb + tb, W[m & np.isin(ch, list(k))][:, :nch]) for k, (Rb, tb) in place.items()]


# ------------------------------------------------------------------ vector features (vec6 / vec1 recipe)
L1 = np.load(LM1, allow_pickle=True)
T10 = L1["lm"][list(L1["clone_id"]).index(REF_ID)].astype(float)
PA_ = list(range(10)) + list(range(20, 30)); PB_ = list(range(10, 20)) + list(range(30, 40))


def kab(X, Y):
    cx, cy = X.mean(0), Y.mean(0); U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy)); d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1, 1, d]) @ Vt


def orient(arc, R, pos):
    v = arc[[20 + p for p in pos]] - arc[pos]; ln = np.linalg.norm(v, axis=1, keepdims=True)
    return np.where(ln >= 0.5, v / np.maximum(ln, 1e-9), 0.0) @ R


def vec_feats(arc, lm):
    """vc (600), uc (60), vg (1200) exactly as vec6 feats() and vec1 vg"""
    Ra, Rb = kab(lm[:5], T10[:5]), kab(lm[5:], T10[5:])
    vc = np.concatenate([((arc[PA_][:, None] - lm[:5][None]) @ Ra).reshape(-1),
                         ((arc[PB_][:, None] - lm[5:][None]) @ Rb).reshape(-1)])
    uc = np.concatenate([orient(arc, Ra, list(range(10))), orient(arc, Rb, list(range(10, 20)))]).reshape(-1)
    Rg = kab(lm, T10)
    vgf = ((arc[:, None] - lm[None]) @ Rg).reshape(-1)
    return vc, uc, vgf


def vec_params():
    th = pd.read_csv(VXP(_REF2.get("n10_thresholds_alt", "")))
    return float(th.lam.iloc[0]), 0.5 * float(th.lam_ori.iloc[0])      # lam, w=0.5 x lam_o (V6 chosen w)


# ------------------------------------------------------------------ chemistry (D3 recipe)
D3 = np.load(VXP("descriptors/out/D3_property_descriptor.npz"), allow_pickle=True)
D2AA = np.load(VXP("descriptors/out/D2_cdr3_atoms.npz"), allow_pickle=True)["aa_table"]
AAP = {r[0]: np.array([float(x) for x in r[1:]]) for r in D2AA}


def arc_resample(pts, vals, n=10):
    if len(pts) == 1:
        return np.repeat(vals[:1], n, axis=0)
    s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))])
    if s[-1] == 0:
        return np.repeat(vals[:1], n, axis=0)
    t = np.linspace(0, s[-1], n)
    return np.stack([np.interp(t, s, vals[:, d]) for d in range(vals.shape[1])], 1)


def chemz_of(path):
    xyz, meta = parse_heavy(path)
    ca = {(a[0], a[1], a[2]): (x, a[3]) for x, a in zip(xyz, meta) if a[4] == "CA" and 105 <= a[1] <= 117}
    blocks = []
    for ch in "AB":
        ks = sorted([k for k in ca if k[0] == ch], key=lambda k: ikey(k[1], k[2]))
        pts = np.array([ca[k][0] for k in ks]); P = np.array([AAP[ca[k][1]] for k in ks])
        blocks.append(arc_resample(pts, P))
    C = np.concatenate(blocks)
    return (C - D3["chem_mean"]) / D3["chem_sd"]


# ------------------------------------------------------------------ sequence baseline
def loop_seqs(path):
    meta = parse_heavy(path)[1]
    R = residues(meta)
    ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
           "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
           "TYR": "Y", "VAL": "V"}
    out = {}
    for ch in "AB":
        for nm, (a, b) in LOOPS.items():
            out[(ch, nm)] = "".join(ONE[r[3]] for r in R if r[0] == ch and a <= r[1] <= b)
    return out


def aligner():
    from Bio import Align
    from Bio.Align import substitution_matrices
    al = Align.PairwiseAligner()
    al.mode = "global"; al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score = -10, -1
    return al


def vregion_identity(pa, pb):
    """positional identity over IMGT 1-104 on both chains (same-V-pair proxy, report-only strata)"""
    def seq(p):
        return {(r[0], r[1], r[2]): r[3] for r in residues(parse_heavy(p)[1]) if r[1] <= 104}
    A, B = seq(pa), seq(pb)
    out = []
    for ch in "AB":
        ka = {k for k in A if k[0] == ch}; kb = {k for k in B if k[0] == ch}
        u = ka | kb
        out.append(sum(1 for k in ka & kb if A[k] == B[k]) / len(u))
    return out
