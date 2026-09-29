"""Pipeline step 1 of the reference-point / vector descriptor: a TCR-intrinsic coordinate system.

Built ONCE on the D2 reference framework (the medoid model's 123 anchor Ca, D2_cdr3_atoms.npz), so
the origin and axes are one fixed point and one fixed set of directions in the D2 frame. Every
model already superposed on D2 is then expressed in these axes by a single rigid transform, which
leaves all pairwise distances unchanged.

  origin  O  centroid of the four CDR3-stem anchors  A104 A118 B104 B118  (Cys104, J-Phe/Trp118)
  z          Va/Vb pseudo-twofold axis: rotation axis of the Kabsch fit of the alpha anchors onto
             the beta anchors at the same IMGT numbers; sign chosen to point from O toward the
             repertoire-mean CDR3 Ca (i.e. up, toward the binding face)
  x          Va -> Vb: centroid of Ca 23+104 (the intradomain disulfide) on beta minus on alpha,
             made perpendicular to z (Rudolph, Stanfield & Wilson 2006 use the same disulfide line)
  y          z cross x   (right-handed)

Checks written:
  - rotation angle of the alpha->beta fit (a true twofold = 180 deg) and screw translation
  - distance of O from the dyad axis
  - bootstrap over paired positions: angular spread of z and x
  - independent validation on the crystal benchmark: each solved TCR-pMHC complex is fitted onto
    the same 123 anchors and the peptide / MHC position read out in the new axes. If z is really
    "toward the pMHC", the peptide must sit on +z, close to the axis.

Usage: python3 h15_intrinsic_frame.py <atoms dir> <tables dir> <out dir>
Numpy + pandas only.
"""
import os, sys, json
import numpy as np, pandas as pd

STD, STB, OUT = sys.argv[1:4]
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0)

def kabsch(P, Q):
    """R, cP, cQ such that (P - cP) @ R + cQ ~ Q  (proper rotation)."""
    cP, cQ = P.mean(0), Q.mean(0)
    H = (P - cP).T @ (Q - cQ)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    return U @ D @ Vt, cP, cQ

def axis_angle(R):
    """Rotation axis (unit) and angle (deg) of a row-vector rotation x @ R."""
    M = R.T                                   # column-vector convention
    w, v = np.linalg.eig(M)
    u = np.real(v[:, np.argmin(np.abs(w - 1))]); u /= np.linalg.norm(u)
    ang = np.degrees(np.arccos(np.clip((np.trace(M) - 1) / 2, -1, 1)))
    return u, ang

def unit(v): return v / np.linalg.norm(v)

# ------------------------------------------------------------------ the D2 reference framework
z = np.load(os.path.join(STD, "D2_cdr3_atoms.npz"), allow_pickle=True)
REF = z["ref_framework"].astype(np.float64)
KEYS = [str(k) for k in z["anchor_keys"]]
K = {k: i for i, k in enumerate(KEYS)}
for k in ("A104", "A118", "B104", "B118", "A23", "B23"):
    assert k in K, "anchor %s missing from D2" % k

O = REF[[K["A104"], K["A118"], K["B104"], K["B118"]]].mean(0)

# alpha/beta pairs at the same IMGT number (both present as anchors)
paired = sorted({int(k[1:]) for k in KEYS if k[0] == "A"} & {int(k[1:]) for k in KEYS if k[0] == "B"})
PA = REF[[K["A%d" % n] for n in paired]]
PB = REF[[K["B%d" % n] for n in paired]]

def dyad(idx):
    R, cA, cB = kabsch(PA[idx], PB[idx])
    u, ang = axis_angle(R)
    return u, ang, R, cA, cB

u, ang, R_ab, cA, cB = dyad(np.arange(len(paired)))
fit_ab = np.sqrt((((PA - cA) @ R_ab + cB - PB) ** 2).sum(1).mean())
# screw translation along the axis, and a point on the axis (midpoint of the two domain centroids
# lies on it for an exact twofold; report its distance anyway)
screw = float(np.dot(cB - cA, u))

# sign: up = toward the CDR3 loops
ca = z["atom"] == "CA"
cdr3_mean = z["xyz"][ca].astype(np.float64).mean(0)
if np.dot(cdr3_mean - O, u) < 0: u = -u
Z = u
d_A = REF[[K["A23"], K["A104"]]].mean(0)
d_B = REF[[K["B23"], K["B104"]]].mean(0)
xv = d_B - d_A
X = unit(xv - np.dot(xv, Z) * Z)
Y = np.cross(Z, X)
AX = np.vstack([X, Y, Z])                     # rows: x, y, z ; local = (p - O) @ AX.T

def dist_point_line(p, a, dirn):
    w = p - a
    return float(np.linalg.norm(w - np.dot(w, dirn) * dirn))
O_off_axis = dist_point_line(O, (cA + cB) / 2, Z)

# alternative "up": framework centroid -> stem centroid (independent of the dyad)
z_alt = unit(O - REF.mean(0))
ang_zalt = float(np.degrees(np.arccos(np.clip(np.dot(z_alt, Z), -1, 1))))

# bootstrap over paired positions
bz, bx = [], []
for _ in range(1000):
    idx = rng.choice(len(paired), len(paired), replace=True)
    ub, *_ = dyad(idx)
    if np.dot(ub, Z) < 0: ub = -ub
    xb = unit(xv - np.dot(xv, ub) * ub)
    bz.append(np.degrees(np.arccos(np.clip(np.dot(ub, Z), -1, 1))))
    bx.append(np.degrees(np.arccos(np.clip(np.dot(xb, X), -1, 1))))
bz, bx = np.array(bz), np.array(bx)

# where the repertoire's CDR3 Ca sit in the new axes (sanity: should be +z)
loc_cdr3 = (z["xyz"][ca].astype(np.float64) - O) @ AX.T
cdr3_z_frac_pos = float((loc_cdr3[:, 2] > 0).mean())

checks = dict(
    n_paired_positions=len(paired),
    alpha_to_beta_rotation_deg=round(float(ang), 2),
    alpha_to_beta_fit_rmsd_A=round(float(fit_ab), 3),
    screw_translation_A=round(screw, 3),
    origin_distance_from_dyad_axis_A=round(O_off_axis, 3),
    angle_z_vs_framework_to_stem_deg=round(ang_zalt, 2),
    bootstrap_z_angle_median_deg=round(float(np.median(bz)), 2),
    bootstrap_z_angle_p95_deg=round(float(np.percentile(bz, 95)), 2),
    bootstrap_x_angle_median_deg=round(float(np.median(bx)), 2),
    bootstrap_x_angle_p95_deg=round(float(np.percentile(bx, 95)), 2),
    angle_x_raw_vs_z_deg=round(float(np.degrees(np.arccos(np.dot(unit(xv), Z)))), 2),
    cdr3_ca_fraction_above_origin=round(cdr3_z_frac_pos, 4),
    cdr3_ca_mean_local_xyz=[round(float(v), 2) for v in loc_cdr3.mean(0)],
    orthonormal_err=float(np.abs(AX @ AX.T - np.eye(3)).max()),
    det=float(np.linalg.det(AX)),
)
print(json.dumps(checks, indent=1))

# ------------------------------------------------------------------ crystal validation
def read_pdb(path):
    rows = []
    for l in open(path):
        if l.startswith("ATOM") and l[16] in " A":
            rows.append((l[21], int(l[22:26]), l[26], l[12:16].strip(), l[17:20],
                         float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return pd.DataFrame(rows, columns=["ch", "num", "ins", "atom", "res", "x", "y", "z"])

bench = pd.read_csv(os.path.join(STB, "B3h_benchmark_structures.csv"))
bench = bench[bench.accepted.astype(bool)]
bdir = os.path.join(STB, "benchmark")
val = []
for _, r in bench.iterrows():
    e = r.entry
    try:
        fx = read_pdb(os.path.join(bdir, "fixed_%s_crystal.pdb" % e))
        raw = read_pdb(os.path.join(bdir, "%s.pdb" % e))
    except FileNotFoundError:
        continue
    fca = fx[(fx.atom == "CA") & (fx.ins == " ")]
    lab = {("A%d" % n): v for n, v in zip(fca[fca.ch == "A"].num, fca[fca.ch == "A"][["x", "y", "z"]].values)}
    lab.update({("B%d" % n): v for n, v in zip(fca[fca.ch == "B"].num, fca[fca.ch == "B"][["x", "y", "z"]].values)})
    common = [k for k in KEYS if k in lab]
    if len(common) < 100:
        val.append(dict(entry=e, note="anchors %d" % len(common))); continue
    P = np.array([lab[k] for k in common]); Q = REF[[K[k] for k in common]]
    R, cP, cQ = kabsch(P, Q)
    fr = float(np.sqrt((((P - cP) @ R + cQ - Q) ** 2).sum(1).mean()))
    tcr = raw[raw.ch.isin([r.auth_A, r.auth_B])][["x", "y", "z"]].values
    rec = dict(entry=e, anchors=len(common), fr_fit=round(fr, 3))
    others = raw[~raw.ch.isin([r.auth_A, r.auth_B])]
    pep, mhc = None, []
    for ch, g in others.groupby("ch"):
        xyz = g[["x", "y", "z"]].values
        # contact with THIS TCR (asymmetric units can hold several complexes)
        dmin = min(np.sqrt(((tcr[i:i + 400, None] - xyz[None]) ** 2).sum(2)).min()
                   for i in range(0, len(tcr), 400))
        if dmin > 5.0: continue
        nres = g[g.atom == "CA"].shape[0]
        if 5 <= nres <= 30: pep = g[g.atom == "CA"][["x", "y", "z"]].values
        elif nres > 60: mhc.append(g[g.atom == "CA"][["x", "y", "z"]].values)
    if pep is None:
        rec["note"] = "no contacting peptide (TCR alone or not pMHC)"; val.append(rec); continue
    tolocal = lambda p: (((p - cP) @ R + cQ) - O) @ AX.T
    pl = tolocal(pep); c = pl.mean(0)
    rec.update(pep_len=len(pep), pep_x=c[0], pep_y=c[1], pep_z=c[2],
               pep_angle_from_z_deg=float(np.degrees(np.arccos(c[2] / np.linalg.norm(c)))),
               pep_off_axis_A=float(np.hypot(c[0], c[1])))
    v = pl[-1] - pl[0]; v[2] = 0            # peptide N->C projected on the xy plane
    rec["crossing_angle_deg"] = float(np.degrees(np.arctan2(v[1], v[0])))
    if mhc:
        m = tolocal(np.vstack(mhc)).mean(0)
        rec.update(mhc_z=m[2], mhc_angle_from_z_deg=float(np.degrees(np.arccos(m[2] / np.linalg.norm(m)))))
    val.append(rec)
V = pd.DataFrame(val)
V.to_csv(os.path.join(OUT, "H15_crystal_validation.csv"), index=False)
C = V.dropna(subset=["pep_z"]) if "pep_z" in V else V.iloc[0:0]
summ = dict(n_accepted=int(len(bench)), n_with_contacting_pMHC=int(len(C)))
if len(C):
    summ.update(
        peptide_above_origin=int((C.pep_z > 0).sum()),
        pep_height_z_median_A=round(float(C.pep_z.median()), 2),
        pep_height_z_range_A=[round(float(C.pep_z.min()), 2), round(float(C.pep_z.max()), 2)],
        pep_angle_from_z_median_deg=round(float(C.pep_angle_from_z_deg.median()), 1),
        pep_angle_from_z_p90_deg=round(float(C.pep_angle_from_z_deg.quantile(.9)), 1),
        pep_off_axis_median_A=round(float(C.pep_off_axis_A.median()), 2),
        crossing_angle_median_deg=round(float(C.crossing_angle_deg.median()), 1),
        crossing_angle_IQR_deg=[round(float(C.crossing_angle_deg.quantile(.25)), 1),
                                round(float(C.crossing_angle_deg.quantile(.75)), 1)],
        crystal_fr_fit_median_A=round(float(C.fr_fit.median()), 3))
print(json.dumps(summ, indent=1))

np.savez(os.path.join(OUT, "H15_intrinsic_frame.npz"), origin=O, axes=AX,
         paired_positions=np.array(paired), anchor_keys=np.array(KEYS),
         note=np.array("local = (p_in_D2_frame - origin) @ axes.T ; rows of axes = x, y, z"))
json.dump(dict(checks=checks, crystal_validation=summ,
               origin=O.tolist(), axes=AX.tolist()),
          open(os.path.join(OUT, "H15_frame_checks.json"), "w"), indent=1)
print("wrote H15_intrinsic_frame.npz, H15_frame_checks.json, H15_crystal_validation.csv")
