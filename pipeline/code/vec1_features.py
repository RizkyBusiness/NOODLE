"""reference step 1 -- vector descriptors of the CDR3 points relative to framework landmarks, in local frames.

Input: landmarks/out/LM1_internal_coords.npz (each model in its own coordinates; arc = 20 Ca then 20 side-chain
centroid points, CDR3a points 0-9 / 20-29, CDR3b 10-19 / 30-39; lm = Ca of IMGT A23 A41 A89 A104 A118 B23 ... B118).

Frames (one Kabsch rotation per molecule onto the reference receptor's landmarks (config reference_receptor; the medoid); row-vector convention
R = U diag(1,1,d) V^T with d = sign det(U V^T), applied as (X - c) @ R):
  vg   global : fitted on the 10 landmarks; vectors from all 40 points to all 10 landmarks        (400 vectors)
  vc   chain  : CDR3a points in the frame of the 5 alpha landmarks, vectors to those 5; same for beta (200)
  vs   stem   : fitted on the 4 CDR3-stem landmarks A104 A118 B104 B118; 40 points x 4 landmarks        (160)
Vector v_pk = (p - L_k) @ R. Geometry distance d_vec = sqrt(mean_pk |v_pk,a - v_pk,b|^2) (Angstrom).
Orientation: u_p = unit(centroid_p - Ca_p) @ R for the 20 positions (zero where |centroid - Ca| < 0.5 A, i.e.
glycine or near-glycine interpolations); d_ori = sqrt(mean_p |u_p,a - u_p,b|^2) (dimensionless, 0..2).
Kabsch needs the landmark set to be non-collinear (rank >= 2); the 2nd singular value is reported.
Writes reference/out/V1_features.npz
"""
import numpy as np, sys
import paths  # package: config (reference receptor, data-set settings)
LM1, OUT = sys.argv[1:3]
L = np.load(LM1, allow_pickle=True)
ids = list(L["clone_id"]); arc = L["arc"].astype(float); lm = L["lm"].astype(float); n = len(ids)
T = lm[ids.index(paths.reference_receptor())]
def kabsch_rows(X, Y):
    cx, cy = X.mean(0), Y.mean(0); U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    d = np.sign(np.linalg.det(U @ Vt)); return U @ np.diag([1, 1, d]) @ Vt, cx, cy
A_IDX, B_IDX = list(range(5)), list(range(5, 10))
PA = list(range(10)) + list(range(20, 30)); PB = list(range(10, 20)) + list(range(30, 40))
def frame(i, ix):
    R, cx, cy = kabsch_rows(lm[i][ix], T[ix]); fit = np.sqrt((((lm[i][ix] - cx) @ R - (T[ix] - cy)) ** 2).sum(1).mean())
    return R, fit
def orient(i, R, pos):
    v = arc[i][[20 + p for p in pos]] - arc[i][pos]; ln = np.linalg.norm(v, axis=1, keepdims=True)
    u = np.where(ln >= 0.5, v / np.maximum(ln, 1e-9), 0.0); return u @ R
out = {k: [] for k in ("vg", "vc", "vs", "ug", "uc", "us", "fit_g", "fit_a", "fit_b", "fit_s", "sv2_g", "sv2_a", "sv2_b", "sv2_s")}
sv2 = lambda X: np.linalg.svd(X - X.mean(0), compute_uv=False)[1]
for i in range(n):
    Rg, fg = frame(i, list(range(10)))
    out["vg"].append(((arc[i][:, None] - lm[i][None]) @ Rg).reshape(-1)); out["ug"].append(orient(i, Rg, list(range(20))).reshape(-1))
    Ra, fa = frame(i, A_IDX); Rb, fb = frame(i, B_IDX)
    va = (arc[i][PA][:, None] - lm[i][A_IDX][None]) @ Ra; vb = (arc[i][PB][:, None] - lm[i][B_IDX][None]) @ Rb
    out["vc"].append(np.concatenate([va.reshape(-1), vb.reshape(-1)]))
    out["uc"].append(np.concatenate([orient(i, Ra, list(range(10))), orient(i, Rb, list(range(10, 20)))]).reshape(-1))
    SI = [3, 4, 8, 9]; Rs, fs = frame(i, SI)
    out["vs"].append(((arc[i][:, None] - lm[i][SI][None]) @ Rs).reshape(-1)); out["us"].append(orient(i, Rs, list(range(20))).reshape(-1))
    for k, v in (("fit_g", fg), ("fit_a", fa), ("fit_b", fb), ("fit_s", fs)): out[k].append(v)
    for k, ix in (("sv2_g", list(range(10))), ("sv2_a", A_IDX), ("sv2_b", B_IDX), ("sv2_s", SI)): out[k].append(sv2(lm[i][ix]))
O = {k: np.array(v, np.float32 if k[0] in "vu" else float) for k, v in out.items()}
np.savez_compressed(OUT, clone_id=np.array(ids), **O, n_vec=np.array([400, 200, 160]))
print({k: v.shape for k, v in O.items()})
for k in ("fit_g", "fit_a", "fit_b", "fit_s", "sv2_g", "sv2_a", "sv2_b", "sv2_s"):
    print("%-6s median %.3f  p99 %.3f  min %.3f" % (k, np.median(O[k]), np.percentile(O[k], 99), O[k].min()))
# ---- identity check (vg): d_vec^2 = RMSD_40(points in the 10-landmark frame)^2 + mean_k |dL_k|^2 (exact, landmark-centred)
rng = np.random.default_rng(3); a, b = rng.integers(0, n, 2000), rng.integers(0, n, 2000)
def placed(i):
    R, cx, cy = kabsch_rows(lm[i], T); return (arc[i] - cx) @ R, (lm[i] - cx) @ R
err = 0.0
for x, y in zip(a, b):
    pa_, la_ = placed(x); pb_, lb_ = placed(y)
    lhs = ((O["vg"][x].astype(float) - O["vg"][y].astype(float)) ** 2).sum() / 400
    rhs = ((pa_ - pb_) ** 2).sum(1).mean() + ((la_ - lb_) ** 2).sum(1).mean()
    err = max(err, abs(lhs - rhs))
print("identity check over 2000 pairs: max |d_vec^2 - (d_geom10^2 + landmark term)| = %.2e A^2" % err)
# ---- rotation invariance: rigidly move a model -> identical features
from scipy.spatial.transform import Rotation
w = 0.0
for i in rng.choice(n, 200, replace=False):
    Q = Rotation.random(random_state=int(rng.integers(1e9))).as_matrix(); t = rng.normal(0, 40, 3)
    a2, l2 = arc[i] @ Q.T + t, lm[i] @ Q.T + t
    R, cx, cy = kabsch_rows(l2, T); v2 = ((a2[:, None] - l2[None]) @ R).reshape(-1)
    w = max(w, np.abs(v2 - O["vg"][i].astype(float)).max())
print("rotation invariance (vg, 200 models): max |feature change| %.1e A (float32 storage)" % w)
