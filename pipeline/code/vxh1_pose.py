"""VXH1 (the voxel pipeline, A12.2 / A12.3): the Va/Vb pose (hinge) descriptor for every structure.

For structure i, with its framework anchors (the landmark file `anch`, 63 alpha / 60 beta; benchmark structures parsed as VX5b, fits on
the anchors present, the full reference set placed):
  1. fit i's alpha anchors onto the reference (the reference receptor (config)) alpha anchors (Kabsch 1976) and apply it to i's beta anchors;
  2. fit the REFERENCE beta anchors onto those placed beta anchors -> rigid pose P_i; place the reference beta anchors: X_i^b;
  3. symmetrically, beta as base, alpha placed: X_i^a;
  4. h_i = [X_i^b ; X_i^a] (123 x 3);  d_hinge(i, j) = sqrt(mean over the 123 points of |h_i - h_j|^2)  (A).
Reference anchors, not i's own, so V-gene framework shape drops out and only the rigid relative pose remains.
Also: rotation vector and angle of P_i (hinge angle to the reference), interdomain translation (displacement of the placed
reference-beta centroid), the same pose from the five landmarks per chain (F1 sensitivity), and TRangle via STCRpy in its
own env (A12.11). Writes out/VXH1_pose.npz, out/VXH1_trangle.csv; checks/VXH1_checks.csv. No state data.
usage: python pipeline/code/vxh1_pose.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, subprocess
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
import vx5b_frame_diagnostic as FB

STPY = __import__("vxpaths").env_python("stcrpy")   # config "envs"
IA = np.array([k[0] == "A" for k in FB.AKEYS])
REF_A, REF_B = FB.T_AN[IA], FB.T_AN[~IA]
LA, LB = np.arange(5), np.arange(5, 10)


def fit(X, Y):
    R, cx, cy = kabsch(X, Y)
    return R, cy - cx @ R


def pose_from(XA, XB, okA, okB, refA, refB):
    """placed reference anchors; okA/okB mark which of the chain's anchors are present"""
    R1, t1 = fit(XA[okA], refA[okA]); Yb = XB @ R1 + t1
    Rp, tp = fit(refB[okB], Yb[okB]); Pb = refB @ Rp + tp
    R2, t2 = fit(XB[okB], refB[okB]); Ya = XA @ R2 + t2
    Rq, tq = fit(refA[okA], Ya[okA]); Pa = refA @ Rq + tq
    return np.vstack([Pb, Pa]), Rp, float(np.linalg.norm(Pb.mean(0) - refB.mean(0)))


def rotvec(R):
    ang = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))
    if ang < 1e-12:
        return np.zeros(3), 0.0
    w = np.array([R[1, 2] - R[2, 1], R[2, 0] - R[0, 2], R[0, 1] - R[1, 0]]) / (2 * np.sin(ang))   # row convention
    return w * ang, float(np.degrees(ang))


def dh(a, b):
    return float(np.sqrt(((a - b) ** 2).sum(-1).mean(-1)))


if __name__ == "__main__":
    C = Checks("VXH1", stop_on_fail=True)
    C.info("status", "A12 VXH1 pose; no state data read")
    ids = FB.IDS; n = len(ids)
    lm3 = load_lm3(); B = lm3["B"].sort_values("entry").reset_index(drop=True)
    FB.PA.lm3_LMPOS = lm3["LMPOS"]
    names, H, HF1, RV, ANG, TR = [], [], [], [], [], []
    allA, allB = np.ones(IA.sum(), bool), np.ones((~IA).sum(), bool)
    for i in range(n):
        an = FB.ANC[i]; h, Rp, tr = pose_from(an[IA], an[~IA], allA, allB, REF_A, REF_B)
        hf, _, _ = pose_from(FB.LMK[i][LA], FB.LMK[i][LB], np.ones(5, bool), np.ones(5, bool), FB.T_LM[LA], FB.T_LM[LB])
        v, a = rotvec(Rp); names.append(ids[i]); H.append(h); HF1.append(hf); RV.append(v); ANG.append(a); TR.append(tr)
    bench = []
    for _, r in B.iterrows():
        for kind, p in (("crystal", os.path.join(FB.PA.BEN, "fixed_%s_crystal.pdb" % r.entry)), ("model", os.path.join(FB.PA.BEN, r.model))):
            lm_, an_, ok_ = FB.bench_points(p, lm3["parse"])
            h, Rp, tr = pose_from(an_[IA], an_[~IA], ok_[IA], ok_[~IA], REF_A, REF_B)
            hf, _, _ = pose_from(lm_[LA], lm_[LB], np.ones(5, bool), np.ones(5, bool), FB.T_LM[LA], FB.T_LM[LB])
            v, a = rotvec(Rp)
            names.append("%s:%s" % (kind, r.entry)); H.append(h); HF1.append(hf); RV.append(v); ANG.append(a); TR.append(tr)
            bench.append((kind, r.entry, r.model, p, Rp))
    H, HF1 = np.array(H), np.array(HF1)
    C.add("no NaN in pose", "", np.isfinite(H).all() and np.isfinite(HF1).all(), "finite")

    # reference maps to itself; invariance to a random rigid transform
    ri = ids.index(REF_ID)
    e = float(np.abs(H[ri] - np.vstack([REF_B, REF_A])).max())
    C.add("reference maps to itself", "%.2e A" % e, e < 1e-6, "< 1e-6 A")
    rng = np.random.default_rng(91); worst = 0.0
    for i in rng.choice(n, 5, replace=False):
        q = rng.normal(size=4); q /= np.linalg.norm(q); a, b, c, d = q
        Q = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                      [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])
        an = FB.ANC[i] @ Q + rng.uniform(-50, 50, 3)
        h2, _, _ = pose_from(an[IA], an[~IA], allA, allB, REF_A, REF_B)
        worst = max(worst, float(np.abs(h2 - H[i]).max()))
    C.add("pose invariant to a random rigid transform (5 models)", "%.2e A" % worst, worst < 1e-6, "< 1e-6 A")

    # crystal vs own model hinge angle reproduces VX5b (route via the reference: median |diff| <= 0.2 deg)
    V5 = pd.read_csv(os.path.join(OUT, "VX5b_hinge.csv")).set_index("entry")
    diffs = []
    for k in range(0, len(bench), 2):
        (_, e_, _, _, Rc), (_, _, _, _, Rm) = bench[k], bench[k + 1]
        ang = float(np.degrees(np.arccos(np.clip((np.trace(Rc.T @ Rm) - 1) / 2, -1, 1))))
        diffs.append(abs(ang - V5.loc[e_, "hinge_deg"]))
    C.add("crystal-vs-model hinge angle reproduces VX5b_hinge.csv", "median |diff| %.3f deg (max %.3f)" % (np.median(diffs), max(diffs)),
          np.median(diffs) <= 0.2, "median <= 0.2 deg")

    # 12 pairs: d_hinge equals an independent loop recomputation
    worst = 0.0
    for a, b in rng.choice(len(H), (12, 2), replace=False):
        s = 0.0
        for p_ in range(H.shape[1]):
            s += sum((H[a][p_][k] - H[b][p_][k]) ** 2 for k in range(3))
        worst = max(worst, abs(dh(H[a], H[b]) - (s / H.shape[1]) ** 0.5))
    C.add("d_hinge == independent loop (12 pairs)", "%.2e" % worst, worst < 1e-9, "< 1e-9")
    C.info("hinge angle to reference, repertoire: median / p95 / max (deg)",
           "%.2f / %.2f / %.2f" % (np.median(ANG[:n]), np.percentile(ANG[:n], 95), max(ANG[:n])))
    C.info("interdomain translation, repertoire: median / p95 (A)", "%.2f / %.2f" % (np.median(TR[:n]), np.percentile(TR[:n], 95)))

    # TRangle in the separate env
    lst = os.path.join(TMP, "VXH1_trangle_paths.txt")
    paths = [VXP("structures/%s.pdb") % c for c in ids] + [b[3] for b in bench]
    open(lst, "w").write("\n".join(paths))
    r = subprocess.run([STPY, os.path.join(os.path.dirname(__file__), "vxh1_trangle_stcrpy.py"), lst,
                        os.path.join(OUT, "VXH1_trangle.csv")], capture_output=True, text=True)
    os.remove(lst)
    C.add("STCRpy TRangle run completed", r.stdout.strip()[-120:] or r.stderr[-300:], r.returncode == 0, "exit 0")
    TG = pd.read_csv(os.path.join(OUT, "VXH1_trangle.csv"))
    TG["name"] = names
    TG["path"] = [os.path.relpath(p, ROOT) for p in TG.path]   # relative to project_root
    TG.to_csv(os.path.join(OUT, "VXH1_trangle.csv"), index=False)
    fail = TG[TG.error.notna() & (TG.error != "")]
    C.info("TRangle failures: repertoire / benchmark", "%d / %d (%s)" % (int(fail.name.isin(ids).sum()), int((~fail.name.isin(ids)).sum()),
           " ".join(fail.name[~fail.name.isin(ids)])[:300]))
    np.savez_compressed(os.path.join(OUT, "VXH1_pose.npz"), names=np.array(names), h=H, h_F1=HF1, rotvec=np.array(RV),
                        hinge_angle_to_ref_deg=np.array(ANG), interdomain_translation_A=np.array(TR), n_repertoire=n)
    C.write()
