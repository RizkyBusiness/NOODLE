"""VXV2 (the voxel pipeline, A12.1 / A12.4): matrices, endpoints and the V1-V3 rules for the descriptor validity panel.

Per descriptor k: C_k (crystal-crystal) and M_k (model-model) over the receptor pairs; S the sequence baseline.
  E1  partial Spearman(C_k, M_k | S); plain Spearman alongside (Mantel p, one-sided, 9,999 label permutations, seed 0)
  E2  partial Spearman(M_grid, C_grid | M_vector, S), and the mirror (primaries only)
  E3  Spearman(M_k, C_l) for every k, l (report only; Mantel p)
  E4  Rel_k = 1 - 2 mean d_k^2(crystal, own model) / mean d_k^2(pilot background pairs)   (A12.8.1; report only)
  plus error over cut (median own distance / pilot 1 % cut) and the count within the cut.
95 % intervals: bootstrap over receptors (2,000 resamples, seed 0; pairs rebuilt from the resampled receptors, pairs of a
receptor with itself dropped), paired between descriptors for Delta. Partial correlations get bootstrap intervals only.
Rules (A12.4): V1 Delta = E1_grid - E1_vector with paired CI, and E2_grid with CI; V2 reinstated if lower CI(Delta) > -0.05
or lower CI(E2_grid) > 0; V3 flag a primary if lower CI(E1) <= 0. Secondary descriptors cannot trigger V2 or V3.
Writes out/VXV2_matrices.npz, VXV2_endpoints.csv, VXV2_crossk.csv, VXV2_strata.csv, VXV2_rules.json; checks/VXV2_checks.csv.
Spearman (rank) partial correlation: Pearson correlation of rank residuals after regressing on the ranks of the controls.
Mantel (1967) matrix permutation test. No state data.
usage: python pipeline/code/vxv2_panel.py
"""
import os, sys, json, time
import numpy as np, pandas as pd
from scipy.stats import rankdata
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
from vxv_common import PRIMARY_GRID, PRIMARY_VEC, CONFIGS, VECTORS

NBOOT, NPERM, SEED = 2000, 9999, 0


def pr(x, y):
    x = x - x.mean(); y = y - y.mean()
    return float((x * y).sum() / np.sqrt((x * x).sum() * (y * y).sum()))


def resid(r, Z):
    A = np.column_stack([np.ones(len(r))] + Z)
    return r - A @ np.linalg.lstsq(A, r, rcond=None)[0]


def pspear(x, y, Z=()):
    rx, ry = rankdata(x), rankdata(y)
    if not Z:
        return pr(rx, ry)
    Zr = [rankdata(z) for z in Z]
    return pr(resid(rx, Zr), resid(ry, Zr))


if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VXV2", stop_on_fail=True)
    C.info("status", "A12 VXV2 panel; no state data read")
    Z = np.load(os.path.join(OUT, "VXV1_descriptors.npz"), allow_pickle=True)
    names = [str(x) for x in Z["names"]]
    C.add("both primaries present", "%s, %s" % (PRIMARY_GRID, PRIMARY_VEC), PRIMARY_GRID in names and PRIMARY_VEC in names, "present")
    n = len(Z["entry"]); iu = np.triu_indices(n, 1)
    S = Z["S"]
    D = {k: dict(C=Z[k + "__C"], M=Z[k + "__M"], own=Z[k + "__own"], cut=float(Z[k + "__cut"]),
                 bg=float(Z[k + "__bg_mean_d2"])) for k in names}
    ok = all(np.allclose(v["C"], v["C"].T) and np.allclose(v["M"], v["M"].T) and not np.diag(v["C"]).any()
             and not np.diag(v["M"]).any() for v in D.values()) and np.allclose(S, S.T)
    C.add("all matrices symmetric with zero diagonal", "%d descriptors + S" % len(D), ok, "exact")
    C.info("receptors / pairs", "%d / %d" % (n, len(iu[0])))
    np.savez(os.path.join(OUT, "VXV2_matrices.npz"), **{k: Z[k] for k in Z.files})

    # ------------------------------------------------------------ point estimates
    up = lambda A, a=iu[0], b=iu[1]: A[a, b]
    s_ = up(S)
    rows = []
    for k in names:
        c, m = up(D[k]["C"]), up(D[k]["M"])
        own = D[k]["own"]
        rows.append(dict(descriptor=k, primary=k in (PRIMARY_GRID, PRIMARY_VEC), E1_partial=pspear(c, m, (s_,)),
                         E1_plain=pspear(c, m), E4_rel=1 - 2 * float((own ** 2).mean()) / D[k]["bg"],
                         error_over_cut=float(np.median(own) / D[k]["cut"]), within_cut=int((own <= D[k]["cut"]).sum()),
                         pilot_cut=D[k]["cut"]))
    E = pd.DataFrame(rows).set_index("descriptor")
    cg, mg = up(D[PRIMARY_GRID]["C"]), up(D[PRIMARY_GRID]["M"]); cv, mv = up(D[PRIMARY_VEC]["C"]), up(D[PRIMARY_VEC]["M"])
    E2g, E2v = pspear(mg, cg, (mv, s_)), pspear(mv, cv, (mg, s_))
    E3 = pd.DataFrame([[pspear(up(D[k]["M"]), up(D[l]["C"])) for l in names] for k in names], index=names, columns=names)

    # ------------------------------------------------------------ bootstrap over receptors
    rng = np.random.default_rng(SEED)
    B = {k: [] for k in names}; B4 = {k: [] for k in names}; Bd, Bg, Bv = [], [], []
    for b in range(NBOOT):
        r = rng.integers(0, n, n)
        a_, b_ = np.triu_indices(n, 1)
        keep = r[a_] != r[b_]
        ia, ib = r[a_][keep], r[b_][keep]
        sb = S[ia, ib]
        for k in names:
            B[k].append(pspear(D[k]["C"][ia, ib], D[k]["M"][ia, ib], (sb,)))
            B4[k].append(1 - 2 * float((D[k]["own"][r] ** 2).mean()) / D[k]["bg"])
        Bd.append(B[PRIMARY_GRID][-1] - B[PRIMARY_VEC][-1])
        cgb, mgb = D[PRIMARY_GRID]["C"][ia, ib], D[PRIMARY_GRID]["M"][ia, ib]
        cvb, mvb = D[PRIMARY_VEC]["C"][ia, ib], D[PRIMARY_VEC]["M"][ia, ib]
        Bg.append(pspear(mgb, cgb, (mvb, sb))); Bv.append(pspear(mvb, cvb, (mgb, sb)))
    ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    for k in names:
        E.loc[k, "E1_lo"], E.loc[k, "E1_hi"] = ci(B[k]); E.loc[k, "E4_lo"], E.loc[k, "E4_hi"] = ci(B4[k])
    C.add("bootstrap resamples all finite", NBOOT, np.isfinite(Bd).all() and np.isfinite(Bg).all(), "all %d" % NBOOT)

    # ------------------------------------------------------------ Mantel permutations (plain correlations only)
    Rc = np.array([rankdata(up(D[k]["C"])) for k in names]); Rc = (Rc - Rc.mean(1, keepdims=True))
    Rc /= np.sqrt((Rc ** 2).sum(1, keepdims=True))
    RM = [np.zeros((n, n)) for _ in names]
    for q, k in enumerate(names):
        RM[q][iu] = rankdata(up(D[k]["M"])); RM[q] = RM[q] + RM[q].T
    obs = np.array([[E3.loc[k, l] for l in names] for k in names])
    ge = np.zeros_like(obs)
    rngp = np.random.default_rng(SEED)
    for _ in range(NPERM):
        p = rngp.permutation(n)
        Mp = np.array([RM[q][np.ix_(p, p)][iu] for q in range(len(names))])
        Mp = Mp - Mp.mean(1, keepdims=True); Mp /= np.sqrt((Mp ** 2).sum(1, keepdims=True))
        ge += (Mp @ Rc.T) >= obs - 1e-12
    Pm = (1 + ge) / (1 + NPERM)
    for q, k in enumerate(names):
        E.loc[k, "E1_plain_mantel_p"] = Pm[q, q]
    X3 = E3.copy(); X3.index.name = "M_k \\ C_l"
    P3 = pd.DataFrame(Pm, index=names, columns=names)
    pd.concat({"spearman": X3, "mantel_p": P3}).to_csv(os.path.join(OUT, "VXV2_crossk.csv"))
    E.to_csv(os.path.join(OUT, "VXV2_endpoints.csv"))

    # ------------------------------------------------------------ strata (report only)
    sp, bd, VID = Z["species"], Z["bound"], Z["vregion_identity"]
    a_, b_ = iu
    mouse = np.array([s == "Mus musculus" for s in sp])
    strata = {"both mouse": mouse[a_] & mouse[b_], "any hybrid human+mouse": ~(mouse[a_] & mouse[b_]),
              "both bound": bd[a_] & bd[b_], "both unbound": ~bd[a_] & ~bd[b_], "bound x unbound": bd[a_] != bd[b_],
              "same V pair (identity >= 0.95 both chains)": (VID[a_, b_] >= 0.95).all(1),
              "different V pair": ~(VID[a_, b_] >= 0.95).all(1)}
    srows = []
    for nm, msk in strata.items():
        for k in names:
            if msk.sum() >= 10:
                srows.append(dict(stratum=nm, pairs=int(msk.sum()), descriptor=k,
                                  E1_plain=pspear(up(D[k]["C"])[msk], up(D[k]["M"])[msk]),
                                  E1_partial=pspear(up(D[k]["C"])[msk], up(D[k]["M"])[msk], (s_[msk],))))
            else:
                srows.append(dict(stratum=nm, pairs=int(msk.sum()), descriptor=k))
    pd.DataFrame(srows).to_csv(os.path.join(OUT, "VXV2_strata.csv"), index=False)

    # ------------------------------------------------------------ rules
    d_pt, d_ci = float(E.loc[PRIMARY_GRID, "E1_partial"] - E.loc[PRIMARY_VEC, "E1_partial"]), ci(Bd)
    g_ci, v_ci = ci(Bg), ci(Bv)
    V2 = bool(d_ci[0] > -0.05 or g_ci[0] > 0)
    V3 = {k: bool(E.loc[k, "E1_lo"] <= 0) for k in (PRIMARY_GRID, PRIMARY_VEC)}
    for k in names:
        C.info("E1 %s" % k, "partial %.3f [%.3f, %.3f]; plain %.3f (Mantel p %.4f); E4 %.3f [%.3f, %.3f]; error/cut %.3f (%d/%d within)"
               % (E.loc[k, "E1_partial"], E.loc[k, "E1_lo"], E.loc[k, "E1_hi"], E.loc[k, "E1_plain"],
                  E.loc[k, "E1_plain_mantel_p"], E.loc[k, "E4_rel"], E.loc[k, "E4_lo"], E.loc[k, "E4_hi"],
                  E.loc[k, "error_over_cut"], E.loc[k, "within_cut"], n))
    C.info("V1: Delta = E1_grid - E1_vector [paired 95 % CI]", "%.3f [%.3f, %.3f]" % (d_pt, *d_ci))
    C.info("V1: E2_grid = partial(M_grid, C_grid | M_vector, S)", "%.3f [%.3f, %.3f]" % (E2g, *g_ci))
    C.info("E2 mirror: partial(M_vector, C_vector | M_grid, S)", "%.3f [%.3f, %.3f]" % (E2v, *v_ci))
    C.info("V2: grid reinstated (lower CI Delta > -0.05 OR lower CI E2_grid > 0)", "%s (%.3f > -0.05: %s; %.3f > 0: %s)"
           % (V2, d_ci[0], d_ci[0] > -0.05, g_ci[0], g_ci[0] > 0))
    C.info("V3: flagged (lower CI of E1 <= 0)", "; ".join("%s: %s" % kv for kv in V3.items()))
    save_json(dict(n_receptors=int(n), n_pairs=int(len(iu[0])), seeds=dict(bootstrap=SEED, mantel=SEED),
                   n_boot=NBOOT, n_perm=NPERM, primary_grid=PRIMARY_GRID, primary_vector=PRIMARY_VEC,
                   E1={k: dict(point=float(E.loc[k, "E1_partial"]), ci=[float(E.loc[k, "E1_lo"]), float(E.loc[k, "E1_hi"])]) for k in names},
                   delta=dict(point=d_pt, ci=list(d_ci)), E2_grid=dict(point=E2g, ci=list(g_ci)),
                   E2_vector=dict(point=E2v, ci=list(v_ci)), V2_reinstated=V2, V3_flag=V3,
                   minutes=round((time.time() - t0) / 60, 1)), os.path.join(OUT, "VXV2_rules.json"))
    C.write()
