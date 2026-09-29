"""VXH2 (the voxel pipeline, A12.3 / A12.4 rule H1): reliability of the hinge descriptor - the hinge's hard stop.

Rel_h = 1 - 2 mean d_hinge^2(crystal, own model) / mean d_hinge^2(background), background = the reference cluster-test procedure's 60,000 random receptor
pairs (default_rng(0)) over the N models; 95 % bootstrap CI over the benchmark pairs (2,000, seed 0). The same for the F1
(five-landmark) pose. H1: pass if the lower bound >= 0.50 (Koo & Li 2016 scale, applied by analogy to this ICC-type ratio).
Report only: the hinge through the VXV panel (E1, E2 against the vector primary, 55 receptors); ICC(2,1) per TRangle
parameter (two-way random, absolute agreement, single measure; Koo & Li 2016) with bootstrap CI and model SD / crystal SD
(compression toward the mean, cf. Bujotzek et al. 2015); d_hinge's own 1 % cut and error over it; Rel_h by release-year
bin, bound state and species. Writes out/VXH2_*.csv / .json; checks/VXH2_checks.csv. No state data.
usage: python pipeline/code/vxh2_reliability.py
"""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NBP = expected("n_benchmark_pairs")
EXP_NBG = expected("n_background_pairs")
from vxv2_panel import pspear

NB, SEED = 2000, 0


def dh(a, b):
    return np.sqrt(((a - b) ** 2).sum(-1).mean(-1))


def icc21(x, y):
    Y = np.column_stack([x, y]); n, k = Y.shape
    gm = Y.mean(); msr = k * ((Y.mean(1) - gm) ** 2).sum() / (n - 1); msc = n * ((Y.mean(0) - gm) ** 2).sum() / (k - 1)
    sse = ((Y - Y.mean(1, keepdims=True) - Y.mean(0, keepdims=True) + gm) ** 2).sum(); mse = sse / ((n - 1) * (k - 1))
    return float((msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n))


if __name__ == "__main__":
    C = Checks("VXH2", stop_on_fail=True)
    C.info("status", "A12 VXH2 hinge reliability; no state data read")
    P = np.load(os.path.join(OUT, "VXH1_pose.npz"), allow_pickle=True)
    names = [str(x) for x in P["names"]]; nr = int(P["n_repertoire"]); pos = {c: k for k, c in enumerate(names)}
    X0 = pd.read_csv(os.path.join(OUT, "VXH0_crystals.csv"))
    rng0 = np.random.default_rng(0)
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    ba, bb = rng0.integers(0, nr, 60000), rng0.integers(0, nr, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    C.add("background pairs (vec2 draw)", len(ba), len(ba) == EXP_NBG, "== %s" % EXP_NBG)
    ic = np.array([pos["crystal:%s" % e] for e in X0.entry]); im = np.array([pos["model:%s" % e] for e in X0.entry])
    res = {}
    for key in ("h", "h_F1"):
        Hh = P[key]
        own = dh(Hh[ic], Hh[im]); bg = dh(Hh[ba], Hh[bb]); mb = float((bg ** 2).mean())
        rel = 1 - 2 * float((own ** 2).mean()) / mb
        rng = np.random.default_rng(SEED)
        bs = [1 - 2 * float((own[rng.integers(0, len(own), len(own))] ** 2).mean()) / mb for _ in range(NB)]
        cut = float(np.percentile(bg, 1))
        res[key] = dict(rel=rel, ci=(float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))), cut=cut,
                        err=float(np.median(own) / cut), within=int((own <= cut).sum()), own=own, mb=mb)
        C.info("Rel %s [95 %% CI]" % ("hinge (anchors, primary)" if key == "h" else "hinge F1 (landmarks, sensitivity)"),
               "%.3f [%.3f, %.3f]; own median %.3f A; bg median %.3f A; 1 %% cut %.3f A; error/cut %.3f (%d/%s within)"
               % (rel, *res[key]["ci"], np.median(own), np.median(bg), cut, res[key]["err"], res[key]["within"], EXP_NBP))
    X0["own_dhinge"] = res["h"]["own"]
    srows = []
    for col in ("year_bin", "bound_pMHC", "species"):
        for v, g in X0.groupby(col):
            o = g.own_dhinge.values
            rng = np.random.default_rng(SEED)
            bs = [1 - 2 * float((o[rng.integers(0, len(o), len(o))] ** 2).mean()) / res["h"]["mb"] for _ in range(NB)]
            srows.append(dict(subset=col, value=str(v), n=len(o), rel=1 - 2 * float((o ** 2).mean()) / res["h"]["mb"],
                              ci_lo=float(np.percentile(bs, 2.5)), ci_hi=float(np.percentile(bs, 97.5))))
    SS = pd.DataFrame(srows); SS.to_csv(os.path.join(OUT, "VXH2_subsets.csv"), index=False)
    for _, r in SS.iterrows():
        C.info("Rel_h subset %s = %s (n %d)" % (r.subset, r.value, r.n), "%.3f [%.3f, %.3f]" % (r.rel, r.ci_lo, r.ci_hi))

    # the hinge through the VXV panel (report only)
    Z = np.load(os.path.join(OUT, "VXV1_descriptors.npz"), allow_pickle=True)
    ent = [str(e) for e in Z["entry"]]; n = len(ent); iu = np.triu_indices(n, 1)
    Hc = P["h"][[pos["crystal:%s" % e] for e in ent]]; Hm = P["h"][[pos["model:%s" % e] for e in ent]]
    Cm = np.sqrt(((Hc[:, None] - Hc[None]) ** 2).sum(-1).mean(-1)); Mm = np.sqrt(((Hm[:, None] - Hm[None]) ** 2).sum(-1).mean(-1))
    S = Z["S"]; Mv, Cv = Z["vec_vcori_geo__M"], Z["vec_vcori_geo__C"]
    e1 = pspear(Cm[iu], Mm[iu], (S[iu],)); e1p = pspear(Cm[iu], Mm[iu]); e2 = pspear(Mm[iu], Cm[iu], (Mv[iu], S[iu]))
    rng = np.random.default_rng(SEED); b1, b2 = [], []
    for _ in range(NB):
        r = rng.integers(0, n, n); a_, b_ = iu; k = r[a_] != r[b_]; ia, ib = r[a_][k], r[b_][k]
        b1.append(pspear(Cm[ia, ib], Mm[ia, ib], (S[ia, ib],))); b2.append(pspear(Mm[ia, ib], Cm[ia, ib], (Mv[ia, ib], S[ia, ib])))
    ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
    C.info("hinge VXV E1 partial | S [CI] (plain)", "%.3f [%.3f, %.3f] (%.3f)" % (e1, *ci(b1), e1p))
    C.info("hinge VXV E2 partial(M_h, C_h | M_vector, S) [CI]", "%.3f [%.3f, %.3f]" % (e2, *ci(b2)))

    # TRangle ICC(2,1), SD ratio (report only)
    T = pd.read_csv(os.path.join(OUT, "VXH1_trangle.csv")).set_index("name")
    trows = []
    for prm in ("BA", "BC1", "AC1", "BC2", "AC2", "dc"):
        xc = T.reindex(["crystal:%s" % e for e in X0.entry])[prm].values; xm = T.reindex(["model:%s" % e for e in X0.entry])[prm].values
        okk = np.isfinite(xc) & np.isfinite(xm); xc, xm = xc[okk], xm[okk]
        rng = np.random.default_rng(SEED); bs = []
        for _ in range(NB):
            k = rng.integers(0, len(xc), len(xc)); bs.append(icc21(xc[k], xm[k]))
        trows.append(dict(param=prm, n=int(okk.sum()), icc21=icc21(xc, xm), ci_lo=float(np.percentile(bs, 2.5)),
                          ci_hi=float(np.percentile(bs, 97.5)), sd_model_over_crystal=float(xm.std(ddof=1) / xc.std(ddof=1)),
                          mean_abs_err=float(np.abs(xc - xm).mean())))
    TT = pd.DataFrame(trows); TT.to_csv(os.path.join(OUT, "VXH2_trangle_icc.csv"), index=False)
    for _, r in TT.iterrows():
        C.info("TRangle %s ICC(2,1) [CI]; SD model/crystal; MAE (n %d)" % (r.param, r.n),
               "%.3f [%.3f, %.3f]; %.3f; %.3f" % (r.icc21, r.ci_lo, r.ci_hi, r.sd_model_over_crystal, r.mean_abs_err))

    H1 = res["h"]["ci"][0] >= 0.50
    C.info("RULE H1: lower 95 % bound of Rel_h >= 0.50", "%.3f -> %s" % (res["h"]["ci"][0], "PASS" if H1 else "FAIL"))
    save_json(dict(rel_h=res["h"]["rel"], rel_h_ci=res["h"]["ci"], rel_hF1=res["h_F1"]["rel"], rel_hF1_ci=res["h_F1"]["ci"],
                   cut_h=res["h"]["cut"], error_over_cut_h=res["h"]["err"], within_cut_h=res["h"]["within"],
                   H1_pass=bool(H1), vxv_E1=[e1, *ci(b1)], vxv_E2=[e2, *ci(b2)], seed=SEED, n_boot=NB),
              os.path.join(OUT, "VXH2_summary.json"))
    C.write()
