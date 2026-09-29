"""VXH6 (the voxel pipeline, A12.3; A12.12): Arm B clustering - the grid alone, so the weight scan does not apply.

Exactly the VX4 procedure (the reference cluster-test procedure's pairs from default_rng(0): the control pairs, 60,000 background pairs over the N
receptors, 40,000 within-length-class draws; cut = 1st percentile rounded to 4 dp; complete linkage on the M
molecules; singletons dropped; ARI at cut x0.98 / x1.02 (Hubert & Arabie 1985)). Crystal-model error over the full cut is
reported on the VX5b F1-loops pairs (same grids) - not a gate for Arm B (A12.0). No state columns are read.
Writes out/VXH6_*; checks/VXH6_checks.csv.
usage: python pipeline/code/vxh6_cluster.py
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
EXP_NCP = expected("n_control_pairs")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NBG = expected("n_background_pairs")
EXP_NWC = expected("n_within_class_pairs")
from vx4_cluster import summ, labels_at, ari

ARM = "prop_voxel_armB"

if __name__ == "__main__":
    C = Checks("VXH6", stop_on_fail=True)
    C.info("status", "A12 VXH6 Arm B clustering (grid alone); no state data read")
    J = json.load(open(os.path.join(OUT, "VXH5_build.json")))
    D = np.load(os.path.join(ROOT, J["receptor_matrix"]), mmap_mode="r")
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; n = len(ids); idx = {c: i for i, c in enumerate(ids)}
    R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
    lclass = (R.cdr3_A.str.len().astype(str) + "-" + R.cdr3_B.str.len().astype(str)).values
    rng = np.random.default_rng(0)
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
    ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    lc = pd.Series(lclass); groups = {k: v.values for k, v in lc.groupby(lc).groups.items()}
    big = [k for k, v in groups.items() if len(v) >= 2]; w = np.array([len(groups[k]) for k in big], float); w /= w.sum()
    wa, wb = [], []
    for k in rng.choice(len(big), 40000, p=w):
        v = groups[big[k]]; i, j = rng.integers(0, len(v), 2)
        if i != j:
            wa.append(v[i]); wb.append(v[j])
    wa, wb = np.array(wa), np.array(wb)
    C.add("pair draw as vec2 / VX4", "%d / %d / %d" % (len(pa), len(ba), len(wa)), len(ba) == EXP_NBG and len(wa) == EXP_NWC,
          "%s / %s / %s" % (EXP_NCP, EXP_NBG, EXP_NWC))
    dp, db, dw = (np.asarray(D[x, y], dtype=np.float64) for x, y in ((pa, pb), (ba, bb), (wa, wb)))
    T = pd.DataFrame([summ(ARM, dp, db, dw)]); T.to_csv(os.path.join(OUT, "VXH6_thresholds.csv"), index=False)
    cut = float(T.bg_p1.iloc[0])
    bgs = np.sort(db); pct = 100 * np.searchsorted(bgs, dp, side="right") / len(bgs)
    C.info("cut / background median", "%.4f / %.4f" % (cut, T.bg_median.iloc[0]))
    C.info("control pairs: share within cut / median bg percentile (p90) / median d over cut",
           "%.4f / %.4f %% (%.3f %%) / %.4f" % (T.sens_at_bg_p1.iloc[0], np.median(pct), np.percentile(pct, 90), np.median(dp) / cut))
    pd.DataFrame({"clone_a": NP.clone_a, "clone_b": NP.clone_b, "d": dp, "over_cut": dp / cut, "bg_pct": pct}) \
      .to_csv(os.path.join(OUT, "VXH6_control_pairs.csv"), index=False)

    Mw = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"))       # molecule rows, the reference cluster-test procedure rule; no state columns
    rep = np.array([idx[c] for c in Mw.clone_id])
    C.add("molecule rows == VXH5 molecule rows", len(rep), list(rep) == J["molecule_rows"], "identical, %s" % EXP_NMOL)
    Dm = np.load(os.path.join(ROOT, J["molecule_condensed"]))
    Z = linkage(Dm, method="complete"); lab = labels_at(Z, cut, len(rep))
    Sq = squareform(Dm)
    worst = max(float(Sq[np.ix_(np.where(lab == c)[0], np.where(lab == c)[0])].max()) for c in range(lab.max() + 1)); del Sq
    C.add("every within-cluster pair within the cut", "max %.4f vs cut %.4f" % (worst, cut), worst <= cut + 1e-9, "<= cut")
    sz = pd.Series(lab[lab >= 0]).value_counts()
    C.info("clusters / molecules clustered / largest / clusters >= 3",
           "%d / %d / %d / %d" % (len(sz), int((lab >= 0).sum()), int(sz.max()), int((sz >= 3).sum())))
    for f_ in (0.98, 1.02):
        l2 = labels_at(Z, cut * f_, len(rep))
        C.info("stability: cut x%.2f -> clusters / ARI vs reported" % f_, "%d / %.4f" % (len(np.unique(l2[l2 >= 0])), ari(lab, l2)))
    out = Mw[["clone_id", "clone_key", "prot_key", "v_A_prot", "v_B_prot", "length_class"]].copy(); out[ARM] = lab
    out.to_csv(os.path.join(OUT, "VXH6_molecule_labels.csv.gz"), index=False)
    V5 = pd.read_csv(os.path.join(OUT, "VX5b_crystal_pairs.csv")); V5 = V5[(V5.frame == "F1") & (V5.atoms == "loops") & (V5.sigma == 2.0)]
    C.info("crystal-model error over the full cut, %s pairs (report only, not a gate)" % EXP_NBP,
           "%.4f (%d/%s within)" % (float(V5.d.median() / cut), int((V5.d <= cut).sum()), EXP_NBP))
    save_json(dict(arm=ARM, cut=cut, clusters=int(len(sz)), clustered=int((lab >= 0).sum()), largest=int(sz.max()),
                   clusters_ge3=int((sz >= 3).sum()), control_within_cut=float(T.sens_at_bg_p1.iloc[0]),
                   control_median_bg_pct=float(np.median(pct)), crystal_error_over_cut=float(V5.d.median() / cut)),
              os.path.join(OUT, "VXH6_summary.json"))
    C.write()
