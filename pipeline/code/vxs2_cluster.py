"""VXS2 (the voxel pipeline, A14.2.2): clustering of the sigma 1.5 A arms D and E - exactly the VXH6 / VXC2 procedure.

the reference cluster-test procedure's pairs from default_rng(0): the control pairs, 60,000 background pairs, 40,000 within-length-class draws; cut = 1st
percentile rounded to 4 dp; complete linkage on the M molecules; singletons dropped; ARI at cut x0.98 / x1.02.
The crystal-model pairs are rebuilt in the arm's configuration (shared atoms, vxv_common.bench_place, the arm's boxes,
7 channels) at sigma 1.5. The recipe is checked first at sigma 2.0: Arm D's configuration must reproduce VX5b's F1-loops
pair distances (1e-4 relative; VX5b placed the benchmark structures by its own fit code, which agrees to ~3e-5) and Arm
E's must reproduce VXC2's (1e-9, same code). The crystal-model error and the controls within cut are REPORTED beside the
sigma 2.0 arm; they are not gates (A14.3). No state columns are read.
Writes out/VXS2_{D,E}_*; checks/VXS2_{D,E}_checks.csv.
usage: python pipeline/code/vxs2_cluster.py D|E
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
from vxv_common import bench_place, build_grid, shared_keys, boxes_for, BEN, H
from vxc_common import load_boxes

ARMS = {"D": dict(col="prop_voxel_armD", atoms="loops", ref="B", ref_lab=("VXH6_molecule_labels.csv.gz", "prop_voxel_armB"), ref_sum="VXH6_summary.json"),
        "E": dict(col="prop_voxel_armE", atoms="cdr3", ref="C", ref_lab=("VXC2_molecule_labels.csv.gz", "prop_voxel_armC"), ref_sum="VXC2_summary.json")}


def pair_d(args):
    pc, pm, bx, atoms, sig, lm3 = args
    sh = shared_keys(pc, pm)
    gc = build_grid((pc, bench_place(pc, "F1", lm3), bx, atoms, 7, sh, None, False, sig))[0]
    gm = build_grid((pm, bench_place(pm, "F1", lm3), bx, atoms, 7, sh, None, False, sig))[0]
    return float(np.sqrt(((gc - gm) ** 2).sum() / H ** 3))


if __name__ == "__main__":
    arm = sys.argv[1]; A = ARMS[arm]; ARM = A["col"]
    C = Checks("VXS2_%s" % arm, stop_on_fail=True)
    C.info("status", "A14 VXS2 Arm %s clustering (sigma 1.5); no state data read" % arm)
    J = json.load(open(os.path.join(OUT, "VXS1_%s_build.json" % arm)))
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
    T = pd.DataFrame([summ(ARM, dp, db, dw)]); T.to_csv(os.path.join(OUT, "VXS2_%s_thresholds.csv" % arm), index=False)
    cut = float(T.bg_p1.iloc[0])
    bgs = np.sort(db); pct = 100 * np.searchsorted(bgs, dp, side="right") / len(bgs)
    C.info("cut / background median", "%.4f / %.4f" % (cut, T.bg_median.iloc[0]))
    C.info("control pairs: share within cut / median bg percentile (p90) / median d over cut",
           "%.4f / %.4f %% (%.3f %%) / %.4f" % (T.sens_at_bg_p1.iloc[0], np.median(pct), np.percentile(pct, 90), np.median(dp) / cut))
    pd.DataFrame({"clone_a": NP.clone_a, "clone_b": NP.clone_b, "d": dp, "over_cut": dp / cut, "bg_pct": pct}) \
      .to_csv(os.path.join(OUT, "VXS2_%s_control_pairs.csv" % arm), index=False)

    Mw = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"))
    rep = np.array([idx[c] for c in Mw.clone_id])
    C.add("molecule rows == VXS1 molecule rows", len(rep), list(rep) == J["molecule_rows"], "identical, %s" % EXP_NMOL)
    Dm = np.load(os.path.join(ROOT, J["molecule_condensed"]))
    Z = linkage(Dm, method="complete"); lab = labels_at(Z, cut, len(rep))
    Sq = squareform(Dm)
    worst = max(float(Sq[np.ix_(np.where(lab == c)[0], np.where(lab == c)[0])].max()) for c in range(lab.max() + 1)); del Sq
    C.add("every within-cluster pair within the cut", "max %.4f vs cut %.4f" % (worst, cut), worst <= cut + 1e-9, "<= cut")
    sz = pd.Series(lab[lab >= 0]).value_counts()
    C.info("clusters / molecules clustered / largest / clusters >= 3",
           "%d / %d / %d / %d" % (len(sz), int((lab >= 0).sum()), int(sz.max()), int((sz >= 3).sum())))
    stab = {}
    for f_ in (0.98, 1.02):
        l2 = labels_at(Z, cut * f_, len(rep)); stab[f_] = (int(len(np.unique(l2[l2 >= 0]))), ari(lab, l2))
        C.info("stability: cut x%.2f -> clusters / ARI vs reported" % f_, "%d / %.4f" % stab[f_])
    out = Mw[["clone_id", "clone_key", "prot_key", "v_A_prot", "v_B_prot", "length_class"]].copy(); out[ARM] = lab
    out.to_csv(os.path.join(OUT, "VXS2_%s_molecule_labels.csv.gz" % arm), index=False)
    LR = pd.read_csv(os.path.join(OUT, A["ref_lab"][0]), usecols=["clone_id", A["ref_lab"][1]])
    assert list(LR.clone_id) == list(out.clone_id)
    C.info("ARI vs Arm %s (same atoms, sigma 2.0)" % A["ref"], "%.4f" % ari(lab, LR[A["ref_lab"][1]].values))

    # the crystal-model pairs: recipe check at sigma 2.0, then sigma 1.5
    lm3 = load_lm3(); BX = boxes_for("F1") if arm == "D" else load_boxes()
    B = lm3["B"].sort_values("entry")
    jobs = lambda sig: [(os.path.join(BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(BEN, r.model), BX, A["atoms"], sig, lm3) for _, r in B.iterrows()]
    d20 = np.array([pair_d(j) for j in jobs(2.0)]); d15 = np.array([pair_d(j) for j in jobs(1.5)])   # serial: lm3 holds exec'd functions
    if arm == "D":
        V5 = pd.read_csv(os.path.join(OUT, "VX5b_crystal_pairs.csv")); V5 = V5[(V5.frame == "F1") & (V5.atoms == "loops") & (V5.sigma == 2.0)].set_index("entry")
        e = float(np.max(np.abs(d20 / V5.loc[list(B.entry), "d"].values - 1)))
        C.add("recipe at sigma 2.0 reproduces VX5b F1-loops crystal pairs (%s)" % EXP_NBP, "%.2e rel" % e, e < 1e-4, "< 1e-4")
    else:
        V2 = pd.read_csv(os.path.join(OUT, "VXC2_crystal_pairs.csv")).set_index("entry")
        e = float(np.max(np.abs(d20 / V2.loc[list(B.entry), "d"].values - 1)))
        C.add("recipe at sigma 2.0 reproduces VXC2 crystal pairs (%s)" % EXP_NBP, "%.2e rel" % e, e < 1e-9, "< 1e-9")
    pd.DataFrame(dict(entry=list(B.entry), model=list(B.model), d=d15, over_cut=d15 / cut, d_sigma20=d20)).to_csv(
        os.path.join(OUT, "VXS2_%s_crystal_pairs.csv" % arm), index=False)
    S0 = json.load(open(os.path.join(OUT, A["ref_sum"])))
    ce = float(np.median(d15) / cut)
    C.info("crystal-model error over the full cut, %s pairs (REPORTED, not a gate, A14.3)" % EXP_NBP, "%.4f (%d/%s within) vs Arm %s %.4f"
           % (ce, int((d15 <= cut).sum()), EXP_NBP, A["ref"], S0["crystal_error_over_cut"]))
    C.info("controls within cut (REPORTED, not a gate)", "%.4f vs Arm %s %.4f" % (T.sens_at_bg_p1.iloc[0], A["ref"], S0["control_within_cut"]))
    save_json(dict(arm=ARM, sigma=1.5, cut=cut, clusters=int(len(sz)), clustered=int((lab >= 0).sum()), largest=int(sz.max()),
                   clusters_ge3=int((sz >= 3).sum()), control_within_cut=float(T.sens_at_bg_p1.iloc[0]),
                   control_median_bg_pct=float(np.median(pct)), crystal_error_over_cut=ce, crystal_within_cut=int((d15 <= cut).sum()),
                   ari_vs_sigma20_arm=float(ari(lab, LR[A["ref_lab"][1]].values)), stability={str(k): v for k, v in stab.items()},
                   reference_arm=A["ref"], reference=dict(cut=S0["cut"], control_within_cut=S0["control_within_cut"],
                   crystal_error_over_cut=S0["crystal_error_over_cut"], clusters=S0["clusters"], clusters_ge3=S0["clusters_ge3"])),
              os.path.join(OUT, "VXS2_%s_summary.json" % arm))
    C.write()
