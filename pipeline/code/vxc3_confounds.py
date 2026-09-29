"""VXC3 (the voxel pipeline, A13.2 / A13.8.2): confounds for Arm C, beside Arm B, the whole-molecule and F0-loops grids and the existing
<reference arm>; and the report-only junctional view. No state data.

Confounds on the reference cluster-test procedure's 60,000 background pairs: Spearman of distance with |delta CDR3 length| (alpha, beta, total);
length-matched vs mismatched medians; same- vs different-V-pair medians; share of single-V-pair clusters; ARI between
partitions (only clone_id and label columns of the existing arm are read).
Junctional view (report only): heavy atoms of the junctional residues (current method's definition, VXC0), Arm C frame /
boxes / channels; distances among the molecules in the view; cut from the reference cluster-test procedure's background pairs with both members in
the view; control pairs likewise; V-pair diagnostic. No clustering decision rests on it.
Writes out/VXC3_*; checks/VXC3_checks.csv.
usage: python pipeline/code/vxc3_confounds.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist, squareform
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column

from vx4_cluster import ari
from vxc_common import load_boxes, build_dataset, gram_distances

if __name__ == "__main__":
    C = Checks("VXC3", stop_on_fail=True)
    C.info("status", "A13 VXC3 confounds + junctional view (report only); no state data read")
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; n = len(ids); idx = {c: i for i, c in enumerate(ids)}
    rng = np.random.default_rng(0); ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
    ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz")).set_index("clone_id").reindex(ids)
    vpair = (ID.v_A_prot.astype(str) + "|" + ID.v_B_prot.astype(str)).values
    la, lb = R.cdr3len_A.values, R.cdr3len_B.values

    LC = pd.read_csv(os.path.join(OUT, "VXC2_molecule_labels.csv.gz"), usecols=["clone_id", "prop_voxel_armC"])
    LB = pd.read_csv(os.path.join(OUT, "VXH6_molecule_labels.csv.gz"), usecols=["clone_id", "prop_voxel_armB"])
    L4 = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"), usecols=["clone_id", "prop_voxel"])
    L6 = pd.read_csv(os.path.join(OUT, "VX6L_molecule_labels.csv.gz"), usecols=["clone_id", "prop_voxel_loops"])
    LE = pd.read_csv(VXP(_N10F), usecols=["clone_id", _N10C])
    M = LC.merge(LB, on="clone_id").merge(L4, on="clone_id").merge(L6, on="clone_id").merge(LE, on="clone_id", how="left")
    M[_N10C] = M[_N10C].fillna(-1).astype(int)
    C.add("partitions aligned on the {:,} molecules".format(EXP_NMOL), len(M), len(M) == EXP_NMOL, "== %s" % EXP_NMOL)
    rep = np.array([idx[c] for c in M.clone_id])

    def diag(name, Dx, lab, a, b):
        d = np.asarray(Dx[a, b], dtype=np.float64); dA, dB = np.abs(la[a] - la[b]), np.abs(lb[a] - lb[b]); sv = vpair[a] == vpair[b]
        o = dict(arm=name, pairs=len(d))
        for k, x in (("alpha", dA), ("beta", dB), ("total", dA + dB)):
            o["spearman_dlen_" + k] = float(spearmanr(d, x).correlation)
        m0 = (dA + dB) == 0
        o["median_len_matched"], o["median_len_mismatched"] = float(np.median(d[m0])), float(np.median(d[~m0]))
        o["median_same_Vpair"], o["median_diff_Vpair"] = float(np.median(d[sv])), float(np.median(d[~sv]))
        if lab is not None:
            vp = pd.Series(vpair[rep]); cl = pd.Series(lab); g = vp[cl >= 0].groupby(cl[cl >= 0]).nunique()
            o["clusters"] = int(len(g)); o["share_single_Vpair_clusters"] = float((g == 1).mean())
            o["ARI_vs_reference"] = ari(lab, M[_N10C].values)
        return o

    rows = []
    for name, f, col in (("Arm C: F1 CDR3, 7 ch", "VXC1_D_receptors.f32.npy", "prop_voxel_armC"),
                         ("Arm B: F1 loops, 7 ch", "VXH5_D_receptors.f32.npy", "prop_voxel_armB"),
                         ("whole molecule (VX4)", "VX3_D_receptors.f32.npy", "prop_voxel"),
                         ("F0 loops (VX6L)", "VX6L_D_receptors.f32.npy", "prop_voxel_loops")):
        rows.append(diag(name, np.load(os.path.join(TMP, f), mmap_mode="r"), M[col].values, ba, bb))
    le = M[_N10C].values; ge = pd.Series(vpair[rep])[le >= 0].groupby(pd.Series(le)[le >= 0]).nunique()
    rows.append(dict(arm="reference method (context)", clusters=int(len(ge)), share_single_Vpair_clusters=float((ge == 1).mean())))
    for a_, b_ in (("prop_voxel_armC", "prop_voxel_armB"), ("prop_voxel_armC", "prop_voxel"), ("prop_voxel_armC", "prop_voxel_loops")):
        C.info("ARI %s vs %s" % (a_, b_), "%.4f" % ari(M[a_].values, M[b_].values))

    # junctional view (report only)
    J = json.load(open(os.path.join(OUT, "VXC0_junctional.json")))
    jid = sorted(J["positions"], key=lambda c: idx[c]); jrows = np.array([idx[c] for c in jid])
    pos = {c: {tuple(p) for p in J["positions"][c]} for c in jid}
    G, na, nb, mass, occ, F = build_dataset(jrows, ids, load_boxes(), "whole", 7, os.path.join(TMP, "VXC3_junctional.f16.npy"), pos)
    C.add("junctional view: typing complete", nb, nb == 0, "== 0")
    rel = max(abs(float(G[k][occ].astype(np.float64).sum()) / na[k] - 1) for k in np.where(na > 0)[0][:200])
    C.add("junctional view: mass conservation, 200 molecules", "%.4f %%" % (100 * rel), rel < 5e-3, "< 0.5 %")
    C.info("junctional heavy atoms per molecule: median (min-max)", "%d (%d-%d)" % (np.median(na), na.min(), na.max()))
    DJ, st = gram_distances(G, F)
    C.add("junctional view: symmetric, no negative D^2", "%.1e / %.2e" % (st["asym"], st["min_off"]), st["asym"] < 1e-6 and st["min_off"] > -1e-6, "exact")
    sub = np.arange(min(200, len(jrows)))
    ref = pdist(np.asarray(G[sub], dtype=np.float64)); rr = float(np.max(np.abs(squareform(DJ[np.ix_(sub, sub)], checks=False) - ref) / np.maximum(ref, 1e-12)))
    C.add("junctional view: blocked Gram == pdist (200)", "%.2e" % rr, rr < 1e-4, "< 1e-4")
    jpos = -np.ones(n, int); jpos[jrows] = np.arange(len(jrows))
    kb = (jpos[ba] >= 0) & (jpos[bb] >= 0)
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b]); kp = (jpos[pa] >= 0) & (jpos[pb] >= 0)
    dbj = DJ[jpos[ba[kb]], jpos[bb[kb]]]; dpj = DJ[jpos[pa[kp]], jpos[pb[kp]]]
    cutj = round(float(np.percentile(dbj, 1)), 4); bgs = np.sort(dbj); pct = 100 * np.searchsorted(bgs, dpj, side="right") / len(bgs)
    C.info("junctional view: background pairs / control pairs in the view", "%d / %d" % (int(kb.sum()), int(kp.sum())))
    C.info("junctional view: cut / background median; controls within cut / median bg pct",
           "%.4f / %.4f; %.4f / %.4f %%" % (cutj, np.median(dbj), float((dpj <= cutj).mean()), float(np.median(pct))))
    oj = diag("junctional only (report)", DJ, None, jpos[ba[kb]], jpos[bb[kb]])
    # diag() indexes the length / V arrays by receptor, so recompute its receptor-level quantities on the view's pairs
    a_, b_ = ba[kb], bb[kb]; dA, dB = np.abs(la[a_] - la[b_]), np.abs(lb[a_] - lb[b_]); sv = vpair[a_] == vpair[b_]
    for k, x in (("alpha", dA), ("beta", dB), ("total", dA + dB)):
        oj["spearman_dlen_" + k] = float(spearmanr(dbj, x).correlation)
    m0 = (dA + dB) == 0
    oj.update(median_len_matched=float(np.median(dbj[m0])), median_len_mismatched=float(np.median(dbj[~m0])),
              median_same_Vpair=float(np.median(dbj[sv])), median_diff_Vpair=float(np.median(dbj[~sv])), cut=cutj,
              control_within_cut=float((dpj <= cutj).mean()))
    rows.append(oj)
    del G; os.remove(os.path.join(TMP, "VXC3_junctional.f16.npy"))
    DG = pd.DataFrame(rows); DG.to_csv(os.path.join(OUT, "VXC3_diagnostics.csv"), index=False)
    for _, o in DG.iterrows():
        C.info("diagnostics: %s" % o.arm, "; ".join("%s=%.4g" % (k, v) for k, v in o.items() if k != "arm" and pd.notna(v)))
    save_json(dict(rows=DG.to_dict("records")), os.path.join(OUT, "VXC3_summary.json"))
    C.write()
