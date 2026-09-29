"""VXT2 (the voxel pipeline, A15.3-4): for every cluster (>= 2 molecules) of Arms B, C, D, E and, for reference, the the reference method primary
(<reference arm>; only its clone_id and label columns are read), is it also a tcrdist3 cluster? Report only; no state data.

Verdict (A15.3, the reference method match() rule): main tcrdist3 cluster = the most common tcrdist3 label among the members (tcrdist3
singletons count for no cluster); frac = its members / cluster size; "also" if frac >= 0.75, "partly" 0.50-0.75, "not"
below. Descriptive (A15.4): median and maximum pairwise TCRdist of the members, share of member pairs within the tcrdist3
cut, random-pair percentile of the median (the reference cluster-test procedure's 60,000 pairs). Member x member TCRdist matrices kept for the reports.
Checks: partitions aligned on the M molecules; for 3 clusters per arm the verdict and pairwise figures recomputed by a
separate route (pairs enumerated from the receptor matrix by clone id); counts of verdicts add up.
Writes out/VXT2_crosscheck_<arm>.csv, out/VXT2_pairs_<arm>.npz, out/VXT2_summary.csv; checks/VXT2_checks.csv.
Runs in the `tcrdist3` (or `docking`) env. usage: python pipeline/code/vxt2_crosscheck.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
from vxpaths import CFG as _CFG
_REF = _CFG.get("reference_method", {})
_refp = lambda k, d="": _REF.get(k, d)
_RD, _REF_KEY = _refp("dir", "reference/out/"), _refp("name", "reference")

ARMS = {"B": ("VXH6_molecule_labels.csv.gz", "prop_voxel_armB", "VXH8_cluster_tests_prop_voxel_armB.csv"),
        "C": ("VXC2_molecule_labels.csv.gz", "prop_voxel_armC", "VXH8_cluster_tests_prop_voxel_armC.csv"),
        "D": ("VXS2_D_molecule_labels.csv.gz", "prop_voxel_armD", "VXS5_cluster_tests_prop_voxel_armD.csv"),
        "E": ("VXS2_E_molecule_labels.csv.gz", "prop_voxel_armE", "VXS5_cluster_tests_prop_voxel_armE.csv"),
        _REF_KEY: (VXP(_RD + _refp("labels")), _refp("label_column"), VXP(_RD + _refp("cluster_tests")))}
verdict = lambda fr: "also" if fr >= 0.75 else ("partly" if fr >= 0.5 else "not")

if __name__ == "__main__":
    C = Checks("VXT2", stop_on_fail=True)
    C.info("status", "A15 VXT2 tcrdist3 cross-check of every cluster; no state data read")
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; idx = {c: i for i, c in enumerate(ids)}; n = len(ids)
    Dt = np.load(os.path.join(TMP, "VXT1_D_receptors.i16.npy")).astype(np.int32)
    cut = float(pd.read_csv(os.path.join(OUT, "VXT1_thresholds.csv")).bg_p1.iloc[0])
    rng = np.random.default_rng(0); ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb
    bg = np.sort(Dt[ba[ok], bb[ok]])
    T = pd.read_csv(os.path.join(OUT, "VXT1_molecule_labels.csv.gz"))
    mids = list(T.clone_id); tl = T.prop_tcrdist3.values; rows_m = np.array([idx[c] for c in mids])
    summ_rows = []
    for arm, (lf, col, tf) in ARMS.items():
        L = pd.read_csv(lf if lf.startswith(("reference/", "/")) else os.path.join(OUT, lf), usecols=["clone_id", col])
        L = L.set_index("clone_id").reindex(mids)
        C.add("arm %s: labels aligned on the {:,} molecules".format(EXP_NMOL) % arm, int(L[col].notna().sum()), L[col].notna().all() and len(L) == EXP_NMOL,
              "%s" % EXP_NMOL)
        lab = L[col].fillna(-1).astype(int).values
        K = pd.read_csv(tf if tf.startswith(("reference/", "/")) else os.path.join(OUT, tf)).set_index("cluster")
        out, mats, offs = [], [], [0]
        for c in sorted(set(lab[lab >= 0])):
            mem = np.where(lab == c)[0]; k = len(mem)
            if k < 2:
                continue
            t = pd.Series(tl[mem]); t = t[t >= 0]
            main, cnt = (int(t.value_counts().index[0]), int(t.value_counts().iloc[0])) if len(t) else (-1, 0)
            fr = cnt / k
            Sm = Dt[np.ix_(rows_m[mem], rows_m[mem])]; iu = np.triu_indices(k, 1); pv = Sm[iu]
            med = float(np.median(pv))
            top = t.value_counts().head(4)
            out.append(dict(cluster=int(c), size=k, tcrdist3_main=main, tcrdist3_frac=round(fr, 4), verdict=verdict(fr),
                            tcrdist3_clusters=int(t.nunique()), tcrdist3_singletons=int(k - len(t)),
                            tcrdist3_top=";".join("%d:%d" % (a, b) for a, b in top.items()),
                            tcrdist_median=med, tcrdist_max=int(pv.max()), share_within_cut=float((pv <= cut).mean()),
                            median_bg_pct=float(100 * np.searchsorted(bg, med, side="right") / len(bg)),
                            q=float(K.loc[c, "q_bh"]) if c in K.index else np.nan))
            mats.append(Sm.astype(np.int16).ravel()); offs.append(offs[-1] + k * k)
        X = pd.DataFrame(out); X.to_csv(os.path.join(OUT, "VXT2_crosscheck_%s.csv" % arm), index=False)
        np.savez_compressed(os.path.join(OUT, "VXT2_pairs_%s.npz" % arm), cluster=X.cluster.values, offset=np.array(offs), d=np.concatenate(mats))
        # independent route for 3 clusters: pairs by clone id from the receptor matrix, verdict by explicit counting
        sub = [X.cluster.iloc[0], X.cluster.iloc[len(X) // 2], X.cluster.iloc[-1]]; badc = 0
        tmap = dict(zip(mids, tl))
        for c in sub:
            cl = list(L.index[L[col] == c]); ds = [int(Dt[idx[a], idx[b]]) for i_, a in enumerate(cl) for b in cl[i_ + 1:]]
            counts = {}
            for m in cl:
                if tmap[m] >= 0:
                    counts[tmap[m]] = counts.get(tmap[m], 0) + 1
            fr = max(counts.values()) / len(cl) if counts else 0.0
            r = X.set_index("cluster").loc[c]
            badc += int(abs(fr - r.tcrdist3_frac) > 1e-4 or verdict(fr) != r.verdict or float(np.median(ds)) != r.tcrdist_median or max(ds) != r.tcrdist_max)
        C.add("arm %s: 3 clusters recomputed by a separate route" % arm, "%d mismatches" % badc, badc == 0, "0")
        ge3 = X[X["size"] >= 3]; sv = X[X.q <= 0.15]
        o = dict(arm=arm, clusters=len(X), clusters_ge3=len(ge3), survivors=len(sv))
        for nm, S in (("all", X), ("ge3", ge3), ("surv", sv)):
            for v in ("also", "partly", "not"):
                o["%s_%s" % (nm, v)] = int((S.verdict == v).sum())
            o["%s_median_share_within" % nm] = float(S.share_within_cut.median()) if len(S) else np.nan
        C.add("arm %s: verdict counts add up" % arm, "%d / %d / %d" % (o["all_also"], o["all_partly"], o["all_not"]),
              o["all_also"] + o["all_partly"] + o["all_not"] == len(X), "== clusters")
        C.info("arm %s: clusters >= 3 also / partly / not a tcrdist3 cluster; survivors also / partly / not" % arm,
               "%d / %d / %d; %d / %d / %d" % (o["ge3_also"], o["ge3_partly"], o["ge3_not"], o["surv_also"], o["surv_partly"], o["surv_not"]))
        summ_rows.append(o)
    pd.DataFrame(summ_rows).to_csv(os.path.join(OUT, "VXT2_summary.csv"), index=False)
    C.write()
