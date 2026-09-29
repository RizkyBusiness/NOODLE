"""VXT1 (the voxel pipeline, A15): paired alpha-beta TCRdist (tcrdist3; Dash et al. 2017, Mayer-Blackwell et al. 2021) for all N
receptors, and an independent tcrdist3 partition of the M molecules with the voxel rule. Report only; no state data.

Input per receptor (landmark-file order): v_A, j_A, cdr3_A, v_B, j_B, cdr3_B from tables//_slim_receptors.csv.gz; genes as *01
alleles of the tcrdist3 mouse database. tcrdist3 defaults: CDR1, CDR2, CDR2.5 (pmhc) weight 1, fixed gap position; CDR3
weight 3, trimmed 3 (N) and 2 (C), best gap position; gap penalty 4; the tcrdist3 substitution distance (BLOSUM62-based,
0-4). Paired distance = alpha + beta.
Partition: cut = 1st percentile of the reference cluster-test procedure's 60,000 random receptor pairs (default_rng(0)), rounded to 4 dp as VX4; complete
linkage on the molecules (VX4 molecule rows); singletons dropped (vx4_cluster.labels_at).
Checks: every V gene in the database; tcrdist3 keeps the input order; zero diagonal, symmetric; 50 random pairs equal an
independent pure-Python implementation of the formula (exactly); the control pairs (one CDR3 amino acid apart) below
the random pairs' 1st percentile (median; share within the cut reported); every within-cluster pair within the cut.
Writes tmp/VXT1_D_receptors.i16.npy, out/VXT1_molecule_labels.csv.gz, out/VXT1_thresholds.csv, out/VXT1_summary.json;
checks/VXT1_checks.csv. Runs in the `tcrdist3` conda env. usage: python pipeline/code/vxt1_tcrdist.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, time
import numpy as np, pandas as pd
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
EXP_NCP = expected("n_control_pairs")
EXP_NBG = expected("n_background_pairs")
EXP_NWC = expected("n_within_class_pairs")
from vx4_cluster import summ, labels_at
import tcrdist
from tcrdist.repertoire import TCRrep
import pwseqdist as pw
from pwseqdist.matrices import tcr_dict_distance_matrix as SUB

if __name__ == "__main__":
    t0 = time.time()
    C = Checks("VXT1", stop_on_fail=True)
    C.info("status", "A15 VXT1 tcrdist3 %s (%s); no state data read" % (getattr(tcrdist, "__version__", "0.3"), os.path.dirname(tcrdist.__file__)))
    ids = [str(c) for c in np.load(LM1, allow_pickle=True)["clone_id"]]; n = len(ids); idx = {c: i for i, c in enumerate(ids)}
    R = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz")).set_index("clone_id").reindex(ids)
    C.add("receptor table complete for the {:,} LM1 receptors".format(EXP_NS), int(R[["v_A", "cdr3_A", "v_B", "cdr3_B"]].notna().all(1).sum()),
          bool(R[["v_A", "cdr3_A", "v_B", "cdr3_B"]].notna().all().all()), "== %s" % EXP_NS)
    db = pd.read_csv(os.path.join(os.path.dirname(tcrdist.__file__), "db", "alphabeta_gammadelta_db.tsv"), sep="\t")
    names = set(db[db.organism == "mouse"].id)
    miss = sorted({g for g in pd.concat([R.v_A, R.v_B, R.j_A, R.j_B]) if g + "*01" not in names})
    C.add("every V and J gene in the tcrdist3 mouse database (*01)", "%d missing %s" % (len(miss), miss[:5]), not miss, "0 missing")
    cell = pd.DataFrame(dict(clone_id=ids, v_a_gene=R.v_A.values + "*01", j_a_gene=R.j_A.values + "*01", cdr3_a_aa=R.cdr3_A.values,
                             v_b_gene=R.v_B.values + "*01", j_b_gene=R.j_B.values + "*01", cdr3_b_aa=R.cdr3_B.values, count=1))
    tr = TCRrep(cell_df=cell, organism="mouse", chains=["alpha", "beta"], db_file="alphabeta_gammadelta_db.tsv",
                deduplicate=False, compute_distances=False, store_all_cdr=False)
    tr.cpus = 8; tr.compute_distances()          # dense on purpose (N > tcrdist3's 10,000 safety default; ~0.4 GB)
    C.add("tcrdist3 keeps the receptor order", len(tr.clone_df), list(tr.clone_df.clone_id) == ids, "identical, %s" % EXP_NS)
    Dt = (tr.pw_alpha.astype(np.int32) + tr.pw_beta.astype(np.int32))
    C.info("distances", "%.1f min; range %d-%d" % ((time.time() - t0) / 60, int(Dt.min()), int(Dt.max())))
    C.add("zero diagonal, symmetric", "%d / %d" % (int(np.abs(np.diag(Dt)).max()), int(np.abs(Dt - Dt.T).max())),
          np.abs(np.diag(Dt)).max() == 0 and (Dt == Dt.T).all(), "exact")
    # independent implementation of the formula (pure Python), from the CDR sequences tcrdist3 looked up
    # substitution distance: the tcrdist3 table for the 20 amino acids; any other symbol (the IMGT gap '.' in the aligned
    # CDR1 / CDR2 / CDR2.5 strings, X) scores 0, exactly as pwseqdist maps unknown symbols (its last row / column is 0)
    sub = lambda x, y: SUB.get((x, y), 0)
    ngap = int(sum(tr.clone_df[c].str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]").sum() for c in tr.clone_df.columns if c.endswith("_aa") and not c.startswith("cdr3")))
    C.info("germline CDR strings containing a non-amino-acid symbol (gap '.', scored 0 by tcrdist3 0.3)", ngap)
    C.add("CDR3 strings contain only the 20 amino acids", "", not tr.clone_df[["cdr3_a_aa", "cdr3_b_aa"]].apply(lambda x: x.str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]")).any().any(), "true")
    def seqd(a, b, ntrim, ctrim, fixed):
        if len(a) == len(b):
            return sum(sub(a[i], b[i]) for i in range(ntrim, len(a) - ctrim))
        s = min(len(a), len(b)); best = None
        rng_ = [min(6, 3 + (s - 5) // 2)] if fixed else None
        if not fixed:
            lo_, hi_ = 5, s - 1 - 4
            while lo_ > hi_:
                lo_ -= 1; hi_ += 1
            rng_ = range(lo_, hi_ + 1)
        for g in rng_:
            d = sum(sub(a[i], b[i]) for i in range(ntrim, g)) + sum(sub(a[len(a) - 1 - i], b[len(b) - 1 - i]) for i in range(ctrim, s - g))
            best = d if best is None else min(best, d)
        return best + 4 * abs(len(a) - len(b))
    def tcrd(i, j):
        r1, r2 = tr.clone_df.iloc[i], tr.clone_df.iloc[j]; t = 0
        for ch in "ab":
            t += 3 * seqd(r1["cdr3_%s_aa" % ch], r2["cdr3_%s_aa" % ch], 3, 2, False)
            for k in ("cdr1", "cdr2", "pmhc"):
                t += seqd(r1["%s_%s_aa" % (k, ch)], r2["%s_%s_aa" % (k, ch)], 0, 0, True)
        return t
    prs = np.random.default_rng(161).integers(0, n, (50, 2))
    bad = sum(int(tcrd(i, j) != Dt[i, j]) for i, j in prs)
    C.add("50 random pairs == independent implementation of the formula", "%d mismatches" % bad, bad == 0, "exact")
    np.save(os.path.join(TMP, "VXT1_D_receptors.i16.npy"), Dt.astype(np.int16))
    # the reference cluster-test procedure's pairs and the voxel cut rule
    rng = np.random.default_rng(0)
    NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
    pa = np.array([idx[c] for c in NP.clone_a]); pb = np.array([idx[c] for c in NP.clone_b])
    ba, bb = rng.integers(0, n, 60000), rng.integers(0, n, 60000); ok = ba != bb; ba, bb = ba[ok], bb[ok]
    lclass = (R.cdr3_A.str.len().astype(str) + "-" + R.cdr3_B.str.len().astype(str)).values
    lc = pd.Series(lclass); groups = {k: v.values for k, v in lc.groupby(lc).groups.items()}
    big = [k for k, v in groups.items() if len(v) >= 2]; w = np.array([len(groups[k]) for k in big], float); w /= w.sum()
    wa, wb = [], []
    for k in rng.choice(len(big), 40000, p=w):
        v = groups[big[k]]; i, j = rng.integers(0, len(v), 2)
        if i != j:
            wa.append(v[i]); wb.append(v[j])
    wa, wb = np.array(wa), np.array(wb)
    C.add("pair draw as vec2 / VX4", "%d / %d / %d" % (len(pa), len(ba), len(wa)), len(ba) == EXP_NBG and len(wa) == EXP_NWC, "%s / %s / %s" % (EXP_NCP, EXP_NBG, EXP_NWC))
    dp, db_, dw = (Dt[x, y].astype(np.float64) for x, y in ((pa, pb), (ba, bb), (wa, wb)))
    T = pd.DataFrame([summ("prop_tcrdist3", dp, db_, dw)]); T.to_csv(os.path.join(OUT, "VXT1_thresholds.csv"), index=False)
    cut = float(T.bg_p1.iloc[0])
    C.info("cut (1st percentile of random pairs) / background median", "%.4f / %.4f" % (cut, T.bg_median.iloc[0]))
    C.add("control pairs (one CDR3 amino acid apart) far below random pairs: median < random 1st percentile", "median %.0f vs %.4f (background median %.0f)"
          % (np.median(dp), cut, np.median(db_)), np.median(dp) < cut, "median < cut")
    C.info("control pairs within the cut (reported)", "%.4f" % T.sens_at_bg_p1.iloc[0])
    Mw = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"))
    rep = np.array([idx[c] for c in Mw.clone_id])
    Dm = squareform(Dt[np.ix_(rep, rep)].astype(np.float64), checks=False)
    Z = linkage(Dm, method="complete"); lab = labels_at(Z, cut, len(rep))
    Sq = squareform(Dm)
    worst = max(float(Sq[np.ix_(np.where(lab == c)[0], np.where(lab == c)[0])].max()) for c in range(lab.max() + 1)); del Sq
    C.add("every within-cluster pair within the cut", "max %.0f vs cut %.4f" % (worst, cut), worst <= cut + 1e-9, "<= cut")
    sz = pd.Series(lab[lab >= 0]).value_counts()
    C.info("clusters / molecules clustered / largest / clusters >= 3", "%d / %d / %d / %d" % (len(sz), int((lab >= 0).sum()), int(sz.max()), int((sz >= 3).sum())))
    out = Mw[["clone_id", "clone_key", "v_A_prot", "v_B_prot", "length_class"]].copy(); out["prop_tcrdist3"] = lab
    out.to_csv(os.path.join(OUT, "VXT1_molecule_labels.csv.gz"), index=False)
    save_json(dict(cut=cut, clusters=int(len(sz)), clustered=int((lab >= 0).sum()), largest=int(sz.max()), clusters_ge3=int((sz >= 3).sum()),
                   control_within_cut=float(T.sens_at_bg_p1.iloc[0]), bg_median=float(T.bg_median.iloc[0]), tcrdist3_version="0.3",
                   matrix=os.path.relpath(VXP("voxel_out/tmp/VXT1_D_receptors.i16.npy"), ROOT), minutes=round((time.time() - t0) / 60, 1)), os.path.join(OUT, "VXT1_summary.json"))
    C.write()
