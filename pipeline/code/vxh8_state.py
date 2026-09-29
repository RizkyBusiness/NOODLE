"""VXH8 (the voxel pipeline, A12.3 / A13.4 / A13.12): state tests, run ONCE each for Arm B and Arm C, mirroring the reference method
(the reference method's own cluster-test and survivor-matching rules).

Per arm, on the M molecules with a state: purity over clusters >= 3; the 200-permutation null ladder with
rng = default_rng(0) re-seeded per arm - free, within mouse, within mouse + V pair (exactly the reference cluster-test procedure / the reference method) - then the
fourth rung A13.4 adds, within mouse + V pair + CDR3 length class, drawn after the first three so their stream is
unchanged; per-cluster one-sided binomial against the base frequency of the cluster's top state (clusters >= 2);
Benjamini-Hochberg over those clusters, survivors at q <= 0.15 (Benjamini & Hochberg 1995).
Overlap (the reference method rule 3.4 match()): for each survivor of the the reference method primary <reference arm>, its robust ones and the
robust N = 10 clusters, the arm's main cluster, fraction kept, and kept-and-significant (>= 75 % and q <= 0.15); and
the same for each of the arm's survivors in the the reference method primary.
Validation first: the same function on the the reference method primary's labels must reproduce its V2_arm_comparison exactly.
Writes out/VXH8_*; checks/VXH8_checks.csv. Arms are reported side by side; neither is preferred after the fact.
usage: python pipeline/code/vxh8_state.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NMOL = expected("n_molecules")
from vxpaths import CFG as _CFG2
_REF2 = _CFG2.get("reference_method", {})
_N10F = _REF2.get("n10_labels", "")   # secondary reference labels file
_N10C = _REF2.get("n10_column", "")   # its cluster-label column


NPERM = 200
# ---- reference method (config block "reference_method"): the established method this pipeline is compared with.
# Keys: dir, labels, label_column, cluster_tests, arm_comparison, survivor_sets {label: file}, name, n10_labels,
# n10_column. Set them to your own files, or leave the block out to skip the comparison.
from vxpaths import CFG as _CFG
_REF = _CFG.get("reference_method", {})
def _refp(key, default=""):
    return _REF.get(key, default)
A30 = VXP(_refp("dir", "reference/out/"))
REF_LABELS, REF_COL = _refp("labels"), _refp("label_column")
REF_TESTS, REF_CMP = _refp("cluster_tests"), _refp("arm_comparison")
REF_NAME = _refp("name", "reference method")


def purity(codes, states, ncl):
    tot = np.bincount(codes, minlength=ncl).astype(float); best = np.zeros(ncl)
    for s in range(states.max() + 1):
        best = np.maximum(best, np.bincount(codes[states == s], minlength=ncl))
    return float((best / tot).mean())


def state_tests(M, aname, fourth=True):
    """the reference method's cluster-test section, verbatim in logic; M has clone_id, state, donor, v_A_prot, v_B_prot, length_class, aname"""
    rng2 = np.random.default_rng(0)
    D = M[M.state.notna()].reset_index(drop=True)
    d = D[D[aname] >= 0]; s_ = d[aname].value_counts(); d = d[d[aname].isin(s_[s_ >= 3].index)]
    cl, _ = pd.factorize(d[aname]); st, _ = pd.factorize(d.state); ncl = cl.max() + 1
    obs = purity(cl, st, ncl)
    r = dict(arm=aname, clusters=ncl, molecules=len(d), observed=round(obs, 4)); nulls = {}
    rungs = [("unstratified", np.zeros(len(d), int)), ("mouse", pd.factorize(d.donor)[0]),
             ("mouse_V", pd.factorize(d.donor.astype(str) + "|" + d.v_A_prot.astype(str) + "|" + d.v_B_prot.astype(str))[0])]
    if fourth:
        rungs.append(("mouse_V_len", pd.factorize(d.donor.astype(str) + "|" + d.v_A_prot.astype(str) + "|" + d.v_B_prot.astype(str)
                                                  + "|" + d.length_class.astype(str))[0]))
    for lab, strata in rungs:
        order = np.argsort(strata, kind="stable"); gs, base = strata[order], st[order]; nl = np.empty(NPERM)
        for kk in range(NPERM):
            i2 = np.lexsort((rng2.random(len(gs)), gs)); s = np.empty_like(st); s[order] = base[i2]; nl[kk] = purity(cl, s, ncl)
        r["null_" + lab] = round(float(nl.mean()), 4); r["delta_" + lab] = round(float(obs - nl.mean()), 4)
        r["z_" + lab] = round(float((obs - nl.mean()) / max(nl.std(), 1e-9)), 1)
        nulls[lab] = dict(mean=float(nl.mean()), sd=float(nl.std()), min=float(nl.min()), max=float(nl.max()))
    d = D[D[aname] >= 0]; s_ = d[aname].value_counts(); d = d[d[aname].isin(s_[s_ >= 2].index)]
    basef = D.state.value_counts(normalize=True); crows = []
    for cid, g in d.groupby(aname):
        top = g.state.value_counts().index[0]; kk = int((g.state == top).sum())
        crows.append(dict(arm=aname, cluster=int(cid), members=len(g), top_state=top, top_n=kk, frac=round(kk / len(g), 3),
                          mice=g.donor.nunique(), v_pairs=g.v_A_prot.astype(str).add("|").add(g.v_B_prot.astype(str)).nunique(),
                          p=stats.binomtest(kk, len(g), float(basef[top]), alternative="greater").pvalue))
    K = pd.DataFrame(crows); pv = K.p.values; o = np.argsort(pv); q = np.empty(len(pv))
    q[o] = np.minimum.accumulate((pv[o] * len(pv) / np.arange(1, len(pv) + 1))[::-1])[::-1]
    K["q_bh"] = q
    return r, K.sort_values("q_bh"), nulls


def match(mem_ids, lab, tst):
    """the reference method's match rule: main cluster of mem in the labelling, fraction kept, its q, kept (>= 75 %) and significant"""
    b = lab.reindex(mem_ids); b = b[b >= 0]
    if not len(b):
        return -1, 0.0, np.nan, False
    main = int(b.value_counts().index[0]); fr = float(b.value_counts().iloc[0] / len(mem_ids))
    q = float(tst.loc[main, "q_bh"]) if main in tst.index else np.nan
    return main, fr, q, bool(fr >= 0.75 and q == q and q <= 0.15)


if __name__ == "__main__":
    C = Checks("VXH8", stop_on_fail=True)
    C.info("status", "A13.12 state tests, once per arm (Arm B, Arm C); first read of state data in this pipeline")
    ST = pd.read_csv(os.path.join(STB, "_receptor_states_all.csv.gz")).drop_duplicates("clone_key")[["clone_key", "state", "donor"]]

    # ------------------------------------------------------------ validation: reproduce the reference method's primary exactly
    AL = pd.read_csv(A30 + REF_LABELS)
    L4 = pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"), usecols=["clone_id", "length_class"])
    AL = AL.merge(L4, on="clone_id", how="left")
    ra, Ka, _ = state_tests(AL, REF_COL, fourth=False)
    ref = pd.read_csv(A30 + REF_CMP).iloc[0]
    diff = max(abs(float(ra[k]) - float(ref[k])) for k in ref.index if k != "arm")
    C.add("mirror reproduces %s primary V2_arm_comparison (all fields)" % REF_NAME, "max |diff| %.2e" % diff, diff < 1e-9, "exact")
    Kref = pd.read_csv(A30 + REF_TESTS)
    C.add("mirror reproduces %s primary survivors (q <= 0.15)" % REF_NAME, "%d vs %d" % ((Ka.q_bh <= 0.15).sum(), (Kref.q_bh <= 0.15).sum()),
          int((Ka.q_bh <= 0.15).sum()) == int((Kref.q_bh <= 0.15).sum()) and len(Ka) == len(Kref), "identical counts")
    _SS = _refp("survivor_sets", {}); A5 = pd.read_csv(A30 + _SS["survivors"]); R8 = pd.read_csv(A30 + _SS["secondary"])
    labA = AL.set_index("clone_id")[REF_COL]; tstA = Ka.set_index("cluster")
    _N10F, _N10C = _refp("n10_labels", ""), _refp("n10_column", "")
    N10 = pd.read_csv(VXP(_N10F), usecols=["clone_id", _N10C]).set_index("clone_id")[_N10C]
    sets = {"%s survivors (%d)" % (REF_NAME, len(A5)): [list(labA.index[labA == c]) for c in A5.cluster],
            "%s robust (%d)" % (REF_NAME, int(A5.robust_n30.sum())): [list(labA.index[labA == c]) for c in A5.cluster[A5.robust_n30]],
            "%s secondary robust (%d)" % (REF_NAME, len(R8)): [list(N10.index[N10 == c]) for c in R8.n10_cluster]}
    C.info("reference sets", "; ".join("%s: %d clusters" % (k, len(v)) for k, v in sets.items()))

    # ------------------------------------------------------------ the two arms, once each
    rows, overl = [], []
    for aname, f in (("prop_voxel_armB", "VXH6_molecule_labels.csv.gz"), ("prop_voxel_armC", "VXC2_molecule_labels.csv.gz")):
        M = pd.read_csv(os.path.join(OUT, f)).merge(ST, on="clone_key", how="left")
        C.add("%s: {:,} molecules, order as VX4".format(EXP_NMOL) % aname, len(M), len(M) == EXP_NMOL and
              list(M.clone_id) == list(pd.read_csv(os.path.join(OUT, "VX4_molecule_labels.csv.gz"), usecols=["clone_id"]).clone_id), "identical")
        r, K, nulls = state_tests(M, aname)
        K.to_csv(os.path.join(OUT, "VXH8_cluster_tests_%s.csv" % aname), index=False)
        r["clusters_tested"] = len(K); r["survivors_q015"] = int((K.q_bh <= 0.15).sum()); rows.append(r)
        C.info("%s: clusters >= 3 / molecules / observed purity" % aname, "%d / %d / %.4f" % (r["clusters"], r["molecules"], r["observed"]))
        for lv in ("unstratified", "mouse", "mouse_V", "mouse_V_len"):
            C.info("%s: rung %s - null mean (sd) / excess / z" % (aname, lv), "%.4f (%.4f) / %+.4f / %.1f"
                   % (nulls[lv]["mean"], nulls[lv]["sd"], r["delta_" + lv], r["z_" + lv]))
        C.info("%s: clusters tested (>= 2) / survivors q <= 0.15" % aname, "%d / %d" % (r["clusters_tested"], r["survivors_q015"]))
        lab = M.set_index("clone_id")[aname]; tst = K.set_index("cluster")
        for sname, groups in sets.items():
            res = [match(g, lab, tst) for g in groups]
            overl.append(dict(arm=aname, reference=sname, clusters=len(groups), recovered=sum(x[3] for x in res),
                              median_frac_kept=float(np.median([x[1] for x in res]))))
            C.info("%s: %s recovered (kept >= 75 %% and q <= 0.15)" % (aname, sname), "%d of %d (median kept %.2f)"
                   % (sum(x[3] for x in res), len(groups), np.median([x[1] for x in res])))
        sv = K[K.q_bh <= 0.15]
        res = [match(list(lab.index[lab == c]), labA, tstA) for c in sv.cluster]
        overl.append(dict(arm=aname, reference="own survivors recovered in %s primary" % REF_NAME, clusters=len(sv),
                          recovered=sum(x[3] for x in res), median_frac_kept=float(np.median([x[1] for x in res])) if res else np.nan))
        C.info("%s: own survivors recovered in the reference method" % aname, "%d of %d" % (sum(x[3] for x in res), len(sv)))
        if len(sv):
            C.info("%s: survivors - members / top state / frac / mice / V pairs / q" % aname,
                   " | ".join("%d:%d %s %.2f m%d v%d q%.3f" % (x.cluster, x.members, x.top_state, x.frac, x.mice, x.v_pairs, x.q_bh)
                              for x in sv.itertuples()))
    ref_row = dict(ra); ref_row["arm"] = "%s (%s, reproduced; 3 rungs)" % (REF_COL, REF_NAME)
    ref_row["clusters_tested"] = len(Ka); ref_row["survivors_q015"] = int((Ka.q_bh <= 0.15).sum())
    pd.DataFrame(rows + [ref_row]).to_csv(os.path.join(OUT, "VXH8_arm_comparison.csv"), index=False)
    pd.DataFrame(overl).to_csv(os.path.join(OUT, "VXH8_overlap.csv"), index=False)
    rows_j = [{k: (v.item() if hasattr(v, "item") else v) for k, v in r_.items()} for r_ in rows]    # numpy scalars -> JSON
    save_json(dict(seed=0, n_perm=NPERM, rungs=["unstratified", "mouse", "mouse_V", "mouse_V_len"], bh_q=0.15, arms=rows_j),
              os.path.join(OUT, "VXH8_summary.json"))
    C.write()
