"""the voxel pipeline report, stage 5: gates on the assembled voxel_report.html (as the the reference method report's RG3 and self-review).

RG3 (data == sources): the JSON is parsed out of the HTML (not from the build files) and compared with <voxel out>/out/:
  cluster rows and survivors (VXH8_cluster_tests_*), summary state figures (VXH8_arm_comparison), cuts (thresholds), labels
  (VXH6 / VXC2 molecule labels), the split (stage 2 files), the worked example against the stored distance matrices.
Self-review (independent route): for 3 clusters per arm (largest survivor, a mid-size and a pair), recompute membership,
the with-state count, top state and share, mice, V pairs and the one-sided binomial p from the raw label and state tables
and the repertoire base rates, and compare with the JSON row and the VXH8 file (q).
Writes checks/RG3[_s15].csv, checks/SELFREVIEW[_s15].csv. Paths come from config/voxel_config.json.: python pipeline/report/code/rv5_gates.py [main|s15].
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", "..", "code"))
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import json, base64, sys, numpy as np, pandas as pd
from scipy import stats
from vxpaths import CFG as _CFGG
RP = VXP("voxel_out/report/"); O = VXP("voxel_out/out/"); rows = []
REPORT = sys.argv[1] if len(sys.argv) > 1 else "main"; assert REPORT in ("main", "s15", "s15v2")
HTML, SUF = {"main": ("voxel_report.html", ""), "s15": ("voxel_report_s15.html", "_s15"), "s15v2": ("voxel_report_s15_v2.html", "_s15v2")}[REPORT]
# per arm: label column, labels file, cluster tests, arm-comparison file, thresholds, distance matrix
SPEC = {"B": ("prop_voxel_armB", "VXH6_molecule_labels.csv.gz", "VXH8_cluster_tests_prop_voxel_armB.csv", "VXH8_arm_comparison.csv", "VXH6_thresholds.csv", "VXH5_D_receptors.f32.npy"),
        "C": ("prop_voxel_armC", "VXC2_molecule_labels.csv.gz", "VXH8_cluster_tests_prop_voxel_armC.csv", "VXH8_arm_comparison.csv", "VXC2_thresholds.csv", "VXC1_D_receptors.f32.npy"),
        "D": ("prop_voxel_armD", "VXS2_D_molecule_labels.csv.gz", "VXS5_cluster_tests_prop_voxel_armD.csv", "VXS5_arm_comparison.csv", "VXS2_D_thresholds.csv", "VXS1_D_D_receptors.f32.npy"),
        "E": ("prop_voxel_armE", "VXS2_E_molecule_labels.csv.gz", "VXS5_cluster_tests_prop_voxel_armE.csv", "VXS5_arm_comparison.csv", "VXS2_E_thresholds.csv", "VXS1_E_D_receptors.f32.npy")}
def chk(gate, name, ok, val=""):
    rows.append(dict(gate=gate, check=name, value=val, passed=bool(ok))); print("[%s] %-70s %s" % ("PASS" if ok else "FAIL", name, val))
h = open(RP + HTML).read()
i = h.index('<script id="data" type="application/json">') + len('<script id="data" type="application/json">')
D = json.loads(h[i:h.index("</script>", i)].replace("<\\/", "</"))
ST = pd.read_csv(VXP("tables/_receptor_states_all.csv.gz")).drop_duplicates("clone_key")[["clone_key", "state", "donor"]]
chk("RG3", "arms in the page == %s report's arms" % REPORT, list(D["arms"]) == (["B", "C"] if REPORT == "main" else ["D", "E", "B", "C"]), ",".join(D["arms"]))
for a in D["arms"]:
    c, lf, kf, sf, tf, dfn = SPEC[a]; S8 = pd.read_csv(O + sf).set_index("arm")
    L = pd.read_csv(O + lf); chk("RG3", "arm %s labels == %s" % (a, lf), list(L[c]) == D["mol"][a])
    K = pd.read_csv(O + kf).set_index("cluster")
    bad = 0
    for r in D["clusters"][a]:
        cid, n, ns, top, fr, mice, vp, p, q = r
        bad += int(n != int((L[c] == cid).sum()))
        if cid in K.index:
            k = K.loc[cid]; bad += int(ns != k.members or top != k.top_state or abs(fr - k.frac) > 1e-12 or mice != k.mice or vp != k.v_pairs or abs(q - k.q_bh) > 1e-15)
    chk("RG3", "arm %s: %d cluster rows == %s and label sizes" % (a, len(D["clusters"][a]), kf), bad == 0, "%d mismatches" % bad)
    sv = K[K.q_bh <= 0.15]
    chk("RG3", "arm %s survivors == %s (set and q)" % (a, kf), sorted(r["cluster"] for r in D["surv"][a]) == sorted(sv.index) and
        all(abs(r["q"] - sv.loc[r["cluster"], "q_bh"]) < 1e-15 for r in D["surv"][a]), "%d" % len(D["surv"][a]))
    sm = next(r for r in D["summary"] if r["arm"] == a); s8 = S8.loc[c]
    keys = [k for k in s8.index if k.startswith(("delta_", "z_", "null_"))]
    chk("RG3", "arm %s summary state figures == %s (%d fields)" % (a, sf, len(keys)), all(abs(sm[k] - float(s8[k])) < 1e-12 for k in keys))
    th = pd.read_csv(O + tf).iloc[0]
    chk("RG3", "arm %s cut == thresholds file" % a, abs(D["thr"][a]["bg_p1"] - float(th.bg_p1)) < 1e-12 and abs(sm["cut"] - float(th.bg_p1)) < 1e-12, "%.4f" % th.bg_p1)
    sp = np.frombuffer(base64.b64decode(D["split"][a]), "<f4").reshape(-1, 14); ref = np.load(RP + "work/split_%s.npz" % a)["split"]
    chk("RG3", "arm %s split in the HTML == stage-2 file" % a, np.allclose(np.nan_to_num(ref, nan=-1), sp, rtol=0, atol=1e-6))
    Dm = np.load(VXP("voxel_out/tmp/") + dfn, mmap_mode="r")
    ids = [str(x) for x in np.load(VXP("landmarks/out/LM1_internal_coords.npz"), allow_pickle=True)["clone_id"]]
    E = D["example"]["arms"][a]
    chk("RG3", "arm %s worked example D == stored matrix; parts sum to D^2" % a,
        abs(E["D"] - float(Dm[ids.index(D["example"]["a"]), ids.index(D["example"]["b"])])) < 1e-4 and abs(sum(E["parts"]) - E["D"] ** 2) < 1e-9, "%.4f" % E["D"])
    # ---- self-review: independent recomputation of three clusters
    M = L.merge(ST, on="clone_key", how="left"); base = M[M.state.notna()].state.value_counts(normalize=True)
    rows_c = {r[0]: r for r in D["clusters"][a]}
    pick = [D["surv"][a][0]["cluster"], sorted(rows_c, key=lambda x: abs(rows_c[x][1] - 8))[0], sorted(rows_c, key=lambda x: (rows_c[x][1] != 2, x))[0]]
    for cid in pick:
        g = M[M[c] == cid]; gs = g[g.state.notna()]; top = gs.state.value_counts().index[0]; kk = int((gs.state == top).sum())
        p = stats.binomtest(kk, len(gs), float(base[top]), alternative="greater").pvalue; r = rows_c[cid]
        ok = (r[1] == len(g) and r[2] == len(gs) and r[3] == top and abs(r[4] - round(kk / len(gs), 3)) < 1e-12 and r[5] == gs.donor.nunique()
              and r[6] == (gs.v_A_prot.astype(str) + "|" + gs.v_B_prot.astype(str)).nunique() and abs(r[7] - p) <= 1e-12 * max(1, p) and abs(r[8] - K.loc[cid, "q_bh"]) < 1e-15)
        chk("SELFREVIEW", "arm %s cluster %d recomputed independently (size %d)" % (a, cid, len(g)), ok, "%s %d/%d p %.3g" % (top, kk, len(gs), p))
# ---- A15 tcrdist3 block: read back against VXT1 / VXT2, and one cluster per arm recomputed from the receptor matrix
TL = pd.read_csv(O + "VXT1_molecule_labels.csv.gz").set_index("clone_id").reindex(D["mol"]["id"]).prop_tcrdist3.astype(int).tolist()
chk("RG3", "tcrdist3 labels in the page == VXT1_molecule_labels", D["mol"]["T"] == TL)
chk("RG3", "tcrdist3 cut == VXT1_thresholds", abs(D["tcrd"]["cut"] - float(pd.read_csv(O + "VXT1_thresholds.csv").bg_p1.iloc[0])) < 1e-12, "%.1f" % D["tcrd"]["cut"])
VS = pd.read_csv(O + "VXT2_summary.csv").set_index("arm")
Dt = np.load(VXP("voxel_out/tmp/VXT1_D_receptors.i16.npy"), mmap_mode="r"); IDS1 = [str(x) for x in np.load(VXP("landmarks/out/LM1_internal_coords.npz"), allow_pickle=True)["clone_id"]]
IX1 = {c: i for i, c in enumerate(IDS1)}; VC = {"also": 0, "partly": 1, "not": 2}
for a in D["arms"]:
    X = pd.read_csv(O + "VXT2_crosscheck_%s.csv" % a); Z = np.load(O + "VXT2_pairs_%s.npz" % a); R_ = D["tcrd"]["arms"][a]["rows"]
    ok = len(R_) == len(X) and all(r[0] == x.cluster and r[1] == x.tcrdist3_main and abs(r[2] - x.tcrdist3_frac) < 1e-12 and r[3] == VC[x.verdict]
                                   and abs(r[6] - x.tcrdist_median) < 1e-12 and r[7] == x.tcrdist_max and abs(r[8] - x.share_within_cut) < 1e-12 for r, x in zip(R_, X.itertuples()))
    chk("RG3", "arm %s: tcrdist3 rows == VXT2_crosscheck (%d clusters)" % (a, len(X)), ok)
    pp = np.frombuffer(base64.b64decode(D["tcrd"]["arms"][a]["d"]), "<i2")
    chk("RG3", "arm %s: tcrdist3 member matrices == VXT2_pairs" % a, np.array_equal(pp, Z["d"]) and D["tcrd"]["arms"][a]["off"] == [int(x) for x in Z["offset"]])
    sm = next(r for r in D["tcrd"]["summ"] if r["arm"] == a)
    chk("RG3", "arm %s: tcrdist3 summary == VXT2_summary" % a, all(sm[k] == VS.loc[a, k] for k in ("ge3_also", "ge3_partly", "ge3_not", "surv_also", "surv_partly", "surv_not")))
    # self-review: the largest survivor, recomputed from the receptor matrix and VXT1 labels
    cid = D["surv"][a][0]["cluster"]; mem = [c for c, l in zip(D["mol"]["id"], D["mol"][a]) if l == cid]
    tl = pd.Series([TL[D["mol"]["id"].index(c)] for c in mem]); tl = tl[tl >= 0]; fr = (tl.value_counts().iloc[0] if len(tl) else 0) / len(mem)
    ds = [int(Dt[IX1[x], IX1[y]]) for i_, x in enumerate(mem) for y in mem[i_ + 1:]]
    r = next(r for r in R_ if r[0] == cid)
    chk("SELFREVIEW", "arm %s cluster %d: tcrdist3 verdict and pairwise TCRdist recomputed" % (a, cid),
        abs(fr - r[2]) < 1e-4 and r[3] == (0 if fr >= 0.75 else 1 if fr >= 0.5 else 2) and float(np.median(ds)) == r[6] and max(ds) == r[7], "%s %.2f med %.0f" % (["also", "partly", "not"][r[3]], fr, np.median(ds)))
# ---- the counts the page text states (D.counts) == this run's files
Cn = D.get("counts", {}); L1c = np.load(VXP("landmarks/out/LM1_internal_coords.npz"), allow_pickle=True)
CPF = {"B": "VXH6_control_pairs.csv", "C": "VXC2_control_pairs.csv", "D": "VXS2_D_control_pairs.csv", "E": "VXS2_E_control_pairs.csv"}
chk("RG3", "counts: models == landmark file rows", Cn.get("n_models") == len(L1c["clone_id"]), "%s" % Cn.get("n_models"))
chk("RG3", "counts: anchors == landmark file anchors", Cn.get("n_anchors") == len(L1c["anch_keys"]), "%s" % Cn.get("n_anchors"))
chk("RG3", "counts: control pairs == every arm's control-pair file", all(Cn.get("n_control_pairs") == len(pd.read_csv(O + CPF[a])) for a in D["arms"]),
    "%s" % Cn.get("n_control_pairs"))
chk("RG3", "counts: crystal-model pairs == VX5_crystal_floor", Cn.get("n_benchmark_pairs") == len(pd.read_csv(O + "VX5_crystal_floor.csv")),
    "%s" % Cn.get("n_benchmark_pairs"))
chk("RG3", "counts: panel receptors == VXV2_rules n_receptors", Cn.get("n_panel_receptors") == json.load(open(O + "VXV2_rules.json"))["n_receptors"],
    "%s" % Cn.get("n_panel_receptors"))
# ---- side-car files and page wiring
nch = D["atoms"]["nchunks"]; import os
chk("RG3", "side-car chunks present (atoms, ca, aln) for all %d chunks" % nch, all(os.path.exists(RP + "report_data/%s_%03d.js" % (k, c)) for k in ("atoms", "ca", "aln") for c in range(nch)))
chk("RG3", "every clustered molecule (any arm) has a 3D chunk", all(D["atoms"]["chunk"][i] >= 0 for a in D["arms"] for i, l in enumerate(D["mol"][a]) if l >= 0))
if REPORT in ("s15", "s15v2"):                   # the A14 block == its source files
    P3 = pd.read_csv(O + "VXS3_panel.csv"); DF = pd.read_csv(O + "VXS3_sigma_diff.csv")
    num = lambda X: X.select_dtypes("number")
    chk("RG3", "s15 validity panel (8 grids) == VXS3_panel.csv", np.allclose(num(pd.DataFrame(D["s15"]["panel"])[P3.columns]).values, num(P3).values, rtol=0, atol=1e-12))
    chk("RG3", "s15 paired sigma differences == VXS3_sigma_diff.csv", np.allclose(num(pd.DataFrame(D["s15"]["diff"])[DF.columns]).values, num(DF).values, rtol=0, atol=1e-12))
    ok = all(abs(D["s15"]["clus"][a][k] - json.load(open(O + "VXS2_%s_summary.json" % a))[k]) < 1e-12 for a in "DE"
             for k in ("cut", "control_within_cut", "crystal_error_over_cut", "ari_vs_sigma20_arm"))
    chk("RG3", "s15 clustering summaries == VXS2 summaries", ok)
    sm = {r["arm"]: r for r in D["summary"]}
    chk("RG3", "s15 crystal error / cut in the arm tiles == VXS2 (D, E)", all(abs(sm[a]["crystal_over_cut"] - D["s15"]["clus"][a]["crystal_error_over_cut"]) < 1e-12 for a in "DE"))
if REPORT == "s15v2":
    order = [h[h.index('href="', i) + 6:h.index('"', h.index('href="', i) + 6)] for i in [m for m in range(len(h)) if h.startswith('<a href="#', m) and h.rfind('<nav id="nav">', 0, m) > h.rfind('</nav>', 0, m)]]
    chk("RG3", "s15v2 tab order Overview, UMAP, Clusters, Methods, Validation, Worked example", order == ["#overview", "#umap", "#clusters", "#methods", "#validation", "#example"], " ".join(order))
    ov = h[h.index("function viewOverview(app)"):h.index("\n}\n", h.index("function viewOverview(app)"))]
    chk("RG3", "s15v2 overview has no sigma-comparison table (heading and per-arm rows of the s15 overview absent)",
        '<h2 style="margin-top:0">σ 1.5 beside σ 2.0</h2>' not in h and "Blur σ (Å)" not in ov and "viewOverview" in h)
for nav in ("#overview", "#methods", "#validation", "#example", "#umap", "#clusters"):
    chk("RG3", "navigation link %s" % nav, 'href="%s"' % nav in h)
stale = [s for s in _CFGG.get("stale_markers", ["PRIM]"]) if s in h]
chk("RG3", "no stale single-method wiring from the source template", not stale, ", ".join(stale))
pd.DataFrame([r for r in rows if r["gate"] == "RG3"]).to_csv(RP + "checks/RG3%s.csv" % SUF, index=False)
pd.DataFrame([r for r in rows if r["gate"] == "SELFREVIEW"]).to_csv(RP + "checks/SELFREVIEW%s.csv" % SUF, index=False)
nf = sum(not r["passed"] for r in rows); print("%d checks, %d failed" % (len(rows), nf))
if nf: raise SystemExit("GATE FAILED")
