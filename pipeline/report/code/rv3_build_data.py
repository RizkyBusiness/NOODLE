"""the voxel pipeline report, stage 3: the data JSON of a report. Configurations (argument, default main):
  main  voxel_report.html      arms B (all four loops) and C (CDR3 only), sigma 2.0 A
  s15   voxel_report_s15.html  arms D and E (the same atoms at sigma 1.5 A, A14) beside B and C

Method-independent parts come from rv0_base_data.py (<voxel out>/report/work/base_data.json, built from primary files): the molecule
table (ids, genes, CDR3s, V proteins, state, mouse, cells), the UMAP cells, transcriptional clusters, gene index, and the
CDR3 germline-origin map. Everything method-specific is read from <voxel out>/out/ (and the stage 1-2 work files):
cluster labels and tests (VXH8), thresholds, stability (VXH6 / VXC2 checks), the exact distance split (stage 2),
3D chunk index (stage 1), the validity panel (VXV), hinge (VXH, exploratory), first design (VX2b, VX5, VX5b), comparison
and predictions (VXH9). Every number that reaches the page is asserted against its source file here.
For s15 also the A14 results (VXS1-VXS5). Writes the voxel pipeline/report/work/voxel_html_data.json (main) or
voxel_s15_html_data.json (s15). Paths come from config/voxel_config.json.: python pipeline/report/code/rv3_build_data.py [main|s15].
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", "..", "code"))
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, base64, re
import numpy as np, pandas as pd
sys.path.insert(0, VXP("voxel_out/code"))
from vxv_common import boxes_for, rep_place, build_grid, H
from vxpaths import CFG as _VXCFG
_REFCOL = _VXCFG.get("reference_method", {}).get("label_column", "")
from vxc_common import load_boxes
_RCFG = _VXCFG.get("report", {})

O = VXP("voxel_out/out/"); RP = VXP("voxel_out/report/"); W = RP + "work/"
# method-independent parts (molecule table, cell map, transcriptional clusters, gene index, CDR3 origin), from rv0_base_data.py
A30 = json.load(open(VXP("voxel_out/report/work/base_data.json")))
REPORT = sys.argv[1] if len(sys.argv) > 1 else "main"; assert REPORT in ("main", "s15")
LOOPS_NAME = "voxel grid of all four loops (CDR1, CDR2, HV4, CDR3), per-chain frames, 7 channels"
CDR3_NAME = "voxel grid of the two CDR3 loops, per-chain frames, 7 channels"
ALL = {"B": dict(col="prop_voxel_armB", name="Arm B: all four loops (CDR1, CDR2, HV4, CDR3), 7 channels", short="all loops", atoms="loops", sigma=2.0,
                 lab="VXH6_molecule_labels.csv.gz", tests="VXH8_cluster_tests_prop_voxel_armB.csv", cmp="VXH8_arm_comparison.csv", ovl="VXH8_overlap.csv",
                 thr="VXH6_thresholds.csv", chk="VXH6_checks.csv", summ="VXH6_summary.json", D=VXP("voxel_out/tmp/VXH5_D_receptors.f32.npy")),
       "C": dict(col="prop_voxel_armC", name="Arm C: CDR3 only, 7 channels", short="CDR3 only", atoms="cdr3", sigma=2.0,
                 lab="VXC2_molecule_labels.csv.gz", tests="VXH8_cluster_tests_prop_voxel_armC.csv", cmp="VXH8_arm_comparison.csv", ovl="VXH8_overlap.csv",
                 thr="VXC2_thresholds.csv", chk="VXC2_checks.csv", summ="VXC2_summary.json", D=VXP("voxel_out/tmp/VXC1_D_receptors.f32.npy")),
       "D": dict(col="prop_voxel_armD", name="Arm D: all four loops (CDR1, CDR2, HV4, CDR3), 7 channels, sigma 1.5 A", short="all loops, sigma 1.5", atoms="loops", sigma=1.5,
                 lab="VXS2_D_molecule_labels.csv.gz", tests="VXS5_cluster_tests_prop_voxel_armD.csv", cmp="VXS5_arm_comparison.csv", ovl="VXS5_overlap.csv",
                 thr="VXS2_D_thresholds.csv", chk="VXS2_D_checks.csv", summ="VXS2_D_summary.json", D=VXP("voxel_out/tmp/VXS1_D_D_receptors.f32.npy")),
       "E": dict(col="prop_voxel_armE", name="Arm E: CDR3 only, 7 channels, sigma 1.5 A", short="CDR3 only, sigma 1.5", atoms="cdr3", sigma=1.5,
                 lab="VXS2_E_molecule_labels.csv.gz", tests="VXS5_cluster_tests_prop_voxel_armE.csv", cmp="VXS5_arm_comparison.csv", ovl="VXS5_overlap.csv",
                 thr="VXS2_E_thresholds.csv", chk="VXS2_E_checks.csv", summ="VXS2_E_summary.json", D=VXP("voxel_out/tmp/VXS1_E_D_receptors.f32.npy"))}
ARM = {a: ALL[a] for a in (("B", "C") if REPORT == "main" else ("D", "E", "B", "C"))}
# page-facing arm table (D.arms); the main report keeps its original labels
JS = {"main": {"B": ("Arm B · all loops", "Arm B: " + LOOPS_NAME), "C": ("Arm C · CDR3 only", "Arm C: " + CDR3_NAME)},
      "s15": {"D": ("Arm D · all loops · σ 1.5", "Arm D: " + LOOPS_NAME + ", σ 1.5 Å"), "E": ("Arm E · CDR3 only · σ 1.5", "Arm E: " + CDR3_NAME + ", σ 1.5 Å"),
              "B": ("Arm B · all loops · σ 2.0", "Arm B: " + LOOPS_NAME + ", σ 2.0 Å"), "C": ("Arm C · CDR3 only · σ 2.0", "Arm C: " + CDR3_NAME + ", σ 2.0 Å")}}[REPORT]
arms_js = {a: dict(lab=a, short=JS[a][0], name=JS[a][1], atoms=P["atoms"], sigma=P["sigma"]) for a, P in ARM.items()}
mol = {k: v for k, v in A30["mol"].items() if k != "a30"}
for a, P in ARM.items():
    L_ = pd.read_csv(O + P["lab"])
    assert list(L_.clone_id) == A30["mol"]["id"], "molecule order differs from the reference report"
    mol[a] = L_[P["col"]].astype(int).tolist()
NMOL = len(mol["id"])
clusters, surv, summary, stab, thr = {}, {}, [], {}, {}
for a, P in ARM.items():
    S8 = pd.read_csv(O + P["cmp"]).set_index("arm"); OV = pd.read_csv(O + P["ovl"])
    K = pd.read_csv(O + P["tests"]).set_index("cluster")
    lab = np.array(mol[a]); size = pd.Series(lab[lab >= 0]).value_counts(); rows = []
    for c, n in size.items():
        k = K.loc[c] if c in K.index else None
        rows.append([int(c), int(n), int(k.members) if k is not None else 0, k.top_state if k is not None else "",
                     float(k.frac) if k is not None else None, int(k.mice) if k is not None else 0, int(k.v_pairs) if k is not None else 0,
                     float(k.p) if k is not None else None, float(k.q_bh) if k is not None else None])
    clusters[a] = sorted(rows, key=lambda r: (-r[1], r[0]))
    sv = K[K.q_bh <= 0.15].reset_index().sort_values("q_bh", kind="stable")
    surv[a] = [dict(cluster=int(r.cluster), members=int(r.members), top=r.top_state, frac=float(r.frac), q=float(r.q_bh), p=float(r.p),
                    vp=int(r.v_pairs), mice=int(r.mice), size=int((lab == r.cluster).sum())) for r in sv.itertuples()]
    s8 = S8.loc[P["col"]]; sm = json.load(open(O + P["summ"])); t = pd.read_csv(O + P["thr"]).iloc[0]
    ck = pd.read_csv(VXP("voxel_out/checks/") + P["chk"]).set_index("check")
    st = [float(ck.loc[k, "value"].split("/")[-1]) for k in ck.index if k.startswith("stability")]
    stab[a] = st
    thr[a] = t.to_dict()
    assert abs(float(t.bg_p1) - sm["cut"]) < 1e-12 and int(s8.survivors_q015) == len(surv[a])
    ov = OV[OV.arm == P["col"]].set_index("reference")
    summary.append(dict(arm=a, name=P["name"], short=P["short"], cut=float(t.bg_p1), sens=float(t.sens_at_bg_p1),
                        crystal_over_cut=sm["crystal_error_over_cut"], clusters=sm["clusters"], clustered=sm["clustered"],
                        clusters_ge3=int(s8.clusters), molecules=int(s8.molecules), purity=float(s8.observed),
                        **{k: float(s8[k]) for k in s8.index if k.startswith(("delta_", "z_", "null_")) and pd.notna(s8[k])},
                        survivors=int(s8.survivors_q015), clusters_tested=int(s8.clusters_tested),
                        overlap={r: dict(clusters=int(v.clusters), recovered=int(v.recovered)) for r, v in ov.iterrows()}))
S8 = pd.read_csv(O + "VXH8_arm_comparison.csv").set_index("arm")
ref30 = S8.loc[[k for k in S8.index if _REFCOL and k.startswith(_REFCOL)][0]]
# ---- exact split (stage 2)
split = {}
for a in ARM:
    z = np.load(W + "split_%s.npz" % a, allow_pickle=True); assert list(z["clone_id"]) == mol["id"]
    s = z["split"].astype("<f4"); assert (np.isnan(s).all(1) == (np.array(mol[a]) < 0)).all()
    split[a] = base64.b64encode(np.nan_to_num(s, nan=-1).astype("<f4").tobytes()).decode()
parts = [str(p) for p in np.load(W + "split_B.npz", allow_pickle=True)["parts"]]
# ---- worked example: the the reference method report's pair, in both arms, split exactly
EX = tuple(_VXCFG.get("report", {}).get("worked_example", []))   # two receptor ids to compare on the worked-example page
if len(EX) != 2:
    raise SystemExit('set report.worked_example (two receptor ids) in the config')
L1 = np.load(VXP("landmarks/out/LM1_internal_coords.npz"), allow_pickle=True); IDS = [str(c) for c in L1["clone_id"]]
BXA = {a: (boxes_for("F1") if P["atoms"] == "loops" else load_boxes()) for a, P in ARM.items()}
example = dict(a=EX[0], b=EX[1], ai=mol["id"].index(EX[0]), bi=mol["id"].index(EX[1]), arms={})
for a, P in ARM.items():
    atoms, bx, Df = P["atoms"], BXA[a], P["D"]
    G = [build_grid((VXP("structures/%s.pdb") % c, rep_place(IDS.index(c), "F1"), bx, atoms, 7, None, None, False, P["sigma"])) for c in EX]
    v = [g[0].astype(np.float16).astype(np.float64) for g in G]; out, s = [], 0
    for k in "AB":
        nv = int(np.prod(bx[k][1]))
        for c in range(7):
            d = v[0][s + c * nv: s + (c + 1) * nv] - v[1][s + c * nv: s + (c + 1) * nv]; out.append(float((d * d).sum() / H ** 3))
        s += 7 * nv
    Dm = np.load(Df, mmap_mode="r"); dst = float(Dm[IDS.index(EX[0]), IDS.index(EX[1])])
    assert abs(np.sqrt(sum(out)) / dst - 1) < 1e-5
    example["arms"][a] = dict(parts=out, D=float(np.sqrt(sum(out))), D_stored=dst, cut=float(thr[a]["bg_p1"]), atoms=[G[0][1], G[1][1]],
                              lab=[mol[a][example["ai"]], mol[a][example["bi"]]])
# ---- validation, hinge, first design, comparison, predictions
E = pd.read_csv(O + "VXV2_endpoints.csv"); V3 = pd.read_csv(O + "VXV3_cdr3.csv"); RU = json.load(open(O + "VXV2_rules.json"))
H2 = json.load(open(O + "VXH2_summary.json")); HS = pd.read_csv(O + "VXH2_subsets.csv"); H3 = json.load(open(O + "VXH3_summary.json"))
TI = pd.read_csv(O + "VXH2_trangle_icc.csv"); HG = pd.read_csv(O + "VX5b_hinge.csv")
SW = pd.read_csv(O + "VX2b_sweep.csv"); V5 = json.load(open(O + "VX5_summary.json")); B5 = pd.read_csv(O + "VX5b_configs.csv")
F1 = np.load(O + "VX1_frames.npz", allow_pickle=True); mm = F1["is_molecule"]
CMP = pd.read_csv(O + "VXH9_comparison.csv"); PRD = pd.read_csv(O + "VXH9_predictions.csv"); DG = pd.read_csv(O + "VXC3_diagnostics.csv")
X0 = pd.read_csv(O + "VXV0_receptors.csv")
# ---- the counts the page text states, from this run's own files (read back by rv5_gates.py, RG3)
_CP = {"B": "VXH6_control_pairs.csv", "C": "VXC2_control_pairs.csv", "D": "VXS2_D_control_pairs.csv", "E": "VXS2_E_control_pairs.csv"}
_XP = {"C": "VXC2_crystal_pairs.csv", "D": "VXS2_D_crystal_pairs.csv", "E": "VXS2_E_crystal_pairs.csv"}
counts = dict(n_models=len(IDS), n_anchors=len(L1["anch_keys"]), n_control_pairs=len(pd.read_csv(O + _CP["B"])),
              n_benchmark_pairs=len(pd.read_csv(O + "VX5_crystal_floor.csv")), n_panel_receptors=int(RU["n_receptors"]))
assert counts["n_models"] == len(F1["clone_id"]) and counts["n_anchors"] == L1["anch"].shape[1], "landmark file vs VX1 frames"
assert all(len(pd.read_csv(O + _CP[a])) == counts["n_control_pairs"] for a in ARM), "control-pair counts differ between arms"
assert all(len(pd.read_csv(O + _XP[a])) == counts["n_benchmark_pairs"] for a in ARM if a in _XP), "crystal-pair counts differ from VX5"
assert counts["n_panel_receptors"] == int(X0.cdr3_complete.sum()), "panel receptors differ from VXV0"
validation = dict(endpoints=E.to_dict("records"), cdr3=V3.to_dict("records"), rules=RU,
                  panel=dict(receptors=int(RU["n_receptors"]), pairs=int(RU["n_pairs"]), mouse=int((X0.species == "Mus musculus").sum()),
                             hybrid=int((X0.species != "Mus musculus").sum()), bound=int(X0.bound_pMHC.sum()), unbound=int((~X0.bound_pMHC).sum())),
                  hinge=dict(rel=H2["rel_h"], ci=H2["rel_h_ci"], relF1=H2["rel_hF1"], ciF1=H2["rel_hF1_ci"], H1=H2["H1_pass"],
                             e1=H2["vxv_E1"], e2=H2["vxv_E2"], subsets=HS.to_dict("records"), trangle=TI.to_dict("records"),
                             cm_median=float(HG.hinge_deg.median()), cm_p95=float(HG.hinge_deg.quantile(.95)),
                             rep_median=float(np.median(F1["vavb_rot_deg"][mm])), rep_p95=float(np.percentile(F1["vavb_rot_deg"][mm], 95)),
                             r2=H3["cv_r2"], r2v=H3["cv_r2_vpair"], r2l=H3["cv_r2_lengths"], rho_vc=H3["spearman_dhinge_vs_vcori"]),
                  sweep=SW[["sigma", "cut_p1", "ratio_shared", "ratio_full"]].to_dict("records"), vx5=V5,
                  frames=B5[~B5.report_only].to_dict("records"), bracket=B5[(B5.frame == "F0") & (B5.atoms == "whole")][["sigma", "crystal_ratio_median"]].to_dict("records"))
frames = dict(lm10=float(np.median(F1["rmsd_lm10"][mm])), anch=float(np.median(F1["rmsd_anch123"][mm])), anch95=float(np.percentile(F1["rmsd_anch123"][mm], 95)),
              fa=float(np.median(F1["fit_alpha5"][mm])), fb=float(np.median(F1["fit_beta5"][mm])), rot=float(np.median(F1["vavb_rot_deg"][mm])))
BX0 = json.load(open(O + "VXC0_boxes.json"))["boxes"]; BX5 = json.load(open(O + "VX5b_boxes.json"))["boxes"]
boxes = {a: {k: (BX5 if P["atoms"] == "loops" else BX0)["F1_%s" % k]["shape"] for k in "AB"} for a, P in ARM.items()}
DI = json.load(open(W + "display_index.json"))
# grids in the page (cluster 3D view): the arms' boxes and blur, and the per-atom channel masses of vxgrid.type_atoms for every
# (residue, atom) name pair in the display files, so the page types atoms with the Python table itself (no JS re-implementation)
import vxgrid as vg
WT = [[vg.type_atoms([("A", 1, "", r, a)])[0][0, :7].tolist() for a in DI["anames"]] for r in DI["rnames"]]
vox = dict(h=H, trunc=3.0, W=WT, cdr3=[105, 117], arms={a: dict(sigma=P["sigma"], atoms=P["atoms"]) for a, P in ARM.items()},
           boxes={a: {k: dict(lo=[float(x) for x in bx[k][0]], shape=[int(x) for x in bx[k][1]]) for k in "AB"} for a, bx in BXA.items()})
assert all(vox["boxes"][a][k]["shape"] == boxes[a][k] for a in ARM for k in "AB")
data = dict(states=A30["states"], mol=mol, clusters=clusters, cells=A30["cells"], tcl=A30["tcl"], genes=A30["genes"],
            atoms=dict(chunk=DI["chunk"], rnames=DI["rnames"], anames=DI["anames"], nchunks=DI["nchunks"]),
            aln=dict(cols=DI["cols"], origin=A30["aln"]["origin"], junction_notes=A30["aln"]["junction_notes"]),
            split=split, parts=parts, example=example, summary=summary, surv=surv, stab=stab, thr=thr, boxes=boxes,
            ref30=dict(z_mouse_V=float(ref30.z_mouse_V), delta_mouse_V=float(ref30.delta_mouse_V), survivors=int(ref30.survivors_q015),
                       z_mouse=float(ref30.z_mouse), clusters_ge3=int(ref30.clusters)),
            comparison=CMP.to_dict("records"), predictions=PRD.to_dict("records"), confounds=DG.to_dict("records"),
            validation=validation, frames=frames, n_state=A30["n_state"], vox=vox, arms=arms_js, report=REPORT,
            dataset_label=_RCFG.get("dataset_label", "data set"), brand_note=_RCFG.get("brand_note", "internal"),
            condition_labels=_RCFG.get("condition_labels", ["group A", "group B"]))
data["counts"] = counts   # models, anchors, control pairs, crystal-model pairs, panel receptors (page text)
# ---- A15 tcrdist3 cross-check (VXT1, VXT2): per arm and cluster, verdict and pairwise TCRdist; the reference method primary for reference
VC = {"also": 0, "partly": 1, "not": 2}
tcrd = dict(cut=json.load(open(O + "VXT1_summary.json"))["cut"], summary=json.load(open(O + "VXT1_summary.json")),
            ref=pd.read_csv(O + "VXT2_summary.csv").set_index("arm").loc[_VXCFG.get("reference_method", {}).get("name", "reference")].to_dict(),
            ref_name=_VXCFG.get("reference_method", {}).get("name", "reference method"),
            summ=pd.read_csv(O + "VXT2_summary.csv").set_index("arm").loc[list(ARM)].reset_index().to_dict("records"), arms={})
for a in ARM:
    X = pd.read_csv(O + "VXT2_crosscheck_%s.csv" % a); Z = np.load(O + "VXT2_pairs_%s.npz" % a)
    assert list(Z["cluster"]) == list(X.cluster) and set(X.cluster) == {r[0] for r in clusters[a] if r[1] >= 2}
    tcrd["arms"][a] = dict(rows=[[int(r.cluster), int(r.tcrdist3_main), float(r.tcrdist3_frac), VC[r.verdict], int(r.tcrdist3_clusters), int(r.tcrdist3_singletons),
                                  float(r.tcrdist_median), int(r.tcrdist_max), float(r.share_within_cut), float(r.median_bg_pct), str(r.tcrdist3_top)] for r in X.itertuples()],
                           off=[int(x) for x in Z["offset"]], d=base64.b64encode(Z["d"].astype("<i2").tobytes()).decode())
data["tcrd"] = tcrd; mol["T"] = pd.read_csv(O + "VXT1_molecule_labels.csv.gz").set_index("clone_id").reindex(mol["id"]).prop_tcrdist3.astype(int).tolist()
if REPORT == "s15":                                      # A14 results (VXS1-VXS5)
    ck = lambda f: pd.read_csv(VXP("voxel_out/checks/%s_checks.csv") % f)
    s15 = dict(panel=pd.read_csv(O + "VXS3_panel.csv").to_dict("records"), diff=pd.read_csv(O + "VXS3_sigma_diff.csv").to_dict("records"),
               vector_E1=json.load(open(O + "VXS3_summary.json"))["vector_E1"],
               clus={a: json.load(open(O + "VXS2_%s_summary.json" % a)) for a in "DE"},
               crystal={a: pd.read_csv(O + "VXS2_%s_crystal_pairs.csv" % a).to_dict("list") for a in "DE"},
               diag=pd.read_csv(O + "VXS4_diagnostics.csv").to_dict("records"), ari=json.load(open(O + "VXS4_summary.json"))["ari"],
               overlap=pd.read_csv(O + "VXS5_overlap.csv").to_dict("records"),
               analytic={a: ck("VXS1_%s" % a).set_index("check").loc["analytic grid-free distance, 12 pairs", "value"] for a in "DE"},
               analytic_ref={a: ck(f).set_index("check").loc["analytic grid-free distance, 12 pairs", "value"] for a, f in (("B", "VXH5"), ("C", "VXC1"))},
               checks={f: dict(n=len(ck(f)), failed=int((~ck(f).passed.astype(bool)).sum())) for f in
                       ("VXS1_D", "VXS1_E", "VXS2_D", "VXS2_E", "VXS3", "VXS4", "VXS5")})
    assert all(v["failed"] == 0 for v in s15["checks"].values())
    data["s15"] = s15
def clean(o):
    if isinstance(o, float) and o != o: return None
    if isinstance(o, dict): return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)): return [clean(v) for v in o]
    if isinstance(o, np.integer): return int(o)
    if isinstance(o, np.bool_): return bool(o)
    if isinstance(o, np.floating): return None if np.isnan(o) else float(o)
    return o
s = json.dumps(clean(data), separators=(",", ":")); open(W + ("voxel_html_data.json" if REPORT == "main" else "voxel_s15_html_data.json"), "w").write(s)
print("%s: json %.2f MB | " % (REPORT, len(s) / 1e6) + " | ".join("%s: clusters %d, survivors %d, example D %.3f" % (a, len(clusters[a]), len(surv[a]), example["arms"][a]["D"]) for a in ARM))
