
// ---------- VALIDATION (validity panel, CDR3-only check, hinge [exploratory], first design) and its plain version
const NAMES = {grid_F1L_occ: "grid · loops · occupancy (primary)", grid_F1L_7ch: "grid · loops · 7 channels (= Arm B)", grid_F0W_occ: "grid · whole molecule, one frame · occupancy",
  grid_F0L_occ: "grid · loops, one frame · occupancy", grid_F1W_occ: "grid · whole molecule, per chain · occupancy", grid_F1C3_occ: "grid · CDR3 · occupancy",
  grid_F1C3_7ch: "grid · CDR3 · 7 channels (= Arm C)", vec_vcori_geo: "reference descriptor · geometry (primary)", vec_vcori_full: "reference descriptor · with chemistry",
  vec_vc: "vector vc · geometry", vec_vg: "vector vg · geometry"};
const ci3 = (x, lo, hi) => `${f(x, 3)} <span class="muted">[${f(lo, 3)}, ${f(hi, 3)}]</span>`;
function viewValidation(app) {
  const V = D.validation, P = V.panel, R = V.rules, Hn = V.hinge;
  sec(app, `<h1>Validation</h1><p class="muted">How much real structure each descriptor recovers from predicted models, checked against lab-solved crystal structures. No transcriptional data are used on this page.</p>
  <h2>1 · Validity panel</h2>
  <p>${P.receptors} receptors with a crystal structure (${P.mouse} mouse, ${P.hybrid} hybrid human + mouse; ${P.bound} bound to pMHC, ${P.unbound} unbound), ${P.pairs.toLocaleString()} receptor pairs. For each descriptor, the crystal–crystal distances (C) and the model–model distances (M) between the same receptors. <b>E1</b> = partial Spearman correlation of C and M after removing what sequence similarity alone explains (BLOSUM62, Henikoff & Henikoff 1992): how much crystal-verified structure the descriptor recovers from models beyond sequence. <b>E4</b> = reliability, 1 − 2 × mean d²(crystal, own model) / mean d²(random model pairs). 95 % intervals by bootstrap over receptors (2,000 resamples); Mantel p (Mantel 1967) for the plain correlations.</p>`);
  const rows = V.endpoints.concat(V.cdr3.map(r => ({descriptor: r.descriptor, E1_partial: r.E1_partial, E1_lo: r.E1_lo, E1_hi: r.E1_hi, E1_plain: r.E1_plain, E4_rel: r.E4_rel, error_over_cut: r.error_over_cut, within_cut: r.within_cut, extra: true})));
  app.append(table([{h: "Descriptor", v: r => NAMES[r.descriptor] || r.descriptor, f: r => (r.descriptor === "grid_F1L_occ" || r.descriptor === "vec_vcori_geo" ? "<b>" : "") + esc(NAMES[r.descriptor] || r.descriptor) + (r.extra ? ' <span class="muted small">(CDR3-only check)</span>' : "")},
    {h: "E1 [95 % CI]", v: r => r.E1_partial, f: r => ci3(r.E1_partial, r.E1_lo, r.E1_hi), num: 1}, {h: "E1 plain", v: r => r.E1_plain, f: r => f(r.E1_plain, 3), num: 1},
    {h: "E4", v: r => r.E4_rel, f: r => f(r.E4_rel, 3), num: 1}, {h: "Crystal err / pilot cut", v: r => r.error_over_cut, f: r => f(r.error_over_cut, 3), num: 1},
    {h: "Within cut", v: r => r.within_cut, f: r => `${r.within_cut} / ${P.receptors}`, num: 1}], rows, {scroll: false}));
  card(app, `<h3 style="margin-top:0">Pre-registered rules (grid primary vs vector primary)</h3><ul>
   <li><b>Delta</b> = E1(grid) − E1(vector) = ${ci3(R.delta.point, R.delta.ci[0], R.delta.ci[1])}.</li>
   <li><b>E2</b>, structure the grid recovers beyond the vector method and sequence = ${ci3(R.E2_grid.point, R.E2_grid.ci[0], R.E2_grid.ci[1])}; the mirror (vector beyond grid) = ${ci3(R.E2_vector.point, R.E2_vector.ci[0], R.E2_vector.ci[1])}.</li>
   <li><b>Rule V2: grid ${R.V2_reinstated ? "reinstated" : "not reinstated"}</b> (reinstated if the lower bound of Delta > −0.05 or the lower bound of E2 > 0). Rule V3: ${Object.entries(R.V3_flag).map(([k, v]) => `${NAMES[k]} ${v ? "flagged" : "not flagged"}`).join("; ")}.</li>
   <li><b>CDR3 only</b> (report only): Delta ${V.cdr3.map(r => `${ci3(r.delta_vs_vector, r.delta_lo, r.delta_hi)} (${NAMES[r.descriptor]})`).join(", ")}.</li></ul>`);
  sec(app, `<h2>2 · The Vα/Vβ hinge — exploratory only</h2>
  <div class="warnbox"><b>Exploratory, labelled descriptor.</b> The hinge failed its pre-registered reliability test and enters no distance, cluster, score or state test.</div>
  <p>Relative orientation of the two domains, built from the reference framework placed on each chain. Reliability ${ci3(Hn.rel, Hn.ci[0], Hn.ci[1])} against the pre-registered bar of 0.50 on the lower bound → <b>${Hn.H1 ? "passed" : "failed"}</b> (five-landmark version ${ci3(Hn.relF1, Hn.ciF1[0], Hn.ciF1[1])}). Crystal vs its own model: median ${f(Hn.cm_median, 1)}° (95th pct ${f(Hn.cm_p95, 1)}°), against a repertoire spread of median ${f(Hn.rep_median, 1)}° (95th pct ${f(Hn.rep_p95, 1)}°) from the reference. Explained by V pair + CDR3 lengths: cross-validated R² ${f(Hn.r2, 2)} (V pair alone ${f(Hn.r2v, 2)}).</p>`);
  app.append(table([{h: "Subset", v: r => r.subset + " = " + r.value}, {h: "n", v: r => r.n, num: 1}, {h: "Reliability [95 % CI]", v: r => r.rel, f: r => ci3(r.rel, r.ci_lo, r.ci_hi), num: 1}], Hn.subsets, {scroll: false}));
  sec(app, `<p class="small muted">TRangle parameters (Dunbar et al. 2014, via STCRpy; ICC(2,1), Koo & Li 2016): ${Hn.trangle.map(t => `${t.param} ${f(t.icc21, 2)}`).join(" · ")} (n = ${Hn.trangle[0].n}). The TCRBuilder2+ training cutoff could not be established; newer crystals (2020 onward) are the least likely to be in its training data.</p>
  <h2>3 · The first design: whole molecule in one frame</h2>
  <p>Blur sweep on a 2,000-molecule pilot (crystal–model error over the 1 % cut, shared atoms; pre-registered rule → σ 2.0 Å):</p>`);
  app.append(table([{h: "σ (Å)", v: r => r.sigma, f: r => f(r.sigma, 1), num: 1}, {h: "Pilot cut", v: r => r.cut_p1, f: r => f(r.cut_p1, 3), num: 1},
    {h: "Error / cut (shared atoms)", v: r => r.ratio_shared, f: r => f(r.ratio_shared, 3), num: 1}, {h: "Error / cut (all atoms)", v: r => r.ratio_full, f: r => f(r.ratio_full, 3), num: 1}], V.sweep, {scroll: false}));
  sec(app, `<p>On the full repertoire the error over cut was <b>${f(V.vx5.ratio_shared, 4)}</b> against the pre-registered limit of 0.60 — the first design failed. A frame diagnostic then showed the frame was not the main cause:</p>`);
  app.append(table([{h: "Frame", v: r => r.frame}, {h: "Atoms", v: r => r.atoms}, {h: "Error / pilot cut", v: r => r.crystal_ratio_median, f: r => f(r.crystal_ratio_median, 3), num: 1},
    {h: "Controls in cut", v: r => r.control_within_cut, f: r => pct1(r.control_within_cut), num: 1}], V.frames, {scroll: false}));
  sec(app, `<p class="small muted">F0 whole molecule, ten landmarks; F1 per chain, five landmarks each (the frame of Arms B and C); F2 per chain, framework anchors; F3 whole molecule, ${D.counts.n_anchors.toLocaleString("en-US")} anchors. Report-only bracketing of the blur: ${V.bracket.map(r => `σ ${f(r.sigma, 1)} → ${f(r.crystal_ratio_median, 3)}`).join(" · ")}.</p>`);
}
function viewValidationExplained(app) {
  const V = D.validation, R = V.rules;
  sec(app, `<h1>Validation explained</h1><p class="muted">The same checks as the Validation page, in plain language.</p>
  <div class="card"><h2 style="margin-top:0">The question</h2><p>The structures in this project are computer predictions. For ${V.panel.receptors} receptors a laboratory-solved crystal structure also exists. So we can ask of any description of receptor shape: when it compares two receptors using their predictions, does it reach the same answer as when it compares their real crystal structures — beyond what their sequences alone would tell us?</p></div>
  <div class="card"><h2 style="margin-top:0">Result</h2><ul class="plain">
   <li>The voxel grid agrees with the crystals more than the existing method does (a score of ${f(V.endpoints.find(r => r.descriptor === "grid_F1L_occ").E1_partial, 2)} against ${f(V.endpoints.find(r => r.descriptor === "vec_vcori_geo").E1_partial, 2)}), and part of what it sees the existing method does not see at all. So the grid was kept.</li>
   <li>Looking only at the CDR3 loops, the grid still does better (by ${f(V.cdr3[0].delta_vs_vector, 2)}–${f(V.cdr3[1].delta_vs_vector, 2)}, confidently above zero).</li>
   <li>But the grid is also noisier: the prediction of the same receptor differs more from its crystal, compared with how different receptors are from each other.</li></ul></div>
  <div class="card"><h2 style="margin-top:0">The hinge</h2><p>We also measured the angle between the two halves of each receptor. Its predictions were not reliable enough (reliability ${f(V.hinge.rel, 2)}, with a lower bound of ${f(V.hinge.ci[0], 2)} against the required 0.50), especially for newer crystal structures that the prediction program probably never saw. It is therefore shown only for interest and is not used anywhere else.</p></div>
  <div class="card"><h2 style="margin-top:0">The first design</h2><p>The first version pictured the whole receptor in one alignment. Its prediction error was just over the pre-agreed limit (${f(V.vx5.ratio_shared, 3)} against 0.60), and changing the alignment did not fix it. That led to the loop-only and CDR3-only versions used in this report.</p></div>`);
}
