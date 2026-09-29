
// ---------- VALIDATION (s15 report: validity panel at sigma 1.5 beside 2.0; crystal-model pairs) and its plain version
const DNAME = r => `${r.atoms === "loops" ? "loops" : "CDR3"} · ${r.channels === 1 ? "occupancy" : "7 channels"} · σ ${f(r.sigma, 1)}` +
  (r.channels === 7 ? ` <span class="muted small">(= Arm ${r.atoms === "loops" ? (r.sigma === 1.5 ? "D" : "B") : (r.sigma === 1.5 ? "E" : "C")})</span>` : "");
const ci3 = (x, lo, hi) => `${f(x, 3)} <span class="muted">[${f(lo, 3)}, ${f(hi, 3)}]</span>`;
function viewValidation(app) {
  const S = D.s15, P = D.validation.panel;
  sec(app, `<h1>Validation — σ 1.5 beside σ 2.0</h1><p class="muted">How much real structure each grid recovers from predicted models, checked against lab-solved crystal structures. No transcriptional data are used on this page.</p>
  <h2>1 · Validity panel</h2>
  <p>The first report's panel (${P.receptors} crystal-benchmarked receptors, ${P.pairs.toLocaleString()} pairs; shared atoms; the same pair masks), rebuilt for the loop and CDR3 grids at both blurs. <b>E1</b> = partial Spearman correlation of crystal–crystal and model–model distances after removing sequence similarity (BLOSUM62): crystal-verified structure recovered beyond sequence. <b>Delta</b> = E1 − E1 of the vector method (${f(S.vector_E1, 3)}). <b>E2</b> = beyond the vector method and sequence. <b>E4</b> = reliability, 1 − 2 × mean d²(crystal, own model) / mean d²(random pairs). 95 % intervals: 2,000 bootstrap resamples of receptors. The σ 2.0 rows reproduce the published panel exactly (checked to 1e-9).</p>`);
  app.append(table([{h: "Grid", v: r => r.descriptor, f: r => DNAME(r)},
    {h: "E1 [95 % CI]", v: r => r.E1_partial, f: r => ci3(r.E1_partial, r.E1_lo, r.E1_hi), num: 1},
    {h: "Delta vs vector [CI]", v: r => r.delta_vs_vector, f: r => ci3(r.delta_vs_vector, r.delta_lo, r.delta_hi), num: 1},
    {h: "E2", v: r => r.E2_grid, f: r => f(r.E2_grid, 3), num: 1},
    {h: "E4 [CI]", v: r => r.E4_rel, f: r => ci3(r.E4_rel, r.E4_lo, r.E4_hi), num: 1},
    {h: "Error / pilot cut", v: r => r.error_over_cut, f: r => f(r.error_over_cut, 3), num: 1},
    {h: "Within cut", v: r => r.within_cut, f: r => `${r.within_cut} / ${P.receptors}`, num: 1}], S.panel, {scroll: false}));
  sec(app, `<h2>2 · Paired difference, σ 1.5 − σ 2.0</h2><p>Same receptors, same bootstrap resamples for both blurs. A positive dE1 whose interval excludes zero would mean the sharper grid recovers more crystal-verified structure.</p>`);
  app.append(table([{h: "Grid", v: r => r.atoms + r.channels, f: r => `${r.atoms === "loops" ? "loops" : "CDR3"} · ${r.channels === 1 ? "occupancy" : "7 channels"}`},
    {h: "dE1 [95 % CI]", v: r => r.dE1, f: r => `${sgn(r.dE1)} <span class="muted">[${f(r.dE1_lo, 3)}, ${f(r.dE1_hi, 3)}]</span>`, num: 1},
    {h: "dE4 [95 % CI]", v: r => r.dE4, f: r => `${sgn(r.dE4)} <span class="muted">[${f(r.dE4_lo, 3)}, ${f(r.dE4_hi, 3)}]</span>`, num: 1},
    {h: "Error / cut, σ 1.5", v: r => r.error_over_cut_s15, f: r => f(r.error_over_cut_s15, 3), num: 1},
    {h: "Error / cut, σ 2.0", v: r => r.error_over_cut_s20, f: r => f(r.error_over_cut_s20, 3), num: 1}], S.diff, {scroll: false}));
  sec(app, `<h2>3 · Crystal vs own model on the full repertoire's cut (${D.counts.n_benchmark_pairs.toLocaleString("en-US")} pairs)</h2>
  <p>The ${D.counts.n_benchmark_pairs.toLocaleString("en-US")} crystal–model pairs of the first design, rebuilt in each arm's configuration (shared atoms). The distance of a crystal to its own model, divided by the arm's cut (1 % of random pairs): lower is better, below 1 means the model falls within the cut of its crystal. <b>Reported, not a gate</b> (A14.3).</p>`);
  const cr = [["D", "B"], ["E", "C"]].map(([a, b]) => { const X = S.crystal[a], cA = SUM(a).cut, cB = SUM(b).cut;
    const r15 = X.d.map(v => v / cA), r20 = X.d_sigma20.map(v => v / cB), med = v => { const s = [...v].sort((p, q) => p - q); return s.length % 2 ? s[(s.length - 1) / 2] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2; };
    return {pair: `${ARMS[a].short} vs ${ARMS[b].short}`, m15: med(r15), m20: med(r20), w15: r15.filter(v => v <= 1).length, w20: r20.filter(v => v <= 1).length,
            worse: r15.filter((v, k) => v > r20[k]).length, n: r15.length}; });
  app.append(table([{h: "Arms", v: r => r.pair}, {h: "Median error / cut, σ 1.5", v: r => r.m15, f: r => f(r.m15, 3), num: 1},
    {h: "Median error / cut, σ 2.0", v: r => r.m20, f: r => f(r.m20, 3), num: 1}, {h: "Within cut, σ 1.5", v: r => r.w15, f: r => `${r.w15} / ${r.n}`, num: 1},
    {h: "Within cut, σ 2.0", v: r => r.w20, f: r => `${r.w20} / ${r.n}`, num: 1}, {h: "Pairs worse at σ 1.5", v: r => r.worse, f: r => `${r.worse} / ${r.n}`, num: 1}], cr, {scroll: false}));
  sec(app, `<p class="small muted">σ 2.0 column: the same ${D.counts.n_benchmark_pairs.toLocaleString("en-US")} pairs rebuilt with the same recipe at σ 2.0, divided by Arm B's or C's cut; the recipe reproduces the first design's pair distances (checked). Control pairs within the cut (${D.counts.n_control_pairs.toLocaleString("en-US")} pairs differing by one CDR3 amino acid): Arm D ${pct1(SUM("D").sens)} (B ${pct1(SUM("B").sens)}), Arm E ${pct1(SUM("E").sens)} (C ${pct1(SUM("C").sens)}).</p>`);
}
function viewValidationExplained(app) {
  const d7L = D.s15.diff.find(r => r.atoms === "loops" && r.channels === 7), d7C = D.s15.diff.find(r => r.atoms === "cdr3" && r.channels === 7);
  const say = r => r.dE1_lo > 0 ? "recovers more real structure than the σ 2.0 grid" : r.dE1_hi < 0 ? "recovers less real structure than the σ 2.0 grid" : "recovers about as much real structure as the σ 2.0 grid — no clear difference";
  sec(app, `<h1>Validation explained</h1><p class="muted">The same checks as the Validation page, in plain language.</p>
  <div class="card"><h2 style="margin-top:0">The question</h2><p>For ${D.validation.panel.receptors} receptors a laboratory-solved crystal structure exists beside the computer prediction. If the sharper grid is picking up real detail, comparing receptors through their predictions should agree better with comparing their crystals. If it is mostly picking up prediction errors, agreement should get worse.</p></div>
  <div class="card"><h2 style="margin-top:0">Result</h2><ul class="plain">
   <li>All four loops: the sharper grid ${say(d7L)} (change ${sgn(d7L.dE1)}, interval ${f(d7L.dE1_lo, 2)} to ${f(d7L.dE1_hi, 2)}).</li>
   <li>CDR3 only: the sharper grid ${say(d7C)} (change ${sgn(d7C.dE1)}, interval ${f(d7C.dE1_lo, 2)} to ${f(d7C.dE1_hi, 2)}).</li>
   <li>Noise: the prediction of the same receptor, compared with its crystal, sits at ${f(SUM("D").crystal_over_cut, 2)} (all loops) and ${f(SUM("E").crystal_over_cut, 2)} (CDR3) of the grouping threshold, against ${f(SUM("B").crystal_over_cut, 2)} and ${f(SUM("C").crystal_over_cut, 2)} at σ 2.0 (lower is better).</li></ul></div>`);
}
