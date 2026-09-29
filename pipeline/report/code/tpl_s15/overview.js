
// ---------- OVERVIEW (s15 report: sigma 1.5 A arms D and E beside their sigma 2.0 counterparts B and C; A14, report only)
const PAIRS = [["D", "B", "all four loops"], ["E", "C", "CDR3 only"]];
// the reference method's survivor set, as named in the state-test overlap table (see vxh8/vxs5 `sets`)
const REF_SET_KEY = Object.keys((D.summary[0] || {}).overlap || {}).find(k => /survivors/.test(k) && !/own/.test(k)) || "";
const REF_SET_LABEL = REF_SET_KEY ? "Reference method: " + REF_SET_KEY : "reference survivors";
function survTable(a, rows) {
  const cols = [
    {h: "Cluster", v: r => r.cluster, f: r => `<a href="#c/${a}/${r.cluster}">${r.cluster}</a>`, num: 1},
    {h: "Molecules", v: r => r.size, num: 1},
    {h: "Top state", v: r => r.top, f: r => `<span class="dot" style="background:var(${stateVar(r.top)})"></span>${r.top} ${pct(r.frac)}`},
    {h: "q", v: r => r.q, f: r => fq(r.q), num: 1}, {h: "Mice", v: r => r.mice, num: 1}, {h: "V pairs", v: r => r.vp, num: 1}];
  return table(cols, rows, {scroll: false});
}
const P15 = (atoms, nch, sig) => D.s15.panel.find(r => r.atoms === atoms && r.channels === nch && r.sigma === sig);
const DIF = (atoms, nch) => D.s15.diff.find(r => r.atoms === atoms && r.channels === nch);
const DIAG = a => D.s15.diag.find(r => r.arm.startsWith("Arm " + a));
const OVL = (a, ref) => D.s15.overlap.find(r => r.arm === "prop_voxel_arm" + a && r.reference === ref);
const ciTxt = (x, lo, hi, d = 3) => `${x >= 0 ? "+" : "−"}${f(Math.abs(x), d)} <span class="muted">[${f(lo, d)}, ${f(hi, d)}]</span>`;
const verdict = (lo, hi, up, down) => lo > 0 ? up : hi < 0 ? down : "no clear change (the interval includes zero)";
function viewOverview(app) {
  sec(app, `<h1>A sharper voxel grid: σ 1.5 Å</h1>
  <p class="muted">${esc(D.dataset_label)} · ${NMOL.toLocaleString()} distinct receptor molecules (TCRBuilder2+ models) · \u03c3 sensitivity arms · ${esc(D.brand_note)}.</p>
  <div class="warnbox"><b>What this report adds.</b> The first report's grids blur each atom with σ 2.0 Å, about twice the width of the usual atom-sized densities, so fine side-chain differences are smoothed out. Here the same two grids are rebuilt with <b>σ 1.5 Å</b> — the sharpest blur the 1 Å voxel grid samples well — and run through exactly the same clustering and state tests: <b>Arm D</b> (all four loops) and <b>Arm E</b> (CDR3 only), shown beside <b>Arm B</b> and <b>Arm C</b> (the same atoms at σ 2.0). σ 1.5 is an intermediate step, not van der Waals-sized (that would be about 0.85 Å on a 0.5 Å grid). D and E are <b>sensitivity arms, report only</b>: they replace nothing, and a higher crystal–model error was expected and accepted in advance.</div>`);
  const rows = [
    ["Blur σ (Å)", a => f(ARMS[a].sigma, 1)],
    ["Cut (1 % of random pairs)", a => f(SUM(a).cut, 4)],
    ["Control pairs within the cut", a => pct1(SUM(a).sens)],
    ["Crystal–model error / cut (" + D.counts.n_benchmark_pairs.toLocaleString("en-US") + " pairs)", a => f(SUM(a).crystal_over_cut, 3)],
    ["E1: crystal structure recovered beyond sequence", a => { const p = P15(ARMS[a].atoms, 7, ARMS[a].sigma); return `${f(p.E1_partial, 3)} <span class="muted">[${f(p.E1_lo, 2)}, ${f(p.E1_hi, 2)}]</span>`; }],
    ["E4: reliability, crystal vs own model", a => { const p = P15(ARMS[a].atoms, 7, ARMS[a].sigma); return `${f(p.E4_rel, 3)} <span class="muted">[${f(p.E4_lo, 2)}, ${f(p.E4_hi, 2)}]</span>`; }],
    ["Clusters ≥ 2 / ≥ 3 (≥ 3 counted on molecules with a state)", a => `${SUM(a).clusters.toLocaleString()} / ${SUM(a).clusters_ge3.toLocaleString()}`],
    ["Molecules clustered", a => SUM(a).clustered.toLocaleString()],
    ["Single-V-pair clusters", a => pct(DIAG(a).share_single_Vpair_clusters)],
    ["Distance vs CDR3 length difference (Spearman)", a => f(DIAG(a).spearman_dlen_total, 3)],
    ["State excess within mouse (z)", a => `${sgn(SUM(a).delta_mouse)} (${f(SUM(a).z_mouse, 1)})`],
    ["… within mouse + V pair (z)", a => `${sgn(SUM(a).delta_mouse_V)} (${f(SUM(a).z_mouse_V, 1)})`],
    ["… + CDR3 length class (z)", a => `${sgn(SUM(a).delta_mouse_V_len)} (${f(SUM(a).z_mouse_V_len, 1)})`],
    ["Significant clusters (q ≤ 0.15)", a => `${SUM(a).survivors} <span class="muted">of ${SUM(a).clusters_tested.toLocaleString()}</span>`],
    [REF_SET_LABEL + " recovered", a => { const o = SUM(a).overlap[REF_SET_KEY]; return o ? o.recovered : "\u2013"; }]];
  const order = ["D", "B", "E", "C"];
  card(app, `<h2 style="margin-top:0">σ 1.5 beside σ 2.0</h2>
   <p class="small muted">Same pairs, seeds, cut rule, linkage and state tests for all four arms. E1 and E4 from the validity panel (${D.counts.n_panel_receptors.toLocaleString("en-US")} crystal-benchmarked receptors, 7-channel grids, 95 % bootstrap intervals). Reference method for comparison: state excess within mouse + V pair z ${f(D.ref30.z_mouse_V, 1)}, ${D.ref30.survivors} significant clusters.</p>
   <div class="tw"><table><thead><tr><th></th>${order.map(a => `<th>${ARMS[a].short}</th>`).join("")}</tr></thead>
   <tbody>${rows.map(([n, g]) => `<tr><td>${n}</td>${order.map(a => `<td class="num">${g(a)}</td>`).join("")}</tr>`).join("")}</tbody></table></div>`);
  const li = PAIRS.map(([a, b, what]) => {
    const d7 = DIF(ARMS[a].atoms, 7), c = D.s15.clus[a], A = SUM(a), B = SUM(b);
    return `<li><b>${what} (${a} vs ${b}).</b> Structure recovered beyond sequence, E1: ${ciTxt(d7.dE1, d7.dE1_lo, d7.dE1_hi)} → ${verdict(d7.dE1_lo, d7.dE1_hi, "the sharper grid recovers more crystal-verified structure", "the sharper grid recovers less")}. Reliability E4: ${ciTxt(d7.dE4, d7.dE4_lo, d7.dE4_hi)}. Crystal–model error / cut ${f(A.crystal_over_cut, 3)} against ${f(B.crystal_over_cut, 3)}; control pairs within the cut ${pct1(A.sens)} against ${pct1(B.sens)}. The partitions agree with ARI ${f(c.ari_vs_sigma20_arm, 2)}. State excess within mouse + V pair z ${f(A.z_mouse_V, 1)} against ${f(B.z_mouse_V, 1)}; ${A.survivors} significant clusters against ${B.survivors}, of which ${OVL(a, "prop_voxel_arm" + b + " survivors recovered here").recovered} of ${b}'s are recovered here.</li>`;
  }).join("");
  card(app, `<h2 style="margin-top:0">What changes at σ 1.5 (read from the numbers above)</h2><ul>${li}
   <li>The paired differences use the same bootstrap resamples of the ${D.counts.n_panel_receptors.toLocaleString("en-US")} receptors for both blurs, so they are more precise than comparing the two intervals by eye. Details on <a href="#validation">Validation</a>.</li>
   <li>The state tests of D and E are an additional look at the same data: q values are corrected within each arm only, and no arm is chosen on them.</li></ul>`);
  tcrdOverview(app);
  for (const a of ["D", "E"]) {
    sec(app, `<h2>Significant clusters · ${ARMS[a].short}</h2><p class="small muted">Sorted by q. Click a cluster to open its page (voxel grids drawn at the arm's own σ).</p>`);
    app.append(survTable(a, D.surv[a]));
  }
}
