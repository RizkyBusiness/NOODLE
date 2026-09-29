
// ---------- OVERVIEW (both arms side by side; neither is preferred)
function survTable(a, rows) {
  const cols = [
    {h: "Cluster", v: r => r.cluster, f: r => `<a href="#c/${a}/${r.cluster}">${r.cluster}</a>`, num: 1},
    {h: "Molecules", v: r => r.size, num: 1},
    {h: "Top state", v: r => r.top, f: r => `<span class="dot" style="background:var(${stateVar(r.top)})"></span>${r.top} ${pct(r.frac)}`},
    {h: "q", v: r => r.q, f: r => fq(r.q), num: 1}, {h: "Mice", v: r => r.mice, num: 1}, {h: "V pairs", v: r => r.vp, num: 1}];
  return table(cols, rows, {scroll: false});
}
function armTiles(a) {
  const P = SUM(a);
  return `<div class="card"><h2 style="margin-top:0">${ARMS[a].short}</h2><p class="small muted">${ARMS[a].name}</p><div class="grid g2">
   <div class="card tile"><div class="lab">Molecules clustered</div><div class="val">${P.clustered.toLocaleString()}</div><div class="sub">in ${P.clusters.toLocaleString()} clusters of ≥ 2</div></div>
   <div class="card tile"><div class="lab">Significant clusters (q ≤ 0.15)</div><div class="val">${P.survivors}</div><div class="sub">of ${P.clusters_tested.toLocaleString()} tested</div></div>
   <div class="card tile"><div class="lab">Excess, same mouse + V pair</div><div class="val">${sgn(P.delta_mouse_V)}</div><div class="sub">z ${f(P.z_mouse_V, 1)} · + CDR3 length class: z ${f(P.z_mouse_V_len, 1)}</div></div>
   <div class="card tile"><div class="lab">Control pairs kept together</div><div class="val">${pct1(P.sens)}</div><div class="sub">${D.counts.n_control_pairs.toLocaleString("en-US")} pairs with one CDR3 amino-acid difference</div></div>
   <div class="card tile"><div class="lab">Crystal–model error / cut</div><div class="val">${f(P.crystal_over_cut, 2)}</div><div class="sub">report only</div></div>
   <div class="card tile"><div class="lab">Single-V-pair clusters</div><div class="val">${pct(D.confounds.find(r => r.arm.startsWith("Arm " + a)).share_single_Vpair_clusters)}</div><div class="sub">share of clusters using a single V-gene pair</div></div></div></div>`;
}
function viewOverview(app) {
  const B = SUM("B"), C = SUM("C"), R = D.ref30;
  sec(app, `<h1>Voxel-grid structure and T-cell state</h1>
  <p class="muted">${esc(D.dataset_label)} · ${NMOL.toLocaleString()} distinct receptor molecules (TCRBuilder2+ models) · voxel-grid descriptor, independent of the reference vector method · ${esc(D.brand_note)}.</p>
  <div class="warnbox"><b>Method: a 3D voxel grid of atoms and chemistry.</b> Each receptor's loop atoms are smeared into a grid of 1 Å cubes (σ 2.0 Å) carrying occupancy and six chemistry channels, with each chain aligned on its own five framework landmarks. Two arms, both reported, neither preferred: <b>Arm B</b> uses all four loops, <b>Arm C</b> the two CDR3 loops only. Details on <a href="#methods">Methods</a>; the benchmark tests on <a href="#validation">Validation</a>.</div>`);
  sec(app, `<div class="grid g2">${armTiles("B")}${armTiles("C")}</div>`);
  const cmp = D.comparison;
  card(app, `<h2 style="margin-top:0">Against the reference method</h2>
   <p class="small muted">Same pairs, same seeds, same state-test machinery (reproduced exactly on the the reference method before the arms were tested). E1 = crystal-verified structure recovered from models beyond sequence, on ${D.counts.n_panel_receptors.toLocaleString("en-US")} crystal-benchmarked receptors (higher is better).</p>`).append(table([
    {h: "Method", v: r => r.arm}, {h: "Scope", v: r => r.scope},
    {h: "Controls in cut", v: r => r.control_within_cut, f: r => pct1(r.control_within_cut), num: 1},
    {h: "Crystal err / cut", v: r => r.crystal_err_over_cut, f: r => f(r.crystal_err_over_cut, 3), num: 1},
    {h: "Single-V clusters", v: r => r.single_Vpair_clusters, f: r => r.single_Vpair_clusters == null ? "–" : pct(r.single_Vpair_clusters), num: 1},
    {h: "E1 (validity)", v: r => r.VXV_E1, f: r => esc(r.VXV_E1 || "–")},
    {h: "Excess, mouse + V (z)", v: r => r.z_mouse_V, f: r => r.z_mouse_V == null ? esc(r.state) : `${sgn(r.excess_mouse_V)} (${f(r.z_mouse_V, 1)})`, num: 1},
    {h: "+ length class (z)", v: r => r.z_mouse_V_len, f: r => r.z_mouse_V_len == null ? "–" : `${sgn(r.excess_mouse_V_len)} (${f(r.z_mouse_V_len, 1)})`, num: 1},
    {h: "Survivors", v: r => r.survivors_q015, f: r => r.survivors_q015 == null ? "–" : r.survivors_q015, num: 1}], cmp, {scroll: false}));
  card(app, `<h2 style="margin-top:0">How to read this</h2><ul>
  <li><b>Does the grid capture real structure?</b> The crystal benchmark on the <a href="#validation">Validation</a> page (E1 beyond sequence, E2 beyond the reference method, E4 reliability, crystal\u2013model error over the cut).</li>
  <li><b>Is it noisier than the reference method?</b> Compare the control pairs kept together and the crystal\u2013model error in the tiles above.</li>
  <li><b>Do the clusters track cell state, and through what?</b> The null ladder: free, within mouse, within mouse + V pair, + CDR3 length class (Methods section 7 and the table above).</li>
  <li><b>Are the clusters V-gene groups?</b> The single-V-pair share per arm; the CDR3-only arm is the control for that.</li>
  <li><b>Are they sequence groups?</b> The tcrdist3 cross-check below and on every cluster page.</li>
  <li><b>How stable are they?</b> The ARI at cut \u00d7 0.98 / 1.02 reported in Methods.</li></ul>`);
  tcrdOverview(app);
  for (const a of ["B", "C"]) {
    sec(app, `<h2>Significant clusters · ${ARMS[a].short}</h2><p class="small muted">Sorted by q. Click a cluster to open its page.</p>`);
    app.append(survTable(a, D.surv[a]));
  }
}
