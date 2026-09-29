
// ---------- CLUSTERS LIST
const LS = {minSize: 3, maxQ: 1, state: "", text: ""};
function viewClusters(app) {
  sec(app, `<h1>Clusters</h1><p class="muted">Complete linkage at the 1 % cut, on distinct molecules. Choose the arm; click a row to open a cluster. Search matches clone id, V/J gene or a CDR3 substring of any member.</p>`);
  const fl = sec(app, `<div class="filters" id="clf">
    <label>Min size <input id="f_min" type="number" min="2" value="${LS.minSize}" style="width:70px"></label>
    <label>Max q <select id="f_q"><option value="1">any</option><option value="0.15">≤ 0.15</option><option value="0.05">≤ 0.05</option></select></label>
    <label>tcrdist3 <select id="f_td"><option value="">any</option><option value="0">also</option><option value="1">partly</option><option value="2">not</option></select></label>
    <label>Top state <select id="f_st"><option value="">any</option>${ST.map(s => `<option>${s}</option>`).join("")}</select></label>
    <label>Search <input id="f_tx" placeholder="receptor id · V gene · CDR3 substring" style="width:240px"></label></div><p class="small muted" id="f_arm"></p><div id="tbl"></div>`);
  $("#clf").prepend(armSelect(() => draw()));
  $("#f_q").value = String(LS.maxQ); $("#f_st").value = LS.state; $("#f_tx").value = LS.text;
  const draw = () => {
    LS.minSize = +$("#f_min").value || 2; LS.maxQ = +$("#f_q").value; LS.state = $("#f_st").value; LS.text = $("#f_tx").value.trim().toUpperCase();
    const a = CUR; $("#f_arm").textContent = ARMS[a].name;
    const td = $("#f_td").value;
    let rows = [...cinfo[a].values()].filter(c => c.size >= LS.minSize && (LS.maxQ >= 1 || (c.q != null && c.q <= LS.maxQ)) && (!LS.state || c.top === LS.state) && (td === "" || (TROW(a, c.id) && TROW(a, c.id).r[3] === +td)));
    if (LS.text) rows = rows.filter(c => members[a].get(c.id).some(i => [M.id[i], M.vA[i], M.jA[i], M.vB[i], M.jB[i], M.c3A[i], M.c3B[i]].some(s => String(s).toUpperCase().includes(LS.text))));
    const tb = $("#tbl"); tb.innerHTML = "";
    tb.append(table([
      {h: "Cluster", v: r => r.id, num: 1}, {h: "Molecules", v: r => r.size, num: 1}, {h: "With state", v: r => r.nstate, num: 1},
      {h: "Top state", v: r => r.top, f: r => r.top ? `<span class="dot" style="background:var(${stateVar(r.top)})"></span>${r.top}` : "–"},
      {h: "Share", v: r => r.frac, f: r => pct(r.frac), num: 1}, {h: "Mice", v: r => r.mice, num: 1}, {h: "V pairs", v: r => r.vp, num: 1},
      {h: "tcrdist3", v: r => TROW(a, r.id) ? TROW(a, r.id).r[2] : null, f: r => TROW(a, r.id) ? `${tBadge(TROW(a, r.id).r[3])} <span class="muted small">${pct(TROW(a, r.id).r[2])}</span>` : "–", num: 1},
      {h: "q (BH)", v: r => r.q, f: r => (r.q != null && r.q <= 0.15 ? "<b>" + fq(r.q) + "</b>" : fq(r.q)), num: 1}],
      rows, {sort: 1, dir: -1, per: 100, click: k => location.hash = `#c/${a}/${k}`, key: r => r.id}));
  };
  fl.querySelectorAll("select,input").forEach(e => { if (e.id !== "armsel") e.oninput = draw; }); draw();
}
// ---------- CLUSTER PAGE
function viewCluster(app, arm, id) {
  if (!ARMS[arm] || !cinfo[arm].has(id)) { sec(app, `<p>No such cluster. <a href="#clusters">Back to the list</a></p>`); return; }
  const c = cinfo[arm].get(id), mem = members[arm].get(id), list = D.clusters[arm], prev = list[c.rank - 1], next = list[c.rank + 1];
  sec(app, `<div class="small"><a href="#clusters">← Clusters</a> · ${prev ? `<a href="#c/${arm}/${prev[0]}">‹ previous (${prev[0]})</a>` : ""} ${next ? ` · <a href="#c/${arm}/${next[0]}">next (${next[0]}) ›</a>` : ""}</div>
   <h1>Cluster ${id} <span class="muted" style="font-weight:400">· ${ARMS[arm].short}</span></h1><p class="muted">${ARMS[arm].name}</p>`);
  const sig = c.q != null && c.q <= 0.15;
  sec(app, `<div class="grid g4">
   <div class="card tile"><div class="lab">Molecules</div><div class="val">${c.size}</div><div class="sub">${c.nstate} with a state · ${mem.reduce((s, i) => s + M.cells[i], 0)} cells</div></div>
   <div class="card tile"><div class="lab">Top state</div><div class="val">${c.top ? `<span class="dot" style="background:var(${stateVar(c.top)})"></span>${c.top} ${pct(c.frac)}` : "–"}</div><div class="sub">binomial vs repertoire base rate</div></div>
   <div class="card tile"><div class="lab">q (BH; ≤ 0.15 = significant)</div><div class="val" style="color:${sig ? "var(--good)" : "inherit"}">${fq(c.q)}</div><div class="sub">p ${fq(c.p)}</div></div>
   <div class="card tile"><div class="lab">Diversity</div><div class="val">${c.vp} V pair${c.vp === 1 ? "" : "s"}</div><div class="sub">${c.mice} mice</div></div></div>`);
  const others = Object.keys(ARMS).filter(a => a !== arm).map(o => [o, [...new Set(mem.map(i => M[o][i]).filter(l => l >= 0))]]);
  const sc = ST.map((s, k) => mem.filter(i => M.st[i] === k).length), ns = sc.reduce((a, b) => a + b, 0) || 1;
  card(app, `<h3 style="margin-top:0">State composition (molecules)</h3><div class="bar">${sc.map((n, k) => n ? `<div title="${ST[k]}: ${n}" style="width:${100 * n / ns}%;background:var(${stateVar(ST[k])})"></div>` : "").join("")}</div>
   <div class="legend">${ST.map((s, k) => `<span><span class="dot" style="background:var(${stateVar(s)})"></span>${s} ${sc[k]}</span>`).join("")}</div>
   ${others.map(([o, ol]) => `<p class="small muted" style="margin:2px 0">In ${ARMS[o].short}, these members fall in ${ol.length ? ol.slice(0, 12).map(l => `<a href="#c/${o}/${l}">${l}</a>`).join(", ") + (ol.length > 12 ? " …" : "") : "no cluster (all singletons)"}.</p>`).join("")}`);
  sec(app, `<div class="grid g2">
    <div class="card"><h3 style="margin-top:0">Cells on the transcriptional UMAP</h3><div class="canvaswrap"><canvas id="um" height="440"></canvas><div class="tip" id="umtip"></div></div>
     <div class="legend">${ST.map(s => `<span><span class="dot" style="background:var(${stateVar(s)})"></span>${s}</span>`).join("")}<span><span class="dot" style="background:var(--cellbg)"></span>other cells</span></div></div>
    <div class="card"><h3 style="margin-top:0">Loops and voxel grids in 3D <span class="small muted">(drag · scroll · hover)</span></h3>
     <div class="filters" id="v3ctl"></div>
     <div id="v3m" style="position:relative;width:100%;height:440px;border:1px solid var(--line);border-radius:8px;overflow:hidden"></div>
     <div class="chk" id="loopchk"></div><div class="legend" id="vleg"></div><div class="small" id="vxnote"></div><div class="small muted" id="v3note"></div>
     <p class="small muted">Coordinates are the method's own frame: each chain fitted on its five framework landmarks onto the same chain of the reference model, then the TCR-intrinsic axes — exactly the placement the grid compares. Framework thin, loops thick, side chains as sticks; colour = V pair.
     <b>Voxels:</b> the members' ${ARMS[arm].atoms === "cdr3" ? "CDR3" : "all-loop"} grids superimposed — each member's grid built from its ${ARMS[arm].atoms === "cdr3" ? "CDR3 heavy atoms (IMGT 105–117; the other loops are drawn for context only)" : "loop heavy atoms (CDR1, CDR2, HV4, CDR3)"} exactly as the arm builds it (σ ${ARMS[arm].sigma} Å Gaussians, ${D.vox.h} Å voxels, the α and β boxes of the arm), then averaged over the shown members, and drawn where the average exceeds the chosen level (a percentage of its peak). Choose one member to see its own grid, or a chemistry channel. The grids are rebuilt in the page from the atoms in report_data/; gate RG4 checks them against the Python build and the stored distances.</p></div></div>`);
  drawUmap($("#um"), $("#umtip"), mem);
  if (window.$3Dmol) { viewer3dmol($("#v3m"), mem, $("#loopchk"), $("#vleg"), $("#v3ctl"), arm); $("#v3note").textContent = "Rendered with 3Dmol.js (Rego & Koes 2015)."; }
  else $("#v3m").innerHTML = '<p class="small muted" style="padding:10px">3Dmol.js unavailable (no WebGL?).</p>';
  tcrdPanel(app, arm, id, mem);
  sec(app, `<h2>Members</h2>`);
  app.append(table([
    {h: "Clone", v: i => M.id[i], cls: "mono"}, {h: "TRAV", v: i => M.vA[i]}, {h: "TRAJ", v: i => M.jA[i]}, {h: "CDR3α", v: i => M.c3A[i], cls: "mono"},
    {h: "TRBV", v: i => M.vB[i]}, {h: "TRBJ", v: i => M.jB[i]}, {h: "CDR3β", v: i => M.c3B[i], cls: "mono"},
    {h: "State", v: i => M.st[i] >= 0 ? ST[M.st[i]] : "", f: i => M.st[i] >= 0 ? `<span class="dot" style="background:var(${stateVar(ST[M.st[i]])})"></span>${ST[M.st[i]]}` : "–"},
    {h: "Mouse", v: i => M.mouse[i]}, {h: "Cells", v: i => M.cells[i], num: 1}],
    mem.slice(), {sort: 3, dir: 1, scroll: false}));
  splitPanel(app, arm, mem);
  alignmentPanel(app, arm, mem);
}
