
// condition labels for the binary cells flag C.gf (config report.condition_labels; [flag true, flag false])
const COND = (D.condition_labels && D.condition_labels.length === 2) ? D.condition_labels : ["group A", "group B"];
// generated from the reference report's UMAP page: PRIM -> CUR (the selected arm), plus the arm selector

// ---------- UMAP WITH CLUSTERS OVERLAID (+ colour by transcriptional state / cluster / gene expression)
const UM = {mode: "sig", minSize: 3, state: "", focus: null, col: "struct", gene: "", rings: true};
function clusterColour(k) {                       // distinct hues by golden-angle spacing; lightness follows the theme
  return `hsl(${Math.round((k * 137.508) % 360)}, 70%, ${isDark() ? 62 : 44}%)`;
}
// transcriptional clusters: each state keeps its hue; where a state holds several clusters, the largest keeps the state
// colour and the others are a lighter / darker step of the same hue (OKLab L +-0.12..0.16; each family validated as a
// monotone single-hue ramp). Cluster identity is carried by the direct labels on the plot, not by colour alone.
const TCL_SHADE = {light: {6: "#b33300", 5: "#ff8f5c", 7: "#0047a1", 10: "#57a4ff"}, dark: {6: "#a82a00", 5: "#ff8656", 7: "#0056af", 10: "#66b4ff"}};
const TCL = new Map(D.tcl.map(t => [t.id, t]));
const tclColour = id => TCL_SHADE[isDark() ? "dark" : "light"][id] || css(stateVar(TCL.get(id).state));
// gene expression: log1p(CP10k), sequential blue ramp from its second step (so the lowest expressing cells stay clear of
// the grey "not detected" cells), capped at the gene's 99th percentile among expressing cells
const EXPR_RAMP = ["#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"];
const GENE_IX = new Map(D.genes.g.map((g, k) => [g.toUpperCase(), k]));
const GENE_VAL = new Map(), GENE_WAIT = {};
window.__dgGene = (k, o) => {
  const idx = b64(o.i, Uint16Array), q = b64(o.q, Uint8Array), v = new Float32Array(D.cells.x.length), s = D.genes.vmax[k] / 255;
  for (let j = 0; j < idx.length; j++) v[idx[j]] = q[j] * s;
  GENE_VAL.set(k, v); (GENE_WAIT[k] || []).forEach(f => f[0](v)); delete GENE_WAIT[k];
};
function loadGene(k) {
  if (GENE_VAL.has(k)) return Promise.resolve(GENE_VAL.get(k));
  return new Promise((res, rej) => { const first = !GENE_WAIT[k]; (GENE_WAIT[k] = GENE_WAIT[k] || []).push([res, rej]);
    if (first) { const s = document.createElement("script"); s.src = "report_data/expr/g_" + String(k).padStart(5, "0") + ".js";
      s.onerror = () => { (GENE_WAIT[k] || []).forEach(f => f[1](k)); delete GENE_WAIT[k]; s.remove(); }; document.head.append(s); } });
}
let GENE_OPTS = null;
function viewUmap(app) {
  if (!GENE_OPTS) GENE_OPTS = D.genes.g.map(g => `<option value="${esc(g)}">`).join("");
  sec(app, `<h1>Transcriptional UMAP with structural clusters</h1>
  <p class="muted">Every cell on the transcriptional UMAP. Colour the cells by structural cluster (of the arm chosen below), by transcriptional state, by transcriptional (Leiden) cluster, or by the expression of any detected gene. In the last three views the cells of the selected structural clusters are ringed. Hover a cell for details; click a structural-cluster cell to open its cluster. Click a structural cluster in the list to highlight it.</p>
  <div class="filters">
   <label>Colour cells by <select id="um_col"><option value="struct">structural clusters</option><option value="state">transcriptional state (6)</option><option value="tcl">transcriptional cluster (10)</option><option value="gene">gene expression</option></select></label>
   <label id="um_gw">Gene <input id="um_gene" list="um_gl" placeholder="e.g. Foxp3, Rorc, Tbx21" style="width:170px" autocomplete="off"><datalist id="um_gl">${GENE_OPTS}</datalist></label>
   <label id="um_rw"><input type="checkbox" id="um_ring"> ring structural-cluster cells</label></div>
  <div class="filters">
   <label>Structural clusters <select id="um_mode"><option value="sig">significant (q ≤ 0.15)</option><option value="all">all</option></select></label>
   <label id="um_minw">Min molecules <input id="um_min" type="number" min="2" value="${UM.minSize}" style="width:70px"></label>
   <label>Top state <select id="um_st"><option value="">any</option>${ST.map(s => `<option>${s}</option>`).join("")}</select></label>
   <span class="small muted" id="um_info"></span></div>
  <div class="card"><div class="canvaswrap"><canvas id="umc" height="640"></canvas><div class="tip" id="umt"></div></div>
   <div class="legend" id="um_key"></div><div class="small muted" id="um_gnote"></div>
   <div class="legend" id="um_leg" style="max-height:170px;overflow:auto"></div>
   <p class="small muted" id="um_foot"></p></div>
  <h2>Structural clusters shown</h2><div class="small muted" id="um_tnote"></div><div id="um_tab"></div>`);
  document.querySelector("#app .filters").prepend(armSelect(() => select()));        // arm selector
  $("#um_mode").value = UM.mode; $("#um_st").value = UM.state; $("#um_col").value = UM.col; $("#um_gene").value = UM.gene; $("#um_ring").checked = UM.rings;
  const cv = $("#umc"), tip = $("#umt"), C = D.cells, NC = C.x.length;
  let x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
  for (let i = 0; i < NC; i++) { x0 = Math.min(x0, C.x[i]); x1 = Math.max(x1, C.x[i]); y0 = Math.min(y0, C.y[i]); y1 = Math.max(y1, C.y[i]); }
  let shown = [], pts = [], X, Y, W, H = 640, gk = null, gv = null, gmsg = "";
  const tclOrder = D.tcl.map(t => t.id).sort((a, b) => ST.indexOf(TCL.get(a).state) - ST.indexOf(TCL.get(b).state) || TCL.get(b).n - TCL.get(a).n);
  function select() {
    UM.mode = $("#um_mode").value; UM.minSize = +$("#um_min").value || 2; UM.state = $("#um_st").value; UM.col = $("#um_col").value; UM.rings = $("#um_ring").checked;
    $("#um_minw").style.display = UM.mode === "all" ? "" : "none";
    $("#um_gw").style.display = UM.col === "gene" ? "" : "none"; $("#um_rw").style.display = UM.col === "struct" ? "none" : "";
    let cs = [...cinfo[CUR].values()];
    if (UM.mode === "sig") cs = cs.filter(c => c.q != null && c.q <= 0.15);
    else cs = cs.filter(c => c.size >= UM.minSize);
    if (UM.state) cs = cs.filter(c => c.top === UM.state);
    cs.sort((a, b) => (a.q ?? 1) - (b.q ?? 1) || b.size - a.size);
    shown = cs.map((c, k) => { const cells = []; members[CUR].get(c.id).forEach(m => (cellsByMol.get(m) || []).forEach(i => cells.push(i)));
      return Object.assign({}, c, {col: clusterColour(k), cells}); });
    if (UM.focus != null && !shown.some(c => c.id === UM.focus)) UM.focus = null;
    pts = []; shown.forEach(c => c.cells.forEach(i => pts.push([i, c])));
    $("#um_info").textContent = `${shown.length} structural clusters · ${pts.length} cells`;
    geneThen();
  }
  function geneThen() {                            // resolve the gene (if any), then redraw everything
    gk = null; gv = null; gmsg = "";
    if (UM.col === "gene") {
      UM.gene = $("#um_gene").value.trim(); const k = GENE_IX.get(UM.gene.toUpperCase());
      if (!UM.gene) gmsg = "Type a gene name (" + D.genes.g.length.toLocaleString("en-US") + " genes detected in ≥ 0.5 % of cells).";
      else if (k == null) gmsg = `“${esc(UM.gene)}” is not among the detected genes (expressed in ≥ 0.5 % of cells).`;
      else if (!GENE_VAL.has(k)) { gmsg = "loading…"; paint();
        loadGene(k).then(() => { if (GENE_IX.get($("#um_gene").value.trim().toUpperCase()) === k) geneThen(); })
          .catch(() => { gmsg = "Expression file not found — keep the report_data/expr/ folder next to this report."; paint(); }); return; }
      else { gk = k; gv = GENE_VAL.get(k); }
    }
    paint();
  }
  function paint() { draw(); key(); legend(); tab(); }
  function sizeCanvas() {
    const dpr = devicePixelRatio || 1; W = cv.clientWidth; cv.width = W * dpr; cv.height = H * dpr; cv.style.height = H + "px";
    const g = cv.getContext("2d"); g.scale(dpr, dpr);
    const pad = 16, s = Math.min((W - 2 * pad) / (x1 - x0), (H - 2 * pad) / (y1 - y0));
    X = x => pad + (x - x0) * s + ((W - 2 * pad) - (x1 - x0) * s) / 2; Y = y => H - pad - (y - y0) * s - ((H - 2 * pad) - (y1 - y0) * s) / 2;
    g.fillStyle = css("--surface"); g.fillRect(0, 0, W, H); return g;
  }
  function label(g, txt, x, y) {                  // direct label with a surface halo (secondary encoding for the categories)
    g.font = "600 12px -apple-system,Segoe UI,Helvetica,Arial,sans-serif"; g.textAlign = "center"; g.textBaseline = "middle";
    g.lineWidth = 4; g.strokeStyle = css("--surface"); g.strokeText(txt, x, y); g.fillStyle = css("--ink"); g.fillText(txt, x, y);
  }
  function centroid(sel) { const xs = [], ys = []; for (let i = 0; i < NC; i++) if (sel(i)) { xs.push(C.x[i]); ys.push(C.y[i]); }
    const med = a => { a.sort((p, q) => p - q); return a[a.length >> 1]; }; return [X(med(xs)), Y(med(ys))]; }
  function draw() {
    const g = sizeCanvas(), bg = css("--cellbg"), labs = [];
    if (UM.col === "struct") {
      g.fillStyle = bg; for (let i = 0; i < NC; i++) g.fillRect(X(C.x[i]) - 0.9, Y(C.y[i]) - 0.9, 1.8, 1.8);
      const surf = css("--surface"), ink = css("--ink");
      const order = UM.focus == null ? pts : pts.filter(p => p[1].id !== UM.focus).concat(pts.filter(p => p[1].id === UM.focus));
      for (const [i, c] of order) {
        const px = X(C.x[i]), py = Y(C.y[i]), dim = UM.focus != null && c.id !== UM.focus, r = c.id === UM.focus ? 5 : 3.6;
        g.globalAlpha = dim ? 0.15 : 1;
        g.beginPath(); g.arc(px, py, r + 1, 0, 7); g.fillStyle = surf; g.fill();
        g.beginPath(); g.arc(px, py, r, 0, 7); g.fillStyle = c.col; g.fill();
      }
      g.globalAlpha = 1; return;
    }
    if (UM.col === "state") {
      const cols = ST.map(s => css(stateVar(s)));
      for (let i = 0; i < NC; i++) { g.fillStyle = cols[C.st[i]]; g.fillRect(X(C.x[i]) - 1.1, Y(C.y[i]) - 1.1, 2.2, 2.2); }
      ST.forEach((s, k) => labs.push([s, ...centroid(i => C.st[i] === k)]));
    } else if (UM.col === "tcl") {
      const cols = new Map(D.tcl.map(t => [t.id, tclColour(t.id)]));
      for (let i = 0; i < NC; i++) { g.fillStyle = cols.get(C.cl[i]); g.fillRect(X(C.x[i]) - 1.1, Y(C.y[i]) - 1.1, 2.2, 2.2); }
      D.tcl.forEach(t => labs.push([`${t.state} · ${t.id}`, ...centroid(i => C.cl[i] === t.id)]));
    } else {
      g.fillStyle = bg; for (let i = 0; i < NC; i++) if (!gv || gv[i] === 0) g.fillRect(X(C.x[i]) - 0.9, Y(C.y[i]) - 0.9, 1.8, 1.8);
      if (gv) { const cap = D.genes.cap[gk] || D.genes.vmax[gk], ord = []; for (let i = 0; i < NC; i++) if (gv[i] > 0) ord.push(i);
        ord.sort((a, b) => gv[a] - gv[b]);                 // highest on top
        for (const i of ord) { g.fillStyle = rampColour(EXPR_RAMP, gv[i] / cap); g.fillRect(X(C.x[i]) - 1.3, Y(C.y[i]) - 1.3, 2.6, 2.6); } }
    }
    if (UM.rings) {                                         // rings: cells of the selected structural clusters
      const ink = css("--ink");
      for (const [i, c] of pts) { const foc = c.id === UM.focus, dim = UM.focus != null && !foc;
        g.globalAlpha = dim ? 0.2 : 1; g.beginPath(); g.arc(X(C.x[i]), Y(C.y[i]), foc ? 5 : 3.4, 0, 7); g.strokeStyle = ink; g.lineWidth = foc ? 1.8 : 1; g.stroke(); }
      g.globalAlpha = 1;
    }
    labs.forEach(([t, x, y]) => label(g, t, x, y));        // labels last, above the rings
  }
  function key() {                                 // colour key for the chosen colouring
    const k = $("#um_key"), n = $("#um_gnote");
    n.innerHTML = UM.col === "gene" ? gmsg : "";
    if (UM.col === "state") k.innerHTML = ST.map(s => `<span><span class="dot" style="background:var(${stateVar(s)})"></span>${s} · ${C.st.filter(x => x === ST.indexOf(s)).length} cells</span>`).join("");
    else if (UM.col === "tcl") k.innerHTML = tclOrder.map(id => { const t = TCL.get(id); return `<span><span class="dot" style="background:${tclColour(id)}"></span>${t.state} · ${id} <span class="muted">(${t.n})</span></span>`; }).join("");
    else if (UM.col === "gene" && gv) {
      const cap = D.genes.cap[gk], R = isDark() ? EXPR_RAMP.slice().reverse() : EXPR_RAMP, n_ = D.genes.n[gk];
      k.innerHTML = `<span><b>${esc(D.genes.g[gk])}</b></span><span>0 <span style="display:inline-block;width:160px;height:10px;border-radius:2px;vertical-align:middle;background:linear-gradient(90deg,${R.join(",")})"></span> ≥ ${f(cap, 2)}</span><span><span class="dot" style="background:var(--cellbg)"></span>not detected</span>`;
      n.innerHTML = `log1p(counts per 10,000) · detected in ${n_.toLocaleString()} of ${NC.toLocaleString()} cells (${pct(n_ / NC)}) · colour saturates at ${f(cap, 2)} (99th percentile of the expressing cells; maximum ${f(D.genes.vmax[gk], 2)}).`;
    } else k.innerHTML = "";
    $("#um_foot").textContent = UM.col === "struct" ? "A molecule may contribute several cells (one clonotype, several cells)."
      : (UM.rings ? "Rings = cells of the structural clusters listed below (bold ring = highlighted cluster). " : "") + (UM.col === "gene" ? "" : "Labels mark the median position of each group" + (UM.col === "tcl" ? "; where a state has several clusters, the largest keeps the state colour and the others are a lighter / darker step of the same hue." : "."));
  }
  function legend() {
    const box = $("#um_leg"), ringStyle = UM.col !== "struct";
    box.style.display = UM.col === "struct" || UM.rings ? "" : "none";
    box.innerHTML = shown.map(c => `<span class="li ${UM.focus != null && UM.focus !== c.id ? "off" : ""}" data-c="${c.id}" style="cursor:pointer"><span class="dot" style="${ringStyle ? "background:transparent;box-shadow:inset 0 0 0 1.5px var(--ink)" : "background:" + c.col}"></span>${c.id} · ${c.top || "–"} · ${c.cells.length} cells</span>`).join("");
    box.querySelectorAll(".li").forEach(e => e.onclick = () => { const c = +e.dataset.c; UM.focus = UM.focus === c ? null : c; draw(); legend(); });
  }
  function tab() {
    const b = $("#um_tab"); b.innerHTML = "";
    const cols = [{h: "", v: r => 0, f: r => `<span class="dot" style="background:${r.col}"></span>`},
      {h: "Cluster", v: r => r.id, f: r => `<a href="#c/${CUR}/${r.id}">${r.id}</a>`, num: 1},
      {h: "Molecules", v: r => r.size, num: 1}, {h: "Cells on UMAP", v: r => r.cells.length, num: 1},
      {h: "Top state", v: r => r.top, f: r => r.top ? `<span class="dot" style="background:var(${stateVar(r.top)})"></span>${r.top} ${pct(r.frac)}` : "–"},
      {h: "Cells in the top state", v: r => r.cells.length ? r.cells.filter(i => ST[C.st[i]] === r.top).length / r.cells.length : null, f: r => r.cells.length ? pct(r.cells.filter(i => ST[C.st[i]] === r.top).length / r.cells.length) : "–", num: 1},
      {h: "Main transcriptional cluster", v: r => mainTcl(r)[0], f: r => { const [id, fr] = mainTcl(r); return id == null ? "–" : `${TCL.get(id).state} · ${id} ${pct(fr)}`; }}];
    if (gv) { const nm = esc(D.genes.g[gk]);
      cols.push({h: `${nm} % expressing`, v: r => r.cells.length ? r.cells.filter(i => gv[i] > 0).length / r.cells.length : null, f: r => r.cells.length ? pct(r.cells.filter(i => gv[i] > 0).length / r.cells.length) : "–", num: 1},
                {h: `${nm} mean`, v: r => r.cells.length ? r.cells.reduce((s, i) => s + gv[i], 0) / r.cells.length : null, f: r => r.cells.length ? f(r.cells.reduce((s, i) => s + gv[i], 0) / r.cells.length, 2) : "–", num: 1});
      $("#um_tnote").innerHTML = `${nm} in all ${NC.toLocaleString()} cells: ${pct(D.genes.n[gk] / NC)} expressing · mean ${f(D.genes.mean[gk], 2)} (log1p CP10k). The cluster columns use the cells of each cluster's molecules.`;
    } else $("#um_tnote").innerHTML = "";
    cols.push({h: "Mice", v: r => r.mice, num: 1}, {h: "V pairs", v: r => r.vp, num: 1}, {h: "q (BH)", v: r => r.q, f: r => fq(r.q), num: 1});
    b.append(table(cols, shown, {per: 100, scroll: false, click: k => location.hash = `#c/${CUR}/${k}`, key: r => r.id}));
  }
  function mainTcl(r) { if (!r._tcl) { const oc = {}; r.cells.forEach(i => oc[C.cl[i]] = (oc[C.cl[i]] || 0) + 1);
      const t = Object.entries(oc).sort((a, b) => b[1] - a[1])[0]; r._tcl = t ? [+t[0], t[1] / r.cells.length] : [null, 0]; } return r._tcl; }
  function nearest(e) {
    const r = cv.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top; let best = null, bd = 49;
    if (UM.col === "struct" || UM.rings)
      for (const p of pts) { if (UM.focus != null && p[1].id !== UM.focus) continue; const d = (X(C.x[p[0]]) - mx) ** 2 + (Y(C.y[p[0]]) - my) ** 2; if (d < bd) { bd = d; best = p; } }
    if (!best && UM.col !== "struct") { let bi = -1; bd = 25;
      for (let i = 0; i < NC; i++) { const d = (X(C.x[i]) - mx) ** 2 + (Y(C.y[i]) - my) ** 2; if (d < bd) { bd = d; bi = i; } }
      if (bi >= 0) best = [bi, null]; }
    return [best, mx, my];
  }
  cv.onmousemove = e => { const [p, mx, my] = nearest(e); if (!p) { tip.style.display = "none"; cv.style.cursor = ""; return; }
    const [i, c] = p, m = C.mi[i]; cv.style.cursor = c ? "pointer" : "";
    const cell = `${ST[C.st[i]]} cell · transcriptional cluster ${TCL.get(C.cl[i]).state} · ${C.cl[i]}` + (gv ? `<br>${esc(D.genes.g[gk])}: ${gv[i] > 0 ? f(gv[i], 2) : "not detected"}` : "");
    const mol = m >= 0 ? `<br><span class="mono">${M.id[m]}</span> · ${M.vA[m]} | ${M.vB[m]} · ${C.gf[i] ? COND[0] : COND[1]} · ${M.mouse[m]}` : "<br>no folded receptor";
    tip.innerHTML = (c ? `<b>structural cluster ${c.id}</b> · ${c.top} ${pct(c.frac)} · q ${fq(c.q)}<br>` : "") + cell + mol;
    tip.style.display = "block"; tip.style.left = Math.min(mx + 12, W - 280) + "px"; tip.style.top = (my + 12) + "px"; };
  cv.onmouseleave = () => tip.style.display = "none";
  cv.onclick = e => { const [p] = nearest(e); if (p && p[1]) location.hash = `#c/${CUR}/${p[1].id}`; };
  ["#um_mode", "#um_min", "#um_st", "#um_col", "#um_ring"].forEach(s => $(s).oninput = select);
  $("#um_gene").oninput = () => { if (GENE_IX.has($("#um_gene").value.trim().toUpperCase()) || !$("#um_gene").value.trim()) geneThen(); };
  $("#um_gene").onchange = geneThen;
  select();
}
