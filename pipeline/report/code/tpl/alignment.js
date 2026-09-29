
// ---------- EXACT DISTANCE SPLIT and IMGT ALIGNMENT (cluster page)
// Split: D^2 is a sum over the two chain boxes and seven channels, so it splits exactly into 14 parts. For each member the
// page shows the mean over the other members of each part (stage 2 of the report build; checked to 1.2e-7 against D).
function splitPanel(app, arm, mem) {
  const S = SPLIT(arm), cut = D.thr[arm].bg_p1;
  const rows = mem.map(i => { const p = Array.from(S.slice(14 * i, 14 * i + 14)); return {i, p, t: p.reduce((s, x) => s + x, 0)}; });
  const avg = Array.from({length: 14}, (_, k) => rows.reduce((s, r) => s + r.p[k], 0) / rows.length), tavg = avg.reduce((s, x) => s + x, 0);
  const scale = Math.max(cut * cut, ...rows.map(r => r.t));
  const bar = (p, t) => `<div class="bar" style="height:14px;width:${Math.max(2, 100 * t / scale)}%">${p.map((v, k) => v > 0 ? `<div title="${k < 7 ? "α" : "β"} ${CHN[k % 7]}: ${f(v, 4)} (${pct(v / t)})" style="width:${100 * v / t}%;background:var(--v${k % 7 + 1});${k === 7 ? "border-left:2px solid var(--surface)" : ""}"></div>` : "").join("")}</div>`;
  const share = (p, t, lo) => pct(p.slice(lo, lo + 7).reduce((s, x) => s + x, 0) / t);
  const box = sec(app, `<h2>Where the members differ: exact split of the distance</h2>
   <p class="small muted">For each member: the mean squared distance to the other members of this cluster, split into its 14 exact parts — α chain (left of the white line) then β chain, each by chemistry channel. Bar length is relative to the cut² (${f(cut * cut, 3)}); every pair in a complete-linkage cluster is within the cut. RMS = √(mean squared distance).</p>
   ${chnLegend()}<div id="sp_tab"></div>`);
  const t = table([
    {h: "Clone", v: r => r.i === -1 ? "" : M.id[r.i], f: r => r.i === -1 ? "<b>cluster mean</b>" : `<span class="mono">${M.id[r.i]}</span>`},
    {h: "RMS to the others", v: r => Math.sqrt(r.t), f: r => f(Math.sqrt(r.t), 3), num: 1},
    {h: "α share", v: r => r.p.slice(0, 7).reduce((s, x) => s + x, 0) / r.t, f: r => share(r.p, r.t, 0), num: 1},
    {h: "Occupancy share", v: r => (r.p[0] + r.p[7]) / r.t, f: r => pct((r.p[0] + r.p[7]) / r.t), num: 1},
    {h: "Split (α | β, by channel)", v: r => r.t, f: r => bar(r.p, r.t)}],
    [{i: -1, p: avg, t: tavg}, ...rows.slice().sort((x, y) => y.t - x.t)], {scroll: false});   // cluster mean first, then members by RMS
  $("#sp_tab", box).append(t);
}
// Alignment: IMGT-numbered V-domain sequences read from the folded models (the reference method/aln/ALN1_sequences.npz), shaded by
// identity to the column's most common residue. The voxel distance has no per-residue decomposition, so no per-position
// structural score is shown; the exact structural split is the panel above.
const ALN_BY_MOL = new Map(), ALN_OK = {}, ALN_WAIT = {};
window.__dgAln = (id, rows) => { rows.forEach(r => ALN_BY_MOL.set(r.m, r)); ALN_OK[id] = true; (ALN_WAIT[id] || []).forEach(f => f[0]()); delete ALN_WAIT[id]; };
function loadAln(mis) {
  const need = [...new Set(mis.map(i => D.atoms.chunk[i]).filter(c => c >= 0 && !ALN_OK[c]))];
  return Promise.all(need.map(c => new Promise((res, rej) => {
    const first = !ALN_WAIT[c]; (ALN_WAIT[c] = ALN_WAIT[c] || []).push([res, rej]);
    if (first) { const s = document.createElement("script"); s.src = "report_data/aln_" + String(c).padStart(3, "0") + ".js";
      s.onerror = () => { (ALN_WAIT[c] || []).forEach(f => f[1](c)); delete ALN_WAIT[c]; s.remove(); }; document.head.append(s); }
  })));
}
const REGIONS = [["FR1", 1, 26], ["CDR1", 27, 38], ["FR2", 39, 55], ["CDR2", 56, 65], ["FR3", 66, 104], ["CDR3", 105, 117], ["FR4", 118, 128]];
const regionOf = n => (REGIONS.find(([, a, b]) => n >= a && n <= b) || ["?"])[0];
const isHV4 = n => n >= 81 && n <= 86;
const segOf = (ch, n) => n <= 104 ? (ch === "A" ? "TRAV" : "TRBV") : n <= 117 ? "junction" : (ch === "A" ? "TRAJ" : "TRBJ");
const ORIG_COL = {V: "--v3", J: "--v7"};
function juncOrigin(mi, ch, cols) {             // origin of each junction residue (IMGT 104-118) from the C1z map (as the vector report)
  const o = D.aln.origin[mi], nv = o[ch === "A" ? 0 : 2], nj = o[ch === "A" ? 1 : 3], J = cols.filter(c => D.aln.cols[c][0] === ch && D.aln.cols[c][1] >= 104 && D.aln.cols[c][1] <= 118), out = new Map();
  J.forEach((c, i) => { const v = nv >= 0 && i < nv, j = nj >= 0 && i >= J.length - nj;
    out.set(c, v && j ? "V/J" : v ? "V" : j ? "J" : (nv < 0 || nj < 0) ? "?" : "N"); });
  return out;
}
const ORIG_TXT = {V: "V-encoded", J: "J-encoded", N: "junctional (N/P/D additions)", "V/J": "germline, V and J extents overlap", "?": "origin unknown (gene with < 10 receptors)"};
const isDark = () => (document.documentElement.dataset.theme || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light")) === "dark";
function rampColour(ramp, x) {
  const R = isDark() ? ramp.slice().reverse() : ramp, t = Math.max(0, Math.min(1, x)) * (R.length - 1), k = Math.min(R.length - 2, Math.floor(t)), u = t - k;
  const h = s => [1, 3, 5].map(i => parseInt(s.slice(i, i + 2), 16)), a = h(R[k]), b = h(R[k + 1]);
  return "#" + a.map((v, i) => Math.round(v + (b[i] - v) * u).toString(16).padStart(2, "0")).join("");
}
function alnRow(mi) {
  const r = ALN_BY_MOL.get(mi); if (!r) return null;
  const cols = b64(r.c, Uint16Array), byCol = new Map();
  for (let k = 0; k < cols.length; k++) byCol.set(cols[k], {aa: r.a[k]});
  const colsArr = Array.from(cols), orig = new Map([...juncOrigin(mi, "A", colsArr), ...juncOrigin(mi, "B", colsArr)]);
  return {mi, byCol, orig};
}
function alignmentPanel(app, arm, mem) {
  const box = sec(app, `<h2>Sequence alignment (IMGT)</h2>
   <p class="small muted">IMGT-numbered V-domain sequences of the members, read from the folded models. Shading: grey = differs from the column's most common residue (identity only, not a structural score); · = no residue at this IMGT position. CDR3 underline: aqua = V-encoded, violet = J-encoded, none = junctional (per-gene map C1z, as in the vector report). ${ARMS[arm].atoms === "cdr3" ? "This arm's grid contains CDR3 only." : "This arm's grid contains CDR1, CDR2, HV4 and CDR3."}</p>
   <div class="filters"><label>Region <select id="al_reg"><option value="cdr">CDR1 · CDR2 · HV4 · CDR3 · J region</option><option value="all">whole V domain</option><option value="cdr3j">CDR3 + J region (IMGT 104–128)</option><option value="cdr3">CDR3 junction (IMGT 104–118)</option></select></label>
   <span class="small muted" id="al_stat">loading sequences…</span></div><div id="al_note"></div><div id="al_res"></div><div class="tip" id="al_tip" style="position:fixed"></div>`);
  $("#al_reg").value = ARMS[arm].atoms === "cdr3" ? "cdr3" : "cdr";
  loadAln(mem).then(() => {
    const rows = mem.map(alnRow).filter(Boolean);
    if (!rows.length) { $("#al_stat").textContent = "no sequence data for these molecules"; return; }
    $("#al_stat").textContent = "";
    const len = r => ["A", "B"].map(ch => [...r.byCol.keys()].filter(c => D.aln.cols[c][0] === ch && D.aln.cols[c][1] >= 105 && D.aln.cols[c][1] <= 117).length).join("-");
    const lens = [...new Set(rows.map(len))], notes = [];
    if (lens.length > 1) notes.push(`<b>CDR3 lengths differ within this cluster</b> (${lens.join(", ")} residues, α-β). The grid keeps loop length: a longer loop fills more voxels.`);
    const jn = rows.filter(r => D.aln.junction_notes.includes(M.id[r.mi])).map(r => M.id[r.mi]);
    if (jn.length) notes.push(`${jn.join(", ")}: the CDR3 recorded by the sequencing annotation differs from the IMGT numbering of the folded sequence; the alignment follows the folded sequence.`);
    $("#al_note").innerHTML = notes.map(s => `<div class="warnbox small" style="margin:6px 0">${s}</div>`).join("");
    const draw = () => {
      const reg = $("#al_reg").value, grey = css("--cellbg"), mute = css("--ink-3");
      const lab = r => `<span class="lc id mono">${M.id[r.mi]} <span class="dot" style="background:var(${M.st[r.mi] >= 0 ? stateVar(ST[M.st[r.mi]]) : "--cellbg"})"></span></span>`;
      const genes = r => `<span class="lc va">${esc(M.vA[r.mi])}</span><span class="lc ja">${esc(M.jA[r.mi])}</span><span class="lc vb">${esc(M.vB[r.mi])}</span><span class="lc jb">${esc(M.jB[r.mi])}</span>`;
      const labHead = `<span class="lc id">clone · state</span><span class="lc va">TRAV</span><span class="lc ja">TRAJ</span><span class="lc vb">TRBV</span><span class="lc jb">TRBJ</span>`;
      const oline = o => o === "V" || o === "J" ? `box-shadow:inset 0 -3px 0 var(${ORIG_COL[o]});` : o === "V/J" ? `box-shadow:inset 0 -3px 0 var(--ink-3);` : "";
      const inReg = c => { const n = D.aln.cols[c][1]; return reg === "all" || (reg === "cdr3" ? n >= 104 && n <= 118 : reg === "cdr3j" ? n >= 104 : ["CDR1", "CDR2", "CDR3", "FR4"].includes(regionOf(n)) || isHV4(n)); };
      const cols = [...new Set(rows.flatMap(r => [...r.byCol.keys()]))].filter(inReg).sort((a, b) => a - b);
      const cons = new Map(cols.map(c => { const t = {}; rows.forEach(r => { const e = r.byCol.get(c); if (e) t[e.aa] = (t[e.aa] || 0) + 1; });
        const top = Object.entries(t).sort((x, y) => y[1] - x[1])[0]; return [c, {aa: top[0], n: top[1]}]; }));
      const runLen = {}; { let key0 = null, start = 0; cols.forEach((c, i) => { const [ch, n] = D.aln.cols[c], k = ch + (isHV4(n) ? "HV4" : regionOf(n));
        if (k !== key0) { key0 = k; start = i; } runLen[start] = (runLen[start] || 0) + 1; }); }
      const segRun = {}; { let key0 = null, start = 0; cols.forEach((c, i) => { const [ch, n] = D.aln.cols[c], k = ch + segOf(ch, n);
        if (k !== key0) { key0 = k; start = i; } segRun[start] = (segRun[start] || 0) + 1; }); }
      let prevCh = null, prevReg = null, prevSeg = null, segStart = 0, sband = "", band = "", nums = "", conr = "", ci = -1, runStart = 0;
      for (const c of cols) { ci++; const [ch, n, ins] = D.aln.cols[c], rg = isHV4(n) ? "HV4" : regionOf(n), sep = prevCh && ch !== prevCh ? '<span class="alsep"></span>' : "";
        const starts = rg !== prevReg || ch !== prevCh; if (starts) runStart = ci;
        const sg = segOf(ch, n), sStarts = sg !== prevSeg || ch !== prevCh; if (sStarts) segStart = ci;
        const sgc = sg.endsWith("V") ? "--v3" : sg.endsWith("J") ? "--v7" : "--line";
        sband += sep + `<span class="ac hd sg" style="box-shadow:inset 0 -3px 0 var(${sgc})">${sStarts && segRun[segStart] >= 3 ? (sg === "junction" ? (ch === "A" ? "CDR3α junction" : "CDR3β junction") : sg) : ""}</span>`; prevSeg = sg;
        band += sep + `<span class="ac hd rg ${rg.startsWith("CDR") ? "cdr" : ""}">${starts && runLen[runStart] >= 3 ? (ch === "A" ? "α " : "β ") + rg : ""}</span>`;
        nums += sep + `<span class="ac hd">${ins ? ins : (n % 5 === 0 ? n : "")}</span>`;
        const k = cons.get(c); conr += sep + `<span class="ac" style="color:${k.n === rows.length ? "inherit" : mute}">${k.aa}</span>`;
        prevCh = ch; prevReg = rg; }
      let h = `<div class="alg"><div class="alrow alh"><span class="allab small muted">gene segment</span>${sband}</div><div class="alrow alh"><span class="allab small muted">IMGT region</span>${band}</div><div class="alrow alh"><span class="allab small muted">IMGT</span>${nums}</div><div class="alrow"><span class="allab small muted">${labHead}</span></div><div class="alrow"><span class="allab"><b>consensus</b></span>${conr}</div>`;
      for (const r of rows) { let s = "", pc = null;
        for (const c of cols) { const [ch] = D.aln.cols[c], e = r.byCol.get(c), sep = pc && ch !== pc ? '<span class="alsep"></span>' : ""; pc = ch;
          if (!e) { s += sep + `<span class="ac" style="color:${mute}">·</span>`; continue; }
          const diff = cons.get(c).aa !== e.aa; s += sep + `<span class="ac" data-r="${r.mi}" data-c="${c}" style="${diff ? `background:${grey};` : ""}${oline(r.orig.get(c))}">${e.aa}</span>`; }
        h += `<div class="alrow"><span class="allab">${lab(r)} ${genes(r)}</span>${s}</div>`; }
      $("#al_res").innerHTML = h + "</div>";
    };
    const tip = $("#al_tip"), byMi = new Map(rows.map(r => [r.mi, r]));
    box.addEventListener("mouseover", e => { const t = e.target.closest("[data-c]"); if (!t) { tip.style.display = "none"; return; }
      const r = byMi.get(+t.dataset.r), c = +t.dataset.c, [ch, n, ins] = D.aln.cols[c], e_ = r.byCol.get(c), sg = segOf(ch, n);
      const gene = sg === "TRAV" ? M.vA[r.mi] : sg === "TRAJ" ? M.jA[r.mi] : sg === "TRBV" ? M.vB[r.mi] : sg === "TRBJ" ? M.jB[r.mi] : "";
      const t_ = {}; rows.forEach(q => { const x = q.byCol.get(c); if (x) t_[x.aa] = (t_[x.aa] || 0) + 1; });
      tip.innerHTML = `<b class="mono">${M.id[r.mi]}</b> · ${ch === "A" ? "α" : "β"} IMGT ${n}${ins} · ${e_.aa} · ${isHV4(n) ? "HV4 (FR3)" : regionOf(n)}${sg === "junction" ? "" : ` · ${sg.endsWith("J") ? "J region" : "V region"} (${esc(gene)})`}${r.orig.has(c) ? `<br>${ORIG_TXT[r.orig.get(c)]}` : ""}<br><span class="muted">${Object.entries(t_).sort((x, y) => y[1] - x[1]).map(([a, k]) => `${a} ×${k}`).join(", ")} in this cluster</span>`;
      tip.style.display = "block"; });
    box.addEventListener("mousemove", e => { tip.style.left = Math.min(e.clientX + 14, innerWidth - 380) + "px"; tip.style.top = (e.clientY + 14) + "px"; });
    box.addEventListener("mouseleave", () => tip.style.display = "none");
    $("#al_reg").onchange = draw; draw();
  }).catch(() => { $("#al_stat").textContent = "sequence files not found — keep the report_data/ folder next to this report"; });
}
