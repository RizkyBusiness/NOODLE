
// ---------- 3D viewer of a cluster: loops (Cα trace, side chains) with the superimposed voxel grids of the arm
// (replaces the source template's 3D viewer; same traces, side chains, hover labels, V-pair legend)
const VOXCOL = {occA: "#5f8fc4", occB: "#c4895f", 1: "#e69f00", 2: "#9467bd", 3: "#56b4e9", 4: "#d55e00", 5: "#0050c8", 6: "#cc3399"};
function viewer3dmol(div, mem, chkBox, legBox, ctl, arm) {
  const V = D.vox;
  const names = ["CDR1α", "CDR2α", "HV4α", "CDR3α", "CDR1β", "CDR2β", "HV4β", "CDR3β"]; const show = [true, true, false, true, true, true, false, true];
  const mols = mem.map(i => ({i, vp: M.vpA[i] + " | " + M.vpB[i]}));
  const vpCount = {}; mols.forEach(m => vpCount[m.vp] = (vpCount[m.vp] || 0) + 1);
  const vps = Object.keys(vpCount).sort((a, b) => vpCount[b] - vpCount[a] || (a < b ? -1 : 1));
  const vcol = {}; vps.forEach((v, k) => vcol[v] = css(k < 7 ? "--v" + (k + 1) : "--vother"));
  const hidden = new Set(); let mode = ARMS[arm].atoms === "cdr3" ? "cdr3" : "sticks", showFw = true;   // CDR3-only arms: CDR3 side chains by default
  const vx = {src: "mean", ch: "0", style: "surface", level: 30, alpha: 60};
  ctl.innerHTML = `<label>View <select id="v3mode"><option value="trace">Cα trace</option><option value="sticks"${mode === "sticks" ? " selected" : ""}>trace + side chains</option><option value="cdr3"${mode === "cdr3" ? " selected" : ""}>CDR3 side chains only</option></select></label>` +
    `<button id="v3reset">reset view</button><button id="v3png">save PNG</button><span class="small muted" id="v3status">loading Cα traces…</span>` +
    `<div class="filters" style="margin-top:6px"><label>Voxels <select id="vxsrc"><option value="mean">cluster mean</option><option value="off">off</option>${mem.map(i => `<option value="${i}">${esc(M.id[i])}</option>`).join("")}</select></label>` +
    `<label>Channel <select id="vxch">${CHN.map((n, c) => `<option value="${c}">${n}</option>`).join("")}<option value="chem">all six chemistry channels</option></select></label>` +
    `<label>Style <select id="vxsty"><option value="surface">surface</option><option value="voxel">voxel cubes</option><option value="mesh">mesh</option></select></label>` +
    `<label>Level <input id="vxlev" type="range" min="5" max="90" step="5" value="${vx.level}" style="width:90px"> <span id="vxlevv">${vx.level} %</span></label>` +
    `<label>Opacity <input id="vxal" type="range" min="10" max="100" step="5" value="${vx.alpha}" style="width:70px"></label></div>`;
  const viewer = $3Dmol.createViewer(div, {backgroundColor: css("--surface")});
  const hex = c => c.startsWith("#") ? c : "#888888";
  let models = [];
  const rangeSel = k => ({chain: LOOPDEF[k][0], resi: [LOOPDEF[k][1] + "-" + LOOPDEF[k][2]]});
  function hover(m, mm, isAtoms) { const lab = `${M.id[m.i]} · ${m.vp}`;
    mm.setHoverable({}, true, (atom, v) => { if (atom.__lab) return; const li = loopIndex(atom.chain, atom.resi);
      atom.__lab = v.addLabel(lab + " · " + (li >= 0 ? names[li] : "framework " + (atom.chain === "A" ? "α" : "β")) + " " + atom.resi + (atom.inscode || "") + (isAtoms ? " " + atom.resn : ""),
        {position: atom, backgroundColor: css("--surface"), fontColor: css("--ink"), borderColor: css("--line"), borderThickness: 1, fontSize: 11, inFront: true}); },
      (atom, v) => { if (atom.__lab) { v.removeLabel(atom.__lab); delete atom.__lab; } }); }
  function style() {
    viewer.setStyle({}, {});
    for (const m of models) {
      if (hidden.has(m.vp)) continue; const col = hex(vcol[m.vp]);
      if (showFw) m.ca.setStyle({}, {cartoon: {style: "trace", color: col, thickness: 0.08, opacity: 0.45}});
      LOOPDEF.forEach((d, k) => { if (!show[k]) return; const cdr3 = k === 3 || k === 7;
        m.ca.setStyle(rangeSel(k), {cartoon: {style: "trace", color: col, thickness: cdr3 ? 0.32 : 0.18}});
        if (m.aa && (mode === "sticks" || (mode === "cdr3" && cdr3)))
          m.aa.setStyle({...rangeSel(k), not: {atom: ["N", "C", "O", "OXT"]}}, {stick: {radius: 0.14, colorscheme: {prop: "elem", map: {C: col, N: "#3050F8", O: "#FF0D0D", S: "#E0C000"}}}});
      });
    }
    viewer.render();
  }
  // ---- voxel grids: built in the page from the side-car atoms (voxAdd), cached per member set
  let grid = null, gridKey = "", gridN = 0;
  function currentGrid() {
    const who = vx.src === "mean" ? models.filter(m => !hidden.has(m.vp)).map(m => m.i) : [+vx.src];
    const key = who.join(",");
    if (key !== gridKey) {
      grid = voxNew(V, arm); gridKey = key; gridN = who.length;
      who.forEach(i => { const E = AT_BY_MOL.get(i); if (E) voxAdd(grid, E, V, arm, 1 / who.length); });
    }
    return grid;
  }
  function drawVox() {
    viewer.removeAllShapes();
    const note = $("#vxnote");
    if (vx.src === "off" || !atomsLoaded || !models.length) { if (note) note.textContent = ""; viewer.render(); return; }
    const g = currentGrid(); if (!gridN) { note.textContent = "no member shown"; viewer.render(); return; }
    const chans = vx.ch === "chem" ? [1, 2, 3, 4, 5, 6] : [+vx.ch], txt = [];
    for (const c of chans) {
      const vols = ["A", "B"].map(k => ({k, ...voxVolume(g[k], V.boxes[arm][k], c, V.h)}));
      const peak = Math.max(vols[0].max, vols[1].max); if (!(peak > 0)) { txt.push(`${CHN[c]}: empty`); continue; }
      const iso = peak * vx.level / 100;
      for (const v of vols) {
        if (v.max <= iso) continue;
        viewer.addIsosurface(v.vol, {isoval: iso, color: c === 0 ? VOXCOL["occ" + v.k] : VOXCOL[c], alpha: vx.alpha / 100, smoothness: vx.style === "voxel" ? 0 : 1,
          voxel: vx.style === "voxel", wireframe: vx.style === "mesh"});
      }
      txt.push(`${CHN[c]}: level ${f(iso, 4)} of peak ${f(peak, 4)}`);
    }
    note.innerHTML = `Voxels: ${vx.src === "mean" ? `mean of the ${gridN} shown member grid${gridN === 1 ? "" : "s"}` : "grid of " + esc(M.id[+vx.src])} · ${txt.join(" · ")} (mass per 1 Å voxel)` +
      (vx.ch === "0" ? ` · <span style="color:${VOXCOL.occA}">■</span> α box <span style="color:${VOXCOL.occB}">■</span> β box` : vx.ch === "chem" ? " · " + [1, 2, 3, 4, 5, 6].map(c => `<span style="color:${VOXCOL[c]}">■</span> ${CHN[c]}`).join(" ") : "");
    viewer.render();
  }
  let atomsLoaded = false;
  function ensureAtoms() {
    if (atomsLoaded) return Promise.resolve(true); $("#v3status").textContent = "loading side chains…";
    return loadAtoms(models.map(m => m.i)).then(() => {
      models.forEach(m => { const p = atomPDBreal(m.i); if (p) { m.aa = viewer.addModel(p, "pdb"); hover(m, m.aa, true); } });
      atomsLoaded = true; $("#v3status").textContent = ""; return true;
    }).catch(() => { $("#v3status").textContent = "side-chain files not found — keep the report_data/ folder next to this report"; return false; });
  }
  const zoomLoops = () => { viewer.zoomTo({or: LOOPDEF.map((d, k) => rangeSel(k))}); viewer.render(); };
  loadCA(mols.map(m => m.i)).then(() => {
    models = mols.map(m => ({...m, ca: viewer.addModel(fullCaPDB(m.i), "pdb"), aa: null})).filter(m => m.ca);
    models.forEach(m => hover(m, m.ca, false)); $("#v3status").textContent = ""; style(); zoomLoops();
    return ensureAtoms().then(ok => { if (!ok) { mode = "trace"; $("#v3mode").value = "trace"; } style(); drawVox(); });
  }).catch(() => {
    $("#v3status").textContent = "3D files not found — keep the report_data/ folder next to this report";
  });
  chkBox.innerHTML = names.map((n, k) => `<label><input type="checkbox" data-k="${k}" ${show[k] ? "checked" : ""}> ${n}</label>`).join("") + `<label><input type="checkbox" id="fw" checked> framework</label>`;
  chkBox.querySelectorAll("input[data-k]").forEach(e => e.onchange = () => { show[+e.dataset.k] = e.checked; style(); });
  $("#fw").onchange = e => { showFw = e.target.checked; style(); };
  legBox.innerHTML = vps.map(v => `<span class="li" data-v="${esc(v)}"><span class="dot" style="background:${vcol[v]}"></span>${esc(v)} (${vpCount[v]})</span>`).join("") + (vps.length > 7 ? `<span class="small muted">grey = V pairs beyond the 7th</span>` : "") + `<span class="small muted">click a V pair to hide it (the cluster-mean voxels follow)</span>`;
  legBox.querySelectorAll(".li").forEach(e => e.onclick = () => { const v = e.dataset.v; if (hidden.has(v)) hidden.delete(v); else hidden.add(v); e.classList.toggle("off"); style(); drawVox(); });
  $("#v3mode").onchange = e => { const want = e.target.value; if (want === "trace") { mode = want; style(); return; }
    ensureAtoms().then(ok => { if (ok) { mode = want; style(); } else { e.target.value = "trace"; mode = "trace"; style(); } }); };
  $("#vxsrc").onchange = e => { vx.src = e.target.value; drawVox(); };
  $("#vxch").onchange = e => { vx.ch = e.target.value; drawVox(); };
  $("#vxsty").onchange = e => { vx.style = e.target.value; drawVox(); };
  $("#vxlev").oninput = e => { vx.level = +e.target.value; $("#vxlevv").textContent = vx.level + " %"; };
  $("#vxlev").onchange = () => drawVox();
  $("#vxal").onchange = e => { vx.alpha = +e.target.value; drawVox(); };
  $("#v3reset").onclick = zoomLoops;
  $("#v3png").onclick = () => { const a = document.createElement("a"); a.href = viewer.pngURI(); a.download = "cluster_loops_voxels.png"; a.click(); };
}
