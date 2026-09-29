
// ---------- WORKED EXAMPLE (a chosen receptor pair, in both arms, split exactly)
function splitBar(parts, scale) {                 // one stacked bar per chain: 7 channel parts of D^2
  const seg = (k, v) => v > 0 ? `<div title="${CHN[k % 7]}: ${f(v, 4)}" style="width:${100 * v / scale}%;background:var(--v${k % 7 + 1})"></div>` : "";
  return ["α", "β"].map((c, ci) => `<div class="xb" style="grid-template-columns:40px 1fr 70px"><span>${c}</span><div class="bar" style="height:16px">${parts.slice(7 * ci, 7 * ci + 7).map((v, j) => seg(7 * ci + j, v)).join("")}</div><span class="xbx">${f(parts.slice(7 * ci, 7 * ci + 7).reduce((s, x) => s + x, 0), 3)}</span></div>`).join("");
}
const chnLegend = () => `<div class="legend">${CHN.map((c, k) => `<span><span class="dot" style="background:var(--v${k + 1})"></span>${c}</span>`).join("")}</div>`;
function viewExample(app) {
  const X = D.example, a = X.ai, b = X.bi;
  sec(app, `<h1>Worked example</h1>
  <p class="muted">The same two receptors as in the vector method's report, compared by the voxel grid in every arm of this report. The distance is D = √(Σ parts); each part is one chain and one chemistry channel, and the parts add up exactly.</p>`);
  app.append(table([{h: "", v: r => r[0]}, {h: "Clone", v: r => r[1], cls: "mono"}, {h: "TRAV · TRAJ", v: r => r[2]}, {h: "CDR3α", v: r => r[3], cls: "mono"},
    {h: "TRBV · TRBJ", v: r => r[4]}, {h: "CDR3β", v: r => r[5], cls: "mono"}, {h: "State", v: r => r[6]}],
    [a, b].map((i, k) => [k ? "b" : "a", M.id[i], M.vA[i] + " · " + M.jA[i], M.c3A[i], M.vB[i] + " · " + M.jB[i], M.c3B[i], M.st[i] >= 0 ? ST[M.st[i]] : "–"]), {scroll: false}));
  for (const arm of Object.keys(ARMS)) {
    const E = X.arms[arm], tot = E.parts.reduce((s, x) => s + x, 0), same = E.lab[0] >= 0 && E.lab[0] === E.lab[1];
    card(app, `<h2 style="margin-top:0">${ARMS[arm].short}</h2>
     <p><b>D = ${f(E.D, 3)}</b> against the cut ${f(E.cut, 4)} → ${E.D <= E.cut ? "within the cut" : "beyond the cut"}. Atoms in the grid: ${E.atoms[0]} and ${E.atoms[1]}. Clusters: a in ${E.lab[0] >= 0 ? `<a href="#c/${arm}/${E.lab[0]}">${E.lab[0]}</a>` : "none (singleton)"}, b in ${E.lab[1] >= 0 ? `<a href="#c/${arm}/${E.lab[1]}">${E.lab[1]}</a>` : "none (singleton)"}${same ? " — the same cluster" : ""}.
     ${E.D <= E.cut && !same ? " Being within the cut does not guarantee the same cluster: complete linkage requires every pair in a cluster to be within the cut." : ""}</p>
     <h3>D² = ${f(tot, 3)}, split by chain and channel</h3>${splitBar(E.parts, tot)}${chnLegend()}
     <p class="small muted">Largest parts: ${E.parts.map((v, k) => [v, k]).sort((p, q) => q[0] - p[0]).slice(0, 3).map(([v, k]) => `${k < 7 ? "α" : "β"} ${CHN[k % 7]} ${pct(v / tot)}`).join(" · ")}. Recomputed from the atoms here and checked against the stored distance (${f(E.D_stored, 4)}).</p>`);
  }
}
