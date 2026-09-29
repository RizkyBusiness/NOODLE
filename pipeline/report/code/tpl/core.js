
// ---------- CORE (the arms of this report are declared in the data, D.arms)
const ARMS = D.arms;                                // per arm: lab, short, name, atoms ("loops" | "cdr3"), sigma (A)
let CUR = Object.keys(ARMS)[0];                     // arm shown on the Clusters and UMAP pages (no arm is preferred)
const members = {}, cinfo = {};
for (const a in ARMS) {
  members[a] = new Map(); const L = M[ARMS[a].lab];
  for (let i = 0; i < NMOL; i++) if (L[i] >= 0) { if (!members[a].has(L[i])) members[a].set(L[i], []); members[a].get(L[i]).push(i); }
  cinfo[a] = new Map();
  D.clusters[a].forEach((r, k) => cinfo[a].set(r[0], {id: r[0], size: r[1], nstate: r[2], top: r[3], frac: r[4], mice: r[5], vp: r[6], p: r[7], q: r[8], rank: k}));
}
const SPLIT_ = {};                                  // decoded on first use: b64 is defined further down (3D section)
const SPLIT = a => SPLIT_[a] || (SPLIT_[a] = b64(D.split[a], Float32Array));   // per molecule 14 parts (alpha ch1-7, beta ch1-7)
const cellsByMol = new Map();
D.cells.mi.forEach((m, i) => { if (m >= 0) { if (!cellsByMol.has(m)) cellsByMol.set(m, []); cellsByMol.get(m).push(i); } });
const SUM = a => D.summary.find(r => r.arm === a);
const CHN = ["occupancy", "hydrophobic", "aromatic", "H-bond donor", "H-bond acceptor", "positive", "negative"];
function armSelect(onchange) {
  const w = document.createElement("label");
  w.innerHTML = `Arm <select id="armsel">${Object.keys(ARMS).map(a => `<option value="${a}" ${a === CUR ? "selected" : ""}>${ARMS[a].short}</option>`).join("")}</select>`;
  w.querySelector("select").onchange = e => { CUR = e.target.value; onchange(); };
  return w;
}
function subnav(app, v, sub) {
  const L = v === "validation" ? [["#validation", "Validation"], ["#validation/explained", "Validation explained"]]
                               : [["#methods", "Methods"], ["#methods/explained", "Methods explained"]];
  sec(app, `<div class="subnav">${L.map(([h, t], k) => `<a href="${h}" class="${(k === 1) === (sub === "explained") ? "on" : ""}">${t}</a>`).join("")}</div>`);
}
function route() {
  const h = (location.hash || "#overview").slice(1); const [v, ...rest] = h.split("/");
  document.querySelectorAll("#nav a").forEach(a => a.classList.toggle("on", a.dataset.v === v || (v === "c" && a.dataset.v === "clusters")));
  const app = $("#app"); app.innerHTML = "";
  if (v === "methods") { subnav(app, v, rest[0]); if (rest[0] === "explained") viewMethodsExplained(app); else viewMethods(app); }
  else if (v === "validation") { subnav(app, v, rest[0]); if (rest[0] === "explained") viewValidationExplained(app); else viewValidation(app); }
  else if (v === "example") viewExample(app);
  else if (v === "clusters") viewClusters(app);
  else if (v === "umap") viewUmap(app);
  else if (v === "c") { if (ARMS[rest[0]]) CUR = rest[0]; viewCluster(app, rest[0], +rest[1]); }
  else viewOverview(app);
  window.scrollTo(0, 0);
}
addEventListener("hashchange", route);
const pct1 = x => x == null ? "–" : (100 * x).toFixed(1) + " %";
const sgn = x => (x >= 0 ? "+" : "−") + f(Math.abs(x), 3);
