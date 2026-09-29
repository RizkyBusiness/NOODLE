
// ---------- tcrdist3 cross-check (A15; report only): is a voxel cluster also a group by sequence (TCRdist)?
const REF_ARM = (D.tcrd && D.tcrd.ref_name) || "reference method";
const TV = ["also", "partly", "not"], TVL = ["also a tcrdist3 cluster", "partly (50–75 %)", "not a tcrdist3 cluster"];
const TVC = ["var(--good)", "var(--warn)", "var(--vother)"];
const TROW_ = {};
const TROW = (a, id) => { if (!TROW_[a]) { TROW_[a] = new Map(D.tcrd.arms[a].rows.map((r, k) => [r[0], {r, k}])); } return TROW_[a].get(id); };
const TPAIR_ = {};
const TPAIRS = a => TPAIR_[a] || (TPAIR_[a] = b64(D.tcrd.arms[a].d, Int16Array));
const tBadge = v => v == null ? "–" : `<span style="color:${TVC[v]};font-weight:600">${TV[v]}</span>`;
function tcrdPanel(app, arm, id, mem) {
  const T = TROW(arm, id), cut = D.tcrd.cut;
  if (!T) { card(app, `<h3 style="margin-top:0">tcrdist3 cross-check</h3><p class="small muted">Only clusters of 2 or more molecules are cross-checked.</p>`); return; }
  const [c, main, fr, v, ncl, nsing, med, mx, share, pctl, top] = T.r;
  const k = mem.length, off = D.tcrd.arms[arm].off[T.k], P = TPAIRS(arm).subarray(off, off + k * k);
  // members in the stored order (ascending molecule index, as mem); sort the heat map by tcrdist3 cluster, then clone id
  const ord = [...Array(k).keys()].sort((p, q) => (M.T[mem[p]] < 0) - (M.T[mem[q]] < 0) || M.T[mem[p]] - M.T[mem[q]] || (M.id[mem[p]] < M.id[mem[q]] ? -1 : 1));
  const HM = Math.min(360, Math.max(150, 12 * k));
  const tops = top ? top.split(";").map(s => s.split(":").map(Number)) : [];
  const d = card(app, `<h3 style="margin-top:0">tcrdist3 cross-check <span class="small muted">(sequence similarity, A15 · report only)</span></h3>
   <p><b style="color:${TVC[v]}">${TVL[v]}</b> — ${pct(fr)} of the ${k} members fall in one tcrdist3 cluster${main >= 0 ? ` (${main})` : ""}; the members span ${ncl} tcrdist3 cluster${ncl === 1 ? "" : "s"}${nsing ? ` and ${nsing} tcrdist3 singleton${nsing === 1 ? "" : "s"}` : ""}${tops.length ? ` (largest: ${tops.map(([a, b]) => `${a} × ${b}`).join(", ")})` : ""}.</p>
   <p class="small">Pairwise TCRdist among members: median ${f(med, 0)} (${pctl < 0.1 ? "fewer than 0.1 %" : f(pctl, 1) + " %"} of random receptor pairs are this close), maximum ${mx}; ${pct(share)} of member pairs within the tcrdist3 cut (${f(cut, 0)}).</p>
   <div style="display:flex;gap:14px;flex-wrap:wrap;align-items:flex-start"><canvas id="tcrdhm" width="${HM}" height="${HM}" style="width:${HM}px;height:${HM}px;max-width:100%;flex:none;border:1px solid var(--line);border-radius:6px"></canvas>
    <div class="small muted" style="max-width:330px">Member × member TCRdist (αβ), members ordered by their tcrdist3 cluster (bars on the edges: one shade per tcrdist3 cluster, none = singleton). Colour: white 0 → dark at the random-pair median (${f(D.tcrd.summary.bg_median, 0)}); cells within the tcrdist3 cut are outlined. Hover for the pair.<div id="tcrdtip" class="mono" style="margin-top:6px;min-height:2.4em"></div></div></div>
   <p class="small muted">"Also a tcrdist3 cluster" = at least 75 % of the members in one cluster of an independent tcrdist3 partition built with the same rule as the voxel clusters (1st percentile of the same 60,000 random pairs, complete linkage). TCRdist compares CDR1, CDR2, CDR2.5 and CDR3 sequences (Dash et al. 2017; tcrdist3, Mayer-Blackwell et al. 2021) — no structure.</p>`);
  const cv = d.querySelector("#tcrdhm"), g = cv.getContext("2d"), W = cv.width, pad = 6, cs = (W - pad) / k, top_ = css("--ink"), bgm = D.tcrd.summary.bg_median;
  const tl = [...new Set(ord.map(p => M.T[mem[p]]).filter(x => x >= 0))], shade = x => x < 0 ? null : `hsl(${(tl.indexOf(x) * 67) % 360},55%,55%)`;
  g.fillStyle = css("--surface"); g.fillRect(0, 0, W, W);
  ord.forEach((p, i) => { const s = shade(M.T[mem[p]]); if (s) { g.fillStyle = s; g.fillRect(pad + i * cs, 0, cs, pad - 1); g.fillRect(0, pad + i * cs, pad - 1, cs); } });
  ord.forEach((p, i) => ord.forEach((q, j) => { const x = P[p * k + q], t = Math.min(1, x / bgm), l = Math.round(97 - 72 * t);
    g.fillStyle = `hsl(215,45%,${l}%)`; g.fillRect(pad + j * cs, pad + i * cs, cs, cs);
    if (i !== j && x <= cut && cs >= 5) { g.strokeStyle = top_; g.lineWidth = 0.6; g.strokeRect(pad + j * cs + 0.5, pad + i * cs + 0.5, cs - 1, cs - 1); } }));
  cv.onmousemove = e => { const r = cv.getBoundingClientRect(), i = Math.floor(((e.clientY - r.top) * W / r.height - pad) / cs), j = Math.floor(((e.clientX - r.left) * W / r.width - pad) / cs);
    if (i < 0 || j < 0 || i >= k || j >= k) return; const p = ord[i], q = ord[j];
    d.querySelector("#tcrdtip").textContent = `${M.id[mem[p]]} × ${M.id[mem[q]]}: TCRdist ${P[p * k + q]} · tcrdist3 clusters ${M.T[mem[p]]} / ${M.T[mem[q]]}`; };
}
function tcrdOverview(app) {
  const S = D.tcrd.summ.concat([Object.assign({arm: REF_ARM}, D.tcrd.ref)]), cut = D.tcrd.cut, T1 = D.tcrd.summary;
  const n3 = r => `${r.ge3_also} / ${r.ge3_partly} / ${r.ge3_not}`, ns = r => `<b>${r.surv_also}</b> / ${r.surv_partly} / ${r.surv_not}`;
  card(app, `<h2 style="margin-top:0">Are the clusters also sequence groups? tcrdist3 cross-check</h2>
   <p class="small muted">Report only (A15). An independent partition of the same ${NMOL.toLocaleString()} molecules by paired αβ TCRdist (tcrdist3; CDR1, CDR2, CDR2.5 and CDR3 sequences, CDR3 weighted 3×; no structure), built with the voxel rule: cut = 1st percentile of the same 60,000 random pairs (TCRdist ${f(cut, 0)}; random-pair median ${f(T1.bg_median, 0)}), complete linkage → ${T1.clusters.toLocaleString()} clusters, ${T1.clustered.toLocaleString()} molecules clustered. A cluster is <b>also</b> a tcrdist3 cluster when ≥ 75 % of its members fall in one tcrdist3 cluster, <b>partly</b> at 50–75 %, <b>not</b> below. the reference method (vector method) for reference.</p>`)
   .append(table([{h: "Arm", v: r => r.arm, f: r => r.arm === "the reference method" ? "reference method" : esc(ARMS[r.arm].short)},
    {h: "Clusters ≥ 3: also / partly / not", v: r => r.ge3_also / r.clusters_ge3, f: n3, num: 1},
    {h: "Share also", v: r => r.ge3_also / r.clusters_ge3, f: r => pct(r.ge3_also / r.clusters_ge3), num: 1},
    {h: "Significant clusters: also / partly / not", v: r => r.surv_also / Math.max(1, r.survivors), f: ns, num: 1},
    {h: "Share also", v: r => r.surv_also / Math.max(1, r.survivors), f: r => r.survivors ? pct(r.surv_also / r.survivors) : "–", num: 1}], S, {scroll: false}));
}
const TCRD_REFS = "Dash P et al. Quantifiable predictive features define epitope-specific T cell receptor repertoires. Nature 2017;547:89–93 (PMID 28636592). Mayer-Blackwell K et al. TCR meta-clonotypes for biomarker discovery with tcrdist3. eLife 2021;10:e68605 (PMID 34845983).";
function tcrdMethods(app) {
  sec(app, `<h2>tcrdist3 cross-check (A15, report only)</h2>
  <p>Paired αβ TCRdist from tcrdist3 0.3 (defaults: substitution distance 0–4 derived from BLOSUM62 over CDR1, CDR2 and CDR2.5 from the V gene, weight 1; CDR3 weight 3, trimmed 3 residues N-terminal and 2 C-terminal, best gap position; 4 per gap; α + β). In tcrdist3 0.3 the IMGT gap symbol in the germline CDR strings scores 0 against any residue. Input: V gene (as the *01 allele of the tcrdist3 mouse database; all genes present) and CDR3 of each receptor. Partition with the voxel rule (1st percentile of the same 60,000 random pairs, complete linkage), and per voxel cluster the share of members in one tcrdist3 cluster (≥ 75 % = also, the reference method's match rule). Checks: 50 random pairs equal an independent implementation of the formula exactly; the control pairs one CDR3 residue apart lie far below the cut (their median below the random pairs' 1st percentile, ${f(D.tcrd.cut, 0)}); every tcrdist3 cluster within its cut; 3 clusters per arm recomputed by a separate route. No state test was run on the tcrdist3 partition.</p>
  <p class="small muted">${TCRD_REFS}</p>`);
}
