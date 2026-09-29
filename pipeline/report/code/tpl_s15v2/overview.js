
// ---------- OVERVIEW (s15 v2: plain-language overview; the sigma comparison table is on the Validation page)
function survTable(a, rows) {
  const cols = [
    {h: "Group", v: r => r.cluster, f: r => `<a href="#c/${a}/${r.cluster}">${r.cluster}</a>`, num: 1},
    {h: "Receptors", v: r => r.size, num: 1},
    {h: "Most common cell type", v: r => r.top, f: r => `<span class="dot" style="background:var(${stateVar(r.top)})"></span>${r.top} ${pct(r.frac)}`},
    {h: "Also a sequence group?", v: r => TROW(a, r.cluster) ? TROW(a, r.cluster).r[3] : null, f: r => TROW(a, r.cluster) ? tBadge(TROW(a, r.cluster).r[3]) : "–"},
    {h: "Mice", v: r => r.mice, num: 1}, {h: "q", v: r => r.q, f: r => fq(r.q), num: 1}];
  return table(cols, rows, {scroll: false});
}
const DIF2 = (atoms, nch) => D.s15.diff.find(r => r.atoms === atoms && r.channels === nch);
const TS = a => D.tcrd.summ.find(r => r.arm === a);
const strength = (z, ref) => z >= ref ? "about as strong as" : z >= 0.75 * ref ? "a little weaker than" : z >= 0.4 * ref ? "clearly weaker than" : "much weaker than";
function armCard(a) {
  const P = SUM(a), T = TS(a), what = ARMS[a].atoms === "loops" ? "all four loops" : "the two CDR3 loops only", sharp = ARMS[a].sigma === 1.5 ? "sharper picture" : "standard picture";
  return `<div class="card"><h3 style="margin-top:0">${ARMS[a].short}</h3><p class="small muted">Shape of ${what}, ${sharp} (σ ${f(ARMS[a].sigma, 1)} Å).</p><div class="grid g4">
   <div class="card tile"><div class="lab">Receptors placed in a group</div><div class="val">${P.clustered.toLocaleString()}</div><div class="sub">of ${NMOL.toLocaleString()} · ${P.clusters_ge3.toLocaleString()} groups of 3 or more</div></div>
   <div class="card tile"><div class="lab">Groups linked to one cell type</div><div class="val">${P.survivors}</div><div class="sub">statistically, after correcting for testing many groups</div></div>
   <div class="card tile"><div class="lab">… of which also sequence groups</div><div class="val">${T.surv_also} of ${T.survivors}</div><div class="sub">receptors that also look alike by sequence</div></div></div></div>`;
}
function viewOverview(app) {
  const R = D.ref30, dL = DIF2("loops", 7), dC = DIF2("cdr3", 7), Tr = D.tcrd.ref, Dd = SUM("D"), Ee = SUM("E"), Bb = SUM("B"), Cc = SUM("C");
  const noGain = r => r.dE1_lo <= 0 && r.dE1_hi >= 0;
  sec(app, `<h1>Do T cells with similarly shaped receptors do similar jobs?</h1>
  <p class="muted">${esc(D.dataset_label)} · ${NMOL.toLocaleString()} different T-cell receptors, each modelled in 3D by computer · ${esc(D.brand_note)}.</p>
  <div class="card"><h2 style="margin-top:0">What this is about</h2>
   <p>Every T cell carries a receptor that recognises a target. The tip of the receptor is made of six flexible loops; the two in the middle (the <b>CDR3 loops</b>) are created fresh in each cell and make most of the contact with the target. From single-cell sequencing we know, for each cell, both its receptor and what kind of T cell it is — for example a regulatory T cell (Treg), a Th1 or Th17 helper cell, or a follicular helper cell (Tfh).</p>
   <p>The question: <b>if two receptors have a similar 3D shape, are the cells that carry them more often the same kind of T cell?</b> If so, receptor shape would carry information about what a T cell does.</p></div>
  <div class="card"><h2 style="margin-top:0">How the shapes were compared</h2><ul class="plain">
   <li>Each receptor's 3D model was turned into a picture made of tiny 1 Å cubes (a "voxel grid"), recording where the loop atoms are and what chemistry they carry (water-repelling, ring-shaped, charged, able to form hydrogen bonds).</li>
   <li>Two receptors are compared cube by cube. Receptors whose pictures are as similar as the closest 1 in 100 random pairs are put in the same <b>group</b>, and every pair inside a group must pass that test.</li>
   <li>Four versions are shown: pictures of <b>all four loops</b> or of the <b>CDR3 loops only</b>, each drawn with a <b>standard blur</b> (σ 2.0 Å, Arms B and C) or a <b>sharper</b> one (σ 1.5 Å, Arms D and E, this report's addition).</li>
   <li>For every group we then asked how often its cells share one cell type, compared with chance, and whether the same receptors would also be grouped by <b>sequence similarity</b> alone (the established TCRdist method).</li></ul></div>`);
  sec(app, `<div class="grid g2">${["D", "E", "B", "C"].map(armCard).join("")}</div>`);
  card(app, `<h2 style="margin-top:0">What we found, in plain words</h2><ul class="plain">
   <li><b>Shape groups do lean towards one cell type — but only a little.</b> Receptors in the same shape group are somewhat more often the same kind of T cell than chance would give, also within the same mouse. Much of that link, however, comes from the gene segments (V genes) the receptors are built from: comparing only receptors with the same V genes, the link is ${strength(Math.max(Dd.z_mouse_V, Ee.z_mouse_V), R.z_mouse_V)} that of the existing vector-based method (${R.survivors} linked groups there).</li>
   <li><b>A sharper picture did not reveal more real detail.</b> Checked against ${D.validation.panel.receptors} receptors whose true 3D structure was solved in the laboratory, the sharper grids matched the real structures ${noGain(dL) && noGain(dC) ? "no better than" : "differently from"} the standard ones, and the computer models looked less like their real structures (more noise). So the extra sharpness mostly adds prediction error; for computer-modelled receptors the standard blur is preferable.</li>
   <li><b>CDR3-only groups are mostly receptors with similar sequences.</b> Of the linked groups found with the CDR3 loops, ${TS("E").surv_also} of ${TS("E").survivors} (sharp) and ${TS("C").surv_also} of ${TS("C").survivors} (standard) are also groups by sequence similarity. They largely confirm what sequence already shows.</li>
   <li><b>All-loop groups often are not.</b> Only ${TS("D").surv_also} of ${TS("D").survivors} (sharp) and ${TS("B").surv_also} of ${TS("B").survivors} (standard) of the linked groups found with all four loops are also sequence groups. These are the candidates for "similar shape without similar sequence" — but the all-loop pictures are also the ones most influenced by V genes (${pct(D.s15.diag.find(r => r.arm.startsWith("Arm D")).share_single_Vpair_clusters)} and ${pct(D.s15.diag.find(r => r.arm.startsWith("Arm B")).share_single_Vpair_clusters)} of their groups use a single V-gene pair), so they need care.</li>
   <li><b>Overall:</b> compare each version's link to cell type with the reference method's in the table above. The voxel grid does capture real structure checked against laboratory structures — including structure the existing method misses — but that extra structure does not translate into a stronger link to what the cells do.</li></ul></div>`);
  card(app, `<h2 style="margin-top:0">How to read this report</h2><ul class="plain">
   <li><a href="#umap"><b>UMAP</b></a> — a map of all cells, where cells with similar gene activity sit close together (coloured by cell type). Choose a version and a group to see where its cells lie, or colour the map by any gene.</li>
   <li><a href="#clusters"><b>Clusters</b></a> — the list of shape groups; filter by size, cell type, significance or whether they are also sequence groups. Each group's page shows:
     <ul><li>its cells on the map and its cell-type make-up;</li>
     <li>the receptors' loops in 3D with their side chains, and the superimposed voxel picture (the average of the members' pictures) — rotate, zoom, pick a chemistry or a single receptor;</li>
     <li>the sequence cross-check: whether the same receptors form a sequence group, with a heat map of how different each pair is by sequence;</li>
     <li>which part of each loop's picture (chain and chemistry) makes members differ, and their aligned sequences.</li></ul></li>
   <li><a href="#methods"><b>Methods</b></a> and <a href="#validation"><b>Validation</b></a> — what was done and the checks against laboratory structures, each with a plain-language version (the detailed sharp-versus-standard checks are on the Validation page).</li>
   <li><a href="#example"><b>Worked example</b></a> — two receptors compared step by step in every version.</li></ul></div>`);
  card(app, `<h2 style="margin-top:0">Words used</h2><ul class="plain">
   <li><b>Receptor / TCR</b> — the T-cell receptor, made of an α and a β chain. <b>CDR1, CDR2, HV4, CDR3</b> — its loops; CDR1, CDR2 and HV4 are fixed by the V gene, CDR3 is made anew in each cell.</li>
   <li><b>V gene</b> — one of the gene segments a receptor chain is built from; receptors with the same V genes share most of their loops.</li>
   <li><b>Voxel</b> — a 1 Å cube of the 3D picture. <b>σ (blur)</b> — how far each atom is smeared: 2.0 Å standard, 1.5 Å sharper.</li>
   <li><b>Group (cluster)</b> — receptors whose pictures are all within the similarity threshold of each other.</li>
   <li><b>Linked / significant</b> — a group whose cells are one cell type more often than chance allows, after correcting for testing many groups (q ≤ 0.15, Benjamini–Hochberg).</li>
   <li><b>Sequence group (TCRdist)</b> — receptors grouped by how similar their loop sequences are, with no 3D information (Dash et al. 2017; tcrdist3).</li>
   <li><b>Crystal structure</b> — a receptor's 3D shape determined in the laboratory, used to check the computer models.</li></ul></div>`);
  card(app, `<h2 style="margin-top:0">Cautions</h2><ul class="plain">
   <li>All shapes are computer predictions; small errors in them become noise in the comparison.</li>
   <li>The links to cell type are statistical and modest; a linked group is a lead to follow up, not a demonstration of function.</li>
   <li>The four versions are several looks at the same cells; each version's significance is corrected within that version only, and no version was chosen as "best" on these results.</li></ul></div>`);
  for (const a of ["D", "E"]) {
    sec(app, `<h2>Linked groups · ${ARMS[a].short}</h2><p class="small muted">Groups whose cells are one cell type more often than chance (q ≤ 0.15), strongest first. Click a group to open its page.</p>`);
    app.append(survTable(a, D.surv[a]));
  }
}
