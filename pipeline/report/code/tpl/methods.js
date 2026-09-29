
// ---------- METHODS (technical) and METHODS EXPLAINED (plain language)
function viewMethods(app) {
  const B = SUM("B"), C = SUM("C"), Fr = D.frames, tb = D.thr, bx = D.boxes;
  sec(app, `<h1>Methods</h1><p class="muted small">NOODLE · method details in the package: references/method.md</p>
  <div class="warnbox"><b>Two arms are reported throughout:</b> Arm B (all four loops) and Arm C (CDR3 only). Both use the same frame, grid, chemistry, distance, threshold rule, clustering and state test; they differ only in which atoms enter the grid.</div>
  <h2>1 · Data and structures</h2><ul>
  <li>Input: paired \u03b1\u03b2 V(D)J sequences with single-cell transcriptomes (${esc(D.dataset_label)}).</li>
  <li>Each receptor folded with TCRBuilder2+ (Quast et al. 2025), ${D.counts.n_models.toLocaleString("en-US")} models; identical V-domain protein sequences are one molecule: ${NMOL.toLocaleString()} molecules, ${D.n_state.toLocaleString()} with a transcriptional state.</li>
  <li>IMGT numbering (Lefranc et al. 2003). Loops: CDR1 27–38, CDR2 56–65, HV4 81–86, CDR3 105–117, insertions included. Heavy atoms only.</li></ul>
  <h2>2 · Frame</h2>
  <p>Each chain is aligned on its own five framework landmarks (Cα of IMGT 23, 41, 89, 104, 118) onto the same chain of a fixed reference receptor by the Kabsch method (Kabsch 1976; row convention), then one fixed rotation into the TCR-intrinsic axes (z toward the pMHC). Median per-chain fits: α ${f(Fr.fa, 2)} Å, β ${f(Fr.fb, 2)} Å. The two chains are placed separately, so the Vα/Vβ angle does not enter the comparison. (A whole-molecule ten-landmark frame was tested first: fit ${f(Fr.lm10, 2)} Å, framework ${f(Fr.anch, 2)} Å; see Validation.)</p>
  <h2>3 · Atom sets and boxes</h2><ul>
  <li><b>Arm B:</b> heavy atoms of CDR1, CDR2, HV4 and CDR3 of both chains (about a quarter of all heavy atoms). Per-chain boxes ${bx.B.A.join(" × ")} and ${bx.B.B.join(" × ")} voxels (7 Å margin).</li>
  <li><b>Arm C:</b> heavy atoms of CDR3 α and β only. Boxes ${bx.C.A.join(" × ")} and ${bx.C.B.join(" × ")} voxels: the full envelope of every CDR3 atom of the repertoire and the benchmark + 7 Å, so no atom is truncated anywhere.</li></ul>
  <h2>4 · Chemistry channels</h2>`);
  app.append(table([{h: "#", v: r => r[0], num: 1}, {h: "Channel", v: r => r[1]}, {h: "Atoms", v: r => r[2]}, {h: "Mass per atom", v: r => r[3]}], [
    [1, "occupancy", "every heavy atom", "1"], [2, "hydrophobic", "carbon not bonded to N or O", "1"], [3, "aromatic", "ring atoms of Phe, Tyr, Trp, His", "1"],
    [4, "H-bond donor", "backbone N except Pro; Arg NE/NH1/NH2, Asn ND2, Gln NE2, Lys NZ, Trp NE1, His ND1/NE2, Ser OG, Thr OG1, Tyr OH", "1"],
    [5, "H-bond acceptor", "all O; His ND1/NE2", "1"], [6, "positive", "Lys NZ; Arg CZ/NE/NH1/NH2; His ring N", "1 / 0.25 / 0.05 (per formal charge)"],
    [7, "negative", "Asp OD1/OD2, Glu OE1/OE2; C-terminal OXT", "0.5 / 1"]], {scroll: false}));
  sec(app, `<p class="small muted">Channels follow the voxelised pharmacophore approach of DeepSite (Jiménez et al. 2017).</p>
  <h2>5 · Grid</h2>
  <p>Each atom is an isotropic Gaussian of width σ = 2.0 Å, point-sampled at the centres of 1.0 Å voxels, truncated at 3σ per axis and renormalised so every atom deposits exactly its mass. Every build was checked against the exact grid-free distance between the Gaussian-smeared atom sets (gridded values within +1.2 to +1.8 %, the effect of the 3σ truncation), against a slow reference loop (≤ 1e-15), and for mass conservation (< 0.5 %). σ = 2.0 Å came from a pre-registered blur sweep on the crystal benchmark.</p>
  <h2>6 · Distance</h2>
  <p>D(a, b)² = Σ over the α and β boxes and the 7 channels of |ρ<sub>a</sub> − ρ<sub>b</sub>|² / h³ (h = 1 Å), computed by a blocked Gram matrix in double precision. Because it is a plain sum, D² splits <b>exactly</b> into 14 non-negative parts (2 chains × 7 channels); the cluster pages show that split, and for every within-cluster pair the parts add up to the stored D² to 1.2 × 10⁻⁷.</p>
  <h2>7 · Threshold and clustering</h2>
  <p>As in the vector method: the cut is the 1st percentile of D over the same 60,000 random receptor pairs (seed 0) — Arm B ${f(tb.B.bg_p1, 4)}, Arm C ${f(tb.C.bg_p1, 4)}; complete linkage on the distinct molecules at the cut; singletons dropped. Arm B: ${B.clusters.toLocaleString()} clusters, ${B.clustered.toLocaleString()} molecules. Arm C: ${C.clusters.toLocaleString()} clusters, ${C.clustered.toLocaleString()} molecules. A ±2 % change of the cut gives ARI ${f(D.stab.B[0], 2)} / ${f(D.stab.B[1], 2)} (B) and ${f(D.stab.C[0], 2)} / ${f(D.stab.C[1], 2)} (C) (Hubert & Arabie 1985).</p>
  <h2>8 · State tests</h2>
  <p>Exactly the machinery of the reference method, reproduced on the reference method's own clusters before these arms were tested: purity over clusters with ≥ 3 molecules; 200 permutations of states, re-seeded per arm, free, within mouse, within mouse and V pair, and — because the grid keeps CDR3 length — within mouse, V pair and CDR3 length class; per-cluster one-sided binomial test against the repertoire base rate; Benjamini–Hochberg at q ≤ 0.15 (Benjamini & Hochberg 1995). Each arm was tested once.</p>`);
  app.append(table([{h: "Null", v: r => r[0]}, {h: "Arm B excess (z)", v: r => r[1]}, {h: "Arm C excess (z)", v: r => r[2]}],
    [["free", "unstratified"], ["within mouse", "mouse"], ["within mouse + V pair", "mouse_V"], ["within mouse + V pair + CDR3 length", "mouse_V_len"]].map(([n, k]) =>
      [n, `${sgn(B["delta_" + k])} (${f(B["z_" + k], 1)})`, `${sgn(C["delta_" + k])} (${f(C["z_" + k], 1)})`]), {scroll: false}));
  sec(app, `<h2>9 · How the design evolved</h2><ul class="plain">
  <li>The first design (whole molecule, one ten-landmark frame) was tested against a pre-registered crystal-error limit; see the Validation page for this run\u2019s value.</li>
  <li>A validity panel then compares the grid and the reference descriptor identically on the crystal-benchmarked receptors (Validation page). Reinterpreting a failed pre-registered test is recorded as a change of criterion (amendment), never as a pass.</li>
  <li>A per-chain Vα/Vβ hinge descriptor failed its reliability gate and is exploratory only. Arms B and C were then fixed before any state data were read.</li></ul>
  <h2>10 · References</h2><ul class="small">
  <li>Kabsch W. Acta Crystallogr A 1976;32:922–3.</li><li>Lefranc M-P et al. Dev Comp Immunol 2003;27:55–77.</li>
  <li>Jiménez J et al. DeepSite. Bioinformatics 2017;33:3036–42.</li><li>Quast NP et al. TCRBuilder2+. Commun Biol 2025;8:362.</li>
  <li>Dunbar J et al. TRangle. PLoS Comput Biol 2014;10:e1003852.</li><li>Koo TK, Li MY. J Chiropr Med 2016;15:155–63.</li>
  <li>Mantel N. Cancer Res 1967;27:209–20.</li><li>Henikoff S, Henikoff JG. PNAS 1992;89:10915–9.</li>
  <li>Benjamini Y, Hochberg Y. J R Stat Soc B 1995;57:289–300.</li><li>Hubert L, Arabie P. J Classif 1985;2:193–218.</li><li>Rego N, Koes D. 3Dmol.js. Bioinformatics 2015;31:1322–4.</li></ul>`);
  tcrdMethods(app);
}
function voxelSketch() {   // generic schematic drawn here (not data): atoms smeared into a grid of cubes
  const cells = []; for (let i = 0; i < 10; i++) for (let j = 0; j < 6; j++) {
    const d = Math.hypot(i - 3.2, j - 2.6), e = Math.hypot(i - 6.6, j - 3.1), v = Math.max(Math.exp(-d * d / 3), 0.8 * Math.exp(-e * e / 2.4));
    cells.push(`<rect x="${40 + 34 * i}" y="${40 + 34 * j}" width="32" height="32" rx="3" fill="var(--accent)" fill-opacity="${(0.06 + 0.8 * v).toFixed(2)}" stroke="var(--line)"/>`); }
  return `<svg viewBox="0 0 420 270" style="width:100%;max-width:420px;height:auto;display:block" role="img" aria-label="Schematic: two atoms smeared into a grid of cubes">
   ${cells.join("")}<circle cx="${40 + 34 * 3.2 + 16}" cy="${40 + 34 * 2.6 + 16}" r="6" fill="var(--ink)"/><circle cx="${40 + 34 * 6.6 + 16}" cy="${40 + 34 * 3.1 + 16}" r="6" fill="var(--ink)"/>
   <text x="210" y="262" text-anchor="middle" font-size="12" fill="var(--ink-2)">each atom (dot) spreads into the nearby cubes; darker = more atom</text></svg>`;
}
function viewMethodsExplained(app) {
  const B = SUM("B"), C = SUM("C");
  sec(app, `<h1>Methods explained</h1><p class="muted">The same method as the technical page, in plain language.</p>
  <div class="card"><h2 style="margin-top:0">A 3D picture made of cubes</h2>
   <p>Every receptor's loops — the parts that touch the target — are turned into a 3D picture made of small cubes, 1 Å on a side. Each atom is blurred into the cubes around it, and each cube records how much atom it holds and what kind of chemistry: water-repelling, aromatic, able to form hydrogen bonds, positively or negatively charged.</p>${voxelSketch()}</div>
  <div class="card"><h2 style="margin-top:0">Two versions</h2><ul class="plain">
   <li><b>Arm B</b> pictures all four loops. CDR1, CDR2 and HV4 are written entirely by the V gene, so receptors sharing V genes look alike here.</li>
   <li><b>Arm C</b> pictures only the two CDR3 loops, the part created when the receptor is made — the same part the existing method describes.</li></ul></div>
  <div class="card"><h2 style="margin-top:0">Lining up and comparing</h2>
   <p>Each half of the receptor is lined up separately on five fixed points of its stable core, so every loop sits in the same place relative to its own chain. Two receptors are then compared cube by cube: the differences are squared and added up. That sum splits exactly into 14 parts — which chain (α or β) and which kind of chemistry — and the cluster pages show those parts.</p></div>
  <div class="card"><h2 style="margin-top:0">Groups and the test against cell behaviour</h2>
   <p>Two receptors count as "the same shape" when they are as close as the closest 1 in 100 random pairs. Receptors closer than that are grouped so that every pair in a group is within the threshold (Arm B: ${B.clusters.toLocaleString()} groups; Arm C: ${C.clusters.toLocaleString()}). We then ask whether cells in the same group behave alike, comparing against shuffles that keep the mouse, the V genes, and finally also the CDR3 length fixed. The last two shuffles matter most: they ask whether shape adds anything beyond the genes.</p></div>
  <div class="card"><h2 style="margin-top:0">What to look at</h2><ul class="plain">
   <li>Whether the cube pictures capture real structure is checked against lab-solved crystal structures on the <a href="#validation">Validation</a> page.</li>
   <li>How much of any link to cell behaviour survives holding the V genes (and CDR3 length) fixed is the null ladder in section 7 above.</li>
   <li>Whether the CDR3-only arm avoids grouping by V gene is the single-V-pair share in the <a href="#overview">Overview</a>.</li></ul></div>
  <div class="card"><h2 style="margin-top:0">For specialists</h2><dl class="gl">
   <dt>cube picture</dt><dd>voxel grid, σ 2.0 Å Gaussians on 1.0 Å voxels, 7 channels</dd>
   <dt>lined up separately</dt><dd>per-chain Kabsch fit on IMGT 23, 41, 89, 104, 118 onto the reference receptor (per-chain frame)</dd>
   <dt>closest 1 in 100</dt><dd>1st percentile of 60,000 random pairs (seed 0)</dd>
   <dt>every pair within</dt><dd>complete linkage</dd>
   <dt>shuffles</dt><dd>200-permutation null ladder; binomial per cluster; Benjamini–Hochberg q ≤ 0.15</dd></dl></div>`);
}
