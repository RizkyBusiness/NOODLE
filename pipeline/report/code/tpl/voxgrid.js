
// ---------- VOXEL GRIDS in the page (cluster 3D view)
// BEGIN voxgrid — pure functions; code/rv6_voxel_gate.py runs this same text in node against the Python build (gate RG4)
// Rebuilds an arm's grid from the side-car atoms exactly as pipeline/code/vxgrid.py (build): each heavy atom an isotropic
// Gaussian (sigma), point-sampled at voxel centres, truncated per axis at |c - x| <= trunc * sigma, each 1-D factor
// renormalised to sum 1 before the box clips it; channel masses from vxgrid.type_atoms (V.W, per residue x atom name).
// Per arm (V.arms[arm]): sigma, and atoms = "loops" (every loop heavy atom: the side-car holds exactly those) or "cdr3"
// (IMGT 105-117 only). Chain A -> box A, B -> box B.
// Layout per box: Float64Array(7 * nx * ny * nz), channel outer, then x, y, z (C order, as the Python grids).
function voxNew(V, arm) {
  const g = {};
  for (const k of ["A", "B"]) { const s = V.boxes[arm][k].shape; g[k] = new Float64Array(7 * s[0] * s[1] * s[2]); }
  return g;
}
function voxAdd(acc, E, V, arm, wgt) {
  const s = V.arms[arm].sigma, loops = V.arms[arm].atoms === "loops", h = V.h, lim = V.trunc * s, K = Math.ceil(2 * V.trunc * s / h) + 2;
  const f = [new Float64Array(K), new Float64Array(K), new Float64Array(K)], i0 = [0, 0, 0];
  const {X, AT, R} = E, sc = E.sc; let ai = E.ao, natoms = 0;
  for (let r = E.ro; r < E.ro + E.nr; r++) {
    const ri = R[5 * r], rn = R[5 * r + 1], ch = R[5 * r + 3] ? "B" : "A", na = R[5 * r + 4];
    const use = loops || (rn >= V.cdr3[0] && rn <= V.cdr3[1]);
    const box = V.boxes[arm][ch], lo = box.lo, sh = box.shape, G = acc[ch], nv = sh[0] * sh[1] * sh[2];
    for (let j = 0; j < na; j++, ai++) {
      if (!use) continue;
      const w = V.W[ri][AT[ai]]; if (!w.some(v => v)) continue;
      natoms++;
      for (let d = 0; d < 3; d++) {
        const u = (X[3 * ai + d] / sc - lo[d]) / h - 0.5, a0 = Math.floor(u - lim / h); let tot = 0;
        for (let q = 0; q < K; q++) { const off = (a0 + q - u) * h; const v = Math.abs(off) > lim ? 0 : Math.exp(-0.5 * (off / s) * (off / s)); f[d][q] = v; tot += v; }
        for (let q = 0; q < K; q++) { const ii = a0 + q; f[d][q] = ii < 0 || ii >= sh[d] ? 0 : f[d][q] / tot; }
        i0[d] = a0;
      }
      for (let c = 0; c < 7; c++) {
        const wc = wgt * w[c]; if (!wc) continue; const base = c * nv;
        for (let p = 0; p < K; p++) { const fx = f[0][p]; if (!fx) continue; const bx = (i0[0] + p) * sh[1];
          for (let q = 0; q < K; q++) { const fy = f[1][q]; if (!fy) continue; const by = (bx + i0[1] + q) * sh[2] + i0[2], fxy = wc * fx * fy;
            for (let t = 0; t < K; t++) { const fz = f[2][t]; if (fz) G[base + by + t] += fxy * fz; } } }
      }
    }
  }
  return natoms;
}
// squared grid distance over both boxes and channels 1-7 (sum of squared differences / h^3, as vxgrid.grid_d2)
function voxD2(g1, g2, V) {
  let s = 0;
  for (const k of ["A", "B"]) { const a = g1[k], b = g2[k]; for (let i = 0; i < a.length; i++) { const d = a[i] - b[i]; s += d * d; } }
  return s / (V.h * V.h * V.h);
}
// END voxgrid
// one channel of one box as a 3Dmol volume in display coordinates: intrinsic (x, y, z) is drawn at (x, z, -y) (P3)
function voxVolume(g, box, c, h) {
  const [nx, ny, nz] = box.shape, nv = nx * ny * nz, lo = box.lo, out = new Float32Array(nv);
  let mx = 0;
  for (let x = 0; x < nx; x++) for (let y = 0; y < ny; y++) for (let z = 0; z < nz; z++) {
    const v = g[c * nv + (x * ny + y) * nz + z]; out[(x * nz + z) * ny + (ny - 1 - y)] = v; if (v > mx) mx = v; }
  return {vol: {size: {x: nx, y: nz, z: ny}, unit: {x: h, y: h, z: h}, origin: {x: lo[0] + 0.5 * h, y: lo[2] + 0.5 * h, z: -(lo[1] + (ny - 0.5) * h)}, data: out, matrix: null}, max: mx};
}
