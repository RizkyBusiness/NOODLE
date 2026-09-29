"""the voxel pipeline report, stage 6: gate RG4 — the voxel grids drawn on the cluster pages equal the arms' grids.

The page rebuilds each member's grid from the side-car atoms (report_data/atoms_NNN.js, unrounded float32) with the voxgrid functions
and draws the cluster mean. Here the voxgrid code and the side-car decoder are cut out of the ASSEMBLED voxel_report.html
and run in node on the same side-car files, for every arm of the report (main: B, C; s15: D, E, B, C), 3 clusters per arm (largest survivor, one of about 8 members, one pair) and
the worked-example pair; the results are compared with the Python build (vxv_common.build_grid, float64) and with the
stored distance matrices (built from float16-stored grids):
  RG4a  atoms per member == build_grid's atom count
  RG4b  mean grid of the cluster: max |page - Python| <= 1e-3 x its peak
  RG4c  every within-cluster pair: |D_page - D_Python| / D_Python <= 1e-3
  RG4d  every within-cluster pair: |D_page - D_stored| / D_stored <= 5e-3
Writes checks/RG4_voxels[_s15].csv (temporary grid files in work/ are deleted). Paths come from config/voxel_config.json.:
python pipeline/report/code/rv6_voxel_gate.py [main|s15].
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", "..", "code"))
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, subprocess
import numpy as np, pandas as pd
sys.path.insert(0, VXP("voxel_out/code"))
from vxv_common import boxes_for, rep_place, build_grid, H
from vxc_common import load_boxes

RP = VXP("voxel_out/report/"); W = RP + "work/"; rows = []
REPORT = sys.argv[1] if len(sys.argv) > 1 else "main"; assert REPORT in ("main", "s15", "s15v2")
HTML, SUF = {"main": ("voxel_report.html", ""), "s15": ("voxel_report_s15.html", "_s15"), "s15v2": ("voxel_report_s15_v2.html", "_s15v2")}[REPORT]
def chk(name, ok, val=""):
    rows.append(dict(gate="RG4", check=name, value=val, passed=bool(ok))); print("[%s] %-78s %s" % ("PASS" if ok else "FAIL", name, val))
h = open(RP + HTML).read()
i = h.index('<script id="data" type="application/json">') + len('<script id="data" type="application/json">')
D = json.loads(h[i:h.index("</script>", i)].replace("<\\/", "</"))
vox = h[h.index("// BEGIN voxgrid"):h.index("// END voxgrid")]
dec = h[h.index("const b64 = (s, T)"):h.index("function loadAtoms(")]
L1 = np.load(VXP("landmarks/out/LM1_internal_coords.npz"), allow_pickle=True); IDS = [str(c) for c in L1["clone_id"]]; IDX = {c: k for k, c in enumerate(IDS)}
ARM = {"B": dict(atoms="loops", boxes=boxes_for("F1"), D=VXP("voxel_out/tmp/VXH5_D_receptors.f32.npy"), sigma=2.0),
       "C": dict(atoms="cdr3", boxes=load_boxes(), D=VXP("voxel_out/tmp/VXC1_D_receptors.f32.npy"), sigma=2.0),
       "D": dict(atoms="loops", boxes=boxes_for("F1"), D=VXP("voxel_out/tmp/VXS1_D_D_receptors.f32.npy"), sigma=1.5),
       "E": dict(atoms="cdr3", boxes=load_boxes(), D=VXP("voxel_out/tmp/VXS1_E_D_receptors.f32.npy"), sigma=1.5)}
assert all(D["vox"]["arms"][a] == dict(sigma=ARM[a]["sigma"], atoms=ARM[a]["atoms"]) for a in D["arms"]), "page arm settings differ"
jobs = []
for a in D["arms"]:
    lab = np.array(D["mol"][a]); rc = {r[0]: r for r in D["clusters"][a]}
    pick = [D["surv"][a][0]["cluster"], sorted(rc, key=lambda x: (abs(rc[x][1] - 8), x))[0], sorted(rc, key=lambda x: (rc[x][1] != 2, x))[0]]
    for cid in pick:
        jobs.append(dict(arm=a, name="cluster %d" % cid, mem=[int(k) for k in np.flatnonzero(lab == cid)]))
    jobs.append(dict(arm=a, name="worked example", mem=[D["mol"]["id"].index(D["example"]["a"]), D["mol"]["id"].index(D["example"]["b"])]))
for k, j in enumerate(jobs):
    j["out"] = W + "rg4_mean_%d.f64" % k
chunks = sorted({D["atoms"]["chunk"][m] for j in jobs for m in j["mem"]})
drv = W + "rg4_driver.js"
open(drv, "w").write("""globalThis.window = globalThis; const fs = require("fs");
%s
%s
const V = %s, jobs = %s;
for (const c of %s) eval(fs.readFileSync("%sreport_data/atoms_" + String(c).padStart(3, "0") + ".js", "utf8"));
const res = jobs.map(j => {
  const g = j.mem.map(i => { const x = voxNew(V, j.arm); const n = voxAdd(x, AT_BY_MOL.get(i), V, j.arm, 1); return [x, n]; });
  const mean = voxNew(V, j.arm); j.mem.forEach(i => voxAdd(mean, AT_BY_MOL.get(i), V, j.arm, 1 / j.mem.length));
  fs.writeFileSync(j.out, Buffer.concat([Buffer.from(mean.A.buffer), Buffer.from(mean.B.buffer)]));
  const d2 = []; for (let p = 0; p < g.length; p++) for (let q = p + 1; q < g.length; q++) d2.push(voxD2(g[p][0], g[q][0], V));
  return {natoms: g.map(x => x[1]), d2};
});
fs.writeFileSync("%srg4_page.json", JSON.stringify(res));
""" % (dec, vox, json.dumps(D["vox"]), json.dumps(jobs), json.dumps(chunks), RP, W))
subprocess.run(["node", drv], check=True)
PG = json.load(open(W + "rg4_page.json"))
worst = dict(b=0.0, c=0.0, d=0.0); na_bad = 0; npairs = 0
for j, r in zip(jobs, PG):
    A = ARM[j["arm"]]; Dm = np.load(A["D"], mmap_mode="r")
    rw = [IDX[D["mol"]["id"][m]] for m in j["mem"]]                    # molecule index -> the landmark file row (placements, stored matrix)
    G = [build_grid((VXP("structures/%s.pdb") % IDS[x], rep_place(x, "F1"), A["boxes"], A["atoms"], 7, None, None, False, A["sigma"])) for x in rw]
    na_bad += sum(int(g[1] != n) for g, n in zip(G, r["natoms"]))
    mean_py = np.mean([g[0] for g in G], 0); mean_pg = np.fromfile(j["out"], "<f8"); os.remove(j["out"])
    eb = float(np.abs(mean_py - mean_pg).max() / mean_py.max())
    iu = [(p, q) for p in range(len(G)) for q in range(p + 1, len(G))]
    dpg = np.sqrt(np.array(r["d2"])); dpy = np.array([np.sqrt(((G[p][0] - G[q][0]) ** 2).sum() / H ** 3) for p, q in iu])
    dst = np.array([float(Dm[rw[p], rw[q]]) for p, q in iu])
    ec, ed = float(np.max(np.abs(dpg - dpy) / dpy)), float(np.max(np.abs(dpg - dst) / dst))
    worst["b"] = max(worst["b"], eb); worst["c"] = max(worst["c"], ec); worst["d"] = max(worst["d"], ed); npairs += len(iu)
    print("  arm %s %-16s %2d members: mean-grid err %.1e of peak | D vs Python %.1e | D vs stored %.1e" % (j["arm"], j["name"], len(j["mem"]), eb, ec, ed))
os.remove(drv); os.remove(W + "rg4_page.json")
nm = sum(len(j["mem"]) for j in jobs)
chk("RG4a atoms per member == build_grid (%d members, %d grids)" % (nm, len(jobs)), na_bad == 0, "%d mismatches" % na_bad)
chk("RG4b cluster-mean grid: max |page - Python| <= 1e-3 x peak", worst["b"] <= 1e-3, "%.2e" % worst["b"])
chk("RG4c within-cluster D, page vs Python float64 build (%d pairs) <= 1e-3 rel" % npairs, worst["c"] <= 1e-3, "%.2e" % worst["c"])
chk("RG4d within-cluster D, page vs stored matrix (%d pairs) <= 5e-3 rel" % npairs, worst["d"] <= 5e-3, "%.2e" % worst["d"])
pd.DataFrame(rows).to_csv(RP + "checks/RG4_voxels%s.csv" % SUF, index=False)
if not all(r["passed"] for r in rows):
    raise SystemExit("RG4 FAILED")
