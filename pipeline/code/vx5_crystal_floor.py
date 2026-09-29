"""VX5 (the voxel pipeline): crystal floor at production settings, against the full-repertoire VX4 cut (DESIGN VX5, 6.5).

The accepted pairs come from lm3 (vxlib.load_lm3, amendment A1.5). Crystal and model grids are built with the
production parameters read from VX2c_build.json through the same frame / box / typing / smearing code as the
repertoire; shared atoms keyed per A3.1 (primary), all atoms reported alongside.
Hard stop: shared-atom crystal-model error over cut > 0.60 -> tool reported not viable, pipeline stops.
usage: python pipeline/code/vx5_crystal_floor.py
"""
import os, sys, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NBP = expected("n_benchmark_pairs")
import vxgrid as vg
import vx2a_pilot as PA

HARD_STOP = 0.60

if __name__ == "__main__":
    C = Checks("VX5")
    P = json.load(open(os.path.join(OUT, "VX2c_build.json")))
    sigma, h = P["sigma"], P["h"]
    same_par = (tuple(P["shape"]) == PA.SHAPE and np.allclose(P["lo"], PA.LO) and P["truncation_sigma"] == 3.0)
    C.add("benchmark built with the production parameters (sigma, voxel, box, truncation)",
          "sigma %.1f, h %.1f, %s" % (sigma, h, P["shape"]), same_par and sigma == 2.0 and h == 1.0, "== VX2c")
    cut = json.load(open(os.path.join(OUT, "VX4_summary.json")))["cut"]
    lm3 = load_lm3(); B = lm3["B"].sort_values("entry").reset_index(drop=True)
    MP = pd.read_csv(os.path.join(OUT, "VX2a_manifest.csv"))
    pil = MP[MP.kind == "crystal_shared"][["entry", "model"]].sort_values("entry").reset_index(drop=True)
    C.add("the %s pairs are lm3's accepted set (and the pilot's)" % EXP_NBP, len(B),
          len(B) == EXP_NBP and B[["entry", "model"]].equals(pil), "identical")

    Gs = np.load(os.path.join(TMP, "VX2b_h1.0_s2.00.f16.npy"), mmap_mode="r")
    rows, bitwise = [], True
    for _, r in B.iterrows():
        pc, pm = os.path.join(PA.BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(PA.BEN, r.model)
        fc, fm = PA.bench_frame(pc, lm3["parse"], lm3["LMPOS"]), PA.bench_frame(pm, lm3["parse"], lm3["LMPOS"])
        keep = {tuple(a) for a in parse_heavy(pc)[1]} & {tuple(a) for a in parse_heavy(pm)[1]}
        g = {}
        for tag, job in (("cs", ("shared", pc, *fc, keep)), ("ms", ("shared", pm, *fm, keep)),
                         ("cf", ("full", pc, *fc, None)), ("mf", ("full", pm, *fm, None))):
            X, W, _, _ = PA.structure(job)
            g[tag] = vg.build(X, W, sigma, h, PA.LO, PA.SHAPE).astype(np.float16)
        # identical to the pilot's sigma-2.0 benchmark grids (same code path, same parameters)
        pr = MP[(MP.entry == r.entry) & (MP.kind == "crystal_shared")].row.iloc[0]
        bitwise &= bool(np.array_equal(g["cs"], Gs[pr]))
        ds, df = np.sqrt(vg.grid_d2(g["cs"], g["ms"], h)), np.sqrt(vg.grid_d2(g["cf"], g["mf"], h))
        rows.append(dict(entry=r.entry, model=r.model, identity=r.identity, shared_atoms=len(keep),
                         d_shared=ds, d_full=df, over_cut_shared=ds / cut, over_cut_full=df / cut))
    C.add("crystal grids identical to the VX2b sigma-2.0 pilot grids (%s crystals, bitwise)" % EXP_NBP, bitwise, bitwise,
          "identical")
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "VX5_crystal_floor.csv"), index=False)
    rs, rf = float(R.d_shared.median() / cut), float(R.d_full.median() / cut)
    S = pd.read_csv(os.path.join(OUT, "VX2b_sweep.csv")).set_index("sigma").loc[2.0]
    C.info("full-repertoire cut (VX4) vs pilot cut (VX2b)", "%.4f vs %.4f (shift %+.2f %%)"
           % (cut, S.cut_p1, 100 * (cut / S.cut_p1 - 1)))
    C.info("median crystal-model distance, shared / all atoms", "%.4f / %.4f" % (R.d_shared.median(), R.d_full.median()))
    C.info("crystal-model error over cut, shared atoms: full build vs pilot", "%.4f vs %.4f (shift %+.4f)"
           % (rs, S.ratio_shared, rs - S.ratio_shared))
    C.info("crystal-model error over cut, all atoms", "%.4f (pilot %.4f)" % (rf, S.ratio_full))
    C.info("pairs within the cut, shared atoms", "%d / %s" % (int((R.over_cut_shared <= 1).sum()), EXP_NBP))
    C.info("existing arms, crystal err/cut (V4_crystal_floor.csv, A6.2)",
           __import__("vxpaths").CFG.get("reference_method", {}).get("floor_context", ""))
    viable = rs <= HARD_STOP
    save_json(dict(cut=cut, ratio_shared=rs, ratio_full=rf, pilot_ratio_shared=float(S.ratio_shared),
                   pilot_cut=float(S.cut_p1), hard_stop=HARD_STOP, viable=bool(viable)),
              os.path.join(OUT, "VX5_summary.json"))
    C.add("HARD STOP (DESIGN 6.5): crystal-model error over cut, shared atoms", "%.4f" % rs, viable, "<= 0.60")
    C.write()
