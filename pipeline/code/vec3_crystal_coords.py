"""reference step 3a (device): CDR3 arc points and the 10 landmark Ca, own coordinates, for the accepted
crystal-model benchmark pairs (same matching and parsing as landmarks/code/lm3_crystal_floor.py).
Writes reference/out/V3_crystal_coords.npz"""
import numpy as np, pandas as pd, os, runpy
import paths
src = open(paths.pkg("pipeline/code/lm3_crystal_floor.py")).read().split("rows = []")[0]   # definitions only (parse, B, LMPOS)
g = {}; exec(src, g)
E, AC, LC, AM, LMm = [], [], [], [], []
for _, r in g["B"].iterrows():
    fc, ac, _ = g["parse"](os.path.join(g["BEN"], "fixed_%s_crystal.pdb" % r.entry))
    fm, am, _ = g["parse"](os.path.join(g["BEN"], r.model))
    E.append(r.entry); AC.append(ac); AM.append(am)
    LC.append([fc[k] for k in g["LMPOS"]]); LMm.append([fm[k] for k in g["LMPOS"]])
np.savez(paths.dst("reference/out/V3_crystal_coords.npz"), entry=np.array(E), arc_c=np.array(AC), lm_c=np.array(LC), arc_m=np.array(AM), lm_m=np.array(LMm))
print(len(E), np.array(AC).shape)
