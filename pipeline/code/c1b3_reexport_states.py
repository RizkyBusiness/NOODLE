"""Re-export every derived per-receptor / per-cell state table from the corrected
tables/cell_states.csv.gz produced by c1b2_relabel_states.py.

Ambient clusters (c9 Lyz1+, c11 Prss2+) are dropped, not relabelled: they are not
T-cell states, so they should not contribute cells to any state, to any receptor's
dominant-state vote, or to the UMAP state colouring.

Rewrites: tables/_receptor_states_all.csv.gz
          tables/umap_coords.csv.gz          (state column only)
[package: the C1g label re-export of the original script is not part of this package]
"""
import os
import pandas as pd
import paths

OUT = "tables"
C = pd.read_csv(paths.src(os.path.join(OUT, "cell_states.csv.gz")))
print("cells total %d, ambient %d" % (len(C), int(C.is_ambient.sum())))
T = C[~C.is_ambient.astype(bool)].copy()
print("cells kept  %d" % len(T))
print(T.state.value_counts().to_string())

# ---- per-receptor dominant state over ALL of that clone's cells
g = (T.dropna(subset=["clone_key"])
       .groupby(["clone_key", "state"]).size().rename("n").reset_index())
tot = g.groupby("clone_key").n.sum().rename("cells")
top = g.sort_values("n", ascending=False).drop_duplicates("clone_key")
don = (T.dropna(subset=["clone_key"]).groupby("clone_key")
         .agg(donor=("donor", lambda s: s.value_counts().index[0]),
              n_donors=("donor", "nunique")))
R = (top.set_index("clone_key").rename(columns={"n": "state_n"})
        .join(tot).join(don))
R["state_frac"] = (R.state_n / R.cells).round(4)
R = R.reset_index()[["clone_key", "state", "state_n", "state_frac",
                     "donor", "n_donors", "cells"]]
R.to_csv(paths.dst(os.path.join(OUT, "_receptor_states_all.csv.gz")), index=False)
print("receptors %d" % len(R))

# ---- umap coords: drop ambient cells, carry corrected state
U = pd.read_csv(paths.src(os.path.join(OUT, "umap_coords.csv.gz")))
U = U.drop(columns=["state"]).merge(C[["cell", "state", "is_ambient"]], on="cell",
                                    how="left")
U = U[~U.is_ambient.astype(bool)].drop(columns=["is_ambient"])
U.to_csv(paths.dst(os.path.join(OUT, "umap_coords.csv.gz")), index=False)
print("umap cells %d" % len(U))

