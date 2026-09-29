"""Stage C helper exports read by the vector method and the report (package; new script).

The project copies of these two files had no producer script (exported during a session). Both are column
subsets of stage B/C outputs, rebuilt here with no new analysis:

  tables/_cell_clone_key.csv.gz   cell, clone_key
        = those two columns of tables/cell_states.csv.gz, same rows, same order
  tables/_slim_receptors.csv.gz   clone_id, clone_key, n_cells, n_donors, v_A, j_A, cdr3_A, v_B, j_B, cdr3_B,
                                      cdr3len_A, cdr3len_B
        = those columns of tables/folding_set.csv.gz, same rows, same order

The earlier copy of _slim_receptors also carried length_class and struct_cluster (labels of a dropped method).
No script in this package reads them, so they are not written.
Run after c1b3_reexport_states.py (cell_states final) and before the stage D/V scripts.
"""
import pandas as pd
import paths

C = pd.read_csv(paths.src("tables/cell_states.csv.gz"), usecols=["cell", "clone_key"], low_memory=False)
C[["cell", "clone_key"]].to_csv(paths.dst("tables/_cell_clone_key.csv.gz"), index=False)
print("_cell_clone_key: %d cells" % len(C))

COLS = ["clone_id", "clone_key", "n_cells", "n_donors", "v_A", "j_A", "cdr3_A", "v_B", "j_B", "cdr3_B", "cdr3len_A", "cdr3len_B"]
F = pd.read_csv(paths.src("tables/folding_set.csv.gz"), usecols=COLS, low_memory=False)[COLS]
F.to_csv(paths.dst("tables/_slim_receptors.csv.gz"), index=False)
print("_slim_receptors: %d receptors" % len(F))
