"""STAGE C prep - the three inputs transcriptome_map.py needs, and why each exists.

Stage A produced state_provisional from a scanpy run whose variable-gene pool was
NOT filtered for receptor genes (code/a3_state_composition.py calls
highly_variable_genes with no exclusion). The packaged mapper excludes
^TR[ABGD][VDJC] and ^IG[HKL][VDJC] from the pool for a specific reason: co-clustered
receptors in this repertoire share TRAV/TRBV 18.8x and 8.2x above chance (B4c), so if
V-segment transcripts drive the embedding then a structural cluster predicts state
partly through the expression of its own V genes and Stage C partly predicts itself.
state_provisional is therefore not usable as the Stage C state input, and this rebuilds
it on a receptor-gene-free embedding.

Writes: raw/gex_main.h5ad          (deposited 10x h5 -> h5ad, the loader takes no .h5)
        raw/gex_lane_map.csv       (aggr barcode suffix 1..8 -> lane; lane is the
                                    technical batch. mouse x lane is fully crossed -
                                    every mouse appears in all 8 lanes - so lane is a
                                    pure batch variable and mouse stays biological)
        tables/paired_cells.csv.gz  (+ cell_gex column; the two assays use
                                    different barcode namespaces - VDJ ids are
                                    "lane|barcode-lane", expression ids are
                                    "barcode-lane" - and the mapper asks for exactly
                                    this column when they differ)

Usage: python c1a_prepare_state_inputs.py
"""
import os

import pandas as pd
import scanpy as sc

import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
H5 = paths.deposit("matrix")
SOURCE_NAME = "gex"          # load_sources prefixes obs_names with "<name>|"
OUT = os.path.join(ROOT, "tables")


def main():
    A = sc.read_10x_h5(H5, gex_only=True)
    A.var_names_make_unique()
    print("GEX h5: %d cells x %d genes" % A.shape)
    rec = A.var_names.str.upper().str.match(r"^TR[ABGD][VDJC]|^IG[HKL][VDJC]")
    print("receptor genes present in the matrix: %d (these are what the mapper drops "
          "from the variable-gene pool)" % int(rec.sum()))

    h5ad = os.path.join(ROOT, "raw", "gex_main.h5ad")
    A.write_h5ad(h5ad)

    # aggr suffix -> lane. suffix i is lane i in this deposit (config _project_notes).
    lanes = pd.DataFrame(dict(library=["lane%d" % i for i in range(1, 9)]))
    lanes.to_csv(os.path.join(ROOT, "raw", "gex_lane_map.csv"), index=False)

    p = os.path.join(OUT, "paired_cells.csv.gz")
    PC = pd.read_csv(p, low_memory=False)
    bc = PC.cell.astype(str).str.split("|").str[-1]
    PC["cell_gex"] = ["%s|%s" % (SOURCE_NAME, b) for b in bc]
    gex_names = pd.Index(["%s|%s" % (SOURCE_NAME, b) for b in A.obs_names])
    hit = PC.cell_gex.isin(set(gex_names)).sum()
    print("paired cells whose expression barcode is in the matrix: %d / %d"
          % (int(hit), len(PC)))
    assert hit > 0.95 * len(PC), "cell_gex join lost %d cells" % (len(PC) - hit)
    PC.to_csv(p, index=False)
    print("wrote %s (+cell_gex), %s, %s" % (p, h5ad, "raw/gex_lane_map.csv"))


if __name__ == "__main__":
    main()
