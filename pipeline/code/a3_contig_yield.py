"""STAGE A3 (contig half) - pairing yield, dual chains, invariant subsets.

Reads the aggregated filtered_contig_annotations.csv.gz. Lane comes from the
cellranger-aggr barcode suffix (-1..-8); the donor/origin columns are both the
literal "pool" and carry no mouse identity, so SPF/GF and mouse must come from
HTO demultiplexing of the h5 (done separately in a3_hto_demux.py).

Writes cells/*.csv and prints the numbers that decide feasibility.

Usage: python a3_contig_yield.py
"""
import os

import pandas as pd

import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
CONTIG = paths.deposit("contigs")
OUT = os.path.join(ROOT, "cells")
os.makedirs(OUT, exist_ok=True)

# canonical mouse invariant usage (skill brief): iNKT TRAV11-TRAJ18,
# MAIT TRAV1-TRAJ33. Match on gene stem so TRAV11D / TRAV11-1 style calls count.
INKT_V, INKT_J = ("TRAV11",), ("TRAJ18",)
MAIT_V, MAIT_J = ("TRAV1",), ("TRAJ33",)


def stem(series):
    """gene call -> comparable stem, dropping allele and D-locus duplicate suffix."""
    return (series.astype(str).str.split("*").str[0]
            .str.replace(r"D$", "", regex=True))


def main():
    df = pd.read_csv(CONTIG)
    print("contigs: %d rows x %d cols" % df.shape)
    assert set(df.chain.unique()) <= {"TRA", "TRB"}, df.chain.unique()

    df["lane"] = df.barcode.str.rsplit("-", n=1).str[-1].astype(int)
    df["cell"] = df.barcode  # aggr barcodes are already lane-unique
    df["v_stem"] = stem(df.v_gene)
    df["j_stem"] = stem(df.j_gene)

    # every row in a filtered_contig file is is_cell/high_confidence/productive
    for c in ("is_cell", "high_confidence", "productive", "full_length"):
        assert df[c].astype(str).str.lower().eq("true").all(), c

    a = df[df.chain == "TRA"]
    b = df[df.chain == "TRB"]
    na = a.groupby("cell").size()
    nb = b.groupby("cell").size()

    cells = pd.DataFrame({"n_TRA": na, "n_TRB": nb}).fillna(0).astype(int)
    cells["lane"] = cells.index.str.rsplit("-", n=1).str[-1].astype(int)
    cells["paired"] = (cells.n_TRA >= 1) & (cells.n_TRB >= 1)
    cells["dual_TRA"] = cells.n_TRA >= 2
    cells["dual_TRB"] = cells.n_TRB >= 2
    cells["dual_any"] = cells.dual_TRA | cells.dual_TRB

    # invariant subsets, judged on the alpha chain of paired cells
    inkt_cells = set(a[a.v_stem.isin(INKT_V) & a.j_stem.isin(INKT_J)].cell)
    mait_cells = set(a[a.v_stem.isin(MAIT_V) & a.j_stem.isin(MAIT_J)].cell)
    cells["iNKT_TCR"] = cells.index.isin(inkt_cells)
    cells["MAIT_TCR"] = cells.index.isin(mait_cells)

    # single-chain-each cells get a clonotype key from observed V/J/CDR3
    a1 = a.sort_values("umis", ascending=False).drop_duplicates("cell")
    b1 = b.sort_values("umis", ascending=False).drop_duplicates("cell")
    key = (a1.set_index("cell")[["v_stem", "j_stem", "cdr3"]]
             .join(b1.set_index("cell")[["v_stem", "j_stem", "cdr3"]],
                   lsuffix="_a", rsuffix="_b", how="inner"))
    key["clone_key"] = (key.v_stem_a + "|" + key.j_stem_a + "|" + key.cdr3_a + "|"
                        + key.v_stem_b + "|" + key.j_stem_b + "|" + key.cdr3_b)
    cells["clone_key"] = key.clone_key

    per_lane = (cells.groupby("lane")
                  .agg(cells_with_any_chain=("paired", "size"),
                       paired=("paired", "sum"),
                       alpha_only=("n_TRB", lambda s: int((s == 0).sum())),
                       beta_only=("n_TRA", lambda s: int((s == 0).sum())),
                       dual_TRA=("dual_TRA", "sum"), dual_TRB=("dual_TRB", "sum"),
                       iNKT=("iNKT_TCR", "sum"), MAIT=("MAIT_TCR", "sum")))
    per_lane["distinct_paired_clonotypes"] = (
        cells[cells.paired].groupby("lane").clone_key.nunique())
    per_lane.to_csv(os.path.join(OUT, "A3_pairing_per_lane.csv"))

    pc = cells[cells.paired]
    totals = pd.Series(dict(
        contigs=len(df), cells_with_any_chain=len(cells),
        productive_TRA_cells=int((cells.n_TRA >= 1).sum()),
        productive_TRB_cells=int((cells.n_TRB >= 1).sum()),
        paired_cells=int(cells.paired.sum()),
        distinct_paired_clonotypes=int(pc.clone_key.nunique()),
        dual_TRA_cells=int(cells.dual_TRA.sum()),
        dual_TRB_cells=int(cells.dual_TRB.sum()),
        dual_any_cells=int(cells.dual_any.sum()),
        dual_among_paired=int((cells.paired & cells.dual_any).sum()),
        iNKT_TCR_paired=int((cells.paired & cells.iNKT_TCR).sum()),
        MAIT_TCR_paired=int((cells.paired & cells.MAIT_TCR).sum()),
    ))
    totals.to_csv(os.path.join(OUT, "A3_pairing_totals.csv"), header=["n"])

    # clone size distribution, for the identical-receptor ceiling
    cs = pc.clone_key.value_counts()
    dist = pd.Series(dict(
        clonotypes=len(cs), singletons=int((cs == 1).sum()),
        largest_clone=int(cs.iloc[0]),
        top_clone_frac=round(cs.iloc[0] / len(pc), 4),
        cells_in_expanded_clones=int(cs[cs >= 2].sum()),
        frac_cells_in_expanded=round(cs[cs >= 2].sum() / len(pc), 4)))
    dist.to_csv(os.path.join(OUT, "A3_clone_size_summary.csv"), header=["value"])

    cells.reset_index(names="cell").to_csv(
        os.path.join(OUT, "cell_chain_index.csv.gz"), index=False)

    # V-gene call coverage against the ANARCI mouse tables is checked in Stage B;
    # here just record the observed gene vocabulary
    vg = (df.groupby(["chain", "v_stem"]).size().rename("contigs")
            .reset_index().sort_values(["chain", "contigs"], ascending=[True, False]))
    vg.to_csv(os.path.join(OUT, "A3_v_gene_usage.csv"), index=False)

    print("\n== totals ==");            print(totals.to_string())
    print("\n== per lane ==");          print(per_lane.to_string())
    print("\n== clone sizes ==");       print(dist.to_string())
    print("\n== distinct V genes observed: TRA %d, TRB %d ==" % (
        vg[vg.chain == "TRA"].v_stem.nunique(), vg[vg.chain == "TRB"].v_stem.nunique()))
    return cells


if __name__ == "__main__":
    main()
