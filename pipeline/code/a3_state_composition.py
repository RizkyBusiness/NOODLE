"""STAGE A3 (state half) - CD4 helper-state composition of the paired-TCR cells.

The authors' own cluster labels exist only inside GSE298371_scRNAseq_EA_Th1.rds
(9.9 Gb), which Stage A4 of the brief says not to download unless there is no
alternative. This therefore produces a PROVISIONAL annotation of our own:
standard scanpy pipeline on the deposited h5, Leiden clusters, then each cluster
labelled by the highest-scoring canonical CD4 program.

These labels are OURS, not the authors'. They are adequate for judging Stage A
feasibility and for the composition table; whether to replace them with the
authors' labels is a decision for the Stage A gate.

Usage: python a3_state_composition.py
"""
import os

import pandas as pd
import scanpy as sc

import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
H5 = paths.deposit("matrix")
OUT = os.path.join(ROOT, "cells")
SEED = 0

# canonical mouse CD4 programs. Tr1 is separated from Th1 by Il10/co-inhibitory
# module on a Tbx21+ background, which is the distinction the paper turns on.
PROGRAMS = {
    "Treg":   ["Foxp3", "Ikzf2", "Il2ra", "Ctla4", "Tnfrsf18", "Folr4", "Izumo1r"],
    "Tr1":    ["Il10", "Tbx21", "Pdcd1", "Lag3", "Havcr2", "Maf", "Ccl5", "Gzmb"],
    "Th1":    ["Tbx21", "Ifng", "Cxcr3", "Klrg1", "Il12rb2"],
    "Th17":   ["Rorc", "Il17a", "Il17f", "Il23r", "Ccr6", "Il1r1"],
    "Tfh":    ["Cxcr5", "Bcl6", "Ascl2", "Izumo1r"],
    "Naive":  ["Sell", "Tcf7", "Lef1", "Ccr7", "Bach2", "Igfbp4"],
    "Cycling": ["Mki67", "Top2a", "Birc5", "Stmn1"],
}


def main():
    ad = sc.read_10x_h5(H5, gex_only=True)
    ad.var_names_make_unique()
    print("GEX: %d cells x %d genes" % ad.shape)

    hto = pd.read_csv(os.path.join(OUT, "hto_assignment.csv.gz")).set_index("cell")
    chain = pd.read_csv(os.path.join(OUT, "cell_chain_index.csv.gz")).set_index("cell")

    ad.obs = ad.obs.join(hto[["hto_class", "mouse", "condition", "lane"]])
    ad.obs = ad.obs.join(chain[["paired", "dual_any", "iNKT_TCR", "MAIT_TCR",
                                "clone_key"]])
    for c in ("paired", "dual_any", "iNKT_TCR", "MAIT_TCR"):
        ad.obs[c] = ad.obs[c].fillna(False).astype(bool)
    ad.obs["clone_key"] = ad.obs.clone_key.fillna("").astype(str)

    ad.var["mt"] = ad.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt"], inplace=True, percent_top=None,
                               log1p=False)
    keep = ((ad.obs.total_counts > 1200) & (ad.obs.n_genes_by_counts > 500)
            & (ad.obs.pct_counts_mt < 5) & (ad.obs.hto_class == "Singlet"))
    print("QC + singlet filter keeps %d / %d cells" % (int(keep.sum()), ad.n_obs))
    ad = ad[keep].copy()

    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=2000)
    sc.pp.pca(ad, n_comps=40, svd_solver="arpack", random_state=SEED)
    sc.pp.neighbors(ad, n_neighbors=15, n_pcs=40, random_state=SEED)
    sc.tl.leiden(ad, resolution=0.4, key_added="leiden", random_state=SEED,
                 flavor="igraph", n_iterations=2, directed=False)
    print("leiden clusters:", ad.obs.leiden.nunique())

    for name, genes in PROGRAMS.items():
        present = [g for g in genes if g in ad.var_names]
        missing = sorted(set(genes) - set(present))
        if missing:
            print("  %-8s missing from matrix: %s" % (name, ",".join(missing)))
        sc.tl.score_genes(ad, present, score_name="score_" + name,
                          random_state=SEED)

    score_cols = ["score_" + k for k in PROGRAMS]
    cl = ad.obs.groupby("leiden", observed=True)[score_cols].mean()
    ad.obs["state_provisional"] = ad.obs.leiden.map(
        cl.idxmax(axis=1).str.replace("score_", "", regex=False))

    cl.round(3).to_csv(os.path.join(OUT, "A3_cluster_program_scores.csv"))

    obs = ad.obs.copy()
    obs.reset_index(names="cell").to_csv(
        os.path.join(OUT, "cell_obs_annotated.csv.gz"), index=False)

    pc = obs[obs.paired]
    comp = (pc.groupby("state_provisional", observed=True)
              .agg(paired_cells=("paired", "size"),
                   clonotypes=("clone_key", "nunique"),
                   iNKT=("iNKT_TCR", "sum"), dual=("dual_any", "sum"))
              .sort_values("paired_cells", ascending=False))
    comp["pct_of_paired"] = (100 * comp.paired_cells / len(pc)).round(1)
    comp.to_csv(os.path.join(OUT, "A3_state_composition_paired.csv"))

    bycond = pd.crosstab(pc.state_provisional, pc.condition)
    bycond.to_csv(os.path.join(OUT, "A3_state_by_condition.csv"))

    cond = (pc.groupby("condition")
              .agg(paired_cells=("paired", "size"),
                   clonotypes=("clone_key", "nunique"),
                   mice=("mouse", "nunique")))
    cond.to_csv(os.path.join(OUT, "A3_paired_by_condition.csv"))

    permouse = (pc.groupby("mouse")
                  .agg(paired_cells=("paired", "size"),
                       clonotypes=("clone_key", "nunique")))
    permouse.to_csv(os.path.join(OUT, "A3_paired_by_mouse.csv"))

    print("\n== cluster program scores (argmax per cluster) ==")
    print(cl.round(2).assign(label=cl.idxmax(axis=1).str.replace("score_", "",
                                                                regex=False)).to_string())
    print("\n== state composition of paired cells (PROVISIONAL labels) ==")
    print(comp.to_string())
    print("\n== state x condition ==");  print(bycond.to_string())
    print("\n== paired by condition =="); print(cond.to_string())
    print("\n== paired by mouse ==");     print(permouse.to_string())
    ad.write(os.path.join(OUT, "gex_annotated.h5ad"))
    return ad


if __name__ == "__main__":
    main()
