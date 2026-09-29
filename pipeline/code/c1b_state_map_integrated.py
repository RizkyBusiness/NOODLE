"""STAGE C prep - repair the integration, then name the states as helper programs.

TWO defects in the run this replaces, both reported rather than worked around:

1. The packaged mapper's harmony call fails on this stack and the failure is caught by
   a broad except, so the run completed on UNINTEGRATED PCs and printed its own
   caveat ("donor structure in the map is then a real confound"). Cause is a version
   incompatibility, not a harmony failure: harmonypy 2.0.0 returns Z_corr/result() as
   (n_cells, n_pcs), while scanpy's wrapper assigns harmony_out.Z_corr.T to obsm, so
   anndata rejects a (n_pcs, n_cells) value. Harmony itself converges. Fixed here by
   calling harmonypy directly with the correct orientation. The skill's script is left
   untouched; this script supersedes its embedding.

2. Stage A named every leiden cluster by argmax of the program scores with no margin
   guard, which is how most of the cells came out labelled
   Tr1: a broadly activated cluster wins that argmax on Tbx21/Maf/Ccl5 alone. Here a
   cluster is named for a program only if the program's mean score is positive AND its
   margin over the runner-up exceeds MARGIN; otherwise it keeps a marker-gene name and
   is excluded from the polarisation arm rather than silently counted as a T-helper
   subset.

The clustering itself is inherited from the packaged mapper's QC, HVG selection
(receptor genes excluded) and PCA, read back from embedding.h5ad - not recomputed.

Usage:  python c1b_state_map_integrated.py [config_main.json]
Writes: tables/cell_states.csv.gz            (the Stage C input, integrated)
        tables/cell_states.unintegrated.csv.gz  (packaged run, preserved)
        tables/embedding_integrated.h5ad
        tables/C1b_cluster_programs.csv, C1b_state_markers.csv,
        C1b_batch_mixing.csv, C1b_vs_cells.csv, figures/C1b_state_map.png
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd

import imgt
import paths  # package: config (reference receptor, data-set settings)

MARGIN = 0.02

# Tr1 is defined by the IL-10 / LAG-3 / CD49b(Itga2) axis (Gagliani 2013, Nat Med;
# Roncarolo 2018, Immunity). An earlier version of this dictionary put Tbx21, Ccl5 and
# Gzmb inside the Tr1 programme, which guaranteed that any T-bet+ cytotoxic effector
# cluster won the Tr1 argmax whether or not it made IL-10 — that is how clusters 0 and 5
# came to be called Tr1 despite Il10 in 1.4% / 0.8% of cells. Those three effector genes
# now sit in Th1-eff, which is what such a cluster actually is.
PROGRAMS = {
    "Treg":   ["Foxp3", "Ikzf2", "Il2ra", "Ctla4", "Tnfrsf18", "Folr4", "Izumo1r"],
    "Tr1":    ["Il10", "Lag3", "Itga2", "Havcr2", "Maf", "Pdcd1"],
    # Naming convention: T-bet+ IFN-gamma+ is Th1, full stop. Whether a Th1 cluster is
    # additionally cytotoxic, tissue-resident or innate-like is a separate question and
    # belongs in a subset annotation, not in the state name (see c1b2_relabel_states.py).
    "Th1":    ["Tbx21", "Ifng", "Cxcr3", "Il12rb2", "Ccl5", "Nkg7"],
    "Th17":   ["Rorc", "Il17a", "Il17f", "Il23r", "Ccr6", "Il1r1"],
    "Tfh":    ["Cxcr5", "Bcl6", "Ascl2", "Izumo1r"],
    "Naive":  ["Sell", "Tcf7", "Lef1", "Ccr7", "Bach2", "Igfbp4"],
    "Cycling": ["Mki67", "Top2a", "Birc5", "Stmn1"],
}

# Clusters carrying non-T transcripts in most of their cells are ambient-RNA artefacts
# (Paneth-cell Lyz1, pancreatic acinar Prss2), not helper states. They are named
# "ambient" and excluded downstream rather than given a T-helper label.
AMBIENT_GENES = ["Lyz1", "Prss2", "Prss1", "Reg3b", "Reg3g", "Defa24"]
AMBIENT_FRAC = 0.50


def same_batch_fraction(A, rep, batch, k):
    """mean fraction of a cell's k nearest neighbours from its own batch.
    1/n_batches is perfect mixing; 1.0 is complete separation."""
    from sklearn.neighbors import NearestNeighbors
    X = A.obsm[rep]
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    idx = nn.kneighbors(X, return_distance=False)[:, 1:]
    b = pd.factorize(A.obs[batch])[0]
    return float((b[idx] == b[:, None]).mean())


def main(config_path="config_main.json"):
    cfg = imgt.load(config_path)
    import scanpy as sc
    import harmonypy
    out, tp = cfg["out_dir"], cfg["transcriptome"]
    sc.settings.verbosity = 1
    fig_dir = os.path.join(out, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    src = os.path.join(out, "cell_states.csv.gz")
    keep = os.path.join(out, "cell_states.unintegrated.csv.gz")
    if os.path.exists(src) and not os.path.exists(keep):
        shutil.copy2(src, keep)
        print("preserved packaged (unintegrated) table as %s" % os.path.basename(keep))
    prev = pd.read_csv(keep, low_memory=False) if os.path.exists(keep) else None

    A = sc.read_h5ad(os.path.join(out, "embedding.h5ad"))
    print("embedding read: %d cells x %d genes | libraries %d"
          % (A.n_obs, A.n_vars, A.obs.library.nunique()))

    # ---- 1. integration, correct orientation
    ho = harmonypy.run_harmony(A.obsm["X_pca"], A.obs[["library"]], ["library"],
                               max_iter_harmony=20, verbose=False)
    Z = np.asarray(ho.result())
    assert Z.shape == A.obsm["X_pca"].shape, "harmony returned %s" % (Z.shape,)
    A.obsm["X_pca_harmony"] = Z
    rep = "X_pca_harmony"

    mix = pd.DataFrame([
        dict(representation="X_pca (unintegrated)", same_library_neighbour_fraction=
             round(same_batch_fraction(A, "X_pca", "library", tp["k"]), 4)),
        dict(representation="X_pca_harmony", same_library_neighbour_fraction=
             round(same_batch_fraction(A, rep, "library", tp["k"]), 4)),
        dict(representation="perfect mixing (1/n_libraries)",
             same_library_neighbour_fraction=round(1 / A.obs.library.nunique(), 4))])
    mix.to_csv(os.path.join(out, "C1b_batch_mixing.csv"), index=False)
    print(mix.to_string(index=False))

    sc.pp.neighbors(A, n_neighbors=tp["k"], use_rep=rep)
    sc.tl.leiden(A, resolution=tp["resolution"], key_added="cluster",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.tsne(A, use_rep=rep, perplexity=tp["perplexity"])
    print("clusters on the integrated representation: %d" % A.obs.cluster.nunique(),
          flush=True)

    # ---- 2. programme scores on log-normalised counts, with a margin guard
    R = A.raw.to_adata()
    R.obs["cluster"] = A.obs.cluster.values
    for name, genes in PROGRAMS.items():
        present = [g for g in genes if g in R.var_names]
        if len(present) < len(genes):
            print("  %-8s markers absent from matrix: %s"
                  % (name, ",".join(sorted(set(genes) - set(present)))))
        sc.tl.score_genes(R, present, score_name="score_" + name,
                          random_state=cfg["controls"]["seed"])
    cols = ["score_" + k for k in PROGRAMS]
    cl = R.obs.groupby("cluster", observed=True)[cols].mean()

    amb_present = [g for g in AMBIENT_GENES if g in R.var_names]
    if amb_present:
        E = pd.DataFrame(
            np.asarray((R[:, amb_present].X > 0).todense())
            if hasattr(R[:, amb_present].X, "todense")
            else (R[:, amb_present].X > 0),
            columns=amb_present)
        E["cluster"] = R.obs.cluster.values
        amb_frac = E.groupby("cluster", observed=True)[amb_present].mean().max(axis=1)
        amb_gene = E.groupby("cluster", observed=True)[amb_present].mean().idxmax(axis=1)
    else:
        amb_frac = pd.Series(0.0, index=cl.index)
        amb_gene = pd.Series("", index=cl.index)
    amb_frac = amb_frac.reindex(cl.index).fillna(0.0)
    amb_gene = amb_gene.reindex(cl.index).fillna("")

    sc.tl.rank_genes_groups(A, "cluster", method="wilcoxon")
    mk = sc.get.rank_genes_groups_df(A, group=None)
    mk.to_csv(os.path.join(out, "C1b_state_markers.csv"), index=False)
    top = (mk[mk.pvals_adj < 0.05].sort_values("scores", ascending=False)
           .groupby("group").names.apply(lambda s: "/".join(s.head(2))))

    srt = np.sort(cl.values, axis=1)
    best = cl.idxmax(axis=1).str.replace("score_", "", regex=False)
    margin = srt[:, -1] - srt[:, -2]
    named, why = {}, {}
    for i, c in enumerate(cl.index):
        if amb_frac.loc[c] >= AMBIENT_FRAC:
            named[c] = "ambient"
            why[c] = ("%s detected in %.1f%% of cells — non-T ambient RNA"
                      % (amb_gene.loc[c], 100 * amb_frac.loc[c]))
        elif srt[i, -1] > 0 and margin[i] > MARGIN:
            named[c] = best.loc[c]
            why[c] = "programme argmax, margin %.3f" % margin[i]
        else:
            named[c] = "unassigned c%s: %s" % (c, top.get(c, "no markers"))
            why[c] = ("no programme separates: top %s score %.3f, margin %.3f"
                      % (best.loc[c], srt[i, -1], margin[i]))
    prog = cl.round(3).copy()
    prog["n_cells"] = R.obs.cluster.value_counts().reindex(cl.index).values
    prog["state"] = [named[c] for c in cl.index]
    prog["ambient_frac"] = amb_frac.reindex(cl.index).round(3).values
    prog["is_ambient"] = [named[c] == "ambient" for c in cl.index]
    prog["assignment"] = [why[c] for c in cl.index]
    prog.to_csv(os.path.join(out, "C1b_cluster_programs.csv"))
    print("\ncluster -> state\n", prog[["n_cells", "state", "assignment"]]
          .to_string())

    A.obs["state"] = [named[c] for c in A.obs.cluster]
    A.obs["state_conf"] = np.nan
    A.obs["state_source"] = "programme-scored cluster (integrated, margin>%.2f)" % MARGIN
    A.obs["condition"] = pd.Series(A.obs.donor.astype(str).values).str.split(
        "_").str[0].values

    T = pd.DataFrame(dict(
        cell=A.obs_names, clone_key=A.obs.clone_key.values, donor=A.obs.donor.values,
        condition=A.obs.condition.values, library=A.obs.library.astype(str).values,
        cluster=A.obs.cluster.astype(str).values, state=A.obs.state.astype(str).values,
        state_conf=A.obs.state_conf.values,
        state_source=A.obs.state_source.astype(str).values,
        x=A.obsm["X_tsne"][:, 0], y=A.obsm["X_tsne"][:, 1]))
    T.to_csv(src, index=False)
    A.write_h5ad(os.path.join(out, "embedding_integrated.h5ad"))
    print("\nstate composition of paired cells\n",
          T.state.value_counts().to_string())
    print("\nstate x condition\n", pd.crosstab(T.state, T.condition).to_string())

    # ---- how much did the receptor-gene leakage and the naming guard change?
    SA = pd.read_csv(os.path.join(cfg["out_dir"].replace("tables", "cells"),
                                  "cell_obs_annotated.csv.gz"),
                     low_memory=False, usecols=["cell", "state_provisional"])
    bc = T.cell.astype(str).str.split("|").str[-1]
    j = pd.DataFrame(dict(gex=bc.values, new=T.state.values)).merge(
        SA.rename(columns={"cell": "gex", "state_provisional": "cells"}),
        on="gex", how="inner")
    ct = pd.crosstab(j.cells, j.new)
    ct.to_csv(os.path.join(out, "C1b_vs_cells.csv"))
    from sklearn.metrics import adjusted_rand_score
    ari = adjusted_rand_score(j.cells, j.new)
    same = float((j.cells == j.new).mean())
    print("\nvs Stage A provisional labels on the same %d cells: identical label "
          "%.1f%%, ARI %.3f" % (len(j), 100 * same, ari))
    if prev is not None:
        pj = pd.DataFrame(dict(cell=T.cell, new=T.state)).merge(
            prev[["cell", "state"]].rename(columns={"state": "unintegrated"}),
            on="cell", how="inner")
        print("vs the unintegrated packaged run: %d cells, ARI %.3f"
              % (len(pj), adjusted_rand_score(pj.unintegrated, pj.new)))

    # ---- figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({
        "figure.dpi": 200, "savefig.dpi": 200, "font.size": 7, "axes.titlesize": 8,
        "axes.labelsize": 7, "legend.fontsize": 6, "xtick.labelsize": 6,
        "ytick.labelsize": 6, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "figure.facecolor": "white", "savefig.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42, "axes.titlelocation": "left"})
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.4))
    order = T.state.value_counts().index.tolist()
    pal = plt.get_cmap("tab10")
    for i, s in enumerate(order):
        m = (T.state == s).values
        axs[0].scatter(T.x[m], T.y[m], s=0.6, lw=0, color=pal(i % 10),
                       rasterized=True)
        cx, cy = float(np.median(T.x[m])), float(np.median(T.y[m]))
        axs[0].annotate("%s (%d)" % (s.split(":")[0], int(m.sum())), (cx, cy),
                        fontsize=5.5, ha="center",
                        bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6))
    for i, c in enumerate(paths.required_list("conditions")):
        m = (T.condition == c).values
        axs[1].scatter(T.x[m], T.y[m], s=0.6, lw=0, label="%s (%d)" % (c, int(m.sum())),
                       color=["#2166ac", "#b2182b"][i], rasterized=True)
    axs[1].legend(loc="upper right", markerscale=8, frameon=False)
    axs[0].set_title("helper state, receptor-gene-free integrated map", pad=8)
    axs[1].set_title("housing condition", pad=8)
    for a in axs:
        a.set_xticks([])
        a.set_yticks([])
        a.set_xlabel("t-SNE 1")
        a.set_ylabel("t-SNE 2")
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, "C1b_state_map.png"), bbox_inches="tight")
    print("wrote %s/cell_states.csv.gz and figures/C1b_state_map.png" % out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config_main.json")
