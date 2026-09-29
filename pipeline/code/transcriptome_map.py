"""STAGE 3a/3b - one embedding over ALL paired cells, and a state name for every cell.

The map is built once, over every cell with a paired receptor - not per cluster,
not per subset. Every later figure and every concordance test reads the same
coordinates and the same state labels, so a difference between two figures is a
difference in the receptors, not in the embedding they were drawn on.

Two deliberate choices:

  receptor genes are excluded from the variable-gene pool. Leaving TRAV/TRBV
    transcripts in makes the embedding partly a receptor embedding, and the
    structure-versus-state test then predicts itself.
  states are named, not numbered. Either transferred from a reference labelling
    by nearest neighbours in the shared PC space, or built from each cluster's
    own top markers. Numbered clusters are not comparable between runs; names
    that carry their markers are.

Usage:  python transcriptome_map.py [config.json]
Writes: <out_dir>/cell_states.csv.gz   (cell, clone_key, donor, state, x, y, ...)
        <out_dir>/gex_qc.csv, <out_dir>/state_markers.csv,
        <out_dir>/embedding.h5ad, <out_dir>/figures/state_map.png
"""
import os
import re
import sys

import numpy as np
import pandas as pd

import imgt


def load_sources(cfg):
    import anndata as ad
    import scanpy as sc
    tp = cfg["transcriptome"]
    if not tp["sources"]:
        raise SystemExit(
            "transcriptome.sources is empty. Each entry is "
            '{"name": ..., "path": ..., "kind": "10x_mtx"|"h5ad", '
            '"lib_map": optional aggregation table mapping barcode suffix -> library}')
    parts = []
    for s in tp["sources"]:
        if s["kind"] == "10x_mtx":
            a = sc.read_10x_mtx(s["path"], var_names="gene_symbols", cache=False)
        elif s["kind"] == "h5ad":
            a = ad.read_h5ad(s["path"])
        else:
            raise SystemExit("unknown source kind %r" % s["kind"])
        a.obs["source"] = s["name"]
        if s.get("lib_map") and os.path.exists(s["lib_map"]):
            lm = pd.read_csv(s["lib_map"])
            col = [c for c in lm.columns if c.lower() in
                   ("library_id", "sample_id", "donor", "library")][0]
            suf = a.obs_names.str.split("-").str[-1].astype(int)
            a.obs["library"] = pd.Categorical(
                [lm[col].iloc[i - 1] if 1 <= i <= len(lm) else "unknown" for i in suf])
        else:
            a.obs["library"] = s["name"]
        a.obs_names = ["%s|%s" % (s["name"], b) for b in a.obs_names]
        parts.append(a)
        print("  %s: %d cells x %d genes" % (s["name"], a.n_obs, a.n_vars), flush=True)
    A = ad.concat(parts, join="outer", label="batch_source", index_unique=None) \
        if len(parts) > 1 else parts[0]
    A.var_names_make_unique()
    return A


def main(config_path="config.json"):
    cfg = imgt.load(config_path)
    import scanpy as sc
    out, tp = cfg["out_dir"], cfg["transcriptome"]
    sc.settings.verbosity = 1
    A = load_sources(cfg)
    qc = [dict(stage="loaded", cells=A.n_obs, genes=A.n_vars)]

    PC = pd.read_csv(os.path.join(out, "paired_cells.csv.gz"), low_memory=False)
    key = "cell_gex" if "cell_gex" in PC.columns else "cell"
    pc = PC.drop_duplicates(key).set_index(key)
    keep = A.obs_names.isin(pc.index)
    if keep.sum() < 0.01 * A.n_obs:
        raise SystemExit(
            "only %d of %d expression barcodes matched paired_cells.%s. The two "
            "assays use different barcode namespaces here - add a 'cell_gex' column "
            "to paired_cells.csv.gz carrying the expression-side barcode."
            % (int(keep.sum()), A.n_obs, key))
    A = A[keep].copy()
    for c in ("clone_key", "donor"):
        A.obs[c] = pc.loc[A.obs_names, c].values
    qc.append(dict(stage="paired receptor", cells=A.n_obs, genes=A.n_vars))
    print("cells with a paired receptor: %d" % A.n_obs, flush=True)

    A.var["mito"] = A.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(A, qc_vars=["mito"], inplace=True, percent_top=None,
                               log1p=False)
    A = A[(A.obs.n_genes_by_counts >= tp["min_genes"])
          & (A.obs.pct_counts_mito <= tp["max_pct_mito"])].copy()
    sc.pp.filter_genes(A, min_cells=tp["min_cells_per_gene"])
    qc.append(dict(stage="QC filtered", cells=A.n_obs, genes=A.n_vars))
    print("after QC: %d cells x %d genes" % (A.n_obs, A.n_vars), flush=True)

    A.layers["counts"] = A.X.copy()
    sc.pp.normalize_total(A, target_sum=1e4)
    sc.pp.log1p(A)
    sc.pp.highly_variable_genes(A, n_top_genes=tp["n_hvg"], batch_key="library"
                                if A.obs.library.nunique() > 1 else None)
    pat = re.compile("|".join(tp["exclude_gene_patterns"]))
    drop = A.var_names.str.upper().str.match(pat)
    n_drop = int((A.var.highly_variable & drop).sum())
    A.var.loc[drop, "highly_variable"] = False
    print("variable genes %d | receptor genes removed from the pool %d"
          % (int(A.var.highly_variable.sum()), n_drop), flush=True)

    A.raw = A
    sc.pp.scale(A, max_value=10, zero_center=True)
    sc.tl.pca(A, n_comps=tp["n_pcs"], svd_solver="arpack", use_highly_variable=True)
    rep = "X_pca"
    if tp["integrate"] == "harmony" and A.obs.library.nunique() > 1:
        try:
            sc.external.pp.harmony_integrate(A, "library", max_iter_harmony=20)
            rep = "X_pca_harmony"
        except Exception as e:
            print("harmony unavailable (%s) - continuing on unintegrated PCs; "
                  "donor structure in the map is then a real confound and the "
                  "cross-donor concordance subset is the only safe headline"
                  % type(e).__name__, flush=True)
    sc.pp.neighbors(A, n_neighbors=tp["k"], use_rep=rep)
    sc.tl.leiden(A, resolution=tp["resolution"], key_added="cluster",
                 flavor="igraph", n_iterations=2, directed=False)
    sc.tl.tsne(A, use_rep=rep, perplexity=tp["perplexity"])
    print("clusters %d" % A.obs.cluster.nunique(), flush=True)

    sc.tl.rank_genes_groups(A, "cluster", method="wilcoxon", use_raw=True)
    mk = sc.get.rank_genes_groups_df(A, None)
    mk.to_csv(os.path.join(out, "state_markers.csv"), index=False)
    top = (mk[mk.pvals_adj < 0.05].sort_values("scores", ascending=False)
           .groupby("group").names.apply(lambda s: "/".join(s.head(2))))
    ref = tp["reference_labels"]
    if ref and os.path.exists(ref):
        R = pd.read_csv(ref)
        lab_col = [c for c in R.columns if c not in ("cell", "barcode")][0]
        cell_col = "cell" if "cell" in R.columns else R.columns[0]
        R = R.set_index(cell_col)[lab_col]
        shared = A.obs_names.intersection(R.index)
        if len(shared) < 20:
            raise SystemExit("reference_labels matched only %d cells" % len(shared))
        from sklearn.neighbors import KNeighborsClassifier
        X = A.obsm[rep]
        m = A.obs_names.isin(shared)
        clf = KNeighborsClassifier(n_neighbors=tp["k"]).fit(
            X[m], R.loc[A.obs_names[m]].values)
        A.obs["state"] = clf.predict(X)
        A.obs["state_conf"] = clf.predict_proba(X).max(1)
        A.obs["state_source"] = np.where(m, "reference", "transferred")
        print("states transferred from %s: %d labels, %d cells seeded"
              % (os.path.basename(ref), A.obs.state.nunique(), int(m.sum())), flush=True)
    else:
        A.obs["state"] = ["c%s: %s" % (c, top.get(c, "unnamed"))
                          for c in A.obs.cluster]
        A.obs["state_conf"] = np.nan
        A.obs["state_source"] = "marker-named cluster"
        print("states named from cluster markers: %d" % A.obs.state.nunique(), flush=True)

    T = pd.DataFrame(dict(
        cell=A.obs_names, clone_key=A.obs.clone_key.values, donor=A.obs.donor.values,
        library=A.obs.library.astype(str).values, cluster=A.obs.cluster.astype(str).values,
        state=A.obs.state.astype(str).values, state_conf=A.obs.state_conf.values,
        state_source=A.obs.state_source.astype(str).values,
        x=A.obsm["X_tsne"][:, 0], y=A.obsm["X_tsne"][:, 1]))
    T.to_csv(os.path.join(out, "cell_states.csv.gz"), index=False)
    pd.DataFrame(qc + [dict(stage="embedded", cells=A.n_obs, genes=A.n_vars)]).to_csv(
        os.path.join(out, "gex_qc.csv"), index=False)
    A.write_h5ad(os.path.join(out, "embedding.h5ad"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    states = sorted(T.state.unique())
    cmap = plt.get_cmap("tab20")
    for i, s in enumerate(states):
        d = T[T.state == s]
        ax.scatter(d.x, d.y, s=2, lw=0, color=cmap(i % 20), label="%s (n=%d)" % (s, len(d)))
        ax.annotate(s.split(":")[0], (d.x.median(), d.y.median()), fontsize=6,
                    ha="center", va="center",
                    bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.5))
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.legend(frameon=False, fontsize=5.5, loc="center left", bbox_to_anchor=(1.0, 0.5),
              markerscale=3, handletextpad=0.2)
    ax.set_title("%s - %d paired cells, %d states" % (cfg["dataset"], len(T), len(states)),
                 fontsize=9, loc="left")
    fig.tight_layout()
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    fig.savefig(os.path.join(out, "figures", "state_map.png"), dpi=200, facecolor="white")
    print("wrote %s/cell_states.csv.gz, embedding.h5ad, figures/state_map.png" % out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.json")
