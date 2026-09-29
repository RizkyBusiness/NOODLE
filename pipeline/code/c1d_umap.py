"""UMAP on the Harmony-integrated PCs; caches coordinates for the plotting scripts.

Reads  tables/embedding_integrated.h5ad  (obsm/X_pca_harmony, obs/cluster, obs/state)
Writes tables/umap_coords.csv.gz
"""
import h5py, numpy as np, pandas as pd, umap, os
import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
H = f"{ROOT}/tables/embedding_integrated.h5ad"
f = h5py.File(H, "r")
X = np.asarray(f["/obsm/X_pca_harmony"][:])
cell = np.array([x.decode() for x in f["/obs/_index/values"][:]])
def cat(p):
    return np.array([x.decode() for x in f[f"/obs/{p}/categories"][:]])[f[f"/obs/{p}/codes"][:]]
cl, st, don, cond = cat("cluster"), cat("state"), cat("donor"), cat("condition")
f.close()
print("running UMAP on %s" % (X.shape,), flush=True)
emb = umap.UMAP(n_neighbors=15, min_dist=0.3, metric="euclidean",
                random_state=0).fit_transform(X)
pd.DataFrame(dict(cell=cell, cluster=cl, state=st, donor=don, condition=cond,
                  umap1=emb[:, 0], umap2=emb[:, 1])).to_csv(
    f"{ROOT}/tables/umap_coords.csv.gz", index=False)
print("wrote tables/umap_coords.csv.gz", flush=True)
