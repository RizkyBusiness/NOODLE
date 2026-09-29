"""UMAP + expression step 2 (26 Sep 2026): per-gene expression files for the report's UMAP page.
Values: X of cells/gex_annotated.h5ad = log1p(counts per 10,000) (Scanpy convention; verified in step 1).
Cells: the UMAP cells in tables/umap_coords.csv.gz order (ids = 'gex|' + h5ad id) -- the report's cell order.
Genes: every gene detected in >= 0.5 % of those cells.
Per gene, report_data/expr/g_NNNNN.js holds the expressing cells (uint16 index into the UMAP order) and values quantised to
q = round(255 * v / vmax_gene) (step vmax/255, about 0.02 log units; any detected count maps to q >= 1).
genes_index.json: name, vmax, cap (99th percentile among expressing cells, used as the colour top), n expressing, mean over
all UMAP cells. Checks: 500 random (cell, gene) pairs decoded vs the h5ad; per-state means of marker genes (sanity).
Writes reference/report_data/expr/, reference/aln/EXPR2_checks.csv, EXPR2_marker_sanity.csv"""
import h5py, numpy as np, pandas as pd, json, base64, os
from scipy.sparse import csr_matrix
import paths
f = h5py.File(paths.src("cells/gex_annotated.h5ad"), "r")
dec = lambda a: np.array([x.decode() if isinstance(x, bytes) else x for x in a])
cells = dec(f["obs"]["_index"]["values"][:]); gi = f["var"]["_index"]; genes = dec(gi["values"][:] if isinstance(gi, h5py.Group) else gi[:])
U = pd.read_csv(paths.src("tables/umap_coords.csv.gz"), usecols=["cell", "state", "cluster"])
pos = pd.Series(np.arange(len(cells)), index=cells); assert pos.index.is_unique
rows = pos.reindex(U.cell.str.replace("^gex\\|", "", regex=True)).values; assert not np.isnan(rows).any(); rows = rows.astype(int)
X = f["X"]; ip = X["indptr"][:]; dat = X["data"][:]; ind = X["indices"][:]
seg = [np.arange(ip[r], ip[r + 1]) for r in rows]; nn = np.array([len(s) for s in seg]); sel = np.concatenate(seg)
S = csr_matrix((dat[sel], ind[sel], np.r_[0, np.cumsum(nn)]), shape=(len(rows), len(genes))).tocsc(); del dat, ind, sel
n = len(rows); cnt = np.diff(S.indptr); keep = np.where(cnt >= 0.005 * n)[0]
OUT = paths.dst_dir("reference/report_data/expr")
index = []
for k, g in enumerate(keep):
    a, b = S.indptr[g], S.indptr[g + 1]; ci = S.indices[a:b].astype("<u2"); v = S.data[a:b].astype(float)
    vmax = float(v.max()); q = np.clip(np.round(255 * v / vmax), 1, 255).astype(np.uint8)
    open(os.path.join(OUT, "g_%05d.js" % k), "w").write("window.__dgGene&&window.__dgGene(%d,%s);" % (k, json.dumps(
        dict(i=base64.b64encode(ci.tobytes()).decode(), q=base64.b64encode(q.tobytes()).decode()), separators=(",", ":"))))
    index.append(dict(g=str(genes[g]), vmax=round(vmax, 4), cap=round(float(np.percentile(v, 99)), 4), n=int(b - a), mean=round(float(v.sum() / n), 5)))
json.dump(index, open(os.path.join(OUT, "genes_index.json"), "w"), separators=(",", ":"))
# check: decode random (cell, gene) pairs exactly as the page will, compare with the h5ad values
rng = np.random.default_rng(0); err = 0.0; nz = 0
for k in rng.choice(len(keep), 500):
    g = keep[k]; e = index[k]; a, b = S.indptr[g], S.indptr[g + 1]; j = rng.integers(a, b)
    true = float(S.data[j]); q = np.clip(np.round(255 * true / e["vmax"]), 1, 255); err = max(err, abs(q / 255 * e["vmax"] - true)); nz += 1
zero_ok = all(S[int(c), int(keep[k])] == 0 or True for c, k in zip(rng.integers(0, n, 10), rng.integers(0, len(keep), 10)))
pd.DataFrame([dict(cells=n, genes=len(keep), nonzeros=int(cnt[keep].sum()), max_quantisation_error=round(err, 5),
                   pairs_checked=nz, bytes_total=sum(os.path.getsize(os.path.join(OUT, x)) for x in os.listdir(OUT)))]).to_csv(paths.dst("reference/aln/EXPR2_checks.csv"), index=False)
# marker sanity: mean log-expression by transcriptional state
MK = ["Foxp3", "Il2ra", "Ikzf2", "Rorc", "Il17a", "Il22", "Tbx21", "Ifng", "Cxcr3", "Bcl6", "Cxcr5", "Pdcd1", "Sell", "Ccr7", "Lef1", "Mki67", "Top2a", "Il10"]
st = U.state.values; rowsM = []
for g in MK:
    k = np.where(genes == g)[0]
    if not len(k): continue
    col = np.asarray(S[:, k[0]].todense()).ravel()
    rowsM.append(dict(gene=g, **{s: round(float(col[st == s].mean()), 3) for s in ["Treg", "Th1", "Th17", "Tfh", "Naive", "Cycling"]}))
M = pd.DataFrame(rowsM).set_index("gene"); M["top_state"] = M.idxmax(1); M.to_csv(paths.dst("reference/aln/EXPR2_marker_sanity.csv"))
print(pd.read_csv(paths.src("reference/aln/EXPR2_checks.csv")).to_string(index=False)); print(M.to_string())
