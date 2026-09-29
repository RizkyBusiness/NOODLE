"""STAGE A3 (HTO half) - mouse identity and SPF/GF condition per cell.

The contig table's donor/origin columns are both the literal "pool" (cells were
pooled during sorting), so mouse identity exists only in the hashtag counts.
This reimplements Seurat HTODemux closely enough for a condition/mouse label:

  CLR-normalise HTO counts across cells, k-means the CLR matrix into
  n_HTO + 1 groups, take for each HTO the group with the lowest mean as that
  HTO's negative population, and set the threshold at the 0.99 empirical
  quantile of the negative population's counts. Cells positive for exactly one
  HTO are singlets; 0 -> negative, >1 -> doublet.

DEVIATION from the authors, on the record: Seurat's HTODemux fits a negative
binomial to the negative population and takes the 0.99 quantile of the FITTED
distribution; this uses the empirical quantile instead. Labels are therefore
close to but not identical to the authors' demultiplexing.

Usage: python a3_hto_demux.py
"""
import os

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans

import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
H5 = paths.deposit("matrix")
HTOF = paths.deposit("hto_features")
# where a hashtag's sample name carries the condition: data-set config "hashtag_names" (no default)
HN = {k: v for k, v in paths.dataset().get("hashtag_names", {}).items() if not k.startswith("_")}
HN_SEP, HN_COND = HN.get("separator"), HN.get("condition_field")
if not HN_SEP or str(HN_SEP).startswith("<") or type(HN_COND) is not int or HN_COND < 0:
    raise SystemExit('data-set config: set "hashtag_names" = {"separator": ..., "condition_field": <0-based index>}: the '
                     'condition is the part condition_field of the hashtag sample name split at the separator (see '
                     'config/upstream.example.json)')
OUT = os.path.join(ROOT, "cells")
os.makedirs(OUT, exist_ok=True)
POS_Q = 0.99
SEED = 0


def main():
    ad = sc.read_10x_h5(H5, gex_only=False)
    ad.var_names_make_unique()
    print("h5: %d cells x %d features" % ad.shape)
    print("feature types:", ad.var.feature_types.value_counts().to_dict())

    hto_names = pd.read_csv(HTOF)
    print("HTO panel:", dict(zip(hto_names.id, hto_names.name)))

    is_ab = ad.var.feature_types != "Gene Expression"
    hto = ad[:, is_ab].copy()
    print("HTO features in matrix:", list(hto.var_names))

    X = np.asarray(hto.X.todense(), dtype=float)          # cells x HTO
    counts = pd.DataFrame(X, index=hto.obs_names, columns=hto.var_names)

    # CLR across cells, per feature (Seurat margin=2 default for HTO)
    def clr(v):
        lv = np.log1p(v)
        return lv - lv[v > 0].mean() if (v > 0).any() else lv
    C = counts.apply(clr, axis=0)

    km = KMeans(n_clusters=C.shape[1] + 1, n_init=10, random_state=SEED)
    grp = km.fit_predict(C.values)

    thresholds, pos = {}, pd.DataFrame(index=counts.index)
    for h in counts.columns:
        means = pd.Series(C[h].values, index=grp).groupby(level=0).mean()
        neg = counts.loc[grp == means.idxmin(), h]
        thr = float(np.quantile(neg.values, POS_Q))
        thresholds[h] = thr
        pos[h] = counts[h] > thr

    npos = pos.sum(axis=1)
    lab = pd.DataFrame(index=counts.index)
    lab["n_positive"] = npos
    lab["hto_class"] = np.where(npos == 1, "Singlet",
                                np.where(npos == 0, "Negative", "Doublet"))
    top = counts.idxmax(axis=1)
    lab["hto_id"] = np.where(npos == 1, pos.idxmax(axis=1), "")
    lab["hto_top_by_count"] = top
    id2name = dict(zip(hto_names.id, hto_names.name))
    # matrix feature names may already be the sample names
    lab["mouse"] = lab.hto_id.map(lambda x: id2name.get(x, x) if x else "")
    lab["condition"] = lab.mouse.str.split(HN_SEP).str[HN_COND]    # config "hashtag_names"
    lab["lane"] = [int(b.rsplit("-", 1)[-1]) for b in lab.index]

    pd.Series(thresholds).to_csv(os.path.join(OUT, "A3_hto_thresholds.csv"),
                                 header=["count_threshold"])
    lab.reset_index(names="cell").to_csv(
        os.path.join(OUT, "hto_assignment.csv.gz"), index=False)

    print("\n== HTO classification ==")
    print(lab.hto_class.value_counts().to_string())
    sing = lab[lab.hto_class == "Singlet"]
    print("\n== singlets per mouse ==")
    print(sing.mouse.value_counts().sort_index().to_string())
    print("\n== singlets: condition x lane ==")
    print(pd.crosstab(sing.condition, sing.lane).to_string())
    return lab


if __name__ == "__main__":
    main()
