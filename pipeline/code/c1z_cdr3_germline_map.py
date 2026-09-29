"""Which CDR3 positions are germline-inherited, and which are junctional?

CDR3 is not one thing. Its N-terminal residues are encoded by the V gene, its C-terminal
residues by the J gene, and only the middle is junctional (D segment, N/P additions) and
therefore genuinely clone-specific. Clustering on "CDR3" is partly clustering on germline.

The boundaries are derived from the data rather than from a reference: for each V gene,
walk in from the start of CDR3 and record the fraction of that gene's receptors carrying
the modal residue; the V-encoded stretch is the run where that fraction stays above
CONS. Same from the C-terminal end for each J gene.

Writes tables/C1z_vj_cdr3_extents.csv     per V and J gene: how many CDR3 positions it fixes
[package: the per-cluster position map of the original script (C1z_position_origin.csv, read from the old C1g
labels) is not part of this package]
"""
import pandas as pd
import numpy as np
import paths

CONS = 0.95
MINN = 10

F = pd.read_csv(paths.src("tables/folding_set.csv.gz"))
ID = pd.read_csv(paths.src("tables/C1w_receptor_identity.csv.gz"))
F = F.merge(ID[["clone_id", "is_representative"]], on="clone_id")
F = F[F.is_representative]

rows = []
for chain, vcol, jcol, ccol in [("A", "v_A", "j_A", "cdr3_A"), ("B", "v_B", "j_B", "cdr3_B")]:
    for gene, grp in F.groupby(vcol):
        if len(grp) < MINN:
            continue
        s = grp[ccol].astype(str).values
        n = 0
        for i in range(min(map(len, s))):
            col = [x[i] for x in s]
            if max(pd.Series(col).value_counts()) / len(col) >= CONS:
                n += 1
            else:
                break
        rows.append(dict(chain=chain, segment="V", gene=gene, n_receptors=len(grp),
                         positions_fixed=n,
                         motif="".join(pd.Series([x[i] for x in s]).value_counts().index[0]
                                       for i in range(n))))
    for gene, grp in F.groupby(jcol):
        if len(grp) < MINN:
            continue
        s = grp[ccol].astype(str).values
        n = 0
        for i in range(1, min(map(len, s)) + 1):
            col = [x[-i] for x in s]
            if max(pd.Series(col).value_counts()) / len(col) >= CONS:
                n += 1
            else:
                break
        rows.append(dict(chain=chain, segment="J", gene=gene, n_receptors=len(grp),
                         positions_fixed=n,
                         motif="".join(pd.Series([x[-i] for x in s]).value_counts().index[0]
                                       for i in range(n, 0, -1))))
E = pd.DataFrame(rows)
E.to_csv(paths.dst("tables/C1z_vj_cdr3_extents.csv"), index=False)
print("V/J germline extents within CDR3 (%.0f%% consensus, genes with >=%d receptors)"
      % (100 * CONS, MINN))
for seg in ("V", "J"):
    d = E[E.segment == seg]
    print("  %s: median %d positions fixed (range %d-%d), n=%d genes"
          % (seg, d.positions_fixed.median(), d.positions_fixed.min(),
             d.positions_fixed.max(), len(d)))

