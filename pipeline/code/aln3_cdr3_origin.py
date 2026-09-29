"""Alignment step 3b (26 Sep 2026): germline origin of each CDR3 junction residue (V-encoded / junctional / J-encoded).

Source (chosen by the user): the project's own data-derived map tables/C1z_vj_cdr3_extents.csv -- for each V and
J gene with >= 10 receptors, how many residues at that end of the CDR3 junction (IMGT 104-118, C...F/W) carry the gene's
modal residue in >= 95 % of its receptors (code/c1z_cdr3_germline_map.py). Approximate: a residue conserved by selection
counts as germline; genes with < 10 receptors have no estimate ('origin unknown').
Mapped onto each molecule's IMGT junction from the folded model (ALN1_sequences.npz): the first nV residues are V-encoded,
the last nJ J-encoded, the rest junctional (N/P/D). Where nV + nJ exceeds the junction length the overlap is flagged.
Writes reference/aln/ALN3_cdr3_origin.csv (clone_id, chain, gene_v, gene_j, junction_len, nV, nJ, overlap)
"""
import numpy as np, pandas as pd
import paths
A = np.load(paths.src("reference/aln/ALN1_sequences.npz"), allow_pickle=True)
E = pd.read_csv(paths.src("tables/C1z_vj_cdr3_extents.csv")).set_index(["chain", "segment", "gene"]).positions_fixed
F = pd.read_csv(paths.src("tables/folding_set.csv.gz"), low_memory=False, usecols=["clone_id", "v_A", "j_A", "v_B", "j_B"]).set_index("clone_id")
ids = A["clone_id"]; mol = A["mol"]; ch = A["chain"]; num = A["num"]
jun = (num >= 104) & (num <= 118)
rows = []
for chn in "AB":
    m = jun & (ch == chn); jl = np.bincount(mol[m], minlength=len(ids))
    for k, c in enumerate(ids):
        gv, gj = F.loc[c, "v_" + chn], F.loc[c, "j_" + chn]
        nv = E.get((chn, "V", gv), np.nan); nj = E.get((chn, "J", gj), np.nan)
        rows.append(dict(clone_id=c, chain=chn, gene_v=gv, gene_j=gj, junction_len=int(jl[k]),
                         nV=nv, nJ=nj, overlap=bool(nv == nv and nj == nj and nv + nj > jl[k])))
R = pd.DataFrame(rows); R.to_csv(paths.dst("reference/aln/ALN3_cdr3_origin.csv"), index=False)
for chn in "AB":
    r = R[R.chain == chn]
    print("%s: chains %d | V extent known %.1f %% | J extent known %.1f %% | nV median %.0f | nJ median %.0f | V/J overlap %d | junctional residues median %.0f"
          % (chn, len(r), 100 * r.nV.notna().mean(), 100 * r.nJ.notna().mean(), r.nV.median(), r.nJ.median(), r.overlap.sum(),
             (r.junction_len - r.nV - r.nJ).clip(lower=0).median()))
