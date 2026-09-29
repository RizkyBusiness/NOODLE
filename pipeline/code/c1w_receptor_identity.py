"""Define a receptor by the molecule, not by the gene label.

This pipeline's clone_key was v_stem|j_stem|cdr3 for both chains. Gene labels turn out to
be a poor identity key here: several mouse V genes are indistinguishable in this assay
(TRBV12-2/TRBV13-2 differ by 0 nucleotides over FR1+FR2+FR3), and several more differ by
1-5 nucleotides while encoding the IDENTICAL V-domain protein (TRAV7-4/7D-4/7N-4,
TRAV10/10D, TRAV7-3/7D-3, TRAV6-5/6N-5, TRAV6-3/6D-3, TRAV13D-3/13N-3). Structures are
folded from the amino-acid sequence, so for every purpose in this project two receptors
with the same V-domain protein are the same receptor, whatever the aligner named them.

Genuinely distinct genes are NOT merged: TRAV16 vs TRAV16N differ by 164 nucleotides and
by protein, TRAV4D-4 vs TRAV5D-4 by 102.

Defines:
  prot_key    the full V-domain amino-acid sequence of both chains
              (fwr1+cdr1+fwr2+cdr2+fwr3+cdr3+fwr4, alpha then beta) -- the identity of the
              molecule that gets folded
  v_A_prot    V-gene names collapsed to one label per distinct V-domain protein
  v_B_prot    the same for beta; used for background/stratum matching

and picks one representative receptor per prot_key (the one contributing the most cells).

Writes tables/C1w_receptor_identity.csv.gz
       tables/C1w_vprot_groups.csv
       tables/C1w_summary.csv
"""
import collections

import pandas as pd

OUT = "tables"
F = pd.read_csv(f"{OUT}/folding_set.csv.gz")
C = pd.read_csv(f"{OUT}/cell_states.csv.gz", low_memory=False)

# ---- 1. collapse V-gene names that encode the same V-domain protein
VPROT = {}
rows = []
for ch, vcol in [("A", "v_A"), ("B", "v_B")]:
    cols = [f"fwr1|{ch}", f"cdr1|{ch}", f"fwr2|{ch}", f"cdr2|{ch}", f"fwr3|{ch}"]
    t = F.groupby(vcol)[cols].agg(lambda s: s.astype(str).value_counts().index[0])
    t["prot"] = t[cols].agg("|".join, axis=1)
    grp = collections.defaultdict(list)
    for g, pr in t.prot.items():
        grp[pr].append(g)
    m = {}
    for names in grp.values():
        label = sorted(names)[0] if len(names) == 1 else "/".join(sorted(names))
        for g in names:
            m[g] = label
        if len(names) > 1:
            rows.append(dict(chain=ch, v_group=label, n_names=len(names),
                             names=";".join(sorted(names)),
                             receptors=int(F[F[vcol].isin(names)].shape[0])))
    VPROT[ch] = m
    F[f"v_{ch}_prot"] = F[vcol].map(m)
    print("chain %s: %d gene names -> %d distinct V-domain proteins" % (ch, len(t), len(grp)))
G = pd.DataFrame(rows).sort_values("receptors", ascending=False)
G.to_csv(f"{OUT}/C1w_vprot_groups.csv", index=False)
print("\ncollapsed groups:")
print(G.to_string(index=False) if len(G) else "  none")

# ---- 2. the molecule itself
SEG = lambda ch: [f"fwr1|{ch}", f"cdr1|{ch}", f"fwr2|{ch}", f"cdr2|{ch}",
                  f"fwr3|{ch}", f"cdr3|{ch}", f"fwr4|{ch}"]
F["vdom_A"] = F[SEG("A")].astype(str).agg("", axis=1) if False else F[SEG("A")].astype(str).agg("".join, axis=1)
F["vdom_B"] = F[SEG("B")].astype(str).agg("".join, axis=1)
F["prot_key"] = F.vdom_A + "_" + F.vdom_B

cells = C.dropna(subset=["clone_key"]).groupby("clone_key").size().rename("cells")
F["cells_own"] = F.clone_key.map(cells).fillna(0).astype(int)

# A merged molecule owns ALL the cells of every row that folded into it, and all their
# mice. Keeping only the representative's cells would silently discard the rest, and the
# receptor's dominant transcriptional state has to be computed over the union.
agg = F.groupby("prot_key").agg(cells=("cells_own", "sum"),
                                n_rows=("clone_id", "size"),
                                clone_keys=("clone_key", lambda s: ";".join(sorted(set(s)))))
F = F.drop(columns=[c for c in ["cells"] if c in F.columns]).merge(
    agg[["cells", "n_rows", "clone_keys"]], on="prot_key", how="left")

F = F.sort_values("cells_own", ascending=False)
rep = F.drop_duplicates("prot_key")[["prot_key", "clone_id"]].rename(
    columns={"clone_id": "representative"})
F = F.merge(rep, on="prot_key", how="left")
F["is_representative"] = F.clone_id == F.representative

n_dup = int((~F.is_representative).sum())
print("\n--- receptor identity ---")
print("folded receptors (rows)                : %d" % len(F))
print("distinct V-domain protein molecules    : %d" % F.prot_key.nunique())
print("rows that are the same molecule twice  : %d" % n_dup)
print("cells re-attributed to the kept row    : %d"
      % int(F.loc[~F.is_representative, "cells_own"].sum()))

keep = ["clone_id", "clone_key", "clone_keys", "prot_key", "representative",
        "is_representative", "v_A", "v_B", "v_A_prot", "v_B_prot", "j_A", "j_B",
        "cdr3_A", "cdr3_B", "cells_own", "cells", "n_rows"]
F[keep].to_csv(f"{OUT}/C1w_receptor_identity.csv.gz", index=False)

# ---- 3. per-molecule transcriptional state, over the union of all merged rows
key2prot = {}
for pk, cks in zip(F.prot_key, F.clone_keys):
    for ck in str(cks).split(";"):
        key2prot[ck] = pk
T = C.dropna(subset=["clone_key"]).copy()
T["prot_key"] = T.clone_key.map(key2prot)
T = T[T.prot_key.notna()]
T = T[~T.get("is_ambient", pd.Series(False, index=T.index)).astype(bool)]
g = T.groupby(["prot_key", "state"]).size().rename("n").reset_index()
tot = g.groupby("prot_key").n.sum().rename("cells_state")
top = g.sort_values("n", ascending=False).drop_duplicates("prot_key")
don = T.groupby("prot_key").agg(donor=("donor", lambda s: s.value_counts().index[0]),
                                n_donors=("donor", "nunique"))
S = (top.set_index("prot_key").rename(columns={"n": "state_n"}).join(tot).join(don))
S["state_frac"] = (S.state_n / S.cells_state).round(4)
S = S.reset_index().merge(F[F.is_representative][["prot_key", "clone_id"]],
                          on="prot_key", how="left")
S.rename(columns={"clone_id": "representative"}, inplace=True)
S.to_csv(f"{OUT}/C1w_molecule_states.csv.gz", index=False)
print("molecules with a transcriptional state : %d of %d"
      % (S.prot_key.notna().sum(), F.prot_key.nunique()))

dup = F[~F.is_representative]
same_as = F.set_index("clone_id").representative
pd.DataFrame([dict(folded_rows=len(F), distinct_molecules=F.prot_key.nunique(),
                   duplicate_rows=n_dup,
                   alpha_names_collapsed=int(G[G.chain == "A"].n_names.sum()) if len(G) else 0,
                   beta_names_collapsed=int(G[G.chain == "B"].n_names.sum()) if len(G) else 0,
                   receptors_in_collapsed_alpha=int(G[G.chain == "A"].receptors.sum()) if len(G) else 0,
                   receptors_in_collapsed_beta=int(G[G.chain == "B"].receptors.sum()) if len(G) else 0)]
             ).to_csv(f"{OUT}/C1w_summary.csv", index=False)
print("\nwrote C1w_receptor_identity.csv.gz / C1w_vprot_groups.csv / C1w_summary.csv")
