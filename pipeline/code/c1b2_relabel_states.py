"""Correct the Stage C state labels on the existing cell_states table.

Two corrections, both at cluster level.

1. The pipeline's "Tr1" programme contained Tbx21, Ccl5 and Gzmb — T-bet-driven effector
   genes — so any T-bet+ cluster won the Tr1 argmax whether or not it made IL-10. Clusters
   0 and 5 were called Tr1 on that basis (Il10 in 1.4% and 0.8% of their cells). They are
   not Tr1. The genuine IL-10+ LAG-3+ CD49b+ cells sit in the Foxp3+ Treg clusters.

   Naming convention adopted: a cluster that is T-bet+ and IFN-gamma+ is Th1. Clusters 0,
   5 and 6 all qualify, so all three are Th1. Whether a Th1 cluster is additionally
   cytotoxic, tissue-resident or innate-like is a separate question and is carried in the
   `subset` column rather than folded into the state name — see C1q_cluster_gene_fraction
   for the evidence behind those descriptors.

2. Clusters 9 and 11 are dominated by non-T ambient RNA (Lyz1 from Paneth cells, Prss2
   from pancreatic acinar cells) and are not T-cell states at all. They are flagged
   is_ambient and excluded downstream rather than relabelled.

Preserves the original call in `state_pipeline` and backs the input up to
tables/cell_states.prerelabel.csv.gz.

Writes tables/cell_states.csv.gz (in place) and tables/C1b2_state_relabel_map.csv
"""
import os
import shutil
import pandas as pd

OUT = "tables"
SRC = os.path.join(OUT, "cell_states.csv.gz")
BAK = os.path.join(OUT, "cell_states.prerelabel.csv.gz")

# cluster -> (state, subset, reason)
RELABEL = {
    "0":  ("Th1", "cytotoxic/resident",
           "Tbx21, Ifng, Ccl5, Nkg7, Cd160, Itgae; "
           "Il12rb2, Il18r1 — T-bet+ effector, low IL-12/IL-18 receptor axis. "
           "Il10, Lag3, Itga2 low — not Tr1"),
    "5":  ("Th1", "innate-like",
           "Tbx21, Ifng; Zbtb16, Klrb1c, Klrd1, Il2rb, Xcl1 "
           "— innate-like signature, but polyclonal TCR (TRAV11/TRAJ18 rare), not iNKT/MAIT. "
           "Il10 low — not Tr1"),
    "6":  ("Th1", "IL-12/IL-18 axis",
           "Tbx21, Ifng, Il12rb2, Il18r1, Il18rap, Ccr5 "
           "— classical Th1; unchanged by this script except for the explicit subset tag"),
    "9":  ("ambient", "ambient",
           "Lyz1+ — Paneth-cell ambient RNA, not a T-cell state"),
    "11": ("ambient", "ambient",
           "Prss2+ — pancreatic acinar ambient RNA"),
}

if not os.path.exists(BAK):
    shutil.copy(SRC, BAK)
    print("backed up -> %s" % BAK)
else:
    print("backup already present, reusing it as the input")
C = pd.read_csv(BAK, low_memory=False)
C["cluster"] = C.cluster.astype(str)

if "state_pipeline" not in C.columns:
    C["state_pipeline"] = C.state
print("before:", C.state_pipeline.value_counts().to_dict())

C["state"] = [RELABEL[c][0] if c in RELABEL else s
              for c, s in zip(C.cluster, C.state_pipeline)]
C["subset"] = [RELABEL[c][1] if c in RELABEL else "" for c in C.cluster]
C["is_ambient"] = C.state == "ambient"

C.to_csv(SRC, index=False)
print("after :", C.state.value_counts().to_dict())
print("subsets:", C[C.subset != ""].subset.value_counts().to_dict())
print("ambient cells flagged: %d" % int(C.is_ambient.sum()))

rows = []
for c, g in C.groupby("cluster"):
    sp = g.state_pipeline.iloc[0]
    st = g.state.iloc[0]
    rows.append(dict(cluster=c, n_cells=len(g), state_pipeline=sp, state=st,
                     subset=g.subset.iloc[0], changed=sp != st,
                     reason=RELABEL[c][2] if c in RELABEL else "unchanged"))
pd.DataFrame(rows).sort_values("cluster", key=lambda s: s.astype(int)).to_csv(
    os.path.join(OUT, "C1b2_state_relabel_map.csv"), index=False)
print("wrote C1b2_state_relabel_map.csv")
