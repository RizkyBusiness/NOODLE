"""the voxel pipeline report, stage 0: the method-independent base data of the report, built from primary files.

Supplies what the report shares with the reference method's report: the molecule table (ids, V / J genes, CDR3s,
V proteins, state, mouse, cells), the UMAP cells (with their transcriptional cluster), the transcriptional clusters, the
gene index for the expression look-up, the CDR3 germline-origin map and the junction-annotation notes. The code is the
reference report's own data code for these blocks (same inputs, same row order, same rounding), so the result equals the
data embedded in that report.
Inputs: the reference method's molecule labels (config reference_method: dir + labels; its optional n10_labels is used
only to check the molecule order), tables/ (_slim_receptors, C1w_receptor_identity, umap_coords, _cell_clone_key),
reference/aln/ (ALN3_cdr3_origin, ALN1_junction_annotation_differences), reference/report_data/expr/genes_index.json.
Writes <voxel out>/report/work/base_data.json with the keys
  states, mol, cells, tcl, genes, n_state, aln = {origin, junction_notes}
read by rv3_build_data.py. Paths come from config/voxel_config.json.: python pipeline/report/code/rv0_base_data.py
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", "..", "code"))
from vxpaths import VXP, CFG  # project paths from voxel_config.json (see config/)
import os, json
import numpy as np, pandas as pd

_REF = CFG.get("reference_method", {})
W = VXP("voxel_out/report/work/"); os.makedirs(W, exist_ok=True)
STATES = ["Treg", "Th1", "Th17", "Tfh", "Naive", "Cycling"]          # cell-state labels, in the page's index order
GF = CFG.get("report", {}).get("flag_condition", "")
if not GF or str(GF).startswith("<"):
    raise SystemExit("voxel config: set report.flag_condition (the raw condition value shown as the flagged group)")                                                             # condition value flagged 1 in cells.gf
# ---- molecule table
M = pd.read_csv(VXP(_REF.get("dir", "reference/out/") + _REF["labels"]))
if _REF.get("n10_labels"):
    M10 = pd.read_csv(VXP(_REF["n10_labels"]), usecols=["clone_id"])
    assert (M10.clone_id.values == M.clone_id.values).all(), "molecule order differs between the reference label files"
S = pd.read_csv(VXP("tables/_slim_receptors.csv.gz")).set_index("clone_id")
ID = pd.read_csv(VXP("tables/C1w_receptor_identity.csv.gz"), usecols=["clone_key", "prot_key"]).drop_duplicates("clone_key")
Mi = M.reset_index().set_index("prot_key")["index"]
sidx = {s: i for i, s in enumerate(STATES)}
mol = dict(id=M.clone_id.tolist(), vA=S.loc[M.clone_id, "v_A"].tolist(), jA=S.loc[M.clone_id, "j_A"].tolist(), c3A=S.loc[M.clone_id, "cdr3_A"].tolist(),
           vB=S.loc[M.clone_id, "v_B"].tolist(), jB=S.loc[M.clone_id, "j_B"].tolist(), c3B=S.loc[M.clone_id, "cdr3_B"].tolist(),
           vpA=M.v_A_prot.astype(str).tolist(), vpB=M.v_B_prot.astype(str).tolist(),
           st=[sidx.get(s, -1) if isinstance(s, str) else -1 for s in M.state], mouse=M.donor.fillna("").tolist(),
           cells=S.loc[M.clone_id, "n_cells"].astype(int).tolist())
# ---- cells on the UMAP
u = pd.read_csv(VXP("tables/umap_coords.csv.gz"), usecols=["cell", "donor", "condition", "umap1", "umap2", "state"])
ck = pd.read_csv(VXP("tables/_cell_clone_key.csv.gz"))
u = u.merge(ck, on="cell", how="left").merge(ID, on="clone_key", how="left"); u = u[u.umap1.notna() & u.state.isin(STATES)]
u["mi"] = u.prot_key.map(Mi).fillna(-1).astype(int)
cells = dict(x=(u.umap1 * 100).round().astype(int).tolist(), y=(u.umap2 * 100).round().astype(int).tolist(), st=[sidx[s] for s in u.state],
             mi=u.mi.tolist(), gf=(u.condition == GF).astype(int).tolist())
# ---- transcriptional clusters and gene index (umap_coords in the same row order as cells, checked)
U = pd.read_csv(VXP("tables/umap_coords.csv.gz"))
C = cells
assert len(U) == len(C["x"]) and (np.round(U.umap1 * 100).values == np.array(C["x"])).all() and (np.round(U.umap2 * 100).values == np.array(C["y"])).all()
assert (np.array(STATES)[np.array(C["st"])] == U.state.values).all()
C["cl"] = U.cluster.astype(int).tolist()
tc = U.groupby("cluster").agg(state=("state", lambda s: s.iloc[0]), n=("state", "size"), nst=("state", "nunique"))
assert (tc.nst == 1).all()                                             # every transcriptional cluster sits inside one state
tcl = [dict(id=int(k), state=r.state, n=int(r.n)) for k, r in tc.iterrows()]
G = json.load(open(VXP("reference/report_data/expr/genes_index.json")))
genes = {k: [x[k] for x in G] for k in ("g", "vmax", "cap", "n", "mean")}
# ---- CDR3 germline origin (per molecule [nV_alpha, nJ_alpha, nV_beta, nJ_beta]; -1 = gene without an estimate) and junction notes
ids = mol["id"]
jd = pd.read_csv(VXP("reference/aln/ALN1_junction_annotation_differences.csv"))
O = pd.read_csv(VXP("reference/aln/ALN3_cdr3_origin.csv")).pivot(index="clone_id", columns="chain", values=["nV", "nJ"]).reindex(ids)
origin = [[int(v) if v == v else -1 for v in (r[("nV", "A")], r[("nJ", "A")], r[("nV", "B")], r[("nJ", "B")])] for _, r in O.iterrows()]
data = dict(states=STATES, mol=mol, cells=cells, tcl=tcl, genes=genes, n_state=int(M.state.notna().sum()),
            aln=dict(origin=origin, junction_notes=jd.clone_id.tolist()))
def clean(o):
    if isinstance(o, float) and o != o: return None
    if isinstance(o, dict): return {k: clean(v) for k, v in o.items()}
    if isinstance(o, list): return [clean(v) for v in o]
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.bool_,)): return bool(o)
    if isinstance(o, (np.floating,)): return None if np.isnan(o) else float(o)
    return o
s = json.dumps(clean(data), separators=(",", ":")); open(W + "base_data.json", "w").write(s)
print("base data: molecules %d | cells %d | transcriptional clusters %d | genes %d | json MB %.2f" % (len(ids), len(C["x"]), len(tcl), len(G), len(s) / 1e6))
