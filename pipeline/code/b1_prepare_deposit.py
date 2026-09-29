"""Stage B1 - build the prepared V(D)J deposits the pipeline will consume.

The packaged pipeline has no notion of invariant-subset receptors or of
dual-chain cells (vdj.max_chains_per_cell is declared in the config schema but
build_vdomains.py never reads it; it keeps the highest-UMI contig per chain and
sets a `multichain` flag instead). Both exclusions therefore have to be applied
to the input table, which is what this script does - explicitly, with counts.

It also injects the two columns the deposit lacks: `mouse` (the donor identity,
recovered from the hashtags; the deposited `donor` column is the literal "pool"
for every contig) and `lane` (from the aggr barcode suffix).

Writes two deposits so the sensitivity analysis is a config switch, not a rerun
of the reasoning:
  raw/prepared_main/contigs.csv.gz    QC + singlet, no iNKT/MAIT, no dual-chain
  raw/prepared_sens/contigs.csv.gz    QC + singlet only (invariant and dual kept)
"""
import os

import pandas as pd

import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
SEG = ["fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]


def main():
    C = pd.read_csv(
        paths.deposit("contigs"),
        low_memory=False)
    obs = pd.read_csv(paths.src("cells/cell_obs_annotated.csv.gz"),
                      low_memory=False)

    C["lane"] = C.barcode.str.rsplit("-", n=1).str[-1].astype(int)
    rows = [dict(step="contigs deposited", contigs=len(C), cells=C.barcode.nunique())]

    # --- donor identity: hashtag singlets that also passed GEX QC
    # cell_obs_annotated.csv.gz contains only cells that already passed GEX QC,
    # so membership in it IS the QC filter; the singlet call is applied here.
    keep = obs[obs.hto_class == "Singlet"]
    m = keep.set_index("cell")[["mouse", "condition", "state_provisional",
                                "iNKT_TCR", "MAIT_TCR", "dual_any"]]
    C = C.merge(m, left_on="barcode", right_index=True, how="inner")
    rows.append(dict(step="cells that are HTO singlets and pass GEX QC",
                     contigs=len(C), cells=C.barcode.nunique()))

    sens = C.copy()

    # --- PRIMARY deposit population: invariant receptors RETAINED (they are the
    # positive control for the geometry method, and exclusion by transcription is
    # deferred until Zbtb16/Klrb1c/Il2rb/Cd244 evidence exists), but the dual-chain
    # exclusion is kept because a dual-chain cell has an ambiguous combining surface.
    prim = C[~C.dual_any].copy()
    rows.append(dict(step="PRIMARY population: invariants retained, dual-chain excluded",
                     contigs=len(prim), cells=prim.barcode.nunique()))

    # --- exclusion 1: invariant-subset receptors
    inv = C.iNKT_TCR | C.MAIT_TCR
    C = C[~inv]
    rows.append(dict(step="after excluding iNKT/MAIT V-J usage",
                     contigs=len(C), cells=C.barcode.nunique()))

    # --- exclusion 2: dual-chain cells (ambiguous combining surface)
    C = C[~C.dual_any]
    rows.append(dict(step="after excluding dual-chain cells",
                     contigs=len(C), cells=C.barcode.nunique()))

    # --- both chains present with complete framework + loop segments
    def finalize(df, name):
        d = df[df[SEG].notna().all(axis=1)]
        n = d.groupby(["barcode", "chain"]).size().unstack(fill_value=0)
        both = n.index[(n.get("TRA", 0) >= 1) & (n.get("TRB", 0) >= 1)]
        d = d[d.barcode.isin(both)]
        out = os.path.join(ROOT, "raw", name)
        os.makedirs(out, exist_ok=True)
        d.to_csv(os.path.join(out, "contigs.csv.gz"), index=False)
        return d

    main_d = finalize(C, "prepared_main")
    rows.append(dict(step="MAIN deposit: complete segments, both chains",
                     contigs=len(main_d), cells=main_d.barcode.nunique()))
    sens_d = finalize(sens, "prepared_sens")
    rows.append(dict(step="SENS deposit: invariant and dual-chain retained",
                     contigs=len(sens_d), cells=sens_d.barcode.nunique()))
    prim_d = finalize(prim, "prepared_primary")
    rows.append(dict(step="PRIMARY deposit: invariants retained, no dual-chain",
                     contigs=len(prim_d), cells=prim_d.barcode.nunique()))

    qc = pd.DataFrame(rows)
    qc.to_csv(os.path.join(ROOT, "tables/B1_deposit_prep_qc.csv"), index=False)
    print(qc.to_string(index=False))

    for nm, d in (("main", main_d), ("sens", sens_d), ("primary", prim_d)):
        cl = d.assign(k=d.v_gene.str.split("*").str[0] + "|"
                      + d.j_gene.str.split("*").str[0] + "|" + d.cdr3)
        key = cl.groupby(["barcode", "chain"]).k.first().unstack()
        key = key.dropna(subset=["TRA", "TRB"])
        ck = key.TRA + "|" + key.TRB
        print("%s: %d cells, %d clonotypes, %d mice, %d lanes"
              % (nm, d.barcode.nunique(), ck.nunique(), d.mouse.nunique(),
                 d.lane.nunique()))


if __name__ == "__main__":
    main()
