"""STAGE 1b - build the paired V-domain sequences on whichever branch stage 1a chose.

observed branch   assemble the donor's own V domain from fwr1-fwr4 + cdr1-cdr3
                  (translating if only nucleotide columns exist), or from the
                  aligned/whole contig sequence, then confirm it with ANARCI and
                  keep the IMGT-numbered 1..128 span.
reference branch  splice: germline V allele positions 1..103 + the observed
                  junction (104..118) + germline J positions 119..128. The
                  framework and the germline loops are then reference sequence,
                  not the donor's - which is what vdomain_provenance.py measures.

Either branch records the constant-region call per chain and carries the
constant sequence through when the deposit provides one. The constant domain is
never folded: TCRBuilder2 (like ABodyBuilder) builds the V module only.

Usage:  python build_vdomains.py [config.json]
Writes: <out_dir>/paired_cells.csv.gz    one row per paired cell
        <out_dir>/vdomains.csv.gz        one row per clonotype, both chain seqs
        <out_dir>/folding_set.csv.gz     the subset that will be folded
        <out_dir>/vdomain_build_qc.csv   cells and clonotypes lost at each step
"""
import json
import os
import sys

import numpy as np
import pandas as pd

import imgt

CODONS = {}
for i, aa in enumerate(
        "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"):
    b = "AACCGGTT"
    CODONS["ACGT"[i // 16] + "ACGT"[(i // 4) % 4] + "ACGT"[i % 4]] = aa


def translate(nt):
    if not isinstance(nt, str):
        return ""
    nt = nt.upper().replace("-", "").replace(".", "")
    return "".join(CODONS.get(nt[i:i + 3], "X") for i in range(0, len(nt) - 2, 3))


def col(df, *names):
    """First present column name among `names`, case-insensitively."""
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


def read_primary(cfg, rep):
    from inspect_deposit import walk_deposit
    want = rep["primary_table"]
    for label, kind, reader in walk_deposit(cfg["deposit_dir"]):
        if label == want:
            sep = next(t["separator"] for t in rep["tables"] if t["file"] == label)
            return pd.read_csv(reader(), sep=sep, low_memory=False)
    raise SystemExit("primary table %r no longer readable" % want)


def chain_of(row, chains, chain_col, vcol):
    """Map a contig to one of the configured chain labels via its locus."""
    v = str(row.get(vcol) or "")
    lab = str(row.get(chain_col) or "")
    for c in chains:
        p = imgt.chain_prefix(None, c)
        if lab.upper().startswith(p) or v.upper().startswith(p):
            return c
    return None


def observed_vdomain(row, df, chains_hint):
    """The donor's own V-domain amino-acid sequence, or '' if not assemblable."""
    parts = []
    for seg in ("fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"):
        aa = row.get(seg) if isinstance(row.get(seg), str) else None
        if aa is None or not str(aa).strip():
            aa = row.get(seg + "_aa") if isinstance(row.get(seg + "_aa"), str) else None
        if aa and set(str(aa).upper()) <= set("ACGTN-."):
            aa = None                      # a nucleotide string in an unsuffixed column
        if not aa:
            nt = row.get(seg + "_nt") or row.get(seg + "_nucleotide")
            aa = translate(nt) if isinstance(nt, str) else ""
        parts.append(str(aa or ""))
    if all(parts):
        return "".join(parts)
    for c in ("sequence_alignment_aa", "sequence_aa"):
        if isinstance(row.get(c), str) and row.get(c).strip():
            return row[c].replace("-", "").replace(".", "")
    for c in ("sequence_alignment", "sequence", "sequence_nt"):
        if isinstance(row.get(c), str) and row.get(c).strip():
            return translate(row[c])
    return ""


def anarci_span(seqs, chains, cfg, ncpu=1):
    """Number sequences with ANARCI and return the IMGT 1..128 span of each.

    Returns {index: (numbered_sequence, {position: residue}, chain_type)}.
    """
    from anarci import run_anarci
    allow = sorted({t for c in chains for t in imgt.SHARED_CHAIN_TYPES.get(c, (c,))})
    inp = [(str(i), s) for i, s in seqs.items() if isinstance(s, str) and len(s) > 60]
    if not inp:
        return {}
    _, numbered, details, _ = run_anarci(
        inp, ncpu=ncpu, scheme="imgt", allowed_species=[cfg["species"]],
        allow=set(allow))
    out = {}
    for (name, _), num, det in zip(inp, numbered, details):
        if not num:
            continue
        pos = {p[0]: a for p, a in num[0][0] if a != "-" and p[1] == " "}
        seq = "".join(a for _, a in num[0][0] if a != "-")
        out[int(name)] = (seq, pos, det[0]["chain_type"])
    return out


def reference_vdomain(v_alleles, j_alleles, v_call, j_call, junction_aa):
    """Germline V 1..103 + junction 104..118 + germline J 119..128."""
    gene_v, gene_j = str(v_call).split("*")[0], str(j_call).split("*")[0]
    if gene_v not in v_alleles or gene_j not in j_alleles or not isinstance(junction_aa, str):
        return "", "", ""
    va, vp = imgt.pick_allele(v_alleles[gene_v], v_call if "*" in str(v_call) else None)
    ja, jp = imgt.pick_allele(j_alleles[gene_j], j_call if "*" in str(j_call) else None)
    fr = "".join(vp[p] for p in sorted(vp) if p < imgt.CYS_104)
    tail = "".join(jp[p] for p in sorted(jp) if p > imgt.ANCHOR_118)
    return fr + junction_aa + tail, va, ja


def main(config_path="config.json"):
    cfg = imgt.load(config_path)
    out, chains = cfg["out_dir"], imgt.chain_labels(cfg)
    rep = json.load(open(os.path.join(out, "deposit_inspection.json")))
    source = rep["vdomain_source"]
    if source == "none":
        raise SystemExit("stage 1a found no V(D)J table; nothing to build")
    print("V-DOMAIN SOURCE: %s" % source.upper(), flush=True)

    df = read_primary(cfg, rep)
    df.columns = [c.lower() for c in df.columns]
    q = cfg["vdj"]
    stages = [dict(step="contigs in table", rows=len(df))]
    for c, want in [("is_cell", True), ("high_confidence", q["require_high_confidence"]),
                    ("full_length", q["require_full_length"]),
                    ("productive", q["require_productive"])]:
        if want and c in df.columns:
            df = df[df[c].astype(str).str.lower().isin(["true", "1", "t", "yes"])]
            stages.append(dict(step="after %s" % c, rows=len(df)))
    ucol = col(df, "umis", "umi_count", "duplicate_count")
    if ucol:
        df = df[pd.to_numeric(df[ucol], errors="coerce").fillna(0) >= q["min_umis"]]
        stages.append(dict(step="after umis>=%d" % q["min_umis"], rows=len(df)))

    vcol = col(df, "v_call", "v_gene")
    jcol = col(df, "j_call", "j_gene")
    ccol = col(df, "c_call", "c_gene")
    dcol = col(df, "d_call", "d_gene")
    j3col = col(df, "junction_aa", "cdr3")
    bcol = col(df, "barcode", "cell_id", "cell")
    chcol = col(df, "chain", "locus") or vcol
    if not all([vcol, jcol, j3col, bcol]):
        raise SystemExit("primary table lacks v/j/junction/barcode columns: %s"
                         % [vcol, jcol, j3col, bcol])
    df["chain_label"] = [chain_of(r, chains, chcol, vcol) for _, r in df.iterrows()]
    df = df[df.chain_label.notna()]
    stages.append(dict(step="assigned to a configured chain", rows=len(df)))

    lib = cfg.get("library_col") and col(df, cfg["library_col"])
    df["lib"] = df[lib] if lib else cfg["dataset"]
    donor = col(df, cfg["donor_col"])
    if donor:
        df["donor"] = df[donor]
    elif cfg.get("sample_metadata"):
        sm = pd.read_csv(cfg["sample_metadata"])
        df = df.merge(sm, left_on="lib", right_on=col(sm, cfg["batch_col"], "lib"), how="left")
        df["donor"] = df[col(df, cfg["donor_col"])]
    else:
        df["donor"] = "donor_unknown"
        print("WARNING: no donor column and no sample_metadata - donor-matched nulls "
              "and the cross-donor pair subset will be unavailable", flush=True)
    df["cell"] = df.lib.astype(str) + "|" + df[bcol].astype(str)

    # one contig per cell per chain, highest UMI count; multi-chain cells flagged
    if ucol:
        df = df.sort_values(ucol, ascending=False)
    n_before = df.groupby(["cell", "chain_label"]).size()
    df["multichain"] = df.cell.isin(n_before[n_before > 1].index.get_level_values(0))
    df = df.drop_duplicates(["cell", "chain_label"])

    wide = {}
    for c in chains:
        d = df[df.chain_label == c].set_index("cell")
        keep = {"v": vcol, "j": jcol, "cdr3": j3col, "c": ccol, "d": dcol}
        w = pd.DataFrame(index=d.index)
        for k, cc in keep.items():
            w["%s_%s" % (k, c)] = d[cc] if cc else None
        w["multichain_%s" % c] = d.multichain
        for extra in ("fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3",
                      "fwr1_nt", "fwr2_nt", "fwr3_nt", "fwr4_nt", "cdr1_nt",
                      "cdr2_nt", "cdr3_nt", "sequence_alignment_aa", "sequence_aa",
                      "sequence_alignment", "sequence", "sequence_nt"):
            if extra in d.columns:
                w["%s|%s" % (extra, c)] = d[extra]
        w["donor"] = d.donor
        w["lib"] = d.lib
        wide[c] = w
    P = wide[chains[0]].join(wide[chains[1]], how="inner", rsuffix="_dup")
    P = P.loc[:, ~P.columns.str.endswith("_dup")].reset_index()
    stages.append(dict(step="cells with both chains", rows=len(P)))
    print("paired cells: %d over %d donors" % (len(P), P.donor.nunique()), flush=True)

    P["clone_key"] = ["|".join(str(r["%s_%s" % (k, c)]) for c in chains
                               for k in ("v", "j", "cdr3")) for _, r in P.iterrows()]
    P.to_csv(os.path.join(out, "paired_cells.csv.gz"), index=False)

    agg = {"n_cells": ("cell", "size"), "n_donors": ("donor", "nunique"),
           "n_libs": ("lib", "nunique")}
    firsts = [c for c in P.columns if c not in ("cell", "clone_key", "donor", "lib")]
    clo = P.groupby("clone_key").agg(**agg, **{c: (c, "first") for c in firsts})
    clo["donors"] = P.groupby("clone_key").donor.apply(lambda s: ",".join(sorted(set(s))))
    clo = clo.reset_index()
    clo["clone_id"] = ["clone%06d" % i for i in range(len(clo))]

    Vg, Jg = imgt.germlines(cfg, "V"), imgt.germlines(cfg, "J")
    for c in chains:
        if source == "observed":
            raw = {i: observed_vdomain({k.split("|")[0]: v for k, v in r.items()
                                        if k.endswith("|" + c) or "|" not in k}, clo, c)
                   for i, r in clo.iterrows()}
            num = anarci_span(pd.Series(raw), (c,), cfg)
            clo["seq_%s" % c] = [num.get(i, ("", {}, ""))[0] for i in clo.index]
            clo["allele_%s" % c] = "observed"
            clo["jallele_%s" % c] = "observed"
        else:
            built = [reference_vdomain(Vg[c], Jg[c], r["v_%s" % c], r["j_%s" % c],
                                       r["cdr3_%s" % c]) for _, r in clo.iterrows()]
            clo["seq_%s" % c] = [b[0] for b in built]
            clo["allele_%s" % c] = [b[1] for b in built]
            clo["jallele_%s" % c] = [b[2] for b in built]
        clo["len_%s" % c] = clo["seq_%s" % c].str.len()
        clo["cdr3len_%s" % c] = clo["cdr3_%s" % c].astype(str).str.len()
    clo["vdomain_source"] = source
    clo.to_csv(os.path.join(out, "vdomains.csv.gz"), index=False)

    ok = np.ones(len(clo), bool)
    for c in chains:
        s = clo["seq_%s" % c].fillna("")
        anchor = clo["cdr3_%s" % c].astype(str)
        ok &= (s.str.len() >= 100) & (s.str.len() <= 145) & ~s.str.contains(r"[X*]")
        ok &= anchor.str.match(r"^C.*[FWC]$").fillna(False)
    fs = clo[ok].copy()
    stages.append(dict(step="clonotypes", rows=len(clo)))
    stages.append(dict(step="folding set (both V domains valid)", rows=len(fs)))
    fs.to_csv(os.path.join(out, "folding_set.csv.gz"), index=False)
    pd.DataFrame(stages).to_csv(os.path.join(out, "vdomain_build_qc.csv"), index=False)

    print("clonotypes %d | folding set %d | source %s"
          % (len(clo), len(fs), source), flush=True)
    print("constant-region calls recorded: %s (never folded)"
          % ", ".join("%s=%d" % (c, int(clo["c_%s" % c].notna().sum())) for c in chains))
    if source == "reference":
        print("NOTE: every framework and germline-loop residue in these models is a "
              "reference allele, not the donor's. Run vdomain_provenance.py next.")
    return fs


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.json")
