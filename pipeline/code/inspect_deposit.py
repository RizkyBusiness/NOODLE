"""STAGE 1a - inspect the V(D)J deposit and decide where the V domain comes from.

This is the branch that decides what every later structural number means, so it
runs before anything else and it never stays silent. It reports the column set
of each V(D)J table, whether a contig FASTA is present, the pipeline version if
the deposit records one, and then chooses:

  observed   the deposit carries the donor's own full-length V region, either as
             fwr1-fwr4 + cdr1-cdr3 columns, or as an aligned/whole contig
             sequence, or as a contig FASTA -> FOLD THE DONOR'S OWN SEQUENCE
  reference  the deposit carries only gene calls plus CDR3 (Cell Ranger 5.x and
             earlier, and many AIRR exports) -> fall back to reference-allele
             reconstruction, and say so everywhere downstream

Handles plain files, gzip, and tar archives - GEO supplementary downloads are
often tars whose members are the real tables, and reading a tar as a table
silently returns padding bytes rather than a header.

Usage:  python inspect_deposit.py [config.json]
Writes: <out_dir>/deposit_inspection.json, <out_dir>/deposit_inspection.md
"""
import gzip
import io
import json
import os
import re
import sys
import tarfile

import imgt

REGION_RE = re.compile(r"^(fwr[1-4]|cdr[1-3])(_nt|_aa)?$", re.I)
FULLSEQ_COLS = ("sequence_alignment_aa", "sequence_alignment", "sequence_aa",
                "sequence", "sequence_nt")
GENE_COLS = ("v_call", "v_gene", "j_call", "j_gene", "d_call", "d_gene",
             "c_call", "c_gene", "junction", "junction_aa", "cdr3", "cdr3_nt")
TABLE_RE = re.compile(r"(contig|airr|rearrangement|clonotype|vdj).*\.(csv|tsv|txt)(\.gz)?$", re.I)
FASTA_RE = re.compile(r"(contig|consensus|vdj).*\.(fa|fasta)(\.gz)?$", re.I)
VERSION_RE = re.compile(r"(cellranger|cell ranger|trust4|mixcr|immcantation|"
                        r"changeo|dandelion)[- _]?v?([0-9]+\.[0-9]+(\.[0-9]+)?)", re.I)


def _open_bytes(path):
    return gzip.open(path, "rb") if path.endswith(".gz") else open(path, "rb")


def walk_deposit(root):
    """Yield (label, kind, reader) for every file and tar member under root.

    reader() returns a binary file object. kind is 'table', 'fasta', 'matrix',
    'json' or 'other', decided on the member name alone.
    """
    def kind_of(name):
        b = os.path.basename(name)
        if TABLE_RE.search(b):
            return "table"
        if FASTA_RE.search(b):
            return "fasta"
        if b.startswith("matrix.mtx") or b.endswith("metrics_summary.csv"):
            return "matrix"
        if b.endswith(".json"):
            return "json"
        return "other"

    for dirpath, _, files in os.walk(root):
        for f in sorted(files):
            p = os.path.join(dirpath, f)
            rel = os.path.relpath(p, root)
            if tarfile.is_tarfile(p):
                with tarfile.open(p) as tf:
                    members = [m for m in tf.getmembers() if m.isfile()]
                for m in members:
                    def rd(path=p, name=m.name):
                        tf2 = tarfile.open(path)
                        fh = tf2.extractfile(name)
                        raw = fh.read()
                        tf2.close()
                        if name.endswith(".gz"):
                            raw = gzip.decompress(raw)
                        return io.BytesIO(raw)
                    yield "%s::%s" % (rel, m.name), kind_of(m.name), rd
                continue
            yield rel, kind_of(f), (lambda path=p: _open_bytes(path))


def header_of(reader, sep=None):
    with reader() as fh:
        first = fh.readline().decode("utf-8", "replace").rstrip("\r\n")
    if sep is None:
        sep = "\t" if first.count("\t") > first.count(",") else ","
    return [c.strip().strip('"') for c in first.split(sep)], sep


def classify(cols):
    low = [c.lower() for c in cols]
    regions = sorted({c for c in low if REGION_RE.match(c)})
    base = {re.sub(r"_(nt|aa)$", "", r) for r in regions}
    have_all_regions = {"fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2"} <= base
    fullseq = [c for c in FULLSEQ_COLS if c in low]
    return dict(n_columns=len(cols), columns=list(cols),
                region_columns=regions, gene_columns=[c for c in GENE_COLS if c in low],
                full_sequence_columns=fullseq,
                has_all_framework_and_loop_columns=bool(have_all_regions),
                has_constant_call=any(c in low for c in ("c_call", "c_gene")))


def main(config_path="config.json"):
    cfg = imgt.load(config_path)
    root, out = cfg["deposit_dir"], cfg["out_dir"]
    if not os.path.isdir(root):
        raise SystemExit("deposit_dir %r does not exist" % root)

    tables, fastas, versions, others = [], [], [], []
    for label, kind, reader in walk_deposit(root):
        if kind == "table":
            try:
                cols, sep = header_of(reader)
            except Exception as e:
                others.append(dict(file=label, error="%s: %s" % (type(e).__name__, e)))
                continue
            rec = dict(file=label, separator=sep)
            rec.update(classify(cols))
            tables.append(rec)
        elif kind == "fasta":
            with reader() as fh:
                head = fh.readline().decode("utf-8", "replace").strip()
            fastas.append(dict(file=label, first_record=head[:120]))
        elif kind in ("matrix", "json"):
            with reader() as fh:
                blob = fh.read(20000).decode("utf-8", "replace")
            for m in VERSION_RE.finditer(blob):
                versions.append(dict(file=label, software=m.group(1),
                                     version=m.group(2)))

    best = None
    for t in tables:
        score = (t["has_all_framework_and_loop_columns"] * 4
                 + bool(t["full_sequence_columns"]) * 2 + len(t["gene_columns"]) / 100)
        if best is None or score > best[0]:
            best = (score, t)
    primary = best[1] if best else None

    if primary is None:
        source, why = "none", "no V(D)J table found under %s" % root
    elif primary["has_all_framework_and_loop_columns"]:
        source, why = "observed", "framework and germline-loop columns present in %s" % primary["file"]
    elif primary["full_sequence_columns"]:
        source, why = "observed", "whole/aligned contig sequence present in %s (%s)" % (
            primary["file"], ", ".join(primary["full_sequence_columns"]))
    elif fastas:
        source, why = "observed", "contig FASTA present (%s)" % fastas[0]["file"]
    else:
        source, why = "reference", (
            "%s carries gene calls plus CDR3 only, and no contig FASTA is present, "
            "so the donor's own framework and germline loops are not in the deposit"
            % primary["file"])

    rep = dict(dataset=cfg["dataset"], species=cfg["species"],
               locus_pair=cfg["locus_pair"], chains=cfg["chains"],
               deposit_dir=root, vdomain_source=source, decision_reason=why,
               pipeline_versions=versions, tables=tables, fasta_files=fastas,
               unreadable=others, primary_table=primary["file"] if primary else None)
    json.dump(rep, open(os.path.join(out, "deposit_inspection.json"), "w"), indent=1)

    ver = "; ".join(sorted({"%s %s" % (v["software"], v["version"]) for v in versions})) or "not recorded"
    md = ["# Deposit inspection - %s" % cfg["dataset"], "",
          "**V-domain source: %s**" % source.upper(), "", why, "",
          "- pipeline version: %s" % ver,
          "- V(D)J tables found: %d" % len(tables),
          "- contig FASTA: %s" % (fastas[0]["file"] if fastas else "ABSENT"),
          "- constant-region call in primary table: %s"
          % ("yes" if primary and primary["has_constant_call"] else "no"), ""]
    if source == "reference":
        md += ["> The folded V domains will NOT be the donors' own sequences. Every",
               "> framework and germline-loop residue comes from a reference allele,",
               "> and roughly half of the compared combining-site positions are",
               "> germline-encoded. Run vdomain_provenance.py and read its report",
               "> before interpreting any structural distance.", ""]
    for t in tables:
        md += ["## %s" % t["file"], "",
               "%d columns, separator %r" % (t["n_columns"], t["separator"]), "",
               "```", ", ".join(t["columns"]), "```", "",
               "- region columns: %s" % (", ".join(t["region_columns"]) or "none"),
               "- whole-sequence columns: %s" % (", ".join(t["full_sequence_columns"]) or "none"),
               "- gene/junction columns: %s" % (", ".join(t["gene_columns"]) or "none"), ""]
    open(os.path.join(out, "deposit_inspection.md"), "w").write("\n".join(md))

    print("V-DOMAIN SOURCE: %s" % source.upper())
    print("  reason: %s" % why)
    print("  pipeline version: %s" % ver)
    print("  tables %d | contig FASTA %s | wrote %s/deposit_inspection.{json,md}"
          % (len(tables), "yes" if fastas else "no", out))
    return rep


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.json")
