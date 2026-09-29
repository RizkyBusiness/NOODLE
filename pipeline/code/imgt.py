"""Run configuration, IMGT conventions and germline access.

Every script in this pipeline reads its species, its two chain labels and its
paths from one JSON config, so the same code runs on an alpha/beta or a
gamma/delta repertoire, human or mouse, without edits.

IMGT unique numbering is shared by every IG and TR V domain, so the position
ranges below are locus- and species-independent. HV4 (CDR4) sits inside FR3 and
contacts peptide-MHC in many solved complexes, which is why it is part of the
combining site here and not part of the framework.
"""
import json
import os

CDR_RANGES = {"CDR1": range(27, 39), "CDR2": range(56, 66), "CDR3": range(105, 118)}
HV4_RANGE = range(81, 87)
SITE_RANGES = dict(CDR_RANGES, HV4=HV4_RANGE)
FR_POS = [p for p in range(1, 129) if not any(p in r for r in CDR_RANGES.values())]
SITE_NAMES = ("CDR1", "CDR2", "CDR3", "HV4")

# IMGT anchors: 104 is the conserved second cysteine, 118 the conserved F/W.
# junction_aa runs 104..118 inclusive, so a reconstructed V domain is
# germline V positions 1..103 + junction + germline J positions 119..128.
CYS_104, ANCHOR_118 = 104, 118

LOCUS_PAIRS = {"AB": ("A", "B"), "GD": ("G", "D")}
# TRAV/DV genes are shared between the TRA and TRD loci, so ANARCI assigns
# those alpha V domains to chain type D (and vice versa). Any chain-type check
# must accept the partner locus or it rejects genuine chains.
SHARED_CHAIN_TYPES = {"A": ("A", "D"), "D": ("D", "A"), "B": ("B",), "G": ("G",)}

DEFAULTS = {
    "dataset": "dataset",
    "species": "human",              # any species in anarci.germlines
    "locus_pair": "AB",              # "AB" or "GD"
    "mhc_class": "unspecified",      # I | II | unspecified - interpretation only
    "deposit_dir": "raw",
    "out_dir": "out",
    "cell_col": "cell",
    "donor_col": "donor",
    "batch_col": "lib",
    "vdj": {"require_productive": True, "require_full_length": True,
            "require_high_confidence": True, "min_umis": 1,
            "max_chains_per_cell": 1},
    "folding": {"builder": "TCRBuilder2", "weights_dir": "weights/tcr",
                "n_shards": 8, "threads": 2, "model_dir": "structures"},
    "clustering": {"linkage": "complete", "min_cluster_size": 2, "n_boot": 20,
                   "boot_frac": 0.8, "n_resample": 10},
    "benchmark": {"percentile": 95, "max_structures": 40, "max_resolution": 3.0,
                  "min_threshold": 0.75},
    "controls": {"n_perm": 5000, "n_tiebreak": 200, "seed": 0,
                 "budgets": [30, 100, 300, 1000, 3000, 10000]},
    "transcriptome": {"n_hvg": 2000, "n_pcs": 50, "k": 15, "resolution": 0.6,
                      "perplexity": 30, "min_genes": 200, "max_pct_mito": 15.0,
                      "min_cells_per_gene": 3, "integrate": "harmony",
                      # gene-name patterns dropped from the variable-gene pool:
                      # a receptor-driven embedding would make the structure /
                      # state comparison circular
                      "exclude_gene_patterns": [r"^TR[ABGD][VDJC]", r"^IG[HKL][VDJC]"],
                      "sources": [],          # [{name, path, kind, lib_map}]
                      "reference_labels": None},
}


def _merge(base, over):
    out = dict(base)
    for k, v in (over or {}).items():
        out[k] = _merge(base[k], v) if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return out


def load(path="config.json"):
    """Read a run config, filling in defaults; validates the locus pair."""
    user = json.load(open(path)) if os.path.exists(path) else {}
    cfg = _merge(DEFAULTS, user)
    if cfg["locus_pair"] not in LOCUS_PAIRS:
        raise ValueError("locus_pair must be one of %s" % list(LOCUS_PAIRS))
    cfg["chains"] = list(LOCUS_PAIRS[cfg["locus_pair"]])
    cfg["config_path"] = path
    os.makedirs(cfg["out_dir"], exist_ok=True)
    return cfg


def chain_labels(cfg):
    """The two model chain ids, first chain first (A,B for alpha/beta)."""
    return tuple(cfg["chains"])


def chain_prefix(cfg, chain):
    """Gene-name prefix for a chain: TRA / TRB / TRG / TRD."""
    return "TR" + chain


def site_loops(cfg):
    """(chain, loop) pairs making up the combining site, in a fixed order."""
    return [(c, n) for c in chain_labels(cfg) for n in SITE_NAMES]


def cdr_loops(cfg):
    """(chain, loop) pairs used for the length-matched geometry descriptor."""
    return [(c, n) for c in chain_labels(cfg) for n in ("CDR1", "CDR2", "CDR3")]


def germlines(cfg, segment="V"):
    """{chain: {gene: {allele: {imgt_position: residue}}}} from ANARCI's tables.

    ANARCI ships IMGT-gapped germline sequences, so the index of a non-gap
    character is its IMGT position - which is what lets a reconstruction splice
    at exact positions instead of guessing offsets.
    """
    from anarci.germlines import all_germlines as G
    sp, out = cfg["species"], {}
    for ch in chain_labels(cfg):
        d = {}
        for ct in SHARED_CHAIN_TYPES.get(ch, (ch,)):
            for allele, seq in G[segment].get(ct, {}).get(sp, {}).items():
                gene = allele.split("*")[0]
                d.setdefault(gene, {})[allele] = {
                    i + 1: c for i, c in enumerate(seq) if c != "-"}
        out[ch] = d
    if not any(out.values()):
        raise ValueError("no %s germlines for species %r chains %s"
                         % (segment, sp, chain_labels(cfg)))
    return out


def pick_allele(alleles, called=None):
    """Choose a germline allele: the exact call if present, else the lowest-numbered.

    Returns (allele_name, positions). Deposits that report bare gene names force
    the fallback, which is the substitution the provenance report has to quantify.
    """
    if called and called in alleles:
        return called, alleles[called]
    name = sorted(alleles)[0]
    return name, alleles[name]
