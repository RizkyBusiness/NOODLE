"""STAGE 2b - derive the clustering threshold from solved structures, not from taste.

A structural clustering threshold is only meaningful next to the error of the
models being clustered, so it is measured here rather than assumed:

  1. find solved structures of the configured locus pair in the PDB;
  2. keep those where ANARCI finds a V domain of each chain type;
  3. renumber the crystal V domains to IMGT;
  4. fold the crystal's OWN sequences with the same builder used for the
     repertoire;
  5. measure combining-site C-alpha RMSD, model against crystal, in the shared
     framework frame.

The threshold is a high percentile of that distribution: below it, two receptors
cannot be told apart from model error. threshold.json also carries the median
and twice the median so the sensitivity of every downstream count to this choice
is visible, and cluster_sites.py refuses to run without this file.

Usage:  python benchmark_threshold.py [config.json]
Writes: <out_dir>/benchmark_structures.csv, <out_dir>/threshold.json,
        <out_dir>/figures/benchmark_threshold.png
"""
import json
import os
import sys
import urllib.request

import numpy as np
import pandas as pd

import imgt
import site_geometry as G

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
PDB_FILE = "https://files.rcsb.org/download/%s.pdb"
LOCUS_TEXT = {"AB": '"T cell receptor" AND alpha AND beta',
              "GD": '"T cell receptor" AND gamma AND delta'}


def _post(url, payload):
    req = urllib.request.Request(
        url, json.dumps(payload).encode(), {"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def search_entries(cfg, rows=300):
    b = cfg["benchmark"]
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "full_text",
         "parameters": {"value": LOCUS_TEXT[cfg["locus_pair"]]}},
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal", "value": b["max_resolution"]}},
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater_or_equal", "value": 2}}]},
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": rows},
                            "results_content_type": ["experimental"]}}
    return [h["identifier"] for h in _post(SEARCH, q).get("result_set", [])]


def entity_sequences(ids):
    """{entry_id: [(auth_chain_ids, one_letter_sequence), ...]} in one GraphQL call."""
    out = {}
    for i in range(0, len(ids), 80):
        chunk = ids[i:i + 80]
        q = ("{entries(entry_ids:%s){rcsb_id polymer_entities{"
             "entity_poly{pdbx_seq_one_letter_code_can}"
             "rcsb_polymer_entity_container_identifiers{auth_asym_ids}}}}"
             % json.dumps(chunk))
        for e in _post(GRAPHQL, {"query": q})["data"]["entries"] or []:
            ent = []
            for pe in e.get("polymer_entities") or []:
                s = (pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
                ch = ((pe.get("rcsb_polymer_entity_container_identifiers") or {})
                      .get("auth_asym_ids") or [])
                if s and ch:
                    ent.append((ch, s.replace("\n", "")))
            out[e["rcsb_id"]] = ent
    return out


def number_chains(seqs, cfg):
    """ANARCI over (label, sequence) pairs -> {label: (chain_type, {pos: aa}, start)}."""
    from anarci import run_anarci
    allow = sorted({t for c in imgt.chain_labels(cfg)
                    for t in imgt.SHARED_CHAIN_TYPES.get(c, (c,))})
    inp = [(k, s) for k, s in seqs if 80 <= len(s) <= 1200]
    if not inp:
        return {}
    _, numbered, details, _ = run_anarci(inp, ncpu=1, scheme="imgt",
                                         allowed_species=[cfg["species"]], allow=set(allow))
    out = {}
    for (k, _), num, det in zip(inp, numbered, details):
        if not num:
            continue
        out[k] = (det[0]["chain_type"],
                  [(p, a) for p, a in num[0][0]],
                  det[0]["query_start"])
    return out


def read_pdb_chain_residues(path):
    """{auth_chain: [(residue_key, one_letter, [atom_lines])]} in file order."""
    out, seen = {}, {}
    for L in open(path):
        if not L.startswith("ATOM") or L[16] not in (" ", "A"):
            continue
        rn = L[17:20].strip()
        if rn not in G.THREE_TO_ONE:
            continue
        ch, key = L[21], (L[22:27])
        d = out.setdefault(ch, [])
        if seen.get(ch) != key:
            d.append((key, G.THREE_TO_ONE[rn], []))
            seen[ch] = key
        d[-1][2].append(L)
    return out


def imgt_renumber(residues, numbering, query_start, label):
    """Rewrite ATOM lines of one crystal chain with IMGT numbers and a new chain id."""
    lines, i = [], query_start
    for (num, ins), aa in numbering:
        if aa == "-":
            continue
        if i >= len(residues):
            break
        key, obs_aa, atoms = residues[i]
        i += 1
        if obs_aa != aa or num > 128:
            continue
        for L in atoms:
            lines.append(L[:21] + label + "%4d%s" % (num, ins if ins != " " else " ")
                         + L[27:])
    return lines


def main(config_path="config.json"):
    cfg = imgt.load(config_path)
    out, chains = cfg["out_dir"], imgt.chain_labels(cfg)
    b, work = cfg["benchmark"], os.path.join(cfg["out_dir"], "benchmark")
    os.makedirs(work, exist_ok=True)
    types = {c: imgt.SHARED_CHAIN_TYPES.get(c, (c,)) for c in chains}

    ids = search_entries(cfg)
    print("PDB candidates: %d" % len(ids), flush=True)
    ents = entity_sequences(ids)

    picked = []
    for eid, entities in ents.items():
        seqs = [("%s|%d" % (eid, k), s) for k, (ch, s) in enumerate(entities)]
        num = number_chains(seqs, cfg)
        assign = {}
        for k, (ch_ids, s) in enumerate(entities):
            key = "%s|%d" % (eid, k)
            if key not in num:
                continue
            ct = num[key][0]
            for c in chains:
                if ct in types[c] and c not in assign:
                    assign[c] = (ch_ids[0], s, num[key])
        if len(assign) == len(chains):
            picked.append((eid, assign))
        if len(picked) >= b["max_structures"]:
            break
    print("entries with a %s/%s V-domain pair: %d"
          % (chains[0], chains[1], len(picked)), flush=True)

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import ib_patches
    ib_patches.apply()
    from ImmuneBuilder import TCRBuilder2
    pred = TCRBuilder2(weights_dir=cfg["folding"]["weights_dir"])

    rows = []
    for eid, assign in picked:
        cryst = os.path.join(work, "%s_crystal.pdb" % eid)
        model = os.path.join(work, "%s_model.pdb" % eid)
        try:
            raw = os.path.join(work, "%s.pdb" % eid)
            if not os.path.exists(raw):
                urllib.request.urlretrieve(PDB_FILE % eid, raw)
            res = read_pdb_chain_residues(raw)
            lines, vseqs = [], {}
            for c in chains:
                auth, _, (ct, numbering, qs) = assign[c]
                if auth not in res:
                    raise KeyError("auth chain %s absent from coordinates" % auth)
                lines += imgt_renumber(res[auth], numbering, qs, c)
                vseqs[c] = "".join(a for _, a in numbering if a != "-")
            open(cryst, "w").writelines(lines + ["END\n"])
            if not os.path.exists(model):
                m = pred.predict(vseqs)
                m.save(model, n_threads=cfg["folding"]["threads"])
            ca_c, ca_m = G.read_ca(cryst, chains), G.read_ca(model, chains)
            fr = G.framework_positions([ca_c, ca_m], chains, min_frac=1.0)
            rmsd, n_site = G.site_rmsd_between(ca_c, ca_m, fr, chains)
            rows.append(dict(entry=eid, site_rmsd=rmsd, n_site_positions=n_site,
                             n_framework=len(fr),
                             **{"len_%s" % c: len(vseqs[c]) for c in chains}))
            print("  %s  site RMSD %.2f A over %d positions"
                  % (eid, rmsd if rmsd else float("nan"), n_site), flush=True)
        except Exception as e:
            rows.append(dict(entry=eid, site_rmsd=None,
                             error="%s: %s" % (type(e).__name__, str(e)[:120])))
    B = pd.DataFrame(rows)
    B.to_csv(os.path.join(out, "benchmark_structures.csv"), index=False)
    ok = B[B.site_rmsd.notna()]
    if len(ok) < 5:
        raise SystemExit("only %d usable benchmark structures; cannot derive a "
                         "threshold. Inspect benchmark_structures.csv" % len(ok))

    p = b["percentile"]
    val = float(np.percentile(ok.site_rmsd, p))
    thr = max(round(val * 20) / 20.0, b["min_threshold"])
    rec = dict(threshold=thr, estimator="p%d of model-vs-crystal combining-site RMSD" % p,
               percentile=p, raw_percentile_value=round(val, 3),
               median=round(float(ok.site_rmsd.median()), 3),
               twice_median=round(2 * float(ok.site_rmsd.median()), 3),
               n_structures=int(len(ok)), n_attempted=int(len(B)),
               min_threshold_floor=b["min_threshold"],
               max_resolution=b["max_resolution"], species=cfg["species"],
               locus_pair=cfg["locus_pair"], builder=cfg["folding"]["builder"],
               entries=sorted(ok.entry))
    json.dump(rec, open(os.path.join(out, "threshold.json"), "w"), indent=1)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.6, 3.0))
    ax.hist(ok.site_rmsd, bins=max(8, len(ok) // 3), color="0.7", edgecolor="white")
    ax.axvline(thr, color="#b2182b", lw=1.6,
               label="threshold %.2f A (p%d)" % (thr, p))
    ax.axvline(ok.site_rmsd.median(), color="#2166ac", lw=1.2, ls="--",
               label="median %.2f A" % ok.site_rmsd.median())
    ax.set_xlabel("combining-site RMSD, model vs crystal (A)")
    ax.set_ylabel("solved structures")
    ax.set_title("Model error floor, n=%d %s structures" % (len(ok), cfg["locus_pair"]),
                 fontsize=9, loc="left")
    ax.legend(frameon=False, fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    fig.savefig(os.path.join(out, "figures", "benchmark_threshold.png"), dpi=200,
                facecolor="white")

    print("THRESHOLD %.2f A  (p%d of %d structures; median %.2f, 2x median %.2f)"
          % (thr, p, len(ok), rec["median"], rec["twice_median"]))
    return rec


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.json")
