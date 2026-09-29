"""How well constrained is the crystal-derived error floor for MOUSE alpha-beta TCR?

benchmark_threshold.py searches the PDB with a full-text query ('"T cell receptor"
AND alpha AND beta') plus resolution <= max_resolution and >= 2 protein entities. It
applies NO species filter, NO non-redundancy filter, and stops at the first
max_structures (40) entries that happen to carry an A/B V-domain pair. The threshold
that governs all downstream clustering is a percentile of model-vs-crystal RMSD over
that set, so if the set is mostly human the floor is a human floor applied to mouse
models.

This script measures the available pool instead of assuming it:
  1. the pipeline's own query, unpaginated, so nothing is truncated at 40;
  2. a broader mouse-targeted query, to bound what the full-text query MISSES;
  3. per entry: resolution, and the source organism of the TCR chains specifically
     (not the whole entry - an entry can pair a mouse TCR with a non-mouse MHC);
  4. ANARCI typing of every polymer entity to find real A/B V-domain pairs;
  5. non-redundancy, both on the exact V-domain sequence pair and on the CDR3 pair.

Writes: tables/B3b_benchmark_structure_survey.csv       one row per usable entry
        tables/B3b_benchmark_survey_summary.csv          the reportable counts
"""
import json
import os
import sys
import urllib.request

import pandas as pd

import imgt  # noqa: E402

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
MOUSE_TAX, HUMAN_TAX = 10090, 9606


def _post(url, payload):
    req = urllib.request.Request(
        url, json.dumps(payload).encode(), {"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.load(r)


def search(nodes, rows=2000):
    q = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
         "return_type": "entry",
         "request_options": {"paginate": {"start": 0, "rows": rows},
                             "results_content_type": ["experimental"]}}
    return [h["identifier"] for h in _post(SEARCH, q).get("result_set", [])]


def res_node(maxres):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                           "operator": "less_or_equal", "value": maxres}}


def entity_info(ids):
    out = {}
    for i in range(0, len(ids), 50):
        chunk = ids[i:i + 50]
        q = ("{entries(entry_ids:%s){rcsb_id "
             "rcsb_entry_info{resolution_combined} "
             "struct{title} "
             "polymer_entities{entity_poly{pdbx_seq_one_letter_code_can} "
             "rcsb_polymer_entity_container_identifiers{auth_asym_ids} "
             "rcsb_entity_source_organism{ncbi_scientific_name ncbi_taxonomy_id}}}}"
             % json.dumps(chunk))
        for e in _post(GRAPHQL, {"query": q})["data"]["entries"] or []:
            rc = (e.get("rcsb_entry_info") or {}).get("resolution_combined") or [None]
            ents = []
            for pe in e.get("polymer_entities") or []:
                s = (pe.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
                ch = ((pe.get("rcsb_polymer_entity_container_identifiers") or {})
                      .get("auth_asym_ids") or [])
                org = pe.get("rcsb_entity_source_organism") or []
                tax = [o.get("ncbi_taxonomy_id") for o in org if o.get("ncbi_taxonomy_id")]
                nm = [o.get("ncbi_scientific_name") for o in org if o.get("ncbi_scientific_name")]
                if s and ch:
                    ents.append(dict(seq=s.replace("\n", ""), chains=ch,
                                     tax=tax[0] if tax else None,
                                     organism=nm[0] if nm else None))
            out[e["rcsb_id"]] = dict(resolution=rc[0], entities=ents,
                                     title=(e.get("struct") or {}).get("title", ""))
    return out


def type_entities(info, cfg, ncpu=4):
    """ANARCI-type every entity across all entries in ONE batched call."""
    from anarci import run_anarci
    chains = imgt.chain_labels(cfg)
    types = {c: imgt.SHARED_CHAIN_TYPES.get(c, (c,)) for c in chains}
    inp, keys = [], []
    for eid, d in info.items():
        for k, e in enumerate(d["entities"]):
            if 80 <= len(e["seq"]) <= 1200:
                inp.append(("%s|%d" % (eid, k), e["seq"]))
                keys.append((eid, k))
    allow = sorted({t for c in chains for t in types[c]})
    # unrestricted on purpose: typing must not depend on a species guess
    _, numbered, details, _ = run_anarci(inp, ncpu=ncpu, scheme="imgt",
                                         allow=set(allow), allowed_species=None)
    got = {}
    for (name, _), num, det in zip(inp, numbered, details):
        if not num:
            continue
        eid, k = name.split("|")
        vseq = "".join(a for _, a in num[0][0] if a != "-")
        nb = {p: a for p, a in num[0][0] if a != "-"}
        cdr3 = "".join(a for (p, _i), a in sorted(nb.items())
                       if p in imgt.CDR_RANGES["CDR3"])
        got[(eid, int(k))] = dict(chain_type=det[0]["chain_type"], vseq=vseq,
                                  cdr3=cdr3, hmm_species=det[0].get("species"))
    return got, types, chains


def main(config_path="config_main.json"):
    cfg = imgt.load(config_path)
    out, b = cfg["out_dir"], cfg["benchmark"]
    maxres = b["max_resolution"]
    locus_text = {"AB": '"T cell receptor" AND alpha AND beta',
                  "GD": '"T cell receptor" AND gamma AND delta'}[cfg["locus_pair"]]

    ids_pipe = search([
        {"type": "terminal", "service": "full_text", "parameters": {"value": locus_text}},
        res_node(maxres),
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "rcsb_entry_info.polymer_entity_count_protein",
                        "operator": "greater_or_equal", "value": 2}}])
    ids_mouse = search([
        {"type": "terminal", "service": "full_text", "parameters": {"value": '"T cell receptor"'}},
        res_node(maxres),
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                        "operator": "exact_match", "value": str(MOUSE_TAX)}}])
    print("pipeline query: %d entries | mouse-targeted query: %d | union %d"
          % (len(ids_pipe), len(ids_mouse), len(set(ids_pipe) | set(ids_mouse))), flush=True)

    info = entity_info(sorted(set(ids_pipe) | set(ids_mouse)))
    got, types, chains = type_entities(info, cfg)

    rows = []
    for eid, d in info.items():
        assign = {}
        for k, e in enumerate(d["entities"]):
            t = got.get((eid, k))
            if not t:
                continue
            for c in chains:
                if t["chain_type"] in types[c] and c not in assign:
                    assign[c] = (e, t)
        if len(assign) != len(chains):
            continue
        tax = {c: assign[c][0]["tax"] for c in chains}
        org = {c: assign[c][0]["organism"] for c in chains}
        allm = all(tax[c] == MOUSE_TAX for c in chains)
        allh = all(tax[c] == HUMAN_TAX for c in chains)
        rows.append(dict(
            pdb_id=eid, resolution=d["resolution"],
            tcr_species=("mouse" if allm else "human" if allh else
                         "|".join(str(org[c]) for c in chains)),
            in_pipeline_query=eid in set(ids_pipe),
            **{"tax_%s" % c: tax[c] for c in chains},
            **{"vseq_%s" % c: assign[c][1]["vseq"] for c in chains},
            **{"cdr3_%s" % c: assign[c][1]["cdr3"] for c in chains},
            title=d["title"][:90]))
    T = pd.DataFrame(rows).sort_values(["tcr_species", "resolution"])
    vcols = ["vseq_%s" % c for c in chains]
    ccols = ["cdr3_%s" % c for c in chains]
    T["redundant_vseq"] = T.duplicated(vcols, keep="first")
    T["redundant_cdr3"] = T.duplicated(ccols, keep="first")
    T.to_csv(os.path.join(out, "B3b_benchmark_structure_survey.csv"), index=False)

    srows = []
    for sp, g in T.groupby("tcr_species"):
        srows.append(dict(tcr_species=sp, entries=len(g),
                          non_redundant_vseq=int((~g.redundant_vseq).sum()),
                          non_redundant_cdr3=int((~g.redundant_cdr3).sum()),
                          best_resolution=g.resolution.min(),
                          median_resolution=round(float(g.resolution.median()), 2),
                          found_by_pipeline_query=int(g.in_pipeline_query.sum())))
    S = pd.DataFrame(srows).sort_values("entries", ascending=False)
    S.to_csv(os.path.join(out, "B3b_benchmark_survey_summary.csv"), index=False)
    print("\n=== usable %s V-domain pairs at <= %.1f A, by TCR-chain species ==="
          % ("/".join(chains), maxres))
    print(S.to_string(index=False))
    mouse = T[T.tcr_species == "mouse"]
    print("\nMOUSE entries (%d; %d non-redundant on the V-domain sequence pair):"
          % (len(mouse), int((~mouse.redundant_vseq).sum())))
    print(mouse[["pdb_id", "resolution", "redundant_vseq", "in_pipeline_query"]]
          .to_string(index=False))
    print("\nmouse PDB IDs: %s" % ", ".join(sorted(mouse.pdb_id)))
    print("non-redundant mouse PDB IDs: %s"
          % ", ".join(sorted(mouse.loc[~mouse.redundant_vseq, "pdb_id"])))
    print("\nmouse entries MISSED by the pipeline's full-text query: %d"
          % int((~mouse.in_pipeline_query).sum()))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config_main.json")
