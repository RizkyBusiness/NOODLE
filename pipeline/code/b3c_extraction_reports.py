"""What Stage B3 extraction implies for clustering, before any threshold is chosen.

Four things, all computed on this repertoire rather than assumed:

  1. loop-length classes. The site descriptor is position-matched, so RMSD is only
     defined between receptors of the same CDR1/2/3 length class. A receptor alone
     in its class can never join a cluster, whatever the threshold - so the
     singleton fraction is a hard ceiling on clustering coverage.
  2. germline sharing among length-comparable pairs. The germline-free control
     discards pairs sharing V genes; its power is set by how many pairs survive,
     and mouse TRBV diversity here is low, so this is measured up front.
  3. the receptors whose site could not be extracted - which loop was missing.
  4. a positive control that does not depend on the deferred invariant set: pairs
     of clonotypes with identical V and J on both chains whose CDR3s differ by a
     single amino acid. These must come out close; if they do not, the descriptor
     is not measuring what it claims. Their RMSDs are reported against a
     within-class background so the comparison is a distribution, not one number.

Writes: tables/B3c_length_classes.csv, B3c_germline_chance_rates.csv,
        B3c_incomplete_sites.csv, B3c_near_identical_pairs.csv.gz,
        B3c_extraction_summary.csv
"""
import itertools
import os
import sys

import numpy as np
import pandas as pd

import imgt  # noqa: E402


def site_rmsd(D, i, j):
    d = D[i] - D[j]
    return float(np.sqrt((d ** 2).sum(1).mean()))


def main(config_path="config_main.json"):
    cfg = imgt.load(config_path)
    out = cfg["out_dir"]
    rng = np.random.default_rng(cfg["controls"]["seed"])
    z = np.load(os.path.join(out, "combining_sites.npz"), allow_pickle=True)
    ids, D, LC = list(z["clone_id"]), z["descriptor"], list(z["length_class"])
    row = {c: i for i, c in enumerate(ids)}
    qc = pd.read_csv(os.path.join(out, "site_extraction_qc.csv"))
    fs = pd.read_csv(os.path.join(out, "folding_set.csv.gz"), low_memory=False)
    fs = fs.set_index("clone_id")
    summary = []

    # ---- 1. length classes -------------------------------------------------
    S = pd.Series(LC, index=ids)
    n = S.value_counts()
    L = pd.DataFrame(dict(length_class=n.index, n_receptors=n.values))
    L["singleton"] = L.n_receptors == 1
    L["pairs_within_class"] = L.n_receptors * (L.n_receptors - 1) // 2
    L = L.sort_values("n_receptors", ascending=False)
    L.to_csv(os.path.join(out, "B3c_length_classes.csv"), index=False)
    n_single = int(L.loc[L.singleton, "n_receptors"].sum())
    tot = len(ids)
    summary += [
        dict(quantity="receptors with a site descriptor", value=tot),
        dict(quantity="loop-length classes", value=len(L)),
        dict(quantity="largest class size", value=int(L.n_receptors.max())),
        dict(quantity="median class size", value=float(L.n_receptors.median())),
        dict(quantity="receptors alone in their length class", value=n_single),
        dict(quantity="ceiling on clustering coverage (%)",
             value=round(100 * (tot - n_single) / tot, 2)),
        dict(quantity="length-comparable pairs in total",
             value=int(L.pairs_within_class.sum())),
    ]

    # ---- 2. germline sharing among length-comparable pairs -----------------
    va = fs["v_A"].astype(str).str.split("*").str[0]
    vb = fs["v_B"].astype(str).str.split("*").str[0]
    df = pd.DataFrame(dict(lc=S, va=va.reindex(S.index), vb=vb.reindex(S.index)))
    pairs = shareA = shareB = shareAB = 0
    for _lc, g in df.groupby("lc"):
        k = len(g)
        if k < 2:
            continue
        pairs += k * (k - 1) // 2
        c = lambda s: int((s * (s - 1) // 2).sum())  # noqa: E731
        shareA += c(g.va.value_counts())
        shareB += c(g.vb.value_counts())
        shareAB += c(g.groupby(["va", "vb"]).size())
    G = pd.DataFrame([
        dict(relation="share TRAV", pairs=shareA, pct=round(100 * shareA / pairs, 3)),
        dict(relation="share TRBV", pairs=shareB, pct=round(100 * shareB / pairs, 3)),
        dict(relation="share both", pairs=shareAB, pct=round(100 * shareAB / pairs, 3)),
        dict(relation="share neither (retained by the germline-free control)",
             pairs=pairs - shareA - shareB + shareAB,
             pct=round(100 * (pairs - shareA - shareB + shareAB) / pairs, 3)),
        dict(relation="all length-comparable pairs", pairs=pairs, pct=100.0)])
    G.to_csv(os.path.join(out, "B3c_germline_chance_rates.csv"), index=False)
    summary += [dict(quantity="observed TRAV genes", value=int(df.va.nunique())),
                dict(quantity="observed TRBV genes", value=int(df.vb.nunique())),
                dict(quantity="length-comparable pairs sharing both V genes (%)",
                     value=float(G.loc[G.relation == "share both", "pct"].iloc[0]))]

    # ---- 3. the incomplete sites -------------------------------------------
    bad = qc[qc.status != "ok"].copy()
    loops = [c for c in qc.columns if c[:1] in ("A", "B") and "_" in c]
    bad["missing_loops"] = bad[loops].isna().apply(
        lambda r: ",".join([c for c in loops if pd.isna(r[c])]), axis=1)
    bad["zero_loops"] = bad[loops].apply(
        lambda r: ",".join([c for c in loops if r[c] == 0]), axis=1)
    keep = ["clone_id", "status", "fr_rmsd", "site_sasa", "missing_loops", "zero_loops"]
    bad[keep + loops].to_csv(os.path.join(out, "B3c_incomplete_sites.csv"), index=False)
    summary += [dict(quantity="receptors with an incomplete site", value=len(bad))]
    print("=== incomplete sites: %d ===" % len(bad))
    if len(bad):
        pat = (bad.missing_loops.replace("", np.nan).fillna(bad.zero_loops)
               .replace("", "none-flagged"))
        print(pat.value_counts().to_string())

    # ---- 4. near-identical CDR3 positive control ---------------------------
    key = ["v_A", "j_A", "v_B", "j_B"]
    F = fs.loc[[i for i in ids if i in fs.index]].copy()
    for c in key:
        F[c] = F[c].astype(str).str.split("*").str[0]
    F["la"], F["lb"] = F.cdr3_A.str.len(), F.cdr3_B.str.len()
    rows = []
    for _, g in F.groupby(key + ["la", "lb"]):
        if len(g) < 2:
            continue
        idx = list(g.index)
        for a, b in itertools.combinations(idx, 2):
            ra, rb = g.loc[a], g.loc[b]
            da = sum(x != y for x, y in zip(ra.cdr3_A, rb.cdr3_A))
            db = sum(x != y for x, y in zip(ra.cdr3_B, rb.cdr3_B))
            if da + db != 1:
                continue
            if S.get(a) != S.get(b):
                continue
            rows.append(dict(clone_a=a, clone_b=b, length_class=S.get(a),
                             chain_differing=("A" if da else "B"),
                             v_A=ra.v_A, j_A=ra.j_A, v_B=ra.v_B, j_B=ra.j_B,
                             cdr3_A_a=ra.cdr3_A, cdr3_A_b=rb.cdr3_A,
                             cdr3_B_a=ra.cdr3_B, cdr3_B_b=rb.cdr3_B,
                             site_rmsd=site_rmsd(D, row[a], row[b])))
    N = pd.DataFrame(rows)
    N.to_csv(os.path.join(out, "B3c_near_identical_pairs.csv.gz"), index=False)

    # background: random length-comparable pairs, same classes as the control pairs
    bg = []
    if len(N):
        want = N.length_class.value_counts() * 20
        for lc, k in want.items():
            pool = [i for i in S.index[S == lc]]
            if len(pool) < 2:
                continue
            for _ in range(int(k)):
                a, b = rng.choice(len(pool), 2, replace=False)
                bg.append(site_rmsd(D, row[pool[a]], row[pool[b]]))
    bg = np.array(bg)
    np.save(os.path.join(out, "B3c_background_rmsd.npy"), bg)
    if len(N):
        summary += [
            dict(quantity="near-identical pairs (1 aa in CDR3, identical V and J)",
                 value=len(N)),
            dict(quantity="distinct clonotypes involved",
                 value=int(pd.unique(N[["clone_a", "clone_b"]].values.ravel()).size)),
            dict(quantity="near-identical pair site RMSD, median (A)",
                 value=round(float(N.site_rmsd.median()), 3)),
            dict(quantity="near-identical pair site RMSD, 90th pct (A)",
                 value=round(float(N.site_rmsd.quantile(0.9)), 3)),
            dict(quantity="near-identical pair site RMSD, max (A)",
                 value=round(float(N.site_rmsd.max()), 3)),
            dict(quantity="background within-class RMSD, median (A)",
                 value=round(float(np.median(bg)), 3)),
            dict(quantity="background within-class RMSD, 5th pct (A)",
                 value=round(float(np.percentile(bg, 5)), 3))]
    summary += [dict(quantity="framework fit RMSD, median (A)",
                     value=round(float(np.median(z["fr_rmsd"])), 3)),
                dict(quantity="framework fit RMSD, 90th pct (A)",
                     value=round(float(np.percentile(z["fr_rmsd"], 90)), 3))]
    Sm = pd.DataFrame(summary)
    Sm.to_csv(os.path.join(out, "B3c_extraction_summary.csv"), index=False)
    print("\n=== length classes ===")
    print(L.head(8).to_string(index=False))
    print("\n=== germline sharing among length-comparable pairs ===")
    print(G.to_string(index=False))
    print("\n=== summary ===")
    print(Sm.to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config_main.json")
