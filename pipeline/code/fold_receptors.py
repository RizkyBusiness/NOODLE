"""STAGE 2a - fold the paired V modules, sharded and resumable.

One process per shard, each writing PDBs into <model_dir> and its own QC row per
receptor. Only a written PDB counts as done, so an interrupted run resumes and
failures are retried rather than silently skipped.

TCRBuilder2 builds alpha/beta V modules only - there is no gamma/delta builder in
ImmuneBuilder 1.x, and the constant domains are never modelled. For a
gamma/delta repertoire, set folding.builder to "external", produce the paired
models with a general-purpose predictor, and drop IMGT-numbered two-chain PDBs
named <clone_id>.pdb into <model_dir>; every later stage is locus-agnostic and
picks up from there.

Usage:  python fold_receptors.py [config.json]              # dispatch all shards
        python fold_receptors.py [config.json] <i> <n>      # one shard directly
Writes: <model_dir>/<clone_id>.pdb, <out_dir>/model_qc.csv
"""
import os
import subprocess
import sys
import time

import pandas as pd

import imgt


def run_shard(cfg, shard, nshards):
    import csv
    import torch
    f = cfg["folding"]
    torch.set_num_threads(f["threads"])
    os.environ["OMP_NUM_THREADS"] = str(f["threads"])
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import ib_patches
    print("ImmuneBuilder patches: %s" % ib_patches.apply(), flush=True)
    from ImmuneBuilder import TCRBuilder2

    chains = imgt.chain_labels(cfg)
    if set(chains) != {"A", "B"}:
        raise SystemExit("TCRBuilder2 builds alpha/beta only; set "
                         "folding.builder='external' for %s/%s" % chains)
    out = f["model_dir"]
    os.makedirs(out, exist_ok=True)
    df = pd.read_csv(os.path.join(cfg["out_dir"], "folding_set.csv.gz"), low_memory=False)
    df = df.iloc[shard::nshards]
    pred = TCRBuilder2(weights_dir=f["weights_dir"])
    qc = os.path.join(cfg["out_dir"], "model_qc_shard%d.csv" % shard)
    new = not os.path.exists(qc)
    fh = open(qc, "a", newline="")
    w = csv.writer(fh)
    if new:
        w.writerow(["clone_id", "status", "seconds", "rmsd_error_mean",
                    "rmsd_error_max", "n_atoms", "message"])
    for k, r in enumerate(df.itertuples()):
        path = os.path.join(out, "%s.pdb" % r.clone_id)
        if os.path.exists(path):
            continue
        t0 = time.time()
        try:
            m = pred.predict({c: getattr(r, "seq_%s" % c) for c in chains})
            m.save(path, n_threads=f["threads"])
            err = m.error_estimates.mean(0).sqrt().cpu().numpy()
            w.writerow([r.clone_id, "ok", round(time.time() - t0, 2),
                        round(float(err.mean()), 4), round(float(err.max()), 4),
                        sum(1 for L in open(path) if L.startswith("ATOM")), ""])
        except Exception as e:
            w.writerow([r.clone_id, "fail", round(time.time() - t0, 2), "", "", "",
                        "%s: %s" % (type(e).__name__, str(e)[:160])])
        fh.flush()
        if k % 25 == 0:
            print("shard %d: %d/%d" % (shard, k + 1, len(df)), flush=True)
    fh.close()
    print("shard %d complete" % shard, flush=True)


def main(config_path="config.json"):
    cfg = imgt.load(config_path)
    f = cfg["folding"]
    if f["builder"] == "external":
        n = len([p for p in os.listdir(f["model_dir"])
                 if p.endswith(".pdb")]) if os.path.isdir(f["model_dir"]) else 0
        print("builder=external: expecting IMGT-numbered PDBs in %s (%d present)"
              % (f["model_dir"], n))
        return
    me = os.path.abspath(__file__)
    procs = [subprocess.Popen([sys.executable, me, config_path, str(i), str(f["n_shards"])])
             for i in range(f["n_shards"])]
    codes = [p.wait() for p in procs]
    parts = [pd.read_csv(os.path.join(cfg["out_dir"], "model_qc_shard%d.csv" % i))
             for i in range(f["n_shards"])
             if os.path.exists(os.path.join(cfg["out_dir"], "model_qc_shard%d.csv" % i))]
    if parts:
        qc = pd.concat(parts).drop_duplicates("clone_id", keep="last")
        qc.to_csv(os.path.join(cfg["out_dir"], "model_qc.csv"), index=False)
        print("models ok %d | failed %d | median error estimate %.2f A"
              % ((qc.status == "ok").sum(), (qc.status != "ok").sum(),
                 qc.rmsd_error_mean.median()))
    print("shard exit codes: %s" % codes)
    if any(codes):
        raise SystemExit("fold_receptors: %d of %d shards failed" % (sum(1 for c in codes if c), len(codes)))


if __name__ == "__main__":
    a = sys.argv[1:]
    if len(a) >= 3:
        run_shard(imgt.load(a[0]), int(a[1]), int(a[2]))
    else:
        main(a[0] if a else "config.json")
