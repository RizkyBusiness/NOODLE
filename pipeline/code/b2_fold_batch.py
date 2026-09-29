"""Stage B2 - fold the repertoire in fixed-size batches.

fold_receptors.py has no batch limit: it folds every row of
<out_dir>/folding_set.csv.gz. It is however resumable by construction - a
receptor counts as done only when its PDB exists, and the per-shard QC rows are
concatenated with drop_duplicates(clone_id, keep="last").

So a batch is just a cumulative prefix of the shuffled full set written over
folding_set.csv.gz: shards skip everything already on disk and fold only the new
rows. After the final batch, folding_set.csv.gz is the complete set, which is
what every downstream stage expects to read.

The order is a seed-0 shuffle of the full set (folding_set.full.csv.gz), so any
interim batch is an unbiased sample of the repertoire rather than a lane-ordered
prefix.

Usage:  python code/b2_fold_batch.py <cumulative_target> [config]
"""
import os
import subprocess
import sys
import time

import pandas as pd

import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
SKILL = paths.LIB  # package: vendored earlier-skill modules (this folder)


def main():
    target = int(sys.argv[1])
    cfgp = sys.argv[2] if len(sys.argv) > 2 else "config_main.json"
    sys.path.insert(0, SKILL)
    import imgt
    cfg = imgt.load(os.path.join(ROOT, cfgp))
    out, mdir = cfg["out_dir"], cfg["folding"]["model_dir"]

    full = pd.read_csv(os.path.join(out, "folding_set.full.csv.gz"), low_memory=False)
    target = min(target, len(full))
    full.iloc[:target].to_csv(os.path.join(out, "folding_set.csv.gz"), index=False)

    before = len([p for p in os.listdir(mdir) if p.endswith(".pdb")]) \
        if os.path.isdir(mdir) else 0
    t0 = time.time()
    rc = subprocess.call([sys.executable, os.path.join(SKILL, "fold_receptors.py"),
                          os.path.join(ROOT, cfgp)], cwd=ROOT)
    dt = time.time() - t0
    after = len([p for p in os.listdir(mdir) if p.endswith(".pdb")])

    # A refinement that raises leaves the UNREFINED coordinates at the final path
    # and extract_sites.py trusts file existence, so every batch is swept for
    # failures before it is called done. See code/b2b_recover_failed.py.
    rc2 = subprocess.call([sys.executable, os.path.join(paths.SCRIPTS, "b2b_recover_failed.py"),
                           cfgp], cwd=ROOT)

    qc = pd.read_csv(os.path.join(out, "model_qc.csv"))
    ok = qc[qc.status.str.startswith("ok")]
    print("\n=== BATCH SUMMARY (target %d of %d) ===" % (target, len(full)))
    print("exit %d | wall %.1f min | new models %d | total models %d"
          % (rc, dt / 60, after - before, after))
    print("qc rows %d | ok %d (of which recovered %d) | failed %d"
          % (len(qc), len(ok), (qc.status == "ok_recovered").sum(),
             (~qc.status.str.startswith("ok")).sum()))
    if after > before:
        print("rate %.2f s/receptor aggregate" % (dt / (after - before)))
    if len(ok):
        print("model error estimate (A): median %.2f  p90 %.2f  max %.2f"
              % (ok.rmsd_error_mean.median(), ok.rmsd_error_mean.quantile(0.9),
                 ok.rmsd_error_mean.max()))
        print("atoms per model: median %d  min %d" % (ok.n_atoms.median(), ok.n_atoms.min()))
    fails = qc[~qc.status.str.startswith("ok")]
    if len(fails):
        print("failure reasons:")
        print(fails.message.str.split(":").str[0].value_counts().to_string())
    print("remaining to fold: %d" % (len(full) - after))
    # a failed model that was neither recovered nor dropped (dropped = moved to structures_failed/ and removed from the
    # folding set by b2b) is an error, as is a crashed fold or recovery step
    keep = set(pd.read_csv(os.path.join(out, "folding_set.full.csv.gz"), usecols=["clone_id"]).clone_id)
    left = sorted(set(fails.clone_id) & keep)
    if rc or rc2 or left:
        sys.exit("b2_fold_batch: fold exit %d, recovery exit %d, %d failed model(s) left%s"
                 % (rc, rc2, len(left), (": " + ", ".join(left[:10])) if left else ""))


if __name__ == "__main__":
    main()
