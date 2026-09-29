"""Stage B2b - recover receptors whose refinement raised, and quarantine the rest.

Why this exists. TCRBuilder2.save() writes the UNREFINED coordinates to the target
path first and only then refines in place. If OpenMM raises (rather than returning
False), the exception escapes save()'s own retry loop, so:

  * an unrefined PDB is left on disk at the final path, and
  * fold_receptors.py records status="fail" for that clone, but
  * extract_sites.py discovers models by file existence alone and never reads
    model_qc.csv - so the unrefined model would enter the analysis unflagged,
    among thousands of refined ones, and
  * re-running the fold would skip it (the file exists) or rewrite the same
    unrefined file.

ImmuneBuilder predicts 4 models per receptor and ranks them. save() aborts on the
first exception instead of falling through to the next-ranked model, so a receptor
whose top-ranked model has bad local geometry is recoverable: refine the 2nd, 3rd
or 4th instead. This script tries all of them, in rank order, and writes the first
that refines cleanly. Anything that fails on all 4 is moved to
structures_failed/ and dropped from folding_set.full.csv.gz, so it is excluded
explicitly and counted rather than silently refolded on the next batch.

Usage:  python code/b2b_recover_failed.py [config]
"""
import csv
import os
import shutil
import sys

import pandas as pd

import paths; ROOT = paths.ROOT  # package: write root from paths.py (was the script's parent folder)
SKILL = paths.LIB  # package: vendored earlier-skill modules (this folder)


def main():
    cfgp = sys.argv[1] if len(sys.argv) > 1 else "config_main.json"
    sys.path.insert(0, SKILL)
    import imgt
    import ib_patches
    ib_patches.apply()
    from ImmuneBuilder import TCRBuilder2
    from ImmuneBuilder.refine import refine
    from ImmuneBuilder.util import add_errors_as_bfactors

    cfg = imgt.load(os.path.join(ROOT, cfgp))
    out, mdir = cfg["out_dir"], cfg["folding"]["model_dir"]
    chains = imgt.chain_labels(cfg)
    qc = pd.read_csv(os.path.join(out, "model_qc.csv"))
    # "ok" and "ok_recovered" are both done; anything else needs a rank retry.
    bad = qc[~qc.status.str.startswith("ok")].clone_id.tolist()
    if not bad:
        print("no failed models to recover")
        return
    fs = pd.read_csv(os.path.join(out, "folding_set.full.csv.gz"), low_memory=False)
    rows = fs.set_index("clone_id")
    pred = TCRBuilder2(weights_dir=cfg["folding"]["weights_dir"])
    quarantine = os.path.join(ROOT, "structures_failed")
    os.makedirs(quarantine, exist_ok=True)

    recovered, dropped, log = [], [], []
    for cid in bad:
        path = os.path.join(mdir, "%s.pdb" % cid)
        r = rows.loc[cid]
        m = pred.predict({c: r["seq_%s" % c] for c in chains})
        done = False
        for rank, idx in enumerate(m.ranking):
            tmp = os.path.join(mdir, "_recover_%s.pdb" % cid)
            m.save_single_unrefined(tmp, index=idx)
            try:
                refine(tmp, tmp, check_for_strained_bonds=True, n_threads=4)
            except Exception as e:
                log.append(dict(clone_id=cid, rank=rank, model_index=int(idx),
                                result="refine_raised",
                                message="%s: %s" % (type(e).__name__, str(e)[:200])))
                os.path.exists(tmp) and os.remove(tmp)
                continue
            add_errors_as_bfactors(
                tmp, m.error_estimates.mean(0).sqrt().cpu().numpy(), header=[m.header])
            os.replace(tmp, path)
            n_atoms = sum(1 for L in open(path) if L.startswith("ATOM"))
            log.append(dict(clone_id=cid, rank=rank, model_index=int(idx),
                            result="recovered", message="n_atoms=%d" % n_atoms))
            recovered.append(cid)
            qc.loc[qc.clone_id == cid, ["status", "message", "n_atoms"]] = [
                "ok_recovered", "refined from rank %d model" % rank, n_atoms]
            # fold_receptors.main() REBUILDS model_qc.csv from the per-shard CSVs
            # on every run (concat then drop_duplicates(clone_id, keep="last")),
            # so editing model_qc.csv alone is undone by the next batch and this
            # receptor gets re-recovered forever. Append the corrected row to the
            # LAST shard file, which keep="last" then honours.
            last = os.path.join(out, "model_qc_shard%d.csv"
                                % (cfg["folding"]["n_shards"] - 1))
            with open(last, "a", newline="") as fh:
                csv.writer(fh).writerow(
                    [cid, "ok_recovered", "", "", "", n_atoms,
                     "refined from rank %d model" % rank])
            done = True
            break
        if not done:
            if os.path.exists(path):
                shutil.move(path, os.path.join(quarantine, "%s.pdb" % cid))
            dropped.append(cid)
            log.append(dict(clone_id=cid, rank=-1, model_index=-1,
                            result="dropped_all_4_ranks_failed", message=""))

    if dropped:
        fs[~fs.clone_id.isin(dropped)].to_csv(
            os.path.join(out, "folding_set.full.csv.gz"), index=False)
        cur = pd.read_csv(os.path.join(out, "folding_set.csv.gz"), low_memory=False)
        cur[~cur.clone_id.isin(dropped)].to_csv(
            os.path.join(out, "folding_set.csv.gz"), index=False)
    qc.to_csv(os.path.join(out, "model_qc.csv"), index=False)
    L = pd.DataFrame(log)
    L.to_csv(os.path.join(out, "B2b_recovery_log.csv"), index=False)
    print(L.to_string(index=False))
    print("\nrecovered %d | dropped %d (moved to structures_failed/)"
          % (len(recovered), len(dropped)))


if __name__ == "__main__":
    main()
