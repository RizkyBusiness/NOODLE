"""Stage B helper: fixed random folding order (package; new script).

b2_fold_batch.py folds the repertoire in batches, each a cumulative prefix of
tables/folding_set.full.csv.gz, so every interim batch is a random sample
of the repertoire rather than a lane-ordered prefix. The project had no script
that wrote that file (it was made during the folding session). This rebuilds it:

  rows of tables/folding_set.csv.gz (as written by build_vdomains.py),
  sorted by clone_key, then shuffled with pandas sample(frac=1, random_state=0).

Sorting first makes the order independent of the row order build_vdomains
writes, so re-running this after folding gives the same file. The rule was
recovered by matching the existing file (identical rows and order for all
receptors); the original one-off command itself is not recorded.

The order affects only which receptors are folded first. The final
folding_set.csv.gz that every later stage reads is the complete set.
Run after build_vdomains.py and before b2_fold_batch.py.
"""
import pandas as pd
import paths

SEED = 0
F = pd.read_csv(paths.src("tables/folding_set.csv.gz"), low_memory=False)
assert F.clone_key.is_unique, "clone_key must be unique to define the order"
S = F.sort_values("clone_key", kind="stable").reset_index(drop=True).sample(frac=1, random_state=SEED)
S.to_csv(paths.dst("tables/folding_set.full.csv.gz"), index=False)
print("folding_set.full: %d receptors, seed %d" % (len(S), SEED))
