#!/bin/bash
# Stage B2c - fold the remaining batches unattended, with a cooldown between each.
#
# Each batch's summary is appended to PROGRESS.md as soon as that batch finishes,
# so the running log stays incremental rather than being reconstructed at the end.
#
# Usage: bash code/b2c_run_remaining.sh <cooldown_seconds> <target> [target ...]
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"   # package: package folder of this script
cd "${VOXEL_UP_WRITE_ROOT:?set VOXEL_UP_WRITE_ROOT to the run folder}" || exit 1
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$PWD/logs/numba_cache}"
export PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}"

COOL=$1; shift
TARGETS=("$@")
N=${#TARGETS[@]}

for i in "${!TARGETS[@]}"; do
  T=${TARGETS[$i]}
  LOG=logs/fold_batch_${T}.log
  python "$HERE/b2_fold_batch.py" "$T" > "$LOG" 2>&1
  RC=$?
  TS=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  {
    echo
    echo "## $TS | Stage B2 | batch to $T complete (exit $RC)"
    grep -vE '^(shard [0-9]+: [0-9]+/|ImmuneBuilder patches)' "$LOG" \
      | sed -n '/=== BATCH SUMMARY/,$p'
    RECOV=$(grep -c 'recovered' "$LOG" 2>/dev/null || true)
    if [ "$RECOV" != "0" ]; then
      echo "recovery sweep output:"
      grep -E '(recovered|dropped|no failed models)' "$LOG" | sed 's/^/  /'
    fi
  } >> PROGRESS.md
  echo "batch to $T done (exit $RC)"
  if [ "$i" -lt "$((N-1))" ]; then
    echo "cooldown ${COOL}s"
    sleep "$COOL"
  fi
done
echo "ALL BATCHES DONE: $(ls structures/*.pdb 2>/dev/null | wc -l | tr -d ' ') models on disk"
