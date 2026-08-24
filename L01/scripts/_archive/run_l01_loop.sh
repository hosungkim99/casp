#!/bin/bash
# L01 stage2: run the pipeline for all 79 binders SEQUENTIALLY in the current shell
# (no SLURM/sbatch). Skips binders that already finished. One log per binder.
# gnina falls back to CPU if no GPU is visible (works, just slower).
#
#   source this env first is handled below.
#   Run in background so it survives disconnect:
#     nohup bash run_l01_loop.sh > l01_loop.out 2>&1 &
#     tail -f l01_loop.out
set -u

CASP17=/path/to/casp17-ligand
source $CASP17/scripts/env_setup.sh

CONS=$CASP17/users/USERNAME/targets/L01/consensus
TEMPLATE=$CASP17/users/USERNAME/config/L01_template.conf
PIPE=$CASP17/users/USERNAME/pipeline/pipelines/run_pipeline.py
LOG=$CASP17/users/USERNAME/logs
mkdir -p "$LOG"

mapfile -t B < <(ls -d $CONS/L01*/ 2>/dev/null | xargs -n1 basename | sort)
echo "[loop] total binders: ${#B[@]}  start=$(date)"

done=0; skip=0; fail=0
for b in "${B[@]}"; do
  if [ -f "$CONS/$b/05_final/selection_summary.csv" ]; then
    echo "[skip] $b (already done)"; skip=$((skip+1)); continue
  fi
  echo "=== $b  start $(date) ==="
  sed "s/__LIG__/$b/g" "$TEMPLATE" > "$CONS/$b/pipeline.conf"
  if python3 "$PIPE" "$CONS/$b/pipeline.conf" --from cache > "$LOG/$b.log" 2>&1; then
    if [ -f "$CONS/$b/05_final/selection_summary.csv" ]; then
      echo "=== $b  OK $(date) ==="; done=$((done+1))
    else
      echo "=== $b  NO-OUTPUT (see $LOG/$b.log) ==="; fail=$((fail+1))
    fi
  else
    echo "=== $b  FAIL (see $LOG/$b.log) ==="; fail=$((fail+1))
  fi
done

echo "[loop] finished  done=$done skip=$skip fail=$fail  end=$(date)"
