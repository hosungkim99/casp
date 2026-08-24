#!/bin/bash
#SBATCH -J s1-all
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --gres=gpu:1
#SBATCH -o /path/to/casp17-ligand/users/USERNAME/logs/s1all_%x_%j.out
#SBATCH -e /path/to/casp17-ligand/users/USERNAME/logs/s1all_%x_%j.err
#
# 전 fragment에 대해 run_stage1_frag.sh(스텝0~5) 실행 (Stage1 배치).
#   GPU 1장/job. 2 GPU면 chunk 0,1 두 job 병렬(run_bt2_batch 방식).
#   idempotent: 05_stage1_binding/binding_row.csv 있으면 그 fragment skip(재실행 안전).
#
# 사용: sbatch run_stage1_all.sbatch <chunk_id> <n_chunks>
#   2 GPU 병렬:  sbatch run_stage1_all.sbatch 0 2 ; sbatch run_stage1_all.sbatch 1 2
#   테스트:      LIMIT=5 sbatch run_stage1_all.sbatch 0 1

set -euo pipefail
CHUNK=${1:-0}; NCHUNK=${2:-1}
LIMIT=${LIMIT:-0}

RUNS_ROOT=${RUNS_ROOT:-/path/to/casp17-ligand/targets/L01/runs}
SMI_DIR=${SMI_DIR:-/path/to/casp17-data/targets_ligand/L01/ligands_smi}
OUT_ROOT=${OUT_ROOT:-/path/to/casp17-ligand/users/USERNAME/targets/L01/outputs}
SC=${SCRIPTS_DIR:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"}   # 이 파일(stage1/)과 같은 폴더

LOGDIR="$OUT_ROOT/logs"
mkdir -p /path/to/casp17-ligand/users/USERNAME/logs "$OUT_ROOT" "$LOGDIR"
source /path/to/casp17-ligand/scripts/env_setup.sh

# GPU: slurm이 준 1장은 job 안에서 0
export GPU=${GPU:-0}

mapfile -t ALL < <(ls -d "$RUNS_ROOT"/L* 2>/dev/null | xargs -n1 basename | sort)
[ "${#ALL[@]}" -gt 0 ] || { echo "error: no L* under $RUNS_ROOT" >&2; exit 2; }
date; echo "▶ 전체 ${#ALL[@]} fragment, chunk $CHUNK/$NCHUNK, GPU=$GPU"

i=0; done=0; skip=0; fail=0
for cid in "${ALL[@]}"; do
    if [ $((i % NCHUNK)) -ne "$CHUNK" ]; then i=$((i+1)); continue; fi
    i=$((i+1))
    if [ -f "$OUT_ROOT/$cid/05_stage1_binding/binding_row.csv" ]; then
        skip=$((skip+1)); continue                       # 이미 됨 → skip
    fi
    if bash "$SC/run_stage1_frag.sh" "$cid" "$RUNS_ROOT" "$SMI_DIR" "$OUT_ROOT" \
         > "$LOGDIR/$cid.log" 2>&1; then
        done=$((done+1))
    else
        echo "  FAIL $cid (로그: $LOGDIR/$cid.log)"; fail=$((fail+1))
    fi
    [ $((done % 20)) -eq 0 ] && [ "$done" -gt 0 ] && echo "  진행 처리 $done, skip $skip, fail $fail (chunk $CHUNK)"
    [ "$LIMIT" -gt 0 ] && [ "$done" -ge "$LIMIT" ] && break
done
date; echo "▶ chunk $CHUNK 완료: 처리 $done, skip $skip, fail $fail"
echo "  취합:  python $SC/3_aggregate.py --outputs $OUT_ROOT --out-dir $OUT_ROOT/stage1"
