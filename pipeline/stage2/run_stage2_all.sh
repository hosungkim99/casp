#!/bin/bash
#SBATCH -J s2-all
#SBATCH -p normal
#SBATCH -c 8
#SBATCH -t 24:00:00
#SBATCH --gres=gpu:1
#SBATCH -o /path/to/casp17-ligand/users/USERNAME/logs/s2all_%x_%j.out
#SBATCH -e /path/to/casp17-ligand/users/USERNAME/logs/s2all_%x_%j.err
#
# 이미 00~04가 만들어진 fragment들에 Stage2(5 final_select + 9 casp_lg) 배치.
#   idempotent: 05_final/selection_summary.csv 있으면 skip.
# 사용: AUTHOR=<그룹코드> sbatch run_stage2_all.sh <chunk_id> <n_chunks>
#   2 GPU:  AUTHOR=xxxx sbatch run_stage2_all.sh 0 2 ; ... 1 2
#   테스트: LIMIT=5 sbatch run_stage2_all.sh 0 1

set -euo pipefail
CHUNK=${1:-0}; NCHUNK=${2:-1}
LIMIT=${LIMIT:-0}

OUT_ROOT=${OUT_ROOT:-/path/to/casp17-ligand/users/USERNAME/targets/L01/outputs}
SC=${SCRIPTS_DIR:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"}   # 이 파일(stage2/)과 같은 폴더
export AUTHOR=${AUTHOR:-CHANGE-ME}
export TASK=${TASK:-P}
export GPU=${GPU:-0}

LOGDIR="$OUT_ROOT/logs_stage2"
mkdir -p /path/to/casp17-ligand/users/USERNAME/logs "$LOGDIR"
source /path/to/casp17-ligand/scripts/env_setup.sh

# 대상 = 04까지 완료된 fragment
mapfile -t ALL < <(ls -d "$OUT_ROOT"/L*/04_ligand_clusters 2>/dev/null \
                   | sed 's#/04_ligand_clusters##' | xargs -n1 basename | sort)
[ "${#ALL[@]}" -gt 0 ] || { echo "error: 04 완료 fragment 없음 under $OUT_ROOT" >&2; exit 2; }
date; echo "▶ Stage2 대상 ${#ALL[@]} fragment, chunk $CHUNK/$NCHUNK, AUTHOR=$AUTHOR, TASK=$TASK"

i=0; done=0; skip=0; fail=0
for cid in "${ALL[@]}"; do
    if [ $((i % NCHUNK)) -ne "$CHUNK" ]; then i=$((i+1)); continue; fi
    i=$((i+1))
    if [ -f "$OUT_ROOT/$cid/05_final/selection_summary.csv" ]; then
        skip=$((skip+1)); continue
    fi
    if bash "$SC/run_stage2_frag.sh" "$cid" "$OUT_ROOT" > "$LOGDIR/$cid.log" 2>&1; then
        done=$((done+1))
    else
        echo "  FAIL $cid (로그: $LOGDIR/$cid.log)"; fail=$((fail+1))
    fi
    [ $((done % 20)) -eq 0 ] && [ "$done" -gt 0 ] && echo "  진행 처리 $done, skip $skip, fail $fail (chunk $CHUNK)"
    [ "$LIMIT" -gt 0 ] && [ "$done" -ge "$LIMIT" ] && break
done
date; echo "▶ chunk $CHUNK 완료: 처리 $done, skip $skip, fail $fail"
echo "  취합:  python $SC/1_aggregate.py --outputs $OUT_ROOT --out-dir $OUT_ROOT/final/stage2"
