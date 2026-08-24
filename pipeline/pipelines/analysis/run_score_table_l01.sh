#!/bin/bash
#SBATCH --job-name=stbl_l01
#SBATCH --output=/path/to/casp17-ligand/users/USERNAME/logs/stbl_l01_%j.out
#SBATCH --error=/path/to/casp17-ligand/users/USERNAME/logs/stbl_l01_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:6
#
# L01 전 fragment(1209개)의 score_table을 **단일 CSV에 누적**(--append).
# 사용:  sbatch run_score_table_l01.sh
#   재개: 이미 CSV에 있는 cid(fragment)는 건너뜀 → 중간에 죽어도 다시 제출하면 이어서.
#   처음부터: OUT 파일을 지우고 제출.
# 환경변수(선택): MINIMIZE(기본 1) · OUT · OUT_ROOT · GPUS · JOBS.
#   BINDERS=<파일>  : 주면 그 파일의 binding=TRUE fragment만 처리(binder만).
#                     파일 형식: 'CASP_ID,...,binding' CSV (마지막 컬럼이 TRUE인 행) 또는 cid 목록.

set -uo pipefail
source /path/to/casp17-ligand/scripts/env_setup.sh
SCI=/path/to/casp17-ligand/.micromamba/envs/vina_gpu/bin/python
# sbatch 스풀 복사 대비: 스크립트 위치를 고정(override 가능)
ANALYSIS_DIR=${ANALYSIS_DIR:-/path/to/casp17-ligand/users/USERNAME/pipeline/pipelines/analysis}
OUT_ROOT=${OUT_ROOT:-/path/to/casp17-ligand/users/USERNAME/targets/L01/outputs}
OUT=${OUT:-/path/to/casp17-ligand/users/USERNAME/targets/L01/score_table.csv}
mkdir -p /path/to/casp17-ligand/users/USERNAME/logs "$(dirname "$OUT")"

# GPU 자동감지 (nvidia-smi 개수). GPUS/JOBS 로 수동 지정 가능.
if [ -z "${GPUS:-}" ]; then
    N=$(nvidia-smi -L 2>/dev/null | wc -l); N=${N:-1}; [ "$N" -lt 1 ] && N=1
    GPUS=$(seq -s, 0 $((N-1)))
fi
JOBS=${JOBS:-$(echo "$GPUS" | tr ',' '\n' | wc -l)}
OPT=(); [ "${MINIMIZE:-1}" = "1" ] && OPT+=(--minimize)

if [ -n "${BINDERS:-}" ]; then
    # binder 목록 파일: 마지막 컬럼이 TRUE인 행의 첫 컬럼(cid)만 (헤더/공백 제외).
    mapfile -t FRAGS < <(awk -F, '{sub(/\r$/,"")} toupper($NF)=="TRUE"{print $1}' "$BINDERS" \
                         | grep -E '^L0' | sort -u)
    echo "[l01] BINDERS=$BINDERS → binder ${#FRAGS[@]}개만 처리"
else
    mapfile -t FRAGS < <(ls -d "$OUT_ROOT"/L0*/ 2>/dev/null | xargs -n1 basename | sort)
fi
[ "${#FRAGS[@]}" -gt 0 ] || { echo "error: 처리할 fragment 없음 (OUT_ROOT=$OUT_ROOT, BINDERS=${BINDERS:-none})" >&2; exit 2; }
date; echo "[l01] fragment ${#FRAGS[@]}, gpus=$GPUS jobs=$JOBS opt=${OPT[*]:-} → $OUT"

i=0; done=0; skip=0; fail=0
for cid in "${FRAGS[@]}"; do
    i=$((i+1))
    # 재개: 이미 누적된 fragment는 건너뜀 (cid 가 각 행 첫 컬럼)
    if [ -f "$OUT" ] && grep -q "^$cid," "$OUT"; then skip=$((skip+1)); continue; fi
    # 클러스터 없는 fragment 는 건너뜀
    [ -s "$OUT_ROOT/$cid/04_ligand_clusters/cluster_members.csv" ] || { skip=$((skip+1)); continue; }
    if $SCI "$ANALYSIS_DIR/12_score_table.py" \
            --outputs "$OUT_ROOT/$cid" --cid "$cid" --out "$OUT" --append \
            --jobs "$JOBS" --gpus "$GPUS" "${OPT[@]}" > "$OUT_ROOT/$cid/12_score_table.log" 2>&1; then
        done=$((done+1))
    else
        echo "  FAIL $cid (로그: $OUT_ROOT/$cid/12_score_table.log)"; fail=$((fail+1))
    fi
    [ $((i % 50)) -eq 0 ] && echo "  진행 $i/${#FRAGS[@]}: done $done skip $skip fail $fail"
done
date; echo "[l01] 완료: done $done, skip $skip, fail $fail → $OUT"
