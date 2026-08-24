#!/bin/bash
# L01 stage2: 79 binder 파이프라인 일괄 실행 (SLURM array).
# 공통 conf 1장(L01_template.conf)의 __LIG__ 를 binder ID로 치환 -> run_pipeline.py --from cache.
# 제출:  sbatch --array=0-78 run_l01_pipeline.sh
# 특정 스텝부터 재실행:  sbatch --array=0-78 run_l01_pipeline.sh --from pocket_validate --force
#SBATCH -J l01-pipe
#SBATCH -p normal
#SBATCH --gres=gpu:1          # gnina GPU. <node-name>(RTX6000 Blackwell) 회피 = RTX4090 명시
#SBATCH --cpus-per-task=8
#SBATCH -t 12:00:00
#SBATCH -o /path/to/casp17-ligand/users/USERNAME/logs/l01pipe_%A_%a.out
#SBATCH -e /path/to/casp17-ligand/users/USERNAME/logs/l01pipe_%A_%a.err
set -e

CASP17=/path/to/casp17-ligand
source $CASP17/scripts/env_setup.sh   # gnina/p2rank/java 등 PATH

# ── 경로(서버 실제 배치) ─────────────────────────────────────────────
CONS=$CASP17/users/USERNAME/targets/L01/consensus
TEMPLATE=$CASP17/users/USERNAME/config/L01_template.conf
# core 파이프라인 지휘자. ★실제 위치 확인해서 맞추기★ (L01/scripts와 별개인 공유 pipeline 폴더)
PIPE=$CASP17/users/USERNAME/pipeline/pipelines/run_pipeline.py
# ─────────────────────────────────────────────────────────────────────

# binder 목록 = consensus 에 합본테이블이 있는 폴더들(=79개)
mapfile -t B < <(ls -d $CONS/L01*/ 2>/dev/null | xargs -n1 basename | sort)
IDX=${SLURM_ARRAY_TASK_ID:-0}
LIG=${B[$IDX]}
[ -z "$LIG" ] && { echo "no binder at index $IDX (총 ${#B[@]}개)"; exit 1; }

CONF=$CONS/$LIG/pipeline.conf
sed "s/__LIG__/$LIG/g" "$TEMPLATE" > "$CONF"

echo "[l01-pipe] idx=$IDX lig=$LIG node=$(hostname) start=$(date)"
python3 "$PIPE" "$CONF" --from cache "$@"
echo "[l01-pipe] lig=$LIG done=$(date)"
