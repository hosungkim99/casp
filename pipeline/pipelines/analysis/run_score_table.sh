#!/bin/bash
#SBATCH --job-name=score_tbl
#SBATCH --output=/path/to/casp17-ligand/users/USERNAME/logs/scoretbl_%j.out
#SBATCH --error=/path/to/casp17-ligand/users/USERNAME/logs/scoretbl_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=12
#SBATCH --gres=gpu:6
# --- 위 자원 헤더는 클러스터 사정에 맞을 때만. 안 맞으면 gpu 수/시간 조정 ---
#
# 12_score_table.py 를 sbatch 로 백그라운드 실행.
# 사용법:
#   sbatch run_score_table.sh <outputs_dir> [cid] [extra args...]
# 예:
#   sbatch run_score_table.sh /path/to/casp17-ligand/.../targets/T2412/outputs T2412
#   sbatch run_score_table.sh /path/to/casp17-ligand/.../targets/T2412/outputs T2412 --minimize
#   MINIMIZE=1 sbatch run_score_table.sh /path/to/casp17-ligand/.../targets/T2412/outputs T2412
#
# 환경변수(선택): GPUS(기본 자동감지), JOBS(기본 GPU수), MINIMIZE(1이면 --minimize),
#                MAXPM(포켓당 멤버 상한), SMILES(SC-RMSD용).

set -euo pipefail

OUTPUTS=${1:?usage: sbatch run_score_table.sh <outputs_dir> [cid] [extra args...]}
CID=${2:-}
shift || true; [ $# -gt 0 ] && shift || true   # $@ = 나머지 extra args

source /path/to/casp17-ligand/scripts/env_setup.sh
SCI=/path/to/casp17-ligand/.micromamba/envs/vina_gpu/bin/python
# ⚠️ sbatch 는 이 스크립트를 /var/spool 로 복사해 실행 → BASH_SOURCE 는 스풀 경로라 쓸 수 없음.
#    12_score_table.py 위치를 고정(ANALYSIS_DIR 로 override 가능).
SC=${ANALYSIS_DIR:-/path/to/casp17-ligand/users/USERNAME/pipeline/pipelines/analysis}
mkdir -p /path/to/casp17-ligand/users/USERNAME/logs

# GPU 자동 감지 (nvidia-smi -L 개수 → 0,1,...,N-1). GPUS 로 수동 지정 가능.
if [ -z "${GPUS:-}" ]; then
    N=$(nvidia-smi -L 2>/dev/null | wc -l); N=${N:-1}; [ "$N" -lt 1 ] && N=1
    GPUS=$(seq -s, 0 $((N-1)))
fi
JOBS=${JOBS:-$(echo "$GPUS" | tr ',' '\n' | wc -l)}   # 기본 병렬수 = GPU 수

OPT=()
[ "${MINIMIZE:-0}" = "1" ] && OPT+=(--minimize)
[ -n "${MAXPM:-}" ] && OPT+=(--max-pocket-members "$MAXPM")
[ -n "${SMILES:-}" ] && OPT+=(--ligand-smiles "$SMILES")
[ -n "$CID" ] && OPT+=(--cid "$CID")

echo "[sbatch] node=$(hostname) start=$(date)"
echo "[sbatch] outputs=$OUTPUTS gpus=$GPUS jobs=$JOBS opt=${OPT[*]:-} extra=$*"
$SCI "$SC/12_score_table.py" --outputs "$OUTPUTS" --gpus "$GPUS" --jobs "$JOBS" "${OPT[@]}" "$@"
echo "[sbatch] done=$(date)"
