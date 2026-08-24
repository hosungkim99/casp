#!/bin/bash
#SBATCH --job-name=lig_pipe
#SBATCH --output=/path/to/casp17-ligand/users/USERNAME/logs/pipe_%j.out
#SBATCH --error=/path/to/casp17-ligand/users/USERNAME/logs/pipe_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:1         # gnina GPU용. 없으면 자동 CPU 폴백 → 일단 끔
# (선택) 특정 노드 회피가 필요하면 주석 해제: --exclude=<node-name>
# --- 위 자원 헤더는 클러스터 사정에 맞을 때만. 안 맞으면 제출 거부됨 ---
#
# 사용법:
#   sbatch run_pipeline.sh <config.conf> [run_pipeline.py 추가인자...]
# 예:
#   sbatch run_pipeline.sh /path/to/casp17-ligand/.../config/T2410.conf
#   sbatch run_pipeline.sh /path/to/casp17-ligand/.../config/T2410.conf --from collect --force

set -e

# 환경 로드 (micromamba/CASP17 등)
source /path/to/casp17-ligand/scripts/env_setup.sh

# 이 스크립트(pipelines/)와 같은 폴더의 run_pipeline.py를 사용 → 재구성 후에도 경로 안전
SC=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

CONF=${1:?usage: sbatch run_pipeline.sh <config.conf> [--from STEP --force ...]}
shift   # 나머지 인자($@)는 run_pipeline.py로 전달

echo "[sbatch] node=$(hostname) conf=$CONF start=$(date)"
python3 "$SC/run_pipeline.py" "$CONF" "$@"
echo "[sbatch] done=$(date)"

# 재실행(이어하기): run_pipeline 은 idempotent → 산출물 있으면 skip.
#   중간 실패 후 재제출하면 끝난 스텝은 건너뛰고 이어서 진행.
#   특정 스텝부터 강제 재실행: sbatch run_pipeline.sh <conf> --from pocket_validate --force
