#!/bin/bash
# run_stage1_frag.sh - fragment 1개에 대해 기존 파이프라인 스텝0~4 실행 (Stage1용).
#
# 스텝0만 boltz 어댑터(1_rank_poses_boltz.py)로 교체, 나머지 0b~4는 기존 스크립트 그대로.
# 출력: <OUT_ROOT>/<CID>/{00_collect,01_protein_clusters,02_pocket_candidates,
#                        03_pocket_validation,04_ligand_clusters}
# 결합확률(9정의)은 이후 2_binding 취합 스텝에서 04 대표 pose로 계산.
#
# 사용: (먼저 source env_setup.sh 로 p2rank/gnina/micromamba 준비)
#   bash run_stage1_frag.sh <CID> [RUNS_ROOT] [SMI_DIR] [OUT_ROOT]
# 예:  bash run_stage1_frag.sh L010001

set -euo pipefail
CID=${1:?usage: bash run_stage1_frag.sh <CID> [RUNS_ROOT] [SMI_DIR] [OUT_ROOT]}
RUNS_ROOT=${2:-/path/to/casp17-ligand/targets/L01/runs}
SMI_DIR=${3:-/path/to/casp17-data/targets_ligand/L01/ligands_smi}
OUT_ROOT=${4:-/path/to/casp17-ligand/users/USERNAME/targets/L01/outputs}

# 스크립트 위치 자동 유도: 이 파일(stage1/)의 상위 = pipeline 루트
ROOT=${SCRIPTS_DIR:-"$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"}
CORE="$ROOT/pipelines/core"   # 0b~4 코어 스텝
ST1="$ROOT/stage1"       # boltz 어댑터, 2_binding
# sci 스텝(gemmi+rdkit+numpy): micromamba가 PATH에 없어서 env python 직접 호출.
# vina_gpu env는 gemmi+rdkit 보유(cluster_bind에서 검증). 다른 env면 PY_SCI로 덮어쓰기.
SCI=${PY_SCI:-/path/to/casp17-ligand/.micromamba/envs/vina_gpu/bin/python}
STD=${PY_STD:-python3}
GPU=${GPU:-0}

FRAG="$RUNS_ROOT/$CID"
O="$OUT_ROOT/$CID"
mkdir -p "$O"

# per-fragment ligand.tsv (4단계 SC-RMSD 대칭순열용 SMILES)
SM=$(awk 'NF{print $1; exit}' "$SMI_DIR/$CID.smi")
printf "ID\tName\tSMILES\tTask\n0\t%s\t%s\tP\n" "$CID" "$SM" > "$O/ligand.tsv"

MT="$O/00_collect/master_table.csv"
CACHE="$O/00_collect/geom_cache.pkl"

echo "=== [0] boltz 어댑터 (수집) ==="
$STD "$ST1/1_rank_poses_boltz.py" --frag-dir "$FRAG" --out "$O/00_collect"

# 0 pose 방어: 팀 cofold가 없는 fragment(예: L010002)는 파이프라인 스킵 → 2_binding이 no_pose 처리
NPOSE=$(($(wc -l < "$MT") - 1))
if [ "$NPOSE" -lt 1 ]; then
    echo "=== [!] $CID: 0 pose — 파이프라인 스킵, no_pose 처리 ==="
    $SCI "$ST1/2_binding.py" --out-dir "$O" --frag-dir "$FRAG" \
         --ligand-tsv "$O/ligand.tsv" --gpu "$GPU" || true
    echo "=== 완료(no_pose): $O ==="
    exit 0
fi

echo "=== [0b] 기하 캐시 ==="
$SCI "$CORE/0b_cache_geometry.py" --table "$MT" --out "$CACHE"

echo "=== [1] 단백질 클러스터 ==="
$SCI "$CORE/1_protein_cluster.py" --table "$MT" --out "$O/01_protein_clusters" \
     --gap-frac 0.20 --cache "$CACHE"

echo "=== [2] 포켓 후보 ==="
$SCI "$CORE/2_pocket_candidates.py" --table "$MT" --out "$O/02_pocket_candidates" \
     --protein-clusters "$O/01_protein_clusters/protein_clusters.csv" \
     --reference "$O/01_protein_clusters/reference.txt" \
     --cache "$CACHE" --threshold 8.0 --topn 10

echo "=== [3] 포켓 검증 (p2rank + gnina) ==="
$SCI "$CORE/3_pocket_validate.py" \
     --candidates "$O/02_pocket_candidates/pocket_candidates.csv" \
     --out "$O/03_pocket_validation" \
     --gnina-cutoff -4.0 --p2rank-dist 6.0 --jobs 1 --gpus "$GPU"

echo "=== [4] 포켓별 포즈 클러스터 (SC-RMSD) ==="
$SCI "$CORE/4_ligand_cluster.py" \
     --members "$O/02_pocket_candidates/members.csv" \
     --validation "$O/03_pocket_validation/pocket_validation.csv" \
     --candidates "$O/02_pocket_candidates/pocket_candidates.csv" \
     --table "$MT" --ligand-tsv "$O/ligand.tsv" \
     --out "$O/04_ligand_clusters" --threshold 2.0

echo "=== [5] Stage1 결합확률 원자료 (2_binding) ==="
$SCI "$ST1/2_binding.py" --out-dir "$O" --frag-dir "$FRAG" \
     --ligand-tsv "$O/ligand.tsv" --gpu "$GPU"

echo ""
echo "=== 완료: $O ==="
ls "$O"
