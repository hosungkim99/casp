#!/bin/bash
# run_stage2_frag.sh - fragment 1개에 대해 Stage2 단계(5 final_select + 9 casp_lg) 실행.
#
# 전제: run_stage1_frag.sh 로 00~04가 이미 있음. 여기서 그 04를 받아 ≤5개 포즈 선정 + 제출포맷.
# 출력: <OUT_ROOT>/<CID>/{05_final, 08_casp_lg}
# 사용: (먼저 source env_setup.sh) bash run_stage2_frag.sh <CID> [OUT_ROOT]

set -euo pipefail
CID=${1:?usage: bash run_stage2_frag.sh <CID> [OUT_ROOT]}
OUT_ROOT=${2:-/path/to/casp17-ligand/users/USERNAME/targets/L01/outputs}
ROOT=${SCRIPTS_DIR:-"$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"}
CORE="$ROOT/pipelines/core"   # 5_final_select, 9_make_casp_lg
SCI=${PY_SCI:-/path/to/casp17-ligand/.micromamba/envs/vina_gpu/bin/python}
GPU=${GPU:-0}
AUTHOR=${AUTHOR:-CHANGE-ME}     # CASP 그룹 코드
TASK=${TASK:-P}                 # P=pose 중심 / PA=affinity 비중↑

O="$OUT_ROOT/$CID"
LC="$O/04_ligand_clusters/ligand_clusters.csv"
[ -f "$LC" ] || { echo "[stage2] $CID: 04 없음 — skip"; exit 0; }

echo "=== [5] 최종 ≤5 선정 ==="
$SCI "$CORE/5_final_select.py" --ligand-clusters "$LC" \
     --validation "$O/03_pocket_validation/pocket_validation.csv" \
     --out "$O/05_final" --top 5 --task "$TASK" --gpus "$GPU"

echo "=== [9] CASP LG 제출포맷 ==="
$SCI "$CORE/9_make_casp_lg.py" --summary "$O/05_final/selection_summary.csv" \
     --out-dir "$O/08_casp_lg" --ligand-tsv "$O/ligand.tsv" \
     --target "$CID" --author "$AUTHOR"

echo "=== 완료(stage2): $O ==="
ls "$O/05_final" "$O/08_casp_lg"
