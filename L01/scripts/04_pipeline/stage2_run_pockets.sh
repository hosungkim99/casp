#!/bin/bash
# ── L01 stage2: pt2 결과에 스텝0b/1/2(pocket 발견)까지 돌리기 (multi-copy 선정 준비) ──
# 스텝3·4(도킹/클러스터)는 multi-copy에서 퇴화하므로 안 돎. 선정은 stage2_select_multicopy.py로.
# 사용:  bash stage2_run_pockets.sh                 # stage2_pt2_out 전체
#        bash stage2_run_pockets.sh L010016 L010589 L010039   # 지정 몇 개만(검증용)
set -uo pipefail

BASE=$CASP17/users/USERNAME/targets/L01
CONTAINER=$BASE/binders                 # per-binder 컨테이너 <id>/{consensus,consensus_s2,inputs,results}
MODELS=${MODELS:-pt2}                   # 수집할 모델(예: MODELS=pt2,bt2,af3 로 3모델)
PIPE=$CASP17/users/USERNAME/pipeline/pipelines/run_pipeline.py
PYSCI=$CONDA_ENVS_PATH/boltz2/bin/python
SCR=$BASE/scripts

if [ $# -gt 0 ]; then BINDERS="$*"; else BINDERS=$(ls "$CONTAINER" | grep -E '^L01[0-9]'); fi   # CSV(L01.smiles..) 제외

for b in $BINDERS; do
    D=$CONTAINER/$b/consensus_s2         # 파이프라인 out_dir(수집·포켓 결과)
    mkdir -p "$D/00_collect"
    # 1) 수집: 모델은 <id>/results/<model>/ 에서 읽고, 테이블은 <id>/consensus_s2/00_collect 에 씀
    python3 "$SCR/03_collect/collect_consensus.py" --only "$b" --models "$MODELS" \
        --results-base "$CONTAINER" --results-subdir results \
        --out-root "$CONTAINER" --out-subdir consensus_s2 || { echo "  collect 실패 $b"; continue; }
    # 2) config (스텝0b/1/2는 ligand_tsv 안 씀 → NA. results_dir는 collect skip이라 미사용)
    cat > "$D/$b.conf" <<EOF
target=$b
results_dir=$CONTAINER
out_dir=$D
ligand_tsv=NA
protein_chain=A
task=P
python_sci=$PYSCI
EOF
    # 3) 스텝 0b(cache) → 1(protein_cluster) → 2(pocket_candidates)
    #    --force: 재수집(경로이동/모델추가)으로 master_table 바뀌었으니 캐시·포켓·멤버 새 경로로 재생성
    python3 "$PIPE" "$D/$b.conf" --only cache            --force || { echo "  cache 실패 $b"; continue; }
    python3 "$PIPE" "$D/$b.conf" --only protein_cluster  --force || { echo "  protein_cluster 실패 $b"; continue; }
    python3 "$PIPE" "$D/$b.conf" --only pocket_candidates --force || { echo "  pocket 실패 $b"; continue; }
    # 4) p2rank 교차검증 (per-pocket, multi-copy 대응 독립 스크립트) — 포켓이 druggable cavity 근처인지
    "$PYSCI" "$SCR/04_pipeline/validate_pockets_p2rank.py" \
        --candidates "$D/02_pocket_candidates/pocket_candidates.csv" \
        --reference  "$D/01_protein_clusters/reference.txt" \
        --out        "$D/03b_p2rank_check" > /dev/null 2>&1 || echo "  p2rank 검증 실패 $b(무시가능)"
    echo "=== $b : pocket 발견 + p2rank 검증 완료 ==="
done
echo "전부 완료. 선정:  python3 $SCR/05_select/stage2_select_multicopy.py --cons $CONTAINER"
