#!/bin/bash
#SBATCH -J protenix-s2
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH -t 24:00:00
#SBATCH -o /path/to/casp17-ligand/users/USERNAME/logs/protenix_s2_%A_%a.out
#SBATCH -e /path/to/casp17-ligand/users/USERNAME/logs/protenix_s2_%A_%a.err
# ── L01 stage2 pt2 재실행 ──────────────────────────────────────────────
# stage1과 달리 input.json을 인라인 생성하지 않고, 미리 만들어 검증한
#   stage2_pt2_inputs/<binder>/input.json  (유기 N copy + Zn) 을 그대로 사용.
# 폴더를 순회하므로 신규 binder(L010123)도 자동 포함(80개).
# 사용:
#   여러 GPU 분산:  sbatch --array=0-5 run_protenix_binders_stage2.sh    (필요시 %3 등으로 동시수 제한)
#   단일 직접실행:  bash run_protenix_binders_stage2.sh 0 1              (GPU 세션)
set -uo pipefail

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    CHUNK=$SLURM_ARRAY_TASK_ID; NCHUNK=${SLURM_ARRAY_TASK_COUNT:-1}
else
    CHUNK=${1:-0}; NCHUNK=${2:-1}
fi

BASE=/path/to/casp17-ligand/users/USERNAME/targets/L01
PROTENIX=/path/to/casp17-ligand/.micromamba/envs/protenix/bin/protenix
CONTAINER=$BASE/binders                 # per-binder 컨테이너
SEEDS=1,2,3,4,5,6                       # 6 seed 유지(모델 밸런스 30 pose)
SAMPLES=5

mkdir -p /path/to/casp17-ligand/users/USERNAME/logs
echo "[node $(hostname)] CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-none}"
mapfile -t CIDS < <(ls -1 "$CONTAINER" | grep -E '^L01[0-9]')   # CSV(L01.smiles..) 제외
echo "binder ${#CIDS[@]}개, chunk $CHUNK/$NCHUNK"

i=0; done=0; skip=0
for cid in "${CIDS[@]}"; do
    if [ $((i % NCHUNK)) -ne "$CHUNK" ]; then i=$((i+1)); continue; fi
    i=$((i+1))
    in="$CONTAINER/$cid/inputs/pt2.json"
    out="$CONTAINER/$cid/results/pt2"
    [ -f "$in" ] || { echo "[$cid] input 없음, skip"; continue; }
    # 완료 판정: 마지막 seed·마지막 sample cif 존재(-e 5 → sample_4가 마지막)
    if compgen -G "$out/$cid/seed_6/predictions/${cid}_sample_4.cif" > /dev/null; then
        echo "[$cid] done, skip"; skip=$((skip+1)); continue
    fi
    mkdir -p "$out"
    echo "[$cid] run"
    "$PROTENIX" pred -i "$in" -o "$out" \
        -s "$SEEDS" -e "$SAMPLES" --use_msa True \
        > "$out/run.log" 2>&1 && done=$((done+1)) || echo "  FAIL $cid (see $out/run.log)"
done
echo "chunk $CHUNK 완료: 실행 $done, skip $skip"
