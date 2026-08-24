#!/bin/bash
#SBATCH -J af3-s2
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH -t 24:00:00
#SBATCH -o /path/to/casp17-ligand/users/USERNAME/logs/af3_s2_%A_%a.out
#SBATCH -e /path/to/casp17-ligand/users/USERNAME/logs/af3_s2_%A_%a.err
# ── L01 stage2 af3(af3_red) 재실행 ─────────────────────────────────────
# 입력 binders/<id>/inputs/af3.json (유기 N copy + Zn, modelSeeds=1..6). 출력 binders/<id>/results/af3/.
# run_alphafold.py가 modelSeeds로 seed_* 생성(1회 실행 = 6 seed).
# 사용:  sbatch --array=0-5%3 run_af3_binders_stage2.sh   / bash run_af3_binders_stage2.sh 0 1
set -uo pipefail

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    CHUNK=$SLURM_ARRAY_TASK_ID; NCHUNK=${SLURM_ARRAY_TASK_COUNT:-1}
else
    CHUNK=${1:-0}; NCHUNK=${2:-1}
fi

BASE=/path/to/casp17-ligand/users/USERNAME/targets/L01
CONTAINER=$BASE/binders
AF3PY=/path/to/casp17-ligand/.micromamba/envs/af3_red/bin/python
AF3=/path/to/casp17-ligand/models/af3_red/run_alphafold.py
DB=/path/to/casp/af3_db/af3
export LD_LIBRARY_PATH=${CUDA_LIB:+$CUDA_LIB:}${LD_LIBRARY_PATH:-}   # CUDA_LIB=클러스터 CUDA lib 경로 (env_setup.sh에서 설정)

mkdir -p /path/to/casp17-ligand/users/USERNAME/logs
mapfile -t CIDS < <(ls -1 "$CONTAINER" | grep -E '^L01[0-9]')
echo "[node $(hostname)] binder ${#CIDS[@]}, chunk $CHUNK/$NCHUNK"

i=0
for cid in "${CIDS[@]}"; do
    if [ $((i % NCHUNK)) -ne "$CHUNK" ]; then i=$((i+1)); continue; fi
    i=$((i+1))
    in="$CONTAINER/$cid/inputs/af3.json"
    out="$CONTAINER/$cid/results/af3"
    [ -f "$in" ] || { echo "[$cid] af3.json 없음, skip"; continue; }
    shopt -s nullglob
    # 1) 이미 정규화됨(canonical) → skip
    if compgen -G "$out/seed_6/sample_4/model.cif" > /dev/null; then
        echo "[$cid] done, skip"; continue
    fi
    # 2) af3_red raw 출력 없으면 실행 (있으면 실행 건너뛰고 정규화만)
    if ! compgen -G "$out/*/seed-*_sample-*/*_model.cif" > /dev/null; then
        mkdir -p "$out"
        "$AF3PY" "$AF3" --json_path "$in" --model_dir "$DB" --db_dir "$DB" \
            --output_dir "$out" --bias_weight 100.0 --bias_sigma 2.0 \
            > "$out/run.log" 2>&1 || { echo "  FAIL $cid (see $out/run.log)"; continue; }
    fi
    # 3) 정규화: af3_red <name>/seed-S_sample-M/*_{model.cif,summary_confidences.json}
    #    → canonical seed_S/sample_M/{model.cif, summary.json} (collect_af3 호환)
    for d in "$out"/*/seed-*_sample-*/; do
        s=$(basename "$d" | sed -E 's/.*seed-([0-9]+).*/\1/')
        m=$(basename "$d" | sed -E 's/.*sample-([0-9]+).*/\1/')
        nd="$out/seed_${s}/sample_${m}"; mkdir -p "$nd"
        cp -f "$d"*_model.cif                "$nd/model.cif"   2>/dev/null
        cp -f "$d"*_summary_confidences.json "$nd/summary.json" 2>/dev/null
    done
    # 원본 af3_red 폴더(<name>/, seed-*_sample-* 보유) 제거 → canonical seed_*/ 만 남김
    for nmdir in "$out"/*/; do
        compgen -G "${nmdir}seed-*_sample-*" > /dev/null && rm -rf "$nmdir"
    done
    echo "[$cid] af3 done"
done
echo "chunk $CHUNK 완료"
