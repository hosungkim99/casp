#!/bin/bash
#SBATCH -J boltz-s2
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH -t 24:00:00
#SBATCH -o /path/to/casp17-ligand/users/USERNAME/logs/boltz_s2_%A_%a.out
#SBATCH -e /path/to/casp17-ligand/users/USERNAME/logs/boltz_s2_%A_%a.err
# ── L01 stage2 bt2(Boltz2) 재실행 ──────────────────────────────────────
# 입력 binders/<id>/inputs/<id>.yaml (유기 N copy + Zn). 출력 binders/<id>/results/bt2/seed_*/.
# 6 seed × 5 sample = 30 pose (pt2/af3와 밸런스).
# 사용:  sbatch --array=0-5%3 run_bt2_binders_stage2.sh    / bash run_bt2_binders_stage2.sh 0 1
set -uo pipefail

if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
    CHUNK=$SLURM_ARRAY_TASK_ID; NCHUNK=${SLURM_ARRAY_TASK_COUNT:-1}
else
    CHUNK=${1:-0}; NCHUNK=${2:-1}
fi

BASE=/path/to/casp17-ligand/users/USERNAME/targets/L01
CONTAINER=$BASE/binders
BOLTZ=/path/to/casp17-ligand/.micromamba/envs/boltz2/bin/boltz
SEEDS="1 2 3 4 5 6"
SAMPLES=5

mkdir -p /path/to/casp17-ligand/users/USERNAME/logs
mapfile -t CIDS < <(ls -1 "$CONTAINER" | grep -E '^L01[0-9]')
echo "[node $(hostname)] binder ${#CIDS[@]}, chunk $CHUNK/$NCHUNK"

i=0
for cid in "${CIDS[@]}"; do
    if [ $((i % NCHUNK)) -ne "$CHUNK" ]; then i=$((i+1)); continue; fi
    i=$((i+1))
    yaml="$CONTAINER/$cid/inputs/$cid.yaml"
    [ -f "$yaml" ] || { echo "[$cid] yaml 없음, skip"; continue; }
    for s in $SEEDS; do
        out="$CONTAINER/$cid/results/bt2/seed_$s"
        # 완료 판정: 마지막 sample(model_4) cif 존재
        if compgen -G "$out/boltz_results_$cid/predictions/$cid/${cid}_model_4.cif" > /dev/null; then
            continue
        fi
        mkdir -p "$out"
        "$BOLTZ" predict "$yaml" --out_dir "$out" --accelerator gpu --no_kernels \
            --diffusion_samples $SAMPLES --seed "$s" > "$out/run.log" 2>&1 \
            || echo "  FAIL $cid seed$s (see $out/run.log)"
        # 중간파일 삭제(predictions만 남김) → 파일 수·메타데이터 부하 대폭 감소
        rm -rf "$out/boltz_results_$cid"/{lightning_logs,msa,processed} 2>/dev/null
    done
    echo "[$cid] bt2 done"
done
echo "chunk $CHUNK 완료"
