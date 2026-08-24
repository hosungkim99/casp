#!/bin/bash
#SBATCH -J oracle-boltz
#SBATCH -p normal
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH -t 08:00:00
#SBATCH -o /path/to/casp17-ligand/users/USERNAME/logs/oracle_boltz_%x_%j.out
#SBATCH -e /path/to/casp17-ligand/users/USERNAME/logs/oracle_boltz_%x_%j.err
#
# oracle 분자 boltz2 affinity cofold (job.sh 방식 재현).
#   분자마다 seed 1..NSEED 로 예측 → runs/<name>/seed_N/boltz_results_<name>/...
#   (fragment 1209개와 동일 레이아웃 → 우리 파이프라인이 그대로 읽음)
#   idempotent: seed 출력 있으면 skip.
# 사용: sbatch run_oracle_boltz.sbatch <chunk_id> <n_chunks>
#   여러 GPU 병렬:  sbatch run_oracle_boltz.sbatch 0 3 ; 1 3 ; 2 3

set -euo pipefail
CHUNK=${1:-0}; NCHUNK=${2:-1}
NSEED=${NSEED:-6}
SAMPLES=${SAMPLES:-5}

ORA=${ORA:-/path/to/casp17-ligand/users/USERNAME/oracle}
YAMLDIR="$ORA/affinity_yaml"
RUNS="$ORA/runs"
mkdir -p /path/to/casp17-ligand/users/USERNAME/logs "$RUNS"

source /path/to/casp17-ligand/scripts/env_setup.sh

mapfile -t YAMLS < <(ls "$YAMLDIR"/*.yaml 2>/dev/null | sort)
[ "${#YAMLS[@]}" -gt 0 ] || { echo "error: no yaml in $YAMLDIR" >&2; exit 2; }
date; echo "▶ oracle ${#YAMLS[@]} 분자 × seed 1..$NSEED, chunk $CHUNK/$NCHUNK"

i=0; done=0; skip=0
for y in "${YAMLS[@]}"; do
    if [ $((i % NCHUNK)) -ne "$CHUNK" ]; then i=$((i+1)); continue; fi
    i=$((i+1))
    name=$(basename "$y" .yaml)
    for s in $(seq 1 "$NSEED"); do
        out="$RUNS/$name/seed_$s"
        if compgen -G "$out/boltz_results_$name/predictions/$name/${name}_model_0.cif" > /dev/null; then
            skip=$((skip+1)); continue
        fi
        mkdir -p "$out"
        micromamba run -p "$CONDA_ENVS_PATH/boltz2" \
            boltz predict "$y" --out_dir "$out" \
            --accelerator gpu --no_kernels --diffusion_samples "$SAMPLES" --seed "$s" \
            > "$out/run.log" 2>&1 && done=$((done+1)) || echo "  FAIL $name seed $s"
    done
    echo "  $name 완료 (처리 $done, skip $skip)"
done
date; echo "▶ chunk $CHUNK 완료: 처리 $done, skip $skip"
echo "  다음: run_stage1_frag.sh 를 각 oracle 분자에 (RUNS_ROOT=$RUNS)"
