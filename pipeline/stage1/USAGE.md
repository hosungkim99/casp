# stage1 — Usage

## 사전 준비물 (필수)
1. **팀 boltz cofold 결과** — `<runs>/<cid>/seed_*/boltz_results_*/predictions/<cid>/` 에 pose cif + confidence json + `affinity_<cid>.json`
2. **리간드 SMILES** — `<smi-dir>/<cid>.smi` (fragment마다 1줄)
3. **환경** — `source /path/to/casp17-ligand/.../scripts/env_setup.sh` (p2rank, gnina singularity)
   - sci python: gemmi+rdkit+numpy 있는 env (예: `.micromamba/envs/vina_gpu/bin/python`). 러너 `PY_SCI`로 지정.
4. **GPU** — 3단계 p2rank + gnina용.
5. common/pipelines가 같은 `pipeline/` 아래 있을 것 (import 부트스트랩이 common/ 자동 탐색).

## 전체 실행 (배치)
```bash
source /path/to/casp17-ligand/.../scripts/env_setup.sh
# GPU 3장 병렬 예시
for i in 0 1 2; do sbatch run_stage1_all.sh $i 3; done
```
경로 기본값이 다르면 env로: `RUNS_ROOT=... SMI_DIR=... OUT_ROOT=... sbatch run_stage1_all.sh 0 1`

## fragment 1개만
```bash
bash run_stage1_frag.sh L010001 <RUNS_ROOT> <SMI_DIR> <OUT_ROOT>
```

## 배치 후 취합
```bash
OUT=/path/to/casp17-ligand/.../users/USERNAME/targets/L01/outputs
python 3_aggregate.py --outputs $OUT --out-dir $OUT/stage1
```

## 최종 제출 파일
```bash
python 4_finalize.py --stage1-dir $OUT/stage1 \
  --definition prob_combined --target L01 --group <그룹번호> --out-dir $OUT/final/stage1
# -> final/stage1/L01LG<그룹>.bind.txt   (이걸 업로드)
```

## 분석 (선택)
```bash
python 5_defs_viz.py --csv $OUT/stage1/binding_scores.csv --out $OUT/stage1/compare.png
python 6_consensus_pocket.py --csv $OUT/stage1/pocket_clusters.csv --out $OUT/stage1/consensus_pocket.txt --min-frac 0.5
```

## 개별 스크립트 실행 예시
각 `.py` 하단 주석에 실행 예시가 있음. 핵심:
```bash
python 1_rank_poses_boltz.py --frag-dir <runs>/<cid> --out <out>/00_collect
python 2_binding.py --out-dir <out>/<cid> --frag-dir <runs>/<cid> --gpu 0
```

## 주의
- **3_aggregate는 전 fragment 완료 후** — 9정의 rank가 라이브러리 전체를 봐야 의미 있음.
- run_stage1_all은 시작 시점의 fragment 목록을 고정 → stage1 다 끝난 뒤 한 번 더 돌리면 누락분 채움(완료분 skip).
- 0 pose fragment(예 L010002)는 자동 no_pose 처리.
