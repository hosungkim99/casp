# stage2 — Usage

## 사전 준비물 (필수)
1. **Stage1 완료** — fragment마다 `<outputs>/<cid>/04_ligand_clusters/ligand_clusters.csv` 존재.
2. **CASP 그룹코드** — `AUTHOR=<그룹번호>` (LG 파일 헤더·파일명).
3. **환경** — `source env_setup.sh` (gnina). sci python(`PY_SCI`).
4. **GPU** — 5_final_select의 gnina 재채점용.
5. (실제 제출은) CASP가 공개한 **결합 확정 fragment 목록** — 그 fragment만 제출.

## 전체 실행 (배치)
```bash
source /path/to/casp17-ligand/.../scripts/env_setup.sh
# 3청크 예시 (그룹코드 넣기)
for i in 0 1 2; do AUTHOR=<그룹번호> sbatch run_stage2_all.sh $i 3; done
```
- Stage1이 다 끝난 뒤 실행(04 있는 fragment만 처리). 중간에 돌렸으면 나중에 재실행(완료분 skip).

## fragment 1개만
```bash
AUTHOR=<그룹번호> bash run_stage2_frag.sh L010001 <OUT_ROOT>
```

## 배치 후 취합
```bash
OUT=/path/to/casp17-ligand/.../users/USERNAME/targets/L01/outputs
python 1_aggregate.py --outputs $OUT --out-dir $OUT/final/stage2
# -> final/stage2/poses.csv + submit/<cid>/<cid>LG.txt
```

## 제출 tarball (CASP 형식)
CASP Stage2는 **개별 포즈 파일**(`L0xxxxxLG<group>_N`)을 `./L01/` 에 모아 tar:
```bash
# (submit 형식을 개별파일로 맞춘 뒤)
tar -czf L01LG<group>.tgz -C <submit_root> L01
```
→ 이 `.tgz` 를 Ligand Series Prediction Upload Form에 업로드.

## Task 옵션
- `TASK=P` (기본) = pose 중심 / `TASK=PA` = affinity 비중↑ (5_final_select 가중치 분기).

## 주의
- LG 파일 = **수용체 PDB + 리간드 MDL 합본**(pdb/sdf 따로 아님).
- Stage2 마감(8/27)이라 여유 있음. binder 목록 공개 후 그 fragment만 제출.
