# stage2 — Scripts

(실행법 [USAGE.md](USAGE.md), 컬럼 정의 [FEATURES.md](FEATURES.md))

## 1_aggregate.py  (전체 취합, stdlib)
- **역할**: 전 fragment 최종 포즈(05_final)와 LG 파일(08_casp_lg)을 모아 제출용으로 정리.
- **input**: `<outputs>/*/05_final/selection_summary.csv`, `<outputs>/*/08_casp_lg/<cid>LG_all_models.txt`
- **output**:
  - `<out-dir>/poses.csv` — 전 리간드 ≤5 포즈 한 표(점수·구조 경로)
  - `<out-dir>/submit/<cid>/<cid>LG.txt` — 각 리간드 제출 포즈 파일
- **비고**: CASP 제출은 개별 포즈 파일(`L0xxxxxLG<group>_N`) + tarball 형식이 정식. 현재는 합본 복사(형식 조정 필요).

## 러너 (오케스트레이터, bash)
- **run_stage2_frag.sh**: fragment 1개 Stage2 = core `5_final_select` + `9_make_casp_lg`.
  - input: `<CID> [OUT_ROOT]` (+ env `AUTHOR`=그룹코드, `TASK`=P/PA)
  - output: `<OUT_ROOT>/<CID>/{05_final, 08_casp_lg}`
  - 전제: `<CID>/04_ligand_clusters/ligand_clusters.csv` 존재(Stage1 04)
- **run_stage2_all.sh**: 04 완료된 전 fragment 배치(idempotent skip on 05_final).
  - input: `<chunk> <n_chunks>` (+ env `AUTHOR`)
  - output: 각 fragment `05_final`/`08_casp_lg` + `logs_stage2/`

## 재사용하는 코어 (pipelines/core)
- **5_final_select.py**: 04 클러스터 후보 → 종합점수 → ≤5개 → `05_final/{model_*.cif, selection_summary.csv, SELECTION_RATIONALE.md}`
- **9_make_casp_lg.py**: selection_summary → CASP LG 포맷 `08_casp_lg/<cid>LG_model*.txt`
