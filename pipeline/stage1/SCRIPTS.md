# stage1 — Scripts

각 스크립트의 역할 / 입력 / 출력. (실행법은 [USAGE.md](USAGE.md), 컬럼 정의는 [FEATURES.md](FEATURES.md))

## 1_rank_poses_boltz.py  (수집 어댑터, stdlib)
- **역할**: boltz2 fragment 출력을 파이프라인이 읽는 `master_table.csv` 로 변환 (af3식 `0_rank_poses`의 boltz판).
- **input**: `<runs>/<cid>/seed_*/boltz_results_*/predictions/<cid>/confidence_<cid>_model_*.json` (+ 같은 폴더 `<cid>_model_*.cif`)
- **output**: `<out>/00_collect/master_table.csv` (30행: rank, seed, sample, ligand_iptm, plddt, cif 경로 …), `progress.txt`

## 2_binding.py  (대표 pose 스코어, sci + gnina)
- **역할**: core 클러스터링 결과에서 지배 포켓·지배 포즈의 **대표 pose 1개**를 뽑아 결합확률 원자료 계산.
- **input**: `<out>/{02_pocket_candidates, 03_pocket_validation, 04_ligand_clusters}/*.csv`, affinity json(`<runs>/<cid>/.../affinity_<cid>.json`), `ligand.tsv`(SMILES→hac)
- **output**: `<out>/05_stage1_binding/binding_row.csv` (1행: boltz/cnn/cnnaff/aff/hac + 포켓/클러스터 요약)

## 3_aggregate.py  (전체 취합, stdlib)
- **역할**: 전 fragment `binding_row.csv` 를 모아 **라이브러리 전체 rank 정규화 → 9정의** 계산.
- **input**: `<outputs>/*/05_stage1_binding/binding_row.csv`
- **output**: `<out-dir>/binding_scores.csv`(9정의 + 원값), `<out-dir>/pocket_clusters.csv`(포켓/클러스터/포켓잔기)

## 4_finalize.py  (최종화, stdlib)
- **역할**: 9정의 중 하나를 골라 CASP 제출 `bind.txt` + 정의선택 이유 txt 작성.
- **input**: `<stage1-dir>/{binding_scores.csv, pocket_clusters.csv}`, `--definition`(예 prob_combined), `--group`
- **output**: `<out-dir>/<target>LG<group>.bind.txt`, `definition_choice.txt`

## 5_defs_viz.py  (분석: 비교차트, sci+matplotlib)
- **역할**: 9정의의 순위상관(Spearman)·상위 겹침(Jaccard)·>0.5 분포를 한 그림으로.
- **input**: `binding_scores.csv`
- **output**: 비교 대시보드 png

## 6_consensus_pocket.py  (분석: 콘센서스 포켓, stdlib)
- **역할**: 전 fragment 포켓 잔기를 집계해 접촉빈도 높은 콘센서스 포켓 도출.
- **input**: `pocket_clusters.csv` (pocket_residues 열)
- **output**: `consensus_pocket.txt` (잔기별 접촉빈도 + ChainResnum 한 줄)

## 러너 (오케스트레이터, bash)
- **run_stage1_frag.sh**: fragment 1개 = 1_rank_poses_boltz → core 0b~4 → 2_binding.
  - input: `<CID> [RUNS_ROOT] [SMI_DIR] [OUT_ROOT]` / output: `<OUT_ROOT>/<CID>/00~05`
- **run_stage1_all.sh**: 전 fragment 배치(idempotent skip) → run_stage1_frag 반복.
  - input: `<chunk> <n_chunks>` (GPU 병렬) / output: 각 fragment 폴더 + `logs/`
