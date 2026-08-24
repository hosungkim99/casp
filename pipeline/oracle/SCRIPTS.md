# oracle — Scripts

(실행법 [USAGE.md](USAGE.md), 컬럼 정의 [FEATURES.md](FEATURES.md))

## make_oracle_yaml.py  (입력 생성, stdlib)
- **역할**: 기존 fragment의 boltz affinity YAML을 틀로, **리간드 SMILES 한 줄만** oracle 분자로 교체.
- BFT1 단백질·MSA·affinity 속성은 그대로 재사용(동일 수용체 → MSA 재계산 불필요).
- **input**: `--template-yaml <기존 L0xxxxx.yaml>`, `--ligands oracle_ligands.tsv`(Name, SMILES)
- **output**: `<out-dir>/<name>.yaml` (분자마다 1개)

## run_oracle_boltz.sh  (boltz cofold, sbatch)
- **역할**: oracle 분자 YAML을 boltz2로 예측 (job.sh 방식 재현, affinity 포함).
- fragment 1209개와 **동일 레이아웃**으로 저장 → stage1 파이프라인이 그대로 읽음.
- **input**: `<chunk> <n_chunks>` (+ env `ORA`=작업폴더, `NSEED`=6, `SAMPLES`=5)
- **output**: `<ORA>/runs/<name>/seed_*/boltz_results_*/predictions/<name>/` (pose cif + confidence + affinity json)

## s1_oracle_auc.py  (AUC 채점, stdlib)
- **역할**: oracle 분자들의 `binding_row.csv` 를 모아 **9정의 각각의 AUC**(양성 vs 음성) 계산.
- AUC = Mann-Whitney U (sklearn 불필요).
- **input**: `--outputs <ORA>/outputs`(각 분자 05_stage1_binding), `--ligands oracle_ligands.tsv`(Name, label)
- **output**: 콘솔에 9정의 AUC 순위 + 최고 정의 (텍스트)

## 재사용하는 것
- **stage1/run_stage1_frag.sh**: oracle 분자에도 그대로 (RUNS_ROOT=oracle/runs, SMI_DIR=oracle/smi, OUT_ROOT=oracle/outputs).
- 이름이 `L*`가 아니어도(ORA_*) CID를 인자로 직접 받아 작동.
