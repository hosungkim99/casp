# oracle — Usage

## 사전 준비물 (필수)
1. **정답 세트** — `oracle_ligands.tsv` (ID, Name, SMILES, label=1양성/0음성) + `smi/<name>.smi`
2. **틀 YAML** — 기존 fragment의 boltz affinity YAML 1개 (BFT1+MSA+affinity 속성). 리간드만 교체됨.
3. **환경** — boltz env (job.sh의 `micromamba run -p $CONDA_ENVS_PATH/boltz2`) + `source env_setup.sh`(파이프라인 단계용 p2rank/gnina)
4. **GPU**

## 전체 순서
```bash
ORA=/path/to/casp17-ligand/.../users/USERNAME/oracle

# 1) 정답 세트 준비 (tsv + smi) — 최초 1회 (분자 SMILES는 검증된 것)
# 2) YAML 생성
python make_oracle_yaml.py \
  --template-yaml /path/to/casp17-ligand/.../targets/L01/affinity_yaml/L010001.yaml \
  --ligands $ORA/oracle_ligands.tsv --out-dir $ORA/affinity_yaml

# 3) boltz cofold (GPU 3장 예시)
for i in 0 1 2; do sbatch run_oracle_boltz.sh $i 3; done

# 4) Stage1 파이프라인 (15분자, stage1 러너 재사용)
source /path/to/casp17-ligand/.../scripts/env_setup.sh
SC=/path/to/casp17-ligand/.../users/USERNAME/pipeline
for y in $ORA/affinity_yaml/*.yaml; do
  name=$(basename "$y" .yaml)
  bash $SC/stage1/run_stage1_frag.sh "$name" "$ORA/runs" "$ORA/smi" "$ORA/outputs"
done

# 5) AUC 채점 (= 최종 답)
python s1_oracle_auc.py --outputs $ORA/outputs --ligands $ORA/oracle_ligands.tsv
```

## 완료 확인
```bash
find $ORA/runs -name "*_model_0.cif" | wc -l          # boltz 완료 (분자×seed)
find $ORA/runs -name "affinity_ORA_*.json" | wc -l     # affinity 나왔나 (prob_boltz 원천)
ls $ORA/outputs/*/05_stage1_binding/binding_row.csv | wc -l   # 파이프라인 완료 (=분자 수)
```

## 병목
- **boltz cofold(3단계)** 가 병목이지만, BFT1 MSA 재사용 + boltz 22초/건이라 15분자×6seed도 ~30분.

## 대안 (빠름/덜 정확)
- cofold 없이 Site B 직접 도킹 → boltz 정의(prob_boltz)는 테스트 못 함 (8/9만).
