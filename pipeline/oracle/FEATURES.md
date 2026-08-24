# oracle — Features (파일/값 정의)

## oracle_ligands.tsv  (정답 세트)
`ID, Name, SMILES, label`
- `label` = **1(양성=결합 억제제) / 0(음성=디코이)** — AUC 채점의 정답지.
- Name 규칙: `ORA_P#_이름`(양성), `ORA_D##_이름`(음성).

## <name>.yaml  (make_oracle_yaml 출력)
boltz2 native affinity 입력. 틀과 동일하되 `ligand.smiles`만 교체:
```yaml
version: 1
sequences:
  - protein: {id: A, sequence: <BFT1>, msa: <재사용 a3m>}
  - ligand:  {id: B, smiles: '<oracle SMILES>'}
properties:
  - affinity: {binder: B}   # ← affinity 예측을 켜는 부분
```

## binding_row.csv  (분자당, stage1 2_binding 출력)
stage1과 동일 컬럼 (boltz/cnn/cnnaff/aff/hac + 포켓/클러스터). → [stage1/FEATURES.md](../stage1/FEATURES.md) 참고.

## AUC  (s1_oracle_auc 출력)
- 정의별 값 = **P(양성점수 > 음성점수)**, 36쌍(양성3×음성12) 비교.
- `1.0` 완벽, `0.5` 랜덤, `<0.5` 거꾸로(나쁨).
- 9정의 각각 AUC → 최고값 = 제출 정의 후보.

## 결과 예시 (L01)
```
prob_boltz    0.833
prob_vina     0.806
prob_combined 0.806
prob_cons3    0.736
prob_LE_vina  0.722
prob_cnn      0.611
prob_cnnaff   0.583
prob_gnina    0.583
prob_LE_caf   0.361   ← 랜덤 이하
```
