# stage2 — Features (생성 파일의 컬럼 정의)

## selection_summary.csv  (core 5_final_select 출력, fragment당)
`model, pocket_id, ligand_cluster_id, size, iptm, lscore, gnina_affinity, cnn_score, composite, pocket_pass, posebusters_valid, source_cif`
- `model` model_1~5 (composite 순) / `size` 그 포즈 클러스터 크기(합의)
- `iptm` interface 신뢰도 / `gnina_affinity` 에너지(음수) / `cnn_score` gnina CNN
- `composite` 종합점수 = 합의·iptm·affinity 가중합 (선정 순위)
- `pocket_pass` 그 포켓 검증 통과 여부 / `source_cif` 원본 pose 경로

## poses.csv  (1_aggregate 출력, 전 리간드 통합)
`cid, model, pocket_id, ligand_cluster_id, size, iptm, gnina_affinity, cnn_score, composite, pocket_pass, source_cif`
- selection_summary를 전 fragment 모은 것 + `cid`. 리간드당 최대 5행.
- **주의**: `gnina_affinity`는 **포즈마다**(model별) 다른 구조의 값. Stage1 `aff`(대표 1개)와는 model_1만 일치.

## LG 파일  (core 9_make_casp_lg 출력)
`<cid>LG_model{N}.txt` (개별) + `<cid>LG_all_models.txt` (합본). CASP LG 포맷:
```
PFRMAT LG
TARGET <cid>
AUTHOR <group>
METHOD ...
MODEL 1
REMARK ...
ATOM ... (수용체 = PDB 형식)
LSCORE <0~1>
LIGAND <id> <name>
<리간드 mol block = MDL 형식>
END
```

## submit 구조  (1_aggregate 출력)
```
final/stage2/
  poses.csv
  submit/<cid>/<cid>LG.txt   (제출용 포즈 파일)
```
- CASP 정식 제출은 `./L01/L0xxxxxLG<group>_N` 개별 파일 → tar → `L01LG<group>.tgz`.
