# L01 전용 파이프라인 — 스크립트 설명 & 실행 방법

CASP17 타겟 **L01 (BFT1, 아연 metalloprotease)** 은 "1 단백질 : **다수** 리간드" 문제입니다.
대회 중반 요건 변경으로 **① 유기 리간드가 구조당 1~4 copy(stoichiometry)**, **② 촉매 아연(Zn) 1개 필수**가 되면서,
"리간드 1개"를 전제한 일반 파이프라인([`../../pipeline/`](../../pipeline/))을 그대로 쓸 수 없어 **stage-2 전용 스크립트(`01_inputs`~`07_validation`)** 를 별도로 구성했습니다.

- **재사용**: 포켓 발견(스텝 0b·1·2)은 일반 파이프라인 `run_pipeline.py`를 그대로 호출 → `04_pipeline`에서 오케스트레이션.
- **우회**: 도킹·배향 클러스터(스텝 3·4)와 신뢰도(5c)는 multi-copy에서 퇴화하므로 사용하지 않고, **구조 단위 선정 + 실험 기반 검증**으로 대체.

> 경로·계정은 특정 환경을 지우고 `/path/to/...`·`$CASP17`·`USERNAME` 으로 일반화했습니다. 실제 실행은 본인 환경에 맞게 채워야 합니다.

---

## 실행 방법 (stage-2 흐름)

각 binder 컨테이너는 `<binders>/<id>/{inputs,results,consensus_s2}` 구조입니다. 아래 순서대로 진행합니다.

```bash
# 0) 환경 (클러스터 경로/계정 채운 뒤)
source scripts/env_setup.sh          # micromamba env(boltz2), gnina/p2rank/PoseBusters 경로

# 1) 입력 생성 — 유기 N copy(stoichiometry) + Zn 1개 반영 (모델별)
python 01_inputs/stage2_pt2_inputs.py  --stage2-csv L01.smiles.stage2.csv --existing-root <binders> ...
python 01_inputs/stage2_bt2_inputs.py  ...
python 01_inputs/stage2_af3_inputs.py  ...

# 2) 추론 — 모델별 co-folding (SLURM 배열, modelSeeds 1..6 → 6 seed x 5 sample = 30 pose/모델)
sbatch --array=0-5%3 02_inference/run_af3_binders_stage2.sh
sbatch             02_inference/run_bt2_binders_stage2.sh
sbatch             02_inference/run_protenix_binders_stage2.sh
python 02_inference/normalize_bt2_unified.py ...     # bt2 출력 레이아웃 통일(collect가 읽게)

# 3~4) 합본 + 포켓 발견 + p2rank 검증 (오케스트레이터가 03/일반스텝0b·1·2/p2rank 순차 실행)
MODELS=pt2,bt2,af3 bash 04_pipeline/stage2_run_pockets.sh          # 전체
#                  bash 04_pipeline/stage2_run_pockets.sh L010016   # 지정 binder만(검증용)

# 5) 최종 pose 선정 (구조 단위)
python 05_select/stage2_select_multicopy.py --cons <binders> --top 5
#  → <binders>/stage2_mc_selection.csv (binder별 MODEL1) + stage2_mc_candidates.csv

# 6) CASP LG 제출 포맷 (수용체 + 유기 N copy + Zn)
CASP_GROUP=YYY CASP_AUTHOR=CHANGE-ME python 06_submission/make_lg_percomplex.py

# 7) 물리 검증 → 정제 → 실험 기반 검증
python 07_validation/prep_validity_stage2.py ...     # PoseBusters 입력(copy별 SDF) + manifest
#      (PoseBusters 검사는 일반 컴포넌트 boost/4c_posebusters.py 가 manifest로 수행)
python 07_validation/refine_stage2.py ...            # 실패 copy만 gnina 국소최소화
python 07_validation/refine_mmff_stage2.py ...       # 남은 covalent 결함 MMFF 교정
python 07_validation/find_clean_alt.py --binder L01xxxx ...   # 그래도 안 되면 대체 pose 교체
python 07_validation/exosite_overlap.py ...          # 실험 exosite 약물자리와 거리
python 07_validation/pose_zn_distance.py ...         # 촉매 Zn 거리(활성부위 vs exosite)
```

---

## 단계별 스크립트 설명

| 단계 | 스크립트 | 하는 일 | 주요 값·설정 | L01 특화 (일반과 차이) |
|---|---|---|---|---|
| **01_inputs**<br>입력 생성 | `stage2_af3_inputs.py`<br>`stage2_bt2_inputs.py`<br>`stage2_pt2_inputs.py` | 모델별 co-folding 입력(JSON/YAML) 재생성. 기존 binder는 단백질·MSA만 재사용하고 리간드를 stage-2 값으로 교체 | 유기 리간드 `count=N`(stoichiometry)<br>이온 `ZN count=1` | **N copy + Zn 추가**가 핵심 변경.<br>신규 binder(L010123)는 기존 input 없어 템플릿 수용체 재사용 |
| **02_inference**<br>추론 | `run_af3_binders_stage2.sh`<br>`run_bt2_binders_stage2.sh`<br>`run_protenix_binders_stage2.sh`<br>`normalize_bt2_unified.py` | 각 모델 co-folding을 SLURM 배열로 실행. `normalize_bt2_unified.py`는 bt2 출력 폴더 구조를 통일 포맷으로 정규화 | `modelSeeds 1..6`<br>→ 6 seed × 5 sample<br>= **30 pose/모델** | stage-1(Boltz2 단독)과 달리 **3모델(af3·bt2·pt2)** 재실행 — 기존 pose엔 copy·Zn이 없었기 때문 |
| **03_collect**<br>합본·랭킹 | `collect_consensus.py` | 3(또는 4)모델의 native 출력을 binder당 `master_table.csv` **한 장**으로 합본. 원본 무수정(경로만 참조) | `ligand_iptm` = chain_pair_iptm off-diag 평균<br>정렬키 = (clash없음 → ligand_iptm → plddt → −gpde)<br>`--max-seeds 6` | 모델 간 **합의(consensus)** 를 한 표에서 랭킹. bt2는 키 매핑(plddt←complex_plddt 등) |
| **04_pipeline**<br>포켓 발견·검증 | `stage2_run_pockets.sh`<br>`validate_pockets_p2rank.py` | 오케스트레이터: 03 합본 → config 작성 → 일반 파이프라인 **스텝 0b(cache)/1(protein_cluster)/2(pocket_candidates)** 호출 → per-pocket p2rank 교차검증 | greedy pocket **8Å**, top **10**<br>p2rank cavity **≤6.0Å** | **copy별 무게중심을 각각 point로**(구조당 N점) → 위치 클러스터링으로 진짜 자리 복원.<br>스텝 3·4는 **호출 안 함**(multi-copy 퇴화) |
| **05_select**<br>최종 선정 | `stage2_select_multicopy.py` | 포켓별 도킹/배향클러스터 대신 **완성 구조를 순위매김** → MODEL1 + top 후보 | 순위키(내림차순) =<br>① clash 없음 ② coverage ③ n_in_top ④ ligand_iptm<br>(N=copy수, topN=dominant pocket) | co-folding이 N copy를 이미 배치 → **"완성된 답"을 고르기**만.<br>일반 파이프라인의 composite 종합점수 대신 **구조 단위** |
| **06_submission**<br>제출 포맷 | `make_lg_percomplex.py` | binder마다 CASP **LG**(Example 6.1) 생성: 수용체 PDB + 유기 리간드 N copy(각 MDL) + Zn 이온 | Zn = 단일원자 MDL + `M CHG +2`<br>env `CASP_GROUP`·`CASP_AUTHOR` | per-complex, **multi-copy + Zn** 포맷. copy마다 `LIGAND N <HTX>` |
| **07_validation**<br>물리·정제 | `prep_validity_stage2.py` | 선정 pose → PoseBusters 입력(**copy별** SDF + 단백질 PDB) + manifest | config="dock"(정답 불필요) | copy 1개가 아니라 **유기 copy 전부**를 각각 SDF로 |
| | `refine_stage2.py` | PoseBusters 실패 copy만 gnina 국소최소화로 정제(실패 copy heavy atom만 cif에 덮어씀) | gnina `--minimize_iters 10`<br>drift **>2.0Å → 원본 유지** | 예측 자체가 틀린 게 아니라 국소 결함만 고침. 나머지 copy·단백질·Zn 불변 |
| | `refine_mmff_stage2.py` | gnina 후에도 남은 covalent 결함(ring/bond)을 RDKit MMFF 제약최소화 | heavy atom 원위치 **±0.5Å(fc=100)**<br>drift **>1.0Å → 원본** | **L01 전용** 2차 정제. 위치 거의 유지하며 결합길이·각·고리 이상화 |
| | `find_clean_alt.py` | gnina·MMFF로 못 고친 clash binder: 후보 중 (전 copy valid + exosite 근접) 대체 pose 첫 통과 채택 | 대체 조건: 전 copy valid **AND** exosite `<8Å` | **L01 전용** 최후 수단(교체). 후보는 `stage2_mc_candidates.csv`(--top 크게) |
| | `exosite_overlap.py` | 선정 pose가 실험 약물결합 exosite와 겹치는지 정량. 실험구조를 우리 프레임으로 정렬 → 리간드↔실험약물 최소거리 | ref 7POL·7POO·7POQ·7POU<br>near **<8Å** (실측 78/80 <3Å) | **정답 없이 신뢰도 확보**의 핵심. co-folding이 exosite에 강하게 합의 |
| | `pose_zn_distance.py` | 각 유기 copy가 촉매 Zn에서 얼마나 떨어졌나 → "활성부위 vs exosite" 판정 | template 불필요(cif에 Zn 포함) | Zn metalloprotease 특성 활용. 70/80이 Zn에서 >15Å → exosite 수렴 확인 |

> **각주**: `refine_mmff_stage2.py`(2차 MMFF 정제)와 `find_clean_alt.py`(대체 pose 교체)는 **L01 전용**입니다. 일반 1:1 파이프라인의 `5b_refine.py`는 gnina 단독이고, 대체 pose 역할은 선정(스텝5) 단계의 medoid 강건화가 담당합니다.

---

## 재사용하는 일반 컴포넌트

L01 전용 스크립트는 다음 일반 파이프라인 요소를 그대로 호출합니다.

| 컴포넌트 | 위치 | L01에서의 역할 |
|---|---|---|
| 포켓 발견 (스텝 0b·1·2) | `pipeline/pipelines/run_pipeline.py` | `04_pipeline/stage2_run_pockets.sh`가 `--only cache/protein_cluster/pocket_candidates`로 호출 |
| PoseBusters 검사 | `pipeline/pipelines/boost/4c_posebusters.py` | `07_validation/prep_validity_stage2.py`가 만든 manifest를 소비해 물리 유효성 판정 |
| 결합 확률(Task A) 오라클 | `pipeline/oracle/` | stage-1에서 "어떤 결합확률 정의가 binder를 가장 잘 구분하나"를 AUC로 채점(`s1_oracle_auc.py`) → 제출용 정의 선택. **stage-2 pose 파이프라인과는 별개 트랙** |

---

## 일반(1:1) 파이프라인과의 차이 요약

| 항목 | 일반 1 단백질:1 리간드 | **L01 (multi-copy + Zn)** |
|---|---|---|
| point 생성 | 구조당 1점(리간드 무게중심) | **구조당 N점**(copy별 centroid) |
| 스텝 3 (pocket 검증) | 사용(p2rank+gnina+템플릿) | **우회** → `exosite_overlap.py`로 대체 |
| 스텝 4 (배향 클러스터) | 사용(SC-RMSD) | **우회**(concat 원자수 ≠ SMILES → SC-RMSD 퇴화) |
| 스텝 5 (선정) | composite 종합점수 → top5 | **구조 단위**(clash→coverage→n_in_top→ligand_iptm) → MODEL1 |
| 스텝 5c (신뢰도) | 사용(HIGH/MEDIUM/LOW) | **스킵** → exosite·Zn거리·PoseBusters로 대체 |
| 정제 | gnina 단독 | gnina → **MMFF → 대체 pose 교체**(L01 전용 2·3차) |
| 제출 | 리간드 1개 LG | 유기 **N copy + Zn** per-complex LG |

---

*상세 딥다이브: [`method_pocket_multicopy.md`](method_pocket_multicopy.md)(copy별 pocket 정의) · 전체 결과 해석은 [`../README.md`](../README.md).*
