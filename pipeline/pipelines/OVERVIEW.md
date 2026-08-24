# CASP 리간드 파이프라인 — 전체 지도 (처음 보는 사람용)

> 이 문서 하나만 읽으면 "이 파이프라인이 뭘 하고, 어떤 스크립트가 어떤 역할이며,
> 어디부터 봐야 하는지"를 알 수 있게 쓴 지도입니다.
> 스크립트별 요약은 [SCRIPTS.md](SCRIPTS.md), 설치/실행은 [USAGE.md](USAGE.md),
> 출력 컬럼 정의는 [FEATURES.md](FEATURES.md),
> 단계별 결과 해석(무엇을 보고 믿을지)은 [INTERPRETATION.md](INTERPRETATION.md)에 있습니다.

## 폴더 구조 (2026-07 재구성)
스크립트는 성격별 하위폴더로 나뉘고, 공용 라이브러리는 상위 `pipeline/common/`에 있음.
import는 각 스크립트 상단 **부트스트랩**이 `common/`을 자동으로 찾아 처리(깊이 무관).
```
pipeline/
├─ common/                 파서·기하·gnina (import 전용)  → ../common/OVERVIEW.md
└─ pipelines/
   ├─ run_pipeline.py      지휘자 (config 읽어 순서대로 자동 실행)
   ├─ core/     0~6,9      메인 흐름 (항상 자동)
   ├─ boost/    4b,4c,5b,5c 보강 (config 플래그로 삽입)
   ├─ analysis/ 7,8,10,11  검증·시각화·보고
   └─ ml/       extract/select/train  ML 실험 (수동)
```
> 아래 본문은 스크립트를 번호로만 부르지만, 실제 경로는 위 하위폴더 기준(`core/2_pocket_candidates.py` 등).

---

## 0. 30초 요약 — 이게 대체 뭘 하는 건가

- **문제**: CASP 대회는 단백질 서열 + 리간드(SMILES) + Task(P=포즈 예측 / A=친화도)를 주고
  "이 리간드가 단백질 어디에 어떻게 붙는지" 3일 안에 맞히라고 한다.
- **재료**: 팀이 먼저 co-folding 모델(AF3, Boltz2 등)을 여러 개·여러 seed로 돌려
  **수천 개의 예측 구조**(단백질+리간드가 붙은 모습)를 만들어 둔다.
- **이 파이프라인이 하는 일**: 그 수천 개를 **자동으로 걸러 제출용 5개**로 압축한다.
- **핵심 아이디어 (consensus, 합의)**: 여러 모델이 **같은 자리에 리간드를 놓으면** 그게 정답일 확률이 높다.
  그래서 raw 점수 하나로 줄 세우지 않고, **"많이 모인 자리 = 좋은 자리"** 로 판단한다.
- **처리 순서 (계층형)**:
  **① 자리(포켓)를 먼저 정한다 → ② 그 자리 안에서 리간드 방향을 정한다 → ③ 물리적으로 말이 되는지 검증 → ④ 최종 5개 선정.**

```
수천 개 예측  ──►  포켓 후보  ──►  검증 통과 포켓  ──►  포켓별 방향  ──►  최종 5개  ──►  제출 파일
 (0단계)          (2단계)        (3단계)            (4단계)          (5단계)       (9단계)
```

---

## 1. 라인 구분 — 한눈에 보기

파이프라인은 성격이 다른 5개 라인 + 공용 라이브러리로 나뉜다.

| 라인 | 스크립트 | 자동 실행? | 한 줄 목적 |
|------|----------|-----------|-----------|
| **① 메인 흐름** | `0`~`6`, `9` | ✅ 항상 자동 | 수천 예측 → 제출용 5개 (핵심 경로) |
| **② 보강** | `4b`,`4c`,`5b`,`5c` | ⚙️ config 플래그 켤 때만 자동 삽입 | 물리 검증·정제·신뢰도 (품질 보강) |
| **③ 검증/비교** | `7`, `10` | 🧪 실험구조/템플릿 있을 때만 | 우리 예측이 실험과 맞는지 대조 |
| **④ 시각화·보고** | `8`, `11` | ✅ 자동 | 그림(`8`) + 정리문서 초안(`11`; 사실만 자동, 판단은 사람이 ✍️) |
| **⑤ ML 실험** | `extract_features`, `select_score_features`, `train_xgb` | ❌ 수동·오프라인 | 여러 타겟 모아 학습용 데이터/모델 (연구용) |
| **공용 라이브러리** | `complex_io.py`, `geom.py`, `scoring.py` | (import 전용) | 모든 스텝이 공유하는 파서·기하·gnina 함수 |
| **legacy/** | `cluster_poses` 등 | 🚫 미사용 | 단일레벨 구버전 (참고용 보관) |

**지휘자**: `run_pipeline.py` 가 config 파일 하나를 읽어 ①②③④를 **순서대로 자동 호출**한다.
(⑤ ML 라인은 파이프라인과 분리된 별도 도구다.)

> 실행 환경 표기: **sci** = boltz2 env(`numpy`+`gemmi`+`rdkit`+`matplotlib`) / **stdlib** = 순수 python /
> **pb** = posebusters env. gnina·p2rank는 singularity/외부 바이너리로 호출.

---

## 2. ① 메인 흐름 (필수 경로) — 상세

각 스텝을 **[받는 것] → [하는 일] → [내놓는 것] → 💡왜** 로 정리.

### `0_rank_poses.py` — 수집·랭킹 (stdlib)
- **[받음]** 팀 결과 폴더 `<model>/seed_*/sample_*/summary.json`(+`model.cif`)
- **[함]** 흩어진 수천 예측을 한 표로 모으고, **interface 신뢰도**(chain_pair_iptm off-diagonal 평균)로 정렬
- **[냄]** `00_collect/master_table.csv` (rank, model, iptm, cif 경로 …), `progress.txt`
- 💡 이후 모든 스텝의 공통 입력. "어떤 예측이 있고 어디 파일에 있나"의 목록.

### `0b_cache_geometry.py` — 기하 캐시 (성능, sci)
- **[받음]** master_table
- **[함]** 모든 cif를 **딱 한 번** 파싱해 (단백질 Cα + 리간드 좌표)만 뽑아 `geom_cache.pkl`로 저장.
- **[냄]** `00_collect/geom_cache.pkl`
- 💡 원래 `1·2·6·8`이 같은 수천 cif를 **각자 재파싱**(gemmi 파싱이 병목)하던 것을 1회로 통일 → 파싱 4회→1회.
  캐시가 없으면 각 스텝이 자동으로 즉시 파싱(하위호환). 좌표만 담아 원본보다 작다.

### `1_protein_cluster.py` — 단백질 형태 분류 (Step 0, sci)
- **[받음]** master_table
- **[함]** 모델마다 단백질 형태가 다를 수 있어(열림/닫힘 등), Cα를 정렬해 **PCA→형태 그룹**으로 라벨링.
  동시에 **정렬 기준 구조(reference)** 를 "PC 공간의 중심 구조"로 선정(=치우치지 않은 대표).
- **[냄]** `01_protein_clusters/protein_clusters.csv`, `reference.txt`, `reps/conf_*.cif`
- 💡 뒤 스텝들이 **모두 같은 기준 구조**에 정렬하도록 통일. (형태 라벨 자체는 참고용)

### `2_pocket_candidates.py` — 포켓 후보 찾기 (Step 1, sci)
- **[받음]** master_table (+ reference.txt)
- **[함]** 모든 예측의 **리간드 무게중심**을 기준 구조 좌표계로 옮긴 뒤, 가까운 것끼리 묶어(그리디 클러스터)
  **가장 붐비는 자리 top 10** 을 포켓 후보로. 작은 자리도 일단 유지.
- **[냄]** `02_pocket_candidates/pocket_candidates.csv`(자리별 크기·참여 모델), `members.csv`(어느 예측이 어느 자리)
- 💡 "리간드가 실제로 자주 놓이는 자리 = 진짜 결합부위 후보". `size`(모인 수)·`n_models`(참여 모델수)가 신뢰의 핵심.

### `3_pocket_validate.py` — 포켓 검증·압축 (Step 2, sci)
- **[받음]** pocket_candidates (+ members.csv)
- **[함]** 후보 자리 대표에 **p2rank**(진짜 캐비티?) + **gnina**(에너지?)를 돌려 `pass/fail` 판정.
  gnina 1회가 **Vina affinity**(`gnina_affinity`)와 **Gnina CNN**(`cnn_score`)을 동시 산출.
  추가로 포켓 **멤버 top-K**(`pocket_members_topk`, 기본 5)를 함께 채점해 **mean/std affinity·cnn** 병기.
- **[냄]** `03_pocket_validation/pocket_validation.csv` (+ `p2rank/`, `rescore/`)
  config `templates` 주면 **template_dist**(템플릿 리간드 ↔ 포켓 center 거리)를 병기하고,
  `template_dist ≤ template_pass_dist`(기본 8Å)인 자리는 **pass로 자동 승격**(pass = 물리 OR 템플릿 지지).
- 💡 "붐빈다고 다 진짜는 아니다" — 물리/구조로 가짜 자리를 거른다. 포켓이 불명확하면(큰 클러스터 없음)
  **평균 score가 낮은(좋은) 포켓**을 참고하고, 리간드 있는 템플릿이 있으면 그 자리를 우선(회의 1·1-1).

### `4_ligand_cluster.py` — 포켓별 리간드 방향 (Step 3, sci)
- **[받음]** members + 통과 포켓 + ligand.tsv(SMILES)
- **[함]** 통과한 자리 **안에서**, 리간드 방향/자세를 **SC-RMSD(대칭 보정 RMSD)** 로 다시 클러스터링 → 자리별 대표 포즈들.
  각 클러스터의 **갯수(size)·평균·분산**(`rmsd_mean`/`rmsd_std`)도 기록(회의 "갯수·평균·분산 분석").
- **[냄]** `04_ligand_clusters/ligand_clusters.csv` (+ cluster_members.csv)
- 💡 같은 자리라도 리간드가 뒤집혀 놓일 수 있다. 방향까지 합의된 대표를 뽑는다(응집도로 신뢰 판단).

### `5_final_select.py` — 최종 선정 (Step 4, sci)
- **[받음]** ligand_clusters + pocket_validation (+ posebusters 있으면)
- **[함]** 각 후보에 gnina를 매기고 **종합점수 = 합의(크기) + Boltz2(iptm) + Vina(affinity) + Gnina(CNNscore)**
  4항 가중합으로 상위 ≤5개 선정(회의 "세 값 조합"). Task별 가중치 분기(**P**=포즈 중심 / **PA**=affinity↑).
  선정 클러스터의 갯수·평균·분산(`rmsd_mean`/`rmsd_std`)도 요약에 병기. PoseBusters 무효 포즈는 제외.
- **[냄]** `05_final/{selection_summary.csv, model_1~5.cif, SELECTION_RATIONALE.md}`
- 💡 파이프라인의 결론. "왜 이 5개인가"를 점수와 함께 남긴다.

### `6_contact_residues.py` — 결합 잔기 (검증, sci)
- **[받음]** master_table (+ reference.txt)
- **[함]** 1순위 포켓 멤버들에서 리간드가 **실제로 닿는 단백질 잔기**를 4Å 이내로 세어 빈도(%) 계산.
- **[냄]** `06_validation/contact_residues.csv`
- 💡 "우리가 예측한 결합부위"의 정의. 빈도 높은 잔기가 많고 고를수록 한 곳으로 수렴(좋음).

### `9_make_casp_lg.py` — 제출 파일 생성 (sci)
- **[받음]** (정제 여부에 따라) selection_summary + ligand.tsv
- **[함]** 최종 5개를 **CASP LG 포맷**(수용체 PDB + 리간드 mol block)으로 변환. 리간드 connectivity는 SMILES에서, 좌표는 예측에서.
- **[냄]** `08_casp_lg/<TARGET>LG_model1~5.txt`, `_all_models.txt`
- 💡 대회에 실제로 제출하는 파일. 마지막 관문.

---

## 3. ② 보강 라인 (옵션 — config로 켜면 자동 삽입)

메인 흐름 품질을 높이는 추가 검증. `run_pipeline.py`가 조건에 맞으면 알맞은 위치에 끼워 넣는다.

| 스크립트 | 켜는 법 | 삽입 위치 | 역할 |
|----------|---------|-----------|------|
| `4b_prep_validity.py` + `4c_posebusters.py` | config `python_pb=` 설정 | final_select **앞** | **①물리 유효성**: 클러스터당 상위 K개 멤버(중심성순)의 SDF/PDB 준비(4b) → PoseBusters 검사(4c). final_select가 무효 제외 + **medoid 강건화**(대표 실패 시 valid 멤버로 교체) |
| `5b_refine.py` | config `refine=true` | final_select **뒤** | **②포즈 정제**: gnina `--minimize`로 국소 최소화. drift가 크면 원본 유지(안전장치). 제출은 정제본 사용 |
| `5c_confidence.py` | 항상 (맨 끝) | 파이프라인 **끝** | **③신뢰도 판정**: 모든 출력을 모아 **HIGH/MEDIUM/LOW** 로 "이 자동선택을 믿어도 되나" 판정 |

> **가장 먼저 볼 결과**: `05c_confidence/confidence_report.md` — HIGH면 그대로 제출, MEDIUM/LOW면 사람 검토.

---

## 4. ③ 검증/비교 라인 (실험 구조·템플릿이 있을 때)

우리 예측이 **실제/유사 구조와 맞는지** 정량 대조. 대회 중엔 유사 구조(템플릿), 대회 후엔 정답으로 쓴다.

### `7_fragment_compare.py` — 실험 fragment 비교 (sci)
- config `experimental_cif=` 있으면 자동. 실험 구조를 model_1에 Cα 정렬(번호 offset 자동) →
  실험 리간드 접촉잔기 ↔ 우리 예측 접촉잔기 **공유 개수**. 오버레이 PDB도 저장.
- ⚠️ **단일 앵커**: model_1(단백질 medoid) 하나에만 맞춰 비교 → 흩어진/다중포켓 타겟엔 대표성이 약함.

### `10_pose_vs_template.py` — 후보별 템플릿 대조 (sci)
- `7`의 단일 앵커 한계를 보완. **최종 model_1~5 각각**을 여러 템플릿과 개별 비교:
  - `overlap_pct` = 우리 리간드 원자 중 실험 리간드 4Å 이내 비율(분모=우리 원자 → 과대평가 방지)
  - `jaccard` = 접촉잔기 겹침 비율, `convergence` = 그 후보 클러스터의 수렴도
- 💡 "L1이 아니라 L3가 실험과 더 맞는" 다중포켓 타겟에서 **진짜 상위 후보**를 골라준다.
- ✅ config `templates=` 지정 시 `run_pipeline.py`가 **자동 실행**(final_select/refine 뒤). 정제 켜져 있으면 정제본을 대상으로.
- ⚠️ 템플릿은 **리간드 포함 RCSB 원본**이어야 함(inputs의 리간드 제거본 불가).

---

## 5. ④ 시각화 라인

### `8_visualize.py` — 종합 그림 (sci)
- **[냄]** `07_viz/visualization.png` (4패널) + `contact_fingerprint.png`
- 4패널: (1) 리간드 PCA by 모델 (2) PCA by 클러스터 (3) 모델별 iptm 분포(=raw 비교 금지 근거) (4) 상위 클러스터의 모델 구성
- 💡 **한 덩어리에 여러 모델 색이 겹치면 = 강한 합의**(신뢰↑). 지문 그래프에서 빨강 = 실험과 공유 잔기.

---

## 6. ⑤ ML 실험 라인 (파이프라인과 분리된 오프라인 연구용)

여러 타겟의 파이프라인 출력을 모아 **"어떤 포즈가 정답이 될지"를 학습**하려는 실험. 수동 실행.

| 스크립트 | 역할 |
|----------|------|
| `extract_features.py` | 한 타겟의 출력 CSV들을 **포즈 단위로 join** → `features.csv`. `--truth` 주면 RMSD-to-정답 + label(≤2Å=1)까지. `--append`로 여러 타겟 누적 |
| `select_score_features.py` | `features.csv`에서 **핵심 피처만** 추려 `features_core.csv` |
| `train_xgb.py` | features.csv로 **XGBoost** 학습, AUC + 피처 중요도 출력. 타겟 단위 train/test 분리(누수 방지) |

> ⚠️ 한계(파일 주석에 명시): 포켓레벨 gnina는 포켓 대표값이 포켓 내 모든 포즈에 브로드캐스트됨(포즈 간 변별 X).

---

## 7. 공용 라이브러리 (import 전용, 단독 실행 X)

상위 `pipeline/common/`에 있음. 각 스텝 상단 부트스트랩이 자동으로 찾으므로 하위폴더 어디에 있어도 됨.
(상세: [../common/OVERVIEW.md](../common/OVERVIEW.md))

### `complex_io.py` — 단백질-리간드 파서 (gemmi)
타겟 종류·체인 수·리간드 개수에 무관하게 동작하는 공용 파서.
- `parse_complex` (Cα + 리간드 검출), `ligand_concat`/`ligand_centroid`, `composition`/`smiles_composition` (원소 조성 매칭),
  `match_ligand_to_smiles`, `read_ligand_tsv`, `write_ligands_pdb`, `reference_cif`(기준 구조 선택)
- ⚠️ 한계: 동일 조성 리간드가 여러 copy면 seqid 순서로 대응(진짜 교환대칭 최적배정 X).

### `geom.py` — 기하/정렬 (numpy; SC-RMSD는 rdkit)
- `kabsch`(회전+이동), `align_to_ref`(Cα 정렬 R,t), `apply_rt`,
  `ligand_automorphisms`(SMILES 대칭 순열), `sc_rmsd`(대칭보정 RMSD)
- ⚠️ 한계: 자기동형 순열은 원자 순서가 SMILES 순서와 같을 때만 정확, 아니면 plain RMSD 폴백.

### `scoring.py` — gnina 점수 (subprocess; 3·5 공유)
- `score_only`/`score_with_fallback`(GPU→CPU 폴백), `run_jobs`(opt-in 병렬 + GPU 라운드로빈), `parse_gpus`
- 💡 원래 `3`·`5`에 중복돼 있던 gnina 호출 코드를 한 곳으로 통일. 병렬은 기본 off(순차)로 공용 서버 배려.

### 캐시 접근 (complex_io 내)
- `load_geom_cache`, `cached_geometry` — `0b`가 만든 캐시를 읽어 재파싱을 건너뛴다(없으면 자동 폴백).

---

## 8. 스크립트별 역할 — 한 줄 사전 (전체 목록)

| 파일 | 라인 | 한 줄 |
|------|------|-------|
| `run_pipeline.py` | 지휘자 | config 읽어 스텝들을 순서대로 자동 실행 (idempotent, `--from`/`--only` 지원) |
| `run_pipeline.sh` | 지휘자 | 위를 slurm(sbatch)로 제출하는 래퍼 |
| `complex_io.py` | 공용 | 복합체 파서 (Cα·리간드·SMILES 매칭) + 기하 캐시 |
| `geom.py` | 공용 | 정렬/RMSD 기하 함수 |
| `scoring.py` | 공용 | gnina 점수 + opt-in 병렬 실행 (3·5 공유) |
| `0_rank_poses.py` | 메인 | 예측 수집 → master_table + 랭킹 |
| `0b_cache_geometry.py` | 메인(성능) | 전 cif 1회 파싱 → geom_cache (1·2·6·8 재파싱 제거) |
| `1_protein_cluster.py` | 메인 | 단백질 형태 클러스터 + 기준 구조 선정 |
| `2_pocket_candidates.py` | 메인 | 리간드 센트로이드 → 포켓 후보 top10 |
| `3_pocket_validate.py` | 메인 | p2rank+gnina로 포켓 검증·압축 |
| `4_ligand_cluster.py` | 메인 | 포켓 내 SC-RMSD로 리간드 방향 클러스터 |
| `5_final_select.py` | 메인 | 종합점수로 최종 5개 선정 |
| `6_contact_residues.py` | 메인 | 1순위 포켓 결합 잔기 빈도 |
| `9_make_casp_lg.py` | 메인 | 제출용 CASP LG 포맷 변환 |
| `4b_prep_validity.py` | 보강① | PoseBusters 입력(SDF/PDB) 준비 |
| `4c_posebusters.py` | 보강① | PoseBusters 물리 유효성 검사 |
| `5b_refine.py` | 보강② | gnina 국소 정제 (drift 안전장치) |
| `5c_confidence.py` | 보강③ | HIGH/MEDIUM/LOW 신뢰도 판정 |
| `7_fragment_compare.py` | 검증 | 실험 fragment와 접촉잔기 비교 (단일 앵커) |
| `10_pose_vs_template.py` | 검증 | 후보별 템플릿 대조 (overlap/jaccard; config `templates=`로 자동) |
| `8_visualize.py` | 시각화 | 4패널 + 접촉지문 PNG. templates 주면 **final pose별×템플릿별 접촉지문 격자**(`contact_fp_*.png`; 클러스터 빈도 x축, 템플릿 리간드 공유잔기 빨강)도 생성 |
| `11_make_summary.py` | 보고 | outputs CSV → 정리문서 초안 `SUMMARY_DRAFT.md`(사실 자동/판단은 ✍️ 자리) |
| `extract_features.py` | ML | 출력 CSV → 포즈 단위 features.csv |
| `select_score_features.py` | ML | 핵심 피처만 추림 |
| `train_xgb.py` | ML | XGBoost 학습/평가 |
| `legacy/*` | 미사용 | 단일레벨 구버전 보관 |

---

## 9. 실행법 / 되돌아가기 (요약)

```bash
# 1) config 작성 (경로 4개 + task 만 채우면 됨)
cp config/TEMPLATE.conf config/<TARGET>.conf

# 2) 환경 로드 후 전체 자동 실행
source /path/to/casp17-ligand/scripts/env_setup.sh
python3 run_pipeline.py config/<TARGET>.conf

# 유용한 옵션
python3 run_pipeline.py config/<T>.conf --dry-run                    # 명령만 확인
python3 run_pipeline.py config/<T>.conf --from pocket_validate --force  # 중간부터 재실행
python3 run_pipeline.py config/<T>.conf --only final_select --force     # 한 스텝만
```

- **idempotent**: 산출물이 이미 있으면 skip → 중간 실패 후 재실행하면 이어서 진행.
- **되돌아가기(피드백 루프)**: 결과가 불만족(예: 통과 포켓 0개)이면 기준을 고치고 `--from <스텝> --force`로 그 단계부터 다시.
- **config 주요 키**:
  - 필수: `target, results_dir, out_dir, ligand_tsv` (`scripts_dir`는 불필요 — 자기 위치 기준 자동 해석)
  - 켜기: `task`(P/PA), `python_pb`(PoseBusters), `refine`(정제), `experimental_cif`(fragment 비교),
    `templates`(후보별 템플릿 대조 = `10`; 파일 목록 콤마구분) 또는
    `templates_dir`(그 폴더의 `*.cif` 전부 자동 사용 — 퍼-타겟 폴더 `templates/{target}/` 권장),
    `split_by_conformation`
    - ⚠️ 템플릿은 **리간드 포함 RCSB 원본**이어야 overlap 계산됨(inputs의 `template_A_*.cif`는 리간드 제거본 → fold RMSD만)
  - 성능: `gnina_jobs`(gnina 병렬 워커 수, 기본 1=순차), `gnina_gpus`(라운드로빈 GPU id 예:`0,1,2`)
  - 강건성: `posebusters_topk`(클러스터당 검사할 상위 멤버 수, 기본 3 — 대표가 물리검증 실패해도
    같은 클러스터의 valid 멤버로 대표를 교체해 큰 합의 클러스터가 통째로 버려지는 것 방지)
  - 데이터 큐레이션: `exclude_models`(쉼표구분 모델 폴더명 제외, 예 `of3_collapse` — 나쁜 변형본이
    클러스터 중심을 오염시킬 때 collect 단계에서 배제. medoid 강건화로도 못 막는 대량 오염에 사용)
  - 자동화: `pocket_auto_relax`(**기본 off/opt-in** — 통과 0개는 "임계값이 빡빡했나 vs 진짜 포켓 없나"를
    사람이 판단해야 하는 신호이므로, 이미 완화가 맞다고 판단했을 때만 켬. 켜면 추가 계산 없이 임계값만 완화 재판정),
    `author`(미지정 시 env `CASP_AUTHOR` 사용, 그것도 없으면 CHANGE-ME)
  - 임계값: `pocket_threshold`, `gnina_cutoff`, `p2rank_dist`, `ligand_threshold`, `top_final`, …
  - 회의 반영(2026-07): `pocket_members_topk`(포켓별 채점 멤버 수→mean/std affinity·cnn, 기본 5, 0=대표만),
    `pocket_sweep`(cutoff 진단 콤마목록 예 `6,8,10,12` — 포켓 불명확 시 threshold 조정 참고),
    `template_pass_dist`(template_dist 이 값 이하면 포켓 pass 자동 승격, 기본 8.0).
    최종 composite = 합의 + **Boltz2(ligand_iptm)** + **Vina(gnina aff)** + **Gnina(CNN)**;
    포켓검증/리간드클러스터에 **평균·분산** 컬럼, 최종 정렬에 **template_dist tie-break** 병기.
    master_table에 Boltz2 native `confidence`(confidence_score) 파싱(정보용).

---

*이 문서는 스크립트 구조 지도입니다. 결과를 "믿을지 말지"의 해석 방법론은 [INTERPRETATION.md](INTERPRETATION.md)를 참고하세요.*
