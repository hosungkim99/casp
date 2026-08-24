# CASP17 Protein–Ligand Structure Prediction Pipeline

> 단백질 서열과 리간드(SMILES)가 주어졌을 때, 팀의 co-folding 모델(AlphaFold3·Boltz-2·Protenix 등)이 만든 **수천 개의 예측 구조를 합의(consensus) 기반으로 자동 선별해 CASP 제출용 5개 포즈로 압축**하는 파이프라인입니다. CASP17 리간드 부문(Ligand category)에 사용했습니다.

![Python](https://img.shields.io/badge/Python-3.13-blue)
![Runtime](https://img.shields.io/badge/Runtime-Linux%20GPU%20%2F%20SLURM-black)
![Domain](https://img.shields.io/badge/Domain-Structural%20Biology-1baf7a)

---

## 한눈에 — 이 파이프라인이 푸는 문제

CASP 대회는 `단백질 서열 + 리간드(SMILES) + Task(P=포즈 예측 / A=친화도)`를 주고 "이 리간드가 단백질의 **어디에, 어떻게** 붙는지"를 3일 안에 맞히라고 요구합니다.

팀은 먼저 co-folding 모델을 여러 개·여러 seed로 돌려 수천 개의 후보 구조를 만듭니다. **이 저장소의 파이프라인은 그 수천 개를 걸러 제출용 5개로 자동 압축합니다.**

핵심 아이디어는 **합의(consensus)** 입니다. 하나의 raw 점수로 줄 세우지 않고, *여러 모델이 리간드를 같은 자리에 놓으면 그 자리가 정답일 확률이 높다* 는 원리로, **"많이 모인 자리 = 좋은 자리"** 를 먼저 정한 뒤 그 안에서 방향·물리 타당성을 검증합니다.

```
수천 개 예측 ──▶ 포켓 후보 ──▶ 검증 통과 포켓 ──▶ 포켓별 리간드 방향 ──▶ 최종 5개 ──▶ 제출 파일
 (수집·랭킹)      (합의 클러스터)   (p2rank·gnina)      (SC-RMSD 클러스터)     (종합점수)    (CASP LG)
```

계층형 순서: **① 자리(포켓)를 먼저 정한다 → ② 그 자리 안에서 리간드 방향을 정한다 → ③ 물리적으로 말이 되는지 검증한다 → ④ 최종 5개를 선정한다.**

---

## 설계에서 특히 신경 쓴 것

- **합의(consensus) 우선** — 모델마다 점수 스케일이 달라 raw 점수 직접 비교는 신뢰할 수 없습니다. 그래서 "여러 모델이 모인 자리"를 신뢰 신호로 삼습니다.
- **계층형 분리** — 포켓(자리) → 리간드 방향(포즈)을 단계로 분리해, 같은 자리에 뒤집혀 놓인 포즈까지 구분합니다.
- **물리 검증 이중화** — p2rank(진짜 캐비티인가) + gnina(에너지가 타당한가) + PoseBusters(충돌·기하 유효성)로 가짜 자리를 걸러냅니다.
- **근거 기반 선택** — 결합확률 정의 후보들을 **oracle 세트로 AUC 채점**해, 취향이 아니라 지표로 제출 정의를 정합니다(→ [`oracle/`](oracle/)).
- **재현성·복구** — 파이프라인은 idempotent(산출물 있으면 skip)라 중간 실패 후 이어서 실행되고, `--from <step> --force`로 특정 단계부터 재실행하는 피드백 루프를 지원합니다.

---

## 저장소 구성

```
.
├─ pipeline/     범용 선별·검증 파이프라인 (모든 타겟 공통)
│   ├─ common/       공용 라이브러리 (파서·기하·gnina 점수)
│   ├─ pipelines/    메인 파이프라인 (지휘자 + core/boost/analysis/ml)
│   ├─ stage1/       결합확률 + 결합 포켓 잔기 산출 → 제출(bind.txt)
│   ├─ stage2/       확정 fragment의 복합체 포즈(≤5) 예측 → 제출(tgz)
│   └─ oracle/       결합확률 정의를 AUC로 검증하는 스크립트
├─ oracle/       oracle 검증 세트 (양성 3 + 디코이 12 분자표) — 위 스크립트의 정답 데이터
└─ L01/
    └─ scripts/  실제 타겟 L01("1 단백질 : 다수 리간드")에 적용한 스크립트 (단계별)
```

`pipeline/` 이 재사용 가능한 본체이고, `L01/scripts/` 는 그 본체를 실제 대회 타겟에 맞게 확장·운용한 예시입니다.

처음 본다면 **[`pipeline/pipelines/OVERVIEW.md`](pipeline/pipelines/OVERVIEW.md)** 부터 읽으면 전체 지도를 볼 수 있습니다.

| 문서 | 내용 |
|---|---|
| [`pipeline/pipelines/OVERVIEW.md`](pipeline/pipelines/OVERVIEW.md) | 전체 지도 — 스크립트별 역할·처리 순서 |
| [`pipeline/pipelines/USAGE.md`](pipeline/pipelines/USAGE.md) | 설치·실행·config 키 |
| [`pipeline/pipelines/FEATURES.md`](pipeline/pipelines/FEATURES.md) | 출력 폴더·컬럼 정의 |
| [`pipeline/pipelines/INTERPRETATION.md`](pipeline/pipelines/INTERPRETATION.md) | 결과를 믿을지 판단하는 기준 |
| [`pipeline/stage1/OVERVIEW.md`](pipeline/stage1/OVERVIEW.md) · [`pipeline/stage2/OVERVIEW.md`](pipeline/stage2/OVERVIEW.md) | Stage1/2 산출물 |
| [`oracle/`](oracle/) · [`L01/scripts/`](L01/scripts/) | oracle 검증 세트 · 실제 타겟 적용 예시 |

---

## 기술 스택

- **언어/환경**: Python 3.13, Bash, Linux GPU 서버, **SLURM** 배치(sbatch)
- **구조/기하**: gemmi(구조 파싱), NumPy, RDKit(SMILES·대칭보정 RMSD), Kabsch 정렬
- **도킹/스코어링**: gnina(CNN score/affinity + Vina), p2rank(포켓 예측), PoseBusters(물리 유효성)
- **구조 예측(입력 생성)**: AlphaFold3, Boltz-2, OpenFold3, Protenix
- **ML 실험**: XGBoost, scikit-learn
- **시각화**: matplotlib (리간드 PCA·접촉 지문)

---

## 실행 방법

경로 4개와 task만 채운 config 하나로 전 단계가 순서대로 자동 실행됩니다.

```bash
# 1) config 작성 (필수 키: target, results_dir, out_dir, ligand_tsv, task)
cp pipeline/pipelines/config/TEMPLATE.conf pipeline/pipelines/config/<TARGET>.conf

# 2) 환경 로드 후 전체 자동 실행
python3 pipeline/pipelines/run_pipeline.py pipeline/pipelines/config/<TARGET>.conf

# 유용한 옵션
python3 pipeline/pipelines/run_pipeline.py <conf> --dry-run                     # 명령만 확인
python3 pipeline/pipelines/run_pipeline.py <conf> --from pocket_validate --force  # 중간부터 재실행
```

가장 먼저 볼 결과는 `05c_confidence/confidence_report.md`의 **HIGH / MEDIUM / LOW** 판정입니다. HIGH면 그대로 제출, MEDIUM/LOW면 사람이 검토합니다.

> 경로 표기: 스크립트·config의 `/path/to/...` 와 `USERNAME` 은 자리표시자입니다. [`env_setup.sh.example`](env_setup.sh.example)를 복사해 자신의 클러스터 경로/계정에 맞게 채운 뒤 `source` 하세요. gnina·p2rank·PoseBusters·co-folding 모델은 pip이 아닌 외부 바이너리/별도 환경이며, 파이썬 의존성은 [`requirements.txt`](requirements.txt) 참고.

---

## 역할 / 기여

CASP17 리간드 팀에서 **예측 구조 선별·검증 파이프라인 전반**을 담당했습니다. 팀이 co-folding 모델로 생성한 수천 개의 후보 구조를 입력으로 받아, 합의 클러스터링 → 포켓 검증(p2rank·gnina) → 포즈 클러스터링 → 종합점수 기반 최종 선정 → CASP 제출 포맷 변환에 이르는 자동화 파이프라인과 그 문서를 설계·구현했으며, 실제 타겟(L01 등)에 맞춰 stage1/stage2 확장 스크립트를 작성하고 운용했습니다. 또한 결합확률 정의를 oracle 세트로 AUC 검증하는 실험을 설계했습니다.

<!-- ▲ 팀 내 정확한 담당 범위(전체 설계 주도 / 특정 단계 담당 등)를 본인 상황에 맞게 한 문장으로 다듬어 주세요. -->

---

## 참고 / 한계

- 이 저장소에는 **코드와 소량의 정의 데이터만** 포함되어 있습니다. 입력 예측 구조·중간 산출물·제출 파일 등 대용량 데이터는 제외했습니다.
- 스크립트·config의 절대경로는 특정 클러스터 환경을 지우고 `/path/to/...`·`USERNAME` 자리표시자로 일반화했습니다.
- CASP 대회 규정상 공개 가능한 범위 내에서 정리한 코드입니다.

---

📧 hosungkim99@gmail.com · 🔗 github.com/hosungkim99
