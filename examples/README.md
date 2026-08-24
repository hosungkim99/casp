# examples — 단일 리간드 end-to-end 예제

[`../pipeline/`](../pipeline/) 본체를 **1 단백질 : 1 리간드** CASP 타겟에 그대로 적용한 예제 모음입니다.
[`../L01/`](../L01/)이 multi-copy + Zn의 복잡한 사례라면, 여기 3개는 파이프라인이 **표준 단일 리간드 문제**에서 어떻게 도는지를 입력→제출까지 한눈에 보여줍니다.

## 예제 3개

| 타겟 | 리간드 | 특징 | 신뢰도 판정 | p2rank 거리 |
|---|---|---|:--:|:--:|
| [**T2411**](T2411/) | U08 (아세트아닐리드) | 작은 방향족 아마이드 | 🟢 **HIGH** | 2.82 Å |
| [**T2412**](T2412/) | DEL (트리아진–피페라진) | 유연한 CF₃ 치환체 | 🟡 **MEDIUM** | 7.14 Å |
| [**T2413**](T2413/) | LIG (글리코사이드) | 다중 OH 대형 분자 | 🟢 **HIGH** | 2.01 Å |

세 예제가 **신뢰도 시스템의 판정 로직**을 잘 보여줍니다: T2412와 T2413은 pose 수렴도가 21%로 같지만, T2412는 포켓 중심이 druggable cavity에서 7.14Å 떨어져(p2rank) **MEDIUM(검토 권장)**, T2413은 2.01Å로 **HIGH(자동선정 신뢰)**. 즉 판정이 하나의 점수가 아니라 합의·검증·수렴을 종합합니다.

## 파이프라인 흐름 (모든 예제 공통)

```
ligand.tsv(SMILES) ──▶ 4모델 co-folding 예측(af3·bt2·of3·pt2, 다중 seed)
   ──▶ 포켓 후보(합의 클러스터) ──▶ p2rank·gnina 검증 ──▶ 리간드 방향 클러스터
   ──▶ 종합점수 최종 5개 선정 ──▶ 신뢰도 판정(HIGH/MEDIUM/LOW) ──▶ CASP LG 제출
```

## 각 폴더에 담긴 것

| 파일 | 내용 |
|---|---|
| `ligand.tsv` | 입력 (ID·이름·SMILES·Task) |
| `confidence_report.md` | **신뢰도 판정** + 근거 (핵심 산출물) |
| `SELECTION_RATIONALE.md` · `selection_summary.csv` | 최종 선정 근거·점수 |
| `confidence.csv` | 판정 수치 원본 |
| `figures/visualization.png` | 리간드 pose 선정 개요 (PCA·iptm·모델구성) |
| `figures/contact_fingerprint.png` | 최종 pose의 단백질 접촉 지문 |
| `model_1.cif` | **최종 1순위 pose 구조** |
| `submission/*LG_model1.txt` | CASP LG 제출 포맷 (1순위) |

> 중간 단계(00–04)·나머지 후보 구조·대용량 로그는 제외했습니다. 서버 경로는 `/path/to/` 로 일반화.
