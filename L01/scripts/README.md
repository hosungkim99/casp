# L01/scripts — 실제 타겟 L01 적용 스크립트

CASP17 타겟 **L01("1 단백질 : 다수 리간드")** 에 [`../../pipeline/`](../../pipeline/) 본체를 적용·확장해 운용한 스크립트입니다. 여러 co-folding 모델로 binder를 예측하고, 합의 포켓 검증을 거쳐 복합체 포즈를 선정·제출하는 Stage2 흐름을 단계별 폴더로 정리했습니다.

## 단계 구성
| 폴더 | 역할 |
|---|---|
| `01_inputs/` | 모델별(AF3·Boltz-2·Protenix) co-folding 입력 생성 |
| `02_inference/` | 각 모델 추론 실행 (SLURM 배치) |
| `03_collect/` | 예측 구조 수집·합의(consensus) 취합 |
| `04_pipeline/` | 포켓 후보 실행·p2rank 검증 |
| `05_select/` | 다중 copy 리간드 포즈 최종 선정 |
| `06_submission/` | CASP LG 제출 포맷 생성 |
| `07_validation/` | 물리 유효성·정제·거리 검증 (exosite/Zn 등) |
| `_archive/` | 초기·보조 분석 스크립트 (참고용 보관) |

> 경로·계정은 특정 환경을 지우고 `/path/to/...`·`USERNAME` 으로 일반화했습니다. 실제 실행은 본인 환경에 맞게 경로를 채워야 합니다.
