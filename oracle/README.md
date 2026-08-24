# oracle — 결합확률 정의 검증 세트

Stage1의 결합확률 정의(9가지 후보) 중 **어느 것이 실제 결합 분자를 가장 잘 골라내는지**를 취향이 아니라 **AUC**로 판정하기 위한 정답 세트입니다. 채점·실행 스크립트는 [`../pipeline/oracle/`](../pipeline/oracle/) 에 있습니다.

## 구성
- `oracle_ligands.tsv` — 검증 분자 15개. `label=1`은 실험적으로 결합이 확인된 **양성 3개**(Flumequine·Hesperetin·Foliosidine, BFT-3 논문 PMC9514063), `label=0`은 일반 약물 **디코이 12개**.
- `ORACLE_검증_정리.md` — 검증 결과 정리.

## 아이디어
각 결합확률 정의로 15개 분자를 점수화하고, **양성이 디코이보다 높은 점수를 받는 비율(AUC)** 로 정의를 비교합니다(양성3 × 디코이12 = 36쌍). AUC가 가장 높은 정의를 제출용으로 채택합니다.

> smi/·affinity_yaml/·af3 input 등 파생 입력은 `oracle_ligands.tsv` 로부터 `pipeline/oracle/make_oracle_yaml.py` 가 생성하므로 저장소에서는 제외했습니다.
