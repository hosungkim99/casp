# oracle — Interpretation (봐야 할 것 + 결론)

## 최종적으로 봐야 할 것
| 목적 | 위치 |
|---|---|
| **9정의 AUC 순위** | `s1_oracle_auc.py` 콘솔 출력 |
| 분자별 원자료 | `oracle/outputs/<name>/05_stage1_binding/binding_row.csv` |
| 결과 정리 | `oracle/ORACLE_검증_정리.md` |

## 결론 (L01 oracle 결과)
- **상위(강함)**: prob_boltz(0.833), prob_vina(0.806), prob_combined(0.806) — **36쌍 중 1쌍 차이 = 통계적 동점**.
- **약함**: prob_cnn/cnnaff/gnina (0.58~0.61) → **GNINA CNN 점수는 구분력 약함**.
- **버려야 함**: prob_LE_caf (0.361) → **랜덤 이하, 사용 금지**.

**핵심 결론**:
1. {GNINA, VINA} 중 **VINA가 신호를 가짐**, GNINA CNN은 약함 → 단순 평균(cons3 0.736)보다 **VINA·boltz 가중**이 나음.
2. **boltz affinity가 오히려 1등** — "boltz affinity 신빙성 낮다"는 우려와 반대 방향(단, 표본 3개라 강하게 주장은 못 함).
3. 서버 제출은 **prob_combined**(상위 동점 + 독립 두 방법 헤지). boltz 신뢰 낮으면 **prob_cons3**.

## 한계 (반드시 같이 볼 것)
- 양성 3개뿐 → AUC가 거침(1/36 단위). 참고치.
- 디코이는 "결합 안 할 것으로 가정"이지 실험 확인 음성 아님.
- 양성은 약물, 라이브러리는 fragment → 화학공간 다름.
- → **대리(proxy) 검증**. 진짜 검증은 Stage2에서 실제 binder 공개 후 재채점.

## AUC를 어떻게 해석하나
- AUC 0.5 근처(cnn/cnnaff/gnina) = 그 정의로는 양성/음성 구분 거의 안 됨 → 신뢰 낮음.
- AUC < 0.5(LE_caf) = 오히려 디코이를 높게 줌 → 크기보정이 이 세트에선 역효과.

## 한 줄 요약
> oracle = "9정의 중 뭐가 맞나"를 AUC로 채점. 결론: **boltz/vina/combined 강함, GNINA CNN 약함, LE_caf 버려라.** 단 표본 작아 proxy — Stage2 실제 정답으로 재검증 예정.
