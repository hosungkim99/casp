# stage1 — Interpretation (봐야 할 파일 + 결론)

## 제출을 위해 최종적으로 봐야 할 파일
| 목적 | 파일 |
|---|---|
| **제출물** | `final/stage1/L01LG<group>.bind.txt` (확률+포켓) |
| **정의 선택 근거** | `final/stage1/definition_choice.txt` |
| 9정의 전체 값 | `stage1/binding_scores.csv` |
| 포켓/수렴도 | `stage1/pocket_clusters.csv` |
| 9정의 비교 | 5_defs_viz 출력 png |
| 콘센서스 포켓 | 6_consensus_pocket 출력 txt |

## 어떤 결합확률 정의를 쓸까 (핵심 결론)
oracle 검증(`oracle/` 참고, 양성3 억제제 vs 음성12 디코이 AUC)에서:
- **상위(강함)**: prob_boltz(0.833), prob_vina(0.806), prob_combined(0.806) — 36쌍 중 1쌍 차이라 **통계적 동점**.
- **약함**: prob_cnn/cnnaff/gnina (0.58~0.61) → **GNINA CNN 점수는 binder/decoy 구분 약함**.
- **버릴 것**: prob_LE_caf (0.361) → **랜덤 이하, 사용 금지**.

**결론**:
- {GNINA, VINA} 중 **VINA가 신호를 가짐**, GNINA CNN은 약함 → "그냥 평균(cons3)"보다 **VINA·boltz 가중**이 유리.
- 서버 제출은 **prob_combined**(상위 동점 + 헤지). boltz 신뢰 낮으면 prob_cons3.
- 단 양성 3개뿐이라 대리(proxy) 검증 — Stage2에서 실제 binder 공개되면 재검증.

## 포켓은 어떻게 보나
- `pocket_clusters.csv`의 `pocket_residues` 콘센서스 = **알로스테릭 Site B (A190/A192/A232/A257/A259/A260/A268 ≈ 98% fragment 공유)**.
- BFT-3 논문(PMC9514063) 알로스테릭 억제제 결합부위와 일치 → 촉매 아연부위 아님.
- `n_pockets_pass=0` fragment(~23%)는 검증 포켓 없는 저신뢰(약한/비결합 후보).

## 수렴도 해석
- **포켓(위치)**: `dom_pocket_frac` 평균 0.94 — 90.7%가 지배 포켓에 pose 80%+ 집중 → 위치는 매우 확실.
- **포즈(방향)**: `dom_pose_frac` 평균 0.73 — 위치보다 덜 수렴(보통 2~3자세). → Stage2에서 포즈 여러 개 제출하는 이유.

## 한 줄 요약
> 제출 = `bind.txt`(prob_combined + Site B 포켓). oracle이 "boltz/vina/combined 강함, GNINA CNN 약함, LE_caf 버려라"를 알려줬고, 포켓은 Site B로 3중 교차검증됨.
