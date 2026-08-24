# stage2 — Interpretation (봐야 할 파일 + 결론)

## 최종적으로 봐야 할 파일
| 목적 | 파일 |
|---|---|
| **제출물** | `L01LG<group>.tgz` (submit/ 를 tar) |
| 리간드별 선정 포즈 요약 | `final/stage2/poses.csv` |
| 개별 포즈(제출) | `final/stage2/submit/<cid>/<cid>LG.txt` |
| 왜 이 포즈인가(fragment별) | `<cid>/05_final/SELECTION_RATIONALE.md` |

## 어떤 포즈를 믿을까
- **composite 순** = 종합점수(합의 클러스터 크기 + iptm + affinity). model_1이 최상위.
- **pocket_pass=True** = 그 포켓이 p2rank+gnina 검증 통과 → 신뢰 높음. `False`면 저신뢰(폴백 선정).
- **여러 model이 같은 포켓·비슷한 자세** = 강한 합의. 서로 다른 포켓이면 hedge(다양성 확보).

## 포즈 개수의 의미
- 리간드당 1~5 포즈. **자세가 여러 종류**(dom_pose_frac 낮음)면 여러 model로 헤지.
- Stage1 분석: 포켓(위치)은 매우 수렴(0.94)하지만 **포즈(방향)는 덜 수렴(0.73)** → 그래서 최대 5개 제출로 방향 불확실성 대비.

## 저신뢰 fragment 판별
- `poses.csv`의 `pocket_pass=False` (또는 Stage1 `n_pockets_pass=0`) = cofold가 Site B로 수렴 못 한 fragment.
- 이런 건 포즈 신뢰 낮음 → 제출하되 상위 후보로 보지 말 것.

## 실제 정답과의 대조 (Stage2 열린 뒤)
- CASP가 **결합 확정 fragment**를 Stage2 입력으로 공개 → 우리 Stage1 예측이 맞았는지 직접 검증 가능.
- 그 목록으로 Stage1 9정의를 **실제 정답 AUC**로 재채점하면 대리 oracle보다 확실.

## 한 줄 요약
> 제출 = binder 각각의 ≤5 포즈(`submit/`) tar. composite·pocket_pass로 신뢰 판단. 위치는 확실(Site B), 방향은 여러 model로 헤지.
