# 선정 신뢰도 리포트

## 판정: **HIGH**
- HIGH: 자동선택 신뢰 가능 / MEDIUM: 검토 권장 / LOW: 사람 개입 필요

## 근거
- 1순위 포켓 P1: 우세도 100%(total 640), 합의 4모델(af3:160;bt2:160;of3:160;pt2:160), af3 포함
- 단일 포켓
- pose 수렴(top 리간드클러스터/포켓) = 99%
- p2rank 거리 = 2.82Å (pass=True), PoseBusters 최종 모두 valid=True
- 정제 최대 drift = 1.01Å (작을수록 안정)

## 최종 5개

| model | pocket | size | iptm | gnina | posebusters |
|---|---|---|---|---|---|
| model_1 | P1 | 631 | 0.94 | -5.98448 | True |
| model_2 | P1 | 3 | 0.779 | -6.16485 | True |
| model_3 | P1 | 6 | 0.871 | -1.95957 | True |

## 해석 가이드
- 합의·검증·유효성이 모두 일치 → 자동 선정 그대로 제출해도 무방.
