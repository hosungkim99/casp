# 선정 신뢰도 리포트

## 판정: **HIGH**
- HIGH: 자동선택 신뢰 가능 / MEDIUM: 검토 권장 / LOW: 사람 개입 필요

## 근거
- 1순위 포켓 P1: 우세도 100%(total 640), 합의 4모델(af3:160;bt2:159;of3:160;pt2:160), af3 포함
- 포켓 분리 P1/P2 = 639.0x
- pose 수렴(top 리간드클러스터/포켓) = 21%
- p2rank 거리 = 2.01Å (pass=True), PoseBusters 최종 모두 valid=True
- 정제 최대 drift = 0.865Å (작을수록 안정)

## 최종 5개

| model | pocket | size | iptm | gnina | posebusters |
|---|---|---|---|---|---|
| model_1 | P1 | 137 | 0.869 | -10.20078 | True |
| model_2 | P1 | 85 | 0.953 | -8.42332 | True |
| model_3 | P1 | 46 | 0.889 | -9.23298 | True |
| model_4 | P1 | 13 | 0.936 | -8.55559 | True |
| model_5 | P1 | 1 | 0.941 | -9.10882 | True |

## 해석 가이드
- 합의·검증·유효성이 모두 일치 → 자동 선정 그대로 제출해도 무방.
