# 선정 신뢰도 리포트

## 판정: **MEDIUM**
- HIGH: 자동선택 신뢰 가능 / MEDIUM: 검토 권장 / LOW: 사람 개입 필요

## 근거
- 1순위 포켓 P1: 우세도 100%(total 640), 합의 4모델(af3:160;bt2:160;of3:160;pt2:158), af3 포함
- 포켓 분리 P1/P2 = 319.0x
- pose 수렴(top 리간드클러스터/포켓) = 21%
- p2rank 거리 = 7.14Å (pass=True), PoseBusters 최종 모두 valid=True
- 정제 최대 drift = 0.808Å (작을수록 안정)

## 최종 5개

| model | pocket | size | iptm | gnina | posebusters |
|---|---|---|---|---|---|
| model_1 | P1 | 132 | 0.96 | -6.57153 | True |
| model_2 | P1 | 63 | 0.921 | -6.29016 | True |
| model_3 | P1 | 59 | 0.947 | -7.36723 | True |
| model_4 | P1 | 2 | 0.921 | -5.84813 | True |
| model_5 | P1 | 7 | 0.87 | -6.57615 | True |

## 해석 가이드
- 부분적으로만 일치 → 1~2순위 포켓 수동 확인, hedge 비중 조정 고려.
