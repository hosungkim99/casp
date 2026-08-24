# 최종 선정 근거 (Task=P)
종합점수 가중치: consensus 0.5, iptm 0.3, affinity 0.2
(P=pose 중심, PA=affinity 비중↑)

| model | pocket | lig_cluster | size | iptm | gnina | composite |
|---|---|---|---|---|---|---|
| model_1 | P1 | L1 | 631 | 0.94 | -5.50445 | 0.9875 |
| model_2 | P1 | L3 | 3 | 0.779 | -5.92895 | 0.2 |
| model_3 | P1 | L2 | 6 | 0.871 | 0.86185 | 0.1738 |

선정 원칙: 포켓 검증 통과 + 합의(클러스터 크기) + interface 신뢰도 + 에너지 종합.
주의: 같은 포켓에서 여러 개가 뽑히면 방향 다양성 확인, 다른 포켓이 섞이면 hedge.
