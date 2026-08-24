# 최종 선정 근거 (Task=P)
종합점수 가중치: consensus 0.5, iptm 0.3, affinity 0.2
(P=pose 중심, PA=affinity 비중↑)

| model | pocket | lig_cluster | size | iptm | gnina | composite |
|---|---|---|---|---|---|---|
| model_1 | P1 | L2 | 137 | 0.869 | -8.49272 | 0.754 |
| model_2 | P1 | L3 | 85 | 0.953 | -7.48088 | 0.73 |
| model_3 | P1 | L4 | 46 | 0.889 | -8.21484 | 0.4578 |
| model_4 | P1 | L9 | 13 | 0.936 | -6.64343 | 0.3909 |
| model_5 | P1 | L37 | 1 | 0.941 | -7.41229 | 0.389 |

선정 원칙: 포켓 검증 통과 + 합의(클러스터 크기) + interface 신뢰도 + 에너지 종합.
주의: 같은 포켓에서 여러 개가 뽑히면 방향 다양성 확인, 다른 포켓이 섞이면 hedge.
