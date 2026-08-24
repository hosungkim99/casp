# 최종 선정 근거 (Task=P)
종합점수 가중치: consensus 0.4, iptm(Boltz2) 0.25, aff(Vina) 0.15, cnn(Gnina) 0.2
(P=pose 중심, PA=affinity 비중↑ / size=클러스터 갯수, rmsd_mean·std=방향 응집도)

| model | pocket | lig_cluster | size | rmsd_mean±std | iptm | vina_aff | cnn | composite |
|---|---|---|---|---|---|---|---|---|
| model_1 | P1 | L1 | 132 | 1.247±0.296 | 0.96 | -6.01321 | 0.88605 | 0.9843 |
| model_2 | P1 | L2 | 63 | 1.0±0.285 | 0.921 | -5.78235 | 0.91092 | 0.7497 |
| model_3 | P1 | L3 | 59 | 1.055±0.276 | 0.947 | -6.57811 | 0.91334 | 0.6563 |
| model_4 | P1 | L44 | 2 | 0.271±0.0 | 0.921 | -4.91356 | 0.89826 | 0.5489 |
| model_5 | P1 | L18 | 7 | 1.173±0.337 | 0.87 | -5.91249 | 0.92197 | 0.5485 |

선정 원칙: 포켓 검증 통과 + 합의(클러스터 크기) + Boltz2(iptm) + Vina(affinity) + Gnina(CNN) 종합.
응집도 참고: size 클수록/rmsd_mean·std 작을수록 방향 합의가 강함(신뢰↑).
주의: 같은 포켓에서 여러 개가 뽑히면 방향 다양성 확인, 다른 포켓이 섞이면 hedge.
