# stage1 — Features (생성 파일의 컬럼 정의)

## 원 스코어 (모든 정의의 재료)
| 기호 | 뜻 | 범위 | 출처 |
|---|---|---|---|
| boltz | boltz-2 `affinity_probability_binary` (6 seed 평균) | 0~1 (보정확률) | affinity json |
| cnn | gnina CNNscore (pose 신뢰도) | 0~1 | gnina |
| cnnaff | gnina CNNaffinity (예측 결합강도, pK류) | ~2~10 | gnina |
| aff | gnina/Vina minimizedAffinity (에너지) | 음수(낮을수록 강함) | gnina |
| hac | heavy atom 수 (fragment 크기) | 정수 | SMILES |

`rank(x)` = 값 x를 라이브러리 전체에서 순위정규화 → 0~1.

## binding_row.csv  (per-fragment 원자료, 2_binding 출력)
`cid, boltz, cnn, cnnaff, aff, hac, n_pockets, n_pockets_pass, pocket_sizes, dom_pocket_size, dom_pocket_frac, n_pose_clusters, pose_sizes, dom_pose_size, dom_pose_frac, n_poses, pocket_residues, rep_cif, note`
- `n_pockets` 포켓 후보 수 / `n_pockets_pass` p2rank+gnina 통과 수
- `pocket_sizes` 포켓 후보 크기들(예 `24;4;2`) / `dom_pocket_frac` 지배 포켓에 모인 비율
- `pose_sizes` 지배 포켓 안 자세 클러스터 크기들 / `dom_pose_frac` 지배 자세 비율
- `pocket_residues` 대표 pose 5Å 잔기 / `rep_cif` 대표 pose 경로 / `note` no_pose 등

## binding_scores.csv  (9정의, 3_aggregate 출력)
`cid` + 9정의 + 원값 `cnn, cnnaff, aff` + `note`

| 정의 | 수식 | 재는 것 |
|---|---|---|
| prob_boltz | boltz (그대로) | 유일한 **절대확률** |
| prob_cnn | rank(cnn) | pose 품질 |
| prob_cnnaff | rank(cnnaff) | 결합 강도 |
| prob_vina | rank(−aff) | 물리 에너지 |
| prob_gnina | 0.5·rank(cnn)+0.5·rank(cnnaff) | gnina 합의 |
| prob_cons3 | (rank(cnn)+rank(cnnaff)+rank(−aff))/3 | 도킹 3신호 합의 |
| prob_LE_caf | rank(cnnaff/hac) | 원자당 강도(크기보정) |
| prob_LE_vina | rank(−aff/hac) | 원자당 물리효율 |
| prob_combined | 0.5·rank(boltz)+0.5·**rank(gnina평균)** | boltz+gnina 독립 합의. gnina평균을 재-rank해 진짜 50:50 blend(타겟무관) |

> 주의: 값 0.9는 "90% 결합"이 아니라 **라이브러리 내 순위**. 절대확률 해석 가능한 건 prob_boltz뿐.

## pocket_clusters.csv  (포켓/클러스터, 3_aggregate 출력)
`cid, n_poses, n_pockets, n_pockets_pass, pocket_sizes, dom_pocket_frac, n_pose_clusters, pose_sizes, dom_pose_frac, pocket_residues, note`
- 포켓 수렴도(`dom_pocket_frac`, 평균 0.94)와 포즈 수렴도(`dom_pose_frac`, 평균 0.73) 구분.
- `pocket_residues` 로 콘센서스 포켓(6_consensus_pocket) 계산.

## bind.txt  (제출, 4_finalize 출력)
탭 구분, fragment마다 한 줄:
```
L010001 <확률>  A190,A192,...   (binder: 확률 + 5Å 포켓 잔기 ChainResnum)
L010002 0                       (non-binder: 0)
```
