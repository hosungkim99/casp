# pipelines — Interpretation (결과를 믿을지 판단하는 법)

각 출력을 "무엇을 보고 믿을지" 기준으로 정리. 컬럼 정의는 [FEATURES.md](FEATURES.md).

## ⭐ 가장 먼저 볼 것 — 05c_confidence/confidence_report.md
- **판정 줄**: HIGH=그대로 제출 / MEDIUM=수동 검토 / LOW=사람 개입.
- 근거 신호: 우세도 ≥30%·다모델 합의·af3 포함 → 사이트 확실 / p2rank ≤6Å·posebusters valid → 물리 타당 /
  fragment 공유 ≥3 → 실험 일치 / refine drift <1Å → 정제 안정.
- MEDIUM/LOW면 어느 신호가 약한지 보고 그 스텝으로 `--from` 복귀.

## 시각화 — 07_viz/
### visualization.png (4패널)
1. **리간드 PCA by 모델**: 색=모델. 여러 모델이 같은 덩어리에 겹치면 **합의**(좋음). 색이 흩어지면 모델별 다른 자리.
2. **리간드 PCA by 클러스터**: 가장 큰 덩어리 = 1순위 포켓.
3. **iptm 분포(boxplot)**: 모델별 스케일 다름 → **raw iptm 직접 비교 금지**의 근거. 그래서 선정은 consensus로.
4. **상위 클러스터 모델 구성(stacked bar)**: 한 막대에 **여러 색=cross-model 합의(신뢰↑)**, 단색=단일모델(의심).
- 판단: (1)(2) 큰 단일 덩어리 + (4) 다모델 = "강한 합의 포켓".

### contact_fingerprint.png
- 막대=잔기별 접촉빈도(%). **높고 고른 막대 다수 = 결합부위 한 곳 수렴**(좋음). **빨강=실험 fragment 공유 잔기**.
- 판단: ~100% 잔기 5개 이상 + 그중 빨강 존재 = 결합부위 확정적.

## 포켓 신뢰 — 02/03
- `size` 큼(합의↑) + `n_models` ≥2 (**af3 포함**) → 1순위 포켓 신뢰. 단일모델 단독이면 의심.
- `p2rank_dist` 작을수록 실제 캐비티. `conf_composition`은 noisy → 무시.

## 최종 선정 — 05(b)/selection_summary.csv
- `composite` 높을수록↑. model_1이 2~5위보다 확연히 높으면 1순위 압도적(좋음).
- 5개가 서로 다른 pocket_id면 = 사이트 hedge(불확실 → 분산 베팅).
- `posebusters_valid`는 반드시 True/NA여야 제출.

## 실험 검증 — 06_validation/fragment_compare.csv
- 진짜 fragment `n_shared` ≥3 = 실험적으로 같은 자리 = 강한 검증. 버퍼(EDO/TLA/UNX)는 0이 정상.

## 후보별 템플릿 대조 — 10_pose_vs_template (다중포켓 타겟)
- `overlap_pct`(우리 원자 중 실험 4Å 이내 비율), `jaccard`(접촉잔기 겹침), `convergence`(그 후보 수렴도).
- "L1이 아니라 L3가 실험과 더 맞는" 경우 진짜 상위 후보를 골라줌.

## 보강 검사
- `posebusters.csv`: 최종 5개 모두 valid면 통과. 무효 많으면 그 포켓/포즈 품질 의심.
- `refine_summary.csv`: `drift`<1 & `used`=refined = 정제 성공. drift 크면 원본 유지(안전).

## 한 줄 요약 — 언제 믿고 제출하나
> **confidence=HIGH** + **1순위 포켓 다모델(af3 포함) 합의** + **최종 5개 posebusters valid** + **fragment n_shared ≥3**
> — 네 가지 충족이면 결과를 신뢰하고 제출.
