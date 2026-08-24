# 파이프라인 분석 & cutoff 정리 (0722 파이프라인 정리)

> 회의 내용에 맞춰 스크립트를 개조하기 위한 **현행 파이프라인 분석 문서**.
> 대상: `pipeline/pipelines/` (core 0~6,9 + boost 4b·4c·5b·5c + analysis 7·8·10·11).
> 지휘자: `run_pipeline.py` (config 1개 읽어 순서대로 자동 실행).

---

## 1. 실행 순서 & 플로우차트

`run_pipeline.py`가 config를 읽어 아래 순서로 스텝을 자동 호출한다.
(`⚙️`=config 플래그로 삽입되는 보강/검증 스텝, `★`=회의 내용과 직접 관련된 핵심 스텝)

```
[00] collect            0_rank_poses.py       수천 예측 수집 → master_table (interface iptm 랭킹)
  │
[00b] cache             0b_cache_geometry.py  전 cif 1회 파싱 → geom_cache.pkl
  │
[01] protein_cluster    1_protein_cluster.py  Cα PCA → 형태(conf) 라벨 + reference 구조 선정   (Step 0)
  │
[02] pocket_candidates  2_pocket_candidates.py 리간드 센트로이드 그리디 클러스터 → 포켓 후보 top10  (Step 1) ★
  │
[03] pocket_validate    3_pocket_validate.py  포켓 대표 1개에 p2rank + gnina → pass/fail 압축      (Step 2) ★
  │
 (⚙️ python_pb):  4b_prep_validity → 4c_posebusters   (PoseBusters 물리 유효성)
  │
[04] ligand_cluster     4_ligand_cluster.py   통과 포켓 내 SC-RMSD 재클러스터 → 방향 대표          (Step 3)
  │
[05] final_select       5_final_select.py     composite(합의+iptm+gnina) 상위 5개 선정             (Step 4) ★
  │
 (⚙️ refine):     5b_refine.py                gnina --minimize 국소정제 (drift 안전장치)
  │
[06] contact_residues   6_contact_residues.py 1순위 포켓 접촉잔기 빈도(4Å)
  │
 (⚙️ experimental_cif): 7_fragment_compare.py 실험 fragment 접촉잔기 비교 (단일 앵커)
 (⚙️ templates):        10_pose_vs_template.py 후보별 템플릿 overlap/jaccard 대조                  ★(template)
  │
[07] visualize          8_visualize.py        4패널 + 접촉지문 PNG
  │
[08] casp_lg            9_make_casp_lg.py      최종 5개 → CASP LG 제출 포맷
  │
 (항상 끝):       5c_confidence.py             HIGH/MEDIUM/LOW 신뢰도 판정
                  11_make_summary.py           SUMMARY_DRAFT.md 초안
```

**핵심 아이디어**: raw 점수 1개로 줄세우지 않고 **"많이 모인 자리 = 좋은 자리"(consensus)**.
계층형: ① 포켓 먼저 → ② 포켓 내 방향 → ③ 물리 검증 → ④ 최종 5개.

---

## 2. Cutoff 정리 (실행 순서대로)

> 값 표기: `키 = 기본값` (config 키 있으면 병기). **★** = 회의에서 손볼 지점.

### [00] collect — `0_rank_poses.py`
- 수치 cutoff 없음. 랭킹 키 = `(has_clash 없음, ligand_iptm, plddt, -gpde)` 내림차순.
- `exclude_models` (opt): 오염 모델 폴더 제외.

### [01] protein_cluster — `1_protein_cluster.py` (Step 0)
- `gap_frac = 0.20` (config `conf_gap_frac`): **PC1에서 (전체범위 × 0.20)보다 큰 틈**에서만 형태 분리. (단일형태 과분할 방지)
- `min_conf_frac = 0.02`: 형태 그룹 최소 크기 = `max(5, 0.02·N)`. 이보다 작아지는 분리는 무시.
- `threshold = 15.0` → **폐기됨(무시)**. 옛 절대임계 잔재.
- reference = **PC공간 중심(medoid)** 구조 (rank1 아님).

### [02] pocket_candidates — `2_pocket_candidates.py` (Step 1) ★
- **`pocket_threshold = 8.0` Å**: 리간드 센트로이드 **그리디 클러스터 거리 임계**. ← 회의 "클러스터링 cutoff distance 조정" 대상.
- `topn_pockets = 10`: 상위 포켓 후보 개수.
- `split_by_conformation` (opt): 형태별 분리 클러스터.

### [03] pocket_validate — `3_pocket_validate.py` (Step 2) ★
- **`gnina_cutoff = -4.0`**: pass 조건 ① gnina affinity(=vina) `<= -4.0`.
- **`p2rank_dist = 6.0` Å**: pass 조건 ② 최근접 p2rank 포켓 거리 `<= 6.0`.
- **pass = (aff ≤ -4.0) AND (p2rank_dist ≤ 6.0)**. → 포켓 **대표 cif 1개**에만 적용.
- `pocket_auto_relax` (opt, 기본 off): 통과 0개면 `cutoff += 1.0`, `dist += 2.0`, 최대 `relax_steps=3`회 완화 재판정(재계산 없음).

### (⚙️) 4b/4c PoseBusters — `python_pb` 설정 시
- `posebusters_topk = 3`: 클러스터당 검사할 상위(중심성) 멤버 수.

### [04] ligand_cluster — `4_ligand_cluster.py` (Step 3)
- **`ligand_threshold = 2.0` Å**: 포켓 내 리간드 **SC-RMSD(대칭보정) 클러스터 임계**.

### [05] final_select — `5_final_select.py` (Step 4) ★
- `top_final = 5`: 최종 선정 개수.
- **가중치(composite)**:
  - Task **P** : `consensus 0.5, iptm 0.3, aff 0.2`
  - Task **PA**: `consensus 0.35, iptm 0.25, aff 0.4`
- `composite = w_c·norm(size) + w_i·norm(iptm) + w_a·norm(-gnina_aff)` (각 항 min-max 정규화).
- gnina 점수는 포켓검증(3)에서 같은 cif면 **재사용**(재계산 생략).

### [06] contact_residues — `6_contact_residues.py`
- `threshold = 2.0` Å: 포켓 멤버 판정 리간드 concat-RMSD.
- `cutoff = 4.0` Å: 리간드-잔기 **접촉 거리**.
- 수렴도 `< 0.4` → 대표성 낮음 경고.

### (⚙️) 보강/검증
- `5b_refine`: `refine_max_drift = 2.0` Å(초과 시 원본 유지), `refine_iters = 25`(config)/10(기본).
- `5c_confidence`: `p2rank_ok = 6.0`, `dominance_high = 0.30`.
  - **HIGH** = n_models≥2 AND af3 참여 AND dominance≥0.30 AND p2d≤6.0 AND posebusters 전부 valid
  - **LOW** = n_models<2 OR dominance<0.15 OR p2d>8.0 OR posebusters 무효
  - 그 외 = MEDIUM.
- `10_pose_vs_template`: overlap_pct(우리원자 중 실험리간드 4Å 이내 비율), jaccard(접촉잔기 겹침).

---

## 3. 회의 내용 vs 현행 — 갭 분석

> 회의 정리:
> **1. 포켓 찾기**: cofolding+클러스터링 / template 찾기.
>   1-1. 포켓 불명확 시(큰 클러스터 없음 / 리간드 template 없음):
>     - 클러스터링 cutoff distance 조정.
>     - **각 포켓 리간드들에 Vina·Gnina 실행 → score 평균이 낮은 곳 = 포켓**.
> **2. 포즈 정하기**: **Boltz2·Vina·Gnina 세 값 조합**해 최고 스코어 선정 (클러스터 리간드 갯수·평균·분산 분석).

| 회의 항목 | 현행 스크립트 | 상태 | 갭 |
|---|---|---|---|
| cofolding + 클러스터링 | `2_pocket_candidates` | ✅ | 없음 |
| template 찾기 | `10_pose_vs_template` | ⚠️ | 존재하나 **검증용**(최종 후보 대조). 포켓 *찾기/결정*에는 미반영 |
| 클러스터링 cutoff distance 조정 | `--threshold`(8.0) | ⚠️ | 수동 1값만. 회의 "조정해본다"=**여러 값 스윕** 필요 |
| 포켓 리간드에 **Vina** 실행 | `gnina_affinity`(=vina항) | ⚠️ | gnina가 vina affinity를 이미 출력. **별도 AutoDock Vina 재도킹은 없음** |
| 포켓 리간드에 **Gnina** 실행 | `cnn_score/cnn_affinity` | ✅ | gnina CNN 점수 이미 계산 |
| **score 평균 낮은 곳 = 포켓** | `3_pocket_validate` | ❌ | 현재 **포켓 대표 1개**만 채점 + **pass/fail 이분법**. 포켓 멤버 **평균**으로 순위/결정하지 않음 |
| 포즈: **Boltz2** 값 | `iptm`(rep_iptm) | ⚠️ | Boltz2 신뢰도를 iptm으로만 대용. 명시적 Boltz2 confidence 항 아님 |
| 포즈: **Vina** 값 | `gnina_affinity` | ✅ | 이미 composite의 `aff` 항 |
| 포즈: **Gnina** 값(CNN) | `cnn_score` (미사용) | ⚠️ | 계산은 하나 composite에 **미반영**(현재 consensus+iptm+aff만) |
| 세 값 **조합** | `composite` 가중합 | ⚠️ | 조합식은 있으나 CNN 미포함, Boltz2 대용 |
| 클러스터 갯수·**평균·분산** 분석 | `size`만 | ❌ | size(갯수)만 사용. **평균·분산 미산출** |

### 핵심 개조 포인트 (요약)
1. **`3_pocket_validate`**: 포켓 대표 1개 → **포켓 멤버 다수 채점 + 평균(±분산)** 기반 포켓 결정으로 확장.
2. **cutoff 스윕**: `2_pocket_candidates`에 여러 threshold 시도/자동 선택 옵션.
3. **`5_final_select`**: composite에 **gnina CNN(gnina 점수)** 항 추가, **Boltz2 confidence** 명시화, **클러스터 평균·분산** 피처 추가.
4. **template**: 포켓 *찾기* 단계에 template 신호 반영(현재는 검증 전용).
5. **Vina 정의 확정**: "Vina" = gnina의 vina affinity항인지 / 별도 AutoDock Vina 재도킹인지 팀 확인 필요.

---

## 4. 설계 결정 & 반영 내역 (2026-07-27)

팀 방향 확정 후 아래와 같이 개조 완료(GPU 실행은 서버에서 별도 진행).

| 결정 | 반영 위치 |
|---|---|
| **Vina** = gnina의 vina affinity항으로 충분(별도 재도킹 X) | 문서/명칭 명확화(`gnina_affinity`=Vina, `cnn_*`=Gnina) |
| **포켓 결정**: pass/fail 자동판정 유지 + **멤버 평균·분산 병기**(사람 판단) | `3_pocket_validate.py`: `mean_affinity`,`std_affinity`,`mean_cnn_score`,`n_scored` 컬럼. `--members`,`--members-topk`(기본5) |
| **cutoff 조정**: 여러 threshold 진단 헬퍼 | `2_pocket_candidates.py`: `--sweep "6,8,10,12"` (config `pocket_sweep`). 기본 동작 불변 |
| **클러스터 갯수·평균·분산 분석** | `4_ligand_cluster.py`: `rmsd_mean`,`rmsd_std` 컬럼(size=갯수는 기존) |
| **포즈 조합식에 Gnina CNN 추가** | `5_final_select.py`: composite 4항(합의+Boltz2 iptm+Vina aff+**Gnina cnn_score**). P=(0.40,0.25,0.15,0.20)/PA=(0.30,0.15,0.30,0.25). `rmsd_mean/std`,`cnn_affinity` 요약 병기 |
| 배선/기본값 | `run_pipeline.py`: `--members`,`--members-topk`,`--sweep` 배선 + config 기본값 `pocket_members_topk=5`,`pocket_sweep=""` |

| **Boltz2 confidence 명시화** | 2·4·5로 `ligand_iptm`(인터페이스 iptm) 스레딩 → composite Boltz2 항이 `ligand_iptm`(폴백 iptm). `rep_ligand_iptm` 컬럼 |
| **template를 포켓 결정에 반영** | `3_pocket_validate.py`: `--templates`/`--reference` → 템플릿 리간드를 reference에 gemmi 서열정렬 후 포켓 center와 거리 `template_dist` 컬럼. `run_pipeline.py`가 config `templates`/`templates_dir` 있으면 자동 배선 |
| **Boltz2 native confidence 파싱** | `0_rank_poses.py`: `confidence_score`→master_table `confidence` 컬럼(AF3엔 빈값). `4→5`로 `rep_confidence` 스레딩→selection_summary `confidence`(정보용, composite엔 미반영—모델간 스케일 상이) |
| **template_dist 자동 편입** | `3_pocket_validate.py`: `template_dist ≤ template_pass_dist`(기본8Å)→`template_support`→**pass 자동 승격**(pass=phys_pass OR template). `5_final_select.py`: composite 동점 시 **template_dist 작은 포켓 우선**(tie-break). config `template_pass_dist` 배선 |

**하위호환**: 모든 소비자(`9_make_casp_lg` 등)가 `csv.DictReader`(이름 기반) → 컬럼 *추가*는 안전.
기존 컬럼 rename/삭제 없음. 오프라인 단위테스트로 sweep·4항 composite·mean/std·template centroid·ligand_iptm 폴백·tie-break·pass 승격 검증 완료.

### 정책 메모
- **Boltz2 confidence**: composite의 Boltz2 항은 `ligand_iptm`(결합면 신뢰도) 유지. `confidence`(Boltz2 native)는
  모델간 스케일이 달라(AF3는 아예 없음) 정보용 컬럼으로만 노출. 필요 시 모델별 정규화 후 편입 검토.
- **template pass 승격**: 리간드 있는 템플릿은 실험적 근거라 물리 마진과 무관히 pass. 다만 템플릿 리간드가
  타겟과 다른 자리를 가리킬 수 있으므로, `template_support` 컬럼으로 승격 여부를 항상 추적(사람이 검토 가능).
