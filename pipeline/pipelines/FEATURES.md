# pipelines — Features (출력 폴더/컬럼 정의)

각 스텝이 번호별 폴더에 남기는 파일과 핵심 컬럼. 해석(뭘 믿을지)은 [INTERPRETATION.md](INTERPRETATION.md).

## 00_collect/ — master_table.csv
`rank, model, seed, sample, iptm, ligand_iptm, cif` — 흩어진 예측을 한 표로. 이후 모든 스텝 입력.
- `iptm` = interface 신뢰도(chain_pair_iptm off-diagonal 평균). ⚠️ **모델 간 raw 비교 금지**(모델별 스케일 다름).

## 01_protein_clusters/ — protein_clusters.csv, reference.txt, reps/
- 구조→`conf_id`(형태 라벨, 참고용/과분할 가능), `reference.txt`=정렬 기준 구조(PC 중심).

## 02_pocket_candidates/ — pocket_candidates.csv, members.csv
- `pocket_id, size, n_models, models` — **size**=모인 구조 수(합의↑), **n_models**=참여 모델 수(**≥2 & af3 포함이 핵심**).
- members.csv = 어느 예측이 어느 포켓에.

## 03_pocket_validation/ — pocket_validation.csv (+p2rank/, rescore/)
- `p2rank_dist`(작을수록 실제 캐비티), `gnina_affinity`(=**Vina** affinity, 낮을수록↑),
  `cnn_score`/`cnn_affinity`(=**Gnina** CNN), `pass`(대표 기준 통과 여부).
- **멤버 평균/분산**(회의 "포켓 리간드 평균 score"): `mean_affinity`±`std_affinity`, `mean_cnn_score`,
  `n_scored`(채점 멤버 수). → 포켓이 불명확할 때 **평균이 낮은(좋은) 포켓**을 사람이 고르는 신호.
  pass/fail 자동판정은 대표 1개 기준 유지(둘 다 병기).
- **`template_dist`**(회의 "template 찾기"): config `templates`/`templates_dir` 주면, 리간드 포함 템플릿을
  reference에 정렬해 얻은 리간드 위치와 **포켓 center 거리**(Å, 작을수록 그 포켓이 진짜일 가능성↑). 없으면 NA.
  - `template_support`(T/F/NA)=template_dist ≤ `template_pass_dist`(기본 8.0Å)면 True. `phys_pass`=물리(gnina+p2rank)만의 통과.
  - **`pass` = phys_pass OR template_support** — 리간드 있는 템플릿이 가리키는 자리는 물리 마진과 무관하게 통과(자동 편입).

## 04_ligand_clusters/ — ligand_clusters.csv
- 포켓별 리간드 방향 클러스터 대표. 1순위 클러스터가 포켓의 큰 비중이면 방향까지 수렴.
- **응집도**(회의 "갯수·평균·분산"): `size`=멤버 갯수(합의), `rmsd_mean`/`rmsd_std`=클러스터 내
  SC-RMSD 평균·분산(작을수록 방향 일치·안정).

## 05_final/ — selection_summary.csv, model_1~5.cif, SELECTION_RATIONALE.md
- `composite`=종합점수(**합의 + Boltz2(ligand_iptm) + Vina(gnina_affinity) + Gnina(cnn_score)** 가중합).
  - Boltz2 항은 `ligand_iptm`(인터페이스 iptm=결합면 신뢰도) 우선, 없으면 global `iptm` 폴백.
- `iptm`/`ligand_iptm`, `confidence`(Boltz2 native, AF3엔 빈값), `rmsd_mean`/`rmsd_std`(응집도),
  `cnn_affinity`, `pocket_pass`, `posebusters_valid`.
- 정렬: composite 내림차순 + 동점 시 `template_dist` 작은 포켓 우선(template 자동 tie-break).
- 가중치: P=(합의0.40, iptm0.25, aff0.15, cnn0.20) / PA=(0.30, 0.15, 0.30, 0.25). RATIONALE.md에 표로 남김.

## 06_validation/
- `contact_residues.csv`: 잔기별 리간드 접촉빈도(%). 높고 고른 잔기 다수 = 결합부위 수렴.
- `fragment_compare.csv`: exp_ligand별 `n_shared`(예측과 공유 잔기). 진짜 fragment 크고 버퍼 0이 정상.
- `overlay_fragment.pdb`: 예측+실험 리간드 겹친 PDB(PyMOL 육안 확인).

## 07_viz/
- `visualization.png`(4패널: 리간드 PCA by 모델/클러스터, iptm 분포, 상위 클러스터 모델 구성).
- `contact_fingerprint.png`(잔기 접촉빈도; 빨강=실험 공유 잔기).

## 08_casp_lg/ — <TARGET>LG_model1~5.txt, _all_models.txt
- 제출 포맷(수용체 PDB + 리간드 mol block). connectivity는 SMILES, 좌표는 예측.

## 보강 폴더 (config 켤 때)
- `04b_posebusters/posebusters.csv`: `valid`(T/F), `failed`(실패 검사명).
- `05b_refined/refine_summary.csv`: `drift`(이동 Å), `used`(refined/orig).
- `05c_confidence/confidence_report.md`: **HIGH/MEDIUM/LOW** 판정 ⭐ 가장 먼저 볼 것.
