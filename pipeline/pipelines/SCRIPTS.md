# pipelines — Scripts (스크립트 사전)

번호 = 실행 순서. 상세 설명·의미는 [OVERVIEW.md](OVERVIEW.md), 해석은 [INTERPRETATION.md](INTERPRETATION.md).

## 지휘자 (pipelines/)
| 파일 | 역할 |
|---|---|
| `run_pipeline.py` | config 하나 읽어 스텝들을 순서대로 자동 호출. idempotent, `--from`/`--only`/`--dry-run`. |
| `run_pipeline.sh` | 위를 slurm(sbatch)로 제출하는 래퍼. |

## core/ — 메인 흐름 (항상 자동)
| 스크립트 | 한 줄 | env |
|---|---|---|
| `0_rank_poses.py` | 수천 예측 수집 → master_table + interface 신뢰도 랭킹 | stdlib |
| `0b_cache_geometry.py` | 전 cif 1회 파싱 → geom_cache (재파싱 제거, 성능) | sci |
| `1_protein_cluster.py` | 단백질 형태 PCA 클러스터 + 기준 구조 선정 | sci |
| `2_pocket_candidates.py` | 리간드 무게중심 클러스터 → 포켓 후보 top10 | sci |
| `3_pocket_validate.py` | p2rank+gnina로 포켓 검증·압축(pass만) | sci |
| `4_ligand_cluster.py` | 포켓 내 SC-RMSD로 리간드 방향 클러스터 | sci |
| `5_final_select.py` | 종합점수(합의·Boltz2 iptm·Vina aff·Gnina cnn)로 최종 ≤5개 | sci |
| `6_contact_residues.py` | 1순위 포켓 결합 잔기 빈도(4Å) | sci |
| `9_make_casp_lg.py` | 최종 5개 → CASP LG 포맷(수용체 PDB+리간드 mol) | sci |

## boost/ — 보강 (config 플래그로 삽입)
| 스크립트 | 켜는 법 | 역할 |
|---|---|---|
| `4b_prep_validity.py` | `python_pb=` | PoseBusters 입력(SDF/PDB) 준비 |
| `4c_posebusters.py` | `python_pb=` | 물리 유효성 검사(충돌·기하) |
| `5b_refine.py` | `refine=true` | gnina `--minimize` 국소 정제(drift 안전장치) |
| `5c_confidence.py` | 항상(맨 끝) | HIGH/MEDIUM/LOW 신뢰도 판정 |

## analysis/ — 검증·시각화·보고
| 스크립트 | 켜는 법 | 역할 |
|---|---|---|
| `7_fragment_compare.py` | `experimental_cif=` | 실험 fragment 접촉잔기 비교(단일 앵커) |
| `10_pose_vs_template.py` | `templates=`/`templates_dir=` | 후보별 템플릿 대조(overlap/jaccard) |
| `8_visualize.py` | 자동 | 4패널 + 접촉지문 PNG |
| `11_make_summary.py` | 자동 | outputs → 정리문서 초안 SUMMARY_DRAFT.md |

## ml/ — ML 실험 (수동·오프라인, 파이프라인과 분리)
| 스크립트 | 역할 |
|---|---|
| `extract_features.py` | 출력 CSV → 포즈 단위 features.csv (`--truth`로 label) |
| `select_score_features.py` | 핵심 피처만 추림 → features_core.csv |
| `train_xgb.py` | XGBoost 학습/평가(AUC·중요도), 타겟 단위 분리 |
