# pipelines — Usage (실행 설명서)

프레임워크(Step 0~4) 계층형 자동화: **포켓 먼저 → 포켓 내 리간드 방향 → 검증 → 최종 5개**.
config에 **경로 4개 + task**만 주면 스텝들이 순서대로 자동 실행되고 결과가 번호별 폴더에 쌓임.

## 1. 준비물
- 팀 cofolding 결과 폴더 (`<model>/seed_*/sample_*/summary.json` + `model.cif`)
- 리간드 목록 tsv (SMILES)
- 환경: `source /path/to/casp17-ligand/scripts/env_setup.sh` (boltz2 env + p2rank/gnina 경로)

## 2. 실행
```bash
cp config/TEMPLATE.conf config/<TARGET>.conf     # 경로 4개 + task 수정
source /path/to/casp17-ligand/scripts/env_setup.sh
python3 run_pipeline.py config/<TARGET>.conf
```

## 3. 유용한 옵션
```bash
python3 run_pipeline.py config/<T>.conf --dry-run                       # 명령만 확인
python3 run_pipeline.py config/<T>.conf --from pocket_validate --force  # 중간부터 재실행
python3 run_pipeline.py config/<T>.conf --only final_select --force     # 한 스텝만
```
- **idempotent**: 산출물 있으면 skip → 중간 실패 후 이어서 진행.
- **되돌아가기**: 결과 불만족(예: 통과 포켓 0개)이면 기준 고치고 `--from <스텝> --force`.

## 4. config 주요 키
- **필수**: `target, results_dir, out_dir, ligand_tsv`
  - (`scripts_dir`는 불필요 — 스텝 경로는 run_pipeline.py 자기 위치 기준으로 자동 해석)
- **켜기**: `task`(P/PA), `python_pb`(PoseBusters), `refine`(정제),
  `experimental_cif`(fragment 비교), `templates`/`templates_dir`(후보별 템플릿 대조=10),
  `split_by_conformation`
  - ⚠️ 템플릿은 **리간드 포함 RCSB 원본**이어야 overlap 계산됨(inputs 리간드 제거본 불가).
- **성능**: `gnina_jobs`(병렬 워커, 기본 1=순차), `gnina_gpus`(라운드로빈 GPU id, 예 `0,1,2`)
- **강건성**: `posebusters_topk`(클러스터당 검사 상위 멤버, 기본 3), `exclude_models`(오염 모델 제외)
- **자동화**: `pocket_auto_relax`(기본 off/opt-in — 통과 0개는 사람이 판단할 신호), `author`
- **회의 반영(2026-07)**: `pocket_members_topk`(포켓별 채점 멤버 수→mean/std affinity·cnn, 기본 5, 0=대표만),
  `pocket_sweep`(cutoff 진단 콤마목록 예 `6,8,10,12` — 포켓 불명확 시 threshold 조정 참고),
  `template_pass_dist`(template_dist 이 값 이하면 포켓 pass 자동 승격, 기본 8.0).
  최종 composite는 4항(합의+Boltz2 ligand_iptm+Vina aff+Gnina cnn), 포켓검증·리간드클러스터에 평균·분산 병기,
  최종 정렬에 template_dist tie-break. master_table에 Boltz2 native `confidence` 파싱(정보용).

## 5. 출력 폴더 → [FEATURES.md](FEATURES.md)
```
00_collect/  01_protein_clusters/  02_pocket_candidates/  03_pocket_validation/
04_ligand_clusters/  05_final/  06_validation/  07_viz/  08_casp_lg/
(+ 04b_posebusters/ 05b_refined/ 05c_confidence/  ← 보강 켤 때)
```

## 6. 범용성 / 한계
- 체인 수·리간드 개수 무관(`common/complex_io`), SMILES 조성 매칭, SC-RMSD(`common/geom`).
- 한계: 동일 조성 리간드 다중 copy는 seqid 순서 고정; 이성질체 구분 X;
  SC-RMSD 대칭순열은 order-match 가정(아니면 plain RMSD 폴백); 템플릿 검색은 로컬(서버 무인터넷).
