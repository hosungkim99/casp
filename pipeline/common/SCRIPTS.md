# common — Scripts (라이브러리 함수 설명)

단독 실행 X. 각 스텝이 `import common.xxx` 로 씀. (사용 예 [USAGE.md](USAGE.md))

## complex_io.py — 복합체 파서 / IO (gemmi)
| 함수 | 역할 |
|---|---|
| `parse_complex(path)` | cif/pdb → 단백질(체인·Cα·서열) + 리간드 리스트로 분해. **모든 스텝의 입구.** |
| `reference_cif(rows, ref_file)` | 후보 중 기준 구조 1개 선택(중심/지정). |
| `ligand_concat(ligands)` | 여러 리간드를 한 덩어리 좌표로 합침(정준 순서). |
| `ligand_centroid(ligands)` | 리간드 무게중심 → **포켓 클러스터링 좌표**. |
| `composition(atoms)` / `smiles_composition(smiles)` | 원소 조성 계산(원자 매칭용). |
| `match_ligand_to_smiles(lig, smiles_list)` | 예측 리간드 ↔ 입력 SMILES 대응(조성 기반). |
| `read_ligand_tsv(path)` | 리간드 목록 tsv 로드. |
| `write_ligands_pdb(ligands, out)` | 리간드만 pdb로. |
| `write_receptor_pdb(cif, out)` | **수용체만(리간드·물 제거) pdb로** — gnina 입력. (5곳 중복 → 여기 통합) |
| `load_geom_cache(path)` / `cached_geometry(path, cache)` | 파싱·정렬 결과 **캐시**(0b 스텝이 만든 것 재사용, 반복 파싱 방지). |

## geom.py — 정렬 / RMSD 기하
| 함수 | 역할 |
|---|---|
| `kabsch(P, Q)` | 두 점군 최적 회전 R,t (SVD). |
| `align_to_ref(ref_ca, ca)` | Cα 기준 구조에 정렬(포켓·포즈 비교 전 좌표 통일). |
| `apply_rt(R, t, coords)` | 회전·평행이동 적용. |
| `ligand_automorphisms(smiles, n_atoms)` | 리간드 **대칭(자기동형)** 순열 목록 → 대칭보정 RMSD용. |
| `sc_rmsd(P, Q, perms)` | **대칭보정 RMSD** — 포즈 클러스터링 거리(같은 모양 다른 원자번호 처리). |

## scoring.py — gnina 점수 (subprocess)
| 함수 | 역할 |
|---|---|
| `score_only(prot, lig, gpu, gpu_id)` | gnina `--score_only --cnn_scoring rescore` → CNNscore/CNNaffinity/Vina 파싱. |
| `score_with_fallback(prot, lig, gpu_id)` | GPU 실패 시 CPU 재시도. |
| `parse_gpus(spec)` / `run_jobs(fns, jobs, gpus)` | **opt-in 병렬** — GPU 여러 장에 점수 작업 분배. |
| `_parse(stdout)` | gnina 출력 파싱(내부용). |
