# L01 Results

파이프라인이 만든 실제 산출물 샘플입니다. (대용량 구조 cif는 제외, 요약 CSV만. 경로는 `/path/to/` 로 일반화.)

| 파일 | 내용 | 핵심 수치 |
|---|---|---|
| `stage2_mc_selection_final.csv` | binder 80개 최종 선정 pose (MODEL1) | 80/80, coverage N/N, clash-free |
| `exosite_overlap.csv` | 선정 pose ↔ 실험 약물(exosite) 최소거리 | **78/80이 <3Å** (실험 결합자리 적중) |
| `pose_validation.csv` | Zn-His 배위 + 실험 Zn 일치 + 리간드↔Zn 거리 | Zn 배치 80/80 정확 |
| `posebusters_mmff.csv` | 물리 검증 (gnina+MMFF 정제 후) | 130/132 valid |
| `posebusters_final.csv` | 물리 검증 (대체 pose 교체 후, 최종) | **132/132 valid** |

PoseBusters valid 진행: 114 (초기) → 125 (gnina, clash) → 130 (MMFF, 기하) → **132 (교체, 최종)**.

*수치·서사는 [../README.md](../README.md) 참조.*
