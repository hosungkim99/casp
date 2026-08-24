#!/usr/bin/env python3
"""
0b_cache_geometry.py - master_table의 모든 cif를 한 번만 파싱해 기하 캐시 저장.

1·2·6·8 스텝은 같은 수천 개 cif를 각자 gemmi로 재파싱했다(gemmi 파싱이 병목 → 실행시간 대부분).
이 스텝이 각 cif를 딱 한 번 파싱해 (단백질 Cα + 리간드 concat 좌표)만 뽑아 저장하고,
이후 스텝들은 complex_io.cached_geometry() 로 재파싱 없이 좌표를 읽는다 → 파싱 4회→1회.

캐시 한 항목: {pk:[(chain,resnum)...], px:(M,3) f64, le:[elem...], lx:(N,3) f64}
 - 좌표만 담으므로 원본 cif 총합보다 작다. (precision은 반드시 float64: 값 불변 = 결과 불변)
 - 파싱 실패/리간드 없는 cif 는 캐시에 넣지 않음 → 그 cif 는 각 스텝이 즉시 파싱으로 폴백(정확성 유지).
boltz2 env (numpy+gemmi). 출력: 00_collect/geom_cache.pkl
"""
import argparse, csv, os, pickle
import numpy as np
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="00의 master_table.csv")
    ap.add_argument("--out", required=True, help="geom_cache.pkl 경로")
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(open(args.table)) if r.get("cif")]
    cache, n_lig = {}, 0
    for i, r in enumerate(rows):
        cif = r["cif"]
        if cif in cache:
            continue
        try:
            ca, ligs = cio.parse_complex(cif)
        except Exception:
            continue
        if not ca:
            continue
        elems, coords = cio.ligand_concat(ligs)
        keys = sorted(ca)
        cache[cif] = {
            "pk": keys,
            # ★ float64 필수: float32로 내리면 값이 미세하게 바뀌어 경계 클러스터의 멤버/대표가
            #    달라지고(→ gnina pass 뒤집힘 → 포켓 선택까지 바뀜) parse_complex 결과와 불일치.
            #    캐시는 "속도만, 값은 그대로"여야 하므로 원본 precision(float64) 유지.
            "px": np.array([ca[k] for k in keys], dtype=np.float64),
            "le": elems,
            "lx": (np.array(coords, dtype=np.float64) if coords else np.zeros((0, 3), np.float64)),
            # copy별 무게중심(원본 프레임, 이온 제외) — 다중 copy 포켓 후보(스텝2)용.
            #   강체정렬은 centroid와 교환가능하므로 원본 프레임 centroid만 저장하면 됨.
            "lc": [[float(x) for x in c] for c in cio.copy_centroids(ligs)],
        }
        if coords:
            n_lig += 1
        if (i + 1) % 500 == 0:
            print(f"  cached {i + 1}/{len(rows)}")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "wb") as f:
        pickle.dump(cache, f, protocol=4)
    print(f"[cache] {len(cache)} structures ({n_lig} with ligand) -> {args.out}")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (boltz2 env python; numpy+gemmi)
#   micromamba run -n boltz2 python 0b_cache_geometry.py \
#       --table $OUT/00_collect/master_table.csv --out $OUT/00_collect/geom_cache.pkl
# (평소엔 run_pipeline.py 가 collect 다음에 자동 호출; 1·2·6·8 이 --cache 로 사용)
