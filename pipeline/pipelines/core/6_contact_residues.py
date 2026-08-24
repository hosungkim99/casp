#!/usr/bin/env python3
"""
contact_residues.py (범용, 워크플로 5번) - 1순위 포켓의 리간드-단백질 per-residue 접촉빈도.
complex_io 사용 → 리간드 개수 무관. boltz2 env (numpy+gemmi).

master_table 1순위 대표와 같은 포즈(리간드 concat RMSD <= threshold)인 멤버를 모아,
각 멤버에서 리간드 4Å 안 단백질 잔기를 세어 접촉빈도(%)를 출력. p2rank '포켓 잔기'와
달리 '실제 리간드가 닿는 잔기'를 준다(= 우리 예측 binding site의 정의).
"""
import argparse, csv, os
import numpy as np
import gemmi
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc


def contact_set(path, cutoff):
    """리간드 heavy atom 4Å 안 단백질 잔기 (chain, resnum) 집합."""
    st = gemmi.read_structure(path); m = st[0]
    ns = gemmi.NeighborSearch(m, st.cell, cutoff + 1).populate()
    hits = set()
    for ch in m:
        for r in ch:
            t = gemmi.find_tabulated_residue(r.name)
            is_aa = bool(t) and t.is_amino_acid()
            is_na = bool(t) and t.is_nucleic_acid()
            if is_aa or is_na or r.name in cio.WATER:
                continue
            for a in r:
                if a.element.name == "H":
                    continue
                for mk in ns.find_atoms(a.pos, "\0", radius=cutoff):
                    rr = mk.to_cra(m).residue
                    cc = mk.to_cra(m).chain
                    tt = gemmi.find_tabulated_residue(rr.name)
                    if tt and tt.is_amino_acid():
                        hits.add((cc.name, rr.seqid.num))
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True, help="contact_residues.csv")
    ap.add_argument("--threshold", type=float, default=2.0, help="포켓 멤버 판정 리간드 RMSD")
    ap.add_argument("--cutoff", type=float, default=4.0, help="접촉 거리")
    ap.add_argument("--reference", default="", help="01의 reference.txt(PC중심 구조; 없으면 rank1 폴백)")
    ap.add_argument("--cache", default="", help="0b geom_cache.pkl(있으면 재파싱 없이 사용)")
    ap.add_argument("--pocket-candidates", default="",
                    help="02 pocket_candidates.csv(주면 1순위 포켓 size로 수렴도=멤버/포켓 병기)")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    pocket_size = 0
    if args.pocket_candidates and os.path.exists(args.pocket_candidates):
        sz = [int(r["size"]) for r in csv.DictReader(open(args.pocket_candidates))
              if str(r.get("size", "")).isdigit()]
        pocket_size = max(sz) if sz else 0   # 1순위(최대) 포켓 크기 = 수렴도 분모

    rows = [r for r in csv.DictReader(open(args.table)) if r.get("cif")]
    rows.sort(key=lambda r: int(r["rank"]))

    cache = cio.load_geom_cache(args.cache)
    ref_ca, ref_elems, ref_coords = cio.cached_geometry(cio.reference_cif(rows, args.reference), cache)
    ref_L = np.array(ref_coords)

    members = []
    for r in rows:
        ca, elems, coords = cio.cached_geometry(r["cif"], cache)
        if not ca or elems != ref_elems:
            continue
        common = [k for k in ref_ca if k in ca]
        if len(common) < 50:
            continue
        R, t = kabsch(np.array([ca[k] for k in common]),
                      np.array([ref_ca[k] for k in common]))
        L = (R @ np.array(coords).T).T + t
        if np.sqrt(((L - ref_L) ** 2).sum(1).mean()) <= args.threshold:
            members.append(r["cif"])
    # 수렴도 병기: 이 접촉지문은 'reference pose 근처(RMSD<=threshold)' 멤버만으로 만들어짐.
    # 그 멤버가 포켓 전체의 소수면(=흩어진 타겟) 지문의 대표성이 낮으니 참고 경고를 남긴다.
    conv = (len(members) / pocket_size) if pocket_size else None
    if conv is not None:
        print(f"[contact_residues] 1순위 포켓 멤버 {len(members)}/{pocket_size} = 수렴도 {conv:.0%}")
        if conv < 0.4:
            print("  [참고] 앵커 pose가 포켓의 소수만 대표(흩어진 타겟) → 이 지문의 대표성 낮음. "
                  "후보별 대조(10_pose_vs_template)를 우선 신뢰하세요.")
    else:
        print(f"[contact_residues] 1순위 포켓 멤버 {len(members)}개")

    freq = {}
    for i, cif in enumerate(members):
        for key in contact_set(cif, args.cutoff):
            freq[key] = freq.get(key, 0) + 1
        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(members)}")
    n = max(len(members), 1)

    # 빈도 내림차순, 동점은 (chain, resnum)로 안정 정렬 → 재실행 간 순서 결정적(set 반복순 의존 제거)
    items = sorted(freq.items(), key=lambda kv: (-kv[1], kv[0][0], kv[0][1]))
    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["chain", "residue", "contact_frequency", "n_members", "pocket_size"])
        for (ch, num), cnt in items:
            w.writerow([ch, num, round(cnt / n, 3), n, pocket_size or ""])
    print(f"[contact_residues] -> {args.out}")
    print("  상위 접촉 잔기:", ", ".join(f"{ch}{num}({cnt/n:.0%})" for (ch, num), cnt in items[:12]))


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (먼저: source $CASP17/scripts/env_setup.sh)
#   SC=$CASP17/users/USERNAME/scripts ; OUT=$CASP17/users/USERNAME/targets/T2383
#   micromamba run -n boltz2 python $SC/6_contact_residues.py \
#       --table $OUT/00_collect/master_table.csv --out $OUT/06_validation/contact_residues.csv --reference $OUT/01_protein_clusters/reference.txt
# (평소엔 run_pipeline.py 가 순서대로 자동 호출)
