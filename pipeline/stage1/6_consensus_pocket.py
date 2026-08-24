#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
6_consensus_pocket.py (Stage1) - fragment별 포켓을 합쳐 콘센서스 포켓 도출.

combined_bind.py가 만든 binding_scores.csv의 pocket_residues(각 fragment
cofold pose 5A 잔기)를 전 fragment에 걸쳐 집계.
  - 잔기별 접촉빈도(몇 % fragment가 접촉) 계산
  - --min-frac 이상 접촉한 잔기 = 콘센서스 포켓
출력: consensus_pocket.txt (빈도표 + CASP ChainResnum 한 줄).
env: 표준 라이브러리만.
"""
import argparse, csv, re
from collections import Counter


def res_sort_key(res):
    """'A188' -> ('A',188) 정렬용."""
    m = re.match(r"([A-Za-z]+)(\d+)", res)
    return (m.group(1), int(m.group(2))) if m else (res, 0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="combined_bind.py 출력 binding_scores.csv")
    ap.add_argument("--out", required=True, help="저장할 consensus_pocket.txt")
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="콘센서스 포함 최소 접촉비율(기본 0.5=50%)")
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(open(args.csv, encoding="utf-8")) if not r.get("note")]
    n = len(rows)
    if n == 0:
        raise SystemExit(f"[consensus] 유효 행 없음: {args.csv}")

    cnt = Counter()
    for r in rows:
        for res in r["pocket_residues"].split(","):
            res = res.strip()
            if res:
                cnt[res] += 1

    ordered = sorted(cnt.items(), key=lambda x: -x[1])
    consensus = [res for res, c in ordered if c / n >= args.min_frac]
    consensus.sort(key=res_sort_key)

    with open(args.out, "w", encoding="utf-8") as f:
        f.write(f"# L01 콘센서스 포켓 (팀 cofold pose {n}개 집계)\n")
        f.write(f"# 접촉비율 >= {args.min_frac:.0%} 잔기 = 콘센서스\n#\n")
        f.write("# 잔기\t접촉fragment수\t비율\n")
        for res, c in sorted(ordered, key=lambda x: res_sort_key(x[0])):
            mark = " *" if c / n >= args.min_frac else ""
            f.write(f"{res}\t{c}\t{100*c/n:.1f}%{mark}\n")
        f.write(f"\n# 콘센서스 포켓({len(consensus)}잔기, ChainResnum):\n")
        f.write(" ".join(consensus) + "\n")

    print(f"[consensus] {n} fragment 집계 -> {args.out}")
    print(f"  콘센서스({args.min_frac:.0%}+, {len(consensus)}잔기): {' '.join(consensus)}")


if __name__ == "__main__":
    main()

# --- 실행 (표준 라이브러리) ---
#   python 6_consensus_pocket.py --csv .../L01/combined_out/binding_scores.csv \
#       --out .../L01/combined_out/consensus_pocket.txt --min-frac 0.5
