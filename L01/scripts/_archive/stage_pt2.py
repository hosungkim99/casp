#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage_pt2.py - Protenix native 출력을 binder별 consensus 폴더로 이동(정리) + master_table 경로 재작성.

이동 전(Protenix native, 이중폴더):
  <pt2-base>/<id>/                 : ERR, input.json, run.log, <id>/
  <pt2-base>/<id>/<id>/seed_1..6/  : predictions/...

이동 후(요청 구조, 이중폴더 폄):
  <cons-root>/<id>/results/pt2/    : seed_1..6, input.json, run.log
  <cons-root>/<id>/00_collect/master_table.csv  : pt2 행의 cif/summary 경로를 새 위치로 재작성

주의:
  - pt2만 처리(users/USERNAME 내부라 이동 OK). bt2는 공용 read-only(targets/L01/runs)라 손대지 않음.
  - master_table의 bt2 행은 그대로 유지, pt2 행 경로만 문자열 치환.
  - idempotent: 이미 옮겼거나 이미 재작성됐으면 건너뜀.
  - 기본 mode=move('변화'). --mode copy 로 원본 보존 복사도 가능.
stdlib만. Linux(서버)에서 실행.
"""
import argparse, csv, glob, os, shutil, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

EMBEDDED_BINDERS = ("L010016 L010039 L010061 L010069 L010078 L010087 L010099 L010128 L010144 "
 "L010167 L010210 L010223 L010281 L010307 L010309 L010316 L010319 L010322 L010327 L010331 "
 "L010337 L010356 L010363 L010365 L010397 L010412 L010432 L010443 L010462 L010537 L010547 "
 "L010552 L010553 L010589 L010594 L010621 L010626 L010630 L010639 L010649 L010662 L010669 "
 "L010685 L010695 L010702 L010712 L010728 L010738 L010761 L010770 L010782 L010801 L010807 "
 "L010886 L010888 L010906 L010912 L010918 L010919 L010930 L010939 L010943 L010984 L010993 "
 "L011019 L011043 L011057 L011070 L011110 L011124 L011140 L011159 L011160 L011166 L011167 "
 "L011177 L011179 L011199 L011207").split()

EXTRA_FILES = ["input.json", "run.log", "ERR"]   # 이동할 outer 파일(있으면)


def place(src, dst, mode):
    """move 또는 copy. dst 이미 있으면 건너뜀."""
    if not os.path.exists(src):
        return "src-missing"
    if os.path.exists(dst):
        return "dst-exists"
    if mode == "copy":
        (shutil.copytree if os.path.isdir(src) else shutil.copy2)(src, dst)
    else:
        shutil.move(src, dst)
    return "ok"


def rewrite_table(table, old_prefix, new_prefix):
    """master_table.csv의 cif/summary 셀에서 old_prefix -> new_prefix 치환(파일 있으면)."""
    if not os.path.exists(table):
        return 0
    rows = list(csv.DictReader(open(table, newline="")))
    if not rows:
        return 0
    cols = list(rows[0].keys())
    n = 0
    for r in rows:
        for c in ("cif", "summary"):
            v = r.get(c) or ""
            if old_prefix in v:
                r[c] = v.replace(old_prefix, new_prefix); n += 1
    if n:
        with open(table, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pt2-base", required=True, help="protenix msa_all 루트 (.../protenix_test/msa_all)")
    ap.add_argument("--cons-root", required=True, help="consensus 루트 (.../L01/consensus)")
    ap.add_argument("--mode", choices=["move", "copy"], default="move")
    ap.add_argument("--binders", default="", help="binder ID 파일(줄당 1개). 미지정시 내장 79개")
    ap.add_argument("--only", default="", help="쉼표구분 binder만 (검증용)")
    args = ap.parse_args()

    binders = list(EMBEDDED_BINDERS)
    if args.binders and os.path.exists(args.binders):
        binders = [x.strip() for x in open(args.binders) if x.strip().startswith("L01")]
    if args.only:
        want = {x.strip() for x in args.only.split(",") if x.strip()}
        binders = [b for b in binders if b in want]

    tot_seed = tot_extra = tot_rewrite = 0
    for b in binders:
        outer = os.path.join(args.pt2_base, b)          # .../msa_all/<id>
        inner = os.path.join(outer, b)                  # .../msa_all/<id>/<id> (seed_* 보유)
        dst = os.path.join(args.cons_root, b, "results", "pt2")
        os.makedirs(dst, exist_ok=True)

        # 1) seed_* 폴더 이동(이중폴더 폄)
        ns = 0
        for sd in sorted(glob.glob(os.path.join(inner, "seed_*"))):
            if place(sd, os.path.join(dst, os.path.basename(sd)), args.mode) == "ok":
                ns += 1
        # 2) outer의 input.json/run.log/ERR 이동
        ne = 0
        for fn in EXTRA_FILES:
            if place(os.path.join(outer, fn), os.path.join(dst, fn), args.mode) == "ok":
                ne += 1
        # 3) master_table의 pt2 경로 재작성 (old: .../<id>/<id>  ->  new: dst)
        old_prefix = inner
        table = os.path.join(args.cons_root, b, "00_collect", "master_table.csv")
        nr = rewrite_table(table, old_prefix, dst)

        tot_seed += ns; tot_extra += ne; tot_rewrite += nr
        print(f"{b}: seed={ns} extra={ne} table_rows_rewritten={nr}")

    print("---")
    print(f"mode={args.mode}  이동seed={tot_seed}  이동파일={tot_extra}  테이블재작성행={tot_rewrite}  (binder {len(binders)}개)")


if __name__ == "__main__":
    main()
