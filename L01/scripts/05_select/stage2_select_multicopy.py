#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage2_select_multicopy.py - L01 stage2 multi-copy binder 최종 pose(구조 단위) 선정.

multi-copy(유기 N copy + Zn)는 co-folding이 N copy를 '한꺼번에' 배치하므로,
포켓별 도킹/배향클러스터(단일-리간드용, multi-copy에서 퇴화) 대신 완성 구조를 순위매긴다:

  순위 키 (내림차순):
    1) clash 없음            (has_clash=False 우선)
    2) coverage              = 그 구조의 copy들이 '상위 N개 dominant pocket'을 몇 개(서로 다른) 커버하나
                               (N = copy 수. dominant = 여러 구조가 수렴한 큰 pocket = 합의)
    3) n_in_top              = top-N pocket 안에 든 copy 수(동률 시 세부)
    4) ligand_iptm           = 모든 copy 결합면 confidence 평균(모델 신뢰)
  MODEL 1 = 최상위 구조(= N copy를 합의 자리에 잘 놓고 신뢰 높은 것), 헤지 ≤ top.

왜 이 설계인가: co-folding이 copy들을 이미 배치했으므로 '완성된 답들'을 고르기만 하면 됨.
  포켓 발견(스텝2, per-copy 고침)은 '합의 자리' 정의에만 쓰고, 도킹/클러스터(스텝3·4)는 우회.

입력(각 binder, 파이프라인 스텝2까지 실행된 상태):
  <cons>/<binder>/02_pocket_candidates/pocket_candidates.csv   (pocket_id, size, center...)
  <cons>/<binder>/02_pocket_candidates/members.csv             (cif, pocket_id : copy-point별 배정)
  <cons>/<binder>/00_collect/master_table.csv                  (cif, ligand_iptm, iptm, has_clash)
출력:
  <cons>/stage2_mc_selection.csv     (binder별 MODEL1 + 근거)
  <cons>/stage2_mc_candidates.csv    (binder별 top 후보)
stdlib만.
"""
import argparse, csv, os, re, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass


def fnum(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def is_clash(x):
    return str(x).strip().lower() in ("true", "1", "yes")


def short_id(path):
    """cif 경로 → 'seed<S>/<basename>' 짧은 라벨(같은 basename 다른 seed 구분용)."""
    m = re.search(r"seed[_-]?(\d+)", path)
    base = os.path.basename(path)
    return f"seed{m.group(1)}/{base}" if m else base


def load_binder(cons, b, sub="consensus_s2"):
    # binders 구조: <cons>/<binder>/consensus_s2/{02_pocket_candidates, 00_collect}
    pc = os.path.join(cons, b, sub, "02_pocket_candidates", "pocket_candidates.csv")
    mem = os.path.join(cons, b, sub, "02_pocket_candidates", "members.csv")
    mt = os.path.join(cons, b, sub, "00_collect", "master_table.csv")
    if not (os.path.exists(pc) and os.path.exists(mem) and os.path.exists(mt)):
        return None
    pockets = list(csv.DictReader(open(pc)))
    members = list(csv.DictReader(open(mem)))
    table = {r["cif"]: r for r in csv.DictReader(open(mt)) if r.get("cif")}
    return pockets, members, table


def select_binder(pockets, members, table, n_override=None):
    # 구조별 copy들의 pocket 목록
    by_cif = {}
    for m in members:
        if m.get("cif"):
            by_cif.setdefault(m["cif"], []).append(m["pocket_id"])
    if not by_cif:
        return None
    # N(copy 수): 지정 없으면 구조당 최대 copy-point 수로 추정
    N = n_override or max(len(v) for v in by_cif.values())
    # 상위 N개 dominant pocket(크기순)
    topN = [p["pocket_id"] for p in
            sorted(pockets, key=lambda p: int(p["size"]), reverse=True)[:N]]
    topset = set(topN)
    rows = []
    for cif, plist in by_cif.items():
        tr = table.get(cif, {})
        lig = fnum(tr.get("ligand_iptm"))
        if lig is None:
            lig = fnum(tr.get("iptm")) or 0.0
        covered = set(plist) & topset
        n_in_top = sum(1 for p in plist if p in topset)
        rows.append(dict(cif=cif, n_copies=len(plist), coverage=len(covered),
                         n_in_top=n_in_top, ligand_iptm=round(lig, 4),
                         has_clash=is_clash(tr.get("has_clash")),
                         pockets=";".join(sorted(plist))))
    rows.sort(key=lambda r: (not r["has_clash"], r["coverage"], r["n_in_top"], r["ligand_iptm"]),
              reverse=True)
    return N, topN, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cons", required=True, help="컨테이너 루트(<binder>/<subdir>/02_.., 00_collect 보유)")
    ap.add_argument("--subdir", default="consensus_s2",
                    help="binder 밑 파이프라인 폴더(binders 구조=consensus_s2). 평면구조면 '' 로.")
    ap.add_argument("--only", default="", help="쉼표구분 binder만")
    ap.add_argument("--top", type=int, default=5, help="헤지 포함 후보 수")
    args = ap.parse_args()

    binders = sorted(d for d in os.listdir(args.cons) if re.match(r"L01\d+$", d))
    if args.only:
        want = {x.strip() for x in args.only.split(",") if x.strip()}
        binders = [b for b in binders if b in want]

    sel, cand, skip = [], [], []
    for b in binders:
        loaded = load_binder(args.cons, b, args.subdir)
        if not loaded:
            skip.append(b); continue
        res = select_binder(*loaded)
        if not res:
            skip.append(b); continue
        N, topN, rows = res
        top = rows[0]
        sel.append(dict(binder=b, MODEL1=short_id(top["cif"]),
                        n_copies=N, coverage=f"{top['coverage']}/{N}",
                        n_in_top=top["n_in_top"], ligand_iptm=top["ligand_iptm"],
                        has_clash=top["has_clash"], cif=top["cif"]))
        for rk, r in enumerate(rows[:args.top], 1):
            cand.append(dict(binder=b, submit_rank=rk, cif=short_id(r["cif"]),
                             cif_full=r["cif"],
                             coverage=f"{r['coverage']}/{N}", n_in_top=r["n_in_top"],
                             ligand_iptm=r["ligand_iptm"], has_clash=r["has_clash"],
                             pockets=r["pockets"]))
        print(f"{b}: N={N} top{N}pockets={topN} -> MODEL1 {os.path.basename(top['cif'])} "
              f"cov {top['coverage']}/{N} n_in_top {top['n_in_top']} "
              f"lig_iptm {top['ligand_iptm']} clash {top['has_clash']}")

    for name, rows in [("stage2_mc_selection.csv", sel), ("stage2_mc_candidates.csv", cand)]:
        if not rows:
            continue
        with open(os.path.join(args.cons, name), "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
    print("---")
    print(f"선정 {len(sel)} / 스킵 {len(skip)} -> {args.cons}/stage2_mc_selection.csv (+ _candidates.csv)")
    if skip:
        print(f"  스킵(스텝2 미완): {skip}")


if __name__ == "__main__":
    main()

# ── 실행 ──
#   python3 stage2_select_multicopy.py --cons $CASP17/users/USERNAME/targets/L01/consensus_s2 --only L010462
