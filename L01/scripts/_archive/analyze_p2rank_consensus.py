#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
79 binder의 기존 p2rank 결과를 잔기 기준으로 합산 -> 전체 p2rank 포켓 순위(1~5).
각 binder는 다른 co-fold 수용체지만 잔기 번호는 공통(BFT1) -> 잔기 overlap으로 global 포켓에 매핑.

입력: --outputs-dir <루트> (각 <id>/03_pocket_validation/p2rank/*/ *predictions.csv 읽음)
출력: 콘솔 표 + csv.
"""
import argparse, csv, os, glob, sys
from collections import defaultdict, Counter
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

POCKETS = {
    "A(188-268)": {188,190,192,230,232,233,254,257,259,260,263,268},
    "B(105-306)": {105,106,109,110,142,151,153,154,273,306,308},
    "C(61-213)":  {61,62,67,68,69,80,82,210,213,341,342},
    "E(194-275)": {194,195,196,197,198,199,200,272,273,274,275},
    "D(8-146)":   {8,33,34,36,81,83,84,86,88,90,119,120,123,141,143,146},
}
EMBEDDED = ("L010016 L010039 L010061 L010069 L010078 L010087 L010099 L010128 L010144 L010167 "
 "L010210 L010223 L010281 L010307 L010309 L010316 L010319 L010322 L010327 L010331 L010337 "
 "L010356 L010363 L010365 L010397 L010412 L010432 L010443 L010462 L010537 L010547 L010552 "
 "L010553 L010589 L010594 L010621 L010626 L010630 L010639 L010649 L010662 L010669 L010685 "
 "L010695 L010702 L010712 L010728 L010738 L010761 L010770 L010782 L010801 L010807 L010886 "
 "L010888 L010906 L010912 L010918 L010919 L010930 L010939 L010943 L010984 L010993 L011019 "
 "L011043 L011057 L011070 L011110 L011124 L011140 L011159 L011160 L011166 L011167 L011177 "
 "L011179 L011199 L011207").split()

def parse_res(s):
    out=set()
    for t in s.split():
        t=t.strip()
        if "_" in t:
            n=t.split("_")[-1]
            if n.isdigit(): out.add(int(n))
    return out

def assign(res):
    ov={k:len(res & v) for k,v in POCKETS.items()}
    b=max(ov,key=ov.get)
    return (b if ov[b]>=3 else "other"), ov[b]

def find_pred(root, binder):
    pats=[os.path.join(root,binder,"03_pocket_validation","p2rank","*","*predictions.csv")]
    for p in pats:
        g=sorted(glob.glob(p))
        if g: return g[0]
    return None

def main():
    ap=argparse.ArgumentParser()
    base=os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--outputs-dir", default=os.path.join(base,"outputs"), dest="outputs_dir")
    ap.add_argument("--out", default=os.path.join(base,"p2rank_consensus.csv"))
    args=ap.parse_args()

    stat=defaultdict(lambda: dict(appear=0, rank1=0, top3=0, scores=[], probs=[], ranks=[]))
    other_res=Counter(); found=0; missing=[]
    per_binder_top1={}
    for b in EMBEDDED:
        f=find_pred(args.outputs_dir, b)
        if not f or not os.path.exists(f): missing.append(b); continue
        rows=list(csv.reader(open(f)))
        if len(rows)<2: missing.append(b); continue
        found+=1
        hdr=[c.strip().lower() for c in rows[0]]
        def col(name): return next((i for i,c in enumerate(hdr) if c.startswith(name)),None)
        ri,si,pi,resi=col("rank"),col("score"),col("probability"),col("residue_ids")
        if None in (ri,si,pi,resi): continue
        for r in rows[1:]:
            if len(r)<=max(x for x in [ri,si,pi,resi] if x is not None): continue
            try: rank=int(r[ri]); score=float(r[si]); prob=float(r[pi])
            except: continue
            g,ov=assign(parse_res(r[resi]))
            if g=="other": other_res.update(parse_res(r[resi]))
            st=stat[g]; st["appear"]+=1; st["scores"].append(score); st["probs"].append(prob); st["ranks"].append(rank)
            if rank==1: st["rank1"]+=1; per_binder_top1[b]=g
            if rank<=3: st["top3"]+=1

    print(f"[입력] p2rank 있는 binder: {found}/{len(EMBEDDED)}" + (f"  (없음 {len(missing)})" if missing else ""))
    # 전체 순위: p2rank #1로 뽑힌 횟수 기준
    order=sorted(stat.items(), key=lambda kv:(-kv[1]["rank1"], -kv[1]["appear"]))
    print(f"\n[전체 p2rank 포켓 순위]  (기준: {found}개 binder 중 p2rank #1로 뽑힌 횟수)")
    print(f"{'순위':<4}{'포켓':<12}{'#1로':<7}{'top3':<7}{'등장':<6}{'평균score':<10}{'평균prob':<9}{'평균rank'}")
    rows_out=[]
    for i,(k,st) in enumerate(order,1):
        n=len(st["scores"]) or 1
        ms=sum(st["scores"])/n; mp=sum(st["probs"])/n; mr=sum(st["ranks"])/n
        print(f"{i:<4}{k:<12}{st['rank1']:<7}{st['top3']:<7}{st['appear']:<6}{ms:<10.2f}{mp:<9.3f}{mr:.1f}")
        rows_out.append(dict(overall_rank=i,pocket=k,as_p2rank_no1=st["rank1"],in_top3=st["top3"],
                             appearances=st["appear"],mean_score=round(ms,2),mean_prob=round(mp,3),mean_rank=round(mr,1)))
    with open(args.out,"w",newline="",encoding="utf-8") as fo:
        w=csv.DictWriter(fo,fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)
    # 참고: co-folding 분포
    print("\n[참고] co-folding이 리간드를 놓은 포켓 분포: A=48 B=11 C=11 E=8 (79개 중)")
    if other_res:
        print("\n['other'(A-E 밖) p2rank 포켓에서 자주 나온 잔기 top10 = 미정의 6번째 cavity 후보]")
        print("  " + ", ".join(f"A{r}({c})" for r,c in other_res.most_common(10)))
    print(f"\n저장: {args.out}")

if __name__=="__main__":
    main()
