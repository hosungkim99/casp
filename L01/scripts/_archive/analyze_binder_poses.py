#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
79 binder 각각의 현 단계 최선 pose + pocket 분석 (로컬/서버 공용).

입력(인자로 경로 지정, 기본은 로컬 레이아웃):
  --bindtxt : 선택 pose의 접촉 잔기 (stage-1 제출 .bind.txt: "<ID>\t<score>\t<res,res,...>")
  --poses   : 후보 pose 표 (poses.csv: cid,model,pocket_id,ligand_cluster_id,size,iptm,
              gnina_affinity,cnn_score,composite,pocket_pass,source_cif)  ※ 서버=전체
  --binders : (선택) binding_truth=TRUE csv. 없으면 내장 79개 목록 사용.
  --out     : 결과 csv 경로

서버 실행 예:
  python analyze_binder_poses.py \
    --bindtxt /path/to/casp17-ligand/.../L01/stage1_out/submit/L01LG<GROUP>.bind.txt \
    --poses   /path/to/casp17-ligand/.../users/USERNAME/targets/L01/outputs/final/stage2/poses.csv \
    --out     /path/to/casp17-ligand/.../L01_binder_pose_analysis.csv
poses.csv가 여전히 일부만 있으면(=consolidated 파일이 미완) 서버에서 재생성 필요.

pocket 정의: L010016 p2rank + fragment 수렴으로 얻은 후보 포켓 4개(원래 번호=BFT1 잔기).
"""
import argparse, csv, os, sys
from collections import defaultdict, Counter
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

POCKETS = {
    "A_consensus(188-268)": {188,190,192,230,232,233,254,257,259,260,263,268},
    "B_p2rank1(105-306)":   {105,106,109,110,142,151,153,154,273,306,308},
    "C_p2rank4(61-213)":    {61,62,67,68,69,80,82,210,213,341,342},
    "D_p2rank2(8-146)":     {8,33,34,36,81,83,84,86,88,90,119,120,123,141,143,146},
}

# CASP 수정본 기준 binder 79개(원래 번호). --binders 미지정 시 사용.
EMBEDDED_BINDERS = ("L010016 L010039 L010061 L010069 L010078 L010087 L010099 L010128 L010144 "
 "L010167 L010210 L010223 L010281 L010307 L010309 L010316 L010319 L010322 L010327 L010331 "
 "L010337 L010356 L010363 L010365 L010397 L010412 L010432 L010443 L010462 L010537 L010547 "
 "L010552 L010553 L010589 L010594 L010621 L010626 L010630 L010639 L010649 L010662 L010669 "
 "L010685 L010695 L010702 L010712 L010728 L010738 L010761 L010770 L010782 L010801 L010807 "
 "L010886 L010888 L010906 L010912 L010918 L010919 L010930 L010939 L010943 L010984 L010993 "
 "L011019 L011043 L011057 L011070 L011110 L011124 L011140 L011159 L011160 L011166 L011167 "
 "L011177 L011179 L011199 L011207").split()


def load_binders(path):
    """형식 무관 파서: eval csv(binding_truth 컬럼) / CASP binders 파일(ID만) / 순수 목록 모두 처리.
    파일 없거나 못 읽으면 내장 79개 사용."""
    import re
    pat=re.compile(r"^L01\d{4,}$")
    if path and os.path.exists(path):
        try:
            rows=list(csv.reader(open(path,encoding="utf-8-sig")))
        except Exception:
            rows=[]
        if rows:
            header=[c.strip().lower() for c in rows[0]]
            bind_i=next((i for i,c in enumerate(header) if "bind" in c), None)
            start=1 if any(("id" in c) or ("smiles" in c) or ("bind" in c) for c in header) else 0
            ids=[]
            for r in rows[start:]:
                if not r: continue
                idc=next((c.strip() for c in r if pat.match(c.strip())), None)
                if not idc: continue
                if bind_i is not None and bind_i < len(r):
                    if r[bind_i].strip().upper() in ("TRUE","1","T"): ids.append(idc)
                else:
                    ids.append(idc)
            if ids: return ids
    return list(EMBEDDED_BINDERS)

def load_bindpocket(path):
    d={}
    for line in open(path,encoding="utf-8"):
        x=line.rstrip("\n").split("\t")
        if len(x)>2 and x[2]:
            d[x[0]]=set(int(t[1:]) for t in x[2].split(",") if t[1:].isdigit())
    return d

def load_poses(poses_csv, outputs_dir, binders):
    """후보 pose 로드. --outputs-dir 있으면 per-ligand selection_summary.csv 우선(서버 완본),
    없으면 consolidated poses.csv(부분) 사용."""
    d=defaultdict(list)
    if outputs_dir and os.path.isdir(outputs_dir):
        for b in binders:
            f=os.path.join(outputs_dir, b, "05_final", "selection_summary.csv")
            if os.path.exists(f):
                try:
                    for r in csv.DictReader(open(f,encoding="utf-8")):
                        if r.get("model"): d[b].append(r)
                except Exception:
                    pass
        if d: return d
    if poses_csv and os.path.exists(poses_csv):
        for r in csv.DictReader(open(poses_csv,encoding="utf-8")):
            d[r["cid"]].append(r)
    return d

def assign_pocket(res):
    if not res: return "no-pocket",0
    ov={k:len(res & s) for k,s in POCKETS.items()}
    best=max(ov,key=ov.get)
    return (best if ov[best]>=3 else "?_other"), ov[best]

def main():
    ap=argparse.ArgumentParser()
    base=os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--bindtxt", default=os.path.join(base,"stage1_out/submit/L01LGCHANGE-ME.bind.txt"))
    ap.add_argument("--poses",   default=os.path.join(base,"outputs/final/stage2/poses.csv"),
                    help="consolidated poses.csv (부분). --outputs-dir 있으면 무시됨")
    ap.add_argument("--outputs-dir", default="", dest="outputs_dir",
                    help="per-ligand 출력 루트 (예: .../targets/L01/outputs). 있으면 <id>/05_final/selection_summary.csv 직접 읽음 → 79개 완본")
    ap.add_argument("--binders", default=os.path.join(base,"stage1_all_ligands_eval.csv"))
    ap.add_argument("--out",     default=os.path.join(base,"stage2_binder_pose_analysis.csv"))
    args=ap.parse_args()

    binders=load_binders(args.binders)
    emb=set(EMBEDDED_BINDERS); got=set(binders)
    if got==emb:
        print("[검증] binders 목록 = 내장 79개(CASP 수정본 예측)와 완전 일치")
    else:
        print(f"[주의] binders 목록이 내장 79개와 다름 → 파일에만: {sorted(got-emb)} / 내장에만: {sorted(emb-got)}")
    bp=load_bindpocket(args.bindtxt)
    poses=load_poses(args.poses, args.outputs_dir, binders)
    src="per-ligand selection_summary.csv" if (args.outputs_dir and os.path.isdir(args.outputs_dir)) else "consolidated poses.csv"
    print(f"[입력] binder {len(binders)}개, bind.txt {len(bp)}개, 후보pose {len(poses)}개 binder분 ({src})")

    rows=[]
    for b in binders:
        res=bp.get(b,set())
        pk,ov=assign_pocket(res)
        cands=sorted(poses.get(b,[]),key=lambda r:-float(r["composite"]))
        if cands:
            sizes=[float(c["size"]) for c in cands]; best=cands[0]
            dom=sizes[0]/sum(sizes) if sum(sizes) else 0
            n=len(cands); bs=int(sizes[0]); comp=float(best["composite"])
            iptm=round(float(best["iptm"]),3); gn=float(best["gnina_affinity"]); cnn=float(best["cnn_score"])
            second=int(sizes[1]) if n>1 else 0; has=True
        else:
            n=bs=comp=iptm=gn=cnn=second=None; dom=None; has=False
        if not has: flag="server-needed"
        elif pk.startswith("A") and dom>=0.6 and bs>=5: flag="HIGH"
        elif pk.startswith("?") or (bs is not None and bs<4) or (dom is not None and dom<0.4): flag="REVIEW"
        else: flag="MED"
        rows.append(dict(id=b,pocket=pk,pocket_overlap=ov,flag=flag,n_cand=n,best_size=bs,
                         second_size=second,dominance=None if dom is None else round(dom,2),
                         composite=comp,iptm=iptm,gnina=gn,cnn=cnn,
                         pocket_residues=",".join(f"A{r}" for r in sorted(res))))
    cols=["id","pocket","pocket_overlap","flag","n_cand","best_size","second_size","dominance",
          "composite","iptm","gnina","cnn","pocket_residues"]
    with open(args.out,"w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader()
        for r in rows: w.writerow(r)

    print(f"\n[포켓 분포 — 79개 전체]")
    for k,c in Counter(r["pocket"] for r in rows).most_common():
        print(f"  {k:24s} {c:2d}  {'#'*c}")
    have=[r for r in rows if r["flag"]!="server-needed"]
    print(f"\n[confidence — 후보데이터 있는 {len(have)}개]")
    for k,c in Counter(r["flag"] for r in have).most_common():
        print(f"  {k:14s} {c}")
    if len(have)<len(rows):
        print(f"  server-needed {len(rows)-len(have)}  (poses.csv에 해당 cid 없음 → 서버 전체 poses.csv 필요)")
    print(f"\n저장: {args.out}")

if __name__=="__main__":
    main()
