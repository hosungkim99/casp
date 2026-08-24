# -*- coding: utf-8 -*-
# Measure whether our co-folded ligand poses occupy the EXPERIMENTAL BFT exosite.
# For each binder's final top pose (05_final/model_1.cif): superpose the drug-bound
# experimental structures (7POU=+hesperetin, 7POL=+flumequine) onto our receptor,
# then measure centroid distance between OUR ligand and the EXPERIMENTAL drug.
#   < 8 A  -> EXOSITE (same pocket, overlaps)
#   8-15 A -> near
#   > 15 A -> OTHER (different site)
# Numbering-free: uses structural superposition (super), so BFT-3<->BFT-1 offset does not matter.
#
# Run inside an open PyMOL:   run /path/to/casp17-local/L01/exosite_overlap.py
# Or headless:                pymol -cq /path/to/casp17-local/L01/exosite_overlap.py
# For the full 79: set CONS to the server/synced consensus root (needs the 05_final/model_1.cif files).
from pymol import cmd
import glob, os, math
from collections import defaultdict

# CONS: consensus root. server면 EXO_CONS 환경변수로 덮어씀.
CONS = os.environ.get("EXO_CONS", r"/path/to/casp17-local/L01/consensus")
# REF: 7pou.cif/7pol.cif 있는 폴더(서버 인터넷 없을 때). 비면 fetch 시도.
REF  = os.environ.get("EXO_REF", "")
OUT  = os.path.join(CONS, "exosite_overlap.csv")
DRUG_PDBS = ["7pou", "7pol"]     # BFT-3 exosite complexes
DRUG_RESN = {"6JP", "7X9"}       # actual exosite drugs only: 6JP=hesperetin(7POU), 7X9=flumequine(7POL)
                                 # (whitelist prevents stray PRO/amino-acid hetero being counted as a drug)

def drug_copies(obj):
    """identities (chain,resn,resi) of each exosite-drug copy (7POU has 2x 6JP; 7POL has 7X9)."""
    seen = set()
    for at in cmd.get_model(f"{obj} and resn 6JP+7X9").atom:
        seen.add((at.chain, at.resn, at.resi))
    return sorted(seen)

# preload experimental drug-bound structures once (local cif if given, else fetch)
for p in DRUG_PDBS:
    local = os.path.join(REF, p + ".cif") if REF else ""
    if local and os.path.exists(local):
        cmd.load(local, p)
    else:
        cmd.fetch(p, p, async_=0)

import re
# 각 binder의 최종 제출 모델 전부(model_*.cif, 최대 5개)를 개별 측정.
tasks = []
for f in glob.glob(os.path.join(CONS, "L01*", "05_final", "model_*.cif")):
    mm = re.search(r"model_(\d+)\.cif$", f)
    if not mm:
        continue
    b = os.path.basename(os.path.dirname(os.path.dirname(f)))
    tasks.append((b, int(mm.group(1)), f))
tasks.sort(key=lambda t: (t[0], t[1]))       # binder -> model idx
print(f"[exosite_overlap] models={len(tasks)} binders={len(set(t[0] for t in tasks))}")

rows = []
for b, idx, mc in tasks:
    cmd.load(mc, "m")
    try:
        lig = cmd.centerofmass("m and not polymer and not solvent")
    except Exception:
        cmd.delete("m"); continue
    best = None
    for p in DRUG_PDBS:
        try:
            cmd.super(f"{p} and polymer", "m and polymer")
        except Exception:
            continue
        for (c, rn, ri) in drug_copies(p):     # nearest of the drug copies
            try:
                dcom = cmd.centerofmass(f"{p} and chain {c} and resn {rn} and resi {ri}")
            except Exception:
                continue
            d = math.dist(lig, dcom)
            if best is None or d < best[0]:
                best = (d, p.upper(), rn)
    if best:
        call = "EXOSITE" if best[0] < 8 else ("near" if best[0] < 15 else "OTHER")
        rows.append((b, idx, round(best[0], 1), best[2], best[1], call))
    cmd.delete("m")

rows.sort(key=lambda r: (r[0], r[1]))          # sort by binder, then model_num
with open(OUT, "w") as f:
    f.write("binder,model_num,dist_A,drug,pdb,call\n")
    for r in rows:
        f.write(",".join(map(str, r)) + "\n")

n = len(rows)
ex = sum(1 for r in rows if r[5] == "EXOSITE")
nr = sum(1 for r in rows if r[5] == "near")
print(f"[exosite_overlap] rows={n}  EXOSITE(<8A)={ex}  near(8-15)={nr}  OTHER(>15)={n-ex-nr}")
for r in rows:
    print(f"  {r[0]}_{r[1]}  {r[2]:>5} A  ({r[3]}/{r[4]})  {r[5]}")
print("[exosite_overlap] saved:", OUT)
