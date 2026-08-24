# -*- coding: utf-8 -*-
# L01 stage2 final pose selection. exosite 판정 = MODEL 리간드 ↔ 자기 수용체 pocket A 거리
# (own-frame, gemmi, super 없음 → contact_residues/7POU-overlap의 신뢰성 문제 회피).
# run:  STAGE2_CONS=<consensus> python stage2_select.py   (gemmi 필요)
import gemmi, csv, os, re, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

ROOT = os.environ.get("STAGE2_CONS",
    "/path/to/casp17-ligand/users/USERNAME/targets/L01/consensus")
POCKA = {188,190,192,230,232,233,254,257,259,260,263,268}
PA_THRESH = 12.0   # 리간드 centroid ↔ pocketA CA centroid < 이 값이면 exosite(pocket A)
EXO = os.path.join(ROOT, "exosite_overlap.csv")

def fnum(x):
    try: return float(x)
    except: return None
def T(x): return str(x).strip().lower() in ("true","1","yes")
def _cen(cs):
    n=len(cs)
    return (sum(c[0] for c in cs)/n, sum(c[1] for c in cs)/n, sum(c[2] for c in cs)/n) if n else None

def pocketA_dist(cif):
    """리간드 무게중심 ↔ 자기 수용체 pocket A(CA) 무게중심 거리. own-frame(중첩 없음)."""
    try:
        st = gemmi.read_structure(cif); st.setup_entities()
    except Exception:
        return None
    lig=[]; pa=[]
    for chain in st[0]:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                if res.seqid.num in POCKA:
                    for at in res:
                        if at.name=="CA": pa.append((at.pos.x,at.pos.y,at.pos.z))
            elif res.name!="HOH":
                for at in res:
                    if at.element.name!="H": lig.append((at.pos.x,at.pos.y,at.pos.z))
    lc=_cen(lig); pc=_cen(pa)
    if not lc or not pc: return None
    return round(sum((a-b)**2 for a,b in zip(lc,pc))**0.5, 1)

# exosite_overlap (info only: 7POU/7POL 거리 — 참고용, 판정엔 안 씀)
exo = {}
if os.path.exists(EXO):
    for r in csv.DictReader(open(EXO)):
        try: exo[(r["binder"], int(r["model_num"]))] = (fnum(r["dist_A"]), r["call"])
        except Exception: pass

def comp(m): return fnum(m["composite"]) or -1
binders = sorted(d for d in os.listdir(ROOT) if re.match(r"L01\d+$", d))
sel=[]; full=[]; no_pipe=[]; exc=[]
for b in binders:
    ss = os.path.join(ROOT, b, "05_final", "selection_summary.csv")
    if not os.path.exists(ss): no_pipe.append(b); continue
    models = list(csv.DictReader(open(ss)))
    if not models: no_pipe.append(b); continue
    for m in models:
        m["_i"] = int(re.search(r"model_(\d+)", m["model"]).group(1))
        m["_pa"] = pocketA_dist(os.path.join(ROOT, b, "05_final", m["model"]+".cif"))
        m["_exo"], m["_call"] = exo.get((b, m["_i"]), (None, "?"))
    def is_exo_pose(m): return m["_pa"] is not None and m["_pa"] < PA_THRESH
    exo_pool = [m for m in models if is_exo_pose(m)]        # 모든 pocket A 포즈(=exosite 최우선)
    if exo_pool:                                           # pocket A 포즈 있으면 그 안에서 PB→composite
        cand = [m for m in exo_pool if T(m["posebusters_valid"])] or exo_pool
    else:                                                  # pocket A 포즈 전무 → 예외
        cand = [m for m in models if T(m["posebusters_valid"])] or models
    ranked = sorted(cand, key=comp, reverse=True)
    p = ranked[0]
    is_exo = is_exo_pose(p)
    site = "A_exosite" if is_exo else "non-A(review)"
    if not is_exo: exc.append((b, p["_pa"]))
    reason = (f"pocket {site}; MODEL1 composite {p['composite']} (size {p['size']}, "
              f"iptm {round(fnum(p['iptm']) or 0,3)}, pocketA_dist {p['_pa']}A, exo7POU {p['_exo']}/{p['_call']}); "
              f"PB {'OK' if T(p['posebusters_valid']) else 'FAIL'}"
              + ("" if is_exo else "  [REVIEW: not in pocket A]"))
    sel.append(dict(binder=b, MODEL1=p["model"], pocket=site, exosite=("Y" if is_exo else "N"),
        composite=p["composite"], size=p["size"], iptm=round(fnum(p["iptm"]) or 0,3),
        posebusters=p["posebusters_valid"], pocketA_dist=p["_pa"], exo7POU_dist=p["_exo"],
        n_models=len(models), cif=os.path.join(ROOT, b, "05_final", p["model"]+".cif"), reason=reason))
    for rk, m in enumerate(ranked, 1):
        full.append(dict(binder=b, submit_rank=rk, model=m["model"], pocket=site,
            size=m["size"], iptm=round(fnum(m["iptm"]) or 0,3), composite=m["composite"],
            posebusters=m["posebusters_valid"], pocketA_dist=m["_pa"], exo7POU_dist=m["_exo"]))

for name, rows in [("stage2_selection.csv", sel), ("stage2_candidates_full.csv", full)]:
    with open(os.path.join(ROOT, name), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

done=len(sel); exoY=sum(1 for r in sel if r["exosite"]=="Y")
print(f"[stage2] 완료 {done}/{len(binders)} (미완 {len(no_pipe)})")
print(f"[stage2] MODEL1 = exosite(pocket A): {exoY}/{done} | 비-A 예외: {len(exc)}")
for b, d in exc: print(f"   예외 {b} -> pocketA_dist {d}A")
print(f"저장: {ROOT}/stage2_selection.csv, stage2_candidates_full.csv")
