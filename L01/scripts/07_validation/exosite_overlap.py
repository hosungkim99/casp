# -*- coding: utf-8 -*-
# exosite_overlap.py — 선정 pose가 BFT의 실험적 allosteric exosite 약물자리와 겹치는지 정량.
#   실험 약물결합 구조(proBFT-3): 7POL(+flumequine=7X9)·7POO/7POQ(+foliosidine=7WK)·7POU(+hesperetin=6JP).
#   각 ref를 우리 구조에 gemmi 단백질정렬 → 실험약물을 우리 프레임으로 → 우리 리간드↔약물 최소거리.
#   작으면(<8Å) 우리 fragment가 실험적으로 검증된 exosite 약물자리에 앉은 것.
# 입력: stage2_mc_selection.csv, ref 폴더(7PO*.pdb/.cif). 출력: exosite_overlap.csv + 요약.
# gemmi만 필요(PyMOL 불필요). 구 exosite_overlap.py(reorg때 삭제) 대체.
import csv, os, sys, math, collections
import gemmi
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

REFIDS = ["7POL", "7POO", "7POQ", "7POU"]
DRUG = {"7X9", "7WK", "6JP"}   # flumequine / foliosidine / hesperetin
NEAR = 8.0


def load_ref(refdir, rid):
    for ext in (".pdb", ".cif"):
        p = os.path.join(refdir, rid + ext)
        if os.path.exists(p):
            st = gemmi.read_structure(p); st.setup_entities()
            drug = [at.pos for ch in st[0] for res in ch if res.name in DRUG
                    for at in res if at.element.name != "H"]
            if drug:
                return st[0][0].get_polymer(), drug
    return None


def organic_ligand_atoms(model):
    out = []
    for ch in model:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if (info and info.is_amino_acid()) or res.name in ("HOH", "WAT", "DOD", "H2O", "ZN"):
                continue
            els = [at.element.name for at in res if at.element.name != "H"]
            if not any(e == "C" for e in els):
                continue
            out += [at.pos for at in res if at.element.name != "H"]
    return out


def main():
    sel = sys.argv[1]; refdir = sys.argv[2]
    refs = []
    for rid in REFIDS:
        r = load_ref(refdir, rid)
        if r:
            refs.append((rid,) + r)
    if not refs:
        sys.exit(f"[에러] ref 약물구조 없음 ({refdir})")
    print(f"참조 약물구조 {len(refs)}개 로드: {[r[0] for r in refs]}")

    out = os.path.join(os.path.dirname(sel), "exosite_overlap.csv")
    rows = []; per = {}; total = 0; err = 0
    for r in csv.DictReader(open(sel, encoding="utf-8")):
        b = r["binder"]; cif = r.get("cif", "")
        if not cif or not os.path.exists(cif):
            continue
        total += 1
        try:
            st = gemmi.read_structure(cif); st.setup_entities()
            op = st[0][0].get_polymer()
            lig = organic_ligand_atoms(st[0])
            if not lig:
                rows.append([b, "NO_LIG", ""]); continue
            best = 1e9; bref = ""
            for rid, rp, drug in refs:
                sup = gemmi.calculate_superposition(op, rp, gemmi.PolymerType.PeptideL,
                                                    gemmi.SupSelect.CaP)
                dt = [sup.transform.apply(dp) for dp in drug]
                for lp in lig:
                    for v in dt:
                        d = math.sqrt((lp.x-v.x)**2 + (lp.y-v.y)**2 + (lp.z-v.z)**2)
                        if d < best:
                            best = d; bref = rid
            rows.append([b, round(best, 2), bref]); per[b] = best
        except Exception as e:
            err += 1; rows.append([b, "ERR", str(e)[:40]])

    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["binder", "min_dist_to_exosite_drug", "nearest_ref"]); w.writerows(rows)

    print(f"\nbinder {total}개 분석(에러 {err}) -> {out}")
    atsite = sum(1 for d in per.values() if d <= NEAR)
    print(f"실험 exosite 약물자리 {NEAR}A 이내(=실험 결합자리 적중): {atsite}/{total}")
    bk = collections.Counter()
    for d in per.values():
        k = ("<3A 정밀중첩" if d < 3 else "3-5A 겹침" if d < 5 else "5-8A 같은자리"
             if d < 8 else "8-15A 근처" if d < 15 else ">15A 딴자리")
        bk[k] += 1
    print("실험약물 최소거리 분포:")
    for k in ["<3A 정밀중첩", "3-5A 겹침", "5-8A 같은자리", "8-15A 근처", ">15A 딴자리"]:
        print(f"  {k}: {bk.get(k, 0)}")


if __name__ == "__main__":
    main()
