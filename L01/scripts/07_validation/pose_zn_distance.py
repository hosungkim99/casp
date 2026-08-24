# -*- coding: utf-8 -*-
# pose_zn_distance.py — 선정 pose의 각 유기 copy가 촉매 Zn에서 얼마나 떨어졌나 측정.
#   BFT1은 아연 metalloprotease → inhibitor는 대개 촉매 Zn 근처 결합. "활성부위 vs exosite" 판정.
#   template 불필요: 선정 cif에 Zn이 이미 들어있음(co-folding에 Zn 넣음).
# 입력: binders/stage2_mc_selection.csv (binder, cif). 출력: pose_zn_distance.csv + 요약.
import csv, os, sys, argparse, collections
import gemmi
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass


def analyze(cif):
    st = gemmi.read_structure(cif); st.setup_entities()
    zns = []; copies = []
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info and info.is_amino_acid():
                continue
            if res.name in ("HOH", "WAT", "DOD", "H2O"):
                continue
            pos = [at.pos for at in res if at.element.name != "H"]
            els = [at.element.name for at in res if at.element.name != "H"]
            if not pos:
                continue
            if res.name == "ZN" or (len(pos) == 1 and els[0] == "Zn"):
                zns.append(pos[0])
            elif any(e == "C" for e in els):
                copies.append(pos)
    return zns, copies


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sel", required=True, help="stage2_mc_selection.csv")
    ap.add_argument("--near", type=float, default=8.0, help="이 거리 이하면 활성부위 근처로 카운트")
    a = ap.parse_args()
    out = os.path.join(os.path.dirname(a.sel), "pose_zn_distance.csv")

    rows = []; per_binder = {}; total = 0; no_zn = 0
    for r in csv.DictReader(open(a.sel, encoding="utf-8")):
        b = r["binder"]; cif = r.get("cif", "")
        if not cif or not os.path.exists(cif):
            continue
        total += 1
        zns, copies = analyze(cif)
        if not zns:
            no_zn += 1; rows.append([b, len(copies), "", "NO_ZN", ""]); continue
        bmin = 1e9
        for i, pos in enumerate(copies):
            mind = min(p.dist(z) for p in pos for z in zns)                 # 원자-원자 최소거리
            cx = sum(p.x for p in pos) / len(pos)
            cy = sum(p.y for p in pos) / len(pos)
            cz = sum(p.z for p in pos) / len(pos)
            cend = min(gemmi.Position(cx, cy, cz).dist(z) for z in zns)     # centroid 거리
            bmin = min(bmin, mind)
            rows.append([b, len(copies), i, round(mind, 2), round(cend, 2)])
        per_binder[b] = bmin

    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["binder", "n_copies", "copy_idx", "min_dist_to_zn", "centroid_dist_to_zn"])
        w.writerows(rows)

    buckets = collections.Counter()
    for d in per_binder.values():
        k = ("<3A 직접배위" if d < 3 else "3-5A" if d < 5 else "5-8A 활성부위"
             if d < 8 else "8-15A" if d < 15 else ">15A exosite")
        buckets[k] += 1
    print(f"binder {total}개 분석 -> {out}")
    near = sum(1 for d in per_binder.values() if d <= a.near)
    print(f"≥1 copy가 촉매 Zn {a.near}A 이내: {near}/{total}")
    if no_zn:
        print(f"Zn 없는 구조: {no_zn}")
    print("촉매 Zn 최소거리 분포 (binder별 가장 가까운 copy):")
    for k in ["<3A 직접배위", "3-5A", "5-8A 활성부위", "8-15A", ">15A exosite"]:
        print(f"  {k}: {buckets.get(k, 0)}")


if __name__ == "__main__":
    main()
