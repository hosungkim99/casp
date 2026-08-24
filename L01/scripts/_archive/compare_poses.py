# -*- coding: utf-8 -*-
# Overlay selected final MODEL 1 poses onto the BFT-3 drug-bound structures (7POU/7POL)
# to eyeball whether our poses sit where the experimental drugs (6JP/7X9) bind (the exosite).
# Reads stage2_selection.csv to get each binder's chosen MODEL 1 cif (e.g. L010016 -> model_4).
# Run inside PyMOL:  run /path/to/casp17-local/L01/compare_poses.py
from pymol import cmd
import csv, os

SEL_CSV = r"/path/to/casp17-local/L01/stage2_selection.csv"
# 비교할 binder (편집 가능). 기본 = exosite 정상 4 + rescue 1 + 포켓B 예외 1(대비).
BINDERS = ["L010016", "L010078", "L010223", "L010695", "L010128", "L010319"]
COLORS  = ["cyan", "yellow", "magenta", "orange", "salmon", "red",
           "green", "slate", "purple", "teal"]

cmd.reinitialize()
cmd.bg_color("white")
cmd.set("ray_shadows", 0)

# experimental drug-bound references (fetch needs internet; else load local 7pou.cif/7pol.cif)
cmd.fetch("7pou", "7pou", async_=0)
cmd.fetch("7pol", "7pol", async_=0)

# reference receptor frame = 7pou (faint cartoon)
cmd.hide("everything")
cmd.show("cartoon", "7pou")
cmd.color("grey80", "7pou")
cmd.set("cartoon_transparency", 0.6, "7pou")

# experimental drugs (exosite marker): 6JP=hesperetin(7POU), 7X9=flumequine(7POL)
cmd.select("drugs", "(7pou or 7pol) and resn 6JP+7X9")
cmd.show("sticks", "drugs")
cmd.set("stick_radius", 0.3, "drugs")
cmd.color("green", "drugs")
cmd.util.cnc("drugs")

# --- BFT1 co-folding templates (aligned onto 7pou frame) ---
INP = r"/path/to/casp17-local/L01/inputs/A1+L010016/af3_msa"
TEMPLATES = {
    "t3p24": r"/path/to/casp17-local/templates/L01/3p24.cif",
    "t8wen": INP + "/template_A_01_8wen.cif",
    "t4on1": INP + "/template_A_02_4on1.cif",
    "t8weo": INP + "/template_A_03_8weo.cif",
}
TCOL = {"t3p24": "wheat", "t8wen": "palegreen", "t4on1": "lightblue", "t8weo": "lightpink"}
loaded_t = []
for name, path in TEMPLATES.items():
    if not os.path.exists(path):
        print(f"[skip tmpl] {name}: {path}"); continue
    cmd.load(path, name)
    try:
        cmd.super(name + " and polymer", "7pou and polymer")
    except Exception as e:
        print(f"[warn] super tmpl {name}: {e}")
    cmd.show("cartoon", name)
    cmd.set("cartoon_transparency", 0.85, name)
    cmd.color(TCOL[name], name)
    loaded_t.append(name)
# template catalytic zinc = active-site marker (should sit far from the green drugs)
if loaded_t:
    cmd.select("tmpl_zn", "(" + " or ".join(loaded_t) + ") and elem Zn")
    cmd.show("spheres", "tmpl_zn")
    cmd.color("purple", "tmpl_zn")
    cmd.set("sphere_scale", 0.4, "tmpl_zn")

# map binder -> chosen MODEL1 cif
pick = {}
if os.path.exists(SEL_CSV):
    for r in csv.DictReader(open(SEL_CSV, encoding="utf-8")):
        pick[r["binder"]] = (r.get("cif", "").replace("\\", "/"), r.get("pocket", ""), r.get("MODEL1", ""))

for i, b in enumerate(BINDERS):
    info = pick.get(b)
    if not info or not info[0] or not os.path.exists(info[0]):
        print(f"[skip] {b}: cif 없음 ({info})"); continue
    cif, pocket, mdl = info
    tag = "tmp_" + b
    cmd.load(cif, tag)
    try:
        cmd.super(tag + " and polymer", "7pou and polymer")   # 각 수용체를 7pou 프레임으로
    except Exception as e:
        print(f"[warn] super 실패 {b}: {e}")
    lig = "lig_" + b
    cmd.create(lig, tag + " and not polymer and not solvent")  # 정렬된 리간드만 복사
    cmd.delete(tag)                                            # 수용체는 버림(화면 정리)
    cmd.show("sticks", lig)
    cmd.set("stick_radius", 0.25, lig)
    cmd.color(COLORS[i % len(COLORS)], lig)
    cmd.util.cnc(lig)
    print(f"  {b}: MODEL1={mdl} pocket={pocket} color={COLORS[i % len(COLORS)]}")

cmd.deselect()
cmd.zoom("drugs", 12)
print("\n[compare] 초록=실험약물(exosite). 각 색=우리 포즈. 초록에 겹치면 exosite, 멀면 다른 포켓(L010319=포켓B 예상).")
print("[tip] disable lig_L010319  /  enable lig_L010319  로 개별 토글. BINDERS 리스트 편집해 다른 binder도 비교.")
