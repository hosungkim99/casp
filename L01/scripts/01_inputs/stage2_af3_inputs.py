#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage2_af3_inputs.py - CASP L01 stage2 AF3(af3_red) 입력 재생성 (유기 N copy + Zn).

기존 af3 입력 <af3-inputs-base>/A1+<id>/af3_msa/input.json 에서 단백질 엔트리
(sequence + unpairedMsaPath/pairedMsaPath + templates)만 그대로 재사용하고,
리간드를 stage2 CSV의 (SMILES × stoichiometry) + Zn 이온으로 교체.

AF3 형식(dialect=alphafold3):
  "sequences": [
    { "protein": { "id":"A", "sequence":..., "unpairedMsaPath":..., "pairedMsaPath":..., "templates":[...] } },
    { "ligand":  { "id": ["B","C",...],  "smiles": "<SMILES>" } },   # id 리스트 = N copy
    { "ligand":  { "id": "F", "ccdCodes": ["ZN"] } }                 # 이온 = ligand + CCD
  ], "modelSeeds": [1..seeds], "dialect":"alphafold3", "version": ...
stdlib만. 서버에서 실행(경로 /path/to 기준).
"""
import argparse, csv, json, os, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass


def read_stage2_csv(path):
    out = {}
    with open(path, newline="") as f:
        rd = csv.reader(f); next(rd, None)
        for row in rd:
            if not row or not row[0].strip().startswith("L01"):
                continue
            try:
                st = int(row[3])
            except (IndexError, ValueError):
                st = 1
            out[row[0].strip()] = (row[2].strip(), st)
    return out


def lig_ids(n):
    """유기 copy용 chain id 리스트 [B, C, ...] (단백질 A 다음부터)."""
    return [chr(ord("B") + i) for i in range(n)]


def build(template_json, name, smiles, stoich, seeds):
    data = json.load(open(template_json))
    # 단백질 등 비-리간드 엔트리 유지(MSA·템플릿 포함), 기존 리간드 제거
    seqs = [s for s in data.get("sequences", []) if "ligand" not in s]
    ids = lig_ids(stoich)
    zn_id = chr(ord("B") + stoich)               # copy 다음 letter
    seqs.append({"ligand": {"id": ids, "smiles": smiles}})
    seqs.append({"ligand": {"id": zn_id, "ccdCodes": ["ZN"]}})
    data["name"] = name
    data["sequences"] = seqs
    data["modelSeeds"] = list(range(1, seeds + 1))
    data.setdefault("dialect", "alphafold3")
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="L01.smiles.stage2.csv")
    ap.add_argument("--af3-inputs-base", required=True,
                    help="기존 af3 입력 루트(.../casp17/targets_ligand/L01/inputs). <base>/A1+<id>/af3_msa/input.json")
    ap.add_argument("--template", default="",
                    help="기존 input.json 없는 binder(L010123)용 폴백 input.json(아무 기존 것)")
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--seeds", type=int, default=6)
    args = ap.parse_args()

    ligs = read_stage2_csv(args.csv)
    n_ok = n_fb = 0
    missing, multi = [], 0
    for b, (smi, st) in sorted(ligs.items()):
        src = os.path.join(args.af3_inputs_base, "A1+" + b, "af3_msa", "input.json")
        tmpl = src if os.path.exists(src) else args.template
        if not tmpl or not os.path.exists(tmpl):
            missing.append(b); continue
        fb = (tmpl == args.template)
        data = build(tmpl, b, smi, st, args.seeds)
        od = os.path.join(args.out_root, b); os.makedirs(od, exist_ok=True)
        json.dump(data, open(os.path.join(od, "af3.json"), "w"), indent=2)
        n_ok += 1; n_fb += fb
        if st > 1:
            multi += 1
        print(f"{b}: af3 ligand x{st} + Zn{' [template 폴백]' if fb else ''}  {smi[:34]}")
    print("---")
    print(f"생성 {n_ok} (multi-copy {multi}, 폴백 {n_fb}) -> {args.out_root}")
    if missing:
        print(f"[주의] 기존 input.json도 --template도 없어 스킵: {missing}")


if __name__ == "__main__":
    main()

# ── 실행(서버) ──
#   SC=$CASP17/users/USERNAME
#   python3 stage2_af3_inputs.py \
#     --csv $SC/inputs/L01.smiles.stage2.csv \
#     --af3-inputs-base /path/to/casp17-data/targets_ligand/L01/inputs \
#     --template /path/to/casp17-data/targets_ligand/L01/inputs/A1+L010016/af3_msa/input.json \
#     --out-root $SC/targets/L01/stage2_af3_inputs
