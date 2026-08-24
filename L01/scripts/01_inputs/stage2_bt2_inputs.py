#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage2_bt2_inputs.py - CASP L01 stage2 Boltz2 입력(yaml) 생성 (유기 N copy + Zn).

Boltz2 yaml 형식:
  version: 1
  sequences:
    - protein:
        id: A
        sequence: <BFT1>
        msa: <a3m 경로 또는 empty>
    - ligand:
        id: [B, C, ...]      # id 리스트 = N copy
        smiles: '<SMILES>'
    - ligand:
        id: <다음 letter>
        ccd: ZN              # 이온 = ligand + ccd
각 binder <out-root>/<id>/<id>.yaml 생성. stdlib만. 서버에서 실행.
"""
import argparse, csv, os, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

# BFT1 서열(모든 binder 공통 수용체) — pt2/af3 입력과 동일.
BFT1 = ("VTASIDLQSVSYTDLATQLNDVSDFGKMIILKDNGFNRQVHVSMDKRTKIQLDNENVRLFNGRDKDSTSFILGDEF"
        "AVLRFYRNGESISYIAYKEAQMMNEIAEFYAAPFKKTRAINEKEAFECIYDSRTRSAGKDIVSVKINIDKAKKILN"
        "LPECDYINDYIKTPQVPHGITESQTRAVPSEPKTVYVICLRENGSTIYPNEVSAQMQDAANSVYAVHGLKRYVNFH"
        "FVLYTTEYSCPSGDAKEGLEGFTASLKSNPKAEGYDDQIYFLIRWGTWDNKILGMSWFNSYNVNTASDFEASGMST"
        "TQLMYPGVMAHELGHILGAEHTDNSKDLMYATFTGYLSHLSEKNMDIIAKNLGWEAADGD")


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


def yaml_for(smiles, stoich, msa):
    ids = ", ".join(chr(ord("B") + i) for i in range(stoich))   # B, C, ...
    zn = chr(ord("B") + stoich)
    return (
        "version: 1\n"
        "sequences:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {BFT1}\n"
        f"      msa: {msa}\n"
        "  - ligand:\n"
        f"      id: [{ids}]\n"
        f"      smiles: '{smiles}'\n"
        "  - ligand:\n"
        f"      id: {zn}\n"
        "      ccd: ZN\n"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="L01.smiles.stage2.csv")
    ap.add_argument("--out-root", required=True)
    ap.add_argument("--msa", default="empty",
                    help="protein msa: a3m 절대경로 또는 'empty'(기본). "
                         "품질 위해 BFT1 a3m 경로 권장(pt2/af3와 동일 MSA).")
    args = ap.parse_args()

    ligs = read_stage2_csv(args.csv)
    n_ok = multi = 0
    for b, (smi, st) in sorted(ligs.items()):
        od = os.path.join(args.out_root, b); os.makedirs(od, exist_ok=True)
        open(os.path.join(od, f"{b}.yaml"), "w").write(yaml_for(smi, st, args.msa))
        n_ok += 1
        if st > 1:
            multi += 1
        print(f"{b}: bt2 ligand x{st} + Zn  msa={args.msa if args.msa=='empty' else 'a3m'}  {smi[:30]}")
    print("---")
    print(f"생성 {n_ok} (multi-copy {multi}) -> {args.out_root}")


if __name__ == "__main__":
    main()

# ── 실행(서버) ──
#   SC=$CASP17/users/USERNAME
#   MSA=/path/to/casp17-data/targets_ligand/L01/inputs/A1+L010016/af3_msa/unpaired_A.a3m
#   python3 stage2_bt2_inputs.py --csv $SC/inputs/L01.smiles.stage2.csv \
#     --out-root $SC/targets/L01/stage2_bt2_inputs --msa $MSA
