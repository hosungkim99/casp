#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
make_oracle_yaml.py (Oracle) - 기존 boltz affinity YAML → oracle 분자별 YAML.

L010001.yaml(BFT1 + MSA + affinity 속성)을 틀로, 리간드 smiles 한 줄만 oracle 분자로 교체.
BFT1·MSA·affinity 설정 전부 그대로 재사용(MSA 재계산 불필요).
출력: <out-dir>/<name>.yaml
env: 표준 라이브러리만.
"""
import argparse, csv, os, re


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template-yaml", required=True, help="틀 YAML (예: L010001.yaml)")
    ap.add_argument("--ligands", required=True, help="oracle_ligands.tsv (Name, SMILES)")
    ap.add_argument("--out-dir", required=True, help="출력 YAML 폴더 (예: .../oracle/affinity_yaml)")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    templ = open(args.template_yaml).read()
    if not re.search(r"smiles:\s*'[^']*'", templ):
        raise SystemExit("template에서 smiles 라인을 못 찾음 (형식 확인)")

    n = 0
    for r in csv.DictReader(open(args.ligands), delimiter="\t"):
        name, smi = r["Name"], r["SMILES"]
        # smiles 값만 교체(백슬래시 이슈 없게 lambda 치환)
        y = re.sub(r"smiles:\s*'[^']*'", lambda m: f"smiles: '{smi}'", templ, count=1)
        open(os.path.join(args.out_dir, f"{name}.yaml"), "w").write(y)
        n += 1
        print(f"  {name}: smiles={smi[:40]}")
    print(f"[make_oracle_yaml] {n}개 YAML -> {args.out_dir}")


if __name__ == "__main__":
    main()

# ── 실행 (서버) ──
#   python make_oracle_yaml.py \
#     --template-yaml /path/to/casp17-ligand/.../targets/L01/affinity_yaml/L010001.yaml \
#     --ligands .../oracle/oracle_ligands.tsv --out-dir .../oracle/affinity_yaml
