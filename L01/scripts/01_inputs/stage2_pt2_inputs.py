#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stage2_pt2_inputs.py - CASP L01 stage2 Protenix(pt2) 입력 재생성.

stage2 업데이트(2026-08-08): 구조마다 유기 리간드 N copy(stoichiometry) + Zn 1개.
각 binder의 기존 input.json에서 단백질(proteinChain+MSA 경로)만 그대로 재사용하고,
리간드를 stage2 CSV의 (canonical_smiles, stoichiometry)로 교체(count 반영) + Zn 이온 추가.

Protenix 입력 형식:
  [ { "name": <binder>,
      "sequences": [
        { "proteinChain": { ...BFT1 서열/MSA... } },
        { "ligand": { "ligand": "<SMILES>", "count": <N> } },   # N = stoichiometry
        { "ion":    { "ion": "ZN",         "count": 1 } }        # 모든 구조 Zn 1개
  ]}]
  ※ 이온 형식은 문서 기준 ion.ion="ZN"(CCD_ 접두어 없음). 설치된 Protenix 버전에서 1개 먼저
    검증 후 전체 실행 권장.

- 기존 79개: <existing-root>/<id>/results/pt2/input.json 의 단백질 엔트리 재사용.
- 신규 binder(L010123): 기존 input.json 없음 → --template(아무 기존 input.json)의 단백질 재사용
  (BFT1 동일 수용체라 MSA 공유 가능).
stdlib만. 서버에서 실행(경로가 /path/to 서버 기준이라 서버에서 만들어야 경로가 유효).
"""
import argparse, csv, json, os, sys
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass


def read_stage2_csv(path):
    """stage2 CSV -> {binder: (smiles, stoich)}.
    컬럼: CASP target, ligand code, canonical_smiles, stoichiometry, (ZN, [Zn++], 1)."""
    out = {}
    with open(path, newline="") as f:
        rd = csv.reader(f)
        next(rd, None)                                   # 헤더
        for row in rd:
            if not row or not row[0].strip().startswith("L01"):
                continue
            binder = row[0].strip()
            smiles = row[2].strip()
            try:
                stoich = int(row[3])
            except (IndexError, ValueError):
                stoich = 1
            out[binder] = (smiles, stoich)
    return out


def protein_entries(input_json_path):
    """input.json에서 리간드/이온이 아닌 엔트리(proteinChain 등)만 반환(깊은복사).
    파일이 비었거나 JSON 파싱 실패/서열 없음이면 None(호출부가 --template로 폴백)."""
    try:
        with open(input_json_path) as f:
            data = json.load(f)
    except (json.JSONDecodeError, OSError, ValueError):
        return None
    obj = data[0] if isinstance(data, list) and data else data
    if not isinstance(obj, dict) or "sequences" not in obj:
        return None
    keep = [s for s in obj["sequences"] if "ligand" not in s and "ion" not in s]
    return json.loads(json.dumps(keep)) if keep else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="L01.smiles.stage2.csv")
    ap.add_argument("--existing-root", required=True,
                    help="consensus 루트(기존 <id>/results/pt2/input.json 의 단백질 재사용)")
    ap.add_argument("--template", default="",
                    help="기존 input.json 없는 binder(L010123 등)용 폴백 input.json(아무 기존 것)")
    ap.add_argument("--out-root", required=True, help="새 input.json 출력 루트(<out>/<id>/input.json)")
    ap.add_argument("--no-zn", dest="add_zn", action="store_false", default=True,
                    help="Zn 이온을 넣지 않음(기본은 넣음)")
    args = ap.parse_args()

    ligs = read_stage2_csv(args.csv)
    templ_prot = (protein_entries(args.template)
                  if args.template and os.path.exists(args.template) else None)

    n_ok = n_fallback = 0
    missing, multi = [], 0
    for binder, (smiles, stoich) in sorted(ligs.items()):
        src = os.path.join(args.existing_root, binder, "results", "pt2", "input.json")
        prot = protein_entries(src) if os.path.exists(src) else None
        used_fallback = False
        if prot is None:                                  # 없음/빈파일/깨짐 → 템플릿 폴백
            if templ_prot is None:
                missing.append(binder); continue
            prot = json.loads(json.dumps(templ_prot)); used_fallback = True; n_fallback += 1

        seqs = list(prot)
        seqs.append({"ligand": {"ligand": smiles, "count": stoich}})
        if args.add_zn:
            seqs.append({"ion": {"ion": "ZN", "count": 1}})

        outdir = os.path.join(args.out_root, binder)
        os.makedirs(outdir, exist_ok=True)
        json.dump([{"name": binder, "sequences": seqs}],
                  open(os.path.join(outdir, "pt2.json"), "w"), indent=2)
        n_ok += 1
        if stoich > 1:
            multi += 1
        tag = " [template 폴백]" if used_fallback else ""
        print(f"{binder}: ligand x{stoich}{' +Zn' if args.add_zn else ''}{tag}  {smiles[:34]}")

    print("---")
    print(f"생성 {n_ok}개 (multi-copy {multi}, template 폴백 {n_fallback}) -> {args.out_root}")
    if missing:
        print(f"[주의] 기존 input.json도 --template도 없어 건너뜀: {missing}")


if __name__ == "__main__":
    main()

# ── 실행(서버) ──
#   SC=$CASP17/users/USERNAME
#   python stage2_pt2_inputs.py \
#     --csv $SC/inputs/L01.smiles.stage2.csv \
#     --existing-root $SC/targets/L01/consensus \
#     --template $SC/targets/L01/consensus/L010016/results/pt2/input.json \
#     --out-root  $SC/targets/L01/stage2_pt2_inputs
