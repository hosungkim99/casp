#!/usr/bin/env python3
"""
1_rank_poses_boltz.py (Step 0, boltz 어댑터) - boltz2 fragment 출력 → master_table.csv.

기존 0_rank_poses.py는 af3식 <model>/seed_*/sample_*/summary.json 레이아웃을 기대한다.
boltz fragment 출력은 레이아웃·키가 다르므로(아래) 이 어댑터로 흡수한다.
  입력: <frag-dir>/seed_*/boltz_results_*/predictions/<cid>/
          confidence_<cid>_model_*.json  +  <cid>_model_*.cif
  매핑: iptm←iptm, ligand_iptm←ligand_iptm, ptm←ptm, plddt←complex_plddt,
        gpde←complex_pde, ranking_score←confidence_score, has_clash←None
출력은 0_rank_poses.py와 동일한 컬럼의 master_table.csv → 1~4단계 그대로 재사용.
stdlib만. 한 fragment(리간드 1개)당 호출.
"""
import argparse, csv, glob, json, os, re, sys
from statistics import mean


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frag-dir", required=True,
                    help="한 fragment의 runs 폴더 (예: .../targets/L01/runs/L010001)")
    ap.add_argument("--out", required=True, help="00_collect 폴더")
    ap.add_argument("--model-name", default="bt2", help="master_table의 model 컬럼값")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    cid = os.path.basename(os.path.normpath(args.frag_dir))

    rows = []
    pat = os.path.join(args.frag_dir, "seed_*", "boltz_results_*", "predictions", cid,
                       f"confidence_{cid}_model_*.json")
    for cj in glob.glob(pat):
        try:
            d = json.load(open(cj))
        except Exception as e:
            sys.stderr.write(f"skip {cj}: {e}\n"); continue
        mm = re.search(r"_model_(\d+)\.json$", cj)
        sample_idx = int(mm.group(1)) if mm else None
        ms = re.search(r"[/\\]seed_(\d+)[/\\]", cj)
        seed = int(ms.group(1)) if ms else None
        cif = os.path.join(os.path.dirname(cj), f"{cid}_model_{sample_idx}.cif")
        rows.append({"model": args.model_name, "seed": seed, "sample_idx": sample_idx,
                     "iptm": d.get("iptm"),
                     "ligand_iptm": d.get("ligand_iptm", d.get("iptm")),
                     "ptm": d.get("ptm"), "plddt": d.get("complex_plddt"),
                     "gpde": d.get("complex_pde"), "has_clash": None,
                     "ranking_score_native": d.get("confidence_score"),
                     "cif": cif if os.path.exists(cif) else "", "summary": cj})

    def key(r):
        li = r["ligand_iptm"] if r["ligand_iptm"] is not None else (r["iptm"] or 0)
        return (0 if r["has_clash"] else 1, li, r["plddt"] or 0,
                -(r["gpde"] if r["gpde"] is not None else 9))

    rows.sort(key=key, reverse=True)
    for i, r in enumerate(rows, 1):
        r["rank"] = i

    cols = ["rank", "model", "seed", "sample_idx", "ligand_iptm", "iptm", "ptm",
            "plddt", "gpde", "has_clash", "ranking_score_native", "cif", "summary"]
    out_csv = os.path.join(args.out, "master_table.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in cols})
    with open(os.path.join(args.out, "progress.txt"), "w") as f:
        f.write(f"cid: {cid}\ntotal samples: {len(rows)}\n  {args.model_name}: {len(rows)}\n")
    print(f"[rank_poses_boltz] {cid}: {len(rows)} poses -> {out_csv}")


if __name__ == "__main__":
    main()

# ── 단독 실행 ──
#   python 1_rank_poses_boltz.py \
#       --frag-dir .../targets/L01/runs/L010001 --out .../L01/outputs/L010001/00_collect
# (이후 1~4단계는 기존 스크립트 그대로 이 out을 입력으로)
