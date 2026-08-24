#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1_aggregate.py (Stage2 취합) - 전 fragment 최종 포즈 → poses.csv + submit/<cid>/.

각 fragment의 05_final/selection_summary.csv(≤5 포즈)와 08_casp_lg/<cid>LG_all_models.txt를 모아:
  1) <out-dir>/poses.csv : 전 리간드 ≤5 포즈 한 표(cid, model, pocket, cluster, 점수)
  2) <out-dir>/submit/<cid>/<cid>LG.txt : 각 리간드 최종 제출 포즈 파일
env: 표준 라이브러리만.
"""
import argparse, csv, glob, os, shutil


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", required=True, help=".../targets/L0x/outputs")
    ap.add_argument("--out-dir", default="", help="기본: <outputs>/final/stage2")
    args = ap.parse_args()
    out_dir = args.out_dir or os.path.join(args.outputs, "final", "stage2")
    submit = os.path.join(out_dir, "submit")
    os.makedirs(submit, exist_ok=True)

    pose_rows = []
    n_cid = 0; n_lg = 0
    for sm in sorted(glob.glob(os.path.join(args.outputs, "*", "05_final", "selection_summary.csv"))):
        cid = sm.split(os.sep)[-3]
        n_cid += 1
        for r in csv.DictReader(open(sm)):
            pose_rows.append(dict(cid=cid, model=r.get("model", ""),
                pocket_id=r.get("pocket_id", ""), ligand_cluster_id=r.get("ligand_cluster_id", ""),
                size=r.get("size", ""), iptm=r.get("iptm", ""),
                gnina_affinity=r.get("gnina_affinity", ""), cnn_score=r.get("cnn_score", ""),
                composite=r.get("composite", ""), pocket_pass=r.get("pocket_pass", ""),
                source_cif=r.get("source_cif", "")))
        # 제출 포즈 파일(LG all_models) → submit/<cid>/
        lg = os.path.join(os.path.dirname(sm).replace("05_final", "08_casp_lg"),
                          f"{cid}LG_all_models.txt")
        if os.path.exists(lg):
            d = os.path.join(submit, cid); os.makedirs(d, exist_ok=True)
            shutil.copy(lg, os.path.join(d, f"{cid}LG.txt")); n_lg += 1

    cols = ["cid", "model", "pocket_id", "ligand_cluster_id", "size", "iptm",
            "gnina_affinity", "cnn_score", "composite", "pocket_pass", "source_cif"]
    pcsv = os.path.join(out_dir, "poses.csv")
    with open(pcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(pose_rows)
    print(f"[1_aggregate] {n_cid} 리간드, 포즈 {len(pose_rows)}행 -> {pcsv}")
    print(f"  제출 포즈파일 {n_lg}개 -> {submit}/<cid>/<cid>LG.txt")


if __name__ == "__main__":
    main()

# ── 실행 ──
#   python 1_aggregate.py --outputs .../targets/L01/outputs --out-dir .../targets/L01/outputs/final/stage2
