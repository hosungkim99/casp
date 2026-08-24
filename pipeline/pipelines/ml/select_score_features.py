#!/usr/bin/env python3
"""features.csv 에서 핵심 피처만 추려 features_core.csv 생성. stdlib만."""
import argparse, csv

ID_COLS = ["target", "model", "seed", "sample_idx", "cif"]          # 식별자(추적/조인용)
FEAT_COLS = ["ligand_iptm", "iptm", "ptm", "plddt", "gpde", "has_clash",
             "pocket_size", "pocket_n_models", "p2rank_dist",
             "n_pockets", "n_ligand_clusters", "single_pocket"]       # 바로 쓸 핵심 피처
                                                                       # (뒤 3개 = 팀 5k 난이도 proxy)
TODO_COLS = ["posebusters_valid", "ligand_cluster_size", "rmsd_to_truth", "label"]  # 채울 placeholder


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.inp)))
    out_cols = ID_COLS + FEAT_COLS + TODO_COLS
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=out_cols)
        w.writeheader()
        for r in rows:
            row = {c: r.get(c, "") for c in ID_COLS + FEAT_COLS}
            for c in TODO_COLS:
                row[c] = ""
            w.writerow(row)
    print(f"[core] {len(rows)}행 -> {args.out}")


if __name__ == "__main__":
    main()

# $PY "$SC/select_score_features.py" --in "$OUT/features.csv" --out "$OUT/features_core.csv"
# head -2 "$OUT/features_core.csv"