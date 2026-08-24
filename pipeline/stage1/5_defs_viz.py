#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
5_defs_viz.py (Stage1) - binding_scores.csv 9정의 비교 대시보드.

combined_bind.py가 만든 binding_scores.csv를 읽어 세 패널 생성:
  1) 정의 간 순위상관 (Spearman 9x9)
  2) 상위 N fragment 겹침 (Jaccard 9x9)
  3) 결합확률 >0.5 fragment 수 (정의별 막대)
팀에게 "어떤 정의를 고를지"를 한눈에 보여주기 위한 그림.
env: numpy + scipy + matplotlib (한글폰트 없으면 자동 기본폰트).
"""
import argparse, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager, rcParams
from scipy.stats import spearmanr

# 9개 정의 (라벨, csv 컬럼)
DEFS = [("boltz", "prob_boltz"), ("gnina", "prob_gnina"), ("cnn", "prob_cnn"),
        ("cnnaff", "prob_cnnaff"), ("vina", "prob_vina"), ("cons3", "prob_cons3"),
        ("LE_caf", "prob_LE_caf"), ("LE_vina", "prob_LE_vina"), ("combined", "prob_combined")]


def set_korean_font():
    """한글 폰트 있으면 사용(맑은고딕/나눔), 없으면 기본폰트로 진행.
    폰트 파일을 직접 등록해 캐시 미스로 □ 깨지는 것 방지."""
    import os
    cands = [r"C:\Windows\Fonts\malgun.ttf",
             "/usr/share/fonts/truetype/nanum/NanumGothic.ttf",
             "/Library/Fonts/AppleGothic.ttf"]
    for path in cands:
        if os.path.exists(path):
            try:
                font_manager.fontManager.addfont(path)
                name = font_manager.FontProperties(fname=path).get_name()
                rcParams["font.family"] = name
                break
            except Exception:
                continue
    else:
        for f in ("Malgun Gothic", "NanumGothic", "NanumBarunGothic", "AppleGothic"):
            try:
                font_manager.findfont(f, fallback_to_default=False)
                rcParams["font.family"] = f
                break
            except Exception:
                continue
    rcParams["axes.unicode_minus"] = False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="combined_bind.py 출력 binding_scores.csv")
    ap.add_argument("--out", required=True, help="저장할 png 경로")
    ap.add_argument("--top", type=int, default=50, help="Jaccard 상위 N (기본 50)")
    args = ap.parse_args()
    set_korean_font()

    labels = [d[0] for d in DEFS]
    rows = []
    with open(args.csv, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r.get("note"):   # no_pose / no_score 제외
                continue
            rows.append(r)
    n = len(rows)
    if n == 0:
        raise SystemExit(f"[5_defs_viz] 유효 행 없음: {args.csv}")

    def col(key):
        return np.array([float(r[key]) if r.get(key) not in ("", None) else np.nan for r in rows])
    mats = {name: col(key) for name, key in DEFS}
    k = len(DEFS)

    # Spearman 순위상관
    S = np.eye(k)
    for i in range(k):
        for j in range(i + 1, k):
            rho, _ = spearmanr(mats[labels[i]], mats[labels[j]], nan_policy="omit")
            S[i, j] = S[j, i] = (rho if rho == rho else 0.0)

    # 상위 N Jaccard 겹침
    topsets = [set(np.argsort(-np.nan_to_num(mats[name], nan=-1))[:args.top].tolist()) for name in labels]
    J = np.eye(k)
    for i in range(k):
        for j in range(i + 1, k):
            a, b = topsets[i], topsets[j]
            J[i, j] = J[j, i] = len(a & b) / len(a | b)

    # >0.5 개수
    cnt = [int(np.nansum(mats[name] > 0.5)) for name in labels]

    fig = plt.figure(figsize=(16, 12.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 0.8], hspace=0.28, wspace=0.22)

    def heat(ax, M, title, cmap, vmin, vmax):
        im = ax.imshow(M, cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_xticks(range(k)); ax.set_yticks(range(k))
        ax.set_xticklabels(labels, rotation=40, ha="right"); ax.set_yticklabels(labels)
        for i in range(k):
            for j in range(k):
                v = M[i, j]
                ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=9,
                        color="white" if abs(v) > 0.6 else "black")
        ax.set_title(title, fontsize=14, pad=10)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    heat(fig.add_subplot(gs[0, 0]), S, "정의 간 순위상관 (Spearman)", "RdBu_r", -1, 1)
    heat(fig.add_subplot(gs[0, 1]), J, f"상위{args.top} fragment 겹침 (Jaccard)", "YlGnBu", 0, 1)

    axb = fig.add_subplot(gs[1, :])
    colors = ["tab:red" if labels[i] == "boltz" else "tab:blue" for i in range(k)]
    bars = axb.bar(labels, cnt, color=colors)
    for b, c in zip(bars, cnt):
        axb.text(b.get_x() + b.get_width() / 2, c + max(cnt) * 0.01, str(c), ha="center", fontsize=11)
    axb.set_ylabel("결합확률 >0.5 개수")
    axb.set_title(f"분포: 결합확률 >0.5 인 fragment 수 (전체 {n})   "
                  f"- boltz만 보정확률(보수적), 나머지는 rank(약 절반)", fontsize=13)
    axb.set_ylim(0, max(cnt) * 1.15)

    fig.suptitle(f"결합확률 9정의 비교  (n={n})", fontsize=18, fontweight="bold")
    fig.savefig(args.out, dpi=160, bbox_inches="tight")
    print(f"[5_defs_viz] saved: {args.out} | n={n}")


if __name__ == "__main__":
    main()

# --- 실행 (numpy+scipy+matplotlib 있는 env; boltz2 env 가능) ---
#   python 5_defs_viz.py --csv .../L01/combined_out/binding_scores.csv \
#       --out .../L01/combined_out/binding_defs_comparison.png
