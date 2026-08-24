#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
s1_oracle_auc.py (Oracle 검증) - 9개 결합확률 정의 중 뭐가 제일 잘 구분하나 AUC로 채점.

oracle 분자(양성=known binder, 음성=decoy)를 Stage1 파이프라인에 태워 나온
05_stage1_binding/binding_row.csv 들을 모아, 9개 정의 각각의 AUC(양성 vs 음성)를 계산.
AUC 1.0=완벽 구분, 0.5=랜덤. 가장 높은 정의 = 제출용으로 추천.
env: 표준 라이브러리만 (sklearn 불필요; Mann-Whitney U로 AUC 정확 계산).
"""
import argparse, csv, glob, os


def rank_norm(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return [0.0] * len(vals)
    order = {x: i for i, x in enumerate(sorted(set(v)))}
    m = (len(order) - 1) or 1
    return [(order[x] / m) if x is not None else 0.0 for x in vals]


def auc(scores, labels):
    """AUC = P(양성점수 > 음성점수). Mann-Whitney U 기반(동점 0.5). 결측은 최저로."""
    s = [(x if x is not None else float("-inf")) for x in scores]
    pos = [s[i] for i in range(len(s)) if labels[i] == 1]
    neg = [s[i] for i in range(len(s)) if labels[i] == 0]
    if not pos or not neg:
        return None
    win = sum((1.0 if a > b else 0.5 if a == b else 0.0) for a in pos for b in neg)
    return win / (len(pos) * len(neg))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", required=True, help="oracle 분자들의 outputs 폴더(상위)")
    ap.add_argument("--ligands", required=True, help="oracle_ligands.tsv (Name, label)")
    args = ap.parse_args()

    lab = {}
    for r in csv.DictReader(open(args.ligands), delimiter="\t"):
        lab[r["Name"]] = int(r["label"])

    rows = []
    for br in sorted(glob.glob(os.path.join(args.outputs, "*", "05_stage1_binding", "binding_row.csv"))):
        try:
            r = next(csv.DictReader(open(br)))
        except StopIteration:
            continue
        cid = r["cid"]
        if cid in lab:
            r["_label"] = lab[cid]
            rows.append(r)
    if not rows:
        raise SystemExit("oracle binding_row.csv 없음 — cofold+파이프라인 먼저 실행")

    labels = [r["_label"] for r in rows]
    boltz = [fnum(r.get("boltz")) for r in rows]
    cnn = [fnum(r.get("cnn")) for r in rows]
    cnnaff = [fnum(r.get("cnnaff")) for r in rows]
    aff = [fnum(r.get("aff")) for r in rows]
    hac = [fnum(r.get("hac")) for r in rows]

    rc = rank_norm(cnn); ra = rank_norm(cnnaff)
    rv = rank_norm([(-a if a is not None else None) for a in aff])
    rb = rank_norm(boltz)
    lc = rank_norm([(cnnaff[i] / hac[i] if (cnnaff[i] is not None and hac[i]) else None) for i in range(len(rows))])
    lv = rank_norm([((-aff[i]) / hac[i] if (aff[i] is not None and hac[i]) else None) for i in range(len(rows))])
    # gnina평균을 재-rank(균일화)해 combined가 boltz와 진짜 50:50이 되게 함(3_aggregate와 동일).
    gnina_avg = [0.5 * rc[i] + 0.5 * ra[i] for i in range(len(rows))]
    gnina_r = rank_norm(gnina_avg)
    defs = {
        "prob_boltz": boltz,
        "prob_cnn": rc, "prob_cnnaff": ra, "prob_vina": rv,
        "prob_gnina": gnina_avg,
        "prob_cons3": [(rc[i] + ra[i] + rv[i]) / 3 for i in range(len(rows))],
        "prob_LE_caf": lc, "prob_LE_vina": lv,
        "prob_combined": [0.5 * rb[i] + 0.5 * gnina_r[i] for i in range(len(rows))],
    }

    npos = sum(labels); nneg = len(labels) - npos
    print(f"[oracle] 양성 {npos}, 음성 {nneg}, 총 {len(rows)}\n")
    res = sorted(((name, auc(sc, labels)) for name, sc in defs.items()),
                 key=lambda x: (x[1] is not None, x[1]), reverse=True)
    print(f"{'정의':14s} AUC   (1.0=완벽, 0.5=랜덤)")
    for name, a in res:
        bar = "#" * int((a or 0) * 20)
        print(f"  {name:14s} {a:.3f}  {bar}" if a is not None else f"  {name:14s}  N/A")
    best = res[0]
    print(f"\n>>> 최고 성능 정의: {best[0]} (AUC {best[1]:.3f})")
    print("    이 정의를 4_finalize --definition 으로 쓰면 oracle 검증된 선택.")


if __name__ == "__main__":
    main()

# ── 실행 (cofold+파이프라인으로 oracle 분자들 05_stage1_binding 만든 뒤) ──
#   python s1_oracle_auc.py --outputs .../oracle/outputs --ligands .../oracle/oracle_ligands.tsv
