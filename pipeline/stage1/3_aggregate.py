#!/usr/bin/env python3
"""
3_aggregate.py (Stage1 취합) - 전 fragment binding_row.csv → 결합확률 9정의 binding_scores.csv.

각 fragment outputs의 05_stage1_binding/binding_row.csv(원자료)만 모아서
라이브러리 전체 rank 정규화 → 결합확률 9정의 계산. (오직 outputs 파일 사용)
정의(모두 0~1):
  prob_boltz   = boltz 보정확률(그대로)
  prob_cnn/cnnaff/vina = rank(CNNscore/CNNaffinity/-Vina)
  prob_gnina   = 0.5*cnn + 0.5*cnnaff
  prob_cons3   = (cnn+cnnaff+vina)/3
  prob_LE_caf/LE_vina = rank(강도/heavy atom) — 크기보정
  prob_combined= w_boltz*rank(boltz) + w_gnina*prob_gnina
출력: binding_scores.csv (9정의 + 클러스터/포켓 열). 이후 s4로 bind.txt.
env: 표준 라이브러리만.
"""
import argparse, csv, glob, os


def rank_norm(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return [0.0] * len(vals)
    order = {x: i for i, x in enumerate(sorted(set(v)))}
    m = (len(order) - 1) or 1
    return [(order[x] / m) if x is not None else 0.0 for x in vals]


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", required=True, help=".../targets/L0x/outputs (fragment 폴더들의 상위)")
    ap.add_argument("--out-dir", default="", help="저장 폴더(기본: <outputs>/stage1)")
    ap.add_argument("--w-boltz", type=float, default=0.5)
    ap.add_argument("--w-gnina", type=float, default=0.5)
    args = ap.parse_args()
    out_dir = args.out_dir or os.path.join(args.outputs, "stage1")
    os.makedirs(out_dir, exist_ok=True)

    rows = []
    for br in sorted(glob.glob(os.path.join(args.outputs, "*", "05_stage1_binding", "binding_row.csv"))):
        try:
            r = next(csv.DictReader(open(br)))
        except StopIteration:
            continue
        rows.append(r)
    if not rows:
        raise SystemExit(f"[aggregate] binding_row.csv 없음: {args.outputs}")

    boltz = [fnum(r.get("boltz")) for r in rows]
    cnn = [fnum(r.get("cnn")) for r in rows]
    cnnaff = [fnum(r.get("cnnaff")) for r in rows]
    aff = [fnum(r.get("aff")) for r in rows]
    hac = [fnum(r.get("hac")) for r in rows]

    rc = rank_norm(cnn)
    ra = rank_norm(cnnaff)
    rv = rank_norm([(-a if a is not None else None) for a in aff])
    rb = rank_norm(boltz)
    lc = rank_norm([(cnnaff[i] / hac[i] if (cnnaff[i] is not None and hac[i]) else None) for i in range(len(rows))])
    lv = rank_norm([((-aff[i]) / hac[i] if (aff[i] is not None and hac[i]) else None) for i in range(len(rows))])
    # gnina = 두 rank(cnn·cnnaff)의 평균. 평균은 균일분포가 아니라 두 점수의 상관에 따라
    # 가운데로 쏠린다(상관 낮을수록 심함). 그대로 combined에 넣으면 boltz(균일 rank)와 명목
    # 50:50이 실효로 비대칭이 되고, 그 정도가 타겟마다(상관마다) 달라진다.
    # → 평균을 다시 rank_norm 하여 균일 스케일로 맞춰, 어느 타겟이든 진짜 50:50 blend가 되게 한다.
    gnina_avg = [0.5 * rc[i] + 0.5 * ra[i] for i in range(len(rows))]
    gnina_r = rank_norm(gnina_avg)
    rnd = lambda x: round(x, 4)

    bind_rows, pock_rows = [], []
    for i, r in enumerate(rows):
        gnina = gnina_avg[i]
        cons3 = (rc[i] + ra[i] + rv[i]) / 3
        combined = args.w_boltz * rb[i] + args.w_gnina * gnina_r[i]
        note = r.get("note") or ("" if (boltz[i] is not None or cnn[i] is not None) else "no_score")
        # ① 결합확률 관련
        bind_rows.append(dict(cid=r["cid"],
            prob_boltz=rnd(boltz[i]) if boltz[i] is not None else 0.0,
            prob_gnina=rnd(gnina), prob_cnn=rnd(rc[i]), prob_cnnaff=rnd(ra[i]),
            prob_vina=rnd(rv[i]), prob_cons3=rnd(cons3),
            prob_LE_caf=rnd(lc[i]), prob_LE_vina=rnd(lv[i]), prob_combined=rnd(combined),
            cnn=r.get("cnn", ""), cnnaff=r.get("cnnaff", ""), aff=r.get("aff", ""), note=note))
        # ② 포켓·포즈(클러스터) 관련
        pock_rows.append(dict(cid=r["cid"],
            n_poses=r.get("n_poses", ""), n_pockets=r.get("n_pockets", ""),
            n_pockets_pass=r.get("n_pockets_pass", ""), pocket_sizes=r.get("pocket_sizes", ""),
            dom_pocket_frac=r.get("dom_pocket_frac", ""),
            n_pose_clusters=r.get("n_pose_clusters", ""), pose_sizes=r.get("pose_sizes", ""),
            dom_pose_frac=r.get("dom_pose_frac", ""),
            pocket_residues=r.get("pocket_residues", ""), note=note))

    bcols = ["cid", "prob_boltz", "prob_gnina", "prob_cnn", "prob_cnnaff", "prob_vina",
             "prob_cons3", "prob_LE_caf", "prob_LE_vina", "prob_combined",
             "cnn", "cnnaff", "aff", "note"]
    pcols = ["cid", "n_poses", "n_pockets", "n_pockets_pass", "pocket_sizes", "dom_pocket_frac",
             "n_pose_clusters", "pose_sizes", "dom_pose_frac", "pocket_residues", "note"]
    bind_p = os.path.join(out_dir, "binding_scores.csv")
    pock_p = os.path.join(out_dir, "pocket_clusters.csv")
    with open(bind_p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=bcols); w.writeheader(); w.writerows(bind_rows)
    with open(pock_p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=pcols); w.writeheader(); w.writerows(pock_rows)
    ok = sum(1 for r in bind_rows if not r["note"])
    print(f"[aggregate] {ok}/{len(bind_rows)} 유효")
    print(f"  결합확률 -> {bind_p}")
    print(f"  포켓/클러스터 -> {pock_p}")
    for k in ("prob_boltz", "prob_combined", "prob_cons3"):
        top = sorted([r for r in bind_rows if not r["note"]], key=lambda r: -r[k])[:3]
        print(f"  {k} 상위:", ", ".join(f"{r['cid']}({r[k]})" for r in top))


if __name__ == "__main__":
    main()

# ── 실행 ──
#   python 3_aggregate.py --outputs .../targets/L01/outputs --out .../targets/L01/binding_scores.csv
#   그다음 split_scores.py / s4_make_bind_txt.py
