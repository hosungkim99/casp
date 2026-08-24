#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4_finalize.py (Stage1 최종화) - 선택한 결합확률 정의로 제출 bind.txt + 선택 이유 txt.

3_aggregate 결과(binding_scores.csv + pocket_clusters.csv)에서:
  - --definition 열을 결합확률로, pocket_clusters의 pocket_residues를 포켓으로 결합
  - CASP 제출 형식 <target>LG<group>.bind.txt 작성
  - 어떤 정의를 왜 골랐는지 definition_choice.txt 작성(팀 검토/수정용 초안)
env: 표준 라이브러리만.
"""
import argparse, csv, os

DEF_DESC = {
    "prob_boltz": "boltz-2 native affinity 보정확률(유일한 절대확률; 보수적)",
    "prob_gnina": "0.5*rank(CNNscore)+0.5*rank(CNNaffinity)",
    "prob_cnn": "rank(CNNscore) — pose 품질",
    "prob_cnnaff": "rank(CNNaffinity) — 결합강도",
    "prob_vina": "rank(-Vina) — 물리에너지",
    "prob_cons3": "(CNNscore+CNNaffinity+(-Vina))/3 — 도킹 3신호 합의",
    "prob_LE_caf": "rank(CNNaffinity/heavy atom) — 크기보정",
    "prob_LE_vina": "rank(-Vina/heavy atom) — 크기보정(물리)",
    "prob_combined": "0.5*rank(boltz)+0.5*prob_gnina — 독립 두 방법(boltz·gnina) 합의",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage1-dir", required=True, help="binding_scores.csv+pocket_clusters.csv 폴더")
    ap.add_argument("--definition", default="prob_combined", help="제출할 결합확률 정의 열")
    ap.add_argument("--target", default="L01")
    ap.add_argument("--group", default="CHANGE-ME", help="CASP 그룹 코드")
    ap.add_argument("--out-dir", default="", help="기본: <stage1-dir>/../final/stage1")
    ap.add_argument("--bind-threshold", type=float, default=0.0,
                    help="이 값 이하는 non-binder(0). 0이면 포켓 있으면 전부 확률 그대로")
    args = ap.parse_args()
    bind_csv = os.path.join(args.stage1_dir, "binding_scores.csv")
    pock_csv = os.path.join(args.stage1_dir, "pocket_clusters.csv")
    out_dir = args.out_dir or os.path.normpath(os.path.join(args.stage1_dir, "..", "final", "stage1"))
    os.makedirs(out_dir, exist_ok=True)

    bind = {r["cid"]: r for r in csv.DictReader(open(bind_csv))}
    pock = {r["cid"]: r for r in csv.DictReader(open(pock_csv))}
    d = args.definition
    if bind and d not in next(iter(bind.values())):
        raise SystemExit(f"정의 열 없음: {d} (있는 열: {list(next(iter(bind.values())).keys())})")

    out_txt = os.path.join(out_dir, f"{args.target}LG{args.group}.bind.txt")
    n_bind = 0
    with open(out_txt, "w") as f:
        for cid in sorted(bind):
            b = bind[cid]; note = b.get("note", "")
            prob = float(b.get(d) or 0)
            res = (pock.get(cid, {}).get("pocket_residues") or "").strip()
            if (not note) and prob > args.bind_threshold and res:
                f.write(f"{cid}\t{prob:.3f}\t{res}\n"); n_bind += 1
            else:
                f.write(f"{cid}\t0\n")
    print(f"[4_finalize] {len(bind)}개 중 binder {n_bind} -> {out_txt}")

    # 선택 이유 초안
    reason = os.path.join(out_dir, "definition_choice.txt")
    with open(reason, "w", encoding="utf-8") as f:
        f.write(f"[Stage1 결합확률 정의 선택]\n\n")
        f.write(f"제출 정의: {d}\n  = {DEF_DESC.get(d, '(설명 없음)')}\n\n")
        f.write("선택 근거:\n")
        f.write("  - boltz와 gnina는 서로 독립(순위상관 낮음)이라, combined는 둘 다에서 높은\n")
        f.write("    fragment만 상위로 올려 한 방법의 오류에 덜 취약(가장 방어적 순위).\n")
        f.write("  - 절대확률로 해석 가능한 건 prob_boltz 뿐(나머지는 라이브러리 내 상대순위).\n")
        f.write("  - 포켓은 팀 cofold 콘센서스 = 알로스테릭 Site B(BFT-3 논문 PMC9514063과 일치).\n\n")
        f.write("대안:\n")
        f.write("  - boltz 신뢰 낮으면 prob_cons3(순수 도킹 합의).\n")
        f.write("  - 크기편향 보려면 prob_LE_caf/prob_LE_vina.\n\n")
        f.write("주의: 정의 간 우열은 known binder(oracle) 검증 후 확정 권장(현재는 잠정 추천).\n")
        f.write(f"\n(전체 9정의 값은 {os.path.basename(bind_csv)}, 포켓/클러스터는 {os.path.basename(pock_csv)} 참고)\n")
    print(f"  선택 이유 초안 -> {reason}")


if __name__ == "__main__":
    main()

# ── 실행 ──
#   python 4_finalize.py --stage1-dir .../outputs/stage1 --definition prob_combined \
#       --target L01 --group <그룹코드> --out-dir .../outputs/final/stage1
