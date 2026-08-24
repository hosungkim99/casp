#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
validate_pockets_p2rank.py - 발견된 포켓을 p2rank cavity로 교차검증 (multi-copy 대응).

기존 step3 버그(포켓 대표구조의 '4-copy 병합중점'을 재서 다 허공 8~15Å) 회피:
  - reference 수용체에 p2rank 1회 실행 → cavity 좌표(reference frame)
  - 각 포켓 center(02_pocket_candidates.csv, 이미 reference frame)를 그 cavity와 비교
  → per-pocket p2rank_dist (포켓별 진짜 값; 병합중점 안 씀, 포켓마다 center 다르니 값도 다름)

self-contained: gemmi + numpy + p2rank(java)만 필요 (pipeline import 불필요).
실행: micromamba run -n boltz2 python validate_pockets_p2rank.py --candidates ... --reference ... --out ...
"""
import argparse, csv, os, subprocess, sys
import numpy as np
import gemmi
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

PRANK = os.environ.get("PRANK",
    "/path/to/casp17-ligand/models/p2rank/p2rank_2.5.1/prank")


def write_receptor_pdb(cif, pdb_out):
    """cif에서 리간드·물 제거한 단백질만 PDB로 (p2rank 입력용)."""
    st = gemmi.read_structure(cif); st.setup_entities()
    st.remove_ligands_and_waters(); st.remove_empty_chains(); st.write_pdb(pdb_out)


def run_prank(pdb, outdir):
    os.makedirs(outdir, exist_ok=True)
    subprocess.run([PRANK, "predict", "-f", pdb, "-o", outdir],
                   check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for root, _, files in os.walk(outdir):
        for f in files:
            if f.endswith("predictions.csv"):
                return os.path.join(root, f)
    return None


def nearest_cavity(pcsv, cen):
    """p2rank predictions.csv에서 cen(3D)에 가장 가까운 cavity → (거리, rank, is_top)."""
    best = None
    with open(pcsv) as fh:
        rd = csv.reader(fh); hdr = [h.strip() for h in next(rd)]
        idx = {h: i for i, h in enumerate(hdr)}
        gi = lambda *ns: next((idx[n] for n in ns if n in idx), None)
        cx, cy, cz, rk = gi("center_x"), gi("center_y"), gi("center_z"), gi("rank")
        if cx is None:
            return None
        for k, row in enumerate(rd):
            if len(row) <= cx:
                continue
            try:
                c = np.array([float(row[cx]), float(row[cy]), float(row[cz])])
            except ValueError:
                continue
            d = float(np.linalg.norm(cen - c))
            if best is None or d < best[0]:
                best = (d, row[rk].strip() if rk is not None else str(k + 1), k == 0)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="02의 pocket_candidates.csv")
    ap.add_argument("--reference", required=True, help="01의 reference.txt (또는 cif 경로)")
    ap.add_argument("--out", required=True, help="출력 폴더")
    ap.add_argument("--pass-dist", type=float, default=6.0, help="이 거리 이하면 druggable cavity 근처(pass)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ref_cif = args.reference
    if ref_cif.endswith(".txt") and os.path.exists(ref_cif):
        txt = open(ref_cif).read().strip().splitlines()
        ref_cif = txt[0].strip() if txt else ""
    if not ref_cif or not os.path.exists(ref_cif):
        sys.exit(f"[에러] reference cif 없음: {ref_cif}")

    prot = os.path.join(args.out, "reference_protein.pdb")
    write_receptor_pdb(ref_cif, prot)
    pcsv = run_prank(prot, os.path.join(args.out, "p2rank"))
    if not pcsv:
        sys.exit("[에러] p2rank 실패 — PRANK 경로/설치 확인 ($PRANK)")
    print(f"[p2rank] reference 수용체 cavity 예측 완료: {os.path.basename(pcsv)}")

    rows = []
    for c in csv.DictReader(open(args.candidates)):
        try:
            cen = np.array([float(c["center_x"]), float(c["center_y"]), float(c["center_z"])])
        except (KeyError, ValueError):
            continue
        nc = nearest_cavity(pcsv, cen)
        d = round(nc[0], 2) if nc else None
        rank = nc[1] if nc else "NA"
        passed = (d is not None) and (d <= args.pass_dist)
        rows.append(dict(pocket_id=c["pocket_id"], size=c["size"],
                         models=c.get("models", ""),
                         center=f"{cen[0]:.1f},{cen[1]:.1f},{cen[2]:.1f}",
                         p2rank_dist=(d if d is not None else "NA"),
                         p2rank_rank=rank, cavity_pass=passed))
        print(f"  P{c['pocket_id']} size {c['size']:>3} ({c.get('models','')}): "
              f"p2rank_dist={d} rank={rank} pass={passed}")

    with open(os.path.join(args.out, "pocket_p2rank.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["pocket_id", "size", "models", "center",
                                          "p2rank_dist", "p2rank_rank", "cavity_pass"])
        w.writeheader(); w.writerows(rows)
    npass = sum(1 for r in rows if r["cavity_pass"])
    print(f"\n[p2rank 검증] {npass}/{len(rows)} 포켓이 cavity 근처(≤{args.pass_dist}Å) "
          f"-> {args.out}/pocket_p2rank.csv")


if __name__ == "__main__":
    main()

# ── 실행 ── (boltz2 env: gemmi+numpy, + p2rank java)
#   D=$CASP17/users/USERNAME/targets/L01/binders/L010462/consensus_s2
#   micromamba run -n boltz2 python validate_pockets_p2rank.py \
#     --candidates $D/02_pocket_candidates/pocket_candidates.csv \
#     --reference  $D/01_protein_clusters/reference.txt \
#     --out        $D/03b_p2rank_check
