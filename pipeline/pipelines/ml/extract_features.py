#!/usr/bin/env python3
"""
extract_features.py - 파이프라인 출력(여러 CSV)을 포즈 단위로 join → ML(XGBoost) 입력 1파일.
정답 구조(--truth)를 주면 각 포즈의 RMSD-to-truth + label(<=2Å=1)까지 계산.

한 행 = 한 포즈(cif). 피처:
  - 신뢰도(per pose): ligand_iptm, iptm, ptm, plddt, gpde, has_clash, ranking_score_native
  - 포켓 배정/포켓레벨: pocket_id, pocket_size, pocket_n_models, p2rank_dist,
    pocket_gnina_affinity, pocket_cnn_score, pocket_pass
  - 난이도 proxy(타겟/포켓레벨): n_pockets(포켓후보수), n_ligand_clusters(그 포켓 방향수),
    single_pocket(포켓1개=1). 팀 5k분석: 포켓수↑·disagreement↑→예측난이도↑, 단일포켓수렴→성공↑.
  - 대표 플래그/유효성: is_pocket_rep, is_ligcluster_rep, ligand_cluster_size,
    posebusters_valid, posebusters_failed
  - 타겟 메타(브로드캐스트): target_verdict, target_dominance, target_separation, target_pose_convergence
  - 라벨(--truth 시): rmsd_to_truth, label
boltz2 env (gemmi+rdkit). 여러 타겟은 --append 로 한 features.csv 에 쌓음.

주의: 포켓레벨 gnina 는 '포켓 대표'의 값이 포켓 내 모든 포즈에 브로드캐스트됨(포켓 내 변별 X).
      포즈마다 gnina 가 필요하면 별도로 25개 전부 scoring 필요(향후).
"""
import argparse, csv, os


def load_csv(path, key=None):
    rows = list(csv.DictReader(open(path))) if os.path.exists(path) else []
    if key:
        return {r[key]: r for r in rows}
    return rows


def g(d, k, default=""):
    return d.get(k, default) if d else default


# ---------- RMSD-to-truth (선택, --truth 줄 때) ----------
def build_ligand_mol(smiles, atoms):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D
    tpl = Chem.RemoveHs(Chem.MolFromSmiles(smiles)); n = tpl.GetNumAtoms()
    if len(atoms) != n:
        raise ValueError(f"원자수 {n} vs {len(atoms)}")
    same = all(tpl.GetAtomWithIdx(i).GetSymbol().capitalize() == el.capitalize()
               for i, (el, _) in enumerate(atoms))
    if same:
        m = Chem.Mol(tpl); conf = Chem.Conformer(n)
        for i, (_, xyz) in enumerate(atoms):
            conf.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
        m.RemoveAllConformers(); m.AddConformer(conf, assignId=True); return m
    from rdkit.Chem import rdDetermineBonds
    rw = Chem.RWMol(); conf = Chem.Conformer(n)
    for i, (el, xyz) in enumerate(atoms):
        rw.AddAtom(Chem.Atom(el.capitalize())); conf.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
    m = rw.GetMol(); m.AddConformer(conf, assignId=True)
    rdDetermineBonds.DetermineConnectivity(m)
    return AllChem.AssignBondOrdersFromTemplate(tpl, m)


def rmsd_to_truth(pose_cif, truth_cif, smiles):
    """단백질 Cα 정렬 후 리간드 RMSD(대칭/원자대응 고려). 실패시 None."""
    import numpy as np
    from rdkit.Chem import rdMolAlign
    from rdkit.Geometry import Point3D
    import complex_io as cio, geom
    ca_p, ligs_p = cio.parse_complex(pose_cif)
    ca_t, ligs_t = cio.parse_complex(truth_cif)
    if not ligs_p or not ligs_t:
        return None
    # 단백질 정렬 (resnum 기준, offset 소폭 탐색)
    by_resnum = lambda ca, off=0: {k[1] + off: v for k, v in ca.items()}
    best = None
    for off in range(-3, 4):
        cp, ct = by_resnum(ca_p, off), by_resnum(ca_t)
        common = sorted(set(cp) & set(ct))
        if len(common) < 50:
            continue
        P = np.array([cp[k] for k in common]); Q = np.array([ct[k] for k in common])
        R, t = geom.kabsch(P, Q)
        r = float(np.sqrt((((R @ P.T).T + t - Q) ** 2).sum(1).mean()))
        if best is None or r < best[0]:
            best = (r, R, t)
    if best is None:
        return None
    _, R, t = best
    try:
        pm = build_ligand_mol(smiles, ligs_p[0]["atoms"])
        tm = build_ligand_mol(smiles, ligs_t[0]["atoms"])
    except Exception:
        return None
    conf = pm.GetConformer()                       # 포즈 리간드를 정답 좌표계로 이동
    for i in range(pm.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        q = R @ np.array([p.x, p.y, p.z]) + t
        conf.SetAtomPosition(i, Point3D(float(q[0]), float(q[1]), float(q[2])))
    try:
        return float(rdMolAlign.CalcRMS(pm, tm))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", required=True, help="한 타겟의 파이프라인 outputs 폴더")
    ap.add_argument("--target", required=True, help="타겟 이름(행에 기록)")
    ap.add_argument("--out", required=True, help="features.csv (없으면 새로, --append면 이어쓰기)")
    ap.add_argument("--truth", default="", help="정답 구조 cif (있으면 RMSD라벨 계산)")
    ap.add_argument("--ligand-tsv", default="", help="RMSD라벨용 SMILES")
    ap.add_argument("--rmsd-cut", type=float, default=2.0)
    ap.add_argument("--append", action="store_true")
    args = ap.parse_args()
    O = args.outputs

    master = load_csv(os.path.join(O, "00_collect/master_table.csv"))
    members = load_csv(os.path.join(O, "02_pocket_candidates/members.csv"), key="cif")
    cand = load_csv(os.path.join(O, "02_pocket_candidates/pocket_candidates.csv"), key="pocket_id")
    val = load_csv(os.path.join(O, "03_pocket_validation/pocket_validation.csv"), key="pocket_id")
    pb = load_csv(os.path.join(O, "04b_posebusters/posebusters.csv"), key="rep_cif")
    ligc = load_csv(os.path.join(O, "04_ligand_clusters/ligand_clusters.csv"), key="rep_cif")
    ligc_list = load_csv(os.path.join(O, "04_ligand_clusters/ligand_clusters.csv"))
    from collections import Counter
    nlc_by_pocket = Counter(r.get("pocket_id", "") for r in ligc_list)  # 포켓별 방향클러스터 수(=disagreement proxy)
    n_pockets = len(cand)   
    conf_rows = load_csv(os.path.join(O, "05c_confidence/confidence.csv"))
    conf = conf_rows[0] if conf_rows else {}
    pocket_reps = {c["rep_cif"] for c in cand.values()}

    smiles = ""
    if args.truth and args.ligand_tsv and os.path.exists(args.ligand_tsv):
        import complex_io as cio
        smis = cio.read_ligand_tsv(args.ligand_tsv)
        smiles = smis[0][2] if smis else ""

    cols = ["target", "model", "seed", "sample_idx", "cif",
            "ligand_iptm", "iptm", "ptm", "plddt", "gpde", "has_clash", "ranking_score_native",
            "pocket_id", "pocket_size", "pocket_n_models", "p2rank_dist",
            "pocket_gnina_affinity", "pocket_cnn_score", "pocket_pass",
            "is_pocket_rep", "is_ligcluster_rep", "ligand_cluster_size",
            "posebusters_valid", "posebusters_failed",
            "target_verdict", "target_dominance", "target_separation", "target_pose_convergence",
            "rmsd_to_truth", "label", "n_pockets", "n_ligand_clusters", "single_pocket"]

    out_rows = []
    for r in master:
        cif = r.get("cif", "")
        if not cif:
            continue
        pid = g(members.get(cif), "pocket_id")
        pc = cand.get(pid, {}); pv = val.get(pid, {})
        pbr = pb.get(cif, {}); lc = ligc.get(cif, {})
        row = {
            "target": args.target, "model": r.get("model"), "seed": r.get("seed"),
            "sample_idx": r.get("sample_idx"), "cif": cif,
            "ligand_iptm": r.get("ligand_iptm"), "iptm": r.get("iptm"), "ptm": r.get("ptm"),
            "plddt": r.get("plddt"), "gpde": r.get("gpde"), "has_clash": r.get("has_clash"),
            "ranking_score_native": r.get("ranking_score_native"),
            "pocket_id": pid, "pocket_size": g(pc, "size"), "pocket_n_models": g(pc, "n_models"),
            "p2rank_dist": g(pv, "p2rank_dist"), "pocket_gnina_affinity": g(pv, "gnina_affinity"),
            "pocket_cnn_score": g(pv, "cnn_score"), "pocket_pass": g(pv, "pass"),
            "n_pockets": n_pockets, "n_ligand_clusters": nlc_by_pocket.get(pid, ""),
            "single_pocket": int(n_pockets == 1),
            "is_pocket_rep": int(cif in pocket_reps),
            "is_ligcluster_rep": int(cif in ligc),
            "ligand_cluster_size": g(lc, "size"),
            "posebusters_valid": g(pbr, "valid"), "posebusters_failed": g(pbr, "failed"),
            "target_verdict": g(conf, "verdict"), "target_dominance": g(conf, "dominance"),
            "target_separation": g(conf, "separation"), "target_pose_convergence": g(conf, "pose_convergence"),
            "rmsd_to_truth": "", "label": "",
        }
        if args.truth and smiles:
            rm = rmsd_to_truth(cif, args.truth, smiles)
            if rm is not None:
                row["rmsd_to_truth"] = round(rm, 3)
                row["label"] = int(rm <= args.rmsd_cut)
        out_rows.append(row)

    write_header = not (args.append and os.path.exists(args.out))
    mode = "a" if args.append and os.path.exists(args.out) else "w"
    with open(args.out, mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        if write_header:
            w.writeheader()
        w.writerows(out_rows)

    npos = sum(1 for r in out_rows if r["label"] == 1)
    print(f"[extract_features] {args.target}: {len(out_rows)}개 포즈 -> {args.out} (mode={mode})")
    if args.truth:
        print(f"  정답(label=1) 포즈: {npos}개 / {len(out_rows)} (recall = 정답 포함 여부)")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (boltz2 env python)
#   $PY $SC/extract_features.py --outputs $OUT --target T2383 --out $OUT/features.csv
#   # 정답 구조 있을 때(distillation): --truth <truth.cif> --ligand-tsv <ligand.tsv>
#   # 여러 타겟 누적: 각 타겟마다 --append 로 한 features.csv 에 쌓기
# $PY "$SC/extract_features.py" --outputs "$OUT" --target T2383 --out "$OUT/features.csv"
# head -2 "$OUT/features.csv"; wc -l "$OUT/features.csv"