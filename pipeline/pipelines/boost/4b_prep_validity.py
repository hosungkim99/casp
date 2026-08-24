#!/usr/bin/env python3
"""
4b_prep_validity.py (보강 ① 준비단계) - PoseBusters 입력 준비. boltz2 env (gemmi+rdkit).
ligand_clusters.csv 의 각 대표 포즈마다:
  - 단백질 PDB 저장
  - 리간드 SDF(SMILES connectivity + 예측 좌표) 저장
  - manifest.csv 에 (pocket,lig,cif,sdf,pdb) 기록
실제 검사는 4c_posebusters.py(casp_eval env)가 manifest를 읽어 수행 → env 분리(gemmi/posebusters).
출력: 04b_posebusters/inputs/* + manifest.csv
"""
import argparse, csv, os
import gemmi
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio


def build_ligand_mol(smiles, atoms):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Geometry import Point3D
    template = Chem.RemoveHs(Chem.MolFromSmiles(smiles))
    n = template.GetNumAtoms()
    if len(atoms) != n:
        raise ValueError(f"원자 수 불일치 {n} vs {len(atoms)}")
    same = all(template.GetAtomWithIdx(i).GetSymbol().capitalize() == el.capitalize()
               for i, (el, _) in enumerate(atoms))
    if same:
        m = Chem.Mol(template); conf = Chem.Conformer(n)
        for i, (_, xyz) in enumerate(atoms):
            conf.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
        m.RemoveAllConformers(); m.AddConformer(conf, assignId=True)
        return m
    from rdkit.Chem import rdDetermineBonds
    rw = Chem.RWMol(); conf = Chem.Conformer(n)
    for i, (el, xyz) in enumerate(atoms):
        rw.AddAtom(Chem.Atom(el.capitalize())); conf.SetAtomPosition(i, Point3D(*[float(v) for v in xyz]))
    m = rw.GetMol(); m.AddConformer(conf, assignId=True)
    rdDetermineBonds.DetermineConnectivity(m)
    return AllChem.AssignBondOrdersFromTemplate(template, m)


def prep_one(cif, tag, smis, inp):
    """한 포즈의 (리간드 SDF, 단백질 PDB, note) 생성. 실패 시 빈 경로+note."""
    from rdkit import Chem
    _, ligs = cio.parse_complex(cif)
    if not ligs:
        return "", "", "no_ligand"
    match = cio.match_ligand_to_smiles(ligs[0], smis)
    if match is None:
        return "", "", "no_smiles_match"
    try:
        mol = build_ligand_mol(match[2], ligs[0]["atoms"])
        sdf = os.path.join(inp, f"{tag}.sdf")
        w = Chem.SDWriter(sdf); w.write(mol); w.close()
        pdb = os.path.join(inp, f"{tag}_protein.pdb")
        cio.write_receptor_pdb(cif, pdb)
        return sdf, pdb, ""
    except Exception as e:
        return "", "", f"prep_err:{str(e)[:60]}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand-clusters", required=True)
    ap.add_argument("--cluster-members", default="",
                    help="4의 cluster_members.csv(있으면 클러스터당 top-K 멤버 준비 → medoid 강건화)")
    ap.add_argument("--topk", type=int, default=3, help="클러스터당 준비할 상위(중심성) 멤버 수")
    ap.add_argument("--ligand-tsv", required=True)
    ap.add_argument("--out", required=True, help="04b_posebusters 폴더")
    args = ap.parse_args()
    inp = os.path.join(args.out, "inputs"); os.makedirs(inp, exist_ok=True)

    smis = cio.read_ligand_tsv(args.ligand_tsv)
    if not smis:
        raise SystemExit("ligand.tsv 비어있음")

    clusters = [(r["pocket_id"], r["ligand_cluster_id"], r["rep_cif"])
                for r in csv.DictReader(open(args.ligand_clusters))]
    # 클러스터별 멤버(중심성 순). 있으면 top-K, 없으면 대표만(폴백)
    members = {}
    if args.cluster_members and os.path.exists(args.cluster_members):
        for r in csv.DictReader(open(args.cluster_members)):
            members.setdefault((r["pocket_id"], r["ligand_cluster_id"]), []).append(
                (int(r["member_rank"]), r["cif"]))
        for k in members:
            members[k].sort()

    man = []
    for pid, lid, rep_cif in clusters:
        picks = members.get((pid, lid), [])[:args.topk] or [(1, rep_cif)]
        for rank, cif in picks:
            tag = f"P{pid}_L{lid}_m{rank}"
            sdf, pdb, note = prep_one(cif, tag, smis, inp)
            man.append({"pocket_id": pid, "ligand_cluster_id": lid, "member_rank": rank,
                        "rep_cif": cif, "sdf": sdf, "protein_pdb": pdb, "note": note})
            print(f"  {tag}: sdf={'O' if sdf else 'X'} pdb={'O' if pdb else 'X'} {note}")

    mcsv = os.path.join(args.out, "manifest.csv")
    with open(mcsv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["pocket_id", "ligand_cluster_id", "member_rank",
                                          "rep_cif", "sdf", "protein_pdb", "note"])
        w.writeheader(); w.writerows(man)
    print(f"[prep_validity] {len(man)}개 준비(클러스터 {len(clusters)}, top-{args.topk}) -> {mcsv}")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (boltz2 env python; gemmi+rdkit)
#   <python_sci> $SC/4b_prep_validity.py \
#       --ligand-clusters $OUT/04_ligand_clusters/ligand_clusters.csv \
#       --ligand-tsv <ligand_tsv> --out $OUT/04b_posebusters
