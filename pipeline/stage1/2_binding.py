#!/usr/bin/env python3
"""
2_binding.py (Stage1 취합 스텝) - 파이프라인 04 대표 pose로 결합확률 원자료 + 포켓.

한 fragment의 outputs(02·03·04)만으로 Stage1용 결합확률 원자료를 뽑는다.
  1) 지배 포켓 = 통과 포켓 중 size 최대 (통과 0이면 전체 후보 폴백)
  2) 지배 포즈 = 그 포켓의 04 리간드 클러스터 중 size 최대 → 대표 pose(rep_cif)
  3) gnina = 03에서 대표가 같으면 재사용, 다르면 대표 pose로 재계산(파이프라인과 동일 방식)
  4) boltz affinity = fragment의 affinity_*.json seed 평균(원자료를 outputs로 캡처)
  5) 포켓 잔기 = 대표 pose 5Å 이내
출력: <out-dir>/05_stage1_binding/binding_row.csv (원자료 1행; 9정의는 취합 aggregate에서 rank).
env: boltz2/vina_gpu(gemmi+rdkit+numpy) + singularity gnina. complex_io/scoring 재사용.
"""
import argparse, csv, glob, json, os
import numpy as np
import gemmi
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio
import common.scoring as scoring
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


def heavy_atoms(smiles):
    """SMILES 최대 유기조각 heavy atom 수(ligand efficiency용). 실패 None."""
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return None
    fr = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    return max(fr, key=lambda x: x.GetNumHeavyAtoms()).GetNumHeavyAtoms() if fr else None


def protein_pdb(cif, out):
    """리간드·물 제거한 단백질 pdb(gnina 수용체)."""
    cio.write_receptor_pdb(cif, out)


def pocket_res(cif, cutoff=5.0):
    """대표 pose의 리간드 5Å 이내 단백질 잔기 → 'A190' 리스트(ChainResnum)."""
    st = gemmi.read_structure(cif); m = st[0]; prot = []; lig = []
    for ch in m:
        for r in ch:
            t = gemmi.find_tabulated_residue(r.name); aa = bool(t and t.is_amino_acid())
            for a in r:
                if a.element.name == "H":
                    continue
                p = [a.pos.x, a.pos.y, a.pos.z]
                if aa:
                    prot.append((ch.name, r.seqid.num, p))
                elif r.name not in ("HOH", "WAT", "DOD"):
                    lig.append(p)
    if not lig or not prot:
        return []
    L = np.array(lig)
    hit = {(c, n) for c, n, p in prot if np.sqrt(((np.array(p) - L) ** 2).sum(1)).min() <= cutoff}
    return [f"{c}{n}" for c, n in sorted(hit, key=lambda x: (x[0], x[1]))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="fragment outputs 폴더 (.../outputs/<cid>)")
    ap.add_argument("--frag-dir", default="", help="runs/<cid> (affinity json용; 미지정시 rep_cif에서 유추)")
    ap.add_argument("--ligand-tsv", default="", help="SMILES(hac용); 기본 <out-dir>/ligand.tsv")
    ap.add_argument("--cutoff", type=float, default=5.0)
    ap.add_argument("--gpu", default=None)
    args = ap.parse_args()
    O = os.path.normpath(args.out_dir)
    cid = os.path.basename(O)
    out05 = os.path.join(O, "05_stage1_binding"); os.makedirs(out05, exist_ok=True)

    val_csv = os.path.join(O, "03_pocket_validation", "pocket_validation.csv")
    lig_csv = os.path.join(O, "04_ligand_clusters", "ligand_clusters.csv")
    mt_csv = os.path.join(O, "00_collect", "master_table.csv")
    pc_csv = os.path.join(O, "02_pocket_candidates", "pocket_candidates.csv")

    row = dict(cid=cid, boltz="", cnn="", cnnaff="", aff="", hac="",
               n_pockets="", n_pockets_pass="", pocket_sizes="",
               dom_pocket_size="", dom_pocket_frac="",
               n_pose_clusters="", pose_sizes="", dom_pose_size="", dom_pose_frac="",
               n_poses="", pocket_residues="", rep_cif="", note="")

    n_poses = 0
    if os.path.exists(mt_csv):
        n_poses = sum(1 for _ in csv.DictReader(open(mt_csv)))
    row["n_poses"] = n_poses

    if n_poses == 0:
        row["note"] = "no_pose"; _write(out05, row); return
    if not (os.path.exists(val_csv) and os.path.exists(lig_csv)):
        row["note"] = "no_pipeline_output"
        _write(out05, row); return

    val = list(csv.DictReader(open(val_csv)))
    pcs = list(csv.DictReader(open(pc_csv))) if os.path.exists(pc_csv) else []
    row["n_pockets"] = len(pcs) or len(val)
    passed = [r for r in val if str(r.get("pass")) in ("True", "1", "true")]
    row["n_pockets_pass"] = len(passed)
    pool = passed or val
    if not pool:
        row["note"] = "no_pocket"; _write(out05, row); return
    # 포켓 후보 크기들(내림차순) — 02 우선, 없으면 03
    row["pocket_sizes"] = ";".join(str(s) for s in sorted(
        (int(r["size"]) for r in (pcs or val)), reverse=True))

    dom = max(pool, key=lambda r: int(r["size"]))       # 지배 포켓
    dom_pid = dom["pocket_id"]
    row["dom_pocket_size"] = dom["size"]
    row["dom_pocket_frac"] = round(int(dom["size"]) / n_poses, 3) if n_poses else ""

    lcs = [r for r in csv.DictReader(open(lig_csv)) if r["pocket_id"] == dom_pid]
    if not lcs:
        row["note"] = "no_pose_cluster"; _write(out05, row); return
    lcs.sort(key=lambda r: -int(r["size"]))
    dom_cl = lcs[0]                                       # 지배 포즈 클러스터
    rep = dom_cl["rep_cif"]
    row["rep_cif"] = rep
    row["n_pose_clusters"] = len(lcs)
    row["pose_sizes"] = ";".join(r["size"] for r in lcs)   # 지배 포켓 안 포즈 클러스터 크기들
    row["dom_pose_size"] = dom_cl["size"]
    row["dom_pose_frac"] = round(int(dom_cl["size"]) / n_poses, 3) if n_poses else ""

    # 포켓 잔기 (대표 pose 5Å)
    try:
        row["pocket_residues"] = ",".join(pocket_res(rep, args.cutoff))
    except Exception:
        pass

    # gnina: 03 대표와 같으면 재사용, 아니면 대표 pose로 재계산(파이프라인과 동일 prep)
    reuse = (dom.get("rep_cif") == rep and dom.get("cnn_score") not in (None, "", "NA"))
    if reuse:
        row.update(cnn=dom.get("cnn_score"), cnnaff=dom.get("cnn_affinity"),
                   aff=dom.get("gnina_affinity"))
    else:
        work = os.path.join(out05, "work"); os.makedirs(work, exist_ok=True)
        prot = os.path.join(work, "prot.pdb"); lig = os.path.join(work, "lig.pdb")
        try:
            protein_pdb(rep, prot)
            _, ligs = cio.parse_complex(rep); cio.write_ligands_pdb(ligs, lig)
            gid = int(args.gpu) if args.gpu not in (None, "", "None") else None
            sc = scoring.score_with_fallback(prot, lig, gid)
            row.update(cnn=sc.get("cnn"), cnnaff=sc.get("cnnaff"), aff=sc.get("aff"))
        except Exception as e:
            row["note"] = f"gnina_fail:{e}"
        for f in (prot, lig):
            if os.path.exists(f):
                os.remove(f)

    # boltz affinity (seed 평균) — 원자료를 outputs로 캡처
    frag = args.frag_dir or rep.split(os.sep + "seed_")[0]
    probs = []
    for aj in glob.glob(os.path.join(frag, "seed_*", "boltz_results_*", "predictions", cid,
                                     f"affinity_{cid}.json")):
        try:
            d = json.load(open(aj))
        except Exception:
            continue
        if isinstance(d.get("affinity_probability_binary"), (int, float)):
            probs.append(d["affinity_probability_binary"])
    row["boltz"] = round(float(np.mean(probs)), 6) if probs else ""

    # hac
    tsv = args.ligand_tsv or os.path.join(O, "ligand.tsv")
    if os.path.exists(tsv):
        smis = cio.read_ligand_tsv(tsv)
        if smis:
            row["hac"] = heavy_atoms(smis[0][2]) or ""

    _write(out05, row)


def _write(out05, row):
    cols = ["cid", "boltz", "cnn", "cnnaff", "aff", "hac",
            "n_pockets", "n_pockets_pass", "pocket_sizes", "dom_pocket_size", "dom_pocket_frac",
            "n_pose_clusters", "pose_sizes", "dom_pose_size", "dom_pose_frac", "n_poses",
            "pocket_residues", "rep_cif", "note"]
    p = os.path.join(out05, "binding_row.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerow(row)
    print(f"[2_binding] {row['cid']}: boltz={row['boltz']} cnn={row['cnn']} "
          f"cnnaff={row['cnnaff']} pocket={row['pocket_residues'][:40]} note={row['note']} -> {p}")


if __name__ == "__main__":
    main()

# ── 단독 실행 ──
#   python 2_binding.py --out-dir .../targets/L01/outputs/L010001 --gpu 0
# (run_stage1_frag.sh 끝에 붙이거나, 배치에서 fragment마다 호출)
