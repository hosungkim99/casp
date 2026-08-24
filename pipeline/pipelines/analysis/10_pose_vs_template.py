#!/usr/bin/env python3
"""
10_pose_vs_template.py - 최종 후보 포즈를 실험 템플릿(리간드 포함)과 "후보별"로 비교.

기존 contact_fingerprint는 단백질 Cα medoid 하나에 앵커 → 흩어진/다중포켓 타겟에서
대표성 없음(예: T2412는 L1이 아니라 L3가 실험과 더 맞음). 이 스크립트는 단일 앵커 대신
최종 model_1~5 각각을 실험 템플릿 리간드와 개별 대조한다.

각 (model, 템플릿)마다:
  - 단백질 CA 정렬(gemmi, 서열기반 → 번호 offset 무관)
  - ligand overlap% : 우리 리간드 원자 중 실험 리간드 4Å 이내 비율 (분모=우리 원자 → '잡아먹기' 함정 회피)
  - contact 잔기 겹침(템플릿 번호) : 실험 리간드 접촉잔기 ∩ 우리 리간드 접촉잔기, Jaccard
  - (선택) ligand_clusters.csv 주면 그 후보의 클러스터 크기/수렴도 병기
다중 리간드 копи(FN9 2벌 등)·이온·펩타이드는 자동 처리(체인 근접 리간드만 사용).
boltz2 env (gemmi+numpy). 출력: pose_vs_template.csv + 콘솔 랭킹.
"""
import argparse, csv, glob, os
import numpy as np
import gemmi


def polymer_chains(model, minlen=40):
    """아미노산 폴리머 체인만."""
    out = []
    for ch in model:
        pol = ch.get_polymer()
        if pol and len(pol) >= minlen:
            out.append(ch)
    return out


def heavy_ligand_atoms(model, near_chain=None, near_cutoff=8.0):
    """비-단백질/비-핵산/비-물, 탄소 포함(유기) 잔기의 heavy atom 좌표.
       near_chain 주면 그 체인 CA 근처(리간드 копи 선택)만."""
    ca = None
    if near_chain is not None:
        ca = np.array([[a.pos.x, a.pos.y, a.pos.z] for r in near_chain
                       for a in r if a.name == "CA"])
    coords = []
    for ch in model:
        for r in ch:
            t = gemmi.find_tabulated_residue(r.name)
            if (t and (t.is_amino_acid() or t.is_nucleic_acid())) or r.name in ("HOH", "WAT"):
                continue
            elems = [a.element.name for a in r if a.element.name != "H"]
            if "C" not in elems:            # 이온 등(탄소 없음) 제외
                continue
            xs = [[a.pos.x, a.pos.y, a.pos.z] for a in r if a.element.name != "H"]
            if ca is not None and len(xs):
                d = np.sqrt(((ca[:, None, :] - np.array(xs)[None, :, :]) ** 2).sum(-1))
                if d.min() > near_cutoff:   # 이 체인에서 먼 копи 제외
                    continue
            coords += xs
    return np.array(coords) if coords else np.zeros((0, 3))


def protein_heavy(chain):
    """체인의 아미노산 heavy atom 좌표 + 잔기번호 라벨."""
    coords, labels = [], []
    for r in chain:
        t = gemmi.find_tabulated_residue(r.name)
        if not (t and t.is_amino_acid()):
            continue
        for a in r:
            if a.element.name == "H":
                continue
            coords.append([a.pos.x, a.pos.y, a.pos.z]); labels.append(r.seqid.num)
    return np.array(coords), np.array(labels)


def contact_residues(prot_xyz, prot_lab, query_xyz, cutoff):
    """query(리간드) cutoff 이내 단백질 잔기번호 집합."""
    if len(query_xyz) == 0 or len(prot_xyz) == 0:
        return set()
    d = np.sqrt(((prot_xyz[:, None, :] - query_xyz[None, :, :]) ** 2).sum(-1))
    return set(int(x) for x in prot_lab[(d <= cutoff).any(1)])


def apply_T(T, xyz):
    out = []
    for p in xyz:
        v = T.mat.multiply(gemmi.Vec3(*p)) + T.vec
        out.append([v.x, v.y, v.z])
    return np.array(out) if out else np.zeros((0, 3))


def best_superpose(fixed_chains, moving_chain):
    """moving을 각 fixed 체인에 정렬 → RMSD 최소인 것 선택 (이량체 대응)."""
    best = None
    for fx in fixed_chains:
        try:
            sup = gemmi.calculate_superposition(fx.get_polymer(), moving_chain.get_polymer(),
                                                gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
        except Exception:
            continue
        if best is None or sup.rmsd < best[0].rmsd:
            best = (sup, fx)
    return best  # (sup, fixed_chain) or None


def load_convergence(path):
    """ligand_clusters.csv → {(pocket,cluster): size}, 포켓별 총합."""
    sizes, tot = {}, {}
    if path and os.path.exists(path):
        for r in csv.DictReader(open(path)):
            k = (r["pocket_id"], r["ligand_cluster_id"]); s = int(r["size"])
            sizes[k] = s; tot[r["pocket_id"]] = tot.get(r["pocket_id"], 0) + s
    return sizes, tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--final-dir", required=True, help="05_final 또는 05b_refined (model_*.cif)")
    ap.add_argument("--templates", required=True, help="쉼표구분 템플릿 cif 경로(리간드 포함 RCSB 원본)")
    ap.add_argument("--out", required=True, help="pose_vs_template.csv")
    ap.add_argument("--ligand-clusters", default="", help="04 ligand_clusters.csv (수렴도 병기, 선택)")
    ap.add_argument("--summary", default="", help="05(b) selection_summary.csv (pocket/cluster 매핑, 선택)")
    ap.add_argument("--cutoff", type=float, default=4.0, help="접촉 거리")
    args = ap.parse_args()

    tmpls = [t.strip() for t in args.templates.split(",") if t.strip()]
    sizes, tot = load_convergence(args.ligand_clusters)
    # model -> (pocket,cluster)
    m2pc = {}
    if args.summary and os.path.exists(args.summary):
        for r in csv.DictReader(open(args.summary)):
            m2pc[r["model"]] = (r.get("pocket_id", ""), r.get("ligand_cluster_id", ""))

    rows = []
    for cif in sorted(glob.glob(os.path.join(args.final_dir, "model_*.cif"))):
        name = os.path.splitext(os.path.basename(cif))[0]
        st = gemmi.read_structure(cif); st.setup_entities(); mm = st[0]
        mpol = polymer_chains(mm)
        if not mpol:
            continue
        our_lig = heavy_ligand_atoms(mm)
        conv = ""
        if name in m2pc and m2pc[name] in sizes:
            pid, cid = m2pc[name]
            conv = f"{sizes[(pid,cid)]}/{tot.get(pid,'?')}"
        for tp in tmpls:
            tn = os.path.splitext(os.path.basename(tp))[0]
            ts = gemmi.read_structure(tp); ts.setup_entities(); tm = ts[0]
            tpol = polymer_chains(tm)
            if not tpol:
                continue
            bs = best_superpose(tpol, mpol[0])
            if bs is None:
                continue
            sup, fx = bs
            pc, pl = protein_heavy(fx)                       # 템플릿 fixed 체인 단백질
            exp_lig = heavy_ligand_atoms(tm, near_chain=fx)  # fixed 체인 근처 실험 리간드
            if len(exp_lig) == 0:
                rows.append(dict(model=name, template=tn, exp_ligand="none/apo",
                                 sup_rmsd=round(sup.rmsd, 2), overlap_pct="", shared="", jaccard=""))
                continue
            our_T = apply_T(sup.transform, our_lig)           # 우리 리간드를 템플릿 좌표계로
            # ligand overlap% (분모=우리 원자)
            if len(our_T):
                d = np.sqrt(((our_T[:, None, :] - exp_lig[None, :, :]) ** 2).sum(-1))
                ov = 100.0 * (d.min(1) <= args.cutoff).mean()
            else:
                ov = 0.0
            setA = contact_residues(pc, pl, exp_lig, args.cutoff)   # 실험 리간드 접촉잔기
            setB = contact_residues(pc, pl, our_T, args.cutoff)     # 우리 리간드 접촉잔기
            shared = sorted(setA & setB)
            jac = len(shared) / len(setA | setB) if (setA | setB) else 0.0
            rows.append(dict(model=name, template=tn,
                             exp_ligand=f"{len(setA)}res", sup_rmsd=round(sup.rmsd, 2),
                             overlap_pct=round(ov, 0), shared=f"{len(shared)}",
                             shared_res=";".join(str(x) for x in shared),
                             jaccard=round(jac, 2), convergence=conv))

    cols = ["model", "template", "sup_rmsd", "overlap_pct", "shared", "jaccard",
            "convergence", "exp_ligand", "shared_res"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)

    print(f"[pose_vs_template] {len(rows)}행 -> {args.out}\n")
    print(f"{'model':9} {'template':7} {'RMSD':5} {'overlap%':8} {'shared':6} {'jaccard':7} {'conv':8}")
    for r in rows:
        print(f"{r['model']:9} {r['template']:7} {str(r['sup_rmsd']):5} "
              f"{str(r.get('overlap_pct','')):8} {str(r.get('shared','')):6} "
              f"{str(r.get('jaccard','')):7} {str(r.get('convergence','')):8}")
    # 후보별 최고 overlap 요약
    bym = {}
    for r in rows:
        if r.get("overlap_pct") == "":
            continue
        bym.setdefault(r["model"], []).append((r["overlap_pct"], r["jaccard"], r["template"]))
    if bym:
        print("\n[후보 랭킹] 최고 overlap 기준 (실험과 가장 일치하는 후보 = 제출 상위 후보):")
        rank = sorted(bym.items(), key=lambda kv: -max(x[0] for x in kv[1]))
        for m, v in rank:
            best = max(v)
            print(f"  {m}: overlap {best[0]}% / jaccard {best[1]} (vs {best[2]})")


if __name__ == "__main__":
    main()

# ── 단독 실행 (boltz2 env; gemmi) ──
#   T=/path/to/casp17-ligand/users/USERNAME/templates
#   OUT=.../users/USERNAME/targets/T2412/outputs
#   micromamba run -n boltz2 python 10_pose_vs_template.py \
#       --final-dir $OUT/05_final \
#       --templates $T/7c9n.cif,$T/6bhd.cif,$T/3dlm.cif \
#       --summary $OUT/05_final/selection_summary.csv \
#       --ligand-clusters $OUT/04_ligand_clusters/ligand_clusters.csv \
#       --out $OUT/06_validation/pose_vs_template.csv
# (템플릿은 리간드 포함 RCSB 원본이어야 함. inputs의 template_A_*.cif는 리간드 제거됨 → 사용 불가)
