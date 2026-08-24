#!/usr/bin/env python3
"""
pocket_validate.py (Step 2) - 포켓 후보 검증 및 압축.
각 포켓 후보 대표에 p2rank(캐비티 존재?) + gnina(에너지)를 적용해 비교·검증하고,
기준 통과(pass) 포켓만 남겨 압축. boltz2 env + p2rank + singularity gnina.
출력: 03_pocket_validation/ {pocket_validation.csv, p2rank/, rescore/}.

autodock은 포켓별 박스 도킹이 필요해 기본 미포함(옵션 훅). gnina affinity가 1차 에너지.
"""
import argparse, csv, os, subprocess, statistics
import numpy as np
import gemmi
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio
import common.scoring as scoring

PRANK = os.environ.get("PRANK", "/path/to/casp17-ligand/models/p2rank/p2rank_2.5.1/prank")


def _fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def protein_pdb_centroid(cif, pdb_out):
    _, ligs = cio.parse_complex(cif)
    cen = cio.ligand_centroid(ligs)
    cio.write_receptor_pdb(cif, pdb_out)
    return np.array(cen) if cen else None


def score_cif(cif, prot, lig, gpu_id):
    """멤버 cif 하나를 gnina 채점 → (aff, cnn). 리간드 없으면 (None, None)."""
    _, ligs = cio.parse_complex(cif)
    if not ligs:
        return None, None
    cio.write_ligands_pdb(ligs, lig)
    cio.write_receptor_pdb(cif, prot)
    sc = scoring.score_with_fallback(prot, lig, gpu_id)
    return sc.get("aff"), sc.get("cnn")


def run_prank(pdb, outdir):
    os.makedirs(outdir, exist_ok=True)
    subprocess.run([PRANK, "predict", "-f", pdb, "-o", outdir],
                   check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    for root, _, files in os.walk(outdir):
        for f in files:
            if f.endswith("predictions.csv"):
                return os.path.join(root, f)
    return None


def nearest_pocket(pcsv, cen):
    best = None
    with open(pcsv) as fh:
        rd = csv.reader(fh); hdr = [h.strip() for h in next(rd)]
        idx = {h: i for i, h in enumerate(hdr)}
        gi = lambda *ns: next((idx[n] for n in ns if n in idx), None)
        cx, cy, cz, rk = gi("center_x"), gi("center_y"), gi("center_z"), gi("rank")
        for k, row in enumerate(rd):
            if cx is None or len(row) <= cx:
                continue
            try:
                c = np.array([float(row[cx]), float(row[cy]), float(row[cz])])
            except ValueError:
                continue
            d = float(np.linalg.norm(cen - c))
            if best is None or d < best[0]:
                best = (d, row[rk].strip() if rk is not None else str(k + 1), k == 0)
    return best  # (dist, rank, is_top)


def _poly_chains(model, minlen=40):
    return [ch for ch in model if ch.get_polymer() and len(ch.get_polymer()) >= minlen]


def template_ligand_centroids(templates, ref_cif):
    """각 템플릿을 reference 구조에 Cα(서열기반) 정렬 → 리간드 무게중심을 reference 좌표계로.
    pocket center_x/y/z 와 같은 프레임의 centroid 목록 반환(포켓별 template_dist 계산용).
    회의 1: 'template 찾기'로 포켓을 결정 — 리간드 있는 템플릿이 가리키는 자리 = 진짜 포켓 후보."""
    cents = []
    if not templates or not ref_cif or not os.path.exists(ref_cif):
        return cents
    try:
        rs = gemmi.read_structure(ref_cif); rs.setup_entities(); rm = rs[0]
    except Exception:
        return cents
    ref_pol = _poly_chains(rm)
    if not ref_pol:
        return cents
    for tp in templates:
        if not tp or not os.path.exists(tp):
            continue
        try:
            ts = gemmi.read_structure(tp); ts.setup_entities(); tm = ts[0]
        except Exception:
            continue
        tpol = _poly_chains(tm)
        if not tpol:
            continue
        # 템플릿 폴리머를 ref 폴리머들에 정렬 → RMSD 최소 변환 채택(이량체/멀티체인 대응)
        best = None
        for rf in ref_pol:
            for tc in tpol:
                try:
                    sup = gemmi.calculate_superposition(rf.get_polymer(), tc.get_polymer(),
                            gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
                except Exception:
                    continue
                if best is None or sup.rmsd < best.rmsd:
                    best = sup
        if best is None:
            continue
        for ch in tm:
            for r in ch:
                t = gemmi.find_tabulated_residue(r.name)
                if (t and (t.is_amino_acid() or t.is_nucleic_acid())) or r.name in cio.WATER:
                    continue
                xs = [[a.pos.x, a.pos.y, a.pos.z] for a in r if a.element.name != "H"]
                if not xs or not any(a.element.name == "C" for a in r):  # 이온(탄소없음) 제외
                    continue
                arr = []
                for p in xs:
                    v = best.transform.mat.multiply(gemmi.Vec3(*p)) + best.transform.vec
                    arr.append([v.x, v.y, v.z])
                cents.append(np.array(arr).mean(0))
    return cents


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", required=True, help="02의 pocket_candidates.csv")
    ap.add_argument("--members", default="", help="02의 members.csv(주면 포켓별 멤버 평균/분산 산출)")
    ap.add_argument("--members-topk", type=int, default=5,
                    help="포켓별 채점할 멤버 수(iptm순, 대표 포함). 0=대표만(멤버평균 생략)")
    ap.add_argument("--out", required=True, help="03_pocket_validation 폴더")
    ap.add_argument("--gnina-cutoff", type=float, default=-4.0, help="pass: gnina affinity <= 이 값")
    ap.add_argument("--p2rank-dist", type=float, default=6.0, help="pass: 최근접 p2rank 포켓 거리 <= 이 값")
    ap.add_argument("--jobs", type=int, default=1, help="병렬 워커 수(공용서버 배려 기본 1=순차)")
    ap.add_argument("--gpus", default="", help="라운드로빈 GPU id 예:'0,1,2'(빈값=미지정)")
    ap.add_argument("--templates", default="",
                    help="쉼표구분 템플릿 cif(리간드 포함 RCSB). 주면 포켓별 template_dist 병기")
    ap.add_argument("--reference", default="",
                    help="01 reference.txt(또는 cif 경로). 템플릿을 이 좌표계로 정렬(center와 동일 프레임)")
    ap.add_argument("--template-pass-dist", type=float, default=8.0,
                    help="template_dist<=이 값이면 template-supported→pass 승격(리간드 있는 템플릿=강한 증거)")
    ap.add_argument("--auto-relax", action="store_true",
                    help="통과 0개면 기준을 단계적으로 완화해 재판정(추가 gnina/p2rank 없음, 임계값만 재적용)")
    ap.add_argument("--relax-steps", type=int, default=3, help="완화 최대 횟수")
    args = ap.parse_args()
    work = args.out
    os.makedirs(os.path.join(work, "p2rank"), exist_ok=True)
    os.makedirs(os.path.join(work, "rescore"), exist_ok=True)

    cands = list(csv.DictReader(open(args.candidates)))

    # 템플릿 리간드 centroid(reference 좌표계) — 포켓별 template_dist 산출용.
    ref_cif = args.reference
    if ref_cif and ref_cif.endswith(".txt") and os.path.exists(ref_cif):
        ref_cif = open(ref_cif).read().strip().splitlines()[0].strip() if os.path.getsize(ref_cif) else ""
    tmpls = [t.strip() for t in args.templates.split(",") if t.strip()]
    tcents = template_ligand_centroids(tmpls, ref_cif)
    if tmpls:
        print(f"[pocket_validate] 템플릿 리간드 centroid {len(tcents)}개 "
              f"(정렬기준 {os.path.basename(ref_cif) if ref_cif else 'NA'})")

    def _template_dist(c):
        """포켓 center(reference 프레임)에서 가장 가까운 템플릿 리간드까지 거리(Å). 없으면 'NA'."""
        if not tcents:
            return "NA"
        try:
            cen = np.array([float(c["center_x"]), float(c["center_y"]), float(c["center_z"])])
        except (KeyError, ValueError):
            return "NA"
        return round(float(min(np.linalg.norm(cen - tc) for tc in tcents)), 2)

    # 포켓별 멤버(iptm 내림차순) — mean/std score 산출용. members.csv 없으면 대표만 채점.
    members_by_pocket = {}
    if args.members and os.path.exists(args.members):
        for m in csv.DictReader(open(args.members)):
            members_by_pocket.setdefault(m["pocket_id"], []).append(m)
        for pid in members_by_pocket:
            members_by_pocket[pid].sort(
                key=lambda m: float(m["iptm"]) if str(m.get("iptm", "")).strip() not in ("", "NA") else -1.0,
                reverse=True)

    def evaluate(c, gpu_id):
        """포켓 대표에 p2rank+gnina(pass/fail) + 멤버 top-K 채점(mean/std). 병렬 워커 단위."""
        pid, cif = c["pocket_id"], c["rep_cif"]
        prot = os.path.join(work, "rescore", f"P{pid}_protein.pdb")
        lig = os.path.join(work, "rescore", f"P{pid}_ligand.pdb")
        cen = protein_pdb_centroid(cif, prot)
        if cen is None:
            print(f"=== pocket {pid}: 리간드 없음"); return None
        _, ligs = cio.parse_complex(cif)
        cio.write_ligands_pdb(ligs, lig)
        pcsv = run_prank(prot, os.path.join(work, "p2rank", f"P{pid}"))
        np_ = nearest_pocket(pcsv, cen) if pcsv else None
        sc = scoring.score_with_fallback(prot, lig, gpu_id)
        aff = sc.get("aff"); p2d = np_[0] if np_ else None
        # pass 판정 = 물리(대표 gnina+p2rank) OR 템플릿 지지(리간드 있는 템플릿이 이 자리를 가리킴).
        phys_pass = (aff is not None and aff <= args.gnina_cutoff) and \
                    (p2d is not None and p2d <= args.p2rank_dist)
        tdist = _template_dist(c)
        tsupport = isinstance(tdist, (int, float)) and tdist <= args.template_pass_dist
        passed = phys_pass or tsupport
        # ── 멤버 평균/분산: 대표(=members 중 하나) 점수를 seed로, 나머지 top-K를 추가 채점 ──
        m_aff = [aff] if aff is not None else []
        m_cnn = [sc.get("cnn")] if sc.get("cnn") is not None else []
        n_scored = 1 if aff is not None else 0
        for j, m in enumerate(members_by_pocket.get(str(pid), [])[:max(args.members_topk, 0)]):
            mcif = m.get("cif", "")
            if not mcif or mcif == cif:            # 대표는 위에서 이미 채점
                continue
            ma, mc = score_cif(mcif, os.path.join(work, "rescore", f"P{pid}_m{j}_prot.pdb"),
                               os.path.join(work, "rescore", f"P{pid}_m{j}_lig.pdb"), gpu_id)
            if ma is not None:
                m_aff.append(ma); n_scored += 1
            if mc is not None:
                m_cnn.append(mc)
        mean_aff = round(statistics.mean(m_aff), 3) if m_aff else "NA"
        std_aff = round(statistics.pstdev(m_aff), 3) if len(m_aff) > 1 else 0.0
        mean_cnn = round(statistics.mean(m_cnn), 3) if m_cnn else "NA"
        treason = " [template-supported]" if (tsupport and not phys_pass) else ""
        print(f"=== pocket {pid} (size {c['size']}, {c['models']}): "
              f"p2rank dist={p2d} vina_aff={aff} cnn={sc.get('cnn')} pass={passed}{treason} "
              f"| mean_aff={mean_aff}±{std_aff} (n={n_scored}) tdist={tdist}")
        return {
            "pocket_id": pid, "size": c["size"], "n_models": c["n_models"],
            "models": c["models"], "conf_composition": c.get("conf_composition", ""),
            "p2rank_dist": round(p2d, 2) if p2d is not None else "NA",
            "p2rank_rank": np_[1] if np_ else "NA",
            "gnina_affinity": aff, "cnn_score": sc.get("cnn"),
            "cnn_affinity": sc.get("cnnaff"),
            "mean_affinity": mean_aff, "std_affinity": std_aff,
            "mean_cnn_score": mean_cnn, "n_scored": n_scored,
            "template_dist": tdist, "template_support": (tsupport if tcents else "NA"),
            "phys_pass": phys_pass, "rep_cif": cif, "pass": passed}

    gpus = scoring.parse_gpus(args.gpus)
    fns = [(lambda gpu, c=c: evaluate(c, gpu)) for c in cands]
    results = [r for r in scoring.run_jobs(fns, args.jobs, gpus) if r is not None]

    npass = sum(1 for r in results if r["pass"])
    # 자동 완화: 통과 0개면 이미 계산된 aff/거리로 임계값만 단계적으로 완화해 재판정(재계산 없음)
    if npass == 0 and args.auto_relax and results:
        cutoff, dist = args.gnina_cutoff, args.p2rank_dist
        for _ in range(args.relax_steps):
            cutoff += 1.0; dist += 2.0
            for r in results:
                aff, p2d = _fnum(r["gnina_affinity"]), _fnum(r["p2rank_dist"])
                r["pass"] = (aff is not None and aff <= cutoff) and (p2d is not None and p2d <= dist)
            npass = sum(1 for r in results if r["pass"])
            if npass > 0:
                print(f"[auto-relax] 기준 완화 gnina<= {cutoff}, p2rank<= {dist} → {npass} 포켓 통과")
                break

    cols = ["pocket_id", "size", "n_models", "models", "conf_composition",
            "p2rank_dist", "p2rank_rank", "gnina_affinity", "cnn_score",
            "cnn_affinity", "mean_affinity", "std_affinity", "mean_cnn_score", "n_scored",
            "template_dist", "template_support", "phys_pass", "pass", "rep_cif"]
    with open(os.path.join(work, "pocket_validation.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(results)
    print(f"\n[pocket_validate] {npass}/{len(results)} 포켓 통과 -> {work}/pocket_validation.csv")
    if npass == 0:
        print("  [주의] 통과 0 — 기준(--gnina-cutoff/--p2rank-dist) 완화 또는 데이터 재검토 필요")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (먼저: source $CASP17/scripts/env_setup.sh)
#   SC=$CASP17/users/USERNAME/scripts ; OUT=$CASP17/users/USERNAME/targets/T2383
#   micromamba run -n boltz2 python $SC/3_pocket_validate.py \
#       --candidates $OUT/02_pocket_candidates/pocket_candidates.csv --out $OUT/03_pocket_validation \
#       --gnina-cutoff -4.0 --p2rank-dist 6.0
# (평소엔 run_pipeline.py 가 순서대로 자동 호출)
