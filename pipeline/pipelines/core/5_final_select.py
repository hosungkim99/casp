#!/usr/bin/env python3
"""
final_select.py (Step 4) - 포켓별 리간드 클러스터 후보들에 에너지(gnina)를 적용해 최종 ≤5개 선정.
종합점수 = 합의(클러스터 크기) + Boltz2(iptm) + Vina(gnina affinity) + Gnina(CNNscore) 가중합.
  (회의 "Boltz2·Vina·Gnina 세 값 조합" 반영. gnina 1회가 vina affinity+CNN점수를 동시 산출)
task에 따라 가중치 분기:
  - P : pose 중심(합의·pose 품질↑)
  - PA: affinity(vina·CNNaffinity) 비중↑
클러스터 응집도(size=갯수, rmsd_mean/rmsd_std=평균·분산)도 요약에 병기(사람 판단용).
boltz2 env + singularity gnina. 출력: 05_final/ {selection_summary.csv, model_*.cif, SELECTION_RATIONALE.md}.
"""
import argparse, csv, os, shutil
import gemmi
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio
import common.scoring as scoring

# 종합점수 가중치. 4항: consensus(합의=클러스터 크기), iptm(Boltz2 신뢰도),
#   aff(Vina affinity, 낮을수록↑ → 정규화 시 부호반전), cnn(Gnina CNNscore, pose 품질 0~1↑).
#   각 task별 합=1.0.
WEIGHTS = {"P":  dict(consensus=0.40, iptm=0.25, aff=0.15, cnn=0.20),
           "PA": dict(consensus=0.30, iptm=0.15, aff=0.30, cnn=0.25)}


def _fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def gnina_score(cif, work, tag, reuse=None, gpu_id=None):
    """gnina 점수. reuse[cif] 에 이미 계산된 affinity 있으면 재사용(포켓검증과 동일 cif → gnina 생략)."""
    if reuse:
        r = reuse.get(cif)
        if r and r.get("aff") is not None:
            return r
    prot = os.path.join(work, f"{tag}_prot.pdb"); lig = os.path.join(work, f"{tag}_lig.pdb")
    _, ligs = cio.parse_complex(cif)
    if not ligs:
        return {}
    cio.write_ligands_pdb(ligs, lig)
    cio.write_receptor_pdb(cif, prot)
    return scoring.score_with_fallback(prot, lig, gpu_id)


def norm(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return lambda x: 0.5
    lo, hi = min(v), max(v)
    return (lambda x: 0.5) if hi == lo else (lambda x: (x - lo) / (hi - lo) if x is not None else 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand-clusters", required=True)
    ap.add_argument("--validation", required=True)
    ap.add_argument("--out", required=True, help="05_final 폴더")
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--protein-clusters", default="",
                    help="01 protein_clusters.csv (있고 형태 2개면 conf별 top개 선정 → CASP v1/v2)")
    ap.add_argument("--min-consensus", type=int, default=2,
                    help="이 크기 미만 클러스터(=합의 없는 singleton)는 그룹 하위로 강등(선정은 되되 상위엔 안 옴)")
    ap.add_argument("--task", default="P", choices=["P", "PA"])
    ap.add_argument("--posebusters", default="", help="04b posebusters.csv (있으면 무효 포즈 제외)")
    ap.add_argument("--jobs", type=int, default=1, help="병렬 워커 수(공용서버 배려 기본 1=순차)")
    ap.add_argument("--gpus", default="", help="라운드로빈 GPU id 예:'0,1,2'(빈값=미지정)")
    args = ap.parse_args()
    work = os.path.join(args.out, "score"); os.makedirs(work, exist_ok=True)

    cands = list(csv.DictReader(open(args.ligand_clusters)))
    val_rows = list(csv.DictReader(open(args.validation)))
    pval = {r["pocket_id"]: r for r in val_rows}
    # 포켓 검증(3)에서 포켓대표 cif 에 이미 계산한 gnina → 동일 cif 후보는 재계산 생략
    reuse = {r["rep_cif"]: {"aff": _fnum(r.get("gnina_affinity")),
                            "cnn": _fnum(r.get("cnn_score")),
                            "cnnaff": _fnum(r.get("cnn_affinity"))}
             for r in val_rows if r.get("rep_cif")}
    if not cands:
        raise SystemExit("ligand_clusters.csv 비어있음")

    # ① PoseBusters 무효 처리 + medoid 강건화
    #   대표(중심성 rank1)가 무효여도 클러스터를 통째로 버리지 않고, 그 클러스터의 valid 멤버 중
    #   가장 중심(rank 최소)인 것으로 대표를 교체(→ 큰 합의 클러스터가 대표 하나 때문에 사라지는 것 방지).
    #   전 멤버 무효인 클러스터만 제외. (posebusters.csv에 member_rank 없으면 기존 동작으로 폴백)
    pbvalid = {}
    if args.posebusters and os.path.exists(args.posebusters):
        by_cl = {}
        for r in csv.DictReader(open(args.posebusters)):
            key = (r["pocket_id"], r["ligand_cluster_id"])
            try:
                rank = int(r.get("member_rank", 1) or 1)
            except ValueError:
                rank = 1
            by_cl.setdefault(key, []).append((rank, r.get("rep_cif", ""), r.get("valid", "")))
        kept, nrescue, ndrop = [], 0, 0
        for c in cands:
            key = (c["pocket_id"], c["ligand_cluster_id"])
            mem = sorted(by_cl.get(key, []))
            valids = [(rk, cif) for rk, cif, v in mem if v == "True"]
            if not mem:
                pbvalid[key] = "NA"; kept.append(c)                 # 미검사 → 보존
            elif valids:
                rk, cif = valids[0]                                 # 가장 중심인 valid 멤버
                if rk != 1 and cif:
                    c["rep_cif"] = cif; nrescue += 1                # 대표 교체(강건화)
                pbvalid[key] = "True"; kept.append(c)
            elif any(v not in ("False",) for _, _, v in mem):
                pbvalid[key] = "NA"; kept.append(c)                 # NA/ERR만 → 보존(보수적)
            else:
                pbvalid[key] = "False"; ndrop += 1                  # 전 멤버 무효 → 제외
        if kept:
            msg = f"[final_select] PoseBusters 무효 {ndrop}개 제외"
            if nrescue:
                msg += f", 대표 교체(medoid 강건화) {nrescue}개"
            print(f"{msg} → {len(kept)} 후보")
            cands = kept
        else:
            print("[final_select] 경고: 전부 무효로 판정 → 필터 미적용(원본 유지)")

    # 각 후보 점수화 (opt-in 병렬; reuse 로 포켓검증 동일 cif 는 gnina 생략)
    gpus = scoring.parse_gpus(args.gpus)
    fns = [(lambda gpu, c=c: gnina_score(
                c["rep_cif"], work, f"P{c['pocket_id']}_L{c['ligand_cluster_id']}", reuse, gpu))
           for c in cands]
    scres = scoring.run_jobs(fns, args.jobs, gpus)
    for c, sc in zip(cands, scres):
        sc = sc or {}
        c["gnina_affinity"] = sc.get("aff"); c["cnn_score"] = sc.get("cnn"); c["cnn_affinity"] = sc.get("cnnaff")
        # Boltz2 신뢰도 = interface iptm(ligand_iptm) 우선, 없으면 global iptm 폴백.
        #   (ligand_iptm = chain_pair off-diagonal 평균 = 리간드-단백질 결합면 신뢰도 → pose에 더 직접적)
        try:
            c["_iptm"] = float(c.get("rep_ligand_iptm") or c.get("rep_iptm") or 0)
        except ValueError:
            c["_iptm"] = 0.0
        c["_size"] = int(c["size"])
        print(f"  P{c['pocket_id']}.L{c['ligand_cluster_id']} size {c['size']} "
              f"iptm {c['_iptm']:.3f} vina_aff {c['gnina_affinity']} cnn {c['cnn_score']}")

    n_size = norm([c["_size"] for c in cands])
    n_iptm = norm([c["_iptm"] for c in cands])
    # Vina affinity는 낮을수록 좋음 → 부호 반전. Gnina CNNscore는 높을수록 좋음(그대로).
    n_aff = norm([-c["gnina_affinity"] if c["gnina_affinity"] is not None else None for c in cands])
    n_cnn = norm([c["cnn_score"] for c in cands])
    w = WEIGHTS[args.task]
    for c in cands:
        aff_v = -c["gnina_affinity"] if c["gnina_affinity"] is not None else None
        c["composite"] = round(
            w["consensus"] * n_size(c["_size"]) + w["iptm"] * n_iptm(c["_iptm"]) +
            w["aff"] * n_aff(aff_v) + w["cnn"] * n_cnn(c["cnn_score"]), 4)

    # 정렬: composite 내림차순, 동점 시 template_dist 작은(템플릿에 가까운) 포켓 우선(자동 tie-break).
    def _tdist(c):
        try:
            return float(pval.get(c["pocket_id"], {}).get("template_dist", "NA"))
        except (TypeError, ValueError):
            return float("inf")   # 템플릿 없음/NA → 동점에서 뒤로
    # ── conf(단백질 형태)별 선택: CASP v1/v2 대응 ──
    # protein_clusters.csv로 각 후보(rep_cif)의 형태를 매핑. 형태가 2개면 각 형태에서 top개씩
    # 뽑아 "우세형태 → model 1..top", "다음형태 → model top+1..2*top" 로 번호부여.
    # 형태정보 없거나 1개면 기존 동작(전체 top개 = model 1..top). → downstream 09가 CASP번호 매핑.
    confmap = {}
    if args.protein_clusters and os.path.exists(args.protein_clusters):
        for pr in csv.DictReader(open(args.protein_clusters)):
            if pr.get("cif"):
                confmap[pr["cif"]] = pr.get("conf_id", "")
    for c in cands:
        c["_conf"] = confmap.get(c["rep_cif"], "")
    confs = [cf for cf in set(c["_conf"] for c in cands) if cf != ""]
    if len(confs) >= 2:
        # 형태별 총 size 합으로 정렬 → 우세형태 먼저(models 1..top). CASP=2형태만 채택.
        conf_order = sorted(confs, key=lambda cf: sum(int(c["size"]) for c in cands
                                                      if c["_conf"] == cf), reverse=True)[:2]
        groups = [(cf, [c for c in cands if c["_conf"] == cf]) for cf in conf_order]
        print(f"[final_select] conf-aware 선택: 형태 {conf_order} "
              f"(형태별 후보수 {[(cf, len(g)) for cf, g in groups]}) → 각 형태 top {args.top}")
    else:
        groups = [("", cands)]                       # 단일형태(기존 동작)

    rows = []
    idx = 0
    for cf, grp in groups:
        # singleton/tiny(합의 없는) 클러스터는 그룹 하위로 강등 → populated pose가 상위 순번을 차지.
        # (size >= min_consensus 인 것 우선, 그 안에서 composite·template 근접 순.)
        grp.sort(key=lambda c: (int(c["size"]) >= args.min_consensus, c["composite"], -_tdist(c)),
                 reverse=True)
        for c in grp[:args.top]:
            idx += 1
            dst = os.path.join(args.out, f"model_{idx}.cif")
            shutil.copy(c["rep_cif"], dst)
            rows.append({"model": f"model_{idx}", "conf_id": cf, "pocket_id": c["pocket_id"],
                         "ligand_cluster_id": c["ligand_cluster_id"], "size": c["size"],
                         "rmsd_mean": c.get("rmsd_mean", ""), "rmsd_std": c.get("rmsd_std", ""),
                         "iptm": c.get("rep_iptm", ""), "ligand_iptm": c.get("rep_ligand_iptm", ""),
                         "confidence": c.get("rep_confidence", ""), "lscore": c.get("rep_iptm", ""),
                         "gnina_affinity": c["gnina_affinity"], "cnn_score": c["cnn_score"],
                         "cnn_affinity": c.get("cnn_affinity", ""),
                         "composite": c["composite"], "pocket_pass": pval.get(c["pocket_id"], {}).get("pass", ""),
                         "posebusters_valid": pbvalid.get((c["pocket_id"], c["ligand_cluster_id"]), "NA"),
                         "source_cif": c["rep_cif"]})

    cols = ["model", "conf_id", "pocket_id", "ligand_cluster_id", "size", "rmsd_mean", "rmsd_std",
            "iptm", "ligand_iptm", "confidence", "lscore", "gnina_affinity", "cnn_score",
            "cnn_affinity", "composite", "pocket_pass", "posebusters_valid", "source_cif"]
    with open(os.path.join(args.out, "selection_summary.csv"), "w", newline="") as f:
        w2 = csv.DictWriter(f, fieldnames=cols); w2.writeheader(); w2.writerows(rows)

    md = [f"# 최종 선정 근거 (Task={args.task})\n",
          f"종합점수 가중치: consensus {w['consensus']}, iptm(Boltz2) {w['iptm']}, "
          f"aff(Vina) {w['aff']}, cnn(Gnina) {w['cnn']}\n",
          "(P=pose 중심, PA=affinity 비중↑ / size=클러스터 갯수, rmsd_mean·std=방향 응집도)\n\n",
          "| model | conf | pocket | lig_cluster | size | rmsd_mean±std | iptm | vina_aff | cnn | composite |\n",
          "|---|---|---|---|---|---|---|---|---|---|\n"]
    for r in rows:
        md.append(f"| {r['model']} | {r.get('conf_id','') or '-'} | P{r['pocket_id']} | L{r['ligand_cluster_id']} | {r['size']} "
                  f"| {r['rmsd_mean']}±{r['rmsd_std']} | {r['iptm']} | {r['gnina_affinity']} "
                  f"| {r['cnn_score']} | {r['composite']} |\n")
    md.append("\n선정 원칙: 포켓 검증 통과 + 합의(클러스터 크기) + Boltz2(iptm) + Vina(affinity) + Gnina(CNN) 종합.\n")
    md.append("응집도 참고: size 클수록/rmsd_mean·std 작을수록 방향 합의가 강함(신뢰↑).\n")
    md.append("주의: 같은 포켓에서 여러 개가 뽑히면 방향 다양성 확인, 다른 포켓이 섞이면 hedge.\n")
    open(os.path.join(args.out, "SELECTION_RATIONALE.md"), "w").write("".join(md))

    print(f"\n[final_select] 최종 {len(rows)}개 -> {args.out}/ (selection_summary.csv, model_*.cif, RATIONALE.md)")
    for r in rows:
        cf = f"conf{r['conf_id']} " if r.get("conf_id") else ""
        print(f"  {r['model']}: {cf}P{r['pocket_id']}.L{r['ligand_cluster_id']} composite {r['composite']}")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (먼저: source $CASP17/scripts/env_setup.sh)
#   SC=$CASP17/users/USERNAME/scripts ; OUT=$CASP17/users/USERNAME/targets/T2383
#   micromamba run -n boltz2 python $SC/5_final_select.py \
#       --ligand-clusters $OUT/04_ligand_clusters/ligand_clusters.csv \
#       --validation $OUT/03_pocket_validation/pocket_validation.csv \
#       --out $OUT/05_final --top 5 --task P        # PA면 --task PA
# (평소엔 run_pipeline.py 가 순서대로 자동 호출)
