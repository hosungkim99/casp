#!/usr/bin/env python3
"""
pocket_candidates.py (Step 1) - 리간드 센트로이드 클러스터링 → 포켓 후보 (상위 N개).
전체 모델의 예측 리간드 위치(무게중심)를 단백질 정렬 좌표계에서 coarse 클러스터링.
작은 클러스터도 유효 사이트일 수 있으므로 top N(기본 10) 유지. boltz2 env (numpy+gemmi).
출력: 02_pocket_candidates/ {pocket_candidates.csv, members.csv}.
"""
import argparse, csv, os
import numpy as np
import gemmi
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio
import common.geom as geom

def greedy_cluster(items, threshold):
    """(row, centroid) 목록을 센트로이드 거리 임계로 그리디 클러스터. [{center, members}] 반환."""
    clusters = []
    for r, c in items:
        for cl in clusters:
            if np.linalg.norm(c - cl["center"]) <= threshold:
                cl["members"].append((r, c)); break
        else:
            clusters.append({"center": c, "members": [(r, c)]})
    for cl in clusters:
        cs = np.array([c for _, c in cl["members"]])
        cl["center"] = cs.mean(0)
    return clusters

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True, help="02_pocket_candidates 폴더")
    ap.add_argument("--protein-clusters", default="", help="01의 protein_clusters.csv(선택; conf 표기)")
    ap.add_argument("--reference", default="", help="01의 reference.txt(PC중심 구조; 없으면 rank1 폴백)")
    ap.add_argument("--threshold", type=float, default=8.0, help="포켓(센트로이드) 클러스터 임계 Å")
    ap.add_argument("--sweep", default="",
                    help="여러 cutoff 진단(콤마구분 예:'6,8,10,12'). 포켓 불명확 시 threshold 조정 참고용 "
                         "(정상 출력은 --threshold 값으로 그대로 생성)")
    ap.add_argument("--topn", type=int, default=10)
    ap.add_argument("--split-by-conf", action="store_true",
                    help="형태(conf)별로 따로 포켓 클러스터링(모양 오염 제거; 기본 off)")
    ap.add_argument("--cache", default="", help="0b geom_cache.pkl(있으면 재파싱 없이 사용)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    conf = {}
    if args.protein_clusters and os.path.exists(args.protein_clusters):
        for r in csv.DictReader(open(args.protein_clusters)):
            conf[r["cif"]] = r.get("conf_id", "")

    rows = [r for r in csv.DictReader(open(args.table)) if r.get("cif")]
    rows.sort(key=lambda r: int(r["rank"]))
    cache = cio.load_geom_cache(args.cache)
    ref_ca, ref_elems, _ = cio.cached_geometry(cio.reference_cif(rows, args.reference), cache)

    data = []  # (row, centroid[np3]) — 다중 copy면 구조당 copy별로 여러 점
    for i, r in enumerate(rows):
        ca, elems, coords = cio.cached_geometry(r["cif"], cache)
        if not ca or elems != ref_elems:
            continue
        rt = geom.align_to_ref(ref_ca, ca)
        if rt is None:
            continue
        # copy별 무게중심을 각각 정렬 → 구조당 copy 수만큼 점 추가(각 copy 자리 보존).
        #   (단일 copy면 1점 = 기존과 동일. 강체정렬은 centroid와 교환가능.)
        cens = cio.cached_copy_centroids(r["cif"], cache)
        if cens:
            for ac in geom.apply_rt(rt[0], rt[1], cens):
                data.append((r, ac))
        else:                                     # 폴백: concat 평균 1점(리간드 있으나 lc 없음)
            data.append((r, geom.apply_rt(rt[0], rt[1], coords).mean(0)))
        if (i + 1) % 500 == 0:
            print(f"  parsed {i+1}/{len(rows)}")
    print(f"[pocket] {len(data)} usable structures")

    # ── cutoff 스윕 진단 (옵션) ──
    # 회의 1-1: 포켓이 불명확할 때(큰 클러스터 없음) threshold 를 "조정해본다".
    # 여러 cutoff 로 클러스터링해 (클러스터 수, 최대 클러스터 크기·비율, 상위3 크기)를 표로 보여줘
    # 사람이 적절한 threshold 를 고르게 돕는다. 정상 출력에는 영향 없음(진단만).
    if args.sweep.strip():
        try:
            ts = [float(x) for x in args.sweep.split(",") if x.strip()]
        except ValueError:
            ts = []
        if ts:
            n = len(data)
            print(f"[sweep] cutoff 진단 (n={n} 구조)")
            print(f"  {'cutoff(Å)':>10} {'#clusters':>10} {'top1':>7} {'top1%':>7} {'top3_sizes':>16}")
            for t in ts:
                cls = greedy_cluster(data, t)
                cls.sort(key=lambda cl: len(cl["members"]), reverse=True)
                sizes = [len(cl["members"]) for cl in cls]
                top1 = sizes[0] if sizes else 0
                top3 = ",".join(str(s) for s in sizes[:3])
                print(f"  {t:>10.1f} {len(cls):>10} {top1:>7} {top1/n*100 if n else 0:>6.0f}% {top3:>16}")
            print("  [해석] top1%가 크고 클러스터 수가 적을수록 포켓이 뚜렷(합의 강함). "
                  "값이 완만하면 threshold를 키워 병합해본다.")

    # 센트로이드 클러스터 (옵션: 형태 conf별로 따로 → 모양 오염 제거)
    use_split = args.split_by_conf and any(conf.get(r["cif"], "") != "" for r, _ in data)
    if use_split:
        conf_ids = sorted(set(conf.get(r["cif"], "") for r, _ in data))
        clusters = []
        for cl_id in conf_ids:
            sub = [(r, c) for r, c in data if conf.get(r["cif"], "") == cl_id]
            for cl in greedy_cluster(sub, args.threshold):
                cl["conf"] = cl_id
                clusters.append(cl)
        print(f"[pocket] split_by_conf: {len(conf_ids)} conformation(s) → 형태별 분리 클러스터")
    else:
        clusters = greedy_cluster(data, args.threshold)
        for cl in clusters:
            cl["conf"] = "all"
    clusters.sort(key=lambda cl: len(cl["members"]), reverse=True)
    clusters = clusters[:args.topn]

    member_rows = []
    cand_rows = []
    for pid, cl in enumerate(clusters, 1):
        mems = cl["members"]
        cs = np.array([c for _, c in mems])
        # medoid 구조(센트로이드 평균에 가장 가까운)
        medoid = mems[int(np.argmin(np.linalg.norm(cs - cl["center"], axis=1)))][0]
        models, confs = {}, {}
        for r, _ in mems:
            models[r["model"]] = models.get(r["model"], 0) + 1
            cid = conf.get(r["cif"], "")
            if cid != "":
                confs[cid] = confs.get(cid, 0) + 1
            member_rows.append({"cif": r["cif"], "pocket_id": pid,
                                "model": r["model"], "iptm": r.get("iptm", ""),
                                "ligand_iptm": r.get("ligand_iptm", "")})
        cand_rows.append({
            "pocket_id": pid, "conf": cl.get("conf", "all"),
            "size": len(mems), "n_models": len(models),
            "models": ";".join(f"{k}:{v}" for k, v in sorted(models.items())),
            "conf_composition": ";".join(f"{k}:{v}" for k, v in sorted(confs.items())) or "NA",
            "center_x": f"{cl['center'][0]:.2f}", "center_y": f"{cl['center'][1]:.2f}",
            "center_z": f"{cl['center'][2]:.2f}",
            "rep_cif": medoid["cif"], "rep_model": medoid["model"],
            "rep_iptm": medoid.get("iptm", ""),
            "rep_ligand_iptm": medoid.get("ligand_iptm", "")})

    with open(os.path.join(args.out, "pocket_candidates.csv"), "w", newline="") as f:
        cols = ["pocket_id", "conf", "size", "n_models", "models", "conf_composition",
                "center_x", "center_y", "center_z", "rep_cif", "rep_model", "rep_iptm",
                "rep_ligand_iptm"]
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(cand_rows)
    with open(os.path.join(args.out, "members.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["cif", "pocket_id", "model", "iptm", "ligand_iptm"])
        w.writeheader(); w.writerows(member_rows)

    print(f"[pocket] {len(clusters)} pocket candidates (top {args.topn}, threshold {args.threshold}Å)")
    for c in cand_rows:
        print(f"  P{c['pocket_id']}: size {c['size']}, #mdl {c['n_models']}, {c['models']}, conf[{c['conf_composition']}]")
    print(f"  -> {args.out}/pocket_candidates.csv (+ members.csv)")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (먼저: source $CASP17/scripts/env_setup.sh)
#   SC=$CASP17/users/USERNAME/scripts ; OUT=$CASP17/users/USERNAME/targets/T2383
#   micromamba run -n boltz2 python $SC/2_pocket_candidates.py \
#       --table $OUT/00_collect/master_table.csv --out $OUT/02_pocket_candidates \
#       --protein-clusters $OUT/01_protein_clusters/protein_clusters.csv --threshold 8.0 --topn 10 --reference $OUT/01_protein_clusters/reference.txt
# (평소엔 run_pipeline.py 가 순서대로 자동 호출)
