#!/usr/bin/env python3
"""
protein_cluster.py (Step 0) - 단백질 구조 클러스터링 → 대표 conformation.
모델에 따라 다른 단백질 형태(예: ConfA/B)가 나오므로, Cα를 정렬해 PCA→그리디 클러스터로
형태를 라벨링하고 대표(medoid)를 저장. 이후 스텝이 conf를 인지/선택할 수 있게 함.
boltz2 env (numpy+gemmi). 출력: 01_protein_clusters/.
"""
import argparse, csv, os, shutil
import numpy as np
import gemmi
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio
import common.geom as geom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True, help="01_protein_clusters 폴더")
    ap.add_argument("--threshold", type=float, default=15.0, help="(폐기) 옛 그리디 절대임계 — 무시됨")
    ap.add_argument("--gap-frac", type=float, default=0.20,
                    help="PC1에서 '전체 범위 x 이 비율'보다 큰 틈에서만 형태 분리(범용; 좌표 단위 무관)")
    ap.add_argument("--min-conf-frac", type=float, default=0.02,
                    help="형태 그룹 최소 크기 비율(이보다 작아지는 분리는 무시 → 이상치 분리 방지)")
    ap.add_argument("--outlier-mad", type=float, default=3.5,
                    help="PCA 전 이상치 제거 강도(중앙값 + k*MAD 이상이면 이상치; 클수록 덜 제거)")
    ap.add_argument("--max-outlier-frac", type=float, default=0.15,
                    help="이상치로 제거할 수 있는 최대 비율(과도제거 방지)")
    ap.add_argument("--sep-ratio", type=float, default=2.0,
                    help="2형태 인정 임계: (두 centroid 거리)/(군집내 산포) 가 이 값 이상일 때만 2형태로 분리")
    ap.add_argument("--max", type=int, default=0, help="속도 위해 상위 N개만(0=전체)")
    ap.add_argument("--cache", default="", help="0b geom_cache.pkl(있으면 재파싱 없이 사용)")
    args = ap.parse_args()
    os.makedirs(os.path.join(args.out, "reps"), exist_ok=True)

    rows = [r for r in csv.DictReader(open(args.table)) if r.get("cif")]
    rows.sort(key=lambda r: int(r["rank"]))
    if args.max:
        rows = rows[:args.max]

    cache = cio.load_geom_cache(args.cache)
    ref_ca, _, _ = cio.cached_geometry(rows[0]["cif"], cache)
    ref_keys = sorted(ref_ca)

    feats, used = [], []
    for i, r in enumerate(rows):
        ca, _, _ = cio.cached_geometry(r["cif"], cache)
        rt = geom.align_to_ref(ref_ca, ca)
        if rt is None:
            continue
        R, t = rt
        # ref_keys 위치의 정렬된 Cα (없으면 제외 표시)
        vec, ok = [], True
        for k in ref_keys:
            if k in ca:
                vec.append(geom.apply_rt(R, t, [ca[k]])[0])
            else:
                ok = False; break
        if not ok:
            continue
        feats.append(np.array(vec).flatten())
        used.append(r)
        if (i + 1) % 500 == 0:
            print(f"  parsed {i+1}/{len(rows)}")

    X = np.array(feats)
    print(f"[protein_cluster] {len(X)} structures, feature dim {X.shape[1] if len(X) else 0}")
    if len(X) < 2:
        open(os.path.join(args.out, "protein_clusters.csv"), "w").write("note,insufficient\n")
        return

    # ── (1) 이상치(깨진/펼쳐진 구조) 제거 후 PCA ──
    # broken 구조가 분산을 장악하면 PC1이 '이상치 축'이 되어 형태분리가 안 됨.
    # robust center(median)로부터의 거리로 이상치를 걸러내고 inlier로만 PCA축을 정의.
    med = np.median(X, axis=0)
    dev = np.linalg.norm(X - med, axis=1)
    mad = np.median(np.abs(dev - np.median(dev))) + 1e-9
    inlier = dev <= np.median(dev) + args.outlier_mad * 1.4826 * mad
    if (~inlier).sum() > args.max_outlier_frac * len(X):       # 과도제거 방지
        inlier = dev <= np.quantile(dev, 1.0 - args.max_outlier_frac)
    n_out = int((~inlier).sum())

    mu = X[inlier].mean(0)
    _, S, Vt = np.linalg.svd(X[inlier] - mu, full_matrices=False)
    K = min(3, Vt.shape[0])
    var = (S ** 2) / (S ** 2).sum() * 100
    PC = (X - mu) @ Vt[:K].T            # 모든 구조를 같은 축에 투영
    PCin = PC[inlier]

    # ── (2) 1 vs 2 형태 판정: inlier에 2-means + 분리도 검정 ──
    # 형태 차이는 PC1이 아니라 PC2/PC3에도 실릴 수 있으므로 top-K 전체를 사용한다.
    # 두 centroid 거리가 군집 내 산포의 sep_ratio배 이상일 때만 2형태로 인정
    # (단일형태 타겟의 과분할 방지). 아니면 1형태(conf 0).
    def _kmeans2(P, iters=100):
        idx = np.argsort(P[:, 0])
        c = np.array([P[idx[len(P) // 10]], P[idx[-len(P) // 10]]], float)
        lab = np.full(len(P), -1)
        for _ in range(iters):
            nl = ((P[:, None, :] - c[None, :, :]) ** 2).sum(2).argmin(1)
            if (nl == lab).all():
                break
            lab = nl
            for k in range(2):
                if (lab == k).any():
                    c[k] = P[lab == k].mean(0)
        return lab, c

    min_size = max(5, int(args.min_conf_frac * len(PCin)))
    lab_in, cen = _kmeans2(PCin)
    sz_in = np.bincount(lab_in, minlength=2)
    sep = float(np.linalg.norm(cen[0] - cen[1]))
    spread = float(np.mean([PCin[lab_in == k].std(0).mean()
                            for k in range(2) if (lab_in == k).any()])) + 1e-9
    two_conf = (sz_in.min() >= min_size) and (sep >= args.sep_ratio * spread)

    if two_conf:                                    # 모든 구조를 가까운 centroid에 배정
        labels = ((PC[:, None, :] - cen[None, :, :]) ** 2).sum(2).argmin(1)
    else:
        labels = np.zeros(len(PC), dtype=int)
    print(f"[protein_cluster] outliers removed={n_out}, sep/spread={sep/spread:.2f} "
          f"(cut={args.sep_ratio}) -> {'2 conformations' if two_conf else '1 conformation'}")

    sizes = np.bincount(labels)
    order = np.argsort(sizes)[::-1]
    relabel = {old: new for new, old in enumerate(order)}

    # ── 정렬 기준 구조: inlier 중 PC중심에 가장 가까운 구조 ──
    inl_idx = np.where(inlier)[0]
    central = used[inl_idx[int(np.argmin(np.linalg.norm(PCin - PCin.mean(0), axis=1)))]]
    with open(os.path.join(args.out, "reference.txt"), "w") as rf:
        rf.write(central["cif"])
    print(f"[protein_cluster] reference(중심) = {os.path.basename(central['cif'])}")
    
    # 대표(medoid) 저장
    rep_cif = {}
    for new, old in enumerate(order):
        idx = np.where(labels == old)[0]
        sub = PC[idx]
        medoid_local = idx[int(np.argmin(np.linalg.norm(sub - sub.mean(0), axis=1)))]
        src = used[medoid_local]["cif"]
        dst = os.path.join(args.out, "reps", f"conf_{new}.cif")
        try:
            shutil.copy(src, dst)
        except Exception:
            pass
        rep_cif[new] = src

    with open(os.path.join(args.out, "protein_clusters.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cif", "model", "seed", "sample", "conf_id", "pc1", "pc2", "pc3"])
        for r, lab, p in zip(used, labels, PC):
            pc = list(p) + [0, 0, 0]
            w.writerow([r["cif"], r["model"], r["seed"], r.get("sample_idx", ""),
                        relabel[lab], f"{pc[0]:.2f}", f"{pc[1]:.2f}", f"{pc[2]:.2f}"])

    print(f"[protein_cluster] {len(sizes)} conformations (PC1 {var[0]:.1f}%, PC2 {var[1]:.1f}%)")
    for new, old in enumerate(order):
        print(f"  conf_{new}: {sizes[old]} structures, rep={os.path.basename(rep_cif[new])}")
    print(f"  -> {args.out}/protein_clusters.csv  (+ reps/)")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (먼저: source $CASP17/scripts/env_setup.sh)
#   SC=$CASP17/users/USERNAME/scripts ; OUT=$CASP17/users/USERNAME/targets/T2383
#   micromamba run -n boltz2 python $SC/1_protein_cluster.py \
#       --table $OUT/00_collect/master_table.csv --out $OUT/01_protein_clusters --gap-frac 0.20
# (평소엔 run_pipeline.py 가 순서대로 자동 호출)
