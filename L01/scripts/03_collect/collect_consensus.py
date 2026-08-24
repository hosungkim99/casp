#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
collect_consensus.py - 4개 co-folding 모델(pt2/bt2/af3/of3)의 native 출력을 읽어
binder 1개당 '모든 모델 합본' master_table.csv 한 장으로 정리(=T2410/00_collect 포맷과 동일).

출력 컬럼(기존 0_rank_poses.py 결과와 완전 동일):
  rank, model, seed, sample_idx, ligand_iptm, iptm, ptm, plddt, gpde,
  has_clash, ranking_score_native, cif, summary
  - model 컬럼으로 pt2/bt2/af3/of3 구분, 전 모델 샘플을 한 표에 랭킹.
  - ligand_iptm = chain_pair_iptm off-diagonal 평균(체인 구성 무관). 없으면 native ligand_iptm/iptm 폴백.
  - 정렬키 = (clash없음, ligand_iptm, plddt, -gpde) 내림차순 → 0_rank_poses.py와 동일.
  - cif/summary = 원본 native 파일 경로를 그대로 참조(복사·symlink 없음, 원본 무수정).

각 binder마다 <out-root>/<binder>/00_collect/master_table.csv 생성(전부 users/USERNAME 내부).
stdlib만. Linux(서버)에서 실행.
"""
import argparse, csv, glob, json, os, re, sys
from statistics import mean
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

EMBEDDED_BINDERS = ("L010016 L010039 L010061 L010069 L010078 L010087 L010099 L010128 L010144 "
 "L010167 L010210 L010223 L010281 L010307 L010309 L010316 L010319 L010322 L010327 L010331 "
 "L010337 L010356 L010363 L010365 L010397 L010412 L010432 L010443 L010462 L010537 L010547 "
 "L010552 L010553 L010589 L010594 L010621 L010626 L010630 L010639 L010649 L010662 L010669 "
 "L010685 L010695 L010702 L010712 L010728 L010738 L010761 L010770 L010782 L010801 L010807 "
 "L010886 L010888 L010906 L010912 L010918 L010919 L010930 L010939 L010943 L010984 L010993 "
 "L011019 L011043 L011057 L011070 L011110 L011124 L011140 L011159 L011160 L011166 L011167 "
 "L011177 L011179 L011199 L011207").split()

COLS = ["rank", "model", "seed", "sample_idx", "ligand_iptm", "iptm", "ptm",
        "plddt", "gpde", "has_clash", "ranking_score_native", "cif", "summary"]


def interface_iptm(cpi):
    """chain_pair_iptm off-diagonal 평균. dict 또는 리스트행렬 모두 허용."""
    if isinstance(cpi, list) and cpi and isinstance(cpi[0], list):
        vals = [v for i, row in enumerate(cpi) for j, v in enumerate(row)
                if i != j and isinstance(v, (int, float))]
        return mean(vals) if vals else None
    if isinstance(cpi, dict):
        vals = [v for a, row in cpi.items() if isinstance(row, dict)
                for b, v in row.items() if a != b and isinstance(v, (int, float))]
        return mean(vals) if vals else None
    return None


def load_existing(path):
    """기존 master_table.csv(예: 팀원 bt2)를 그대로 읽어 행(문자열 dict) 리스트로. rank는 무시(재계산)."""
    rows = []
    if path and os.path.exists(path):
        with open(path, newline="") as f:
            for r in csv.DictReader(f):
                d = {k: r.get(k) for k in COLS}
                d["rank"] = None
                rows.append(d)
    return rows


def _seed_of(path):
    m = re.search(r"[/\\]seed[_-]?(\d+)[/\\]", path)
    return int(m.group(1)) if m else None


# ── 모델별 수집: (row dict) 리스트 반환. cif/summary는 원본 native 경로 ────────────
def collect_pt2(binder, base):
    """Protenix. 평탄화 leaf <base>/<id>/[results]/pt2/seed_<S>/sample_<M>/{summary.json, model.cif}.
    summary.json 키가 AF3식(iptm/chain_pair_iptm/ptm/plddt/gpde/has_clash/ranking_score)이라 af3 리더 재사용."""
    return _collect_af3like(binder, base, "pt2")


def collect_bt2(binder, base):
    """Boltz2. 평탄화 leaf <base>/<id>/[results]/bt2/seed_<S>/sample_<M>/{summary.json, model.cif}.
    summary.json = boltz raw confidence(그대로). 매핑: plddt<-complex_plddt, gpde<-complex_pde,
    ranking_score_native<-confidence_score, has_clash 없음.
    ligand_iptm: (pair_)chain_pair_iptm off-diag, 없으면 native ligand_iptm, 그것도 없으면 iptm."""
    rows = []
    root = os.path.join(base, binder, RESULTS_SUBDIR, "bt2")
    def _seed(p):
        m = re.search(r"seed[_-]?(\d+)", p); return int(m.group(1)) if m else 0
    sjs = glob.glob(os.path.join(root, "**", "summary.json"), recursive=True)
    if MAX_SEEDS > 0:
        keep = sorted(set(_seed(p) for p in sjs))[:MAX_SEEDS]
        sjs = [p for p in sjs if _seed(p) in keep]
    for sj in sorted(sjs, key=lambda p: (_seed(p), p)):
        seed = _seed(sj)
        mm = re.search(r"sample_(\d+)", sj); midx = int(mm.group(1)) if mm else None
        cif = os.path.join(os.path.dirname(sj), "model.cif")
        if midx is None or not os.path.exists(cif):
            continue
        try:
            d = json.load(open(sj))
        except Exception as e:
            sys.stderr.write(f"[bt2] skip {sj}: {e}\n"); continue
        lig = interface_iptm(d.get("chain_pair_iptm") or d.get("pair_chains_iptm"))
        if lig is None:
            lig = d.get("ligand_iptm", d.get("iptm"))
        rows.append({"model": "bt2", "seed": seed, "sample_idx": midx,
                     "ligand_iptm": lig, "iptm": d.get("iptm"), "ptm": d.get("ptm"),
                     "plddt": d.get("complex_plddt"), "gpde": d.get("complex_pde"),
                     "has_clash": None,
                     "ranking_score_native": d.get("confidence_score"),
                     "cif": cif, "summary": sj})
    return rows


MAX_SEEDS = 6   # af3/of3 seed 상한(bt2/pt2 6 seed와 밸런스). main에서 --max-seeds로 덮음.
RESULTS_SUBDIR = ""   # 모델 읽기경로에 끼울 중간폴더. binders 구조면 "results"(→ <base>/<id>/results/<model>/).


def _collect_af3like(binder, base, model):
    """AF3/OF3. 통합 레이아웃 <base>/<id>/<model>/ 아래 재귀검색.
    seed_<S>/sample_<M>/{summary.json, model.cif}. seed는 경로 seed_<S>, MAX_SEEDS 상한(앞 N seed).
    summary.json은 pt2와 동일한 AF3식 키(iptm/chain_pair_iptm/ptm/plddt/gpde/has_clash/ranking_score)."""
    rows = []
    root = os.path.join(base, binder, RESULTS_SUBDIR, model)
    def _seed(p):
        m = re.search(r"seed[_-]?(\d+)", p); return int(m.group(1)) if m else 0
    sjs = glob.glob(os.path.join(root, "**", "summary.json"), recursive=True)
    if MAX_SEEDS > 0:
        keep = sorted(set(_seed(p) for p in sjs))[:MAX_SEEDS]
        sjs = [p for p in sjs if _seed(p) in keep]
    for sj in sorted(sjs, key=lambda p: (_seed(p), p)):
        seed = _seed(sj)
        mm = re.search(r"sample_(\d+)", sj); midx = int(mm.group(1)) if mm else None
        cif = os.path.join(os.path.dirname(sj), "model.cif")
        if midx is None or not os.path.exists(cif):
            continue
        try:
            d = json.load(open(sj))
        except Exception as e:
            sys.stderr.write(f"[{model}] skip {sj}: {e}\n"); continue
        rows.append({"model": model, "seed": seed, "sample_idx": midx,
                     "ligand_iptm": interface_iptm(d.get("chain_pair_iptm")),
                     "iptm": d.get("iptm"), "ptm": d.get("ptm"),
                     "plddt": d.get("plddt"), "gpde": d.get("gpde"),
                     "has_clash": d.get("has_clash"),
                     "ranking_score_native": d.get("ranking_score"),
                     "cif": cif, "summary": sj})
    return rows


def collect_af3(binder, base):
    return _collect_af3like(binder, base, "af3")


def collect_of3(binder, base):
    return _collect_af3like(binder, base, "of3")


COLLECTORS = {"pt2": collect_pt2, "bt2": collect_bt2, "af3": collect_af3, "of3": collect_of3}


def _num(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def _is_clash(x):
    if isinstance(x, bool): return x
    if x is None: return False
    return str(x).strip().lower() in ("true", "1", "yes")


def rank_key(r):
    """기존행(문자열)·신규행(native) 혼용 안전. 정렬: clash없음 > ligand_iptm > plddt > -gpde."""
    li = _num(r.get("ligand_iptm"))
    if li is None:
        li = _num(r.get("iptm")) or 0
    gp = _num(r.get("gpde"))
    return (0 if _is_clash(r.get("has_clash")) else 1, li, _num(r.get("plddt")) or 0,
            -(gp if gp is not None else 9))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-root", required=True,
                    help="출력 루트(users/USERNAME 내부). <out-root>/<binder>/00_collect/master_table.csv 생성. "
                         "--existing-root와 같게 주면 그 파일에 그대로 덮어써(=기존 테이블에 pt2 추가)")
    ap.add_argument("--existing-root", default="",
                    help="기존 master_table.csv 루트(예: .../L01/outputs). <root>/<binder>/00_collect/master_table.csv "
                         "를 읽어 그 위에 새 모델 행을 더함(팀원 bt2 재활용). 새로 긁는 --models와 겹치는 모델은 자동 대체.")
    ap.add_argument("--models", default="pt2", help="쉼표구분, '추가로 긁을' 모델 (pt2,bt2,af3,of3)")
    ap.add_argument("--results-base", default="",
                    help="통합 결과 루트. 주면 전 모델이 <base>/<id>/[results-subdir]/<model>/ 에서 읽음.")
    ap.add_argument("--results-subdir", default="",
                    help="모델 읽기경로 중간폴더(binders 구조면 'results' → <base>/<id>/results/<model>/)")
    ap.add_argument("--out-subdir", default="",
                    help="테이블 쓰기경로 중간폴더(binders 구조면 'consensus_s2' → <out-root>/<id>/consensus_s2/00_collect)")
    ap.add_argument("--pt2-base", default="", help="(개별지정) pt2 base. <base>/<id>/pt2/ 재귀")
    ap.add_argument("--bt2-runs", default="", help="(개별지정) bt2 base. <base>/<id>/bt2/ 재귀")
    ap.add_argument("--af3-base", default="", help="(개별지정) af3 base. <base>/<id>/af3/ 재귀")
    ap.add_argument("--of3-base", default="", help="(개별지정) of3 base. <base>/<id>/of3/ 재귀")
    ap.add_argument("--max-seeds", type=int, default=6,
                    help="af3/of3 seed 상한(0=전량). bt2/pt2 6 seed와 밸런스 위해 기본 6")
    ap.add_argument("--binders", default="", help="binder ID 파일(줄당 1개). 미지정시 내장 79개")
    ap.add_argument("--only", default="", help="쉼표구분 binder만 (검증용, 예: L010016)")
    args = ap.parse_args()

    binders = list(EMBEDDED_BINDERS)
    if args.binders and os.path.exists(args.binders):
        binders = [x.strip() for x in open(args.binders) if x.strip().startswith("L01")]
    if args.only:
        # --only는 준 binder를 그대로 사용(임베드 79 목록과 무관 → 신규 L010123 등도 처리)
        binders = [x.strip() for x in args.only.split(",") if x.strip()]

    global MAX_SEEDS, RESULTS_SUBDIR
    MAX_SEEDS = args.max_seeds
    RESULTS_SUBDIR = args.results_subdir
    models = [m.strip() for m in args.models.split(",") if m.strip()]
    # --results-base 주면 전 모델 base로 사용(<base>/<id>/<model>/). 개별 --*-base가 있으면 그게 우선.
    rb = args.results_base
    bases = {"pt2": args.pt2_base or rb, "bt2": args.bt2_runs or rb,
             "af3": args.af3_base or rb, "of3": args.of3_base or rb}
    for m in models:
        if not bases[m]:
            sys.stderr.write(f"[경고] --{m}-* base 미지정 → {m} 건너뜀\n")

    collect_set = {m for m in models if bases[m]}
    grand = {m: 0 for m in models}
    for b in binders:
        rows = []
        if args.existing_root:
            # 기존 테이블 재활용. 단, 지금 새로 긁는 모델과 겹치면(재실행) 기존 것 버리고 새 걸로.
            ex = load_existing(os.path.join(args.existing_root, b, "00_collect", "master_table.csv"))
            rows.extend([r for r in ex if r.get("model") not in collect_set])
        for m in models:
            if bases[m]:
                mr = COLLECTORS[m](b, bases[m])
                rows.extend(mr); grand[m] += len(mr)
        rows.sort(key=rank_key, reverse=True)
        for i, r in enumerate(rows, 1):
            r["rank"] = i
        out_dir = os.path.join(args.out_root, b, args.out_subdir, "00_collect")
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "master_table.csv"), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS); w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k) for k in COLS})
        cnt = {}
        for r in rows:
            cnt[r["model"]] = cnt.get(r["model"], 0) + 1
        print(f"{b}: total={len(rows)}  " + " ".join(f"{m}={c}" for m, c in sorted(cnt.items())))

    print("---")
    print("새로 긁은 모델 총 샘플: " + " ".join(f"{m}={grand[m]}" for m in models) + f"   (binder {len(binders)}개)")


if __name__ == "__main__":
    main()
