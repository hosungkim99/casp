#!/usr/bin/env python3
"""
12_score_table.py - 클러스터별 raw 스코어 테이블 생성 (medoid + average + 분산).

목적: 흩어진 pipeline 출력(클러스터·포켓·포즈·gnina)을 한데 모아 보기 편한 표 하나로.
**팀 회의록(CASP17-ligand 기본 프로토콜)에 맞춘 3종 row**:
  - pocket_average : 회의록 1-1 — "각 **포켓**에 속한 리간드들에 Vina·Gnina → 평균 낮은 곳이 포켓".
                     후보 포켓 **전체**(02_pocket_candidates/members.csv)의 멤버 평균 + 분산.
  - medoid         : 포켓 내 포즈 클러스터의 **구조 대표**(모양상 중심) pose 실제 값.
                     ※ 회의록 2의 "가장 스코어가 좋은 것"(=top_score)과 다름(composite 준비 후 별도 행 추가 예정).
  - average        : 회의록 2 — 그 포즈 클러스터 전 멤버 평균 (+ 분산 std). "개수·평균·분산" 대응.

컬럼 배치: [1] 분산(평균행만) → [2] 공통(모든행) → [3] pocket_residues·source_cif(medoid행·긴 데이터, 맨 뒤).

점수 방향: gnina_affinity(=Vina) 낮을수록 좋음 / cnn_score·cnn_affinity 높을수록 좋음.
gnina/cnn 은 pipeline 이 포켓 대표에만 계산 → 여기서 **전 멤버에 gnina 직접 실행**(--minimize 옵션).
boltz_affinity: bt2 포즈에만 존재(그 예측의 affinity json). af3/of3/pt2 는 NA. composite 에선 제외.
"""
import sys, os, csv, argparse, glob, json, statistics as st

_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import numpy as np
import gemmi
import common.complex_io as cio
import common.geom as geom
import common.scoring as scoring

WATER = getattr(cio, "WATER", {"HOH", "WAT", "H2O"})

COLS = [  # 출력 컬럼 순서 (= 빈 행 템플릿의 기준)
    # 식별/구조
    "pocket_id", "cluster_id", "cluster_size", "model_composition", "model_status",
    # raw 결합 점수 + 각 분산
    "gnina_affinity", "cnn_score", "cnn_affinity",
    "gnina_affinity_std", "cnn_score_std", "cnn_affinity_std", "rmsd",
    # 신뢰도 · 파생
    "ligand_iptm", "boltz_affinity", "composite",
    # 검증 · 긴 데이터
    "p2rank_dist", "posebusters_pass", "pocket_residues", "source_cif"]


# ---- 작은 유틸 ---------------------------------------------------------------
def load_csv(path):
    return list(csv.DictReader(open(path))) if os.path.exists(path) else []


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def mean(xs):
    xs = [x for x in xs if x is not None]
    return round(sum(xs) / len(xs), 4) if xs else ""


def std(xs):
    xs = [x for x in xs if x is not None]
    return round(st.pstdev(xs), 4) if len(xs) >= 2 else (0.0 if xs else "")


def row(**kw):
    """모든 컬럼을 ''로 채운 뒤 준 값만 덮어쓴 행 dict (빈 필드 반복 제거)."""
    return {**{c: "" for c in COLS}, **kw}


def make_norm(vals):
    """표 전체 기준 min-max 0~1 정규화 함수. None은 통과."""
    v = [x for x in vals if x is not None]
    if not v:
        return lambda x: None
    lo, hi = min(v), max(v)
    return lambda x: (None if x is None else ((x - lo) / (hi - lo) if hi > lo else 0.5))


# ---- 포즈 단위 계산 ----------------------------------------------------------
def _prot_lig_atoms(cif):
    """cif → (단백질 잔기별 heavy atom 좌표 dict, 리간드 heavy atom 좌표 list)."""
    st_ = gemmi.read_structure(cif)
    if len(st_) == 0:
        return {}, []
    prot, lig = {}, []
    for ch in st_[0]:
        for res in ch:
            t = gemmi.find_tabulated_residue(res.name)
            if t and t.is_amino_acid():
                c = [[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H"]
                if c:
                    prot[(ch.name, res.seqid.num)] = np.array(c)
            elif (t and t.is_nucleic_acid()) or res.name in WATER:
                continue
            else:
                lig += [[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H"]
    return prot, lig


def pocket_res(cif, cutoff):
    """리간드 heavy atom 에서 단백질 heavy atom(all-atom)이 cutoff 이내인 잔기 → 'A_190'."""
    prot, lig = _prot_lig_atoms(cif)
    if not prot or not lig:
        return ""
    L = np.array(lig)
    hit = [f"{ch}_{num}" for (ch, num), c in prot.items()
           if np.sqrt(((c[:, None, :] - L[None, :, :]) ** 2).sum(2)).min() <= cutoff]
    return ",".join(sorted(hit))


def gnina_one(cif, work, gpu_id, minimize=False):
    """멤버 cif 하나에 gnina → {aff,cnn,cnnaff}. minimize=True면 국소 최소화 후 값(clash 완화)."""
    tag = os.path.splitext(os.path.basename(cif))[0]
    prot, lig = os.path.join(work, f"{tag}_prot.pdb"), os.path.join(work, f"{tag}_lig.pdb")
    _, ligs = cio.parse_complex(cif)
    if not ligs:
        return {}
    cio.write_ligands_pdb(ligs, lig)
    cio.write_receptor_pdb(cif, prot)
    return scoring.score_with_fallback(prot, lig, gpu_id, minimize)


def aligned_lig(cif, ref_ca):
    """cif 리간드 좌표를 ref_ca(=medoid Cα)에 단백질 정렬해서 반환. 정렬 실패 시 원좌표."""
    ca, ligs = cio.parse_complex(cif)
    if not ligs:
        return None
    _, coords = cio.ligand_concat(ligs)   # (elements, coords) 튜플 → 좌표만
    L = np.asarray(coords)
    rt = geom.align_to_ref(ref_ca, ca) if ref_ca else None
    return geom.apply_rt(rt[0], rt[1], L) if rt else L


def rmsd_to_medoid(medoid, mems, smiles):
    """멤버들을 medoid 에 단백질 정렬 후 SC-RMSD 평균(구조 분산). 실패 멤버는 skip."""
    try:
        ref_ca, _ = cio.parse_complex(medoid)
        ref = aligned_lig(medoid, None)
        perms = geom.ligand_automorphisms(smiles, len(ref)) if smiles else None
        ds = []
        for m in mems[1:]:
            L = aligned_lig(m["cif"], ref_ca)
            if L is not None and len(L) == len(ref):
                ds.append(geom.sc_rmsd(L, ref, perms))
        return round(sum(ds) / len(ds), 3) if ds else 0.0
    except Exception:
        return ""


def boltz_aff(cif, model, key, gpat):
    """bt2 포즈면 affinity json 을 찾아 값 반환. 아니면 ''."""
    if model != "bt2":
        return ""
    for base in (os.path.dirname(cif), os.path.dirname(os.path.dirname(cif))):
        for jf in glob.glob(os.path.join(base, gpat)):
            try:
                j = json.load(open(jf))
                if key in j:
                    return j[key]
            except Exception:
                pass
    return ""


# ---- 클러스터/포켓 단위 집계 -------------------------------------------------
def model_comp(cifs, mt):
    """멤버들의 모델 구성 문자열: 'bt2:3;pt2:2'."""
    comp = {}
    for c in cifs:
        mo = mt.get(c, {}).get("model", "?")
        comp[mo] = comp.get(mo, 0) + 1
    return ";".join(f"{k}:{v}" for k, v in sorted(comp.items()))


def pb_frac(cifs, pb):
    """멤버 중 posebusters valid 비율 (검사된 것 기준). 없으면 ''."""
    v = [pb[c].get("valid") for c in cifs if c in pb]
    return round(sum(1 for x in v if x in ("True", "1", True)) / len(v), 2) if v else ""


def avg_row(cifs, ctx, **fixed):
    """포켓/클러스터 멤버들의 평균+분산 행 (pocket_average·average 공용)."""
    gsc, mt, pb, pose_comp, args = ctx["gsc"], ctx["mt"], ctx["pb"], ctx["pose_comp"], ctx["args"]
    g = lambda k: [fnum(gsc.get(c, {}).get(k)) for c in cifs]
    aff, cnn, caf = g("aff"), g("cnn"), g("cnnaff")
    lip = [fnum(mt.get(c, {}).get("ligand_iptm")) for c in cifs]
    blz = [fnum(boltz_aff(c, mt.get(c, {}).get("model", ""),
                          args.boltz_affinity_key, args.boltz_affinity_glob)) for c in cifs]
    return row(model_composition=model_comp(cifs, mt),
               ligand_iptm=mean(lip), boltz_affinity=mean(blz),
               gnina_affinity=mean(aff), cnn_score=mean(cnn), cnn_affinity=mean(caf),
               composite=mean([pose_comp.get(c) for c in cifs]), posebusters_pass=pb_frac(cifs, pb),
               gnina_affinity_std=std(aff), cnn_score_std=std(cnn), cnn_affinity_std=std(caf),
               **fixed)


def medoid_row(medoid, mcomp, ctx, **fixed):
    """medoid(구조 대표 pose) 실제값 행."""
    gsc, mt, pb, pose_comp, args = ctx["gsc"], ctx["mt"], ctx["pb"], ctx["pose_comp"], ctx["args"]
    mg, pbm = mt.get(medoid, {}), pb.get(medoid, {})
    return row(model_composition=mcomp, model_status="medoid",
               ligand_iptm=mg.get("ligand_iptm", ""),
               boltz_affinity=boltz_aff(medoid, mg.get("model", ""),
                                        args.boltz_affinity_key, args.boltz_affinity_glob),
               gnina_affinity=gsc[medoid].get("aff", ""), cnn_score=gsc[medoid].get("cnn", ""),
               cnn_affinity=gsc[medoid].get("cnnaff", ""), composite=pose_comp.get(medoid, ""),
               posebusters_pass=pbm.get("valid", ""),
               pocket_residues=pocket_res(medoid, args.pocket_cutoff), source_cif=medoid, **fixed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", required=True, help="한 타겟 outputs 폴더")
    ap.add_argument("--cid", default="", help="타겟 id (기본: outputs 상위 폴더명)")
    ap.add_argument("--ligand-smiles", default="", help="SC-RMSD 대칭보정용 SMILES(없으면 plain)")
    ap.add_argument("--out", default="", help="출력 csv (기본: <outputs>/score_table.csv)")
    ap.add_argument("--append", action="store_true",
                    help="누적 모드: --out 파일에 이어붙임(맨 앞 cid 컬럼 추가, 헤더는 새 파일일 때만).")
    ap.add_argument("--pocket-cutoff", type=float, default=5.0, help="pocket_residues all-atom 거리")
    ap.add_argument("--cluster-cutoff", default="", help="클러스터링 cutoff (메타로 기록만)")
    ap.add_argument("--max-pocket-members", type=int, default=0,
                    help="포켓당 gnina 채점 멤버 상한(0=전체). 계산량 조절용")
    ap.add_argument("--minimize", action="store_true",
                    help="gnina --minimize (국소 최소화 후 채점 → clash 오염 완화). 계산 시간↑")
    ap.add_argument("--jobs", type=int, default=1, help="gnina 병렬 워커(기본 1=순차)")
    ap.add_argument("--gpus", default="", help="라운드로빈 GPU id 예: 0,1")
    ap.add_argument("--boltz-affinity-key", default="affinity_probability_binary",
                    help="bt2 affinity json 에서 읽을 필드")
    ap.add_argument("--boltz-affinity-glob", default="affinity*.json",
                    help="cif 폴더(및 상위)에서 찾을 affinity json 패턴")
    args = ap.parse_args()

    O = args.outputs
    cid = args.cid or os.path.basename(os.path.dirname(O.rstrip("/")))
    out = args.out or os.path.join(O, "score_table.csv")

    members = load_csv(os.path.join(O, "04_ligand_clusters", "cluster_members.csv"))
    if not members:
        raise SystemExit(f"[score_table] cluster_members.csv 없음: {O} (통과 포켓/클러스터 필요)")
    lclust = {(r["pocket_id"], r["ligand_cluster_id"]): r
              for r in load_csv(os.path.join(O, "04_ligand_clusters", "ligand_clusters.csv"))}
    mt = {r["cif"]: r for r in load_csv(os.path.join(O, "00_collect", "master_table.csv")) if r.get("cif")}
    pval = {r["pocket_id"]: r for r in load_csv(os.path.join(O, "03_pocket_validation", "pocket_validation.csv"))}
    pb = {r.get("rep_cif", ""): r for r in load_csv(os.path.join(O, "04b_posebusters", "posebusters.csv"))}

    # 회의록 1-1 대응: 후보 포켓 전체의 멤버 (포켓 단위 평균 비교용)
    by_pocket = {}
    for r in load_csv(os.path.join(O, "02_pocket_candidates", "members.csv")):
        by_pocket.setdefault(r["pocket_id"], []).append(r)
    if args.max_pocket_members > 0:
        by_pocket = {k: v[:args.max_pocket_members] for k, v in by_pocket.items()}

    # 1) gnina 대상 = 후보 포켓 전 멤버 ∪ 포즈클러스터 멤버 (중복 제거 후 1회씩)
    cifs = sorted({m["cif"] for m in members} | {r["cif"] for v in by_pocket.values() for r in v})
    work = os.path.join(O, "12_score_table_work")
    os.makedirs(work, exist_ok=True)
    gpus = scoring.parse_gpus(args.gpus)
    print(f"[score_table] {cid}: 멤버 {len(members)}, gnina 대상 cif {len(cifs)} (jobs={args.jobs}, gpus={gpus})")
    fns = [(lambda c: (lambda gpu_id: gnina_one(c, work, gpu_id, args.minimize)))(c) for c in cifs]
    gsc = {c: (r or {}) for c, r in zip(cifs, scoring.run_jobs(fns, args.jobs, gpus))}

    # 2) composite용 포즈별 점수: 표 전체 min-max 0~1 + 부호 통일(전부 "높을수록 좋음"). boltz 제외.
    fv = make_norm([(-fnum(gsc[c].get("aff"))) if fnum(gsc[c].get("aff")) is not None else None for c in cifs])
    fc = make_norm([fnum(gsc[c].get("cnn")) for c in cifs])
    fa = make_norm([fnum(gsc[c].get("cnnaff")) for c in cifs])
    pose_comp = {}
    for c in cifs:
        a, s2, a2 = fnum(gsc[c].get("aff")), fnum(gsc[c].get("cnn")), fnum(gsc[c].get("cnnaff"))
        parts = [p for p in (fv(-a if a is not None else None), fc(s2), fa(a2)) if p is not None]
        pose_comp[c] = round(sum(parts) / len(parts), 4) if parts else None

    ctx = dict(gsc=gsc, mt=mt, pb=pb, pose_comp=pose_comp, args=args)
    rows = []

    # (A) 포켓 단위 (회의록 1-1): 후보 포켓 전 멤버 평균
    for pid, pm in sorted(by_pocket.items(), key=lambda kv: -len(kv[1])):
        cs = [r["cif"] for r in pm]
        rows.append(avg_row(cs, ctx, pocket_id=pid, cluster_id="", cluster_size=len(cs),
                            model_status="pocket_average",
                            p2rank_dist=pval.get(pid, {}).get("p2rank_dist", "")))

    # (B) 포즈 클러스터 단위 (회의록 2): medoid(대표) + average(평균+분산)
    by_cl = {}
    for m in members:
        by_cl.setdefault((m["pocket_id"], m["ligand_cluster_id"]), []).append(m)
    for (pid, lcid), mems in sorted(by_cl.items(), key=lambda kv: -len(kv[1])):
        mems.sort(key=lambda m: int(m.get("member_rank", 999)))
        cs = [m["cif"] for m in mems]
        csize = int(lclust.get((pid, lcid), {}).get("size") or len(mems))
        p2d = pval.get(pid, {}).get("p2rank_dist", "")
        mcomp = model_comp(cs, mt)
        base = dict(pocket_id=pid, cluster_id=lcid, cluster_size=csize, p2rank_dist=p2d)
        rows.append(medoid_row(cs[0], mcomp, ctx, **base))
        rows.append(avg_row(cs, ctx, model_status="average",
                            rmsd=rmsd_to_medoid(cs[0], mems, args.ligand_smiles), **base))

    if args.append:
        # 누적 모드: cid 컬럼을 맨 앞에 붙이고, 기존 파일에 이어붙임(헤더는 새 파일일 때만).
        cols = ["cid"] + COLS
        for r in rows:
            r["cid"] = cid
        new_file = not os.path.exists(out) or os.path.getsize(out) == 0
        os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
        with open(out, "a", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            if new_file:
                w.writeheader()
            w.writerows(rows)
        print(f"[score_table] += {out}  (cid={cid}, +{len(rows)}행)")
    else:
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=COLS)
            w.writeheader()
            w.writerows(rows)
        print(f"[score_table] -> {out}  ({len(rows)}행, 클러스터 {len(by_cl)})")


if __name__ == "__main__":
    main()
