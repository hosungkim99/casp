#!/usr/bin/env python3
"""
visualize.py (범용, 워크플로 7번) - 선별 결과 시각화. complex_io 사용 → 리간드 개수 무관.
boltz2 env (numpy+gemmi+matplotlib). 그래프 텍스트 영어, 주석 한글.

출력 1) <out>: 4패널 (리간드 PCA by model / by cluster / 모델별 iptm 분포 / 상위클러스터 구성)
출력 2) --contacts 주면 contact fingerprint PNG (per-residue 접촉빈도; fragment 공유잔기 강조)
출력 3) --templates+--summary+--cluster-members 주면 final pose별 × 템플릿별 접촉지문 격자 PNG
        (contact_fp_<model>.png; x=그 pose 클러스터 접촉빈도, 템플릿 리간드 공유잔기=빨강)
"""
import argparse, csv, os, glob
import numpy as np
import gemmi
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys, os  # 부트스트랩: 상위 pipeline/ 를 path에 추가 → common 패키지 import 가능
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio

PALETTE = ["tab:blue", "tab:green", "tab:orange", "tab:red", "tab:purple", "tab:brown"]


def kabsch(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, Qc - R @ Pc


def color_for(models):
    uniq = sorted(set(models))
    return {m: PALETTE[i % len(PALETTE)] for i, m in enumerate(uniq)}


# ── 출력 3: 템플릿 접촉지문용 헬퍼 (final pose별 × 템플릿별) ──
def _main_chain(model, minlen=40):
    best = None
    for ch in model:
        pol = ch.get_polymer()
        if pol and len(pol) >= minlen and (best is None or len(pol) > len(best.get_polymer())):
            best = ch
    return best


def _polymer_chains(model, minlen=40):
    return [ch for ch in model if ch.get_polymer() and len(ch.get_polymer()) >= minlen]


def _heavy_ligand_atoms(model, near_chain=None, near_cutoff=8.0):
    ca = None
    if near_chain is not None:
        ca = np.array([[a.pos.x, a.pos.y, a.pos.z] for r in near_chain for a in r if a.name == "CA"])
    coords = []
    for ch in model:
        for r in ch:
            t = gemmi.find_tabulated_residue(r.name)
            if (t and (t.is_amino_acid() or t.is_nucleic_acid())) or r.name in ("HOH", "WAT"):
                continue
            heavy = [a for a in r if a.element.name != "H"]
            if "C" not in [a.element.name for a in heavy]:
                continue                     # 이온 등(탄소 없음) 제외
            xs = [[a.pos.x, a.pos.y, a.pos.z] for a in heavy]
            if ca is not None and len(xs) and len(ca):
                d = np.sqrt(((ca[:, None, :] - np.array(xs)[None, :, :]) ** 2).sum(-1))
                if d.min() > near_cutoff:
                    continue
            coords += xs
    return np.array(coords) if coords else np.zeros((0, 3))


def _protein_heavy(chain):
    coords, labels = [], []
    for r in chain:
        t = gemmi.find_tabulated_residue(r.name)
        if not (t and t.is_amino_acid()):
            continue
        for a in r:
            if a.element.name == "H":
                continue
            coords.append([a.pos.x, a.pos.y, a.pos.z]); labels.append(r.seqid.num)
    return np.array(coords) if coords else np.zeros((0, 3)), np.array(labels)


def _apply_T(T, xyz):
    out = []
    for p in xyz:
        v = T.mat.multiply(gemmi.Vec3(*p)) + T.vec
        out.append([v.x, v.y, v.z])
    return np.array(out) if out else np.zeros((0, 3))


def _contact_resnums(prot_xyz, prot_lab, query_xyz, cutoff):
    if len(query_xyz) == 0 or len(prot_xyz) == 0:
        return set()
    d = np.sqrt(((prot_xyz[:, None, :] - query_xyz[None, :, :]) ** 2).sum(-1))
    return set(int(x) for x in prot_lab[(d <= cutoff).any(1)])


def _cluster_contact_freq(cifs, cutoff):
    """클러스터 멤버들의 per-residue 접촉빈도 {resnum: frac}, 멤버 수."""
    freq, n = {}, 0
    for cif in cifs:
        try:
            st = gemmi.read_structure(cif); st.setup_entities(); m = st[0]
        except Exception:
            continue
        ns = gemmi.NeighborSearch(m, st.cell, cutoff + 1).populate()
        hit = set()
        for ch in m:
            for r in ch:
                t = gemmi.find_tabulated_residue(r.name)
                if (t and (t.is_amino_acid() or t.is_nucleic_acid())) or r.name in ("HOH", "WAT"):
                    continue
                if "C" not in [a.element.name for a in r if a.element.name != "H"]:
                    continue
                for a in r:
                    if a.element.name == "H":
                        continue
                    for mk in ns.find_atoms(a.pos, "\0", radius=cutoff):
                        rr = mk.to_cra(m).residue
                        tt = gemmi.find_tabulated_residue(rr.name)
                        if tt and tt.is_amino_acid():
                            hit.add(rr.seqid.num)
        for rn in hit:
            freq[rn] = freq.get(rn, 0) + 1
        n += 1
    n = max(n, 1)
    return {rn: c / n for rn, c in freq.items()}, n


def _template_contact_set(model_cif, tmpl_cif, cutoff):
    """템플릿을 우리 모델에 중첩 → 템플릿 리간드 접촉 '우리 잔기번호' 집합. 리간드 없으면 None."""
    ms = gemmi.read_structure(model_cif); ms.setup_entities()
    ts = gemmi.read_structure(tmpl_cif); ts.setup_entities()
    mm, tm = ms[0], ts[0]
    our = _main_chain(mm)
    if our is None:
        return None
    prot_xyz, prot_lab = _protein_heavy(our)
    best = None
    for tch in _polymer_chains(tm):
        try:
            sup = gemmi.calculate_superposition(our.get_polymer(), tch.get_polymer(),
                                                gemmi.PolymerType.PeptideL, gemmi.SupSelect.CaP)
        except Exception:
            continue
        if best is None or sup.rmsd < best[0].rmsd:
            best = (sup, tch)
    if best is None:
        return None
    sup, tch = best
    lig = _heavy_ligand_atoms(tm, near_chain=tch)
    if len(lig) == 0:
        return None                          # apo/이온만
    return _contact_resnums(prot_xyz, prot_lab, _apply_T(sup.transform, lig), cutoff)


def template_fingerprints(args):
    """출력 3: final pose별 PNG(리간드 있는 템플릿마다 subplot). 인자 다 있을 때만 호출."""
    tmpls = [t.strip() for t in args.templates.split(",") if t.strip()]
    m2pc = {r["model"]: (r.get("pocket_id", ""), r.get("ligand_cluster_id", ""))
            for r in csv.DictReader(open(args.summary))}
    members = {}
    for r in csv.DictReader(open(args.cluster_members)):
        members.setdefault((r["pocket_id"], r["ligand_cluster_id"]), []).append(
            (int(r["member_rank"]), r["cif"]))
    for k in members:
        members[k].sort()
    outdir = os.path.dirname(args.out) or "."
    made = 0
    for cif in sorted(glob.glob(os.path.join(args.final_dir, "model_*.cif"))):
        name = os.path.splitext(os.path.basename(cif))[0]
        key = m2pc.get(name)
        if not key or key not in members:
            continue
        mem = [c for _, c in members[key]]
        if args.max_members:
            mem = mem[:args.max_members]
        freq, nmem = _cluster_contact_freq(mem, args.fp_cutoff)
        data = sorted([(rn, f) for rn, f in freq.items() if f >= args.fp_threshold])
        if not data:
            continue
        xs = [rn for rn, _ in data]
        panels = []
        for tp in tmpls:
            tset = _template_contact_set(cif, tp, args.fp_cutoff)
            if tset is not None:
                panels.append((os.path.splitext(os.path.basename(tp))[0], tset))
        if not panels:
            continue
        ncol = min(2, len(panels)); nrow = (len(panels) + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(7.5 * ncol, 4.2 * nrow), squeeze=False)
        fig.suptitle(f"{name} contact fingerprint  (cluster n={nmem}; red = shared with template ligand)",
                     fontsize=13)
        ys = [f * 100 for _, f in data]
        for i, (tn, tset) in enumerate(panels):
            ax = axes[i // ncol][i % ncol]
            cols = ["crimson" if rn in tset else "steelblue" for rn in xs]
            nsh = sum(1 for rn in xs if rn in tset)   # 공유(우리 접촉 잔기 중 템플릿도 접촉)
            nour, ntmpl = len(xs), len(tset)          # 우리 / 템플릿 접촉 잔기 수
            recall = (nsh / ntmpl) if ntmpl else 0.0  # 템플릿 대비(재현율)
            jac = (nsh / (nour + ntmpl - nsh)) if (nour + ntmpl - nsh) else 0.0  # 합집합 대비(공정 비교)
            ax.bar([str(x) for x in xs], ys, color=cols)
            ax.set_title(f"vs {tn}: shared {nsh}  (our {nour} / tmpl {ntmpl}, "
                         f"recall {recall:.0%}, jaccard {jac:.2f})", fontsize=10)
            ax.set_ylabel("contact frequency (%)"); ax.set_xlabel("residue number"); ax.set_ylim(0, 105)
            for lab in ax.get_xticklabels():
                lab.set_rotation(60); lab.set_fontsize(7)
        for j in range(len(panels), nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        plt.tight_layout()
        fp = os.path.join(outdir, f"contact_fp_{name}.png")
        plt.savefig(fp, dpi=130); plt.close(fig)
        print(f"[viz] saved: {fp} ({len(panels)} templates)")
        made += 1
    print(f"[viz] template fingerprints: {made} PNG")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", required=True, help="4패널 PNG")
    ap.add_argument("--contacts", default="", help="contact_residues.csv (있으면 fingerprint PNG)")
    ap.add_argument("--fragment", default="", help="fragment_compare.csv (공유잔기 강조)")
    ap.add_argument("--threshold", type=float, default=2.0)
    ap.add_argument("--reference", default="", help="01의 reference.txt(PC중심 구조; 없으면 rank1 폴백)")
    ap.add_argument("--cache", default="", help="0b geom_cache.pkl(있으면 재파싱 없이 사용)")
    ap.add_argument("--protein-clusters", default="",
                    help="01 protein_clusters.csv (있으면 단백질 형태 PCA = v1/v2 그래프도 출력)")
    # 출력 3(템플릿 접촉지문): 아래 4개 다 주면 final pose별 격자 PNG 생성
    ap.add_argument("--templates", default="", help="쉼표구분 템플릿 cif(리간드 포함)")
    ap.add_argument("--summary", default="", help="05(b) selection_summary.csv (model→pocket/cluster)")
    ap.add_argument("--cluster-members", default="", help="04 cluster_members.csv")
    ap.add_argument("--final-dir", default="", help="05_final 또는 05b_refined (model_*.cif)")
    ap.add_argument("--fp-cutoff", type=float, default=4.0, help="접촉지문 거리")
    ap.add_argument("--fp-threshold", type=float, default=0.2, help="접촉지문 x축 최소 빈도")
    ap.add_argument("--max-members", type=int, default=0, help="빈도 계산 클러스터 멤버 상한(0=전체)")
    args = ap.parse_args()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)

    rows = [r for r in csv.DictReader(open(args.table)) if r.get("cif")]
    rows.sort(key=lambda r: int(r["rank"]))
    cache = cio.load_geom_cache(args.cache)
    ref_ca, ref_elems, _ = cio.cached_geometry(cio.reference_cif(rows, args.reference), cache)

    vecs, models, iptms, cifs = [], [], [], []
    for i, r in enumerate(rows):
        ca, elems, coords = cio.cached_geometry(r["cif"], cache)
        if not ca or elems != ref_elems:
            continue
        common = [k for k in ref_ca if k in ca]
        if len(common) < 50:
            continue
        R, t = kabsch(np.array([ca[k] for k in common]), np.array([ref_ca[k] for k in common]))
        L = (R @ np.array(coords).T).T + t
        vecs.append(L.flatten()); models.append(r["model"]); cifs.append(r["cif"])
        try:
            iptms.append((r["model"], float(r["iptm"])))
        except (TypeError, ValueError):
            pass
        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{len(rows)}")
    X = np.array(vecs); models = np.array(models)
    print(f"[viz] {len(X)} structures")
    MCOL = color_for(models)

    Xc = X - X.mean(0)
    _, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    PC = Xc @ Vt[:2].T
    var = (S ** 2) / (S ** 2).sum() * 100

    reps, labels = [], []
    for v in X:
        for ci, rep in enumerate(reps):
            if np.sqrt(((v - rep) ** 2).reshape(-1, 3).sum(1).mean()) <= args.threshold:
                labels.append(ci); break
        else:
            reps.append(v); labels.append(len(reps) - 1)
    labels = np.array(labels); sizes = np.bincount(labels)
    top = np.argsort(sizes)[::-1][:6]

    fig, ax = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle("Ligand pose selection overview", fontsize=15)
    for mdl in sorted(set(models)):
        sel = models == mdl
        ax[0, 0].scatter(PC[sel, 0], PC[sel, 1], s=10, alpha=0.5, c=MCOL[mdl], label=mdl)
    ax[0, 0].set_title(f"(1) Ligand PCA by model (PC1 {var[0]:.1f}%, PC2 {var[1]:.1f}%)")
    ax[0, 0].legend()
    # (2) 하이브리드: 전역 PCA 유지 + pocket_id로 색칠(04/05 선정과 일치) +
    #     최종 선정 모델 별표 강조 + 검증탈락/이상치는 회색.
    pocket_of, clust_of, rank_of = {}, {}, {}   # cif -> pocket_id / ligand_cluster / member_rank
    if args.cluster_members and os.path.exists(args.cluster_members):
        for cr in csv.DictReader(open(args.cluster_members)):
            c = cr.get("cif")
            if c:
                pocket_of[c] = cr.get("pocket_id", "")
                clust_of[c] = cr.get("ligand_cluster_id", "")
                rank_of[c] = cr.get("member_rank", "")
    selected_cl = set()  # 선정된 (pocket, ligand_cluster) — source_cif 대신 클러스터로 매칭
    if args.summary and os.path.exists(args.summary):     #   (refine on/off 무관하게 성립)
        for sr in csv.DictReader(open(args.summary)):
            selected_cl.add((sr.get("pocket_id", ""), sr.get("ligand_cluster_id", "")))
    pk_list = [pocket_of.get(c, "") for c in cifs]           # PC 순서와 동일
    gray = np.array([p == "" for p in pk_list])
    if gray.any():
        ax[0, 1].scatter(PC[gray, 0], PC[gray, 1], s=8, alpha=0.25,
                         c="lightgray", label="rejected/other")
    uniq_pk = sorted({p for p in pk_list if p},
                     key=lambda s: int(s) if s.isdigit() else 999)
    for j, pk in enumerate(uniq_pk):
        m = np.array([p == pk for p in pk_list])
        ax[0, 1].scatter(PC[m, 0], PC[m, 1], s=14, alpha=0.6,
                         color=PALETTE[j % len(PALETTE)], label=f"pocket {pk} (n={int(m.sum())})")
    # 별표: 선정 클러스터의 대표(member_rank 1) pose (refine 여부와 무관하게 앙상블에서 찾음)
    selmask = np.array([(pocket_of.get(c, ""), clust_of.get(c, "")) in selected_cl
                        and rank_of.get(c, "") == "1" for c in cifs])
    if selmask.any():
        ax[0, 1].scatter(PC[selmask, 0], PC[selmask, 1], s=90, marker="*",
                         facecolors="none", edgecolors="black", linewidths=1.2,
                         label=f"final models (n={int(selmask.sum())})")
    ax[0, 1].set_title("(2) Ligand PCA by pocket  (final models = *)")
    ax[0, 1].legend(fontsize=7)
    by = {}
    for mdl, v in iptms:
        by.setdefault(mdl, []).append(v)
    names = sorted(by)
    ax[1, 0].boxplot([by[n] for n in names], labels=names)
    ax[1, 0].set_title("(3) iptm distribution by model (calibration)"); ax[1, 0].set_ylabel("iptm")
    topc = np.argsort(sizes)[::-1][:8]
    bottom = np.zeros(len(topc))
    for mdl in sorted(set(models)):
        counts = [int(np.sum((labels == ci) & (models == mdl))) for ci in topc]
        ax[1, 1].bar(range(len(topc)), counts, bottom=bottom, color=MCOL[mdl], label=mdl)
        bottom += counts
    ax[1, 1].set_xticks(range(len(topc)))
    ax[1, 1].set_xticklabels([f"c{ci+1}" for ci in topc])
    ax[1, 1].set_title("(4) Model composition of top clusters"); ax[1, 1].legend()
    plt.tight_layout(); plt.savefig(args.out, dpi=130)
    print(f"[viz] saved: {args.out}")

    # ── 단백질 형태(conf) PCA : CASP v1/v2 시각화 (01 protein_clusters.csv 기반) ──
    # 리간드 PCA(위 4패널)는 리간드 위치만 봐서 v1/v2가 안 보임. 단백질 Cα PCA는 형태를 직접 보여줌.
    if args.protein_clusters and os.path.exists(args.protein_clusters):
        prows = [r for r in csv.DictReader(open(args.protein_clusters)) if r.get("pc1")]
        if prows:
            Pp = np.array([[float(r["pc1"]), float(r["pc2"]), float(r.get("pc3", 0) or 0)]
                           for r in prows])
            cids = [r.get("conf_id", "0") for r in prows]
            uconf = sorted(set(cids))
            pal = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
            fig2, ax2 = plt.subplots(1, 2, figsize=(14, 6))
            for pj, (a, b, lbl) in enumerate([(0, 1, "PC1 vs PC2"), (0, 2, "PC1 vs PC3")]):
                for k, cf in enumerate(uconf):
                    m = [i for i, c in enumerate(cids) if c == cf]
                    ax2[pj].scatter(Pp[m, a], Pp[m, b], s=12, alpha=0.6,
                                    c=pal[k % len(pal)], label=f"conf {cf} (n={len(m)})")
                ax2[pj].set_xlabel(f"protein PC{a+1}"); ax2[pj].set_ylabel(f"protein PC{b+1}")
                ax2[pj].set_title(lbl); ax2[pj].legend()
            fig2.suptitle(f"PROTEIN-conformation PCA -> v1/v2  ({len(uconf)} conformation(s))",
                          fontweight="bold")
            fig2.tight_layout()
            pconf = os.path.join(os.path.dirname(args.out) or ".", "protein_conformation_pca.png")
            fig2.savefig(pconf, dpi=130)
            print(f"[viz] saved: {pconf}")

    # contact fingerprint (선택)
    if args.contacts and os.path.exists(args.contacts):
        crows = list(csv.DictReader(open(args.contacts)))
        frag_res = set()
        if args.fragment and os.path.exists(args.fragment):
            for fr in csv.DictReader(open(args.fragment)):
                for x in (fr.get("shared_with_prediction", "") or "").split(";"):
                    if x.strip().isdigit():
                        frag_res.add(int(x))
        data = []
        for r in crows:
            try:
                data.append((int(r["residue"]), float(r["contact_frequency"])))
            except (ValueError, KeyError):
                pass
        data = [d for d in data if d[1] >= 0.2]
        data.sort(key=lambda x: x[0])
        if data:
            xs = [str(k) for k, _ in data]; ys = [v * 100 for _, v in data]
            cols = ["crimson" if k in frag_res else "steelblue" for k, _ in data]
            plt.figure(figsize=(14, 6))
            plt.bar(xs, ys, color=cols)
            plt.ylabel("contact frequency (%)"); plt.xlabel("residue number")
            plt.title("Ligand-contact fingerprint  |  red = shared with experimental fragment")
            plt.xticks(rotation=60); plt.tight_layout()
            fp = os.path.join(os.path.dirname(args.out), "contact_fingerprint.png")
            plt.savefig(fp, dpi=130)
            print(f"[viz] saved: {fp}")

    # 출력 3: final pose별 × 템플릿별 접촉지문 (인자 다 주어졌을 때만; 템플릿 없으면 skip)
    if args.templates and args.summary and args.cluster_members and args.final_dir \
            and os.path.exists(args.summary) and os.path.exists(args.cluster_members):
        template_fingerprints(args)


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (먼저: source $CASP17/scripts/env_setup.sh)
#   SC=$CASP17/users/USERNAME/scripts ; OUT=$CASP17/users/USERNAME/targets/T2383
#   micromamba run -n boltz2 python $SC/8_visualize.py \
#       --table $OUT/00_collect/master_table.csv --out $OUT/07_viz/visualization.png \
#       --contacts $OUT/06_validation/contact_residues.csv \
#       --reference $OUT/01_protein_clusters/reference.txt \
#       [--fragment $OUT/06_validation/fragment_compare.csv]
# (평소엔 run_pipeline.py 가 순서대로 자동 호출)
