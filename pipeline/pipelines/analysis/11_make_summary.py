#!/usr/bin/env python3
"""
11_make_summary.py - 파이프라인 출력(CSV들)에서 타겟 정리 문서의 "사실 부분"만 자동 초안 생성.

사람이 손으로 표를 옮겨 적던 부분(포켓 후보/최종 poses/템플릿 overlap/contact residue)을
outputs CSV에서 그대로 뽑아 Markdown 초안으로 만든다. 해석·PyMOL 캡처·제출전략처럼
'사람의 판단'이 필요한 부분은 채우지 않고 `✍️` 자리표시자로만 남긴다(전사는 자동, 판단은 사람).

Markdown 으로 내보내므로 Google Docs→docx 변환 시 생기던 탭 깨짐이 없다.
stdlib 만 사용(gemmi/numpy 불필요) → 어떤 python 으로도 실행 가능.
출력: 09_summary/SUMMARY_DRAFT.md
"""
import argparse, csv, os


def load(path):
    return list(csv.DictReader(open(path, encoding="utf-8"))) if os.path.exists(path) else []


def read_smiles_tsv(path):
    """ligand.tsv → [(ID, Name, SMILES)]. complex_io(gemmi 의존) 대신 직접 읽음."""
    out = []
    if path and os.path.exists(path):
        with open(path, newline="", encoding="utf-8") as h:
            for r in csv.DictReader(h, delimiter="\t"):
                if r.get("SMILES"):
                    out.append((r.get("ID", ""), r.get("Name", ""), r["SMILES"]))
    return out


def read_progress(path):
    """00_collect/progress.txt → (total, {model:count})."""
    total, per = None, {}
    if os.path.exists(path):
        for ln in open(path, encoding="utf-8"):
            ln = ln.rstrip("\n")
            if ln.lower().startswith("total"):
                try:
                    total = int(ln.split(":")[1])
                except (IndexError, ValueError):
                    pass
            elif ":" in ln and ln.startswith(" "):
                k, v = ln.strip().split(":", 1)
                try:
                    per[k.strip()] = int(v)
                except ValueError:
                    pass
    return total, per


def md_table(rows, cols, headers=None):
    headers = headers or cols
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join(["---"] * len(cols)) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def fnum(x, d=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


HUMAN = "> ✍️ **[사람이 판단해 채우기]** "


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outputs", required=True, help="한 타겟의 outputs 폴더")
    ap.add_argument("--target", default="TARGET")
    ap.add_argument("--task", default="P")
    ap.add_argument("--ligand-tsv", default="")
    ap.add_argument("--out", default="", help="SUMMARY_DRAFT.md (기본 outputs/09_summary/)")
    args = ap.parse_args()
    O = args.outputs
    out = args.out or os.path.join(O, "09_summary", "SUMMARY_DRAFT.md")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    # ── 출력 CSV 로드 ──
    cand = load(os.path.join(O, "02_pocket_candidates/pocket_candidates.csv"))
    val = load(os.path.join(O, "03_pocket_validation/pocket_validation.csv"))
    ligc = load(os.path.join(O, "04_ligand_clusters/ligand_clusters.csv"))
    sel = load(os.path.join(O, "05b_refined/selection_summary.csv")) or \
          load(os.path.join(O, "05_final/selection_summary.csv"))
    conf_rows = load(os.path.join(O, "05c_confidence/confidence.csv"))
    conf = conf_rows[0] if conf_rows else {}
    contacts = load(os.path.join(O, "06_validation/contact_residues.csv"))
    pvt = load(os.path.join(O, "06_validation/pose_vs_template.csv"))
    frag = load(os.path.join(O, "06_validation/fragment_compare.csv"))
    total, per_model = read_progress(os.path.join(O, "00_collect/progress.txt"))
    smis = read_smiles_tsv(args.ligand_tsv)

    # 1순위(focus) 포켓
    focus = conf.get("focus_pocket") or (sel[0].get("pocket_id") if sel else "")
    verdict = conf.get("verdict", "(신뢰도 리포트 없음)")

    L = []  # 문서 라인
    L.append(f"# {args.target} 정리 — 파이프라인 & 템플릿 비교 (자동 초안)\n")
    L.append("> 이 문서는 `11_make_summary.py`가 파이프라인 출력에서 **사실만 자동 채운 초안**입니다.\n"
             "> 표/수치는 자동, `✍️` 표시는 **사람이 해석·판단해 채우는 자리**입니다. "
             "PyMOL 캡처·최종 순위 결정·제출전략은 사람 몫입니다.\n")

    # ── 0. 타겟 소개 ──
    L.append("\n## 0. 타겟 소개\n")
    L.append(f"- **Target**: {args.target}  |  **Task**: {args.task}")
    if smis:
        L.append(f"- **리간드**: {smis[0][1] or smis[0][0] or '?'} (SMILES는 맨 끝 [참고])")
    if total is not None:
        pm = ", ".join(f"{k}:{v}" for k, v in sorted(per_model.items()))
        L.append(f"- **샘플 수**: 총 {total}개 ({pm})")
    L.append(f"- **★ 최종 신뢰도 판정(5c)**: {verdict}")
    L.append(f"{HUMAN}단백질 이름/설명, 실험구조 존재 여부, 사용한 템플릿 목록\n")

    # ── 1-1. 포켓 후보 ──
    L.append("\n## 1. 파이프라인 고찰\n")
    L.append("### 1-1. 포켓 후보 (pocket_validation.csv)\n")
    L.append("- `p2rank_dist`: 예측 포켓 중심~최근접 캐비티 거리(작을수록 진짜 캐비티, pass ≤ 6Å 기본)")
    L.append("- `gnina_affinity`: 결합 에너지(음수일수록 안정, pass ≤ -4 기본) / `pass`: 둘 다 만족\n")
    src = val or cand
    if src:
        L.append(md_table(src,
                 ["pocket_id", "size", "n_models", "models", "p2rank_dist", "gnina_affinity", "pass"],
                 ["pocket", "size", "n_models", "models", "p2rank_dist", "gnina_aff", "pass"]))
    else:
        L.append("_(pocket_validation.csv 없음)_")
    if conf:
        L.append(f"\n- 1순위 포켓 **P{focus}**: 우세도 {conf.get('dominance','?')}, "
                 f"합의 {conf.get('n_models','?')}모델, af3 {'포함' if conf.get('af3_present')=='True' else '없음'}, "
                 f"포켓분리 {conf.get('separation','?')}x")
    L.append(f"\n{HUMAN}포켓 위치가 한 곳으로 몰렸는지, pass 여부의 의미(위치 문제 vs 검증기준 문제)\n")

    # ── 1-2. 최종 poses ──
    L.append("\n### 1-2. 최종 제출 poses (selection_summary.csv)\n")
    if sel:
        L.append(md_table(sel,
                 ["model", "pocket_id", "ligand_cluster_id", "size", "iptm",
                  "gnina_affinity", "composite", "pocket_pass", "posebusters_valid"],
                 ["model", "pocket", "lig_cl", "size", "iptm", "gnina_aff",
                  "composite", "pocket_pass", "posebusters"]))
    else:
        L.append("_(selection_summary.csv 없음)_")
    L.append(f"\n{HUMAN}5개가 같은 포켓인지/방향만 다른지, composite 1위가 압도적인지\n")

    # ── 1-3. 방향 수렴 ──
    L.append("\n### 1-3. 방향 수렴 (ligand_clusters.csv)\n")
    focus_lc = [r for r in ligc if r.get("pocket_id") == str(focus)]
    if focus_lc:
        sizes = sorted((int(r.get("size", 0)) for r in focus_lc), reverse=True)
        # 분모 = 포켓 전체 size (5c/문서 정의와 동일). 없으면 클러스터 합으로 폴백.
        psize = next((int(r.get("size", 0)) for r in (val or cand)
                      if r.get("pocket_id") == str(focus)), 0) or sum(sizes)
        conv = (sizes[0] / psize) if psize else 0
        L.append(f"- P{focus} 방향 클러스터 {len(focus_lc)}개, 최대 {sizes[0]}/{psize} = **{conv:.0%} 수렴** "
                 f"(높을수록 방향까지 일치).")
    elif conf.get("pose_convergence"):
        L.append(f"- pose 수렴(5c): {float(conf['pose_convergence']):.0%}")
    else:
        L.append("_(ligand_clusters.csv 없음)_")
    L.append(f"\n{HUMAN}수렴도가 다른 타겟 대비 높/낮은지, 낮다면 원인(리간드 유연/대칭 등)\n")

    # ── 2. 템플릿 비교 ──
    L.append("\n## 2. 템플릿 비교 (pose_vs_template.csv)\n")
    if pvt:
        templates = sorted(set(r.get("template", "") for r in pvt))
        L.append(f"- 템플릿: {', '.join(t for t in templates if t)}")
        L.append("- `overlap_pct`: 우리 리간드 원자 중 실험 리간드 4Å 이내 비율 / "
                 "`shared`·`jaccard`: 접촉잔기 겹침 / `sup_rmsd`: CA 정렬 RMSD\n")
        L.append(md_table(pvt,
                 ["model", "template", "sup_rmsd", "overlap_pct", "shared", "jaccard", "convergence"],
                 ["model", "template", "sup_rmsd", "overlap%", "shared", "jaccard", "conv"]))
        # composite 순위 vs 최고 overlap (docx 4-3 핵심 표)
        comp = {r["model"]: r.get("composite", "") for r in sel}
        bym = {}
        for r in pvt:
            ov = fnum(r.get("overlap_pct"))
            if ov is None:
                continue
            bym.setdefault(r["model"], []).append((ov, r.get("template", "")))
        if bym:
            L.append("\n**composite 순위 vs 실험 overlap (자동 랭킹이 실험 최적 포즈를 1위로 뽑았나?)**\n")
            rank = sorted(bym.items(), key=lambda kv: -max(x[0] for x in kv[1]))
            rows = []
            for m, v in rank:
                best = max(v)
                rows.append({"model": m, "composite": comp.get(m, ""),
                             "best_overlap": f"{best[0]:.0f}%", "vs_template": best[1]})
            L.append(md_table(rows, ["model", "composite", "best_overlap", "vs_template"],
                              ["model", "composite", "best overlap", "vs template"]))
    else:
        L.append("_(pose_vs_template.csv 없음 — config `templates=` 설정 시 자동 생성)_")
    L.append(f"\n{HUMAN}PyMOL 육안 비교(캡처 첨부), overlap이 낮/높은 이유, 위치는 맞는데 pose가 어긋나는지\n")

    # ── 3. contact residue ──
    L.append("\n## 3. contact residue (contact_residues.csv)\n")
    if contacts:
        top = sorted(contacts, key=lambda r: -fnum(r.get("contact_frequency"), 0))[:15]
        s = ", ".join(f"{r.get('chain','')}{r.get('residue','')}"
                      f"({float(r.get('contact_frequency',0))*100:.0f}%)" for r in top)
        L.append(f"- 우리 예측(P{focus}) 고빈도 접촉잔기: {s}")
    else:
        L.append("_(contact_residues.csv 없음)_")
    if frag:
        L.append("\n- 실험 fragment 공유(fragment_compare.csv):")
        L.append(md_table(frag, ["exp_ligand", "n_shared", "n_exp"],
                          ["exp_ligand", "shared", "n_exp"]))
    L.append(f"\n{HUMAN}실험 리간드 접촉잔기와 공유 개수(같은 포켓 확정 근거), 위 2장 overlap과 일치하는지\n")

    # ── 4. 제출 전략 (사람) ──
    L.append("\n## 4. 제출 전략\n")
    L.append(f"{HUMAN}자동 composite 순위를 그대로 낼지, overlap/contact 근거로 순위를 손볼지(인간 선택). "
             "그 근거와 주의점(예: 실험 리간드 ≠ 타겟 리간드라면 얼마나 신뢰?)\n")

    # ── 5. 요약 (사람) ──
    L.append("\n## 5. 요약\n")
    L.append(f"{HUMAN}consensus/검증/템플릿/overlap/contact를 종합한 한 문단 결론과 신뢰도(HIGH/MEDIUM/LOW) 근거\n")

    # ── 참고: 리간드 ──
    if smis:
        L.append("\n## [참고] 리간드\n")
        for sid, name, smi in smis:
            L.append(f"- **{name or sid or 'ligand'}**: `{smi}`")

    open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
    print(f"[summary] draft -> {out}")
    print("  (facts auto-filled; fill the marked sections by human judgment)")


if __name__ == "__main__":
    main()

# ── 단독 실행 ── (아무 python; stdlib)
#   python3 11_make_summary.py --outputs $OUT --target T2412 --task P --ligand-tsv <ligand_tsv>
#   # -> $OUT/09_summary/SUMMARY_DRAFT.md (Markdown; docx 변환 없이 그대로 열람/붙여넣기)
# (평소엔 run_pipeline.py 가 맨 끝에 자동 호출)
