#!/usr/bin/env python3
"""
run_pipeline.py - 계층형 CASP 리간드 선별 파이프라인 (교수님 Step 0~4 프레임워크).

config(경로 몇 개 + task)만 주면 각 스텝을 순서대로 자동 실행하고, 각 스텝 결과를
번호별 폴더에 저장(중간 결과 열람 가능). 산출물 있으면 skip(idempotent).

흐름:
  00 collect(rank) → 01 protein_cluster(Step0) → 02 pocket_candidates(Step1)
  → 03 pocket_validate(Step2) → 04 ligand_cluster(Step3) → 05 final_select(Step4)
  → 06 contact/fragment 검증 → 07 visualize → 08 casp_lg

사용:
  source /path/to/casp17-ligand/scripts/env_setup.sh
  python3 run_pipeline.py config/T2383.conf
  python3 run_pipeline.py config/T2383.conf --from pocket_validate --force
"""
import argparse, os, shlex, subprocess, sys, csv, glob

# 스텝 스크립트(core/boost/analysis)는 항상 이 파일과 같은 폴더에 붙어 다님.
# → 스크립트 위치 = 이 파일 위치. config에 scripts_dir를 둘 필요 없음(안 읽음).
_SELF_DIR = os.path.dirname(os.path.abspath(__file__))


def load_conf(path):
    c = {}
    for line in open(path, encoding="utf-8"):
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            k, v = line.split("=", 1); c[k.strip()] = v.strip()
    # scripts_dir는 더 이상 필수 아님 — 스텝 경로는 _SELF_DIR로 자동 해석.
    for k in ["target", "results_dir", "out_dir", "ligand_tsv"]:
        if not c.get(k):
            sys.exit(f"[config] 필수 키 누락/빈값: {k}")
    d = dict(protein_chain="", task="P", conf_gap_frac="0.20", split_by_conformation="",
             exclude_models="",
             pocket_threshold="8.0", pocket_sweep="",
             topn_pockets="10", ligand_threshold="2.0", gnina_cutoff="-4.0",
             p2rank_dist="6.0", pocket_members_topk="5", template_pass_dist="8.0",
             top_final="5", diversity_topk="10",
             author=os.environ.get("CASP_AUTHOR", "CHANGE-ME"),
             python_stdlib="python3", python_sci="micromamba run -n boltz2 python",
             experimental_cif="", python_pb="", templates="", templates_dir="",
             gnina_jobs="1", gnina_gpus="", pocket_auto_relax="", posebusters_topk="3",
             refine="", refine_iters="25", refine_max_drift="2.0")
    for k, v in d.items():
        c.setdefault(k, v)
    return c


def build_steps(c):
    O = c["out_dir"]
    P = lambda *a: os.path.join(O, *a)
    std, sci = c["python_stdlib"], c["python_sci"]
    D = dict(collect="00_collect", prot="01_protein_clusters", pock="02_pocket_candidates",
             val="03_pocket_validation", lig="04_ligand_clusters", pb="04b_posebusters",
             fin="05_final", ref="05b_refined", v6="06_validation", v7="07_viz",
             lg="08_casp_lg", conf="05c_confidence", sum9="09_summary")

    rank_args = ["--results", c["results_dir"], "--out", P(D["collect"])]
    if c["protein_chain"]:
        rank_args += ["--protein-chain", c["protein_chain"]]
    if c.get("exclude_models"):
        rank_args += ["--exclude-models", c["exclude_models"]]
    MT = P(D["collect"], "master_table.csv")
    CACHE = P(D["collect"], "geom_cache.pkl")   # 0b가 만든 기하 캐시(1·2·6·8 공유)

    pock_args = ["--table", MT, "--out", P(D["pock"]),
                 "--protein-clusters", P(D["prot"], "protein_clusters.csv"),
                 "--reference", P(D["prot"], "reference.txt"),
                 "--cache", CACHE,
                 "--threshold", c["pocket_threshold"], "--topn", c["topn_pockets"]]
    if str(c.get("split_by_conformation", "")).lower() in ("1", "true", "yes", "on"):
        pock_args.append("--split-by-conf")
    if c.get("pocket_sweep", "").strip():
        pock_args += ["--sweep", c["pocket_sweep"]]
    pv_relax = (["--auto-relax"]
                if str(c.get("pocket_auto_relax", "")).lower() in ("1", "true", "yes", "on") else [])

    # 템플릿 목록: templates=(명시) 우선, 없으면 templates_dir 안의 *.cif 전부.
    #   포켓검증(3, template_dist)·후보대조(10)·시각화(8)가 공유.
    tmpl = c.get("templates", "").strip()
    if not tmpl and c.get("templates_dir", "").strip():
        tmpl = ",".join(sorted(glob.glob(os.path.join(c["templates_dir"], "*.cif"))))
    pv_tmpl = (["--templates", tmpl, "--reference", P(D["prot"], "reference.txt"),
                "--template-pass-dist", c["template_pass_dist"]] if tmpl else [])

    steps = [
        dict(name="collect", env=std, script="core/0_rank_poses.py", args=rank_args, out=MT),
        # 기하 캐시: 모든 cif 1회 파싱 → 1·2·6·8 재파싱 제거(gemmi 파싱 병목 해소)
        dict(name="cache", env=sci, script="core/0b_cache_geometry.py",
             args=["--table", MT, "--out", CACHE], out=CACHE),
        # Step 0
        dict(name="protein_cluster", env=sci, script="core/1_protein_cluster.py",
             args=["--table", MT, "--out", P(D["prot"]), "--gap-frac", c["conf_gap_frac"],
                   "--cache", CACHE],
             out=P(D["prot"], "protein_clusters.csv")),
        # Step 1
        dict(name="pocket_candidates", env=sci, script="core/2_pocket_candidates.py",
             args=pock_args,
             out=P(D["pock"], "pocket_candidates.csv")),
        # Step 2
        dict(name="pocket_validate", env=sci, script="core/3_pocket_validate.py",
             args=["--candidates", P(D["pock"], "pocket_candidates.csv"), "--out", P(D["val"]),
                   "--members", P(D["pock"], "members.csv"),
                   "--members-topk", c["pocket_members_topk"],
                   "--gnina-cutoff", c["gnina_cutoff"], "--p2rank-dist", c["p2rank_dist"],
                   "--jobs", c["gnina_jobs"], "--gpus", c["gnina_gpus"]] + pv_relax + pv_tmpl,
             out=P(D["val"], "pocket_validation.csv")),
        # Step 3
        dict(name="ligand_cluster", env=sci, script="core/4_ligand_cluster.py",
             args=["--members", P(D["pock"], "members.csv"),
                   "--validation", P(D["val"], "pocket_validation.csv"),
                   "--candidates", P(D["pock"], "pocket_candidates.csv"),
                   "--table", MT, "--ligand-tsv", c["ligand_tsv"],
                   "--out", P(D["lig"]), "--threshold", c["ligand_threshold"]],
             out=P(D["lig"], "ligand_clusters.csv")),
        # Step 4
        dict(name="final_select", env=sci, script="core/5_final_select.py",
             args=["--ligand-clusters", P(D["lig"], "ligand_clusters.csv"),
                   "--validation", P(D["val"], "pocket_validation.csv"),
                   "--out", P(D["fin"]), "--top", c["top_final"], "--task", c["task"],
                   "--protein-clusters", P(D["prot"], "protein_clusters.csv"),
                   "--jobs", c["gnina_jobs"], "--gpus", c["gnina_gpus"]],
             out=P(D["fin"], "selection_summary.csv")),
        # 검증: contact residues
        dict(name="contact_residues", env=sci, script="core/6_contact_residues.py",
             args=["--table", MT, "--out", P(D["v6"], "contact_residues.csv"),
                   "--reference", P(D["prot"], "reference.txt"), "--cache", CACHE,
                   "--pocket-candidates", P(D["pock"], "pocket_candidates.csv")],
             out=P(D["v6"], "contact_residues.csv")),
    ]

    # ① PoseBusters 유효성 (python_pb 설정 시에만): prep(boltz2/gemmi) → bust(casp_eval/posebusters)
    if c.get("python_pb"):
        pb_csv = P(D["pb"], "posebusters.csv")
        prep = dict(name="prep_validity", env=sci, script="boost/4b_prep_validity.py",
                    args=["--ligand-clusters", P(D["lig"], "ligand_clusters.csv"),
                          "--cluster-members", P(D["lig"], "cluster_members.csv"),
                          "--topk", c["posebusters_topk"],
                          "--ligand-tsv", c["ligand_tsv"], "--out", P(D["pb"])],
                    out=P(D["pb"], "manifest.csv"))
        bust = dict(name="posebusters", env=c["python_pb"], script="boost/4c_posebusters.py",
                    args=["--manifest", P(D["pb"], "manifest.csv"), "--out", pb_csv],
                    out=pb_csv)
        fi = next(i for i, s in enumerate(steps) if s["name"] == "final_select")
        steps.insert(fi, bust)
        steps.insert(fi, prep)    # 순서: prep → bust → final_select
        for s in steps:
            if s["name"] == "final_select":
                s["args"] += ["--posebusters", pb_csv]
    # ② 포즈 정제 (refine=true): final_select 다음에 끼우고, casp_lg는 정제본 사용
    casp_summary = P(D["fin"], "selection_summary.csv")
    final_dir = P(D["fin"])                       # 최종 model_*.cif 위치(정제 시 갱신)
    if str(c.get("refine", "")).lower() in ("1", "true", "yes", "on"):
        refine_step = dict(name="refine", env=sci, script="boost/5b_refine.py",
            args=["--final-dir", P(D["fin"]), "--ligand-tsv", c["ligand_tsv"],
                  "--out", P(D["ref"]), "--max-drift", c["refine_max_drift"], "--iters", c["refine_iters"]],
            out=P(D["ref"], "selection_summary.csv"))
        fi = next(i for i, s in enumerate(steps) if s["name"] == "contact_residues")
        steps.insert(fi, refine_step)
        casp_summary = P(D["ref"], "selection_summary.csv")
        final_dir = P(D["ref"])
    # 검증: 실험 fragment 비교 (experimental_cif 있을 때만)
    if c["experimental_cif"]:
        steps.append(dict(name="fragment_compare", env=sci, script="analysis/7_fragment_compare.py",
             args=["--model1", P(D["fin"], "model_1.cif"), "--experimental", c["experimental_cif"],
                   "--our-contacts", P(D["v6"], "contact_residues.csv"),
                   "--out", P(D["v6"], "fragment_compare.csv"),
                   "--overlay", P(D["v6"], "overlay_fragment.pdb")],
             out=P(D["v6"], "fragment_compare.csv")))
    # 검증: 후보별 템플릿 대조 — 다중포켓 타겟에서 진짜 상위 후보 식별
    # (tmpl 은 위에서 이미 해석 = templates 명시 우선, 없으면 templates_dir/*.cif).
    # 비면 스텝 생략(템플릿 아직 없을 때 나머지는 정상 실행).
    if tmpl:
        steps.append(dict(name="pose_vs_template", env=sci, script="analysis/10_pose_vs_template.py",
             args=["--final-dir", final_dir, "--templates", tmpl,
                   "--summary", casp_summary,
                   "--ligand-clusters", P(D["lig"], "ligand_clusters.csv"),
                   "--out", P(D["v6"], "pose_vs_template.csv")],
             out=P(D["v6"], "pose_vs_template.csv")))
    # 시각화 (8_visualize: 4패널 + 접촉지문. templates 있으면 final pose별×템플릿별 접촉지문도 생성)
    viz = ["--table", MT, "--out", P(D["v7"], "visualization.png"),
           "--contacts", P(D["v6"], "contact_residues.csv"),
           "--protein-clusters", P(D["prot"], "protein_clusters.csv"),
           "--reference", P(D["prot"], "reference.txt"), "--cache", CACHE]
    if c["experimental_cif"]:
        viz += ["--fragment", P(D["v6"], "fragment_compare.csv")]
    if tmpl:
        viz += ["--templates", tmpl, "--summary", casp_summary,
                "--cluster-members", P(D["lig"], "cluster_members.csv"),
                "--final-dir", final_dir]
    steps.append(dict(name="visualize", env=sci, script="analysis/8_visualize.py", args=viz,
             out=P(D["v7"], "visualization.png")))
    # CASP LG
    steps.append(dict(name="casp_lg", env=sci, script="core/9_make_casp_lg.py",
             args=["--summary", casp_summary, "--out-dir", P(D["lg"]),
                   "--ligand-tsv", c["ligand_tsv"], "--target", c["target"], "--author", c["author"]],
             out=P(D["lg"], f"{c['target']}LG_model1.txt")))
    # ③ 선정 신뢰도 리포트 (맨 끝, 모든 출력 집계)
    steps.append(dict(name="confidence", env=sci, script="boost/5c_confidence.py",
             args=["--outputs", O], out=P(D["conf"], "confidence_report.md")))
    # 정리 문서 초안: 출력 CSV의 '사실'만 자동 채운 Markdown(해석/전략은 사람이 ✍️). stdlib.
    steps.append(dict(name="summary", env=std, script="analysis/11_make_summary.py",
             args=["--outputs", O, "--target", c["target"], "--task", c["task"],
                   "--ligand-tsv", c["ligand_tsv"], "--out", P(D["sum9"], "SUMMARY_DRAFT.md")],
             out=P(D["sum9"], "SUMMARY_DRAFT.md")))
    return steps


def run_step(step, scripts_dir, force, dry):
    cmd = shlex.split(step["env"]) + [os.path.join(scripts_dir, step["script"])] + step["args"]
    if (not force) and step["out"] and os.path.exists(step["out"]):
        print(f"[skip] {step['name']} -> {step['out']}"); return True
    print(f"\n[run ] {step['name']}: {' '.join(cmd)}")
    if dry:
        return True
    if subprocess.run(cmd).returncode != 0:
        print(f"[FAIL] {step['name']} — 중단 (이전 단계로 돌아가려면 --from 사용)"); return False
    if step["out"] and not os.path.exists(step["out"]):
        print(f"[WARN] {step['name']} 산출물 없음: {step['out']}")
    return True


def review(c):
    O = c["out_dir"]
    print("\n" + "=" * 64 + "\n  REVIEW (제출 전 사람 확인 / 불만족 시 --from 으로 복귀)\n" + "=" * 64)
    for f, label in [("02_pocket_candidates/pocket_candidates.csv", "포켓 후보(Step1)"),
                     ("03_pocket_validation/pocket_validation.csv", "포켓 검증(Step2)"),
                     ("05_final/selection_summary.csv", "최종 선정(Step4)")]:
        p = os.path.join(O, f)
        if os.path.exists(p):
            rows = list(csv.DictReader(open(p)))
            print(f"\n[{label}] {p}")
            for r in rows[:6]:
                print("  " + ", ".join(f"{k}={v}" for k, v in list(r.items())[:6]))
    print("\n[확인] task=", c["task"], "| AUTHOR=", c["author"], "(CHANGE-ME면 교체) | 리간드 ID/개수")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only")
    ap.add_argument("--from", dest="from_step")
    args = ap.parse_args()

    c = load_conf(args.config)
    steps = build_steps(c)
    names = [s["name"] for s in steps]
    print(f"[pipeline] target={c['target']} task={c['task']}")
    print("  스텝: " + " -> ".join(names))
    if args.only:
        steps = [s for s in steps if s["name"] == args.only] or sys.exit(f"--only: 없음 {names}")
    elif args.from_step:
        if args.from_step not in names:
            sys.exit(f"--from: 없음 {names}")
        steps = steps[names.index(args.from_step):]

    os.makedirs(c["out_dir"], exist_ok=True)
    for s in steps:
        if not run_step(s, _SELF_DIR, args.force, args.dry_run):
            sys.exit(1)
    if not args.dry_run:
        review(c)
    print("\n[pipeline] 완료")


if __name__ == "__main__":
    main()

# ── 단독 실행 (전체 파이프라인 지휘자) ── (먼저: source $CASP17/scripts/env_setup.sh)
#   SC=$CASP17/users/USERNAME/scripts
#   python3 $SC/run_pipeline.py $SC/config/T2383.conf                       # 전체 자동
#   python3 $SC/run_pipeline.py $SC/config/T2383.conf --dry-run             # 명령만 확인
#   python3 $SC/run_pipeline.py $SC/config/T2383.conf --from pocket_validate --force   # 중간부터
