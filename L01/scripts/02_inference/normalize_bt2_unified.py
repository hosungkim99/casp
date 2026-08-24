# -*- coding: utf-8 -*-
# normalize_bt2_unified.py — 내 boltz bt2 raw → 팀원 unified leaf 포맷으로 재배치.
#
# raw : results/bt2/seed_S/boltz_results_<id>/predictions/<id>/
#         {<id>_model_M.cif, confidence_<id>_model_M.json, pae/pde/plddt_<id>_model_M.npz}
# out : results/bt2/seed_S/sample_r/
#         { model.cif, summary.json, meta.json, arrays/{pae,pde,plddt}.npz }
#   - 5 sample을 ranking_score(=boltz confidence_score) 내림차순 정렬 → sample_0..4 재번호
#   - meta.json 에 native_sample_idx(원래 boltz model_M) 기록
#   - 검증(모든 sample model.cif 존재) 후 옛 boltz_results_<id>/ 삭제. idempotent.
#
# ⚠️ arrays/ 는 팀원의 logits(pae_logits/contact_probs/plddt_token; 추론시점 patch writer 산물)와
#    다릅니다. 나는 boltz 최종 npz(pae/pde/plddt)만 있어 그걸 넣습니다(아무도 소비 안 하므로 무해).
#
# 사용:
#   test1개(원본 보존): python normalize_bt2_unified.py --only L010016 --keep-raw
#   전체 롤아웃       : python normalize_bt2_unified.py            # $CONTAINER 전체, 옛 nesting 삭제
import os, sys, json, glob, shutil, argparse
import gemmi
try: sys.stdout.reconfigure(encoding="utf-8")
except Exception: pass

JOB_ID = os.environ.get("BT2_JOB_ID", "L01_A1")
WRITER = "post-hoc unified (normalize_bt2_unified.py)"   # 팀원 실물과 출처 구분(정직). 필요시 교체.


def cif_counts(cif):
    """n_atoms(전체 원자), n_tokens(단백질=잔기당1, 그 외=원자당1)."""
    st = gemmi.read_structure(cif); st.setup_entities()
    n_atoms = n_tokens = 0
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            natr = len(res)
            n_atoms += natr
            n_tokens += 1 if (info and info.is_amino_acid()) else natr
    return n_atoms, n_tokens


def _letters1(d):
    """{'0':v,...} → {'A':v,...} (정수 순)."""
    return {chr(65 + int(k)): (round(v, 3) if isinstance(v, (int, float)) else v)
            for k, v in sorted(d.items(), key=lambda kv: int(kv[0]))}


def _letters2(dd):
    out = {}
    for i, row in sorted(dd.items(), key=lambda kv: int(kv[0])):
        out[chr(65 + int(i))] = {chr(65 + int(j)): round(v, 3)
                                 for j, v in sorted(row.items(), key=lambda kv: int(kv[0]))}
    return out


def make_summary(conf, cif, seed, sample_idx):
    n_atoms, n_tokens = cif_counts(cif)
    cp = _letters1(conf.get("chains_ptm", {}))
    ppi = _letters2(conf.get("pair_chains_iptm", {}))
    ci = {}                                              # chain_iptm[X] = 비대각 평균(best-effort)
    for x, row in ppi.items():
        offs = [v for y, v in row.items() if y != x]
        ci[x] = round(sum(offs) / len(offs), 3) if offs else 0.0
    return {
        "model_name": "bt2",
        "seed": int(seed),
        "sample_idx": int(sample_idx),
        "n_tokens": n_tokens,
        "n_atoms": n_atoms,
        "plddt": round(conf.get("complex_plddt", 0) * 100, 3),
        "ranking_score": round(conf.get("confidence_score", 0), 3),
        "ptm": round(conf.get("ptm", 0), 3),
        "iptm": round(conf.get("iptm", 0), 3),
        "gpde": round(conf.get("complex_pde", 0), 3),
        "has_clash": None,
        "chain_ptm": cp,
        "chain_pair_iptm": ppi,
        "chain_iptm": ci,
        "inputs": {"variant": "af3_msa"},
        # 내 파이프라인 선정용 원본값 보존(팀원 스키마에 없지만 무해). collect가 우선 사용.
        "_boltz": {"ligand_iptm": round(conf.get("ligand_iptm", conf.get("iptm", 0)), 4),
                   "protein_iptm": round(conf.get("protein_iptm", 0), 4),
                   "complex_plddt": round(conf.get("complex_plddt", 0), 4)},
    }


def convert_seed(pred_dir, seed_dir, seed, binder):
    confs = sorted(glob.glob(os.path.join(pred_dir, f"confidence_{binder}_model_*.json")))
    if not confs:
        return 0
    items = []
    for cj in confs:
        m = int(cj.rsplit("_model_", 1)[1].split(".")[0])
        conf = json.load(open(cj))
        cif = os.path.join(pred_dir, f"{binder}_model_{m}.cif")
        if os.path.exists(cif):
            items.append((m, conf, cif))
    items.sort(key=lambda t: t[1].get("confidence_score", 0), reverse=True)   # ranking desc
    for rank, (m, conf, cif) in enumerate(items):
        sd = os.path.join(seed_dir, f"sample_{rank}")
        ar = os.path.join(sd, "arrays"); os.makedirs(ar, exist_ok=True)
        shutil.copy2(cif, os.path.join(sd, "model.cif"))
        with open(os.path.join(sd, "summary.json"), "w") as f:
            json.dump(make_summary(conf, cif, seed, rank), f, indent=2)
        with open(os.path.join(sd, "meta.json"), "w") as f:
            json.dump({"writer": WRITER, "model": "bt2", "job_id": JOB_ID,
                       "seed": int(seed), "sample_idx": rank, "native_sample_idx": m},
                      f, indent=2)
        for kind in ("pae", "pde", "plddt"):
            src = os.path.join(pred_dir, f"{kind}_{binder}_model_{m}.npz")
            if os.path.exists(src):
                shutil.copy2(src, os.path.join(ar, f"{kind}.npz"))
    return len(items)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--container", default=os.environ.get("CONTAINER"),
                    help="binders 컨테이너 (기본 $CONTAINER)")
    ap.add_argument("--only", nargs="*", help="특정 binder만")
    ap.add_argument("--keep-raw", action="store_true", help="옛 boltz_results_ 삭제 안 함(검증용)")
    args = ap.parse_args()
    C = args.container
    if not C or not os.path.isdir(C):
        sys.exit(f"[에러] container 경로 없음: {C} (CONTAINER env 또는 --container)")

    binders = args.only or sorted(b for b in os.listdir(C)
                                  if b.startswith("L01") and b[3:4].isdigit())
    tot = 0
    for b in binders:
        bt2 = os.path.join(C, b, "results", "bt2")
        if not os.path.isdir(bt2):
            continue
        for seed_dir in sorted(glob.glob(os.path.join(bt2, "seed_*"))):
            seed = os.path.basename(seed_dir).split("_")[1]
            pred = glob.glob(os.path.join(seed_dir, f"boltz_results_{b}", "predictions", b))
            if not pred:
                if glob.glob(os.path.join(seed_dir, "sample_0", "model.cif")):
                    continue                                    # 이미 변환됨
                print(f"  [{b} seed{seed}] raw 없음, skip"); continue
            n = convert_seed(pred[0], seed_dir, seed, b)
            ok = n > 0 and all(os.path.exists(os.path.join(seed_dir, f"sample_{r}", "model.cif"))
                               for r in range(n))
            if ok and not args.keep_raw:
                shutil.rmtree(os.path.join(seed_dir, f"boltz_results_{b}"), ignore_errors=True)
            tot += n
        print(f"[{b}] bt2 unified 변환 완료")
    print(f"총 {tot} sample 변환 -> seed_S/sample_r/{{model.cif, summary.json, meta.json, arrays/}}")


if __name__ == "__main__":
    main()
