#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
L01 stage-1 결합확률 성능 검증 (정답 T/F 공개 후).

입력:
  - 제출 파일: <ID>\t<score>\t<pocket>  (기본 stage1_out/submit/L01LGCHANGE-ME.bind.txt)
  - 정답 파일: CSV  헤더 'CASP ID,canonical_smiles,binding'  (binding = TRUE/FALSE)
평가: ID로 join → binding(정답) vs score(제출).
지표: 개수/유병률, AUC-ROC(Mann-Whitney), AUC-PR(average precision),
      enrichment factor(1/5/10%), precision@K·recall@K(K=양성수, 2×양성수),
      binder vs non-binder 점수 분포, calibration(확률구간별 실제 양성비율).
순수 표준 라이브러리만 사용(numpy 불필요).
주의: 정답 파일의 canonical_smiles 열은 사용하지 않음(ID+binding만). L010357부터
      canonical_smiles가 .smi와 +1 어긋난 이슈가 있으나 binding 열은 ID 정렬로 가정.
"""
import argparse, csv, os, sys

# Windows 콘솔(cp949)에서도 한글이 깨지지 않도록 UTF-8 강제 (Linux는 무해)
try:
    sys.stdout.reconfigure(encoding="utf-8")
    sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass


def read_submission(path):
    sub = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            cid = parts[0].strip()
            try:
                score = float(parts[1]) if len(parts) > 1 and parts[1] != "" else 0.0
            except ValueError:
                score = 0.0
            sub[cid] = score
    return sub


def read_answer(path):
    """CASP ID -> bool(binding). 헤더 컬럼명 유연 처리."""
    ans = {}
    with open(path, encoding="utf-8-sig", newline="") as f:
        rdr = csv.reader(f)
        header = next(rdr)
        cols = [h.strip().lower() for h in header]
        # ID 열 / binding 열 위치 탐색
        id_i = next((i for i, c in enumerate(cols) if "id" in c), 0)
        bind_i = next((i for i, c in enumerate(cols) if "bind" in c), len(cols) - 1)
        for row in rdr:
            if not row or len(row) <= max(id_i, bind_i):
                continue
            cid = row[id_i].strip()
            b = row[bind_i].strip().upper()
            ans[cid] = (b == "TRUE" or b == "1" or b == "T")
    return ans


def auc_roc(scores, labels):
    """Mann-Whitney U 기반 AUC (tie 평균순위 보정)."""
    pairs = sorted(zip(scores, labels))
    n = len(pairs)
    ranks = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j < n and pairs[j][0] == pairs[i][0]:
            j += 1
        avg = (i + 1 + j) / 2.0  # 1-based 평균순위
        for k in range(i, j):
            ranks[k] = avg
        i = j
    npos = sum(1 for _, l in pairs if l)
    nneg = n - npos
    if npos == 0 or nneg == 0:
        return float("nan")
    sum_ranks_pos = sum(r for r, (_, l) in zip(ranks, pairs) if l)
    u = sum_ranks_pos - npos * (npos + 1) / 2.0
    return u / (npos * nneg)


def average_precision(scores, labels):
    """AUC-PR = precision-recall 곡선 아래 면적(step, tie는 함께 처리)."""
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    P = sum(1 for l in labels if l)
    if P == 0:
        return float("nan")
    ap = 0.0
    tp = 0
    prev_recall = 0.0
    i = 0
    n = len(order)
    while i < n:
        j = i
        while j < n and scores[order[j]] == scores[order[i]]:
            j += 1
        for k in range(i, j):
            if labels[order[k]]:
                tp += 1
        prec = tp / j
        recall = tp / P
        ap += prec * (recall - prev_recall)
        prev_recall = recall
        i = j
    return ap


def enrichment_factor(scores, labels, frac):
    n = len(scores)
    P = sum(1 for l in labels if l)
    if P == 0:
        return float("nan")
    k = max(1, int(round(n * frac)))
    order = sorted(range(n), key=lambda i: -scores[i])
    hits = sum(1 for i in order[:k] if labels[i])
    return (hits / k) / (P / n)


def precision_recall_at_k(scores, labels, k):
    n = len(scores)
    P = sum(1 for l in labels if l)
    order = sorted(range(n), key=lambda i: -scores[i])
    hits = sum(1 for i in order[:k] if labels[i])
    return hits / k, (hits / P if P else float("nan")), hits


def summarize(vals):
    v = sorted(vals)
    n = len(v)
    if n == 0:
        return (float("nan"),) * 3
    mean = sum(v) / n
    med = v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2
    return mean, med, (v[0], v[-1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--submission", default="stage1_out/submit/L01LGCHANGE-ME.bind.txt")
    ap.add_argument("--answer", default="", help="정답 CSV 경로 (미지정시 L01/에서 자동탐색)")
    args = ap.parse_args()

    ans_path = args.answer
    if not ans_path:
        for cand in ("stage1_answer_official.csv", "L01_answer.csv", "answer.csv",
                     "stage1_answer.csv"):
            if os.path.exists(cand):
                ans_path = cand
                break
    if not ans_path or not os.path.exists(ans_path):
        sys.exit("[에러] 정답 CSV를 찾지 못함. --answer 로 경로 지정하거나 "
                 "공식 파일을 L01/stage1_answer_official.csv 로 저장하세요.")

    sub = read_submission(args.submission)
    ans = read_answer(ans_path)
    print(f"[입력] 제출 {len(sub)}개, 정답 {len(ans)}개  (정답 파일: {ans_path})")

    # ID join (제출에 있는 ID만 평가)
    ids = [c for c in sub if c in ans]
    missing = [c for c in sub if c not in ans]
    scores = [sub[c] for c in ids]
    labels = [ans[c] for c in ids]
    P = sum(labels)
    N = len(labels) - P
    print(f"[join] 평가대상 {len(ids)}개  (정답에 없어 제외 {len(missing)}개)")
    print(f"[정답] binder(TRUE) {P}개, non-binder {N}개, 유병률 {P/len(ids)*100:.2f}%")
    if P != 80:
        print(f"  [주의] binder가 80개가 아님({P}개) -> 정답 파일/파싱 재확인 필요")

    print("\n===== 판별 성능 =====")
    print(f"AUC-ROC        : {auc_roc(scores, labels):.4f}   (0.5=무작위, 1.0=완벽)")
    print(f"AUC-PR (AP)    : {average_precision(scores, labels):.4f}   (기준선=유병률 {P/len(ids):.4f})")
    for frac in (0.01, 0.05, 0.10):
        print(f"Enrichment@{int(frac*100):>2}%  : {enrichment_factor(scores, labels, frac):.2f}x")
    for k in (P, 2 * P):
        pr, rc, hits = precision_recall_at_k(scores, labels, k)
        print(f"top-{k:<4}       : precision {pr:.3f}, recall {rc:.3f}  (binder {hits}/{P})")

    print("\n===== 점수 분포 (binder vs non-binder) =====")
    b_scores = [s for s, l in zip(scores, labels) if l]
    n_scores = [s for s, l in zip(scores, labels) if not l]
    bm, bmed, brng = summarize(b_scores)
    nm, nmed, nrng = summarize(n_scores)
    print(f"binder     : mean {bm:.3f}  median {bmed:.3f}  range [{brng[0]:.3f},{brng[1]:.3f}]")
    print(f"non-binder : mean {nm:.3f}  median {nmed:.3f}  range [{nrng[0]:.3f},{nrng[1]:.3f}]")
    print(f"평균차이(binder-nonbinder): {bm-nm:+.3f}")

    print("\n===== binder들의 순위 (내림차순 점수, 1=최고) =====")
    order = sorted(range(len(ids)), key=lambda i: -scores[i])
    rank_of = {ids[idx]: r + 1 for r, idx in enumerate(order)}
    branks = sorted(rank_of[c] for c in ids if ans[c])
    print(f"binder 순위 최고 {branks[0]}, 최저 {branks[-1]}, 중앙값 {branks[len(branks)//2]}")
    for cut in (80, 160, 240, 605):
        cnt = sum(1 for r in branks if r <= cut)
        print(f"  top-{cut:<4} 안의 binder: {cnt}/{P}")

    print("\n===== calibration (확률구간별 실제 binder 비율) =====")
    bins = [(i/10, (i+1)/10) for i in range(10)]
    for lo, hi in bins:
        grp = [l for s, l in zip(scores, labels) if (lo <= s < hi or (hi == 1.0 and s == 1.0))]
        if grp:
            print(f"  [{lo:.1f},{hi:.1f}) n={len(grp):<4} 실제 binder비율 {sum(grp)/len(grp)*100:5.1f}%")


if __name__ == "__main__":
    main()
