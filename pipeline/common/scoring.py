#!/usr/bin/env python3
"""
scoring.py - gnina 점수 공용 유틸 (3·5 스텝이 공유) + 선택적 병렬 실행. import 전용.

기존엔 3_pocket_validate 와 5_final_select 가 각자 gnina 호출 코드를 중복 보유했고,
후보마다 singularity 프로세스를 순차로 띄웠다(프로세스 생성 오버헤드 + GPU 1장만 사용).
이 모듈로:
  - gnina --score_only 파싱을 한 곳으로 통일(중복 제거)
  - run_jobs 로 여러 후보를 워커 풀 병렬 실행 + GPU 라운드로빈(CUDA_VISIBLE_DEVICES)
공용 GPU 서버 배려로 병렬은 **opt-in**: jobs<=1 이면 기존과 동일한 순차 동작(라운드로빈만 적용).
"""
import os, subprocess
from concurrent.futures import ThreadPoolExecutor

GNINA_SIF = os.environ.get("GNINA_SIF", "/path/to/casp17-ligand/models/gnina/gnina.sif")


def _parse(stdout):
    """gnina stdout → {aff,cnn,cnnaff}. 콜론 유무 무관, minimize 시 minimizedAffinity 우선."""
    sc = {}
    for ln in stdout.splitlines():
        parts = ln.split()
        if not parts:
            continue
        tok = parts[0].rstrip(":")
        try:
            val = float(parts[1])
        except (IndexError, ValueError):
            continue
        if tok == "minimizedAffinity":
            sc["aff"] = val                      # 최소화 후 값(우선)
        elif tok == "Affinity" and "aff" not in sc:
            sc["aff"] = val                      # score_only 값(minimize 없을 때)
        elif tok == "CNNscore":
            sc["cnn"] = val
        elif tok == "CNNaffinity":
            sc["cnnaff"] = val
    return sc


def score_only(prot, lig, gpu=True, gpu_id=None, minimize=False):
    """gnina 점수 → {aff,cnn,cnnaff}. minimize=True면 --minimize(국소 최소화 후 값), 아니면 --score_only.
    gpu_id 주면 그 GPU로 고정."""
    env = dict(os.environ)
    if gpu and gpu_id is not None:
        # singularity --nv 는 호스트 GPU를 모두 노출 → 특정 GPU 고정은 아래 두 env 필요
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        env["SINGULARITYENV_CUDA_VISIBLE_DEVICES"] = str(gpu_id)
    # --minimize 는 최소화된 포즈를 --out 에 써야 실제 최소화가 이뤄짐(/dev/null 이면 최소화 전 값 출력).
    # → minimize 시 lig 옆에 실제 sdf 로 출력. score_only 는 포즈 안 쓰므로 /dev/null.
    mode = ["--minimize"] if minimize else ["--score_only"]
    outp = (os.path.splitext(lig)[0] + "_min.sdf") if minimize else "/dev/null"
    cmd = ["singularity", "exec"] + (["--nv"] if gpu else [])
    cmd += ["--bind", "/path/to", GNINA_SIF, "gnina", "--receptor", prot, "--ligand", lig]
    cmd += mode + ["--cnn_scoring", "rescore", "--out", outp]
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return _parse(r.stdout)


def score_with_fallback(prot, lig, gpu_id=None, minimize=False):
    """GPU 시도 → 실패(빈 결과)면 CPU 폴백. 기존 `gnina(...) or gnina(...,False)` 동작 동일."""
    return score_only(prot, lig, True, gpu_id, minimize) or score_only(prot, lig, False, None, minimize)


def parse_gpus(spec):
    """"0,1,2" → [0,1,2]. 빈 값이면 [None](= CUDA_VISIBLE_DEVICES 미설정, 기존 동작)."""
    ids = [x.strip() for x in str(spec or "").split(",") if x.strip() != ""]
    return [int(x) for x in ids] if ids else [None]


def run_jobs(fns, jobs=1, gpus=None):
    """
    fns: [callable(gpu_id) -> result]. 결과를 입력 순서대로 반환.
    jobs<=1: 순차(기존 동작, GPU 라운드로빈만). jobs>1: ThreadPoolExecutor 로 병렬
    (subprocess 는 GIL 을 놓으므로 스레드로 충분). gpus: 라운드로빈할 GPU id 리스트.
    """
    gpus = gpus or [None]
    if jobs <= 1:
        return [fn(gpus[i % len(gpus)]) for i, fn in enumerate(fns)]
    out = [None] * len(fns)
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = {ex.submit(fn, gpus[i % len(gpus)]): i for i, fn in enumerate(fns)}
        for fut in futs:
            out[futs[fut]] = fut.result()
    return out

# 이 파일은 import 전용 공용 라이브러리입니다 (단독 실행 X).
