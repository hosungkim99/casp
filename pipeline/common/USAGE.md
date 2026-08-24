# common — Usage (import 사용법)

라이브러리라 실행하지 않음. 새 스텝을 만들거나 고칠 때 **import 해서** 씀.

## 1. 부트스트랩 (스크립트 맨 위, 필수)
하위폴더 어디에 있어도 common/ 을 찾게:
```python
import sys, os
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)
sys.path.insert(0, _d)

import common.complex_io as cio
import common.geom as geom
import common.scoring as scoring
```

## 2. 자주 쓰는 패턴
```python
# 복합체 읽기 → 리간드 무게중심(포켓 좌표)
prot, ligs = cio.parse_complex("model_0.cif")
c = cio.ligand_centroid(ligs)

# 캐시된 기하 재사용 (0b 스텝이 미리 만듦)
cache = cio.load_geom_cache(os.path.join(out, "00b_geom_cache.pkl"))
g = cio.cached_geometry(path, cache)

# 대칭보정 RMSD (포즈 클러스터 거리)
perms = geom.ligand_automorphisms(smiles, n_atoms)
d = geom.sc_rmsd(P, Q, perms)

# gnina 점수 (수용체 pdb 먼저 뽑고)
cio.write_receptor_pdb("model_0.cif", "rec.pdb")
s = scoring.score_only("rec.pdb", "lig.sdf", gpu=True, gpu_id=0)
# s → CNNscore / CNNaffinity / Vina
```

## 3. 규칙
- **새 중복 함수 만들지 말 것** — 파서·정렬·gnina는 여기에만 둠(과거 5곳 중복 → 통합함).
- 새 공용 기능이 생기면 이 3파일 중 맞는 곳에 추가하고 [SCRIPTS.md](SCRIPTS.md) 갱신.
- common 자체는 numpy/gemmi 외 파이프라인 스텝을 import 하지 않음(순환 방지).
