# common — Overview

## 목적
모든 스텝(pipelines/stage1/stage2)이 공유하는 **공용 라이브러리**. 단독 실행 X, import 전용.
중복 코드를 한 곳으로 모아 파이프라인 전체가 같은 파서·기하·gnina 함수를 쓰게 한다.

## 구성
| 파일 | 역할 |
|---|---|
| `complex_io.py` | 단백질-리간드 복합체 파서 (gemmi) + 기하 캐시 + 수용체 pdb 쓰기 |
| `geom.py` | 정렬/RMSD 기하 (kabsch, SC-RMSD) |
| `scoring.py` | gnina 점수 (subprocess) + opt-in 병렬 |
| `__init__.py` | 패키지 표시(빈 파일) |

## import 방법 (중요)
스텝들이 하위폴더(pipelines/core 등)에 있어도 되게, 각 소비 스크립트 상단에 **부트스트랩**이 있음:
```python
import sys, os
_d = os.path.dirname(os.path.abspath(__file__))
while _d != os.path.dirname(_d) and not os.path.isdir(os.path.join(_d, "common")):
    _d = os.path.dirname(_d)   # common/ 가진 상위(pipeline 루트)까지 위로 탐색
sys.path.insert(0, _d)
import common.complex_io as cio
```
→ **깊이 무관**. common/ 을 가진 폴더(pipeline 루트)를 자동으로 찾아 path에 추가.

## 특징
- 타겟 종류·체인 수·리간드 개수에 **무관**하게 동작 (범용 파서).
- 함수별 상세는 [SCRIPTS.md](SCRIPTS.md), 사용 예 [USAGE.md](USAGE.md).
