# O2 터널 애니메이션 - "원점에서 날아오는" 이슈 해결

## 최종 확정된 원인 (usda 측정)

```
startTimeCode=1.0  endTimeCode=6.0  timeCodesPerSecond=2.0
-> 전체 2.500초,  타임코드 1칸 = 0.500초
xformOp:translate : 5개 prim 전부 (0,0,0) 정적 -> 진짜 이동은 없음
```

**타임코드 1칸이 0.5초다.** ARKit 은 샘플 사이를 선형 보간하므로
`scale 0.0001 -> 1` 이 0.5초에 걸쳐 일어난다. 오리진이 13.7~24 단위
벗어나 있어 그 거리를 0.5초 동안 이동한다. 이것이 "스윽 날아오는"의 정체다.

프레임 수를 6개로 줄여 "순간이동"을 만들려던 전략이 오히려 보간 시간을
늘렸다. 전환을 빠르게 하려면 프레임이 **더 많아야** 한다.

→ `09_retime_fix.py` 로 240프레임 / 24fps 로 되돌린다.
  타임코드 1칸이 0.042초가 되어 보간이 팝으로 보인다 (12배 단축).

### Blender 5.3 에서 확인된 제약

- 애니메이션 `visibility` 를 내보내지 않는다 (`visibility.timeSamples` 0개).
- `export_textures`, `visible_objects_only` 옵션이 없다.
- **숨겨진 프레임의 scale 샘플을 쓰지 않는다.** `hide_viewport` 키를 넣으면
  오히려 scale 키프레임이 유실된다. 블렌더에서 6키였던 것이 USD 에서 1~3개로
  줄고, 마지막 프레임에만 보이던 F6 은 샘플 1개(1.0)만 남아 **영구 표시**됐다.
  → `07_visibility_fix.py` 는 이 버전에서 쓰면 안 된다. `09` 가 그 키를 지운다.

## 처음 측정한 내용 (참고)

`00_list_objects.py` 로 측정한 결과:

| 오브젝트 | verts | 반경 | 붕괴 시 중심 이동 | 단백질 반경(51.4) 대비 |
|---|---|---|---|---|
| O2_Trajectory_F2 | 968 | 10.90 | 13.74 | 26.7% |
| O2_Trajectory_F3 | 2420 | 15.42 | 15.88 | 30.9% |
| O2_Trajectory_F4 | 2904 | 14.49 | 16.08 | 31.3% |
| O2_Trajectory_F5 | 3872 | 28.36 | 19.57 | 38.1% |
| O2_Trajectory_F6 | 5324 | 23.29 | 23.99 | 46.7% |

두 가지가 겹쳐 있다.

1. **오리진이 모두 월드 원점에 있다.** 붕괴 시 시각적 중심이 단백질 반경의
   27~47% 를 이동한다. 이게 보이는 "날아옴"이다.
2. **각 오브젝트가 분자 하나가 아니다.** 버텍스 수의 최대공약수가 484 이고
   각각 484 단위의 2, 5, 6, 8, 11 개 배수(합계 32개)다. 반경 10.9~28.4 는
   터널 반경 27.3 과 맞먹는다. 즉 **여러 구체가 합쳐진 궤적 묶음**이다.

따라서 **오리진을 고쳐도 해결되지 않는다.** 오리진을 그룹 중심으로 옮기면
이동이 13.7 → 최대 10.9 로 줄어들 뿐, scale 수축이 묶음 안의 구체들을 그룹
중심으로 모으는 움직임은 남는다. **이 씬에서는 scale 트릭 자체가 부적합하다.**

→ `07_visibility_fix.py` 로 `visibility` 방식을 쓴다.

## 일반론: 오리진과 스케일

블렌더에서 **스케일은 항상 오브젝트의 오리진(Origin)을 향해** 수축/팽창한다.
인계서 1절의 "단백질 바운딩 박스 중심을 씬 원점으로 정렬"을 **버텍스를 옮기는
방식**(Edit Mode 이동 / mesh data transform / 원본에서 Separate)으로 했다면,
산소 오브젝트의 **오리진은 (0,0,0)에 남고 버텍스만 터널 위치로 이동**한 상태가 된다.

이 상태에서 `scale 0.001 → 1` 을 주면

- `0.001` → 지오메트리가 **월드 원점으로 뭉개짐**
- `1` → 원점에서 제자리로 **늘어나며 날아오는 것처럼 보임**

Location 키프레임이 하나도 없어도 증상이 동일하다. ARKit 의 선형 보간은 이걸
부드럽게 만들어 더 눈에 띄게 할 뿐이고, 범인이 아니다. 그래서 `0.001` 트릭과
`CONSTANT` 보간으로는 해결되지 않았다.

## 블렌더 Text Editor 복붙 워크플로 (권장)

`o2_all_in_one.py` **한 개만** 붙여넣고 맨 위 `MODE` 만 바꿔가며 실행하면 된다.

0. 먼저 `00_list_objects.py` 를 붙여넣고 실행 → 산소 이름과 `PATTERNS` 값을 알려준다
1. Scripting 탭 → **New** → `o2_all_in_one.py` 전체 붙여넣기
2. 맨 위 `PATTERNS` 를 0번이 알려준 값으로 수정
3. `MODE` 를 바꿔가며 **Alt+P** (Run Script)

```
MODE = "diagnose"   # 원인 확인 (아무것도 바꾸지 않음)
MODE = "fix"        # 수정 (★ .blend 백업 후!)
MODE = "export"     # usda + usdz 내보내기
MODE = "verify"     # 내보낸 usda 검증
```

### ★ 중요: 결과를 어디서 보는가

Text Editor 에서 실행한 `print()` 는 **블렌더 UI 에 안 보인다.** 그래서 이 스크립트는
결과를 두 곳에 기록한다.

- Text Editor 상단 **텍스트 블록 드롭다운 → `O2_REPORT`** 선택 ← 가장 편함
- `.blend` 와 같은 폴더의 **`O2_REPORT.txt`**

시스템 콘솔을 열고 싶다면: Windows 는 `Window ▸ Toggle System Console`,
macOS/Linux 는 터미널에서 블렌더를 직접 실행해야 한다
(`/Applications/Blender.app/Contents/MacOS/Blender`).

### 복붙하지 않는 방법

파일을 로컬에 받아두면 Text Editor 의 `Text ▸ Open` 으로 바로 열 수 있다.

```bash
git clone https://github.com/mingyu221/office_claude
# 또는 이미 있으면
git pull
```

## 개별 스크립트 (분리 버전)

통합 스크립트와 기능이 같다. 터미널에서 돌리거나 단계별로 관리하고 싶을 때 사용.



| 스크립트 | 실행 위치 | 역할 |
|---|---|---|
| `00_list_objects.py` | 블렌더 Scripting 탭 | 읽기 전용. 오브젝트 이름 파악 + `PATTERNS` 자동 추출 |
| `01_diagnose.py` | 블렌더 Scripting 탭 | 읽기 전용. 오리진↔지오메트리 거리로 원인 확정 |
| `02_fix_origin_and_rebuild.py` | 블렌더 Scripting 탭 | 오리진 교정 + 키프레임 재생성 (**백업 필수**) |
| `03_export_usdz.py` | 블렌더 Scripting 탭 | `.usda` + `.usdz` 내보내기 (ARKit 옵션) |
| `04_verify_usda.py` | 일반 터미널 | 내보낸 `.usda` 안에 문제가 남았는지 검증 |
| `05_inject_visibility.py` | 일반 터미널 | USD `visibility` 주입 (07 의 폴백) |
| `07_visibility_fix.py` | 블렌더 Scripting 탭 | **이 씬의 해법.** scale 대신 visibility + 내보내기 + 검증 |

`05` 는 블렌더 안에서 돌아가지 않는다. 블렌더는 USD 를 C++ 라이브러리로만 품고 있고
파이썬 `pxr` 모듈은 노출하지 않으므로, 일반 python3 에서 실행해야 한다.

각 스크립트 상단의 `PATTERNS` 를 실제 산소 오브젝트 이름에 맞게 수정할 것.

```bash
# 04, 05 는 블렌더 없이 돌아감
pip install usd-core          # 05 에만 필요
python3 04_verify_usda.py o2_tunnel.usda
python3 05_inject_visibility.py o2_tunnel.usda
```

## 02 가 하는 일

1. 기존 `scale` 커브에서 **"보이는 구간"을 먼저 읽어 타이밍을 그대로 보존**
2. `location` / `rotation` / `delta` F-Curve **전부 제거** — ARKit 이 보간할 채널 자체를 없앰
3. 부모 해제(트랜스폼 유지) → `matrix_parent_inverse` 간섭 제거
4. delta transform 을 베이스 트랜스폼에 흡수
5. 회전·스케일을 메쉬에 적용 → 베이스를 `rot=0, scale=1` 로 정규화
6. **오리진을 지오메트리 중심으로 이동** ← 핵심
7. 키프레임을 `CONSTANT` + **1프레임 전환**으로 재생성
   (선형 보간기에서도 전환이 1/24초 안에 끝나므로 '팝'으로 보인다)

## 왜 Plan B 가 가장 확실한가

USD 의 `visibility` 는 숫자가 아니라 token(`inherited` / `invisible`)이다.
token 은 보간이 정의되지 않으므로 **RealityKit 이 중간값을 만들 수 없다.**
`FREEZE_SCALE=True` 면 scale 애니메이션까지 제거하고 `visibility` 만 남기므로,
ARKit 이 보간할 트랜스폼 채널이 **하나도 남지 않는다.**

검증된 출력 형태:

```usda
def Xform "O2_02"
{
    token visibility.timeSamples = {
        1: "invisible",
        49: "inherited",
        97: "invisible",
    }
    float3 xformOp:scale = (1, 1, 1)        # 정적
    double3 xformOp:translate = (4, -1, 2)  # 정적 - 보간 불가능
}
```

## 오리진 이탈 판정 기준

절대 거리가 아니라 **오브젝트 자기 크기 대비 비율**로 판정한다. 산소 구체 반경이
0.07 인데 오리진이 6 떨어져 있으면 반경의 약 83배이므로 확정이고, 반경의 0.5배
미만이면 오리진이 지오메트리 안에 있어 제자리에서 수축한다.

| 비율 | 판정 |
|---|---|
| 반경의 2배 이상 | ★ 범인 확정 (오리진이 지오메트리 바깥) |
| 반경의 0.5 ~ 2배 | 주의 (가장자리 근처) |
| 반경의 0.5배 미만 | OK |

## 알려진 함정

- `OFF_SCALE` 은 **절대 0 으로 두지 말 것** (인계서 3-2: ARKit 좌표 리셋 버그). `0.001` 유지.
- 내보낼 때 `xform_op_mode='TRS'`. `'MAT'`(단일 4×4 행렬)로 나가면 RealityKit 이
  행렬을 TRS 로 분해해 보간하는데, `scale 0.001` 처럼 거의 축퇴된 행렬에서 결과가 튈 수 있다.
- `VIS_KEYS=True` 를 쓸 거면 반드시 `visible_objects_only=False` 로 내보낼 것.
  아니면 숨겨진 상태의 오브젝트가 통째로 빠진다 (인계서 3-1 과 같은 증상).
- `02` 는 되돌리기 어렵다. **반드시 .blend 사본에서 실행.**
