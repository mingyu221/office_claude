# O2 터널 애니메이션 - "원점에서 날아오는" 이슈 해결

## 원인

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

## 실행 순서

| 스크립트 | 실행 위치 | 역할 |
|---|---|---|
| `01_diagnose.py` | 블렌더 Scripting 탭 | 읽기 전용. 오리진↔지오메트리 거리로 원인 확정 |
| `02_fix_origin_and_rebuild.py` | 블렌더 Scripting 탭 | 오리진 교정 + 키프레임 재생성 (**백업 필수**) |
| `03_export_usdz.py` | 블렌더 Scripting 탭 | `.usda` + `.usdz` 내보내기 (ARKit 옵션) |
| `04_verify_usda.py` | 일반 터미널 | 내보낸 `.usda` 안에 문제가 남았는지 검증 |
| `05_inject_visibility.py` | 일반 터미널 | [Plan B] USD `visibility` 로 확실하게 팝 처리 |

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

## 알려진 함정

- `OFF_SCALE` 은 **절대 0 으로 두지 말 것** (인계서 3-2: ARKit 좌표 리셋 버그). `0.001` 유지.
- 내보낼 때 `xform_op_mode='TRS'`. `'MAT'`(단일 4×4 행렬)로 나가면 RealityKit 이
  행렬을 TRS 로 분해해 보간하는데, `scale 0.001` 처럼 거의 축퇴된 행렬에서 결과가 튈 수 있다.
- `VIS_KEYS=True` 를 쓸 거면 반드시 `visible_objects_only=False` 로 내보낼 것.
  아니면 숨겨진 상태의 오브젝트가 통째로 빠진다 (인계서 3-1 과 같은 증상).
- `02` 는 되돌리기 어렵다. **반드시 .blend 사본에서 실행.**
