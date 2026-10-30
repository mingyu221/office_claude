# =============================================================================
# 02_fix_origin_and_rebuild.py
#
#  "산소가 원점에서 날아오는" 이슈 수정.
#  핵심: 각 산소 오브젝트의 오리진을 자기 지오메트리 중심으로 옮긴다.
#        -> scale 0.001 이 '월드 원점'이 아니라 '제자리'에서 수축하게 됨.
#
#  추가로 하는 일:
#    - 기존 scale 키프레임에서 "보이는 구간"을 자동으로 읽어와 타이밍을 그대로 보존
#    - location / rotation / delta F-Curve 를 전부 제거 (ARKit 이 보간할 채널 자체를 없앰)
#    - 부모 해제(트랜스폼 유지) -> parent inverse matrix 간섭 제거
#    - delta transform 을 베이스 트랜스폼에 흡수
#    - 회전/스케일을 메쉬에 적용 -> 베이스 트랜스폼을 rot=0, scale=1 로 정규화
#    - 키프레임을 CONSTANT + 1프레임 전환으로 재생성 (선형 보간기에서도 '팝'으로 보임)
#
#  ※ 반드시 .blend 를 백업하거나 "Save As" 로 사본에서 돌리세요. 되돌리기 어렵습니다.
# =============================================================================
import bpy
from mathutils import Vector

# ---------------------------- 설정 -------------------------------------------
PATTERNS   = ("O2", "OXY", "OXYGEN", "산소")   # 산소 오브젝트 이름 패턴
OFF_SCALE  = 0.001      # 숨길 때의 스케일. 0 은 ARKit 좌표 리셋 버그 유발하므로 금지
UNPARENT   = True       # 부모 해제(트랜스폼 유지). parent inverse 간섭 제거용
VIS_KEYS   = False      # True 면 hide_viewport/hide_render 도 키프레임 (아래 주석 참고)
DRY_RUN    = False      # True 면 아무것도 바꾸지 않고 계획만 출력
# -----------------------------------------------------------------------------
# VIS_KEYS 주의: 블렌더 USD 익스포터가 애니메이션 visibility 를 써주면 ARKit 에서
# 가장 확실한 해법(token 이라 보간 자체가 불가능)이지만, 익스포트 옵션
# "Visible Objects Only" 가 켜져 있으면 숨겨진 상태의 오브젝트가 통째로
# 빠질 수 있습니다. 쓸 거면 반드시 visible_objects_only=False 로 내보내세요.


def iter_fcurves(obj):
    """블렌더 버전(4.4+ slotted action 포함) 상관없이 F-Curve 를 꺼낸다."""
    ad = obj.animation_data
    if not ad or not ad.action:
        return []
    act = ad.action
    try:
        legacy = list(act.fcurves)
        if legacy:
            return legacy
    except Exception:
        pass
    out = []
    try:
        for layer in getattr(act, "layers", []):
            for strip in layer.strips:
                cbag = strip.channelbag(ad.action_slot)
                if cbag:
                    out.extend(cbag.fcurves)
    except Exception:
        pass
    return out


def read_on_blocks(obj, f0, f1):
    """기존 scale 커브에서 '보이는' 프레임 구간 [(start, end), ...] 과 ON 스케일값을 추출."""
    curves = [fc for fc in iter_fcurves(obj) if fc.data_path == "scale"]
    if not curves:
        return None, None
    fc = min(curves, key=lambda c: c.array_index)
    kvals = [k.co[1] for k in fc.keyframe_points]
    if not kvals:
        return None, None
    on_val = max(kvals)
    if on_val <= 1e-9:
        return None, None
    thr = on_val * 0.5

    blocks, start = [], None
    for f in range(f0, f1 + 1):
        on = fc.evaluate(f) > thr
        if on and start is None:
            start = f
        elif not on and start is not None:
            blocks.append((start, f - 1))
            start = None
    if start is not None:
        blocks.append((start, f1))
    return blocks, on_val


def prepare_object(obj):
    """선택/활성화 가능한 상태로 만든다."""
    obj.hide_viewport = False
    obj.hide_select = False
    try:
        obj.hide_set(False)
    except RuntimeError:
        return False          # 뷰 레이어에서 제외된 컬렉션에 있음
    obj.select_set(True)
    bpy.context.view_layer.objects.active = obj
    return True


# ------------------------------- 1) 대상 수집 --------------------------------
scn = bpy.context.scene
F0, F1 = scn.frame_start, scn.frame_end

targets = sorted(
    [o for o in bpy.data.objects
     if o.type == 'MESH' and any(p.upper() in o.name.upper() for p in PATTERNS)],
    key=lambda o: o.name,
)

print("\n" + "=" * 86)
print(f"대상 {len(targets)}개 / 프레임 {F0}-{F1} / OFF_SCALE={OFF_SCALE}"
      f" / DRY_RUN={DRY_RUN}")
if not targets:
    raise SystemExit("!! 이름 패턴에 맞는 메쉬가 없습니다. PATTERNS 를 수정하세요.")

if bpy.context.object and bpy.context.object.mode != 'OBJECT':
    bpy.ops.object.mode_set(mode='OBJECT')

# ------------------- 2) 기존 타이밍 먼저 읽기 (손대기 전에!) ------------------
plan = {}
for o in targets:
    blocks, on_val = read_on_blocks(o, F0, F1)
    plan[o.name] = {"blocks": blocks, "on": on_val}
    print(f"  [{o.name}] 기존 ON 구간 = {blocks}   ON 스케일 = {on_val}")

missing = [n for n, p in plan.items() if not p["blocks"]]
if missing:
    # scale 키가 없는 오브젝트는 이름 순서대로 240프레임을 균등 분할
    n_all = len(targets)
    span = (F1 - F0 + 1) / max(n_all, 1)
    for i, o in enumerate(targets):
        if o.name in missing:
            a = int(round(F0 + i * span))
            b = int(round(F0 + (i + 1) * span)) - 1
            plan[o.name] = {"blocks": [(a, min(b, F1))], "on": 1.0}
            print(f"  [{o.name}] scale 키 없음 -> 균등 분할로 {plan[o.name]['blocks']} 배정")

if DRY_RUN:
    print("=" * 86)
    print("DRY_RUN=True 이므로 여기서 종료합니다. 계획이 맞으면 False 로 바꾸고 다시 실행하세요.")
    raise SystemExit(0)

# ------------------------- 3) 트랜스폼 정리 ----------------------------------
prefs = bpy.context.preferences.edit
saved_interp = prefs.keyframe_new_interpolation_type
prefs.keyframe_new_interpolation_type = 'CONSTANT'   # 새 키는 전부 계단식

try:
    for o in targets:
        print("-" * 86)
        print(f"[{o.name}] 정리 시작")

        # 애니메이션 전부 제거 (origin_set 이 0.001 스케일 상태에서 돌지 않도록)
        o.animation_data_clear()

        bpy.ops.object.select_all(action='DESELECT')
        if not prepare_object(o):
            print("  !! 선택 불가(제외된 컬렉션?) - 건너뜀")
            continue

        # 메쉬 데이터 독립 (공유 상태면 origin_set / transform_apply 가 실패)
        if o.data.users > 1:
            o.data = o.data.copy()
            print("  - 메쉬 데이터 독립화")

        # 델타 트랜스폼 흡수 + 부모 해제 (월드 위치는 유지)
        mw = o.matrix_world.copy()
        o.delta_location = (0.0, 0.0, 0.0)
        o.delta_scale = (1.0, 1.0, 1.0)
        o.delta_rotation_euler = (0.0, 0.0, 0.0)
        o.delta_rotation_quaternion = (1.0, 0.0, 0.0, 0.0)
        if UNPARENT and o.parent:
            print(f"  - 부모 해제: {o.parent.name}")
            o.parent = None
            o.matrix_parent_inverse.identity()
        o.matrix_world = mw

        # ON 스케일을 메쉬에 굽고, 회전도 적용 -> 베이스를 rot=0 / scale=1 로 정규화
        on_val = plan[o.name]["on"] or 1.0
        o.scale = (on_val, on_val, on_val)
        try:
            bpy.ops.object.transform_apply(location=False, rotation=True, scale=True)
            print(f"  - 회전/스케일({on_val}) 메쉬에 적용 -> base rot=0, scale=1")
        except RuntimeError as e:
            print(f"  !! transform_apply 실패: {e}")

        # ★ 핵심: 오리진을 지오메트리 중심으로
        before = Vector((0.0, 0.0, 0.0))
        for v in o.data.vertices:
            before += v.co
        before /= max(len(o.data.vertices), 1)
        try:
            bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')
            print(f"  - ★ 오리진 이동: 로컬 median {tuple(round(c,3) for c in before)}"
                  f" -> (0,0,0),  location = {tuple(round(c,3) for c in o.location)}")
        except RuntimeError as e:
            print(f"  !! origin_set 실패: {e}")

    # ------------------------ 4) 키프레임 재생성 -----------------------------
    for o in targets:
        blocks = plan[o.name]["blocks"]
        off = OFF_SCALE

        keys = {}
        for a, b in blocks:
            keys[a] = 1.0
            keys[b] = 1.0
        for a, b in blocks:                      # ON 구간 바로 앞/뒤 1프레임에 OFF
            for f in (a - 1, b + 1):
                if F0 <= f <= F1 and f not in keys:
                    keys[f] = off
        for f in (F0, F1):                       # 양 끝 보장
            if f not in keys:
                keys[f] = off

        o.animation_data_clear()
        for f in sorted(keys):
            v = keys[f]
            o.scale = (v, v, v)
            o.keyframe_insert(data_path="scale", frame=f)
            if VIS_KEYS:
                hidden = (v <= off * 1.5)
                o.hide_viewport = hidden
                o.hide_render = hidden
                o.keyframe_insert(data_path="hide_viewport", frame=f)
                o.keyframe_insert(data_path="hide_render", frame=f)

        # location / rotation 커브는 '만들지 않는다' -> USD 에 정적 값으로만 기록됨
        for fc in iter_fcurves(o):
            try:
                if fc.data_path != "scale" and not fc.data_path.startswith("hide_"):
                    print(f"  !! 예상치 못한 커브 제거: {o.name}.{fc.data_path}")
                    o.animation_data.action.fcurves.remove(fc)
                else:
                    for k in fc.keyframe_points:
                        k.interpolation = 'CONSTANT'
                    fc.update()
            except Exception:
                pass

        print(f"  [{o.name}] 키 재생성: {sorted(keys.items())}")
finally:
    prefs.keyframe_new_interpolation_type = saved_interp

# --------------------------- 5) 최종 확인 ------------------------------------
print("=" * 86)
bad = []
for o in targets:
    lm = Vector((0.0, 0.0, 0.0))
    for v in o.data.vertices:
        lm += v.co
    lm /= max(len(o.data.vertices), 1)
    radius = max((v.co - lm).length for v in o.data.vertices) if o.data.vertices else 0.0
    tol = max(1e-4, radius * 0.01)       # 반경에 비례한 허용 오차
    paths = sorted({fc.data_path for fc in iter_fcurves(o)})
    ok = lm.length <= tol and all(p == "scale" or p.startswith("hide_") for p in paths)
    print(f"  [{o.name}] 로컬 median |{lm.length:.6f}| (허용 {tol:.6f})"
          f"   애니 채널={paths}   {'OK' if ok else '<<< 확인 필요'}")
    if not ok:
        bad.append(o.name)

print("-" * 86)
if bad:
    print(f"!! 아직 문제 있는 오브젝트: {bad}")
else:
    print("■ 완료. 뷰포트에서 타임라인을 스크럽해 보세요:")
    print("   - 산소가 '제자리에서' 톡 나타나야 합니다 (원점으로 끌려가면 안 됨)")
    print("   - 각 산소를 선택해 오리진(주황 점)이 산소 구체 안에 있는지 확인")
    print("   그 다음 03_export_usdz.py 로 내보내고 04_verify_usda.py 로 검증하세요.")
print("=" * 86 + "\n")
