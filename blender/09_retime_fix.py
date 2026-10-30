# =============================================================================
# 09_retime_fix.py  -  '0.5초 보간'을 '1/24초 팝'으로 바꾼다 (최종 해법)
#
#  측정으로 확정된 원인:
#    timeCodesPerSecond = 2.0, 프레임 1-6  ->  타임코드 1칸 = 0.5초
#    ARKit 은 샘플 사이를 선형 보간하므로 scale 0.0001 -> 1 이 0.5초에 걸쳐
#    일어난다. 오리진이 13.7~24 단위 벗어나 있어 그 거리를 0.5초 동안 이동한다.
#    프레임 수를 줄여 '순간이동'을 만들려던 전략이 오히려 보간 시간을 늘렸다.
#
#  해법: 240프레임 / 24fps 로 되돌리고 전환을 1프레임으로 만든다.
#        타임코드 1칸 = 0.042초 -> 보간이 눈에 '팝'으로 보인다.
#
#  이 스크립트가 하는 일
#    1. 07 이 넣은 hide_viewport / hide_render 키 제거 (★ 중요)
#       Blender 5.3 은 이를 USD visibility 로 내보내지 않으면서, 숨겨진
#       프레임의 scale 샘플을 아예 쓰지 않는다. 그래서 F6 이 1샘플(1.0)만
#       남아 영구 표시됐다.
#    2. 씬을 240프레임 / 24fps 로 설정
#    3. 기존 순서(F2->F3->F4->F5->F6)를 읽어 240프레임을 균등 분배
#    4. scale 키를 CONSTANT + 1프레임 전환으로 재생성
#    5. (선택) 오리진을 각 오브젝트의 지오메트리 중심으로 교정
#    6. 내보내고, 결과 usda 의 샘플 수와 보간 시간을 직접 측정해 보고
#
#  ※ 반드시 .blend 백업 후 실행하세요.
# =============================================================================
import os
import re
import bpy
from mathutils import Vector

# ---------------------------- 설정 -------------------------------------------
PATTERNS      = ("O2_Trajectory_F",)
TARGET_FRAMES = 240       # 24fps 에서 10초
TARGET_FPS    = 24
OFF_SCALE     = 0.001     # 0 은 ARKit 좌표 리셋 버그 -> 금지
FIX_ORIGIN    = True      # 오리진을 지오메트리 중심으로 (보험. 전환이 1프레임이라 없어도 되지만)
DO_EXPORT     = True
OUT_NAME      = "o2_tunnel_retimed"
# -----------------------------------------------------------------------------


def fcurves(obj):
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


def first_on_frame(obj, f0, f1):
    """현재 scale 커브에서 '처음 보이는' 프레임. 순서 정렬에만 사용."""
    curves = [fc for fc in fcurves(obj) if fc.data_path == "scale"]
    if not curves:
        return None
    fc = min(curves, key=lambda c: c.array_index)
    vals = [k.co[1] for k in fc.keyframe_points]
    if not vals or max(vals) <= 1e-9:
        return None
    thr = max(vals) * 0.5
    for f in range(f0, f1 + 1):
        if fc.evaluate(f) > thr:
            return f
    return None


def view3d():
    for win in bpy.context.window_manager.windows:
        for area in win.screen.areas:
            if area.type != 'VIEW_3D':
                continue
            region = next((r for r in area.regions if r.type == 'WINDOW'), None)
            if region:
                return {"window": win, "screen": win.screen,
                        "area": area, "region": region}
    return None


_OV = view3d()


def op(func, **kw):
    if _OV and hasattr(bpy.context, "temp_override"):
        with bpy.context.temp_override(**_OV):
            return func(**kw)
    return func(**kw)


scn = bpy.context.scene
OLD0, OLD1, OLDFPS = scn.frame_start, scn.frame_end, scn.render.fps
objs = [o for o in bpy.data.objects
        if o.type == 'MESH' and any(p in o.name for p in PATTERNS)]

print("\n" + "=" * 92)
print(f"Blender {bpy.app.version_string}")
print(f"현재: 프레임 {OLD0}-{OLD1} @ {OLDFPS}fps"
      f"  -> 타임코드 1칸 = {1.0/OLDFPS:.3f}초")
print(f"목표: 프레임 1-{TARGET_FRAMES} @ {TARGET_FPS}fps"
      f"  -> 타임코드 1칸 = {1.0/TARGET_FPS:.3f}초")
print("=" * 92)
if not objs:
    print("!! 패턴에 맞는 메쉬가 없습니다.")
    raise SystemExit

# ---- 1) 기존 순서 읽기 (프레임 범위를 바꾸기 전에!) --------------------------
order = []
for o in objs:
    f = first_on_frame(o, OLD0, OLD1)
    order.append((f if f is not None else 10 ** 6, o.name, o))
order.sort()
print("기존 순서 (현재 scale 커브에서 읽음):")
for f, name, _ in order:
    print(f"    프레임 {f if f < 10**6 else '?':>4}  {name}")
seq = [o for _, _, o in order]

# ---- 2) 씬 리타이밍 ---------------------------------------------------------
scn.frame_start = 1
scn.frame_end = TARGET_FRAMES
scn.render.fps = TARGET_FPS
scn.render.fps_base = 1.0
F0, F1 = 1, TARGET_FRAMES

# ---- 3) 구간 균등 분배 ------------------------------------------------------
n = len(seq)
span = TARGET_FRAMES / n
windows = {}
for i, o in enumerate(seq):
    a = int(round(1 + i * span))
    b = int(round(1 + (i + 1) * span)) - 1
    windows[o.name] = (a, min(b, F1))
print("-" * 92)
print(f"새 구간 ({n}개 × 약 {span:.0f}프레임 = {span/TARGET_FPS:.2f}초씩):")
for o in seq:
    a, b = windows[o.name]
    print(f"    {o.name:<26} 프레임 {a:>4}-{b:<4}"
          f"  ({(b-a+1)/TARGET_FPS:.2f}초 동안 표시)")

# ---- 4) 정리 + 키 재생성 ----------------------------------------------------
prefs = bpy.context.preferences.edit
saved = prefs.keyframe_new_interpolation_type
prefs.keyframe_new_interpolation_type = 'CONSTANT'

try:
    if bpy.context.object and bpy.context.object.mode != 'OBJECT':
        op(bpy.ops.object.mode_set, mode='OBJECT')

    for o in seq:
        print("-" * 92)
        print(f"[{o.name}]")

        # ★ hide 키 제거 + 강제 표시. 이게 scale 샘플이 날아간 원인이었다.
        o.animation_data_clear()
        o.hide_viewport = False
        o.hide_render = False
        o.hide_select = False
        try:
            o.hide_set(False)
        except RuntimeError:
            pass
        print("  - hide 애니메이션 제거, 강제 표시")

        if FIX_ORIGIN:
            op(bpy.ops.object.select_all, action='DESELECT')
            o.select_set(True)
            bpy.context.view_layer.objects.active = o
            if o.data.users > 1:
                o.data = o.data.copy()
                print("  - 메쉬 데이터 독립화")
            o.scale = (1.0, 1.0, 1.0)
            med = Vector((0.0, 0.0, 0.0))
            for v in o.data.vertices:
                med += v.co
            med /= max(len(o.data.vertices), 1)
            try:
                op(bpy.ops.object.origin_set, type='ORIGIN_GEOMETRY', center='MEDIAN')
                print(f"  - 오리진 이동: 로컬 median |{med.length:.3f}| -> 0"
                      f"   location={tuple(round(c,2) for c in o.location)}")
            except RuntimeError as exc:
                print(f"  !! origin_set 실패: {exc}")

        a, b = windows[o.name]
        keys = {F0: OFF_SCALE, F1: OFF_SCALE}
        if a - 1 >= F0:
            keys[a - 1] = OFF_SCALE
        keys[a] = 1.0
        keys[b] = 1.0
        if b + 1 <= F1:
            keys[b + 1] = OFF_SCALE
        keys[a] = 1.0          # 경계 충돌 시 ON 우선 (앞 구간의 b+1 보다 나중에 씀)
        keys[b] = 1.0

        o.animation_data_clear()
        for f in sorted(keys):
            v = keys[f]
            o.scale = (v, v, v)
            o.keyframe_insert(data_path="scale", frame=f)

        for fc in fcurves(o):
            for k in fc.keyframe_points:
                k.interpolation = 'CONSTANT'
            fc.update()

        paths = sorted({fc.data_path for fc in fcurves(o)})
        print(f"  - scale 키 {len(keys)}개: {sorted(keys.items())}")
        print(f"  - 애니 채널: {paths}")
finally:
    prefs.keyframe_new_interpolation_type = saved
    scn.frame_set(1)

print("=" * 92)
print("■ 리타이밍 완료. 뷰포트에서 타임라인을 스크럽해 확인하세요.")
print(f"  전환이 1프레임({1.0/TARGET_FPS:.3f}초)이라 ARKit 에서 팝으로 보입니다.")

# ---- 5) 내보내기 ------------------------------------------------------------
if not DO_EXPORT:
    print("=" * 92)
    raise SystemExit

folder = bpy.path.abspath("//")
if not folder:
    print("!! .blend 를 먼저 저장하세요.")
    print("=" * 92)
    raise SystemExit

wanted = {
    "export_animation": True,
    "xform_op_mode": 'TRS',
    "selected_objects_only": False,
    "export_materials": True,
    "generate_preview_surface": True,
    "relative_paths": True,
    "evaluation_mode": 'RENDER',
    "triangulate_meshes": True,
    "use_instancing": False,
}
props = bpy.ops.wm.usd_export.get_rna_type().properties
opts = {k: v for k, v in wanted.items() if k in props}
dropped = sorted(set(wanted) - set(opts))

print("-" * 92)
if dropped:
    print(f"  이 버전이 지원하지 않아 생략된 옵션: {dropped}")
usda = os.path.join(folder, OUT_NAME + ".usda")
for ext in (".usda", ".usdz"):
    path = os.path.join(folder, OUT_NAME + ext)
    try:
        bpy.ops.wm.usd_export(filepath=path, **opts)
        print(f"  -> {path}")
    except Exception as exc:
        print(f"  !! {path} 실패: {exc}")

# ---- 6) 결과 측정 -----------------------------------------------------------
print("-" * 92)
if not os.path.exists(usda):
    print("!! usda 가 없어 측정을 건너뜁니다.")
    print("=" * 92)
    raise SystemExit

text = open(usda, encoding="utf-8", errors="replace").read()
NUM = r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'


def meta(k):
    m = re.search(rf'^\s*{k}\s*=\s*"?({NUM}|[A-Za-z]+)"?\s*$', text, re.M)
    return m.group(1) if m else None


def block_at(pos):
    i = text.find("{", pos)
    if i < 0:
        return ""
    depth, j = 0, i
    while j < len(text):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[i + 1:j]
        j += 1
    return text[i:]


prims = [(m.start(), m.group(1))
         for m in re.finditer(r'^\s*def\s+\w+\s+"([^"]+)"', text, re.M)]


def owner(pos):
    name = "<stage>"
    for off, nm in prims:
        if off <= pos:
            name = nm
        else:
            break
    return name


tps = float(meta("timeCodesPerSecond") or meta("framesPerSecond") or 0)
t0, t1 = float(meta("startTimeCode") or 0), float(meta("endTimeCode") or 0)
print(f"측정: {usda}")
print(f"  startTimeCode={t0}  endTimeCode={t1}  timeCodesPerSecond={tps}")
ok_timing = False
if tps > 0:
    step = 1.0 / tps
    print(f"  -> 전체 {(t1-t0)/tps:.2f}초,  타임코드 1칸 = {step:.3f}초")
    ok_timing = step <= 0.1
    print(f"  {'✓' if ok_timing else '★'} ARKit 보간 시간 = {step:.3f}초"
          f" {'(팝으로 보임)' if ok_timing else '(아직 눈에 보임)'}")

print("  per-object scale 샘플:")
ok_samples = True
for m in re.finditer(r'xformOp:scale\.timeSamples\s*=', text):
    body = block_at(m.end())
    frames = [float(x) for x in re.findall(rf'({NUM})\s*:\s*\(', body)]
    vals = [float(mm.group(1)) for mm in
            re.finditer(rf'{NUM}\s*:\s*\(\s*({NUM})', body)]
    lo, hi = (min(vals), max(vals)) if vals else (0, 0)
    good = len(frames) >= 4 and lo < 0.5 < hi
    ok_samples &= good
    print(f"    {'✓' if good else '★'} {owner(m.start()):<26}"
          f" 샘플 {len(frames):>3}개  {lo:.4f}~{hi:.4f}"
          f"  프레임 {[int(f) for f in frames[:8]]}")

n_vis = len(re.findall(r'\bvisibility\.timeSamples\s*=', text))
moving = 0
for m in re.finditer(r'xformOp:translate\.timeSamples\s*=', text):
    vals = set(re.findall(
        rf'\(\s*{NUM}\s*,\s*{NUM}\s*,\s*{NUM}\s*\)', block_at(m.end())))
    if len(vals) > 1:
        moving += 1
print(f"  translate 가 움직이는 prim: {moving}개")
print(f"  visibility 타임샘플: {n_vis}개")

print("=" * 92)
if ok_timing and ok_samples and moving == 0:
    print("✓ 통과. 보간 시간이 1프레임으로 줄고, 각 오브젝트가 ON/OFF 샘플을 모두")
    print("  갖고 있으며, 위치는 정적입니다.")
    print(f"  -> {OUT_NAME}.usdz 를 기기로 보내 확인하세요.")
else:
    if not ok_timing:
        print("★ 보간 시간이 아직 깁니다. TARGET_FPS 를 올리거나 TARGET_FRAMES 를 늘리세요.")
    if not ok_samples:
        print("★ 일부 오브젝트의 scale 샘플이 부족합니다 (ON 또는 OFF 가 누락).")
        print("  숨김 상태가 남아 있는지 확인하세요. 이 스크립트는 hide 키를 지웁니다.")
    if moving:
        print(f"★ translate 가 움직이는 prim 이 {moving}개 있습니다.")
print("=" * 92)
