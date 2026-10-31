# =============================================================================
# 13_full_fix.py  -  09(리타이밍) + 11(중심 교정) 을 한 번에. 블렌더 재시작 불필요.
#
#  09 와 11 을 따로 돌리면, 중간에 .blend 를 저장하지 않고 블렌더를 다시 열면
#  앞 작업이 사라진다. 실제로 그렇게 돼서 centered 파일이 6프레임 2fps 로 나왔다.
#  이 스크립트는 한 번의 실행으로 전부 처리한다.
#
#  하는 일
#   1. 기존 순서를 현재 scale 커브에서 읽는다 (F2->F3->F4->F5->F6)
#   2. 씬을 240프레임 / 24fps 로 리타이밍 (10초)
#   3. hide 애니메이션 제거 (Blender 5.3 은 이걸 두면 scale 샘플이 유실된다)
#   4. 각 O2 오리진을 자기 지오메트리 중심으로
#   5. scale 키를 CONSTANT + 1프레임 전환으로 재생성
#   6. 씬 전체를 평행이동해 메쉬 bbox 중심을 원점으로 (AR 회전축)
#   7. 메쉬만 선택해 내보낸다 (카메라/라이트 배제 -> bbox 가 안 늘어남)
#   8. 결과를 직접 측정해 보고
#
#  ※ .blend 백업 후 실행. 끝나면 File > Save 로 저장하세요.
#  ※ 그 다음 터미널에서: python3 -m pip install usd-core
#                        python3 freeze.py o2_tunnel_final.usda
# =============================================================================
import os
import re
import bpy
from mathutils import Vector

# ---------------------------- 설정 -------------------------------------------
PATTERNS      = ("O2_Trajectory_F",)
RECENTER_ON   = ("Fixed_Protein",)
TARGET_FRAMES = 240
TARGET_FPS    = 24
OFF_SCALE     = 0.001
FIX_ORIGIN    = True
OUT_NAME      = "o2_tunnel_final"
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


def world_bbox(objs):
    lo = Vector((1e30, 1e30, 1e30))
    hi = Vector((-1e30, -1e30, -1e30))
    found = False
    for o in objs:
        if o.type != 'MESH' or not o.data.vertices:
            continue
        mw = o.matrix_world
        try:
            import numpy as np
            n = len(o.data.vertices)
            flat = np.empty(n * 3, dtype=np.float32)
            o.data.vertices.foreach_get("co", flat)
            m = np.array([list(r) for r in mw])
            w = flat.reshape(n, 3) @ m[:3, :3].T + m[:3, 3]
            mn, mx = w.min(axis=0), w.max(axis=0)
            for i in range(3):
                lo[i] = min(lo[i], float(mn[i]))
                hi[i] = max(hi[i], float(mx[i]))
        except Exception:
            for v in o.data.vertices:
                w = mw @ v.co
                for i in range(3):
                    lo[i] = min(lo[i], w[i])
                    hi[i] = max(hi[i], w[i])
        found = True
    return (lo, hi) if found else (None, None)


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


def fmt(v, n=3):
    return "(" + ", ".join(f"{c:.{n}f}" for c in v) + ")"


scn = bpy.context.scene
OLD0, OLD1, OLDFPS = scn.frame_start, scn.frame_end, scn.render.fps
meshes = [o for o in bpy.data.objects if o.type == 'MESH']
o2s = [o for o in meshes if any(p in o.name for p in PATTERNS)]
prot = [o for o in meshes if any(p in o.name for p in RECENTER_ON)]
cams = [o for o in bpy.data.objects if o.type == 'CAMERA']
lights = [o for o in bpy.data.objects if o.type == 'LIGHT']

print("\n" + "=" * 92)
print(f"Blender {bpy.app.version_string}")
print(f"현재: 프레임 {OLD0}-{OLD1} @ {OLDFPS}fps -> 타임코드 1칸 {1.0/OLDFPS:.3f}초")
print(f"목표: 프레임 1-{TARGET_FRAMES} @ {TARGET_FPS}fps"
      f" -> 타임코드 1칸 {1.0/TARGET_FPS:.3f}초")
print(f"메쉬 {len(meshes)}개 / O2 {len(o2s)}개 / 단백질 {len(prot)}개"
      f" / 카메라 {len(cams)}개 / 라이트 {len(lights)}개")
print("=" * 92)
if not o2s:
    print(f"!! PATTERNS={PATTERNS} 에 맞는 메쉬가 없습니다.")
    raise SystemExit
if not prot:
    print(f"!! RECENTER_ON={RECENTER_ON} 에 맞는 메쉬가 없습니다. 메쉬 목록:")
    for o in meshes:
        print(f"    - {o.name}")
    raise SystemExit

# ---- 1) 기존 순서 (프레임 범위 바꾸기 전에) ---------------------------------
order = sorted((first_on_frame(o, OLD0, OLD1) or 10 ** 6, o.name, o) for o in o2s)
print("기존 순서:")
for f, name, _ in order:
    print(f"    프레임 {f if f < 10**6 else '?':>4}  {name}")
seq = [o for _, _, o in order]

# ---- 2) 리타이밍 -----------------------------------------------------------
scn.frame_start, scn.frame_end = 1, TARGET_FRAMES
scn.render.fps, scn.render.fps_base = TARGET_FPS, 1.0
F0, F1 = 1, TARGET_FRAMES

n = len(seq)
span = TARGET_FRAMES / n
windows = {}
for i, o in enumerate(seq):
    a = int(round(1 + i * span))
    b = min(int(round(1 + (i + 1) * span)) - 1, F1)
    windows[o.name] = (a, b)
print("-" * 92)
print(f"새 구간 ({n}개 × 약 {span:.0f}프레임 = {span/TARGET_FPS:.2f}초씩):")
for o in seq:
    a, b = windows[o.name]
    print(f"    {o.name:<26} {a:>4}-{b:<4} ({(b-a+1)/TARGET_FPS:.2f}초)")

# ---- 3~5) hide 제거 + 오리진 + 키 재생성 ------------------------------------
prefs = bpy.context.preferences.edit
saved = prefs.keyframe_new_interpolation_type
prefs.keyframe_new_interpolation_type = 'CONSTANT'
try:
    if bpy.context.object and bpy.context.object.mode != 'OBJECT':
        op(bpy.ops.object.mode_set, mode='OBJECT')

    for o in seq:
        print("-" * 92)
        print(f"[{o.name}]")
        o.animation_data_clear()
        o.hide_viewport = o.hide_render = o.hide_select = False
        try:
            o.hide_set(False)
        except RuntimeError:
            pass

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
                op(bpy.ops.object.origin_set,
                   type='ORIGIN_GEOMETRY', center='MEDIAN')
                print(f"  - 오리진 이동: 로컬 median |{med.length:.3f}| -> 0")
            except RuntimeError as exc:
                print(f"  !! origin_set 실패: {exc}")

        a, b = windows[o.name]
        keys = {F0: OFF_SCALE, F1: OFF_SCALE}
        if a - 1 >= F0:
            keys[a - 1] = OFF_SCALE
        if b + 1 <= F1:
            keys[b + 1] = OFF_SCALE
        keys[a] = 1.0
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
        print(f"  - scale 키 {len(keys)}개: {sorted(keys.items())}")
        print(f"  - 애니 채널: {sorted({fc.data_path for fc in fcurves(o)})}")
finally:
    prefs.keyframe_new_interpolation_type = saved

# ---- 6) 중심 교정 ----------------------------------------------------------
bpy.context.view_layer.update()          # matrix_world 를 최신으로
scn.frame_set(1)

e_lo, e_hi = world_bbox(meshes)
e_center = (e_lo + e_hi) / 2
p_lo, p_hi = world_bbox(prot)
p_center = (p_lo + p_hi) / 2
prot_half = (p_hi - p_lo).length / 2

print("=" * 92)
print("중심 교정 전")
print(f"  단백질 bbox 중심     {fmt(p_center)}   크기 {fmt(p_hi - p_lo)}")
print(f"  내보낼 메쉬 bbox 중심 {fmt(e_center)}")
if cams:
    c_lo = Vector([min(o.matrix_world.translation[i] for o in cams)
                   for i in range(3)])
    c_hi = Vector([max(o.matrix_world.translation[i] for o in cams)
                   for i in range(3)])
    w_lo = Vector([min(e_lo[i], c_lo[i]) for i in range(3)])
    w_hi = Vector([max(e_hi[i], c_hi[i]) for i in range(3)])
    print(f"  카메라 {len(cams)}개 위치 {fmt(c_lo)}")
    print(f"  ★ 카메라까지 포함하면 중심이 {fmt((w_lo + w_hi) / 2)} 로 끌린다")

delta = -e_center
roots = [o for o in bpy.data.objects if o.parent is None]
for o in roots:
    for fc in [f for f in fcurves(o) if f.data_path == "location"]:
        for k in fc.keyframe_points:
            d = delta[fc.array_index]
            k.co[1] += d
            k.handle_left[1] += d
            k.handle_right[1] += d
        fc.update()
    o.location = o.location + delta

bpy.context.view_layer.update()          # ★ 없으면 낡은 행렬을 읽는다
scn.frame_set(1)

e_lo2, e_hi2 = world_bbox(meshes)
e_center2 = (e_lo2 + e_hi2) / 2
pivot_err = e_center2.length
print("-" * 92)
print(f"평행이동 delta = {fmt(delta)}  (부모 없는 {len(roots)}개 이동)")
print("중심 교정 후")
print(f"  내보낼 메쉬 bbox 중심 {fmt(e_center2)}  |{pivot_err:.6f}|  <- AR 회전축")
print(f"  내보낼 메쉬 bbox 크기 {fmt(e_hi2 - e_lo2)}")

# ---- 7) 내보내기 (메쉬만) ---------------------------------------------------
folder = bpy.path.abspath("//")
if not folder:
    print("!! .blend 를 먼저 저장하세요.")
    print("=" * 92)
    raise SystemExit

for o in bpy.data.objects:
    try:
        o.select_set(False)
    except RuntimeError:
        pass
picked = []
for o in bpy.data.objects:
    if o.type not in ('MESH', 'EMPTY'):
        continue
    try:
        o.hide_select = False
        o.select_set(True)
        picked.append(o.name)
    except RuntimeError:
        pass
if picked:
    bpy.context.view_layer.objects.active = bpy.data.objects[picked[0]]

wanted = {
    "export_animation": True,
    "xform_op_mode": 'TRS',
    "selected_objects_only": True,
    "export_materials": True,
    "generate_preview_surface": True,
    "relative_paths": True,
    "evaluation_mode": 'RENDER',
    "triangulate_meshes": True,
    "use_instancing": False,
    "export_cameras": False,
    "export_lights": False,
}
props = bpy.ops.wm.usd_export.get_rna_type().properties
opts = {k: v for k, v in wanted.items() if k in props}
dropped = sorted(set(wanted) - set(opts))

print("-" * 92)
print(f"내보낼 오브젝트 {len(picked)}개 (카메라 {len(cams)}개 / 라이트"
      f" {len(lights)}개 제외)")
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

# ---- 8) 측정 ---------------------------------------------------------------
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


prims = [(m.start(), m.group(1), m.group(2)) for m in
         re.finditer(r'^\s*def\s+(\w+)\s+"([^"]+)"', text, re.M)]


def owner(pos):
    name = "<stage>"
    for off, _, nm in prims:
        if off <= pos:
            name = nm
        else:
            break
    return name


def lerp(sm, f):
    ks = sorted(sm)
    if f <= ks[0]:
        return sm[ks[0]]
    if f >= ks[-1]:
        return sm[ks[-1]]
    for x, y in zip(ks, ks[1:]):
        if x <= f <= y:
            t = (f - x) / (y - x) if y > x else 0.0
            return sm[x] + (sm[y] - sm[x]) * t
    return sm[ks[-1]]


tps = float(meta("timeCodesPerSecond") or meta("framesPerSecond") or 0)
t0, t1 = float(meta("startTimeCode") or 0), float(meta("endTimeCode") or 0)
ok = True

print(f"측정: {usda}")
print(f"  startTimeCode={t0} endTimeCode={t1} timeCodesPerSecond={tps}")
if tps > 0:
    step = 1.0 / tps
    good = step <= 0.1
    ok &= good
    print(f"  {'✓' if good else '★'} 전체 {(t1-t0)/tps:.2f}초,"
          f" 타임코드 1칸 {step:.3f}초")

per_obj = {}
for m in re.finditer(r'xformOp:scale\.timeSamples\s*=', text):
    body = block_at(m.end())
    sm = {}
    for mm in re.finditer(rf'({NUM})\s*:\s*\(\s*({NUM})', body):
        sm[float(mm.group(1))] = float(mm.group(2))
    if sm:
        per_obj[owner(m.start())] = sm

print("  per-object scale:")
for name, sm in sorted(per_obj.items()):
    lo, hi = min(sm.values()), max(sm.values())
    ks = sorted(sm)
    slow = [(x, y) for x, y in zip(ks, ks[1:])
            if (sm[x] > 0.5) != (sm[y] > 0.5) and y - x != 1]
    on = [f for f in range(int(t0), int(t1) + 1) if lerp(sm, f) > 0.5]
    good = lo < 0.5 < hi and not slow and bool(on)
    ok &= good
    print(f"    {'✓' if good else '★'} {name:<26} 샘플 {len(sm):>3}개"
          f"  {lo:.4f}~{hi:.4f}  ON {min(on) if on else '-'}-{max(on) if on else '-'}")
    if slow:
        print(f"      ★ 전환이 1프레임이 아님: {slow}")

if per_obj:
    bad_cov = [f for f in range(int(t0), int(t1) + 1)
               if len([1 for sm in per_obj.values() if lerp(sm, f) > 0.5]) != 1]
    if bad_cov:
        ok = False
        print(f"    ★ 동시 표시/미표시 프레임 {len(bad_cov)}개 {bad_cov[:6]}")
    else:
        print(f"    ✓ 모든 프레임({int(t0)}-{int(t1)})에 정확히 1개만 표시")

kinds = {}
for _, k, nm in prims:
    kinds.setdefault(k, []).append(nm)
leaked = [f"{k}({len(kinds[k])})" for k in kinds
          if k in ("Camera", "DistantLight", "SphereLight", "RectLight",
                   "DiskLight", "CylinderLight", "DomeLight")]
if leaked:
    print(f"  ★ 카메라/라이트 prim 이 남아 있음: {leaked}")
    ok = False
else:
    print("  ✓ 카메라/라이트 prim 없음")

tol = max(1e-3, prot_half * 0.05)
if pivot_err <= tol:
    print(f"  ✓ AR 회전축이 원점 근처 |{pivot_err:.6f}| <= {tol:.3f}")
else:
    print(f"  ★ AR 회전축이 원점에서 {pivot_err:.3f} (허용 {tol:.3f})")
    ok = False

moving = sum(1 for m in re.finditer(r'xformOp:translate\.timeSamples\s*=', text)
             if len(set(re.findall(rf'\(\s*{NUM}\s*,\s*{NUM}\s*,\s*{NUM}\s*\)',
                                   block_at(m.end())))) > 1)
print(f"  translate 가 움직이는 prim: {moving}개")
ok &= moving == 0

print("=" * 92)
if ok:
    print("✓ 전부 통과.")
    print("  ★ 지금 바로 File > Save 로 .blend 를 저장하세요 (안 하면 또 사라집니다).")
    print("  그 다음 터미널에서:")
    print("    python3 -m pip install usd-core")
    print(f"    python3 freeze.py {OUT_NAME}.usda")
    print(f"    -> {OUT_NAME}_frozen.usdz 를 기기로 보내 확인")
else:
    print("★ 위 경고를 확인하세요.")
print("=" * 92)
