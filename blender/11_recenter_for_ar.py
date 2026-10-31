# =============================================================================
# 11_recenter_for_ar.py  -  단백질 중심을 원점으로 (AR 조작 회전축 교정)
#
#  증상: AR Quick Look 에서 손가락으로 돌리면 회전축이 단백질 바깥 먼 곳에 생긴다.
#
#  원인 2가지
#   1. 카메라가 usdz 에 포함됐다. Viewpoint 6개가 (0, -279.84, 0) 에 있어서
#      씬 바운딩 박스가 y [-279.8, 51.4] 로 늘어나고 중심이 y = -114.2 로 끌렸다.
#      AR Quick Look 은 씬 바운딩 박스 중심을 조작 회전축으로 쓴다.
#   2. 단백질의 '바운딩 박스 중심'이 원점이 아니다. 오브젝트 오리진은 버텍스
#      median 에 가깝지만, median 과 bbox 중심은 다르다 (밀집 영역으로 치우친다).
#
#  이 스크립트
#   - 단백질 메쉬의 월드 바운딩 박스 중심을 구해 씬 전체를 그만큼 평행이동
#   - 카메라/라이트를 제외하고 메쉬만 선택해 내보낸다 (selected_objects_only=True)
#   - 내보낸 usda 에 Camera/Light prim 이 없는지, 중심이 맞는지 검증
#
#  ※ .blend 백업 후 실행하세요. 끝나면 freeze 스크립트를 '새 usda' 에 다시 돌리세요.
# =============================================================================
import os
import re
import bpy
from mathutils import Vector

# ---------------------------- 설정 -------------------------------------------
RECENTER_ON      = ("Fixed_Protein",)   # 기준 메쉬 패턴
RECENTER_MODE    = "scene"   # "scene": 내보낼 메쉬 전체 bbox 중심을 원점으로 (권장)
                             # "protein": 기준 메쉬 bbox 중심을 원점으로
# AR Quick Look 이 조작 회전축으로 쓰는 것은 '씬 전체' 바운딩 박스 중심이다.
# 단백질이 씬을 지배하면 둘이 거의 같지만, 터널이나 산소 궤적이 한쪽으로
# 삐져나오면 어긋난다. 그래서 기본값은 "scene" 이다.
# 스크립트가 두 중심을 모두 출력하므로 차이를 보고 바꿔도 된다.
EXCLUDE_CAMERAS  = True
EXCLUDE_LIGHTS   = True
DO_EXPORT        = True
OUT_NAME         = "o2_tunnel_centered"
# -----------------------------------------------------------------------------


def world_bbox(objs):
    """메쉬 버텍스 기준 월드 바운딩 박스 (min, max). 없으면 (None, None)."""
    lo = Vector(( 1e30,  1e30,  1e30))
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
            co = flat.reshape(n, 3)
            m = np.array([list(r) for r in mw])          # 4x4
            w = co @ m[:3, :3].T + m[:3, 3]
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


def fmt(v, n=3):
    return "(" + ", ".join(f"{c:.{n}f}" for c in v) + ")"


scn = bpy.context.scene
meshes = [o for o in bpy.data.objects if o.type == 'MESH']
targets = [o for o in meshes if any(p in o.name for p in RECENTER_ON)]

print("\n" + "=" * 92)
print(f"Blender {bpy.app.version_string}")
print(f"메쉬 {len(meshes)}개, 기준 오브젝트 {len(targets)}개"
      f" {[o.name for o in targets]}")
print("=" * 92)
if not targets:
    print(f"!! RECENTER_ON={RECENTER_ON} 에 맞는 메쉬가 없습니다. 메쉬 목록:")
    for o in meshes:
        print(f"    - {o.name}")
    raise SystemExit

# ---- 1) 기준 bbox 와 씬 bbox (이동 전) --------------------------------------
bpy.context.view_layer.update()      # matrix_world 를 최신으로
p_lo, p_hi = world_bbox(targets)
p_center = (p_lo + p_hi) / 2
s_lo, s_hi = world_bbox(meshes)
s_center = (s_lo + s_hi) / 2

print("이동 전")
print(f"  기준(단백질) bbox  min {fmt(p_lo)}  max {fmt(p_hi)}")
print(f"                     중심 {fmt(p_center)}   크기 {fmt(p_hi - p_lo)}")
print(f"  전체 메쉬 bbox     중심 {fmt(s_center)}   크기 {fmt(s_hi - s_lo)}")

cams = [o for o in bpy.data.objects if o.type == 'CAMERA']
lights = [o for o in bpy.data.objects if o.type == 'LIGHT']
if cams:
    c_lo = Vector((min(o.matrix_world.translation[i] for o in cams) for i in range(3)))
    c_hi = Vector((max(o.matrix_world.translation[i] for o in cams) for i in range(3)))
    with_cam_lo = Vector((min(s_lo[i], c_lo[i]) for i in range(3)))
    with_cam_hi = Vector((max(s_hi[i], c_hi[i]) for i in range(3)))
    print(f"  카메라 {len(cams)}개 위치 {fmt(c_lo)} ~ {fmt(c_hi)}")
    print(f"  ★ 카메라 포함 bbox 중심 {fmt((with_cam_lo + with_cam_hi) / 2)}"
          f"  <- AR 이 이 지점을 회전축으로 씀")

# ---- 2) 씬 전체 평행이동 ----------------------------------------------------
# 내보낼 메쉬만으로 bbox 를 다시 잡는다 (카메라/라이트는 어차피 빠진다)
export_meshes = meshes
e_lo, e_hi = world_bbox(export_meshes)
e_center = (e_lo + e_hi) / 2
gap = (e_center - p_center).length
print(f"  내보낼 메쉬 전체 bbox 중심 {fmt(e_center)}")
print(f"  단백질 중심과의 차이 {gap:.3f}"
      f" (단백질 반경의 {gap / max((p_hi - p_lo).length / 2, 1e-9) * 100:.1f}%)")

if RECENTER_MODE == "scene":
    delta = -e_center
    basis = "내보낼 메쉬 전체 bbox 중심 (AR 회전축과 동일)"
else:
    delta = -p_center
    basis = "단백질 bbox 중심"

print("-" * 92)
print(f"기준: {basis}")
print(f"평행이동 delta = {fmt(delta)}")

roots = [o for o in bpy.data.objects if o.parent is None]
moved = 0
for o in roots:
    loc_curves = [fc for fc in fcurves(o) if fc.data_path == "location"]
    if loc_curves:
        # location 에 애니메이션이 있으면 키 값 자체를 옮긴다
        for fc in loc_curves:
            for k in fc.keyframe_points:
                k.co[1] += delta[fc.array_index]
                k.handle_left[1] += delta[fc.array_index]
                k.handle_right[1] += delta[fc.array_index]
            fc.update()
        print(f"  {o.name}: location 키 {sum(len(f.keyframe_points) for f in loc_curves)}개 이동")
    o.location = o.location + delta
    moved += 1
print(f"  부모 없는 오브젝트 {moved}개 이동 (자식은 따라감)")

# ---- 3) 이동 후 확인 --------------------------------------------------------
# ★ location 을 바꿔도 matrix_world 는 depsgraph 가 다시 돌기 전까지 낡은 값이다.
#   이 호출이 없으면 이동 전과 똑같은 bbox 가 나와 교정이 안 된 것처럼 보인다.
bpy.context.view_layer.update()
scn.frame_set(scn.frame_current)
p_lo2, p_hi2 = world_bbox(targets)
p_center2 = (p_lo2 + p_hi2) / 2
s_lo2, s_hi2 = world_bbox(meshes)
s_center2 = (s_lo2 + s_hi2) / 2
print("-" * 92)
print("이동 후")
print(f"  단백질 bbox 중심     {fmt(p_center2)}  |{p_center2.length:.6f}|")
print(f"  전체 메쉬 bbox 중심  {fmt(s_center2)}  |{s_center2.length:.6f}|"
      f"   <- AR 회전축")
print(f"  전체 메쉬 bbox 크기  {fmt(s_hi2 - s_lo2)}")
pivot_err = s_center2.length
prot_half = (p_hi2 - p_lo2).length / 2
if pivot_err < 1e-3:
    print("  ✓ AR 회전축이 정확히 원점")
elif prot_half > 1e-9 and pivot_err / prot_half < 0.05:
    print(f"  ○ AR 회전축이 원점에서 {pivot_err:.3f}"
          f" (단백질 반경의 {pivot_err/prot_half*100:.1f}%) — 체감되지 않을 수준")
else:
    print(f"  ★ AR 회전축이 원점에서 {pivot_err:.3f}"
          f" (단백질 반경의 {pivot_err/prot_half*100:.1f}%)")
    print('    RECENTER_MODE = "scene" 으로 바꿔 다시 실행하세요.')

# ---- 4) 내보내기 (카메라/라이트 제외) ---------------------------------------
if not DO_EXPORT:
    print("=" * 92)
    raise SystemExit

folder = bpy.path.abspath("//")
if not folder:
    print("!! .blend 를 먼저 저장하세요.")
    print("=" * 92)
    raise SystemExit

# 메쉬만 선택 -> selected_objects_only 로 카메라/라이트를 확실히 배제
skip = set()
if EXCLUDE_CAMERAS:
    skip.add('CAMERA')
if EXCLUDE_LIGHTS:
    skip.add('LIGHT')

for o in bpy.data.objects:
    try:
        o.select_set(False)
    except RuntimeError:
        pass
picked = []
for o in bpy.data.objects:
    if o.type in skip:
        continue
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

print("-" * 92)
print(f"내보낼 오브젝트 {len(picked)}개 (제외: {sorted(skip) or '없음'})")
print(f"  제외된 것: 카메라 {len(cams)}개, 라이트 {len(lights)}개")

wanted = {
    "export_animation": True,
    "xform_op_mode": 'TRS',
    "selected_objects_only": True,      # ★ 선택한 메쉬만 -> 카메라/라이트 배제
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

# ---- 5) 검증 ---------------------------------------------------------------
print("-" * 92)
if not os.path.exists(usda):
    print("!! usda 가 없어 검증을 건너뜁니다.")
    print("=" * 92)
    raise SystemExit

text = open(usda, encoding="utf-8", errors="replace").read()
kinds = {}
for m in re.finditer(r'^\s*def\s+(\w+)\s+"([^"]+)"', text, re.M):
    kinds.setdefault(m.group(1), []).append(m.group(2))
print("검증: usda 안의 prim 종류")
for k in sorted(kinds):
    names = kinds[k]
    print(f"  {k:<18} {len(names):>3}개  {names[:4]}{' ...' if len(names) > 4 else ''}")

ok = True
leaked = []
for k in kinds:
    if k in ("Camera", "DistantLight", "SphereLight", "RectLight",
             "DiskLight", "CylinderLight", "DomeLight"):
        leaked.append(f"{k}({len(kinds[k])})")
if leaked:
    print(f"  ★ 카메라/라이트가 usda 에 남아 있습니다: {leaked}")
    ok = False
else:
    print("  ✓ 카메라/라이트 prim 없음 -> 바운딩 박스를 늘리지 않음")

tol = max(1e-3, prot_half * 0.05)
if pivot_err <= tol:
    print(f"  ✓ AR 회전축(씬 bbox 중심)이 원점 근처 |{pivot_err:.6f}| <= {tol:.3f}")
else:
    print(f"  ★ AR 회전축이 원점에서 {pivot_err:.3f} 떨어져 있음 (허용 {tol:.3f})")
    ok = False

print("=" * 92)
if ok:
    half = (s_hi2 - s_lo2) / 2
    print("✓ 통과. 씬 bbox 중심이 원점이고 카메라/라이트가 빠졌습니다.")
    print(f"  씬 bbox 는 원점 기준 ±{fmt(half, 1)} -> AR 회전축이 단백질 중앙에 생깁니다.")
    print()
    print("  다음: 새로 나온 usda 에 freeze 를 다시 돌리세요.")
    print(f"    python3 freeze.py {OUT_NAME}.usda")
    print(f"    -> {OUT_NAME}_frozen.usdz 를 기기로 보내 확인")
else:
    print("★ 위 경고를 확인하세요.")
print("=" * 92)
