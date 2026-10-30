# =============================================================================
# 01_diagnose.py  -  O2 터널 애니메이션 "원점에서 날아오는" 이슈 진단
#
# 사용법: 블렌더에서 o2_tunnel_video.blend 열고 Scripting 탭 -> 이 파일 열고 실행
#         결과는 Blender 콘솔(Window > Toggle System Console) 또는 Info 영역에 출력
#
# 아무것도 바꾸지 않는 읽기 전용 스크립트입니다.
# =============================================================================
import bpy
from mathutils import Vector

# --- 산소 오브젝트를 고르는 이름 패턴. 실제 이름에 맞게 수정하세요 --------------
PATTERNS = ("O2", "OXY", "OXYGEN", "산소")
# -----------------------------------------------------------------------------


def pick_targets():
    out = []
    for o in bpy.data.objects:
        up = o.name.upper()
        if any(p.upper() in up for p in PATTERNS):
            out.append(o)
    return sorted(out, key=lambda o: o.name)


def local_median(obj):
    """오브젝트 로컬 좌표계에서 버텍스 median. 오리진이 지오메트리 중심이면 (0,0,0)에 가까움."""
    if obj.type != 'MESH' or not obj.data.vertices:
        return None
    acc = Vector((0.0, 0.0, 0.0))
    for v in obj.data.vertices:
        acc += v.co
    return acc / len(obj.data.vertices)


def is_identity(m, eps=1e-6):
    return all(abs(m[i][j] - (1.0 if i == j else 0.0)) < eps
               for i in range(4) for j in range(4))


def fmt(v, n=4):
    return "(" + ", ".join(f"{c:.{n}f}" for c in v) + ")"


targets = pick_targets()
scn = bpy.context.scene

print("\n" + "=" * 86)
print(f"SCENE: frame {scn.frame_start}-{scn.frame_end} @ {scn.render.fps}fps"
      f"   (current frame: {scn.frame_current})")
print(f"대상 오브젝트 {len(targets)}개  패턴={PATTERNS}")
if not targets:
    print("!! 이름 패턴에 맞는 오브젝트가 없습니다. PATTERNS 를 수정하세요.")
    print("   씬의 전체 오브젝트 목록:")
    for o in bpy.data.objects:
        print(f"     - {o.name}  ({o.type})")

verdict = []

for o in targets:
    print("-" * 86)
    print(f"[{o.name}]  type={o.type}  mesh_users={o.data.users if o.data else '-'}")

    lm = local_median(o)
    origin_w = o.matrix_world.translation.copy()

    print(f"  location (local)     : {fmt(o.location)}")
    print(f"  origin   (world)     : {fmt(origin_w)}")
    print(f"  scale                : {fmt(o.scale)}")
    print(f"  rotation_mode        : {o.rotation_mode}")

    if lm is not None:
        off = lm.length
        print(f"  버텍스 median(로컬)  : {fmt(lm)}   |거리| = {off:.4f}")
        if off > 1e-3:
            print(f"  >>> ★ 범인 확정: 오리진이 지오메트리에서 {off:.3f} 만큼 떨어져 있음.")
            print(f"  >>>   scale 0.001 을 주면 지오메트리가 오리진(월드 {fmt(origin_w, 2)})으로 수축함.")
            verdict.append((o.name, off))
        else:
            print(f"  >>> OK: 오리진이 지오메트리 중심에 붙어 있음 (제자리 수축)")

    # 델타 트랜스폼
    dl, ds = o.delta_location, o.delta_scale
    if dl.length > 1e-6 or any(abs(c - 1.0) > 1e-6 for c in ds):
        print(f"  !! delta_location    : {fmt(dl)}")
        print(f"  !! delta_scale       : {fmt(ds)}")

    # 부모 / 부모 역행렬
    if o.parent:
        print(f"  parent               : {o.parent.name} (type={o.parent_type})")
        if not is_identity(o.matrix_parent_inverse):
            print(f"  !! matrix_parent_inverse 가 단위행렬이 아님")
        if o.parent.animation_data and o.parent.animation_data.action:
            print(f"  !! 부모에 애니메이션이 있음: {o.parent.animation_data.action.name}")
    else:
        print(f"  parent               : None")

    if o.constraints:
        print(f"  !! constraints       : {[c.type for c in o.constraints]}")
    if o.modifiers:
        print(f"  modifiers            : {[m.type for m in o.modifiers]}")

    # 애니메이션 커브
    ad = o.animation_data
    if ad and ad.nla_tracks:
        print(f"  !! NLA tracks        : {[t.name for t in ad.nla_tracks]}")
    if ad and ad.action:
        print(f"  action               : {ad.action.name}")
        for fc in ad.action.fcurves:
            interp = sorted({k.interpolation for k in fc.keyframe_points})
            frames = [int(round(k.co[0])) for k in fc.keyframe_points]
            vals = [round(k.co[1], 4) for k in fc.keyframe_points]
            flag = ""
            if fc.data_path in ("location", "rotation_euler", "rotation_quaternion",
                                "delta_location", "delta_rotation_euler"):
                flag = "   <<< ★ 위치/회전 커브가 남아 있음! ARKit 이 보간할 수 있음"
            print(f"    {fc.data_path}[{fc.array_index}]  keys={len(frames)}"
                  f"  interp={interp}{flag}")
            print(f"       frames = {frames[:14]}{' ...' if len(frames) > 14 else ''}")
            print(f"       values = {vals[:14]}{' ...' if len(vals) > 14 else ''}")
    else:
        print(f"  action               : None")

print("=" * 86)
if verdict:
    print("■ 진단 결과: 오리진 문제입니다. 아래 오브젝트들의 오리진이 지오메트리와 분리돼 있습니다.")
    for n, d in verdict:
        print(f"    - {n}: 거리 {d:.3f}")
    print("  -> 02_fix_origin_and_rebuild.py 를 실행하세요.")
else:
    print("■ 오리진은 정상입니다. 위에 '★' 로 표시된 다른 항목(위치 커브 / 부모 역행렬 /")
    print("  델타 트랜스폼 / 컨스트레인트)을 확인하고, 그래도 없으면 03_verify_usda.py 로")
    print("  실제로 내보낸 USD 안에 translate 타임샘플이 들어가는지 확인하세요.")
print("=" * 86 + "\n")
