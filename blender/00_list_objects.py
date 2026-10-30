# =============================================================================
# 00_list_objects.py  -  산소 오브젝트 이름부터 파악하기
#   블렌더 Scripting 탭 -> New -> 붙여넣고 Alt+P. 터미널에 출력됩니다.
#   아무것도 바꾸지 않는 읽기 전용 스크립트입니다.
# =============================================================================
import bpy
from mathutils import Vector


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


def stats(obj):
    """(버텍스 수, 지오메트리 반경, 로컬 median 크기)"""
    if obj.type != 'MESH' or not obj.data.vertices:
        return 0, 0.0, 0.0
    verts = obj.data.vertices
    med = Vector((0.0, 0.0, 0.0))
    for v in verts:
        med += v.co
    med /= len(verts)
    radius = max((v.co - med).length for v in verts)
    return len(verts), radius, med.length


rows = []
for o in bpy.data.objects:
    n_v, radius, med_off = stats(o)
    paths = {}
    for fc in fcurves(o):
        paths[fc.data_path] = paths.get(fc.data_path, 0) + len(fc.keyframe_points)
    try:
        colls = ",".join(c.name for c in o.users_collection) or "-"
    except Exception:
        colls = "-"
    loc = o.matrix_world.translation
    rows.append({
        "obj": o, "name": o.name, "type": o.type, "verts": n_v,
        "radius": radius, "med_off": med_off, "paths": paths,
        "colls": colls, "loc": (loc[0], loc[1], loc[2]),
    })

# scale 애니메이션이 있는 메쉬 = 산소 후보
cand = [r for r in rows if r["type"] == 'MESH' and "scale" in r["paths"]]

print("\n" + "=" * 118)
print(f"전체 오브젝트 {len(rows)}개   |   scale 애니메이션 있는 메쉬 {len(cand)}개")
print("=" * 118)
# 헤더는 ASCII 로. 한글은 터미널에서 2칸 폭이라 정렬이 깨진다.
print(f"{'':2}{'NAME':<30} {'TYPE':<8} {'VERTS':>7} {'RADIUS':>8} {'TRAVEL':>9} "
      f"{'WORLD POSITION':<29} {'ANIM':<22} {'COLLECTION'}")
print("   (RADIUS = 오브젝트가 퍼져 있는 반경,"
      "  TRAVEL = scale 붕괴 시 시각적 중심이 이동하는 절대 거리)")
print("-" * 118)

for r in sorted(rows, key=lambda x: (x["type"] != 'MESH', -len(x["paths"]), x["name"])):
    mark = "★" if r in cand else " "
    anim = ",".join(f"{k}:{v}" for k, v in sorted(r["paths"].items())) or "-"
    loc = "(" + ",".join(f"{c:8.2f}" for c in r["loc"]) + ")"
    off = (f"{r['med_off']:9.2f}" if r["type"] == 'MESH' and r["verts"]
           else "        -")
    print(f"{mark} {r['name']:<30} {r['type']:<8} {r['verts']:>7} {r['radius']:>8.3f} "
          f"{off:>9} {loc:<29} {anim:<22} {r['colls']}")

print("=" * 118)
if cand:
    names = [r["name"] for r in cand]
    scene_scale = max((x["radius"] for x in rows if x["type"] == 'MESH'), default=0.0)
    print(f"■ 산소 후보 {len(names)}개 (★ 표시, scale 애니메이션 보유)")
    print(f"   씬 크기 기준 = 가장 큰 메쉬 반경 {scene_scale:.2f}")
    print("-" * 118)

    moving, spread_out = [], []
    for r in cand:
        travel = r["med_off"]
        notes = []
        if scene_scale > 1e-9 and travel > scene_scale * 0.02:
            notes.append(f"★ 붕괴 시 {travel:.2f} 이동 (씬의 {travel/scene_scale*100:.1f}%)")
            moving.append(r["name"])
        else:
            notes.append("오리진 OK")
        if scene_scale > 1e-9 and r["radius"] > scene_scale * 0.10:
            notes.append(f"★ 넓게 퍼진 묶음 메쉬 (씬의 {r['radius']/scene_scale*100:.1f}%)")
            spread_out.append(r["name"])
        print(f"    - {r['name']:<28} radius={r['radius']:<8.2f} travel={travel:<8.2f}"
              f" -> {' / '.join(notes)}")

    # 공통 접두사 자동 추출
    prefix = names[0]
    for n in names[1:]:
        while prefix and not n.startswith(prefix):
            prefix = prefix[:-1]
    print("-" * 118)
    if len(prefix) >= 2:
        print(f'■ 이걸 그대로 쓰세요:   PATTERNS = ("{prefix}",)')
    else:
        print("■ 공통 접두사가 없습니다. 이름을 전부 나열하세요:")
        print("   PATTERNS = (" + ", ".join(f'"{n}"' for n in names) + ",)")

    print("-" * 118)
    if spread_out:
        print(f"■ 진단: {len(spread_out)}개가 '넓게 퍼진 묶음 메쉬' 입니다.")
        print("  분자 하나가 아니라 여러 구체가 합쳐진 덩어리이므로, 오리진을 고쳐도")
        print("  scale 수축 시 내부 구체들이 그룹 중심으로 모여듭니다. scale 방식 부적합.")
        print("  -> 07_visibility_fix.py 를 쓰세요 (scale 대신 visibility).")
    elif moving:
        print(f"■ 진단: {len(moving)}개의 오리진이 지오메트리에서 벗어나 있습니다.")
        print('  -> o2_all_in_one.py 를 MODE = "fix" 로 돌리세요.')
    else:
        print("■ 진단: 오리진 정상이고 묶음도 아닙니다. 원인이 다른 곳입니다.")
        print('  -> o2_all_in_one.py 를 MODE = "export" -> "verify" 로 USD 를 보세요.')
else:
    print("■ scale 애니메이션이 있는 메쉬가 없습니다.")
    print("  위 목록에서 산소로 보이는 이름(버텍스 적고 반경 작은 것들)을 찾아 알려주세요.")
print("=" * 118 + "\n")
