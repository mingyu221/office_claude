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
print(f"{'':2}{'NAME':<30} {'TYPE':<8} {'VERTS':>7} {'RADIUS':>8} {'OFFSET':>9} "
      f"{'WORLD POSITION':<29} {'ANIM':<22} {'COLLECTION'}")
print("   (OFFSET = 오리진이 지오메트리에서 벗어난 정도. 자기 반경의 몇 배인지)")
print("-" * 118)

for r in sorted(rows, key=lambda x: (x["type"] != 'MESH', -len(x["paths"]), x["name"])):
    mark = "★" if r in cand else " "
    anim = ",".join(f"{k}:{v}" for k, v in sorted(r["paths"].items())) or "-"
    loc = "(" + ",".join(f"{c:8.2f}" for c in r["loc"]) + ")"
    ratio = (r["med_off"] / r["radius"]) if r["radius"] > 1e-9 else 0.0
    off = f"{ratio:8.2f}x" if r["type"] == 'MESH' and r["verts"] else "        -"
    print(f"{mark} {r['name']:<30} {r['type']:<8} {r['verts']:>7} {r['radius']:>8.3f} "
          f"{off:>9} {loc:<29} {anim:<22} {r['colls']}")

print("=" * 118)
if cand:
    names = [r["name"] for r in cand]
    print(f"■ 산소 후보 {len(names)}개 (★ 표시, scale 애니메이션 보유):")
    for r in cand:
        ratio = (r["med_off"] / r["radius"]) if r["radius"] > 1e-9 else 0.0
        verdict = ("★ 오리진 이탈 (범인)" if ratio >= 2.0
                   else "주의 (가장자리)" if ratio >= 0.5 else "오리진 OK")
        print(f"    - {r['name']:<30} radius={r['radius']:<8.3f}"
              f" offset={ratio:>8.2f}x   -> {verdict}")

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

    bad = [r for r in cand
           if r["radius"] > 1e-9 and r["med_off"] / r["radius"] >= 2.0]
    print("-" * 118)
    if bad:
        print(f"■ 진단: {len(bad)}/{len(cand)} 개의 오리진이 지오메트리 바깥입니다."
              " -> 오리진 문제 확정.")
        print('   o2_all_in_one.py 를 MODE = "fix" 로 돌리세요.')
    else:
        print("■ 진단: 오리진은 모두 정상입니다. 원인이 다른 곳입니다.")
        print('   o2_all_in_one.py 를 MODE = "export" -> "verify" 로 돌려 USD 내용을 보세요.')
else:
    print("■ scale 애니메이션이 있는 메쉬가 없습니다.")
    print("  위 목록에서 산소로 보이는 이름(버텍스 적고 반경 작은 것들)을 찾아 알려주세요.")
print("=" * 118 + "\n")
