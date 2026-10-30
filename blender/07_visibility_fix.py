# =============================================================================
# 07_visibility_fix.py  -  scale 대신 visibility 로 팝 처리 (이 씬에 맞는 해법)
#
#  왜 scale 이 아니라 visibility 인가:
#    O2_Trajectory_Fn 은 분자 하나가 아니라 484버텍스 단위 구체 여러 개가
#    합쳐진 '궤적 묶음'이다 (반경 10.9~28.4, 터널 반경 27.3 과 맞먹음).
#    scale 0.001 은 오리진이 어디 있든 묶음 안의 구체들을 한 점으로 모은다.
#    오리진을 그룹 중심으로 옮겨도 이동이 줄어들 뿐 사라지지 않는다.
#    visibility 는 지오메트리를 전혀 건드리지 않으므로 움직임이 원리적으로 없다.
#
#  ■ 이 스크립트는 '추가'만 한다 (매우 안전)
#    - 기존 scale 애니메이션, 오리진, 부모, 지오메트리를 건드리지 않는다
#    - hide_viewport / hide_render 키프레임만 추가한다
#    - 되돌리려면 그 두 채널의 F-Curve 만 지우면 된다
#    - scale 과 visibility 가 공존해도 무해하다: 축소된 구간은 항상 '숨김'이므로
#      보이는 순간의 스케일은 늘 1 이다. pxr 폴백 경로도 그대로 살아 있다.
#
#  ■ 사용법
#    블렌더 Scripting 탭 -> New -> 붙여넣고 Alt+P (터미널에 출력됨)
# =============================================================================
import os
import re
import bpy

# ---------------------------- 설정 -------------------------------------------
PATTERNS  = ("O2_Trajectory_F",)      # 00_list_objects.py 가 알려준 값
DO_EXPORT = True                      # 끝나면 usda + usdz 자동 내보내기
OUT_NAME  = "o2_tunnel_visibility"
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


def on_blocks(obj, f0, f1):
    """기존 scale 커브에서 '보이는' 프레임 구간 [(start, end), ...] 추출."""
    curves = [fc for fc in fcurves(obj) if fc.data_path == "scale"]
    if not curves:
        return None
    fc = min(curves, key=lambda c: c.array_index)
    vals = [k.co[1] for k in fc.keyframe_points]
    if not vals or max(vals) <= 1e-9:
        return None
    thr = max(vals) * 0.5

    blocks, start = [], None
    for f in range(f0, f1 + 1):
        if fc.evaluate(f) > thr:
            if start is None:
                start = f
        elif start is not None:
            blocks.append((start, f - 1))
            start = None
    if start is not None:
        blocks.append((start, f1))
    return blocks


scn = bpy.context.scene
F0, F1 = scn.frame_start, scn.frame_end
objs = sorted([o for o in bpy.data.objects
               if o.type == 'MESH' and any(p in o.name for p in PATTERNS)],
              key=lambda o: o.name)

print("\n" + "=" * 86)
print(f"Blender {bpy.app.version_string}   프레임 {F0}-{F1}   대상 {len(objs)}개")
print("=" * 86)
if not objs:
    print("!! 패턴에 맞는 메쉬가 없습니다. PATTERNS 를 수정하세요.")
    raise SystemExit

prefs = bpy.context.preferences.edit
saved = prefs.keyframe_new_interpolation_type
prefs.keyframe_new_interpolation_type = 'CONSTANT'
added = 0

try:
    for o in objs:
        blocks = on_blocks(o, F0, F1)
        if not blocks:
            print(f"  [{o.name}] scale 커브가 없어 구간을 알 수 없음 - 건너뜀")
            continue

        # 기존 hide 커브가 있으면 먼저 제거 (재실행 가능하게)
        for fc in list(fcurves(o)):
            if fc.data_path.startswith("hide_"):
                try:
                    o.animation_data.action.fcurves.remove(fc)
                except Exception:
                    pass

        # frame -> hidden?  (True = 숨김)
        keys = {}
        for a, b in blocks:
            keys[a] = False
            if b + 1 <= F1:
                keys[b + 1] = True
        if F0 not in keys:
            keys[F0] = True

        for f in sorted(keys):
            o.hide_viewport = keys[f]
            o.hide_render = keys[f]
            o.keyframe_insert(data_path="hide_viewport", frame=f)
            o.keyframe_insert(data_path="hide_render", frame=f)

        for fc in fcurves(o):
            if fc.data_path.startswith("hide_"):
                for k in fc.keyframe_points:
                    k.interpolation = 'CONSTANT'
                fc.update()

        shown = ", ".join(f"{a}-{b}" for a, b in blocks)
        print(f"  [{o.name}] 보이는 구간 {shown}  ->  visibility 키 {len(keys)}개 추가")
        added += 1
finally:
    prefs.keyframe_new_interpolation_type = saved
    scn.frame_set(scn.frame_current)     # 현재 프레임 상태 다시 평가

print("-" * 86)
print(f"■ {added}/{len(objs)} 개에 visibility 키 추가 완료.")
print("  뷰포트에서 타임라인을 스크럽하면 산소가 켜졌다 꺼졌다 해야 합니다.")

# ----------------------------- 내보내기 --------------------------------------
if not DO_EXPORT:
    print("=" * 86)
    raise SystemExit

folder = bpy.path.abspath("//")
if not folder:
    print("!! .blend 를 먼저 저장하세요. 내보내기를 건너뜁니다.")
    print("=" * 86)
    raise SystemExit

wanted = {
    "export_animation": True,
    "xform_op_mode": 'TRS',
    "visible_objects_only": False,     # ★ 숨김 오브젝트가 빠지지 않게 반드시 False
    "selected_objects_only": False,
    "export_materials": True,
    "generate_preview_surface": True,
    "export_textures": True,
    "relative_paths": True,
    "evaluation_mode": 'RENDER',
    "triangulate_meshes": True,
    "use_instancing": False,
}
props = bpy.ops.wm.usd_export.get_rna_type().properties
opts = {k: v for k, v in wanted.items() if k in props}
dropped = sorted(set(wanted) - set(opts))

print("-" * 86)
if dropped:
    print(f"  이 버전이 지원하지 않아 생략된 옵션: {dropped}")
usda_path = os.path.join(folder, OUT_NAME + ".usda")
for ext in (".usda", ".usdz"):
    path = os.path.join(folder, OUT_NAME + ext)
    try:
        bpy.ops.wm.usd_export(filepath=path, **opts)
        print(f"  -> {path}")
    except Exception as exc:
        print(f"  !! {path} 실패: {exc}")

# ------------------------------ 검증 -----------------------------------------
print("-" * 86)
if not os.path.exists(usda_path):
    print("!! usda 가 없어 검증을 건너뜁니다.")
    print("=" * 86)
    raise SystemExit

with open(usda_path, "r", encoding="utf-8", errors="replace") as fh:
    text = fh.read()

n_vis = len(re.findall(r'\bvisibility\.timeSamples\s*=', text))
n_tr_ts = len(re.findall(r'xformOp:translate\.timeSamples\s*=', text))
n_mat_ts = len(re.findall(r'xformOp:transform\.timeSamples\s*=', text))
n_sc_ts = len(re.findall(r'xformOp:scale\.timeSamples\s*=', text))

print(f"검증: {usda_path}")
print(f"  visibility 타임샘플      : {n_vis}개")
print(f"  translate 타임샘플       : {n_tr_ts}개")
print(f"  scale 타임샘플           : {n_sc_ts}개")
print(f"  단일 4x4 행렬 타임샘플   : {n_mat_ts}개")
print("-" * 86)

if n_vis > 0:
    print("✓ 블렌더가 애니메이션 visibility 를 내보냈습니다.")
    print("  visibility 는 token 이라 ARKit 이 보간할 수 없습니다 -> 깔끔히 팝됩니다.")
    print(f"  -> {OUT_NAME}.usdz 를 기기로 보내 확인하세요.")
    if n_sc_ts:
        print("  (scale 타임샘플이 함께 있어도 무해합니다. 축소 구간은 항상 숨김이라")
        print("   보이는 순간의 스케일은 늘 1 입니다.)")
else:
    print("✗ 블렌더가 visibility 를 내보내지 않았습니다 (이 버전은 미지원).")
    print("  -> 폴백: 내보낸 usda 를 블렌더 밖에서 후처리하세요.")
    print("       pip install usd-core")
    print(f"       python3 05_inject_visibility.py '{usda_path}'")
    print("     05 가 scale 타임샘플을 읽어 visibility 로 바꾸고 scale 을 고정합니다.")
    print("     (그래서 이 스크립트는 scale 애니메이션을 일부러 지우지 않았습니다.)")

if n_mat_ts:
    print("!! 단일 4x4 행렬로 나갔습니다. xform_op_mode='TRS' 를 확인하세요.")
print("=" * 86)
