# =============================================================================
#  o2_all_in_one.py  -  O2 터널 애니메이션 "원점에서 날아오는" 이슈
#                       진단 / 수정 / 내보내기 / 검증 통합 스크립트
#
#  ■ 사용법 (블렌더 Text Editor 복붙용)
#     1. 블렌더 Scripting 탭 -> New -> 이 전체를 붙여넣기
#     2. 아래 MODE 를 바꿔가며 Run Script (Alt+P) 를 누른다
#          MODE = "diagnose"  ->  원인 확인 (아무것도 바꾸지 않음)
#          MODE = "fix"       ->  수정 (★ 반드시 .blend 백업 후!)
#          MODE = "export"    ->  usda + usdz 내보내기
#          MODE = "verify"    ->  내보낸 usda 검증
#     3. 결과는 Text Editor 의 텍스트 블록 드롭다운에서 "O2_REPORT" 를 고르면 보인다
#        (시스템 콘솔을 열 필요 없음). .blend 폴더에 O2_REPORT.txt 로도 저장됨.
#
#  ■ 원인 요약
#     블렌더에서 스케일은 항상 오브젝트의 '오리진'을 향해 수축한다.
#     오리진이 월드 원점에 남아 있고 버텍스만 터널 위치로 옮겨진 상태라면,
#     Location 키프레임이 하나도 없어도 scale 0.001 -> 1 구간에서
#     "원점에서 제자리로 날아오는" 것처럼 보인다.
# =============================================================================
import os
import re
import tempfile
from collections import defaultdict

import bpy
from mathutils import Vector

# ============================== 설정 =========================================
MODE = "diagnose"          # "diagnose" | "fix" | "export" | "verify"

PATTERNS   = ("O2", "OXY", "OXYGEN", "산소")   # 산소 오브젝트 이름 패턴
OFF_SCALE  = 0.001         # 숨길 때 스케일. 0 은 ARKit 좌표 리셋 버그 유발 -> 금지
UNPARENT   = True          # 부모 해제(트랜스폼 유지). parent inverse 간섭 제거
VIS_KEYS   = False         # hide_viewport/hide_render 도 키프레임 (아래 주의 참고)
OUT_NAME   = "o2_tunnel"   # 내보낼 파일 이름 (.usda / .usdz)
# =============================================================================
# VIS_KEYS = True 로 쓸 거면 반드시 export 옵션 visible_objects_only=False 여야 한다.
# 아니면 숨겨진 상태의 오브젝트가 통째로 빠져서 "일부 산소만 보이는" 증상이 재발한다.


# ----------------------------- 로그 유틸 -------------------------------------
_LOG = []


def log(*args):
    s = " ".join(str(a) for a in args)
    _LOG.append(s)
    print(s)


def flush_log():
    """Text Editor 의 'O2_REPORT' 블록과 .blend 옆 O2_REPORT.txt 에 기록."""
    body = "\n".join(_LOG)

    txt = bpy.data.texts.get("O2_REPORT")
    if txt is None:
        txt = bpy.data.texts.new("O2_REPORT")
    txt.clear()
    txt.write(body)

    folder = bpy.path.abspath("//") or tempfile.gettempdir()
    try:
        path = os.path.join(folder, "O2_REPORT.txt")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(body)
        print(f"\n[report] Text Editor 의 'O2_REPORT' 블록 또는 {path}")
    except Exception as exc:
        print(f"\n[report] 파일 저장 실패({exc}). Text Editor 의 'O2_REPORT' 블록을 보세요.")


# -------------------------- 컨텍스트 / 오퍼레이터 ----------------------------
def _find_view3d():
    """Text Editor 에서 실행해도 오퍼레이터가 돌도록 VIEW_3D 컨텍스트를 찾는다."""
    for win in bpy.context.window_manager.windows:
        for area in win.screen.areas:
            if area.type != 'VIEW_3D':
                continue
            region = next((r for r in area.regions if r.type == 'WINDOW'), None)
            if region is None:
                continue
            return {"window": win, "screen": win.screen, "area": area, "region": region}
    return None


_OVERRIDE = _find_view3d()


def op(func, **kwargs):
    """VIEW_3D 컨텍스트를 덮어쓰고 오퍼레이터 실행."""
    if _OVERRIDE and hasattr(bpy.context, "temp_override"):
        with bpy.context.temp_override(**_OVERRIDE):
            return func(**kwargs)
    return func(**kwargs)


# ------------------------------- 공통 헬퍼 -----------------------------------
def targets():
    return sorted(
        [o for o in bpy.data.objects
         if o.type == 'MESH' and any(p.upper() in o.name.upper() for p in PATTERNS)],
        key=lambda o: o.name,
    )


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


def local_median(obj):
    """로컬 좌표계 버텍스 median. 오리진이 지오메트리 중심이면 (0,0,0)."""
    if obj.type != 'MESH' or not obj.data.vertices:
        return None
    acc = Vector((0.0, 0.0, 0.0))
    for v in obj.data.vertices:
        acc += v.co
    return acc / len(obj.data.vertices)


def geo_stats(obj):
    """(로컬 median, median 기준 반경).

    오리진 이탈은 '절대 거리'가 아니라 '오브젝트 자기 크기 대비'로 판정해야 한다.
    산소 구체 반경이 0.07 인데 오리진이 6 떨어져 있으면 반경의 88배 -> 확정.
    반경의 0.5배 미만이면 오리진이 지오메트리 안에 있으므로 제자리에서 수축한다.
    """
    med = local_median(obj)
    if med is None:
        return None, 0.0
    radius = 0.0
    for v in obj.data.vertices:
        d = (v.co - med).length
        if d > radius:
            radius = d
    return med, radius


def is_identity(m, eps=1e-6):
    return all(abs(m[i][j] - (1.0 if i == j else 0.0)) < eps
               for i in range(4) for j in range(4))


def fmt(vec, n=4):
    return "(" + ", ".join(f"{c:.{n}f}" for c in vec) + ")"


def read_on_blocks(obj, f0, f1):
    """기존 scale 커브에서 '보이는' 구간 [(start, end), ...] 과 ON 스케일값 추출."""
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

    blocks, start, frame = [], None, f0
    for frame in range(f0, f1 + 1):
        if fc.evaluate(frame) > thr:
            if start is None:
                start = frame
        elif start is not None:
            blocks.append((start, frame - 1))
            start = None
    if start is not None:
        blocks.append((start, f1))
    return blocks, on_val


def build_plan(objs, f0, f1):
    """오브젝트별 ON 구간 계획. 기존 키가 없으면 이름 순서대로 균등 분할."""
    plan = {}
    for o in objs:
        blocks, on_val = read_on_blocks(o, f0, f1)
        plan[o.name] = {"blocks": blocks, "on": on_val}
        log(f"  [{o.name}] 기존 ON 구간 = {blocks}   ON 스케일 = {on_val}")

    missing = [n for n, p in plan.items() if not p["blocks"]]
    if missing:
        span = (f1 - f0 + 1) / max(len(objs), 1)
        for i, o in enumerate(objs):
            if o.name in missing:
                a = int(round(f0 + i * span))
                b = min(int(round(f0 + (i + 1) * span)) - 1, f1)
                plan[o.name] = {"blocks": [(a, b)], "on": 1.0}
                log(f"  [{o.name}] scale 키 없음 -> 균등 분할로 {plan[o.name]['blocks']} 배정")
    return plan


# =============================== diagnose ====================================
def mode_diagnose():
    scn = bpy.context.scene
    objs = targets()

    log("=" * 86)
    log(f"MODE=diagnose   Blender {bpy.app.version_string}")
    log(f"SCENE: frame {scn.frame_start}-{scn.frame_end} @ {scn.render.fps}fps"
        f"  (current {scn.frame_current})")
    log(f"대상 {len(objs)}개  패턴={PATTERNS}")

    if not objs:
        log("!! 패턴에 맞는 메쉬가 없습니다. PATTERNS 를 아래 목록 보고 수정하세요.")
        for o in bpy.data.objects:
            log(f"     - {o.name}  ({o.type})")
        log("=" * 86)
        return

    guilty, suspect = [], []
    for o in objs:
        log("-" * 86)
        log(f"[{o.name}]  mesh_users={o.data.users}")
        lm, radius = geo_stats(o)
        origin_w = o.matrix_world.translation.copy()

        log(f"  location (local)    : {fmt(o.location)}")
        log(f"  origin   (world)    : {fmt(origin_w)}")
        log(f"  scale               : {fmt(o.scale)}")

        if lm is not None:
            dist = lm.length
            ratio = dist / radius if radius > 1e-9 else 0.0
            log(f"  버텍스 median(로컬) : {fmt(lm)}")
            log(f"  오브젝트 자체 반경  : {radius:.4f}")
            log(f"  >>> origin -> geo 거리 : {dist:.4f}  (자기 반경의 {ratio:.2f}배)")
            if ratio >= 2.0:
                log(f"  >>> ★ 범인 확정. 오리진이 지오메트리 '바깥'에 있음"
                    f" (거리 {dist:.3f} > 반경 {radius:.3f}).")
                log(f"  >>>   scale {OFF_SCALE} 을 주면 지오메트리가 오리진"
                    f"(월드 {fmt(origin_w, 2)})으로 수축함.")
                guilty.append((o.name, dist, ratio))
            elif ratio >= 0.5:
                log(f"  >>> 주의. 오리진이 지오메트리 가장자리 근처 (반경의 {ratio:.2f}배).")
                suspect.append((o.name, dist, ratio))
            else:
                log("  >>> OK. 오리진이 지오메트리 안에 있음 (제자리 수축).")

        if o.delta_location.length > 1e-6 or any(abs(c - 1.0) > 1e-6 for c in o.delta_scale):
            log(f"  !! delta_location   : {fmt(o.delta_location)}")
            log(f"  !! delta_scale      : {fmt(o.delta_scale)}")

        if o.parent:
            log(f"  parent              : {o.parent.name} (type={o.parent_type})")
            if not is_identity(o.matrix_parent_inverse):
                log("  !! matrix_parent_inverse 가 단위행렬이 아님")
            if o.parent.animation_data and o.parent.animation_data.action:
                log(f"  !! 부모에 애니메이션: {o.parent.animation_data.action.name}")
        else:
            log("  parent              : None")

        if o.constraints:
            log(f"  !! constraints      : {[c.type for c in o.constraints]}")

        ad = o.animation_data
        if ad and ad.nla_tracks:
            log(f"  !! NLA tracks       : {[t.name for t in ad.nla_tracks]}")
        curves = iter_fcurves(o)
        if curves:
            log(f"  action              : {ad.action.name}")
            for fc in curves:
                interp = sorted({k.interpolation for k in fc.keyframe_points})
                frames = [int(round(k.co[0])) for k in fc.keyframe_points]
                vals = [round(k.co[1], 4) for k in fc.keyframe_points]
                flag = ""
                if fc.data_path in ("location", "rotation_euler", "rotation_quaternion",
                                    "delta_location", "delta_rotation_euler"):
                    flag = "   <<< ★ 위치/회전 커브 존재! ARKit 이 보간할 수 있음"
                log(f"    {fc.data_path}[{fc.array_index}] keys={len(frames)}"
                    f" interp={interp}{flag}")
                log(f"       frames = {frames[:14]}{' ...' if len(frames) > 14 else ''}")
                log(f"       values = {vals[:14]}{' ...' if len(vals) > 14 else ''}")
        else:
            log("  action              : None")

    log("=" * 86)
    if guilty:
        log("■ 진단: 오리진 문제입니다. 아래 오브젝트의 오리진이 지오메트리 바깥에 있습니다.")
        for name, dist, ratio in guilty:
            log(f"    - {name}: 거리 {dist:.3f} (자기 반경의 {ratio:.1f}배)")
        log('  -> .blend 를 백업하고 MODE = "fix" 로 바꿔 다시 실행하세요.')
    elif suspect:
        log("■ 오리진이 지오메트리 가장자리 근처입니다. 수정해 두는 편이 안전합니다:")
        for name, dist, ratio in suspect:
            log(f"    - {name}: 거리 {dist:.3f} (자기 반경의 {ratio:.2f}배)")
        log('  -> MODE = "fix" 를 권합니다.')
    else:
        log("■ 오리진은 정상입니다. 위에서 '★' 또는 '!!' 로 표시된 항목을 확인하세요.")
        log('  아무것도 없으면 MODE = "export" -> "verify" 로 실제 USD 내용을 보세요.')
    log("=" * 86)


# ================================= fix =======================================
def mode_fix():
    scn = bpy.context.scene
    f0, f1 = scn.frame_start, scn.frame_end
    objs = targets()

    log("=" * 86)
    log(f"MODE=fix   Blender {bpy.app.version_string}")
    log(f"대상 {len(objs)}개 / 프레임 {f0}-{f1} / OFF_SCALE={OFF_SCALE}"
        f" / UNPARENT={UNPARENT} / VIS_KEYS={VIS_KEYS}")
    if not objs:
        log("!! 패턴에 맞는 메쉬가 없습니다. PATTERNS 를 수정하세요.")
        log("=" * 86)
        return

    if bpy.context.object and bpy.context.object.mode != 'OBJECT':
        op(bpy.ops.object.mode_set, mode='OBJECT')

    # 1) 기존 타이밍을 손대기 전에 먼저 읽는다
    log("-" * 86)
    log("기존 타이밍 수집")
    plan = build_plan(objs, f0, f1)

    prefs = bpy.context.preferences.edit
    saved = prefs.keyframe_new_interpolation_type
    prefs.keyframe_new_interpolation_type = 'CONSTANT'

    try:
        # 2) 트랜스폼 정리
        for o in objs:
            log("-" * 86)
            log(f"[{o.name}] 정리")
            o.animation_data_clear()          # origin_set 이 0.001 상태에서 돌지 않게

            op(bpy.ops.object.select_all, action='DESELECT')
            o.hide_viewport = False
            o.hide_select = False
            try:
                o.hide_set(False)
            except RuntimeError:
                log("  !! 뷰 레이어에서 제외된 컬렉션에 있어 선택 불가 - 건너뜀")
                continue
            o.select_set(True)
            bpy.context.view_layer.objects.active = o

            if o.data.users > 1:
                o.data = o.data.copy()
                log("  - 메쉬 데이터 독립화")

            # 델타 흡수 + 부모 해제 (월드 위치 유지)
            mw = o.matrix_world.copy()
            o.delta_location = (0.0, 0.0, 0.0)
            o.delta_scale = (1.0, 1.0, 1.0)
            o.delta_rotation_euler = (0.0, 0.0, 0.0)
            o.delta_rotation_quaternion = (1.0, 0.0, 0.0, 0.0)
            if UNPARENT and o.parent:
                log(f"  - 부모 해제: {o.parent.name}")
                o.parent = None
                o.matrix_parent_inverse.identity()
            o.matrix_world = mw

            # ON 스케일/회전을 메쉬에 굽고 베이스를 rot=0, scale=1 로 정규화
            on_val = plan[o.name]["on"] or 1.0
            o.scale = (on_val, on_val, on_val)
            try:
                op(bpy.ops.object.transform_apply, location=False, rotation=True, scale=True)
                log(f"  - 회전/스케일({on_val}) 메쉬에 적용 -> base rot=0, scale=1")
            except RuntimeError as exc:
                log(f"  !! transform_apply 실패: {exc}")

            # ★ 핵심: 오리진을 지오메트리 중심으로
            before = local_median(o) or Vector((0.0, 0.0, 0.0))
            try:
                op(bpy.ops.object.origin_set, type='ORIGIN_GEOMETRY', center='MEDIAN')
                log(f"  - ★ 오리진 이동: 로컬 median {fmt(before, 3)} -> (0,0,0)")
                log(f"       location = {fmt(o.location, 3)}")
            except RuntimeError as exc:
                log(f"  !! origin_set 실패: {exc}")

        # 3) 키프레임 재생성
        for o in objs:
            blocks = plan[o.name]["blocks"]
            keys = {}
            for a, b in blocks:
                keys[a] = 1.0
                keys[b] = 1.0
            for a, b in blocks:                       # ON 구간 앞/뒤 1프레임에 OFF
                for f in (a - 1, b + 1):
                    if f0 <= f <= f1 and f not in keys:
                        keys[f] = OFF_SCALE
            for f in (f0, f1):                        # 양 끝 보장
                if f not in keys:
                    keys[f] = OFF_SCALE

            o.animation_data_clear()
            for f in sorted(keys):
                v = keys[f]
                o.scale = (v, v, v)
                o.keyframe_insert(data_path="scale", frame=f)
                if VIS_KEYS:
                    hidden = v <= OFF_SCALE * 1.5
                    o.hide_viewport = hidden
                    o.hide_render = hidden
                    o.keyframe_insert(data_path="hide_viewport", frame=f)
                    o.keyframe_insert(data_path="hide_render", frame=f)

            # location/rotation 커브는 애초에 만들지 않는다 -> USD 에 정적 값으로 기록
            for fc in list(iter_fcurves(o)):
                try:
                    if fc.data_path != "scale" and not fc.data_path.startswith("hide_"):
                        log(f"  !! 예상치 못한 커브 제거: {o.name}.{fc.data_path}")
                        o.animation_data.action.fcurves.remove(fc)
                    else:
                        for k in fc.keyframe_points:
                            k.interpolation = 'CONSTANT'
                        fc.update()
                except Exception:
                    pass

            log(f"  [{o.name}] 키 재생성: {sorted(keys.items())}")
    finally:
        prefs.keyframe_new_interpolation_type = saved

    # 4) 최종 확인
    log("=" * 86)
    bad = []
    for o in objs:
        lm, radius = geo_stats(o)
        lm = lm if lm is not None else Vector((0.0, 0.0, 0.0))
        tol = max(1e-4, radius * 0.01)       # 반경에 비례한 허용 오차
        paths = sorted({fc.data_path for fc in iter_fcurves(o)})
        ok = lm.length <= tol and all(p == "scale" or p.startswith("hide_") for p in paths)
        log(f"  [{o.name}] 로컬 median |{lm.length:.6f}| (허용 {tol:.6f})"
            f"  애니 채널={paths}  {'OK' if ok else '<<< 확인 필요'}")
        if not ok:
            bad.append(o.name)

    log("-" * 86)
    if bad:
        log(f"!! 아직 문제 있는 오브젝트: {bad}")
    else:
        log("■ 완료. 뷰포트에서 타임라인을 스크럽해 확인하세요:")
        log("   - 산소가 '제자리에서' 톡 나타나야 함 (원점으로 끌려가면 안 됨)")
        log("   - 각 산소를 선택해 오리진(주황 점)이 구체 안에 있는지 확인")
        log('   그 다음 MODE = "export" 로 내보내세요.')
    log("=" * 86)


# ================================ export =====================================
def mode_export():
    # xform_op_mode='TRS' 가 중요. 'MAT'(단일 4x4 행렬)로 나가면 RealityKit 이
    # 행렬을 TRS 로 분해해 보간하는데, scale 0.001 처럼 거의 축퇴된 행렬에서 튄다.
    wanted = {
        "export_animation": True,
        "xform_op_mode": 'TRS',
        "visible_objects_only": False,
        "selected_objects_only": False,
        "export_materials": True,
        "generate_preview_surface": True,
        "export_textures": True,
        "relative_paths": True,
        "export_normals": True,
        "export_uvmaps": True,
        "evaluation_mode": 'RENDER',
        "triangulate_meshes": True,
        "export_armatures": False,
        "export_shapekeys": False,
        "export_blendshapes": False,
        "export_hair": False,
        "use_instancing": False,
        "export_lights": False,
        "export_cameras": False,
    }
    props = bpy.ops.wm.usd_export.get_rna_type().properties
    supported = {k: v for k, v in wanted.items() if k in props}
    dropped = sorted(set(wanted) - set(supported))

    folder = bpy.path.abspath("//")
    if not folder:
        folder = tempfile.gettempdir()
        log("!! .blend 가 저장되지 않아 임시 폴더로 내보냅니다. 먼저 저장을 권합니다.")

    log("=" * 86)
    log(f"MODE=export   Blender {bpy.app.version_string}")
    if dropped:
        log(f"  이 버전이 지원하지 않아 생략된 옵션: {dropped}")
    for k in sorted(supported):
        log(f"  {k:26s} = {supported[k]!r}")
    log("-" * 86)

    for ext in (".usda", ".usdz"):
        path = os.path.join(folder, OUT_NAME + ext)
        try:
            res = bpy.ops.wm.usd_export(filepath=path, **supported)
            log(f"  -> {path}   {res}")
        except Exception as exc:
            log(f"  !! {path} 실패: {exc}")

    log("-" * 86)
    log('■ 다음: MODE = "verify" 로 바꿔 실행하면 내보낸 usda 를 검사합니다.')
    log("=" * 86)


# ================================ verify =====================================
def mode_verify():
    folder = bpy.path.abspath("//") or tempfile.gettempdir()
    path = os.path.join(folder, OUT_NAME + ".usda")

    log("=" * 86)
    log("MODE=verify")
    if not os.path.exists(path):
        log(f"!! 파일이 없습니다: {path}")
        log('   먼저 MODE = "export" 를 실행하세요.')
        log("=" * 86)
        return

    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        text = fh.read()

    prims = [(m.start(), m.group(1))
             for m in re.finditer(r'^\s*def\s+\w+\s+"([^"]+)"', text, re.M)]

    def owner(pos):
        name = "<stage>"
        for off, n in prims:
            if off <= pos:
                name = n
            else:
                break
        return name

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

    num = r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'
    vec3 = re.compile(rf'\(\s*({num})\s*,\s*({num})\s*,\s*({num})\s*\)')

    ts_translate = list(re.finditer(r'xformOp:translate\.timeSamples\s*=', text))
    static_translate = len(re.findall(r'xformOp:translate\s*=', text))

    moving = {}
    for m in ts_translate:
        vals = {tuple(round(float(x), 6) for x in g)
                for g in vec3.findall(block_at(m.end()))}
        if len(vals) > 1:
            srt = sorted(vals)
            moving[owner(m.start())] = (srt[0], srt[-1], len(vals))

    ts_matrix = {owner(m.start())
                 for m in re.finditer(r'xformOp:transform\.timeSamples\s*=', text)}

    zero_scale = defaultdict(list)
    for m in re.finditer(r'xformOp:scale\.timeSamples\s*=', text):
        body = block_at(m.end())
        for fm in re.finditer(rf'({num})\s*:\s*\(\s*({num})\s*,', body):
            if abs(float(fm.group(2))) < 1e-9:
                zero_scale[owner(m.start())].append(fm.group(1))

    ts_vis = {owner(m.start())
              for m in re.finditer(r'\bvisibility\.timeSamples\s*=', text)}

    log(f"FILE: {path}")
    log(f"  prim 수: {len(prims)}")
    log(f"  xformOp:translate -> 정적 {static_translate}개 / 타임샘플 {len(ts_translate)}개")
    log(f"  xformOp:scale 타임샘플: {len(re.findall(r'xformOp:scale.timeSamples', text))}개")
    log(f"  visibility 타임샘플: {len(ts_vis)}개 {sorted(ts_vis) if ts_vis else ''}")
    log("-" * 86)

    problems = []
    if moving:
        problems.append("translate 가 실제로 움직임")
        log("★ [문제] 아래 prim 의 위치가 애니메이션 중입니다. '날아오는' 직접 원인:")
        for who, (lo, hi, n) in moving.items():
            log(f"    - {who}: 서로 다른 값 {n}개, {lo} ... {hi}")
        log('    -> MODE="fix" 가 location F-Curve 를 지웠는지 다시 확인하세요.')
    elif ts_translate:
        log("○ translate 에 타임샘플은 있지만 값이 전부 같습니다 (움직이지 않음). 무해.")
    else:
        log("✓ translate 는 정적 값으로만 기록됨. 위치 보간 자체가 불가능. (최상)")

    if ts_matrix:
        problems.append("단일 4x4 행렬로 내보내짐")
        log(f"★ [문제] xformOp:transform(단일 행렬) 타임샘플: {sorted(ts_matrix)}")
        log("    -> xform_op_mode='TRS' 로 내보내세요.")
    else:
        log("✓ 단일 4x4 행렬 애니메이션 없음 (TRS 분리 상태).")

    if zero_scale:
        problems.append("scale 0 프레임 존재")
        log("★ [문제] scale 이 완전한 0 인 프레임 (ARKit 좌표 리셋 버그):")
        for who, frames in zero_scale.items():
            log(f"    - {who}: 프레임 {frames[:10]}{' ...' if len(frames) > 10 else ''}")
        log("    -> OFF_SCALE 을 0.001 같은 작은 양수로 두세요.")
    else:
        log("✓ scale 0 프레임 없음.")

    if ts_vis:
        log("✓ visibility 애니메이션 있음. token 이라 ARKit 이 보간 불가 -> 깔끔히 팝됨.")
    else:
        log("○ visibility 애니메이션 없음. scale 트릭만으로 동작합니다.")
        log("    -> 더 확실히 하려면 VIS_KEYS=True 로 fix 재실행, 또는")
        log("       05_inject_visibility.py (블렌더 밖, pip install usd-core) 사용.")

    log("-" * 86)
    if problems:
        log(f"판정: 문제 {len(problems)}건 -> {problems}")
    else:
        log("판정: 통과. 같은 설정으로 나온 .usdz 를 기기에서 확인하세요.")
    log("=" * 86)


# ================================= main ======================================
_MODES = {
    "diagnose": mode_diagnose,
    "fix": mode_fix,
    "export": mode_export,
    "verify": mode_verify,
}

try:
    handler = _MODES.get(MODE)
    if handler is None:
        log(f"!! MODE 가 잘못됐습니다: {MODE!r}. {sorted(_MODES)} 중 하나여야 합니다.")
    else:
        handler()
except Exception:
    import traceback
    log("!! 예외 발생")
    log(traceback.format_exc())
finally:
    flush_log()
