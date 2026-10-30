# =============================================================================
# 03_export_usdz.py  -  ARKit 용 USD/USDZ 내보내기 (블렌더 버전 무관 옵션 필터링)
#
#  .usda(텍스트) 와 .usdz 를 같은 이름으로 둘 다 내보냅니다.
#  .usda 는 04_verify_usda.py 로 "translate 에 타임샘플이 없는지" 검증하는 용도.
# =============================================================================
import os
import bpy

# ---------------------------- 설정 -------------------------------------------
OUT_DIR  = bpy.path.abspath("//")          # .blend 와 같은 폴더
OUT_NAME = "o2_tunnel"
# -----------------------------------------------------------------------------

# ARKit 친화 옵션.
#  xform_op_mode='TRS' 가 중요합니다. 'MAT'(단일 4x4 행렬)로 나가면 RealityKit 이
#  행렬을 TRS 로 분해해서 보간하는데, scale 0.001 같은 거의 축퇴된 행렬에서
#  분해 결과가 튈 수 있습니다. TRS 면 translate/rotate/scale 이 분리돼서 나갑니다.
WANTED = {
    "export_animation":        True,
    "xform_op_mode":           'TRS',
    "visible_objects_only":    False,   # 숨김 상태 오브젝트가 빠지는 사고 방지
    "selected_objects_only":   False,
    "export_materials":        True,
    "generate_preview_surface": True,   # ARKit 은 UsdPreviewSurface 를 읽음
    "export_textures":         True,
    "relative_paths":          True,
    "export_normals":          True,
    "export_uvmaps":           True,
    "evaluation_mode":         'RENDER',
    "triangulate_meshes":      True,
    "export_armatures":        False,
    "export_shapekeys":        False,
    "export_blendshapes":      False,
    "export_hair":             False,
    "use_instancing":          False,   # 인스턴싱은 ARKit 에서 트랜스폼 꼬임 소지
    "export_lights":           False,
    "export_cameras":          False,
}

props = bpy.ops.wm.usd_export.get_rna_type().properties
supported = {k: v for k, v in WANTED.items() if k in props}
dropped = sorted(set(WANTED) - set(supported))

print("\n" + "=" * 86)
print(f"Blender {bpy.app.version_string} / USD export")
if dropped:
    print(f"  이 버전에서 지원하지 않아 생략된 옵션: {dropped}")
for k in sorted(supported):
    print(f"  {k:26s} = {supported[k]!r}")

os.makedirs(OUT_DIR, exist_ok=True)
for ext in (".usda", ".usdz"):
    path = os.path.join(OUT_DIR, OUT_NAME + ext)
    res = bpy.ops.wm.usd_export(filepath=path, **supported)
    print(f"  -> {path}  {res}")

print("-" * 86)
print("다음: 04_verify_usda.py 로 .usda 를 검증하세요.")
print(f"  blender --background --python-expr \"\" ... 없이, 그냥 터미널에서:")
print(f"  python3 04_verify_usda.py \"{os.path.join(OUT_DIR, OUT_NAME + '.usda')}\"")
print("=" * 86 + "\n")
