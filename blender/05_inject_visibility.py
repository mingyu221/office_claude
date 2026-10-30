#!/usr/bin/env python3
# =============================================================================
# 05_inject_visibility.py  -  [Plan B] USD 의 visibility 로 확실하게 팝 처리
#
#  왜 이게 가장 확실한가:
#    USD 의 visibility 는 숫자가 아니라 token("inherited"/"invisible") 입니다.
#    token 은 보간 자체가 정의되지 않으므로 ARKit/RealityKit 이 중간값을 만들 수
#    없습니다. 즉 "서서히 나타남 / 날아옴"이 원리적으로 불가능합니다.
#    FREEZE_SCALE=True 로 두면 scale 애니메이션도 통째로 없애고 visibility 만
#    남기므로, ARKit 이 보간할 트랜스폼 채널이 하나도 남지 않습니다.
#
#  필요한 것:  pip install usd-core        (블렌더 내장 파이썬이 아니라 일반 python3)
#  사용법:     python3 05_inject_visibility.py o2_tunnel.usda [출력.usdz]
# =============================================================================
import os
import sys

try:
    from pxr import Usd, UsdGeom, UsdUtils, Sdf
except ImportError:
    sys.exit("pxr 모듈이 없습니다.  pip install usd-core  후 다시 실행하세요.")

# ---------------------------- 설정 -------------------------------------------
PATTERNS     = ("O2", "OXY", "OXYGEN")
FREEZE_SCALE = True    # True: scale 애니메이션을 제거하고 ON 값으로 고정, visibility 만 사용
ON_THRESHOLD = 0.5     # scale 이 (최대값 * 이 비율)보다 크면 '보이는' 프레임으로 판정
# -----------------------------------------------------------------------------

if len(sys.argv) < 2:
    sys.exit("usage: python3 05_inject_visibility.py <in.usda|in.usdz> [out.usdz]")

src = sys.argv[1]
out_usdz = sys.argv[2] if len(sys.argv) > 2 else None

stage = Usd.Stage.Open(src)
if not stage:
    sys.exit(f"열 수 없습니다: {src}")

start = stage.GetStartTimeCode()
end = stage.GetEndTimeCode()
print("=" * 86)
print(f"IN : {src}")
print(f"     timeCode {start} .. {end}   upAxis={UsdGeom.GetStageUpAxis(stage)}")
print(f"     FREEZE_SCALE={FREEZE_SCALE}")
print("-" * 86)

touched = 0
for prim in stage.Traverse():
    name = prim.GetName().upper()
    if not any(p in name for p in PATTERNS):
        continue

    xf = UsdGeom.Xformable(prim)
    if not xf:
        continue

    scale_op = None
    for op in xf.GetOrderedXformOps():
        if op.GetOpType() == UsdGeom.XformOp.TypeScale:
            scale_op = op
            break
    if scale_op is None:
        print(f"  [{prim.GetPath()}] scale op 없음 - 건너뜀")
        continue

    times = scale_op.GetTimeSamples()
    if not times:
        print(f"  [{prim.GetPath()}] scale 타임샘플 없음 - 건너뜀")
        continue

    samples = [(t, scale_op.Get(t)) for t in times]
    on_val = max(max(abs(c) for c in v) for _, v in samples)
    thr = on_val * ON_THRESHOLD

    imageable = UsdGeom.Imageable(prim)
    vis_attr = imageable.CreateVisibilityAttr()
    vis_attr.Clear()

    # 첫 샘플 전 구간은 숨김으로 시작
    prev = None
    written = []
    if times[0] > start:
        vis_attr.Set(UsdGeom.Tokens.invisible, start)
        written.append((start, "invisible"))
        prev = False

    for t, v in samples:
        on = max(abs(c) for c in v) > thr
        if on != prev:
            tok = UsdGeom.Tokens.inherited if on else UsdGeom.Tokens.invisible
            vis_attr.Set(tok, t)
            written.append((t, "inherited" if on else "invisible"))
            prev = on

    if FREEZE_SCALE:
        # scale 애니메이션 제거 -> ON 값으로 고정. 보간할 채널이 사라진다.
        attr = scale_op.GetAttr()
        for t in times:
            attr.ClearAtTime(t)
        full = next(v for _, v in samples if max(abs(c) for c in v) > thr)
        attr.Set(full)

    print(f"  [{prim.GetPath()}] ON 스케일={on_val:.4f}  "
          f"visibility 키 {len(written)}개: "
          f"{[(round(t,1), s) for t, s in written]}")
    if FREEZE_SCALE:
        print(f"      scale 애니메이션 제거 -> 고정 {tuple(round(c,4) for c in full)}")
    touched += 1

print("-" * 86)
if not touched:
    sys.exit("!! 처리한 prim 이 없습니다. PATTERNS 를 실제 prim 이름에 맞게 수정하세요.\n"
             "   prim 이름 확인:  python3 -c \"from pxr import Usd;"
             " [print(p.GetPath()) for p in Usd.Stage.Open('%s').Traverse()]\"" % src)

base, _ = os.path.splitext(src)
fixed_usda = base + "_vis.usda"
stage.GetRootLayer().Export(fixed_usda)
print(f"OUT: {fixed_usda}  (prim {touched}개 수정)")

if out_usdz is None:
    out_usdz = base + "_vis.usdz"
if UsdUtils.CreateNewUsdzPackage(Sdf.AssetPath(fixed_usda), out_usdz):
    print(f"OUT: {out_usdz}")
    print("     -> 이 .usdz 를 아이폰/아이패드로 보내서 확인하세요.")
else:
    print(f"!! usdz 패키징 실패. {fixed_usda} 를 Reality Converter 로 변환하세요.")
print("=" * 86)
