#!/usr/bin/env python3
# =============================================================================
# 12_verify_center.py  -  내보낸 USD 의 '실제' 월드 바운딩 박스를 측정
#
#   AR Quick Look 이 조작 회전축으로 쓰는 값을 그대로 계산한다.
#   UsdGeom.BBoxCache 는 ARKit 과 동일한 방식으로 월드 bbox 를 구한다.
#
#   필요: pip install usd-core
#   사용: python3 12_verify_center.py o2_tunnel_centered.usda
# =============================================================================
import sys

try:
    from pxr import Usd, UsdGeom, Gf
except ImportError:
    sys.exit("pxr 가 없습니다.  pip3 install usd-core  후 다시 실행하세요.")

if len(sys.argv) < 2:
    sys.exit("usage: python3 12_verify_center.py <file.usda|file.usdz>")

path = sys.argv[1]
stage = Usd.Stage.Open(path)
if not stage:
    sys.exit(f"열 수 없습니다: {path}")

t0, t1 = stage.GetStartTimeCode(), stage.GetEndTimeCode()
print("=" * 90)
print(f"FILE: {path}")
print(f"      timeCode {t0}-{t1}   upAxis={UsdGeom.GetStageUpAxis(stage)}")
print("-" * 90)


def rng(prim, tc):
    cache = UsdGeom.BBoxCache(Usd.TimeCode(tc), [UsdGeom.Tokens.default_])
    return cache.ComputeWorldBound(prim).ComputeAlignedRange()


def show(label, r):
    if r.IsEmpty():
        print(f"  {label:<34} (빈 bbox)")
        return None
    mn, mx, mid = r.GetMin(), r.GetMax(), r.GetMidpoint()
    sz = mx - mn
    print(f"  {label:<34} 중심 ({mid[0]:8.3f},{mid[1]:8.3f},{mid[2]:8.3f})"
          f"  크기 ({sz[0]:7.2f},{sz[1]:7.2f},{sz[2]:7.2f})")
    return mid


# 씬 전체 (= AR 회전축). 시간에 따라 변할 수 있으니 여러 시점에서 본다
print("[1] 씬 전체 월드 bbox  <- AR Quick Look 이 이 중심을 회전축으로 씀")
mids = []
root = stage.GetPseudoRoot()
for tc in sorted({t0, t0 + (t1 - t0) * 0.25, (t0 + t1) / 2,
                  t0 + (t1 - t0) * 0.75, t1}):
    m = show(f"timeCode {tc:.0f}", rng(root, tc))
    if m is not None:
        mids.append((tc, Gf.Vec3d(m)))

print()
print("[2] 개별 prim (timeCode 1 기준)")
for prim in stage.Traverse():
    if not prim.IsA(UsdGeom.Xformable):
        continue
    if prim.GetParent() != stage.GetPseudoRoot() and prim.GetParent().GetName() != "root":
        continue
    if prim.GetName() == "root":
        continue
    show(prim.GetName(), rng(prim, t0))

print("-" * 90)
if not mids:
    sys.exit("!! bbox 를 계산할 수 없습니다 (지오메트리가 없나요?)")

worst_tc, worst = max(mids, key=lambda x: Gf.Vec3d(x[1]).GetLength())
worst_len = Gf.Vec3d(worst).GetLength()
# 단백질 반경을 척도로
prot = next((p for p in stage.Traverse()
             if "Protein" in p.GetName() and p.IsA(UsdGeom.Xformable)), None)
scale_ref = 1.0
if prot:
    pr = rng(prot, t0)
    if not pr.IsEmpty():
        scale_ref = (pr.GetMax() - pr.GetMin()).GetLength() / 2

print(f"[3] 판정")
print(f"  회전축이 원점에서 가장 많이 벗어난 시점: timeCode {worst_tc:.0f}")
print(f"  그때의 중심 ({worst[0]:.3f}, {worst[1]:.3f}, {worst[2]:.3f})"
      f"   거리 {worst_len:.3f}")
print(f"  단백질 반경 기준 {worst_len / max(scale_ref, 1e-9) * 100:.1f}%")
print("=" * 90)
if worst_len <= max(1e-3, scale_ref * 0.05):
    print("✓ 통과. AR 회전축이 단백질 중앙에 있습니다.")
else:
    print("★ 회전축이 아직 벗어나 있습니다.")
    print("  11_recenter_for_ar.py 를 최신 버전으로 다시 받아 블렌더에서 재실행하세요")
    print("  (view_layer.update() 누락 버그를 고쳤습니다).")
print("=" * 90)
