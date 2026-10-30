#!/usr/bin/env python3
# =============================================================================
# 04_verify_usda.py  -  내보낸 .usda 안을 들여다보고 "날아오는 원인"이 남았는지 검증
#
#  사용법:  python3 04_verify_usda.py /path/to/o2_tunnel.usda
#  (블렌더 없이 그냥 파이썬으로 돌아갑니다)
#
#  검사 항목
#   1. xformOp:translate 에 .timeSamples 가 붙어 있는가        -> 있으면 위치가 애니메이션 중
#   2. translate 타임샘플 값이 실제로 변하는가                  -> 변하면 그게 '날아오는' 원인
#   3. xformOp:transform(단일 4x4 행렬)으로 나갔는가            -> TRS 로 바꾸는 게 안전
#   4. xformOp:scale 타임샘플이 0 인 프레임이 있는가            -> ARKit 좌표 리셋 버그
#   5. visibility 타임샘플이 있는가                            -> 있으면 가장 확실한 팝 처리
# =============================================================================
import re
import sys
from collections import defaultdict

if len(sys.argv) < 2:
    sys.exit("usage: python3 04_verify_usda.py <file.usda>")

path = sys.argv[1]
with open(path, "r", encoding="utf-8", errors="replace") as fh:
    text = fh.read()

prim_re = re.compile(r'^\s*def\s+\w+\s+"([^"]+)"', re.M)
prims = [(m.start(), m.group(1)) for m in prim_re.finditer(text)]


def owner(pos):
    """해당 위치가 속한 가장 가까운 상위 prim 이름."""
    name = "<stage>"
    for off, n in prims:
        if off <= pos:
            name = n
        else:
            break
    return name


def block_at(pos):
    """attr.timeSamples = { ... } 의 중괄호 블록 본문을 돌려준다."""
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
vec3_re = re.compile(rf'\(\s*({num})\s*,\s*({num})\s*,\s*({num})\s*\)')

problems, notes = [], []

# --- 1,2) translate -----------------------------------------------------------
ts_translate = list(re.finditer(r'xformOp:translate\.timeSamples\s*=', text))
static_translate = len(re.findall(r'xformOp:translate\s*=', text))

moving = defaultdict(list)
for m in ts_translate:
    who = owner(m.start())
    vals = {tuple(round(float(x), 6) for x in g) for g in vec3_re.findall(block_at(m.end()))}
    if len(vals) > 1:
        xs = sorted(vals)
        moving[who] = [xs[0], xs[-1], len(vals)]

# --- 3) 단일 행렬 ------------------------------------------------------------
ts_matrix = {owner(m.start()) for m in
             re.finditer(r'xformOp:transform\.timeSamples\s*=', text)}

# --- 4) scale 0 --------------------------------------------------------------
zero_scale = defaultdict(list)
for m in re.finditer(r'xformOp:scale\.timeSamples\s*=', text):
    who = owner(m.start())
    body = block_at(m.end())
    for fm in re.finditer(rf'({num})\s*:\s*\(\s*({num})\s*,', body):
        if abs(float(fm.group(2))) < 1e-9:
            zero_scale[who].append(fm.group(1))

# --- 5) visibility -----------------------------------------------------------
ts_vis = {owner(m.start()) for m in
          re.finditer(r'\bvisibility\.timeSamples\s*=', text)}

# ------------------------------- 리포트 --------------------------------------
print("=" * 86)
print(f"FILE: {path}")
print(f"  prim 수: {len(prims)}")
print(f"  xformOp:translate  -> 정적 {static_translate}개 / 타임샘플 {len(ts_translate)}개")
print(f"  xformOp:scale 타임샘플: {len(re.findall(r'xformOp:scale.timeSamples', text))}개")
print(f"  visibility 타임샘플: {len(ts_vis)}개 {sorted(ts_vis) if ts_vis else ''}")
print("-" * 86)

if moving:
    problems.append("translate 가 실제로 움직이고 있음")
    print("★ [문제] 아래 prim 의 위치가 애니메이션 중입니다. 이게 '날아오는' 직접 원인입니다:")
    for who, (lo, hi, n) in moving.items():
        print(f"    - {who}: {n}개의 서로 다른 값, {lo} ... {hi}")
    print("    -> 02_fix_origin_and_rebuild.py 가 location F-Curve 를 정말 지웠는지 확인하세요.")
elif ts_translate:
    print("○ translate 에 타임샘플은 있지만 값이 전부 동일합니다 (움직이지 않음). 무해합니다.")
else:
    print("✓ translate 는 정적 값으로만 기록됐습니다. 위치 보간 자체가 불가능합니다. (최상)")

if ts_matrix:
    problems.append("단일 4x4 행렬로 내보내짐")
    print(f"★ [문제] xformOp:transform(단일 행렬)에 타임샘플이 있음: {sorted(ts_matrix)}")
    print("    -> 내보낼 때 xform_op_mode='TRS' 로 바꾸세요. RealityKit 의 행렬 분해를 피합니다.")
else:
    print("✓ 단일 4x4 행렬 애니메이션 없음 (TRS 분리 상태).")

if zero_scale:
    problems.append("scale 0 프레임 존재")
    print("★ [문제] scale 이 완전한 0 인 프레임이 있습니다 (ARKit 좌표 리셋 버그):")
    for who, frames in zero_scale.items():
        print(f"    - {who}: 프레임 {frames[:10]}{' ...' if len(frames) > 10 else ''}")
    print("    -> OFF_SCALE 을 0.001 같은 작은 양수로 두세요.")
else:
    print("✓ scale 0 프레임 없음.")

if ts_vis:
    print("✓ visibility 애니메이션이 있습니다. token 이라 ARKit 이 보간할 수 없어 깔끔히 팝됩니다.")
else:
    notes.append("visibility 애니메이션 없음 (scale 트릭에만 의존)")
    print("○ visibility 애니메이션 없음. scale 트릭만으로 동작합니다.")
    print("    -> 더 확실하게 하려면 02 스크립트의 VIS_KEYS=True, 또는 05_inject_visibility.py 사용.")

print("-" * 86)
if problems:
    print(f"판정: 문제 {len(problems)}건 -> {problems}")
    sys.exit(1)
print("판정: 통과. 이 .usda 와 같은 설정으로 나온 .usdz 를 기기에서 확인하세요.")
if notes:
    print(f"참고: {notes}")
print("=" * 86)
