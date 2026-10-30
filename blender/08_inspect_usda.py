#!/usr/bin/env python3
# =============================================================================
# 08_inspect_usda.py  -  내보낸 usda 에서 '날아오는' 메커니즘을 확정한다
#
#   python3 08_inspect_usda.py /Users/mingyukim/Desktop/vision/frame/o2_tunnel_visibility.usda
#
#   블렌더 없이, pip 설치 없이 그냥 python3 로 돌아갑니다.
#
#   확인하는 것
#     1. 스테이지 타이밍 (startTimeCode/endTimeCode/timeCodesPerSecond)
#        -> 프레임 1개가 실제 몇 초인지. ARKit 이 보간하는 시간이 이 값이다.
#     2. xformOp:translate 타임샘플 값이 실제로 변하는가, 변하면 몇 단위인가
#        -> 변하면 그게 '날아오는' 직접 원인
#     3. xformOp:scale 타임샘플 범위
#     4. visibility 타임샘플 유무
# =============================================================================
import re
import sys
from collections import OrderedDict

if len(sys.argv) < 2:
    sys.exit("usage: python3 08_inspect_usda.py <file.usda>")

path = sys.argv[1]
with open(path, "r", encoding="utf-8", errors="replace") as fh:
    text = fh.read()

NUM = r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'


def meta(key, default=None):
    m = re.search(rf'^\s*{key}\s*=\s*"?({NUM}|[A-Za-z]+)"?\s*$', text, re.M)
    return m.group(1) if m else default


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


def samples(attr):
    """{prim: OrderedDict{frame: (v0,v1,v2)}}"""
    out = {}
    for m in re.finditer(rf'{re.escape(attr)}\.timeSamples\s*=', text):
        body = block_at(m.end())
        od = OrderedDict()
        for fm in re.finditer(
                rf'({NUM})\s*:\s*\(\s*({NUM})\s*,\s*({NUM})\s*,\s*({NUM})\s*\)', body):
            od[float(fm.group(1))] = tuple(float(fm.group(i)) for i in (2, 3, 4))
        out[owner(m.start())] = od
    return out


t0 = float(meta("startTimeCode", 0) or 0)
t1 = float(meta("endTimeCode", 0) or 0)
tps = float(meta("timeCodesPerSecond", 0) or 0)
fps = float(meta("framesPerSecond", 0) or 0)
rate = tps or fps

print("=" * 92)
print(f"FILE: {path}")
print("-" * 92)
print("[1] 스테이지 타이밍")
print(f"    startTimeCode = {t0}   endTimeCode = {t1}")
print(f"    timeCodesPerSecond = {tps or '(없음)'}   framesPerSecond = {fps or '(없음)'}")
print(f"    upAxis = {meta('upAxis', '(없음)')}   metersPerUnit = {meta('metersPerUnit', '(없음)')}")
if rate > 0:
    dur = (t1 - t0) / rate
    step = 1.0 / rate
    print(f"    -> 전체 길이 {dur:.3f}초,  타임코드 1칸 = {step:.3f}초")
    print(f"    ★ ARKit 이 샘플 사이를 선형 보간하는 시간 = {step:.3f}초")
    if step > 0.1:
        print(f"    ★★ {step:.3f}초는 눈에 보이는 길이입니다."
              " 이게 '스윽 날아오는'의 정체입니다.")
    else:
        print(f"    ○ {step:.3f}초면 사실상 '팝'으로 보입니다.")
else:
    print("    !! 재생 속도 메타데이터가 없습니다.")

tr = samples("xformOp:translate")
sc = samples("xformOp:scale")
vis = {owner(m.start()) for m in re.finditer(r'\bvisibility\.timeSamples\s*=', text)}
static_tr = {m.group(1) for m in
             re.finditer(rf'xformOp:translate\s*=\s*\(\s*{NUM}', text)}

print("-" * 92)
print("[2] xformOp:translate 타임샘플")
moving = {}
if not tr:
    print("    (없음) -> 위치는 정적. 보간 불가능. 최상.")
for name, od in tr.items():
    vals = list(od.values())
    uniq = set(vals)
    if len(uniq) <= 1:
        print(f"    ○ {name:<28} 샘플 {len(vals)}개, 값 전부 동일 {vals[0] if vals else ''}"
              " -> 무해")
        continue
    span = max(
        sum((a - b) ** 2 for a, b in zip(p, q)) ** 0.5
        for p in uniq for q in uniq)
    moving[name] = span
    first, last = vals[0], vals[-1]
    print(f"    ★ {name:<28} 샘플 {len(vals)}개, 서로 다른 값 {len(uniq)}개")
    print(f"      {first}  ->  {last}")
    print(f"      최대 이동 거리 {span:.3f} 단위  <<< '날아오는' 직접 원인")

print("-" * 92)
print("[3] xformOp:scale 타임샘플")
for name, od in sc.items():
    mags = [max(abs(c) for c in v) for v in od.values()]
    if not mags:
        continue
    zero = [f for f, v in od.items() if max(abs(c) for c in v) < 1e-9]
    print(f"    {name:<28} 샘플 {len(mags)}개, {min(mags):.4f} ~ {max(mags):.4f}"
          + (f"   !! scale 0 프레임 {zero}" if zero else ""))

print("-" * 92)
print(f"[4] visibility 타임샘플 : {len(vis)}개 {sorted(vis) if vis else '(없음)'}")

print("=" * 92)
print("결론")
if vis:
    print("  visibility 가 이미 애니메이션돼 있습니다. token 이라 보간 불가 -> 팝.")
else:
    print("  visibility 애니메이션이 없습니다. ARKit 이 translate/scale 을 보간합니다.")
if moving:
    print(f"  ★ translate 가 실제로 움직입니다 ({len(moving)}개,"
          f" 최대 {max(moving.values()):.2f} 단위).")
    print("    오리진 문제가 아니라 USD 안에 '진짜 이동'이 기록돼 있습니다.")
elif tr:
    print("  ○ translate 타임샘플은 있지만 값은 고정입니다. 이동의 원인이 아닙니다.")
if rate > 0 and 1.0 / rate > 0.1 and not vis:
    print(f"  ★ 타임코드 1칸이 {1.0/rate:.3f}초입니다. scale 0.001 -> 1 보간이")
    print(f"    {1.0/rate:.3f}초에 걸쳐 일어나므로, 오리진 오프셋만큼 '천천히' 날아옵니다.")
print("=" * 92)
