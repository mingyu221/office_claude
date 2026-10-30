#!/usr/bin/env python3
# =============================================================================
# 10_freeze_to_visibility.py
#
#   scale 애니메이션을 USD visibility 로 완전히 치환한다.
#   이동이 '줄어드는' 게 아니라 '원리적으로 불가능'해진다.
#
#   USD 의 visibility 는 숫자가 아니라 token("inherited"/"invisible") 이다.
#   token 은 보간이 정의되지 않으므로 ARKit/RealityKit 이 중간값을 만들 수 없다.
#   scale 을 정적 값으로 고정하면 보간할 트랜스폼 채널이 하나도 남지 않는다.
#
#   ■ 의존성 없음. python3 만 있으면 된다 (pip install 불필요).
#     usda 는 텍스트이고, usdz 는 무압축 zip 이라 표준 라이브러리로 처리한다.
#
#   ■ 사용법
#       python3 10_freeze_to_visibility.py o2_tunnel_retimed.usda
#     결과: o2_tunnel_retimed_frozen.usda, o2_tunnel_retimed_frozen.usdz
# =============================================================================
import os
import re
import struct
import sys
import zlib

# ---------------------------- 설정 -------------------------------------------
PATTERNS = ("O2_Trajectory_F", "O2", "OXY")   # 대상 prim 이름 패턴
# -----------------------------------------------------------------------------

NUM = r'-?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?'


# ============================ usdz 패키징 ====================================
def write_usdz(path, files):
    """USDZ = 무압축(ZIP_STORED) zip + 각 파일 데이터가 64바이트 정렬.

    files: [(arcname, bytes)] — 첫 항목이 루트 레이어여야 한다.
    정렬은 local file header 의 extra field 를 0으로 채워 맞춘다.
    """
    ALIGN = 64
    entries = []
    with open(path, "wb") as fh:
        for arcname, data in files:
            name = arcname.encode("utf-8")
            header_len = 30 + len(name)
            pad = (-(fh.tell() + header_len)) % ALIGN
            offset = fh.tell()
            crc = zlib.crc32(data) & 0xFFFFFFFF
            fh.write(struct.pack(
                "<IHHHHHIIIHH",
                0x04034B50,          # local file header signature
                20,                  # version needed
                0,                   # flags
                0,                   # method: 0 = stored
                0, 0,                # mod time, mod date
                crc, len(data), len(data),
                len(name), pad))
            fh.write(name)
            fh.write(b"\x00" * pad)
            if fh.tell() % ALIGN != 0:
                raise RuntimeError("정렬 실패")
            fh.write(data)
            entries.append((name, crc, len(data), offset))

        cd_start = fh.tell()
        for name, crc, size, offset in entries:
            fh.write(struct.pack(
                "<IHHHHHHIIIHHHHHII",
                0x02014B50,          # central directory signature
                20, 20,              # version made by, version needed
                0, 0,                # flags, method
                0, 0,                # mod time, mod date
                crc, size, size,
                len(name), 0, 0,     # namelen, extralen, commentlen
                0, 0,                # disk start, internal attrs
                0,                   # external attrs
                offset))
            fh.write(name)
        cd_end = fh.tell()

        fh.write(struct.pack(
            "<IHHHHIIH",
            0x06054B50,              # end of central directory
            0, 0,
            len(entries), len(entries),
            cd_end - cd_start, cd_start, 0))
    return path


# ============================== 메인 ========================================
def main():
    if len(sys.argv) < 2:
        sys.exit("usage: python3 10_freeze_to_visibility.py <file.usda>")

    src = sys.argv[1]
    with open(src, "r", encoding="utf-8") as fh:
        text = fh.read()

    def meta(key):
        m = re.search(rf'^\s*{key}\s*=\s*"?({NUM}|[A-Za-z]+)"?\s*$', text, re.M)
        return m.group(1) if m else None

    t0 = float(meta("startTimeCode") or 0)
    t1 = float(meta("endTimeCode") or 0)
    tps = float(meta("timeCodesPerSecond") or meta("framesPerSecond") or 0)

    print("=" * 88)
    print(f"IN : {src}")
    print(f"     timeCode {t0}-{t1}   timeCodesPerSecond={tps}")

    # 외부 에셋 참조가 있으면 usdz 에 같이 담아야 한다
    assets = [a for a in re.findall(r'@([^@\n]+)@', text)
              if not a.startswith("./") or os.path.splitext(a)[1]]
    assets = [a for a in assets if os.path.splitext(a)[1].lower() not in ("", ".usda", ".usdc", ".usd")]
    if assets:
        print(f"     ! 외부 에셋 참조 {len(assets)}개: {assets[:5]}")

    prims = [(m.start(), m.group(1))
             for m in re.finditer(r'^\s*def\s+\w+\s+"([^"]+)"', text, re.M)]

    def owner(pos):
        name = "<stage>"
        for off, nm in prims:
            if off <= pos:
                name = nm
            else:
                break
        return name

    def block_span(pos):
        """{...} 블록의 (시작인덱스, 끝인덱스+1, 본문)"""
        i = text.find("{", pos)
        if i < 0:
            return None
        depth, j = 0, i
        while j < len(text):
            if text[j] == "{":
                depth += 1
            elif text[j] == "}":
                depth -= 1
                if depth == 0:
                    return i, j + 1, text[i + 1:j]
            j += 1
        return None

    def lerp(sm, f):
        ks = sorted(sm)
        if f <= ks[0]:
            return sm[ks[0]]
        if f >= ks[-1]:
            return sm[ks[-1]]
        for a, b in zip(ks, ks[1:]):
            if a <= f <= b:
                t = (f - a) / (b - a) if b > a else 0.0
                return sm[a] + (sm[b] - sm[a]) * t
        return sm[ks[-1]]

    # scale 타임샘플 블록을 찾아 뒤에서부터 치환 (앞에서 바꾸면 인덱스가 밀린다)
    jobs = []
    for m in re.finditer(r'(?P<indent>[ \t]*)(?P<decl>\w+3?\s+xformOp:scale)\.timeSamples\s*=',
                         text):
        name = owner(m.start())
        if not any(p in name for p in PATTERNS):
            continue
        span = block_span(m.end())
        if span is None:
            continue
        b0, b1, body = span
        sm = {}
        for mm in re.finditer(rf'({NUM})\s*:\s*\(\s*({NUM})\s*,\s*({NUM})\s*,\s*({NUM})\s*\)',
                              body):
            sm[float(mm.group(1))] = tuple(float(mm.group(i)) for i in (2, 3, 4))
        if not sm:
            continue
        jobs.append({
            "name": name, "start": m.start(), "end": b1,
            "indent": m.group("indent"), "decl": m.group("decl"), "sm": sm,
        })

    if not jobs:
        sys.exit("!! 대상 prim 의 xformOp:scale.timeSamples 를 찾지 못했습니다.\n"
                 f"   PATTERNS={PATTERNS} 를 실제 prim 이름에 맞게 수정하세요.")

    print("-" * 88)
    patched = text
    report = []
    for job in sorted(jobs, key=lambda j: j["start"], reverse=True):
        sm = job["sm"]
        mag = {f: max(abs(c) for c in v) for f, v in sm.items()}
        on_val = max(mag.values())
        thr = on_val * 0.5

        # ON 구간 복원 (선형 보간 기준) -> visibility 전환 시점 산출
        frames = range(int(t0), int(t1) + 1)
        vis, prev = [], None
        for f in frames:
            cur = lerp(mag, f) > thr
            if cur != prev:
                vis.append((f, "inherited" if cur else "invisible"))
                prev = cur
        if not vis:
            vis = [(int(t0), "invisible")]

        full = next(v for f, v in sm.items() if mag[f] > thr)
        ind = job["indent"]
        lines = [f'{ind}{job["decl"]} = ({full[0]:g}, {full[1]:g}, {full[2]:g})',
                 f'{ind}token visibility.timeSamples = {{']
        for f, tok in vis:
            lines.append(f'{ind}    {f}: "{tok}",')
        lines.append(f'{ind}}}')
        patched = patched[:job["start"]] + "\n".join(lines) + patched[job["end"]:]

        on_frames = [f for f in frames if lerp(mag, f) > thr]
        report.append((job["name"], len(sm), vis,
                       (min(on_frames), max(on_frames)) if on_frames else None,
                       full))

    for name, n_sm, vis, win, full in sorted(report):
        print(f"  {name}")
        print(f"    scale 타임샘플 {n_sm}개 -> 정적 ({full[0]:g}, {full[1]:g}, {full[2]:g})")
        print(f"    ON 구간 {win[0]}-{win[1]}" if win else "    ON 구간 없음")
        print(f"    visibility {len(vis)}키: {vis}")

    base, _ = os.path.splitext(src)
    out_usda = base + "_frozen.usda"
    with open(out_usda, "w", encoding="utf-8") as fh:
        fh.write(patched)
    print("-" * 88)
    print(f"OUT: {out_usda}")

    # usdz 패키징
    files = [(os.path.basename(out_usda), patched.encode("utf-8"))]
    missing = []
    for a in assets:
        p = a if os.path.isabs(a) else os.path.join(os.path.dirname(src), a)
        if os.path.exists(p):
            with open(p, "rb") as fh:
                files.append((os.path.basename(a), fh.read()))
        else:
            missing.append(a)
    if missing:
        print(f"  ! usdz 에 담지 못한 에셋: {missing}")

    out_usdz = base + "_frozen.usdz"
    write_usdz(out_usdz, files)
    size = os.path.getsize(out_usdz)
    print(f"OUT: {out_usdz}  ({size:,} 바이트, 파일 {len(files)}개)")

    # ---------------------------- 자체 검증 ----------------------------------
    print("-" * 88)
    ok = True
    chk = patched
    n_sc_ts = len(re.findall(r'xformOp:scale\.timeSamples', chk))
    n_tr_ts = len(re.findall(r'xformOp:translate\.timeSamples', chk))
    n_mat_ts = len(re.findall(r'xformOp:transform\.timeSamples', chk))
    n_vis = len(re.findall(r'\bvisibility\.timeSamples', chk))
    print(f"검증  scale 타임샘플 {n_sc_ts}개 (0이어야 함)")
    print(f"      translate 타임샘플 {n_tr_ts}개")
    print(f"      단일 4x4 행렬 타임샘플 {n_mat_ts}개 (0이어야 함)")
    print(f"      visibility 타임샘플 {n_vis}개 (대상 {len(jobs)}개와 같아야 함)")
    for label, cond in (("scale 애니메이션 잔존", n_sc_ts == 0),
                        ("행렬 애니메이션 잔존", n_mat_ts == 0),
                        ("visibility 개수 불일치", n_vis == len(jobs))):
        if not cond:
            print(f"      ★ {label}")
            ok = False

    # usdz 규격 재확인: local file header 를 직접 훑어 무압축 + 64바이트 정렬 검사
    import zipfile
    try:
        with open(out_usdz, "rb") as fh:
            blob = fh.read()
        pos, seen = 0, []
        while pos + 30 <= len(blob) and blob[pos:pos + 4] == b"PK\x03\x04":
            (_, _, _, method, _, _, _, csize, usize,
             nlen, elen) = struct.unpack("<IHHHHHIIIHH", blob[pos:pos + 30])
            nm = blob[pos + 30:pos + 30 + nlen].decode("utf-8", "replace")
            data_off = pos + 30 + nlen + elen
            seen.append((nm, method, data_off, usize))
            pos = data_off + csize
        if not seen:
            print("      ★ local file header 를 찾지 못함")
            ok = False
        for nm, method, off, usize in seen:
            a, st = off % 64 == 0, method == 0
            print(f"      usdz: {nm}  dataOffset={off}"
                  f" {'64정렬 OK' if a else '★정렬 실패'},"
                  f" {'무압축 OK' if st else '★압축됨'}, {usize:,}바이트")
            if not (a and st):
                ok = False
        with zipfile.ZipFile(out_usdz) as zf:
            bad = zf.testzip()
            if bad:
                print(f"      ★ zip 손상: {bad}")
                ok = False
    except Exception as exc:
        print(f"      ★ usdz 검사 실패: {exc}")
        ok = False

    print("=" * 88)
    if ok:
        print("✓ 통과. 트랜스폼은 전부 정적이고 visibility 만 애니메이션됩니다.")
        print("  visibility 는 token 이라 ARKit 이 보간할 수 없습니다.")
        print("  -> 이동이 '줄어든' 게 아니라 '불가능'합니다.")
        print(f"  {os.path.basename(out_usdz)} 를 기기로 보내 확인하세요.")
    else:
        print("★ 위 경고를 확인하세요.")
    print("=" * 88)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
