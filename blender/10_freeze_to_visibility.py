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
FREEZE_ALL_XFORM = True   # scale 외 translate/rotate 타임샘플도 정적 값으로 고정
PREFER_USDC = True        # pxr 이 있으면 바이너리 usdc 로 담아 파일 크기를 줄인다
SIZE_WARN_MB = 25         # usdz 가 이보다 크면 경고 (AR Quick Look 권장 상한 근처)
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

    # scale 외의 xformOp 타임샘플(translate/rotate/orient)도 정적 값으로 고정한다.
    # 값이 고정이어도 타임샘플이 남아 있으면 '보간할 채널이 없다'는 보장이 아니다.
    # 값 리터럴은 arity 를 모르므로 원문을 그대로 뽑아 쓴다 (3요소/4요소/스칼라 공통).
    freeze_jobs = []
    if FREEZE_ALL_XFORM:
        other = r'(?P<indent>[ \t]*)(?P<decl>\w+\d*\s+xformOp:(?:translate|rotate\w*|orient|transform))'
        for m in re.finditer(other + r'\.timeSamples\s*=', text):
            name = owner(m.start())
            if not any(p in name for p in PATTERNS):
                continue
            span = block_span(m.end())
            if span is None:
                continue
            b0, b1, body = span
            entries = {}
            for line in body.splitlines():
                lm = re.match(rf'\s*({NUM})\s*:\s*(.+?),?\s*$', line)
                if lm:
                    entries[float(lm.group(1))] = lm.group(2).rstrip(",").strip()
            if not entries:
                continue
            freeze_jobs.append({
                "name": name, "start": m.start(), "end": b1,
                "indent": m.group("indent"), "decl": m.group("decl"),
                "entries": entries,
            })

    print("-" * 88)
    report = []
    frozen_report = []
    frames = range(int(t0), int(t1) + 1)

    # ON 구간 시작 프레임 (translate/orient 를 어느 프레임 값으로 고정할지 결정)
    on_start = {}
    for job in jobs:
        mag = {f: max(abs(c) for c in v) for f, v in job["sm"].items()}
        thr = max(mag.values()) * 0.5
        hits = [f for f in frames if lerp(mag, f) > thr]
        if hits:
            on_start[job["name"]] = hits[0]

    # ★ 치환 텍스트를 먼저 전부 만들고, 단 한 번의 내림차순 패스로 적용한다.
    #   두 번에 나눠 적용하면 앞선 치환이 뒤쪽 job 의 인덱스를 밀어 usda 가 깨진다.
    patches = []

    for job in jobs:                                   # scale -> 정적 + visibility
        sm = job["sm"]
        mag = {f: max(abs(c) for c in v) for f, v in sm.items()}
        thr = max(mag.values()) * 0.5

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
        patches.append((job["start"], job["end"], "\n".join(lines)))

        on_frames = [f for f in frames if lerp(mag, f) > thr]
        report.append((job["name"], len(sm), vis,
                       (min(on_frames), max(on_frames)) if on_frames else None,
                       full))

    for job in freeze_jobs:                            # translate/orient -> 정적
        ent = job["entries"]
        pick = on_start.get(job["name"])
        if pick is None or pick not in ent:
            pick = min(ent)
        raw = ent[pick]
        patches.append((job["start"], job["end"],
                        f'{job["indent"]}{job["decl"]} = {raw}'))
        frozen_report.append((job["name"], job["decl"].split()[-1],
                              len(ent), len(set(ent.values())), raw))

    # 구간이 겹치면 치환이 서로를 깨뜨린다. 겹침을 먼저 확인한다.
    spans = sorted((a, b) for a, b, _ in patches)
    for (a1, b1), (a2, b2) in zip(spans, spans[1:]):
        if a2 < b1:
            sys.exit(f"!! 치환 구간이 겹칩니다 ({a1}-{b1} vs {a2}-{b2}). 중단합니다.")

    patched = text
    for start, end, repl in sorted(patches, key=lambda x: x[0], reverse=True):
        patched = patched[:start] + repl + patched[end:]

    for name, n_sm, vis, win, full in sorted(report):
        print(f"  {name}")
        print(f"    scale 타임샘플 {n_sm}개 -> 정적 ({full[0]:g}, {full[1]:g}, {full[2]:g})")
        print(f"    ON 구간 {win[0]}-{win[1]}" if win else "    ON 구간 없음")
        print(f"    visibility {len(vis)}키: {vis}")

    if frozen_report:
        print("-" * 88)
        print("  scale 외 xformOp 고정:")
        for name, opname, n, n_uniq, raw in sorted(frozen_report):
            note = "값 고정이었음" if n_uniq == 1 else f"★ 값이 {n_uniq}종류였음"
            print(f"    {name:<26} {opname:<22} 샘플 {n}개 -> {raw}  ({note})")

    base, _ = os.path.splitext(src)
    out_usda = base + "_frozen.usda"
    with open(out_usda, "w", encoding="utf-8") as fh:
        fh.write(patched)
    print("-" * 88)
    print(f"OUT: {out_usda}")

    # ---- 루트 레이어: 가능하면 바이너리 usdc 로 (ASCII usda 는 수십 배 크다) ----
    root_name = os.path.basename(out_usda)
    root_bytes = patched.encode("utf-8")
    usdc_note = ""
    if PREFER_USDC:
        try:
            # ★ Sdf.Layer 로 변환해야 한다. Usd.Stage.Open().Export() 는 리졸버를
            #   거치면서 @./textures/x.exr@ 를 절대 경로로 바꿔버린다. 그러면
            #   usdz 안의 상대 경로 참조가 깨져 기기에서 텍스처가 안 보인다.
            from pxr import Sdf                       # 있으면 쓰고, 없으면 그냥 usda
            out_usdc = base + "_frozen.usdc"
            lay = Sdf.Layer.FindOrOpen(out_usda)
            if lay and lay.Export(out_usdc):
                with open(out_usdc, "rb") as fh:
                    root_bytes = fh.read()
                root_name = os.path.basename(out_usdc)
                usdc_note = (f"  (usda {len(patched.encode('utf-8')):,} 바이트"
                             f" -> usdc {len(root_bytes):,} 바이트,"
                             f" {len(patched.encode('utf-8'))/max(len(root_bytes),1):.1f}배 축소)")
                print(f"OUT: {out_usdc}{usdc_note}")
        except ImportError:
            usdc_note = ("  ! pxr 없음 -> ASCII usda 를 담습니다. 파일이 큽니다.\n"
                         "    작게 만들려면: pip install usd-core 후 다시 실행")
        except Exception as exc:
            usdc_note = f"  ! usdc 변환 실패({exc}) -> usda 를 담습니다."

    # ---- usdz 패키징 ----
    # 에셋은 usda 가 참조하는 '상대 경로 그대로' 담아야 참조가 깨지지 않는다.
    files = [(root_name, root_bytes)]
    missing = []
    for a in assets:
        arc = a[2:] if a.startswith("./") else a
        if os.path.isabs(arc):
            arc = os.path.basename(arc)
        full = a if os.path.isabs(a) else os.path.join(os.path.dirname(src), a)
        if os.path.exists(full):
            with open(full, "rb") as fh:
                files.append((arc, fh.read()))
        else:
            missing.append(a)
    if missing:
        print(f"  ! usdz 에 담지 못한 에셋: {missing}")
        print("    (해당 머티리얼의 텍스처가 기기에서 안 보일 수 있습니다)")

    out_usdz = base + "_frozen.usdz"
    write_usdz(out_usdz, files)
    size = os.path.getsize(out_usdz)
    print(f"OUT: {out_usdz}  ({size:,} 바이트 = {size/1048576:.1f} MB,"
          f" 파일 {len(files)}개)")
    if usdc_note and not usdc_note.startswith("  ("):
        print(usdc_note)

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
    checks = [("scale 애니메이션 잔존", n_sc_ts == 0),
              ("행렬 애니메이션 잔존", n_mat_ts == 0),
              ("visibility 개수 불일치", n_vis == len(jobs))]
    if FREEZE_ALL_XFORM:
        checks.append(("translate 타임샘플 잔존", n_tr_ts == 0))
    for label, cond in checks:
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

    if size > SIZE_WARN_MB * 1048576:
        print(f"      ★ usdz 가 {size/1048576:.1f} MB 입니다"
              f" (권장 {SIZE_WARN_MB} MB 이하).")
        print("        AR Quick Look 에서 로딩이 느리거나 실패할 수 있습니다.")
        if "pxr 없음" in usdc_note:
            print("        pip install usd-core 후 다시 실행하면 크게 줄어듭니다.")
        else:
            print("        메쉬 폴리곤을 줄이거나(데시메이트) 단백질 표현을 단순화하세요.")

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
