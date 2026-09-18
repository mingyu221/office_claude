# Folddisco 민감도 — 독립 노트북

금속 모티프 질의(`A112,A114`)가 **Ch CooC1 의 직교체인 Y19 CooC(`AHZ96930.1`)를 못 찾는다.**
같은 단백질을 Walker A 질의는 찾고, Foldseek 은 tm 0.938 로 찾는다. 구조는 인덱스에 있다.

그러면 지금까지 나온 모든 "0" 은 음성이 아니라 **민감도 밖**일 수 있다.
이 노트북은 그 민감도를 실측한다. S1 → S2 → S3 순.

- **S2** 가 갈림길이다. Y19 CooC 가 *자기 인덱스에서 자기 질의로* 잡히는가?
  - 잡히면 → 문제는 **종간 발산**이다. S3 의 `-d`/`-a` 확장으로 푼다.
  - 안 잡히면 → 그 단백질의 **Cys 배치 자체가 질의 기하와 다르다**. 확장으로 안 풀린다.
- **S3** 은 Ch CooC1 질의의 허용폭을 넓히며 Y19 CooC 가 나타나는 지점을 찾는다.

---

## CELL S1 — 설정

```python
# =============================================================================
# CELL S1 | 독립 실행용 설정
# =============================================================================
import re, subprocess, itertools, math
from pathlib import Path
import pandas as pd

BASE  = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS = Path("/mnt/af2results/mingyu")
TBL   = BASE/"result"/"table"
OUT   = BASE/"result"/"folddisco"/"sensitivity"; OUT.mkdir(parents=True, exist_ok=True)
THREADS = 8

FD  = TOOLS/"folddisco"/"bin"/"folddisco"
IDX = {"BL21":   TOOLS/"database"/"bacteriaDB"/"folddisco"/"UP000503272",
       "Y19":    TOOLS/"database"/"bacteriaDB"/"folddisco"/"UP000034085",
       "MG1655": TOOLS/"database"/"folddisco"/"ecoli_folddisco"/"e_coli_folddisco"}
STRUCT = {"BL21": TOOLS/"database"/"bacteriaDB"/"structures_UP000503272",
          "Y19":  TOOLS/"database"/"bacteriaDB"/"structures_UP000034085"}
COOC1_PDB = TOOLS/"workspace"/"seek_ni_insertase"/"input"/"3kji.pdb"

# Y19 의 CooC — structure_accessions.xlsx 의 Y19 시트에서 확인한 값
Y19_COOC_GB   = "AHZ96930.1"
Y19_COOC_UNI  = "A0A059VK30"
Y19_COOC_FILE = STRUCT["Y19"]/f"AF-{Y19_COOC_UNI}-F1-model_v6.cif"

def sh(cmd, quiet=False):
    r = subprocess.run(cmd, shell=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if not quiet: print(r.stdout.rstrip())
    return r.stdout

print("=" * 96); print("### 경로"); print("=" * 96)
for n, p in [("folddisco", FD), ("3kji", COOC1_PDB), ("Y19 CooC 구조", Y19_COOC_FILE),
             ("Y19 index", IDX["Y19"]), ("BL21 index", IDX["BL21"]),
             ("MG1655 index", IDX["MG1655"]), ("출력", OUT)]:
    mark = "O" if (p.exists() or p.parent.exists()) else "X"
    print(f"  {mark}  {n:16s} {p}")
if not Y19_COOC_FILE.exists():
    print("\n⚠ Y19 CooC 구조 파일이 없다. 확장자·버전이 다를 수 있다:")
    for f in sorted(STRUCT['Y19'].glob(f"*{Y19_COOC_UNI}*"))[:5]: print("   ", f.name)
```

---

## CELL S0 — 3kji 의 조립체를 제대로 읽는다

```python
# =============================================================================
# CELL S0 | REMARK 350 을 조립체 단위로 읽는다
#   앞서 쓴 파서는 BIOMT 행을 연산자 번호로만 모았다. 조립체가 둘이면
#   BIOMOLECULE 1 의 BIOMT1 과 BIOMOLECULE 2 의 BIOMT1 이 같은 키로 덮어써지고,
#   APPLY THE FOLLOWING TO CHAINS 도 무시했다. 그래서 B–B2 가 2.13 Å 으로 나온
#   것이 정당한 조립체인지, 다른 조립체의 연산자를 B 에 잘못 건 결과인지
#   구분할 수 없었다.
#
#   여기서는 조립체마다 따로 만들고, 사슬 적용 범위를 지키고, 각각에서
#   Cys 쌍 거리를 재서 '어느 조립체가 Ni 자리를 만드는가' 를 직접 본다.
# =============================================================================
import itertools, math
RES = (112, 114)
SITE_MAX = 12.0

def parse_remark350(path):
    """[{id, chains:[...], ops:[(R,t)...]}] 로 조립체를 읽는다."""
    asm, cur, ops = [], None, {}
    def flush():
        nonlocal cur, ops
        if cur is not None:
            o = []
            for k in sorted(ops, key=lambda s: int(s)):
                m = ops[k]
                if len(m) == 3:
                    o.append(([m[i][:3] for i in (1, 2, 3)], [m[i][3] for i in (1, 2, 3)]))
            cur["ops"] = o
            asm.append(cur)
        cur, ops = None, {}
    for l in open(path, errors="ignore"):
        if not l.startswith("REMARK 350"): continue
        s = l[10:].strip()
        if s.startswith("BIOMOLECULE:"):
            flush(); cur = {"id": s.split(":")[1].strip(), "chains": [], "ops": []}
        elif "APPLY THE FOLLOWING TO CHAINS" in s and cur is not None:
            cur["chains"] += [c.strip() for c in s.split(":")[1].split(",") if c.strip()]
        elif "AND CHAINS:" in s and cur is not None:
            cur["chains"] += [c.strip() for c in s.split(":")[1].split(",") if c.strip()]
        elif s.startswith("BIOMT"):
            k = int(s[5:6]); idx = s[6:10].strip()
            ops.setdefault(idx, {})[k] = [float(x) for x in s[10:].split()[:4]]
    flush()
    return asm

def load_atoms(path):
    out = []
    for l in open(path, errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")): continue
        rs = l[22:26].strip()
        if rs.lstrip("-").isdigit():
            out.append((l[21], int(rs), l[17:20].strip(), l.rstrip("\n")))
    return out

def xform(atoms, chains, R, t, tag):
    out = []
    for c, r, n, l in atoms:
        if chains and c not in chains: continue
        x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        nx = R[0][0]*x+R[0][1]*y+R[0][2]*z+t[0]
        ny = R[1][0]*x+R[1][1]*y+R[1][2]*z+t[1]
        nz = R[2][0]*x+R[2][1]*y+R[2][2]*z+t[2]
        out.append((f"{c}{tag}", r, n, l[:30] + f"{nx:8.3f}{ny:8.3f}{nz:8.3f}" + l[54:]))
    return out

def sg(atoms, ch, rs):
    for c, r, _n, l in atoms:
        if c == ch and r == rs and l[12:16].strip() == "SG":
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return None

def center(atoms, ch):
    p = [sg(atoms, ch, r) for r in RES]
    if any(x is None for x in p): return None
    return tuple(sum(v[i] for v in p)/2 for i in range(3))

A0  = load_atoms(COOC1_PDB)
ASM = parse_remark350(COOC1_PDB)
print("=" * 96); print(f"### {COOC1_PDB.name} 의 조립체 {len(ASM)}개"); print("=" * 96)

BUILT = {}
for a in ASM:
    parts = []
    for i, (R, t) in enumerate(a["ops"], 1):
        ident = (all(abs(R[x][y]-(1 if x == y else 0)) < 1e-6 for x in range(3) for y in range(3))
                 and all(abs(v) < 1e-6 for v in t))
        parts += xform(A0, a["chains"], R, t, "" if ident else str(i))
    BUILT[a["id"]] = parts
    chs = sorted({c for c, _r, _n, _l in parts})
    print(f"\n  조립체 {a['id']}  적용 사슬 {a['chains']}  연산 {len(a['ops'])}개")
    print(f"    생성 사슬: {chs}")
    cen = {c: center(parts, c) for c in chs}
    cen = {k: v for k, v in cen.items() if v}
    if len(cen) < 2:
        print(f"    Cys{RES[0]}/{RES[1]} 자리를 가진 사슬이 {len(cen)}개 — 자리를 못 만든다")
        continue
    for x, y in itertools.combinations(sorted(cen), 2):
        d = math.dist(cen[x], cen[y])
        print(f"    {x}–{y}  자리 중심 {d:6.2f} Å" + ("   ← 한 자리" if d <= SITE_MAX else ""))

print("\n" + "=" * 96); print("### 판정"); print("=" * 96)
GOOD = []
for aid, parts in BUILT.items():
    chs = sorted({c for c, _r, _n, _l in parts})
    cen = {c: center(parts, c) for c in chs}
    cen = {k: v for k, v in cen.items() if v}
    for x, y in itertools.combinations(sorted(cen), 2):
        d = math.dist(cen[x], cen[y])
        if d <= SITE_MAX: GOOD.append((aid, x, y, d))
if GOOD:
    for aid, x, y, d in sorted(GOOD, key=lambda v: v[3]):
        print(f"  조립체 {aid}: {x}+{y} 가 {d:.2f} Å 로 한 자리를 만든다")
    aid, x, y, d = sorted(GOOD, key=lambda v: v[3])[0]
    print(f"\n  ★ 템플릿은 조립체 {aid} 의 {x}+{y} 로 만든다.")
    print("    S2 이후 셀에서 이 조합을 쓸 것. (앞서 만든 merged.pdb 가 이것과")
    print("    같은 조합이면 결과는 그대로다 — 우연히 맞았던 것이고, 이제 근거가 있다.)")
    ASM_BEST = (aid, x, y, BUILT[aid])
else:
    print("  어느 조립체에서도 두 Cys 쌍이 한 자리를 이루지 않는다.")
    print("  잔기 번호나 구조 형태를 다시 봐야 한다.")
    ASM_BEST = None
```

---

## CELL S2 — Y19 CooC 를 직접 뜯어본다

```python
# =============================================================================
# CELL S2 | 이 단백질의 Cys 가 실제로 어떻게 놓여 있는가, 그리고
#          자기 인덱스에서 자기 질의로 잡히는가
#   (1) 구조에서 모든 CYS 의 SG 좌표를 읽어 쌍 거리를 잰다
#   (2) Ch CooC1 의 질의 기하(A112–A114 = 3.58 Å)와 견줄 만한 쌍이 있는가
#   (3) 있으면 그 쌍으로 Y19 인덱스에 질의한다 — 자기 자신은 반드시 나와야 한다
#
#   (3) 이 실패하면 확장으로 풀 수 있는 문제가 아니다. 그 단백질의 Cys 배치가
#   질의 기하와 애초에 다르다는 뜻이고, 그건 AlphaFold 가 열린(apo) 형태를
#   예측했을 때 생긴다 — 3kji 의 ASU 에서 사슬 간이 60 Å 이던 것과 같은 상황이다.
# =============================================================================
assert Y19_COOC_FILE.exists(), f"구조 없음: {Y19_COOC_FILE}"

def cys_sg(path):
    """CYS 의 SG 좌표. mmCIF / PDB 둘 다 읽는다."""
    out = {}
    txt = path.read_text(errors="ignore").splitlines()
    if path.suffix == ".cif":
        cols, inloop = [], False
        for l in txt:
            s = l.strip()
            if s.startswith("_atom_site."): cols.append(s.split(".")[1]); inloop = True; continue
            if inloop and (s.startswith("#") or not s): inloop = False; continue
            if inloop and (s.startswith("ATOM") or s.startswith("HETATM")):
                f = s.split()
                d = dict(zip(cols, f))
                if d.get("label_comp_id") == "CYS" and d.get("label_atom_id") == "SG":
                    out[int(d["label_seq_id"])] = (float(d["Cartn_x"]), float(d["Cartn_y"]),
                                                   float(d["Cartn_z"]))
    else:
        for l in txt:
            if l.startswith("ATOM") and l[17:20].strip() == "CYS" and l[12:16].strip() == "SG":
                out[int(l[22:26])] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return out

SG = cys_sg(Y19_COOC_FILE)
print("=" * 96); print(f"### (1) {Y19_COOC_GB} 의 CYS"); print("=" * 96)
print(f"  CYS {len(SG)}개: {sorted(SG)}")
if len(SG) < 2:
    print("  ⚠ Cys 가 2개 미만이다. 이 단백질에는 Cys 쌍 자리가 없다.")

REF = 3.58        # Ch CooC1 A112–A114 의 SG–SG 거리
print(f"\n  SG–SG 쌍 거리 (Ch CooC1 기준 {REF} Å)")
pairs = []
for a, b in itertools.combinations(sorted(SG), 2):
    d = math.dist(SG[a], SG[b])
    if d <= 12:
        pairs.append((abs(d - REF), d, a, b))
        print(f"    A{a:<5d} A{b:<5d} {d:6.2f}   Δ(기준) {d-REF:+6.2f}"
              + ("   ← 기준과 가장 가까움" if False else ""))
pairs.sort()
if not pairs:
    print("    12 Å 안에 있는 Cys 쌍이 없다 — 이 모델에는 금속 자리가 열려 있다.")
print("\n" + "=" * 96); print("### (2) 자기 질의 — 자기 자신을 찾는가"); print("=" * 96)
SELF_OK = None
if pairs:
    _, d0, i, j = pairs[0]
    q = f"A{i},A{j}"
    print(f"  질의 {q}  (SG–SG {d0:.2f} Å, 기준 {REF} Å 과 Δ {d0-REF:+.2f})")
    o = OUT/f"self_y19cooc_{i}_{j}.tsv"
    if not (o.exists() and o.stat().st_size):
        sh(f'"{FD}" query -p "{Y19_COOC_FILE}" -q {q} -i "{IDX["Y19"]}" -t {THREADS} '
           f'--per-structure --header --sort-by idf --top 5000 -o "{o}"')
    if o.exists() and o.stat().st_size:
        dd = pd.read_csv(o, sep="\t"); dd.columns = [c.strip().lstrip("#") for c in dd.columns]
        me = dd[dd.tid.astype(str).str.contains(Y19_COOC_UNI, na=False)]
        SELF_OK = len(me) > 0
        print(f"  히트 {len(dd)}개, 그중 자기 자신 {len(me)}개")
        if len(me): print(me.to_string(index=False))
    else:
        SELF_OK = False; print("  결과 없음")

print("\n" + "=" * 96); print("### (3) 갈림길"); print("=" * 96)
if SELF_OK:
    print("  자기 질의로는 잡힌다. 즉 이 단백질의 Cys 자리는 인덱스에 제대로 있고")
    print("  folddisco 도 그것을 색인했다.")
    print("  ★ 그렇다면 Ch CooC1 질의가 이것을 놓친 이유는 종간 발산이다.")
    print("    S3 에서 -d / -a 를 넓혀 어디서 잡히기 시작하는지 찾는다.")
elif SELF_OK is False:
    print("  자기 질의로도 안 잡힌다. -d/-a 확장으로 풀 문제가 아니다.")
    print("  ★ 이 모델의 Cys 배치가 질의 기하와 다르다 — AlphaFold 가 열린(apo)")
    print("    형태를 예측했을 가능성이 크다. 위 (1) 의 쌍 거리를 볼 것.")
    print("    그 경우 할 일은 ADP 결합 형태나 이량체 조립체를 질의로 쓰는 것이지")
    print("    허용폭을 넓히는 것이 아니다.")
else:
    print("  Cys 쌍이 없어 판정 불가.")
```

---

## CELL S3 — `-d` / `-a` 민감도 스윕

```python
# =============================================================================
# CELL S3 | Ch CooC1 질의의 허용폭을 넓히며 Y19 CooC 가 나타나는 지점을 찾는다
#   -d 거리 허용폭 (기본 0.5 Å) / -a 각도 허용폭 (기본 5.0°)
#   찾는 것은 두 가지다.
#     (a) AHZ96930.1 이 처음 잡히는 설정 — 이 방법의 실제 민감도
#     (b) 그때 세 균주에서 몇 개가 딸려 나오는가 — 그 민감도의 비용
#   (a) 를 모르면 어떤 0 도 음성으로 읽을 수 없다.
# =============================================================================
GRID = [(0.5, 5.0), (1.0, 10.0), (1.5, 15.0), (2.0, 20.0), (3.0, 30.0)]
QRES = "A112,A114"          # Ch CooC1 의 Ni 배위 Cys 쌍 (3kji 체인 A)
assert COOC1_PDB.exists(), f"3kji 없음: {COOC1_PDB}"

rows = []
for d_, a_ in GRID:
    for s, idx in IDX.items():
        o = OUT/f"sweep_d{d_}_a{a_}_{s}.tsv"
        if not (o.exists() and o.stat().st_size):
            sh(f'"{FD}" query -p "{COOC1_PDB}" -q {QRES} -i "{idx}" -t {THREADS} '
               f'-d {d_} -a {a_} --per-structure --header --sort-by idf --top 20000 -o "{o}"',
               quiet=True)
        n, found = 0, False
        if o.exists() and o.stat().st_size:
            dd = pd.read_csv(o, sep="\t"); dd.columns = [c.strip().lstrip("#") for c in dd.columns]
            n = len(dd)
            if s == "Y19":
                found = dd.tid.astype(str).str.contains(Y19_COOC_UNI, na=False).any()
        rows.append({"-d": d_, "-a": a_, "strain": s, "hits": n,
                     "Y19_CooC": ("O" if found else ".") if s == "Y19" else ""})
        print(f"  d={d_:<4} a={a_:<5} {s:8s} {n:6d}개"
              + (f"   CooC {'O' if found else 'X'}" if s == "Y19" else ""))
S = pd.DataFrame(rows)
P = S.pivot_table(index=["-d", "-a"], columns="strain", values="hits").reindex(
        columns=["BL21", "MG1655", "Y19"])
P["Y19_CooC"] = S[S.strain == "Y19"].set_index(["-d", "-a"])["Y19_CooC"]
print("\n" + "=" * 96); print("### 스윕 결과"); print("=" * 96)
print(P.to_string())
P.to_csv(TBL/"folddisco_sensitivity_sweep.csv", encoding="utf-8-sig")

hit = P[P.Y19_CooC == "O"]
print("\n" + "=" * 96); print("### 판정"); print("=" * 96)
if len(hit):
    d_, a_ = hit.index[0]
    base = P.loc[(0.5, 5.0)]
    print(f"  ★ Y19 CooC 가 처음 잡히는 설정: -d {d_} -a {a_}")
    print(f"    그때 히트 수 BL21 {int(hit.iloc[0]['BL21'])} / "
          f"MG1655 {int(hit.iloc[0]['MG1655'])} / Y19 {int(hit.iloc[0]['Y19'])}")
    print(f"    기본(0.5/5.0) 대비 {int(hit.iloc[0]['BL21'])/max(int(base['BL21']),1):.1f}배")
    print("\n  이것이 이 방법의 실제 민감도다. Track C 를 다시 돌린다면 이 설정이어야")
    print("  하고, 기존 80개 목록은 '직교체도 못 잡는 설정'에서 나온 것임을 명시할 것.")
    print("  히트가 수백~수천으로 늘면 그다음은 배위 등급·pLDDT 로 줄이는 단계다.")
else:
    print("  ★ 넓힌 범위 전체에서 Y19 CooC 가 안 잡힌다.")
    print("    허용폭 문제가 아니다. S2 의 (1) 쌍 거리를 볼 것 —")
    print("    이 모델의 Cys 배치가 질의 기하와 다르다는 뜻이다.")
    print("    그 경우 Folddisco 로는 이 자리를 찾을 수 없고, Foldseek(폴드) 축으로")
    print("    가는 것이 맞다. 실제로 대조군을 통과한 것은 Foldseek 뿐이었다.")
```

---

## CELL S4 — 원인 확인 + 질의 프레임 교체

```python
# =============================================================================
# CELL S4 | S2·S3 가 말해준 것
#   Y19 CooC 의 CXC 는 제자리에 있다 (A115,A117 — Ch 의 112/114 에서 +3).
#   자기 질의로 rmsd 0.000 에 잡힌다. 서열도 폴드도 문제가 아니다.
#   다른 것은 거리다.  Ch CooC1 3.58 Å  vs  Y19 CooC 모델 8.16 Å.
#   AlphaFold 가 그 자리를 '열린(apo)' 상태로 예측했다.
#
#   이것이 뜻하는 바는 Track C 전체에 걸린다.
#     질의 = 결정구조의 닫힌 기하
#     대상 = 전부 AlphaFold 모델 (금속 없는 상태를 예측한다)
#   닫힌 기하로 열린 모델들을 훑고 있었다. 그러면 80개 히트는
#   'Ni 자리 기하를 가진 단백질' 이 아니라
#   'AlphaFold 가 우연히 Cys 쌍을 닫아 예측한 단백질' 이다.
#   프로테옴의 0.5~0.7% 만 나온 것도, 균주별로 들쭉날쭉한 것도 이걸로 설명된다.
#
#   (a) 거리가 원인인지 확인한다. 8.16 Å 이 걸리려면 -d >= 4.6 이어야 한다.
#       S3 은 3.0 까지만 봤다. 거기서 잡히면 설명이 맞고, 안 잡히면 틀린 것이다.
#   (b) 질의를 AlphaFold 프레임으로 바꾼다. Y19 CooC 의 AF 모델을 질의로 쓰면
#       질의와 대상이 같은 종류의 구조가 되어 상태 불일치가 사라진다.
# =============================================================================
Y19_CYS = (115, 117)          # S2 에서 확인한 CXC
OBS_D   = 8.16                # 그 쌍의 SG–SG (AF 모델)
REF_D   = 3.58                # Ch CooC1 (3kji)

# ---------- (a) -d 를 실제로 필요한 데까지 넓힌다 ----------
print("=" * 96); print("### (a) 거리가 원인인가 — -d 를 4.6 이상으로"); print("=" * 96)
print(f"  필요 조건: |{OBS_D} - {REF_D}| = {OBS_D-REF_D:.2f} Å 이내를 봐줘야 한다")
rows = []
for d_, a_ in [(4.0, 40.0), (5.0, 45.0), (6.0, 60.0)]:
    o = OUT/f"wide_d{d_}_a{a_}_Y19.tsv"
    if not (o.exists() and o.stat().st_size):
        sh(f'"{FD}" query -p "{COOC1_PDB}" -q {QRES} -i "{IDX["Y19"]}" -t {THREADS} '
           f'-d {d_} -a {a_} --per-structure --header --sort-by idf --top 50000 -o "{o}"',
           quiet=True)
    n, found = 0, False
    if o.exists() and o.stat().st_size:
        dd = pd.read_csv(o, sep="\t"); dd.columns = [c.strip().lstrip("#") for c in dd.columns]
        n = len(dd); found = dd.tid.astype(str).str.contains(Y19_COOC_UNI, na=False).any()
    rows.append({"-d": d_, "-a": a_, "Y19 히트": n, "CooC": "O" if found else "."})
    print(f"  d={d_:<4} a={a_:<5} Y19 {n:7d}개   CooC {'O' if found else 'X'}")
W = pd.DataFrame(rows)
if (W.CooC == "O").any():
    d0 = W[W.CooC == "O"].iloc[0]
    print(f"\n  ★ -d {d0['-d']} 에서 잡힌다. 거리가 원인이라는 설명이 맞다.")
    print(f"    다만 그때 히트가 {int(d0['Y19 히트'])}개다 — 이 허용폭은 선별력이 없다.")
    print("    즉 '넓히면 된다' 가 아니라 '이 질의 방식이 맞지 않는다' 가 결론이다.")
else:
    print("\n  ★ 여기서도 안 잡힌다. 거리만의 문제가 아니다 —")
    print("    각도나 다른 기하 특징도 함께 어긋나 있다는 뜻이다.")

# ---------- (b) 질의를 AlphaFold 프레임으로 ----------
print("\n" + "=" * 96); print("### (b) AlphaFold 모델을 질의로 — 상태를 맞춘다"); print("=" * 96)
print(f"  질의: {Y19_COOC_GB} 의 A{Y19_CYS[0]},A{Y19_CYS[1]}  (SG–SG {OBS_D} Å)")
print("  대상도 전부 AlphaFold 모델이므로, 같은 상태끼리 비교하게 된다.")
QA = f"A{Y19_CYS[0]},A{Y19_CYS[1]}"
rows = []
for s, idx in IDX.items():
    o = OUT/f"afframe_{s}.tsv"
    if not (o.exists() and o.stat().st_size):
        sh(f'"{FD}" query -p "{Y19_COOC_FILE}" -q {QA} -i "{idx}" -t {THREADS} '
           f'--per-structure --header --sort-by idf --top 20000 -o "{o}"', quiet=True)
    n = 0; d = None
    if o.exists() and o.stat().st_size:
        d = pd.read_csv(o, sep="\t"); d.columns = [c.strip().lstrip("#") for c in d.columns]
        n = len(d)
    rows.append({"strain": s, "hits": n})
    print(f"  {s:8s} {n:6d}개")
AF = pd.DataFrame(rows)
AF.to_csv(TBL/"folddisco_afframe_query.csv", index=False, encoding="utf-8-sig")

# 기존 질의(닫힌 기하)와 겹치는가 — 겹치지 않으면 서로 다른 집합을 보고 있었다는 뜻
print("\n  기존 질의(3kji 닫힌 기하)와의 비교")
FDD = BASE/"result"/"folddisco"
for s in IDX:
    old = FDD/f"metal_run01_{s}.tsv"; new = OUT/f"afframe_{s}.tsv"
    if not (old.exists() and new.exists() and new.stat().st_size): continue
    a = pd.read_csv(old, sep="\t"); a.columns = [c.strip().lstrip("#") for c in a.columns]
    b = pd.read_csv(new, sep="\t"); b.columns = [c.strip().lstrip("#") for c in b.columns]
    A_ = {Path(str(x)).stem for x in a.tid}; B_ = {Path(str(x)).stem for x in b.tid}
    print(f"  {s:8s} 닫힌 {len(A_):5d} / 열린 {len(B_):5d} / 겹침 {len(A_&B_):5d}")

print("\n" + "=" * 96); print("### 읽는 법"); print("=" * 96)
print("  겹침이 거의 없다면, 두 질의는 전혀 다른 단백질 집합을 보고 있었다는 뜻이다.")
print("  그러면 '어느 쪽이 맞느냐' 가 아니라 '무엇을 찾고 있었느냐' 의 문제다.")
print("    닫힌 기하 = 금속이 든 상태의 자리        (결정구조에서 나온 것)")
print("    열린 기하 = AlphaFold 가 예측한 apo 자리 (대상과 같은 상태)")
print("  Ni 을 실제로 잡고 있는 단백질을 찾는 것이 목적이라면 전자가 옳다.")
print("  그러나 대상이 전부 AF 모델이면 전자로는 대조군조차 못 찾는다 — 방금 봤다.")
print("  ★ 그래서 Folddisco 를 이 목적에 쓰려면 대상 구조를 바꿔야 한다")
print("    (금속 결합 상태로 재예측하거나, Boltz 로 Ni 을 넣고 그 구조를 색인하거나).")
print("    그전까지 이 축의 '없음' 은 전부 미검출이고, 판단은 Foldseek 쪽에 둔다.")
```
