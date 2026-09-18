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
