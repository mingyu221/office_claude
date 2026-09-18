# Folddisco — 질의를 예측구조 프레임으로 바꾼다

## 왜

지금까지의 Track C 는 **질의와 대상의 상태가 달랐다.**

| | 무엇이었나 | 상태 |
|---|---|---|
| 질의 | 3kji 결정구조 (Ch CooC1) | **금속이 들어가 닫힌 holo 형태** |
| 대상 | 세 균주 AlphaFold 모델 전부 | **금속이 없는 열린 apo 형태** |

금속이 붙으면 배위 잔기가 안쪽으로 모인다. 그 닫힌 기하로 열린 모델을 훑었으니
안 잡히는 게 당연하다. 실측으로도 확인됐다 — Y19 CooC 의 CXC 는 제자리에 있고
자기 질의로 rmsd 0.000 에 잡히는데 **SG–SG 가 8.16 Å**, Ch CooC1 결정은 **3.58 Å** 다.
`-d` 를 3.0 까지 넓혀도 4.58 Å 을 못 메운다.

## 무엇을 할 것인가

**Ch CooC1 의 AlphaFold 예측 모델**을 질의로 쓴다. 그러면 질의와 대상이
같은 예측 프레임에 놓인다.

검증은 대조군 하나로 한다 — **Y19 의 CooC(`A0A059VK30`)가 잡히는가.**
같은 과의 직교체이고 Foldseek 이 tm 0.938 로 잡는 단백질이다. 이게 안 잡히면
설정을 바꿔도 Track C 는 못 쓴다.

F1 → F2 → F3 → F4 → F5 순으로 실행.

---

## CELL F1 — 설정

```python
# =============================================================================
# CELL F1 | 독립 실행용 설정
# =============================================================================
import re, json, math, subprocess, itertools, urllib.request
from pathlib import Path
import pandas as pd

BASE  = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS = Path("/mnt/af2results/mingyu")
TBL   = BASE/"result"/"table"
OUT   = BASE/"result"/"folddisco"/"afframe"; OUT.mkdir(parents=True, exist_ok=True)
THREADS = 8

FD  = TOOLS/"folddisco"/"bin"/"folddisco"
FS  = TOOLS/"foldseek"/"bin"/"foldseek"
IDX = {"BL21":   TOOLS/"database"/"bacteriaDB"/"folddisco"/"UP000503272",
       "Y19":    TOOLS/"database"/"bacteriaDB"/"folddisco"/"UP000034085",
       "MG1655": TOOLS/"database"/"folddisco"/"ecoli_folddisco"/"e_coli_folddisco"}
STRUCT = {"BL21": TOOLS/"database"/"bacteriaDB"/"structures_UP000503272",
          "Y19":  TOOLS/"database"/"bacteriaDB"/"structures_UP000034085"}

CRYSTAL   = TOOLS/"workspace"/"seek_ni_insertase"/"input"/"3kji.pdb"
PDB_ID    = "3KJI"
CRY_CHAIN = "A"
CRY_RES   = (112, 114)         # 결정구조에서의 Ni 배위 Cys 쌍

# 대조군 — Y19 의 CooC
Y19_GB, Y19_UNI = "AHZ96930.1", "A0A059VK30"
Y19_FILE = STRUCT["Y19"]/f"AF-{Y19_UNI}-F1-model_v6.cif"

# AlphaFold 모델을 손으로 넣고 싶을 때 여기에 경로를 준다 (없으면 F2 가 받아온다)
AF_MANUAL = None

def sh(cmd, quiet=False):
    r = subprocess.run(cmd, shell=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if not quiet: print(r.stdout.rstrip())
    return r.stdout

def get(url, timeout=30):
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read()

print("=" * 96); print("### 경로"); print("=" * 96)
for n, p in [("folddisco", FD), ("foldseek", FS), ("3kji", CRYSTAL),
             ("Y19 CooC", Y19_FILE), ("Y19 index", IDX["Y19"]),
             ("BL21 index", IDX["BL21"]), ("MG1655 index", IDX["MG1655"]), ("출력", OUT)]:
    print(f"  {'O' if (p.exists() or p.parent.exists()) else 'X'}  {n:14s} {p}")
```

---

## CELL F2 — Ch CooC1 의 AlphaFold 모델을 확보한다

```python
# =============================================================================
# CELL F2 | 결정구조와 같은 단백질의 예측 모델을 가져온다
#
#   경로는 셋. 위에서부터 시도한다.
#     (1) AF_MANUAL 에 경로를 줬으면 그걸 쓴다
#     (2) RCSB API 로 3KJI → UniProt 을 얻고 AlphaFold DB 에서 받는다
#     (3) 둘 다 안 되면 결정구조 서열을 FASTA 로 떨궈 둔다.
#         ColabFold 로 직접 접어 AF_MANUAL 에 경로를 주고 다시 돌리면 된다.
#
#   ★ 중요: 목표는 '예측 모델' 이면 된다는 것이 아니다. **대상 DB 와 같은 방식으로
#     예측된 모델** 이어야 한다. 대상은 AlphaFold2 계열이므로 AFDB 모델이 1순위,
#     ColabFold(역시 AF2) 가 2순위다.
# =============================================================================
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}

def chain_seq(path, chain=None):
    """(서열, [잔기번호...]) — CA 원자 순서 그대로."""
    seq, nums, seen = [], [], set()
    for l in open(path, errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")): continue
        if l[12:16].strip() != "CA": continue
        ch, rs, rn = l[21], l[22:26].strip(), l[17:20].strip()
        if chain and ch != chain: continue
        if not rs.lstrip("-").isdigit() or (ch, rs) in seen: continue
        seen.add((ch, rs)); seq.append(AA3.get(rn, "X")); nums.append(int(rs))
    return "".join(seq), nums

CRY_SEQ, CRY_NUM = chain_seq(CRYSTAL, CRY_CHAIN)
print(f"  결정구조 {CRY_CHAIN} 사슬 {len(CRY_SEQ)} 잔기, 번호 {CRY_NUM[0]}–{CRY_NUM[-1]}")
for r in CRY_RES:
    i = CRY_NUM.index(r) if r in CRY_NUM else None
    print(f"    {r}번 = {CRY_SEQ[i] if i is not None else '?'}"
          + (f"   주변 ±7: {CRY_SEQ[max(0,i-7):i+8]}" if i is not None else ""))

AF_PATH, AF_SRC, UNI = None, None, None
if AF_MANUAL:
    AF_PATH, AF_SRC = Path(AF_MANUAL), "수동 지정"
else:
    try:
        u = (f"https://data.rcsb.org/rest/v1/core/polymer_entity/{PDB_ID}/1")
        j = json.loads(get(u))
        ids = j.get("rcsb_polymer_entity_container_identifiers", {}) \
               .get("reference_sequence_identifiers", [])
        UNI = next((r["database_accession"] for r in ids
                    if r.get("database_name", "").lower().startswith("uniprot")), None)
        print(f"\n  RCSB → UniProt: {UNI}")
    except Exception as e:
        print(f"\n  RCSB 조회 실패: {e}")
    if UNI:
        for v in (4, 6):
            dst = OUT/f"AF-{UNI}-F1-model_v{v}.pdb"
            if dst.exists() and dst.stat().st_size:
                AF_PATH, AF_SRC = dst, f"AFDB v{v} (캐시)"; break
            try:
                dst.write_bytes(get(f"https://alphafold.ebi.ac.uk/files/AF-{UNI}-F1-model_v{v}.pdb"))
                AF_PATH, AF_SRC = dst, f"AFDB v{v}"; break
            except Exception as e:
                print(f"  AFDB v{v} 실패: {e}")

if AF_PATH and AF_PATH.exists():
    print(f"\n  ★ 예측 모델: {AF_PATH.name}  ({AF_SRC})")
else:
    fa = OUT/"cooc1_crystal.fasta"
    fa.write_text(f">CooC1_{PDB_ID}_{CRY_CHAIN}\n{CRY_SEQ}\n")
    print(f"\n  ⚠ 예측 모델을 못 받았다. 결정구조 서열을 떨궜다: {fa}")
    print("    ColabFold 로 접은 뒤 F1 의 AF_MANUAL 에 그 경로를 주고 다시 실행할 것.")
    print("    예)  colabfold_batch --num-recycle 3 "
          f"{fa} {OUT/'colabfold'}")
```

---

## CELL F3 — 두 프레임의 기하를 나란히 잰다

```python
# =============================================================================
# CELL F3 | 예측 모델에서 같은 Cys 쌍을 찾고, 결정구조와 거리를 비교한다
#
#   번호가 같으리라 가정하지 않는다. 결정구조의 Cys 주변 서열로 예측 모델을
#   정렬해 대응 잔기를 찾는다. 예측 모델에는 결정에서 안 보이던 앞뒤 잔기가
#   들어 있어 번호가 밀리는 일이 흔하다.
#
#   보려는 것은 하나다 — **예측 모델의 SG–SG 가 결정보다 벌어져 있는가.**
#   벌어져 있다면 'AF 는 이 자리를 열린 상태로 예측한다' 가 Ch CooC1 자신에게서도
#   확인되는 것이고, Y19 CooC 의 8.16 Å 이 이상값이 아니라 규칙이 된다.
# =============================================================================
assert AF_PATH and AF_PATH.exists(), "F2 에서 예측 모델을 먼저 확보할 것"

def sg_of(path, chain, rs):
    for l in open(path, errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")): continue
        if l[12:16].strip() != "SG": continue
        if chain and l[21] != chain: continue
        if l[22:26].strip() == str(rs):
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return None

AF_SEQ, AF_NUM = chain_seq(AF_PATH)
print(f"  예측 모델 {len(AF_SEQ)} 잔기, 번호 {AF_NUM[0]}–{AF_NUM[-1]}")

# --- 결정구조 Cys 주변 창으로 대응 위치를 찾는다 ---
def map_res(r):
    """결정구조 잔기번호 r 에 대응하는 예측 모델의 잔기번호."""
    if r not in CRY_NUM: return None, "결정구조에 없음"
    i = CRY_NUM.index(r)
    for w in (9, 7, 5, 3, 2):
        lo  = max(0, i - w)
        win = CRY_SEQ[lo:i+w+1]
        off = i - lo               # ★ 창 안에서 r 의 위치. 글자를 찾으면 안 된다 —
                                   #   창에 같은 아미노산이 둘이면 첫 번째로 끌려간다
                                   #   (C112/C114 가 정확히 그 경우였다)
        j = AF_SEQ.find(win)
        if j >= 0 and AF_SEQ.count(win) == 1:
            return AF_NUM[j+off], f"창 ±{w} 로 유일 매칭"
    return None, "유일 매칭 실패"

AF_RES, notes = [], []
for r in CRY_RES:
    m, why = map_res(r); AF_RES.append(m); notes.append(why)
    aa = AF_SEQ[AF_NUM.index(m)] if m in AF_NUM else "?"
    cr = CRY_SEQ[CRY_NUM.index(r)]
    ok = "OK" if aa == cr else f"★ 불일치 (결정 {cr} vs 예측 {aa})"
    print(f"    결정 {r}({cr}) → 예측 {m}({aa})   {why}   {ok}")
    if m in AF_NUM:
        k = AF_NUM.index(m)
        print(f"        예측 주변 ±7: {AF_SEQ[max(0,k-7):k+8]}")

# --- 매핑 검산 — 여기서 걸러야 퇴화 질의로 0개가 나오는 일이 없다 ---
if len(set(AF_RES)) < len(AF_RES):
    raise SystemExit(f"매핑이 겹쳤다: {AF_RES}. 같은 잔기를 두 번 질의하게 된다.")
for r, m in zip(CRY_RES, AF_RES):
    if m is not None and AF_SEQ[AF_NUM.index(m)] != CRY_SEQ[CRY_NUM.index(r)]:
        raise SystemExit(f"결정 {r} 과 예측 {m} 의 아미노산이 다르다. 매핑을 확인할 것.")

# --- 매핑이 안 되면 모델의 CXC 후보를 전부 보여준다 ---
if any(m is None for m in AF_RES):
    print("\n  매핑 실패. 예측 모델의 CYS 쌍을 전부 잰다:")
    cys = [n for n, a in zip(AF_NUM, AF_SEQ) if a == "C"]
    print(f"    CYS {len(cys)}개: {cys}")
    for a, b in itertools.combinations(cys, 2):
        pa, pb = sg_of(AF_PATH, None, a), sg_of(AF_PATH, None, b)
        if pa and pb and abs(a-b) <= 6:
            print(f"      {a}–{b}  SG–SG {math.dist(pa, pb):6.2f} Å")
    raise SystemExit("AF_RES 를 손으로 정하고 다시 실행할 것")

# --- 거리 비교 ---
print("\n" + "=" * 96); print("### 같은 Cys 쌍, 두 프레임"); print("=" * 96)
d_cry = math.dist(sg_of(CRYSTAL, CRY_CHAIN, CRY_RES[0]), sg_of(CRYSTAL, CRY_CHAIN, CRY_RES[1]))
d_af  = math.dist(sg_of(AF_PATH, None, AF_RES[0]),       sg_of(AF_PATH, None, AF_RES[1]))
print(f"  결정구조 {CRY_CHAIN}{CRY_RES[0]}–{CRY_CHAIN}{CRY_RES[1]}   SG–SG {d_cry:6.2f} Å   (금속 결합, 닫힘)")
print(f"  예측 모델 {AF_RES[0]}–{AF_RES[1]}   SG–SG {d_af:6.2f} Å")
print(f"  차이 {d_af - d_cry:+6.2f} Å")

# Y19 CooC 도 같이 (S2 에서 8.16 Å 으로 나왔던 값)
if Y19_FILE.exists():
    ys, yn = chain_seq(Y19_FILE)
    ycys = [n for n, a in zip(yn, ys) if a == "C"]
    pairs = [(a, b) for a, b in itertools.combinations(ycys, 2) if abs(a-b) == 2]
    for a, b in pairs:
        pa, pb = sg_of(Y19_FILE, None, a), sg_of(Y19_FILE, None, b)
        if pa and pb:
            D_Y19 = math.dist(pa, pb)
            print(f"  Y19 CooC {a}–{b}          SG–SG {D_Y19:6.2f} Å")

# 질의가 대조군에 닿으려면 허용폭이 얼마여야 하는가 — 실측에서 바로 나온다
GAP = abs(D_Y19 - d_af) if "D_Y19" in dir() else None
if GAP:
    print(f"\n  ※ 예측 프레임 질의와 Y19 CooC 의 간격: {GAP:.2f} Å")
    print(f"    (결정 프레임이었을 때는 {abs(D_Y19 - d_cry):.2f} Å 이었다)")
    print(f"    → 이 질의로 대조군을 잡으려면 -d 를 최소 {GAP:.1f} 이상 줘야 한다.")

print("\n  읽는 법")
print("    예측 모델도 벌어져 있으면 → AF 가 이 자리를 열린 상태로 예측한다는 것이")
print("      Ch CooC1 자신에게서 확인된 것이다. F4 의 질의가 맞는 방향이다.")
print("    예측 모델이 결정만큼 닫혀 있으면 → 프레임 문제가 아니라 Y19 쪽 모델의")
print("      개별 문제다. 그때는 이 노트북으로 안 풀린다.")
```

---

## CELL F4 — 예측 프레임 질의 + 대조군

```python
# =============================================================================
# CELL F4 | 같은 질의를 두 프레임으로 돌려 나란히 본다
#   판정은 하나다 — 예측 프레임 질의가 Y19 CooC 를 잡는가.
# =============================================================================
QCRY = f"{CRY_CHAIN}{CRY_RES[0]},{CRY_CHAIN}{CRY_RES[1]}"
AF_CH = open(AF_PATH, errors="ignore").readline()  # 사슬 ID 확인용
_ch = None
for l in open(AF_PATH, errors="ignore"):
    if l.startswith("ATOM"): _ch = l[21]; break
QAF  = f"{_ch}{AF_RES[0]},{_ch}{AF_RES[1]}"
assert AF_RES[0] != AF_RES[1], f"질의 잔기가 같다: {AF_RES}. F3 의 매핑을 다시 볼 것"
print(f"  결정 프레임 질의: {CRYSTAL.name}  -q {QCRY}")
print(f"  예측 프레임 질의: {AF_PATH.name}  -q {QAF}")

def query(tag, qpath, qres, strain, d=0.5, a=5.0):
    # 파일명에 질의 잔기를 넣는다. 안 넣으면 질의를 고쳐도 옛 결과를 재사용한다
    # (A112,A112 퇴화 질의의 0개 결과가 그대로 살아남는 사고가 있었다)
    o = OUT/f"{tag}_{qres.replace(',', '-')}_d{d}_a{a}_{strain}.tsv"
    if not (o.exists() and o.stat().st_size):
        sh(f'"{FD}" query -p "{qpath}" -q {qres} -i "{IDX[strain]}" -t {THREADS} '
           f'-d {d} -a {a} --per-structure --header --sort-by idf --top 20000 -o "{o}"',
           quiet=True)
    if not (o.exists() and o.stat().st_size): return set(), o
    dd = pd.read_csv(o, sep="\t"); dd.columns = [c.strip().lstrip("#") for c in dd.columns]
    return set(dd.tid.astype(str)), o

# F3 에서 잰 간격 너머까지 반드시 훑는다. 앞서 1.5 에서 끊어 놓고
# "안 잡힌다" 로 읽은 적이 있다 — 실제 간격은 2.86 Å 이었다.
GRID = [(0.5, 5.0), (1.0, 10.0), (1.5, 15.0), (2.0, 20.0), (3.0, 30.0), (4.0, 40.0)]
try:
    if GAP and max(d for d, _ in GRID) < GAP + 1.0:
        GRID.append((round(GAP + 1.0, 1), round((GAP + 1.0) * 10, 1)))
        print(f"  간격 {GAP:.2f} Å 를 넘기려고 -d {GRID[-1][0]} 를 추가했다")
except NameError:
    pass
rows, HITS = [], {}
for d_, a_ in GRID:
    for tag, qp, qr in [("cry", CRYSTAL, QCRY), ("af", AF_PATH, QAF)]:
        for s in ["BL21", "MG1655", "Y19"]:
            h, o = query(tag, qp, qr, s, d_, a_)
            HITS[(tag, d_, a_, s)] = h
            ok = any(Y19_UNI in t for t in h) if s == "Y19" else None
            rows.append({"프레임": "결정" if tag == "cry" else "예측",
                         "-d": d_, "-a": a_, "strain": s, "hits": len(h),
                         "Y19_CooC": ("O" if ok else "X") if s == "Y19" else ""})
            print(f"  {'결정' if tag=='cry' else '예측'}  d={d_:<4} a={a_:<5} {s:8s} {len(h):6d}개"
                  + (f"   CooC {'O' if ok else 'X'}" if s == "Y19" else ""))

R = pd.DataFrame(rows)
print("\n" + "=" * 96); print("### 두 프레임 비교"); print("=" * 96)
print(R.pivot_table(index=["-d", "-a", "strain"], columns="프레임",
                    values="hits").to_string())
R.to_csv(TBL/"folddisco_afframe_compare.csv", index=False, encoding="utf-8-sig")

print("\n" + "=" * 96); print("### 대조군 판정"); print("=" * 96)
PASS = [(d_, a_) for d_, a_ in GRID
        if any(Y19_UNI in t for t in HITS[("af", d_, a_, "Y19")])]
if PASS:
    d_, a_ = PASS[0]
    print(f"  ★ 예측 프레임 질의가 Y19 CooC 를 잡는다 — 처음 잡히는 설정 -d {d_} -a {a_}")
    nb, nm, ny = (len(HITS[('af', d_, a_, s)]) for s in ("BL21", "MG1655", "Y19"))
    print(f"    그때 히트 BL21 {nb} / MG1655 {nm} / Y19 {ny}")
    base = len(HITS[("af", 0.5, 5.0, "BL21")]) or 1
    print(f"    기본(0.5/5.0) 대비 {nb/base:.1f}배")
    if d_ >= 2.0:
        print("\n    ★ 단, 이 허용폭은 매우 느슨하다. -d 3.0 이면 사실상")
        print("      'Cys 두 개가 2.4~8.4 Å 안에 있다' 를 묻는 것이다. 대조군을")
        print("      통과해도 모티프로서의 변별력은 거의 없다고 봐야 한다.")
        print("      F5 의 균주 비율이 세 균주 모두 같게 나오면 그 확인이다.")
    print("    → F5 로 기존 목록과 비교할 것.")
    BEST = (d_, a_)
else:
    print("  넓힌 범위 전체에서도 Y19 CooC 가 안 잡힌다.")
    print("    프레임을 맞춰도 안 되면 원인은 리간드 상태가 아니다.")
    print("    예측 모델끼리도 이 자리의 기하가 보존되지 않는다는 뜻이고,")
    print("    그러면 2잔기 거리 모티프로는 종간 비교가 원리상 성립하지 않는다.")
    print("    → F6 의 Foldseek(폴드 전체) 축으로 간다.")
    BEST = None
```

---

## CELL F5 — 살아났다면, 기존 목록과 무엇이 다른가

```python
# =============================================================================
# CELL F5 | 예측 프레임 히트 목록을 저장하고 결정 프레임 80개와 비교한다
#   묻는 것 셋.
#     (1) 기존 80개 중 몇 개가 살아남는가 — 기존 목록이 인공물이었는지
#     (2) 새로 들어온 것은 무엇인가 — 놓치고 있던 것
#     (3) 균주 차이가 생기는가 — 기존엔 비율이 세 균주 같았다
# =============================================================================
if BEST is None:
    print("  F4 대조군이 실패했다. 이 셀은 의미가 없다.")
else:
    d_, a_ = BEST
    xw = pd.read_csv(TBL/"id_crosswalk_struct_to_genbank.csv") \
         if (TBL/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
    XW = ({Path(str(k)).stem: str(v).split(",")[0]
           for k, v in xw.set_index("tid")["protein"].to_dict().items()}
          if len(xw) else {})
    def pid(t): return XW.get(Path(str(t)).stem, Path(str(t)).stem)

    PROT = {"BL21": 4110, "MG1655": 4300, "Y19": 5325}
    print("=" * 96); print(f"### 예측 프레임 (-d {d_} -a {a_})"); print("=" * 96)
    NEW = {}
    for s in ["BL21", "MG1655", "Y19"]:
        af  = {pid(t) for t in HITS[("af",  d_, a_, s)]}
        cry = {pid(t) for t in HITS[("cry", 0.5, 5.0, s)]}
        NEW[s] = af
        print(f"\n  [{s}]  예측 {len(af)}개 ({len(af)/PROT[s]*100:.2f}%)  "
              f"vs 결정 기본 {len(cry)}개")
        print(f"    양쪽 다: {len(af & cry)}   예측에만: {len(af - cry)}   "
              f"결정에만: {len(cry - af)}")
        if af - cry:
            print("    새로 들어온 것(앞 12개): " + ", ".join(sorted(af - cry)[:12]))

    pd.DataFrame([{"strain": s, "protein": p} for s, v in NEW.items() for p in sorted(v)]) \
      .to_csv(TBL/"folddisco_afframe_hits.csv", index=False, encoding="utf-8-sig")
    print(f"\n저장: {TBL/'folddisco_afframe_hits.csv'}")

    print("\n" + "=" * 96); print("### 균주 비율"); print("=" * 96)
    for s in ["BL21", "MG1655", "Y19"]:
        print(f"  {s:8s} {len(NEW[s]):5d} / {PROT[s]}  = {len(NEW[s])/PROT[s]*100:5.2f}%")
    print("\n  기존 결정 프레임에서는 세 균주 비율이 같았다(0.5~0.7%). 그건 질의가")
    print("  'CooC 같은 자리' 가 아니라 'AF 가 우연히 닫아 예측한 자리' 를 세고")
    print("  있었기 때문이다. 예측 프레임에서도 비율이 같으면 균주 변별력은")
    print("  여전히 없는 것이고, 갈리면 그때 처음으로 의미 있는 차이다.")
```

---

## CELL F6 — Foldseek 으로 같은 질문을 던진다

```python
# =============================================================================
# CELL F6 | 모티프가 안 되면 폴드로 간다
#
#   Folddisco 는 잔기 2~4개의 국소 기하를 본다. 그래서 리간드 유무로 그 자리가
#   열리고 닫히는 것에 그대로 휘둘린다. Foldseek 은 구조 전체를 정렬하므로
#   그 자리가 조금 벌어져도 폴드는 같게 잡는다.
#   실제로 대조군을 통과한 것은 Foldseek 뿐이었다 (CooC1 → Y19 CooC tm 0.938).
#
#   여기서 두 가지를 같이 본다.
#     (1) 결정 프레임 / 예측 프레임 질의가 Foldseek 에서도 갈리는가.
#         안 갈리면 — 그것이 Foldseek 을 쓰는 이유의 증거다.
#     (2) 세 균주에서 CooC 과(科)에 드는 것이 무엇인가. 이게 Track C 를
#         대신할 후보 목록이다.
#
#   ※ MG1655 는 구조 디렉터리가 없고 folddisco 인덱스만 있다. 없으면 건너뛰고
#     F7 에서 서열로 묻는다.
# =============================================================================
FFMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits,prob,alntmscore,lddt"
TM_FOLD, TM_NEAR, TM_CORE = 0.50, 0.70, 0.90
FSD = OUT/"foldseek"; FSD.mkdir(parents=True, exist_ok=True)

# 결정구조는 A 사슬만 떼어 질의로 쓴다 (B 사슬·리간드가 섞이면 정렬이 흐려진다)
CRY_A = OUT/"3kji_chainA.pdb"
if not CRY_A.exists():
    with open(CRYSTAL, errors="ignore") as fh, open(CRY_A, "w") as o:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")) and l[21] == CRY_CHAIN: o.write(l)
        o.write("END\n")
print(f"  결정 질의 {CRY_A.name} / 예측 질의 {AF_PATH.name}")

STRUCT_ALL = dict(STRUCT)
_mg = TOOLS/"database"/"protein_list"/"structures_MG1655"
for cand in [_mg, TOOLS/"database"/"structures_mg1655",
             TOOLS/"database"/"bacteriaDB"/"structures_MG1655"]:
    if cand.is_dir(): STRUCT_ALL["MG1655"] = cand; break
print(f"  구조 DB: {', '.join(STRUCT_ALL)}"
      + ("" if "MG1655" in STRUCT_ALL else "   (MG1655 구조 없음 → F7 에서 서열로)"))

def fs(tag, qpath, strain):
    o = FSD/f"{tag}_vs_{strain}.m8"
    if not (o.exists() and o.stat().st_size):
        sh(f'"{FS}" easy-search "{qpath}" "{STRUCT_ALL[strain]}" "{o}" "{FSD}/tmp_{strain}" '
           f'--format-output "{FFMT}" -e 10 --max-seqs 2000 --exact-tmscore 1 '
           f'--threads {THREADS}', quiet=True)
    if not (o.exists() and o.stat().st_size): return pd.DataFrame()
    d = pd.read_csv(o, sep="\t", names=FFMT.split(","))
    d["stem"] = d.target.apply(lambda x: Path(str(x)).stem)
    return d.sort_values("alntmscore", ascending=False).drop_duplicates("stem")

xw = pd.read_csv(TBL/"id_crosswalk_struct_to_genbank.csv") \
     if (TBL/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
XW = ({Path(str(k)).stem: str(v).split(",")[0]
       for k, v in xw.set_index("tid")["protein"].to_dict().items()} if len(xw) else {})
def pid(stem): return XW.get(stem, stem)

RES_FS, SETS = [], {}
for tag, qp in [("cry", CRY_A), ("af", AF_PATH)]:
    for s in STRUCT_ALL:
        d = fs(tag, qp, s)
        if not len(d):
            print(f"  {tag:3s} → {s:8s} 결과 없음"); continue
        fold = d[d.alntmscore >= TM_FOLD]
        SETS[(tag, s)] = {pid(x) for x in fold.stem}
        top = d.iloc[0]
        ctrl = ""
        if s == "Y19":
            hit = d[d.stem.str.contains(Y19_UNI, na=False)]
            ctrl = (f"   대조군 {Y19_GB} tm {hit.iloc[0].alntmscore:.3f}"
                    if len(hit) else "   ★ 대조군 없음")
        print(f"  {tag:3s} → {s:8s} 폴드(≥{TM_FOLD}) {len(fold):4d}개   "
              f"최고 {pid(top.stem)} tm {top.alntmscore:.3f}{ctrl}")
        for _, r in fold.iterrows():
            RES_FS.append({"frame": "결정" if tag == "cry" else "예측", "strain": s,
                           "protein": pid(r.stem), "tm": round(float(r.alntmscore), 3),
                           "fident": round(float(r.fident), 3), "qcov": round(float(r.qcov), 3),
                           "tier": ("핵심" if r.alntmscore >= TM_CORE else
                                    "주변" if r.alntmscore >= TM_NEAR else "폴드만")})

F = pd.DataFrame(RES_FS)
F.to_csv(TBL/"cooc_fold_afframe.csv", index=False, encoding="utf-8-sig")
print(f"\n저장: {TBL/'cooc_fold_afframe.csv'}")

print("\n" + "=" * 96); print("### (1) 두 프레임이 Foldseek 에서 갈리는가"); print("=" * 96)
for s in STRUCT_ALL:
    a, b = SETS.get(("cry", s), set()), SETS.get(("af", s), set())
    if not a and not b: continue
    j = len(a & b) / max(len(a | b), 1)
    print(f"  {s:8s} 결정 {len(a):4d}  예측 {len(b):4d}  교집합 {len(a&b):4d}  자카드 {j:.3f}")
print("\n  자카드가 높으면 → Foldseek 은 리간드 유무에 휘둘리지 않는다.")
print("    Folddisco 가 0 이었던 것과 대비되는 지점이고, 축을 바꾸는 근거가 된다.")

print("\n" + "=" * 96); print("### (2) CooC 과 구성원 (예측 프레임, tier 별)"); print("=" * 96)
for s in STRUCT_ALL:
    sub = F[(F.frame == "예측") & (F.strain == s)]
    if not len(sub): continue
    print(f"\n  [{s}]  " + "  ".join(f"{t} {len(sub[sub.tier==t])}"
                                     for t in ["핵심", "주변", "폴드만"]))
    print(sub.sort_values("tm", ascending=False)
             .head(12)[["protein", "tm", "fident", "tier"]].to_string(index=False))
```

---

## CELL F7 — 찾아낸 과 구성원에게 MG1655 를 묻는다

```python
# =============================================================================
# CELL F7 | Foldseek 이 뽑은 BL21 구성원이 MG1655 에도 있는가
#   MG1655 는 구조 DB 가 없으므로 구조로는 못 묻는다. 서열로 묻는다.
#   프로테옴에 없으면 게놈까지 — 앞서 어노테이션 누락에 두 번 속았다.
# =============================================================================
import tempfile
FAA_BL21 = TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"bl21_db_match_qjz.faa"
FAA_MG   = TOOLS/"database"/"protein_list"/"mg1655_protein.faa"
GEN_MG   = BASE/"input"/"external"/"MG1655_U00096.3.fna"
SRCH     = BASE/"result"/"search"; SRCH.mkdir(parents=True, exist_ok=True)
MFMT = "query,target,fident,alnlen,evalue,bits,qlen,tlen,qcov,tcov,qstart,qend"

def read_faa(p):
    s, k = {}, None
    for l in open(p):
        if l.startswith(">"): k = l[1:].split()[0]; s[k] = []
        elif k: s[k].append(l.strip())
    return {k: "".join(v) for k, v in s.items()}

SEQB = read_faa(FAA_BL21) if FAA_BL21.exists() else {}
MEM  = sorted({p for p in F[(F.frame == "예측") & (F.strain == "BL21")].protein
               if p in SEQB})
print(f"  BL21 과 구성원 {len(MEM)}개에게 묻는다 (FASTA 에 있는 것만)")

QF = SRCH/"coocfold_bl21.faa"
QF.write_text("".join(f">{k}\n{SEQB[k]}\n" for k in MEM))

def msearch(q, tgt, out, extra=""):
    if not (out.exists() and out.stat().st_size):
        with tempfile.TemporaryDirectory() as td:
            sh(f'mmseqs easy-search "{q}" "{tgt}" "{out}" "{td}" '
               f'--format-output "{MFMT}" -s 7.5 -e 1e-3 --max-seqs 5 '
               f'--threads {THREADS} {extra}', quiet=True)
    if not (out.exists() and out.stat().st_size): return pd.DataFrame()
    return pd.read_csv(out, sep="\t", names=MFMT.split(","))

d = msearch(QF, FAA_MG, SRCH/"coocfold_vs_MG1655.m8")
PH = {}
if len(d):
    d = d[d.qcov >= 0.70].sort_values("bits", ascending=False).drop_duplicates("query")
    PH = {str(r.query): (float(r.fident), str(r.target)) for _, r in d.iterrows()}
print(f"  프로테옴 히트(qcov≥0.70) {len(PH)} / {len(MEM)}")

NOHIT = [m for m in MEM if m not in PH]
GH = {}
if NOHIT:
    gf = SRCH/"coocfold_nohit.faa"
    gf.write_text("".join(f">{k}\n{SEQB[k]}\n" for k in NOHIT))
    g = msearch(gf, GEN_MG, SRCH/"coocfold_vs_MG1655genome.m8", "--search-type 2")
    if len(g):
        # --search-type 2 의 alnlen 은 염기 단위다. 커버리지는 qstart/qend 로 낸다
        g["qcov2"] = (g.qend - g.qstart + 1) / g.qlen
        g = g[(g.fident >= 0.80) & (g.qcov2 >= 0.70)] \
              .sort_values("bits", ascending=False).drop_duplicates("query")
        GH = {str(r.query): float(r.fident) for _, r in g.iterrows()}
print(f"  프로테옴 미검출 {len(NOHIT)} 중 게놈에 있는 것 {len(GH)}")

TMX = F[(F.frame == "예측") & (F.strain == "BL21")].set_index("protein")["tm"].to_dict()
rows = []
for m in MEM:
    if m in PH:   v, why = ".", f"MG1655 에 있음 {PH[m][1]} fident {PH[m][0]:.3f}"
    elif m in GH: v, why = "~", f"게놈에만 있음 fident {GH[m]:.3f}"
    else:         v, why = "O", "프로테옴·게놈 둘 다 없음"
    rows.append({"protein": m, "tm_CooC": TMX.get(m), "MG1655": v, "근거": why})
C = pd.DataFrame(rows).sort_values(["MG1655", "tm_CooC"], ascending=[True, False])
C.to_csv(TBL/"cooc_fold_mg1655.csv", index=False, encoding="utf-8-sig")

print("\n" + "=" * 96); print("### CooC 폴드 × MG1655 부재"); print("=" * 96)
print(C.to_string(index=False))
print(f"\n저장: {TBL/'cooc_fold_mg1655.csv'}")
print("\n  O   폴드는 CooC 인데 MG1655 에 없다  ← 가장 강한 조합")
print("  ~   어노테이션 누락")
print("  .   MG1655 에도 있다")
print("\n  ※ BL21 에는 tm 0.9 이상의 진짜 CooC 이 없다는 것이 이미 확인됐다")
print("    (최고 ApbC tm 0.690). 그러니 여기서 O 가 나와도 'CooC 자체' 가 아니라")
print("    'CooC 과에 드는 BL21 단백질' 이다. 그 구분을 표에 유지할 것.")
```
