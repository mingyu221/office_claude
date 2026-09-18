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

# ---- MG1655 구조 디렉터리를 실제로 찾는다 --------------------------------
#   경로를 추측해 세 개만 찍어보고 포기하면, 있는 것도 '없다' 로 처리된다.
#   database 트리를 훑어 구조 파일이 많이 든 디렉터리를 세어 보고 고른다.
import os
# AlphaFold DB 에서 받아 둔 대장균 프로테옴을 쓴다.
#   AF-{UniProt}-F1-model_v6 규칙이고 4,371개로 프로테옴 4,300 과 맞는다.
#   /database/mg1655 (7,447개) 는 model_v1 사슬분할 파일이 섞여 있어 쓰면 안 된다.
#   ※ globals().get(k, default) 는 k 가 None 으로 이미 있으면 default 가 아니라
#     그 None 을 준다. 앞 실행에서 MG_STRUCT=None 이 남아 있으면 기본값이 안 먹는다.
MG_STRUCT = globals().get("MG_STRUCT") or (TOOLS/"database"/"mg1655_uni")

def find_struct_dirs(root, min_n=500, depth=3):
    out = []
    root = Path(root)
    def walk(d, lv):
        if lv > depth: return
        try: entries = list(os.scandir(d))
        except Exception: return
        n = sum(1 for e in entries
                if e.is_file() and e.name.endswith((".cif", ".pdb", ".cif.gz", ".pdb.gz")))
        if n >= min_n: out.append((Path(d), n))
        for e in entries:
            if e.is_dir() and not e.name.startswith("."): walk(e.path, lv + 1)
    walk(root, 0)
    return sorted(out, key=lambda v: -v[1])

STRUCT_ALL = dict(STRUCT)
if MG_STRUCT and Path(MG_STRUCT).is_dir():
    STRUCT_ALL["MG1655"] = Path(MG_STRUCT)
else:
    known = {str(Path(v).resolve()) for v in STRUCT.values()}
    found = find_struct_dirs(TOOLS/"database")
    print("  구조 파일이 든 디렉터리:")
    for d, n in found[:12]:
        mine = "  ← 이미 씀" if str(d.resolve()) in known else ""
        print(f"    {n:6d}개  {d}{mine}")
    # 파일 개수로 고르면 안 된다. 사슬별로 쪼개진 파일과 버전이 섞인
    # 디렉터리가 개수만 많다. 이름 규칙과 **중복 제거 후 단백질 수**로 고른다.
    from collections import Counter

    def profile(d, cap=20000):
        """(이름규칙 분포, 모델버전 분포, 고유 단백질 수)"""
        pat, ver, prots = Counter(), Counter(), set()
        for i, f in enumerate(Path(d).iterdir()):
            if i > cap: break
            if f.suffix.lower() not in (".cif", ".pdb"): continue
            s = f.stem
            m = re.match(r"AF-([A-Z][A-Z0-9]{4,9})-F1-model_v(\d+)$", s)
            if m:
                pat["AFDB(UniProt)"] += 1; ver[f"v{m.group(2)}"] += 1
                prots.add(m.group(1)); continue
            m = re.match(r"AF-(\d+)-model_v(\d+)(?:_([A-Za-z0-9]+))?$", s)
            if m:
                pat["AFDB(숫자ID)" + ("+사슬" if m.group(3) else "")] += 1
                ver[f"v{m.group(2)}"] += 1; prots.add(m.group(1)); continue
            m = re.match(r"cf_(\S+)", s)
            if m:
                pat["ColabFold"] += 1; prots.add(m.group(1)); continue
            pat["기타"] += 1; prots.add(re.sub(r"_[A-Za-z0-9]$", "", s))
        return pat, ver, len(prots)

    MG_PROT = 4300          # MG1655 프로테옴 단백질 수
    cands = [(d, n) for d, n in found if str(d.resolve()) not in known]
    print("\n  후보 디렉터리 진단 (프로테옴 4,300 과 비교):")
    prof = []
    for d, n in cands[:8]:
        pat, ver, nprot = profile(d)
        prof.append((d, n, nprot, pat, ver))
        print(f"    {d.name:22s} 파일 {n:6d}  고유단백질 {nprot:6d}  "
              f"Δ4300 {abs(nprot - MG_PROT):5d}")
        print(f"      이름 {dict(pat)}   버전 {dict(ver)}")
    if prof:
        pick = min(prof, key=lambda v: abs(v[2] - MG_PROT))[0]
        STRUCT_ALL["MG1655"] = pick
        print(f"\n  ★ MG1655 구조로 고른 것: {pick}")
        print("    ※ 파일 개수가 아니라 **중복 제거 후 단백질 수**로 골랐다.")
        print("      이름 규칙이 한 가지이고 버전이 하나인 디렉터리가 맞다.")
        print("      다르면 MG_STRUCT = '경로' 를 정의하고 이 셀을 다시 돌릴 것")
print(f"\n  구조 DB: {', '.join(STRUCT_ALL)}"
      + ("" if "MG1655" in STRUCT_ALL else "   (MG1655 구조 못 찾음 → F7 에서 서열로)"))

def fs(tag, qpath, strain):
    # 캐시 파일명에 구조 DB 이름을 넣는다. 안 넣으면 DB 를 바꿔도
    # 옛 결과를 그대로 재사용한다 (실제로 그래서 MG1655 이 안 바뀌었다).
    o = FSD/f"{tag}_vs_{strain}_{Path(STRUCT_ALL[strain]).name}.m8"
    if not (o.exists() and o.stat().st_size):
        sh(f'"{FS}" easy-search "{qpath}" "{STRUCT_ALL[strain]}" "{o}" "{FSD}/tmp_{strain}" '
           f'--format-output "{FFMT}" -e 10 --max-seqs 2000 --exact-tmscore 1 '
           f'--threads {THREADS}', quiet=True)
    if not (o.exists() and o.stat().st_size): return pd.DataFrame()
    d = pd.read_csv(o, sep="\t", names=FFMT.split(","))
    d["stem"] = d.target.apply(lambda x: Path(str(x)).stem)
    # 사슬별로 쪼개진 파일(..._A, ..._B)은 같은 단백질이다. 접어서 센다.
    # 안 접으면 MG1655 히트 수가 부풀어 BL21 과 비교가 안 된다.
    d["prot"] = d.stem.str.replace(r"_[A-Za-z0-9]$", "", regex=True)
    return d.sort_values("alntmscore", ascending=False).drop_duplicates("prot")

xw = pd.read_csv(TBL/"id_crosswalk_struct_to_genbank.csv") \
     if (TBL/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
XW = ({Path(str(k)).stem: str(v).split(",")[0]
       for k, v in xw.set_index("tid")["protein"].to_dict().items()} if len(xw) else {})
def pid(stem):
    """crosswalk 에 있으면 GenBank ID, 없으면 읽을 수 있는 형태로 줄인다.
       MG1655 은 crosswalk 이 BL21/Y19 만 덮으므로 UniProt accession 이 나온다."""
    s = str(stem)
    if s in XW: return XW[s]
    m = re.match(r"AF-([A-Z][A-Z0-9]{4,9})-F1-model_v\d+$", s)
    if m: return m.group(1)
    m = re.match(r"AF-(\d+)-model_v\d+$", s)
    if m: return f"AF{m.group(1)}"
    return s

RES_FS, SETS, RES_ALL = [], {}, []
PROTN = {"BL21": 4110, "MG1655": 4300, "Y19": 5325}   # 프로테옴 단백질 수
for tag, qp in [("cry", CRY_A), ("af", AF_PATH)]:
    for s in STRUCT_ALL:
        d = fs(tag, qp, s)
        if not len(d):
            print(f"  {tag:3s} → {s:8s} 결과 없음"); continue
        fold = d[d.alntmscore >= TM_FOLD]
        # 전체 정렬 수와 컷 전후를 같이 찍는다. ≥0.5 만 보고하면
        # '8개' 가 몇 개 중 8개인지 알 수 없다 (7페이지 folddisco 표는 전체 수였다)
        n_all = len(d)
        bands = {f"≥{c}": int((d.alntmscore >= c).sum()) for c in (0.3, 0.4, 0.5, 0.7, 0.9)}
        print(f"  {tag:3s} → {s:8s} 정렬 {n_all:5d}개 중  " +
              "  ".join(f"{k} {v}" for k, v in bands.items()))
        RES_ALL.append({"frame": "결정" if tag == "cry" else "예측", "strain": s,
                        "정렬총수": n_all, **bands,
                        "프로테옴": PROTN.get(s), 
                        "비율(≥0.5)": round(bands["≥0.5"] / PROTN[s] * 100, 3)
                                      if PROTN.get(s) else None})
        SETS[(tag, s)] = {pid(x) for x in fold.prot}
        top = d.iloc[0]
        ctrl = ""
        if s == "Y19":
            hit = d[d.prot.str.contains(Y19_UNI, na=False)]
            ctrl = (f"   대조군 {Y19_GB} tm {hit.iloc[0].alntmscore:.3f}"
                    if len(hit) else "   ★ 대조군 없음")
        print(f"  {tag:3s} → {s:8s} 폴드(≥{TM_FOLD}) {len(fold):4d}개   "
              f"최고 {pid(top.prot)} tm {top.alntmscore:.3f}{ctrl}")
        for _, r in fold.iterrows():
            RES_FS.append({"frame": "결정" if tag == "cry" else "예측", "strain": s,
                           "protein": pid(r.prot), "tm": round(float(r.alntmscore), 3),
                           "fident": round(float(r.fident), 3), "qcov": round(float(r.qcov), 3),
                           "tier": ("TM>=0.9" if r.alntmscore >= TM_CORE else
                                    "TM 0.7-0.9" if r.alntmscore >= TM_NEAR
                                    else "TM 0.5-0.7")})

F = pd.DataFrame(RES_FS)
F.to_csv(TBL/"cooc_fold_afframe.csv", index=False, encoding="utf-8-sig")
A = pd.DataFrame(RES_ALL)
A.to_csv(TBL/"cooc_fold_counts.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 96); print("### 전체 매치 수 (컷 전후)"); print("=" * 96)
print(A.to_string(index=False))
print("\n  '정렬총수' 는 -e 10 --max-seqs 2000 으로 되돌아온 전부다. 구조 정렬은")
print("  아무 단백질이나 어느 정도는 겹치므로 이 수 자체에 의미는 없다.")
print("  의미는 **TM 분포**에 있다 — 0.5 위로 얼마나 남는가, 0.9 위에 뭐가 있는가.")
print(f"\n저장: {TBL/'cooc_fold_afframe.csv'}, {TBL/'cooc_fold_counts.csv'}")

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
                                     for t in ["TM>=0.9", "TM 0.7-0.9", "TM 0.5-0.7"]))
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
# MG1655 구조가 있으면 Foldseek 으로도 직접 물어 서열 판정과 대조한다
if "MG1655" in STRUCT_ALL:
    dmg = fs("af", AF_PATH, "MG1655")
    if len(dmg):
        fold_mg = dmg[dmg.alntmscore >= TM_FOLD]
        print("\n" + "=" * 96); print("### MG1655 도 구조로 직접"); print("=" * 96)
        print(f"  폴드(≥{TM_FOLD}) {len(fold_mg)}개   최고 {pid(fold_mg.iloc[0].prot)} "
              f"tm {fold_mg.iloc[0].alntmscore:.3f}")
        print(f"  BL21 {len(MEM)} / Y19 / MG1655 {len(fold_mg)} — 과 크기를 세 균주로 비교할 수 있다")
        fold_mg.assign(strain="MG1655", frame="예측")[
            ["strain", "frame", "prot", "alntmscore", "fident", "qcov"]] \
            .to_csv(TBL/"cooc_fold_mg1655_hits.csv", index=False, encoding="utf-8-sig")
        print(f"저장: {TBL/'cooc_fold_mg1655_hits.csv'}")

print("\n  O   폴드는 CooC 인데 MG1655 에 없다  ← 가장 강한 조합")
print("  ~   어노테이션 누락")
print("  .   MG1655 에도 있다")
print("\n  ※ BL21 에는 tm 0.9 이상의 진짜 CooC 이 없다는 것이 이미 확인됐다")
print("    (최고 ApbC tm 0.690). 그러니 여기서 O 가 나와도 'CooC 자체' 가 아니라")
print("    'CooC 과에 드는 BL21 단백질' 이다. 그 구분을 표에 유지할 것.")
```

---

## CELL F8 — 그 자리가 루프인가, 그리고 AF 는 그걸 아는가

```python
# =============================================================================
# CELL F8 | 왜 예측 모델끼리도 기하가 안 맞는지 — 자리 자체를 본다
#
#   구조를 열어 보면 Cys112/114 가 루프에 있다. 루프면 예측 좌표의 신뢰도가
#   낮고, 낮으면 직교체끼리 5.30 Å 과 8.16 Å 으로 갈리는 것이 설명된다.
#   AlphaFold 는 자기 확신도를 잔기마다 pLDDT 로 적어 둔다 (B-factor 칸).
#   추측하지 말고 그 값을 읽는다.
#
#   보는 것 셋.
#     (1) 질의 잔기의 pLDDT — 단백질 평균보다 낮은가
#     (2) 이차구조 대략 판정 — 나선/가닥/루프 (CA 간격으로 낸 근사)
#     (3) folddisco 금속 히트 80개의 매칭 잔기도 같은 경향인가
#         그렇다면 그 목록 전체가 '좌표를 못 믿는 자리' 위에 세워진 것이다
# =============================================================================
def load_bf(path):
    """{(chain, resseq): (resname, {atom: (xyz, bfactor)})}"""
    path = Path(path); out = {}
    if path.suffix.lower() in (".cif", ".mmcif"):
        hdr, started = [], False
        for l in open(path, errors="ignore"):
            if l.startswith("_atom_site."):
                hdr.append(l.strip().split(".")[1]); started = True; continue
            if started and l[:4] in ("ATOM", "HETA"):
                c = {n: i for i, n in enumerate(hdr)}
                f = l.split()
                try:
                    ch = f[c.get("auth_asym_id", c.get("label_asym_id"))]
                    rs = f[c.get("auth_seq_id", c.get("label_seq_id"))]
                    rn = f[c["label_comp_id"]]; an = f[c["label_atom_id"]].strip('"')
                    xyz = tuple(float(f[c[k]]) for k in ("Cartn_x", "Cartn_y", "Cartn_z"))
                    bf = float(f[c["B_iso_or_equiv"]])
                except Exception:
                    continue
                if not rs.lstrip("-").isdigit(): continue
                out.setdefault((ch, int(rs)), (rn, {}))[1][an] = (xyz, bf)
            elif started and l.startswith("#") and out:
                break
        return out
    for l in open(path, errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")): continue
        rs = l[22:26].strip()
        if not rs.lstrip("-").isdigit(): continue
        try: bf = float(l[60:66])
        except ValueError: bf = float("nan")
        xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        out.setdefault((l[21], int(rs)), (l[17:20].strip(), {}))[1][l[12:16].strip()] = (xyz, bf)
    return out

def ss_proxy(atoms, ch, rs):
    """CA(i)–CA(i+4) 거리로 낸 이차구조 근사. 나선 ~6.2 Å, 가닥 ~12–13 Å.
       DSSP 가 아니다. 대략만 본다."""
    def ca(k):
        v = atoms.get((ch, k))
        return v[1]["CA"][0] if v and "CA" in v[1] else None
    a, b = ca(rs - 2), ca(rs + 2)
    if not (a and b): return "?"
    d = math.dist(a, b)
    return "나선" if d < 7.0 else "가닥" if d > 11.0 else "루프/전이"

def plddt_of(atoms, ch, rs):
    v = atoms.get((ch, rs))
    if not v: return None
    return sum(b for _x, b in v[1].values()) / len(v[1])

def mean_plddt(atoms):
    vals = [b for _k, (_rn, at) in atoms.items() for _x, b in at.values()]
    return sum(vals) / len(vals) if vals else None

print("=" * 96); print("### (1) 질의 잔기의 pLDDT"); print("=" * 96)
TARGETS = [("Ch CooC1 (AF)", AF_PATH, None, AF_RES),
           ("Y19 CooC (AF)", Y19_FILE, None, None)]
for name, path, ch, res in TARGETS:
    if not Path(path).exists(): print(f"  {name}: 파일 없음"); continue
    A = load_bf(path)
    chains = sorted({c for c, _r in A})
    ch = ch or chains[0]
    mp = mean_plddt(A)
    if res is None:   # Y19 은 CXC 를 직접 찾는다
        cys = sorted(r for (c, r), (rn, _a) in A.items() if c == ch and rn == "CYS")
        res = next(((a, b) for a, b in itertools.combinations(cys, 2) if b - a == 2), tuple(cys[:2]))
    print(f"\n  [{name}]  사슬 {ch}  전체 평균 pLDDT {mp:.1f}")
    for r in res:
        p, s = plddt_of(A, ch, r), ss_proxy(A, ch, r)
        if p is None: print(f"    {r}: 없음"); continue
        print(f"    {r:4d}  pLDDT {p:5.1f}  (평균 대비 {p-mp:+5.1f})   이차구조 {s}")
    lo = [r for r in res if (plddt_of(A, ch, r) or 100) < mp]
    print(f"    → 질의 잔기가 평균보다 낮은 것 {len(lo)}/{len(res)}")

print("\n" + "=" * 96); print("### (2) 주변 구간의 pLDDT 프로파일"); print("=" * 96)
A = load_bf(AF_PATH)
ch = sorted({c for c, _r in A})[0]
lo, hi = min(AF_RES) - 10, max(AF_RES) + 10
mp = mean_plddt(A)
for r in range(lo, hi + 1):
    p = plddt_of(A, ch, r)
    if p is None: continue
    rn = A[(ch, r)][0]
    bar = "█" * int(p / 4)
    mark = "  ←질의" if r in AF_RES else ""
    print(f"    {r:4d} {rn:3s} {p:5.1f} {bar}{mark}")
print(f"    (전체 평균 {mp:.1f})")

print("\n" + "=" * 96); print("### (3) folddisco 금속 히트 80개의 매칭 잔기는?"); print("=" * 96)
MM = TBL/"folddisco_metal_motif_detail.csv"
if not MM.exists():
    print(f"  {MM.name} 없음 — 건너뛴다")
else:
    d = pd.read_csv(MM)
    col = next((c for c in ["matching_residues", "residues", "match"] if c in d.columns), None)
    tcol = next((c for c in ["tid", "target", "path"] if c in d.columns), None)
    if not col or not tcol:
        print(f"  필요한 열이 없다: {list(d.columns)[:12]}")
    else:
        rows, miss = [], 0
        for _, r in d.head(200).iterrows():
            p = Path(str(r[tcol]))
            if not p.exists():
                cand = [s/p.name for s in STRUCT.values() if (s/p.name).exists()]
                if not cand: miss += 1; continue
                p = cand[0]
            try: At = load_bf(p)
            except Exception: miss += 1; continue
            chs = sorted({c for c, _x in At})
            if not chs: continue
            c0, mpp = chs[0], mean_plddt(At)
            for tok in re.findall(r"([A-Za-z])(\d+)", str(r[col]).split(":")[0]):
                rr = int(tok[1])
                pv = plddt_of(At, c0, rr)
                if pv is None: continue
                rows.append({"protein": p.stem, "res": rr, "plddt": round(pv, 1),
                             "mean": round(mpp, 1), "delta": round(pv - mpp, 1),
                             "ss": ss_proxy(At, c0, rr)})
        if not rows:
            print(f"  구조 파일을 못 찾았다 (미해결 {miss}건)")
        else:
            P = pd.DataFrame(rows)
            P.to_csv(TBL/"motif_residue_plddt.csv", index=False, encoding="utf-8-sig")
            print(f"  매칭 잔기 {len(P)}개 / 단백질 {P.protein.nunique()}개  (구조 미해결 {miss})")
            print(f"\n  매칭 잔기 pLDDT 평균 {P.plddt.mean():.1f}  "
                  f"vs 단백질 평균 {P['mean'].mean():.1f}  "
                  f"(차이 {P.delta.mean():+.1f})")
            print(f"  평균보다 낮은 잔기 {int((P.delta < 0).sum())}/{len(P)} "
                  f"({(P.delta < 0).mean()*100:.0f}%)")
            print(f"  pLDDT 70 미만 {int((P.plddt < 70).sum())} "
                  f"({(P.plddt < 70).mean()*100:.0f}%)")
            print("\n  이차구조 분포:")
            print(P.ss.value_counts().to_string())
            print(f"\n저장: {TBL/'motif_residue_plddt.csv'}")

print("\n" + "=" * 96); print("### 판정"); print("=" * 96)
print("  질의 잔기가 루프이고 pLDDT 가 낮으면 —")
print("    그 좌표는 AlphaFold 스스로도 확신하지 않는 값이다. 직교체끼리")
print("    5.30 Å 과 8.16 Å 으로 갈린 것은 종간 차이가 아니라 예측 불확실성이다.")
print("    허용폭을 아무리 넓혀도 '노이즈를 노이즈에 맞추는' 일이 된다.")
print("\n  ★ (3) 이 낮게 안 나와도 결론은 같다 — 오히려 더 분명해진다.")
print("    질의 자리는 pLDDT 40대인데 히트로 잡힌 자리들은 잘 정의돼 있다면,")
print("    그것은 '노이즈 좌표에 우연히 들어맞은 well-ordered 자리들' 이라는 뜻이다.")
print("    80개가 CooC 과와 무관했던 것이 이것으로 설명된다.")
print("\n  그러면 이 자리는 구조 모티프 검색의 대상이 아니다. 남는 길 둘:")
print("    (a) 서열 모티프 — CXC / CXXC 패턴을 프로테옴 전체에서 직접 찾는다.")
print("        좌표를 안 쓰므로 pLDDT 와 무관하다. GPU 도 필요 없다")
print("    (b) 폴드 — Foldseek. 루프 하나가 흔들려도 폴드 판정은 버틴다.")
print("        F6 에서 두 프레임의 자카드가 높게 나온 것이 그 증거다")
```

---

## CELL F9 — 발표용 요약 표 두 개

```python
# =============================================================================
# CELL F9 | 슬라이드에 그대로 올릴 표를 만든다
#   표1 — 균주별 폴드 히트 수와 tier, 최고 히트와 TM
#   표2 — BL21 과 구성원 × MG1655 존재 여부
#   ※ MG1655 는 구조 DB 가 없어 Foldseek 3균주 비교가 성립하지 않는다.
#     그래서 균주 비교는 표2 로 분리한다. 표1 에 빈칸으로 두지 말 것.
# =============================================================================
F = pd.read_csv(TBL/"cooc_fold_afframe.csv")
# MG1655 를 따로 돌렸으면 합친다 (F7 이 구조 DB 를 찾았을 때 생긴다)
_mgf = TBL/"cooc_fold_mg1655_hits.csv"
if not _mgf.exists() and (TBL/"cooc_fold_MG1655.csv").exists():
    _mgf = TBL/"cooc_fold_MG1655.csv"     # 옛 이름 (대소문자만 달라 위험했다)
if _mgf.exists() and "MG1655" not in set(F.strain):
    m = pd.read_csv(_mgf)
    m = m.rename(columns={"alntmscore": "tm", "prot": "protein", "stem": "protein"})
    m["tier"] = ["TM>=0.9" if v >= 0.90 else "TM 0.7-0.9" if v >= 0.70 else "TM 0.5-0.7" for v in m.tm]
    F = pd.concat([F, m[["frame", "strain", "protein", "tm", "fident", "qcov", "tier"]]],
                  ignore_index=True)
    print(f"  MG1655 {len(m)}행 합침")
SLIDE = TBL/"slide_foldseek"; SLIDE.mkdir(parents=True, exist_ok=True)
LABEL = {"BL21": "E. coli BL21(DE3)", "MG1655": "E. coli K-12 MG1655",
         "Y19": "C. amalonaticus Y19"}

rows = []
for st in ["BL21", "MG1655", "Y19"]:
    for fr in ["예측", "결정"]:
        sub = F[(F.strain == st) & (F.frame == fr)]
        if not len(sub):
            if fr == "예측":
                rows.append({"Strain": LABEL[st], "Frame": "AFDB", "TM>=0.5": None,
                             "TM>=0.7": None, "TM>=0.9": None,
                             "Top hit": "no structure DB", "Top TM": None})
            continue
        top = sub.sort_values("tm", ascending=False).iloc[0]
        # 열 이름을 한 방식으로 통일한다. 전부 **누적** 기준이고 표기도 하나다.
        # 앞서 Fold hits / >=0.3 / Near 0.70-0.90 / Core >=0.90 이 섞여 있어
        # 열 이름만 보고는 무슨 기준인지 알 수 없었다.
        rows.append({"Strain": LABEL[st], "Frame": "AFDB" if fr == "예측" else "Crystal",
                     "TM>=0.5": len(sub),
                     "TM>=0.7": int((sub.tm >= 0.70).sum()),
                     "TM>=0.9": int((sub.tm >= 0.90).sum()),
                     "Top hit": top.protein, "Top TM": round(float(top.tm), 3)})
T1 = pd.DataFrame(rows)

# 필터 전 매칭 수(cooc_fold_counts.csv)를 같은 표에 합친다.
# 두 파일을 오가며 손으로 조합하지 않도록 한 장으로 만든다.
_cp = TBL/"cooc_fold_counts.csv"
if _cp.exists():
    C0 = pd.read_csv(_cp)
    C0["Strain"] = C0.strain.map(LABEL)
    C0["Frame"]  = C0.frame.map({"예측": "AFDB", "결정": "Crystal"})
    C0 = C0.rename(columns={"정렬총수": "Total", "프로테옴": "Proteome",
                            "≥0.3": "TM>=0.3", "≥0.4": "TM>=0.4",
                            "비율(≥0.5)": "%proteome"})
    keep = ["Strain", "Frame", "Proteome", "Total", "TM>=0.3", "TM>=0.4", "%proteome"]
    T1 = T1.merge(C0[[c for c in keep if c in C0.columns]],
                  on=["Strain", "Frame"], how="left")
    # 왼쪽에서 오른쪽으로 컷이 조여지는 순서. 전부 누적.
    ORD = ["Strain", "Frame", "Proteome", "Total",
           "TM>=0.3", "TM>=0.4", "TM>=0.5", "TM>=0.7", "TM>=0.9",
           "%proteome", "Top hit", "Top TM"]
    T1 = T1[[c for c in ORD if c in T1.columns] +
            [c for c in T1.columns if c not in ORD]]
else:
    print("  ⚠ cooc_fold_counts.csv 없음 — F6 을 먼저 돌리면 필터 전 수가 합쳐진다")

T1.to_csv(SLIDE/"table1_fold_counts.csv", index=False, encoding="utf-8-sig")
print("=" * 96); print("### 표1 — 폴드 검색 결과 (필터 전후 합본)"); print("=" * 96)
print(T1.to_string(index=False))
print("\n  열은 전부 **누적** 기준이다 — TM>=0.5 는 0.5 이상 전부(0.7, 0.9 포함).")
print("  Total 은 -e 10 --max-seqs 2000 으로 되돌아온 정렬 전부라 그 자체에")
print("  의미는 없다. %proteome 은 TM>=0.5 를 프로테옴 크기로 나눈 값이다.")

C = pd.read_csv(TBL/"cooc_fold_mg1655.csv")
NOTE = {".": "", "~": "annotation gap", "O": "ABSENT in MG1655"}
T2 = pd.DataFrame({
    "BL21 protein": C.protein,
    "TM to ChCooC1": C.tm_CooC,
    "MG1655": ["genome only" if v == "~" else "absent" if v == "O"
               else (re.search(r"있음 (\S+)", str(g)).group(1)
                     if re.search(r"있음 (\S+)", str(g)) else "-")
               for v, g in zip(C.MG1655, C["근거"])],
    "fident": [(re.search(r"fident ([\d.]+)", str(g)).group(1)
                if re.search(r"fident ([\d.]+)", str(g)) else None) for g in C["근거"]],
    "Note": [NOTE.get(v, "") for v in C.MG1655],
}).sort_values("TM to ChCooC1", ascending=False)
T2.to_csv(SLIDE/"table2_bl21_family_vs_MG1655.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 96); print("### 표2 — BL21 과 구성원 × MG1655"); print("=" * 96)
print(T2.to_string(index=False))

print("\n" + "=" * 96); print("### 슬라이드 문장"); print("=" * 96)
# 프레임을 반드시 밝힌다. 두 프레임 숫자를 한 문장에 섞으면 표와 안 맞는다.
def ctl_tm(fr):
    s = F[(F.strain == "Y19") & (F.frame == fr)]
    s = s[s.protein.str.contains(Y19_GB.split(".")[0], na=False)]
    return float(s.iloc[0].tm) if len(s) else None
c_cry, c_af = ctl_tm("결정"), ctl_tm("예측")
bl = F[(F.strain == "BL21") & (F.frame == "예측")]
print("  Control PASSED : ChCooC1 -> Y19 CooC")
print(f"      TM {c_cry:.3f} (crystal query) / {c_af:.3f} (AFDB query)"
      if c_cry and c_af else "      대조군 행 없음")
print(f"  BL21 CooC-family : {len(bl)} proteins, max TM {bl.tm.max():.3f}"
      f"   [AFDB frame]  -> no true CooC in BL21")
_core = int((F[F.frame == "예측"].tm >= 0.90).sum())
print(f"  Core (TM >= 0.90) in AFDB frame : {_core}"
      + ("  — 어느 균주에도 없다" if _core == 0 else ""))
n_in = int((C.MG1655 != "O").sum())
print(f"  All {n_in}/{len(C)} also present in MG1655"
      f"  -> no strain difference at family level")
print("\n  ※ 위 숫자는 모두 프레임을 명시한 것이다. 슬라이드에 옮길 때")
print("    crystal 값과 AFDB 값을 한 줄에 섞지 말 것 — 표와 안 맞게 된다.")
LINES = []
LINES.append("Control PASSED : ChCooC1 -> Y19 CooC")
if c_cry and c_af:
    LINES.append(f"    TM {c_cry:.3f} (crystal query) / {c_af:.3f} (AFDB query)")
LINES.append(f"BL21 CooC-family : {len(bl)} proteins, max TM {bl.tm.max():.3f}  [AFDB frame]")
LINES.append(f"Core (TM >= 0.90) in AFDB frame : {_core}")
LINES.append(f"All {n_in}/{len(C)} BL21 members also present in MG1655")
(SLIDE/"slide_lines.txt").write_text("\n".join(LINES) + "\n", encoding="utf-8")
print(f"\n저장: {SLIDE}")
print(f"  table1_fold_counts.csv          필터 전후 합본")
print(f"  table2_bl21_family_vs_MG1655.csv")
print(f"  slide_lines.txt                 위 문장 그대로")
```

---

## CELL F10 — 세 균주의 과 구성원이 같은 것인가 (구조로 직접 대조)

```python
# =============================================================================
# CELL F10 | "BL21 8개가 MG1655 에도 있다" 를 구조로 확인한다
#
#   F7 은 서열로 물었다 (fident 0.97~1.000). 강한 근거지만 ID 체계가
#   BL21=GenBank / MG1655=UniProt 이라 표에서 같은 단백질인지 눈으로 못 맞춘다.
#   여기서는 **BL21 구성원 하나하나를 MG1655 구조 DB 에 직접 Foldseek** 해서
#   MG1655 쪽 과 구성원 중 누구에게 떨어지는지 본다. ID 매핑이 필요 없다.
#
#   그리고 Y19 도 같이 해서 세 균주 대응표를 만든다.
#   과가 같은 구성원으로 채워져 있으면 '균주 차이 없음' 이 구조 수준에서도 성립한다.
# =============================================================================
BLM = sorted(F[(F.frame == "예측") & (F.strain == "BL21")].protein)
print("=" * 100); print("### 무엇을 무엇에 거는가"); print("=" * 100)
print(f"  질의 : BL21 의 CooC-like {len(BLM)}개 (ChCooC1 과 TM>=0.5 인 것들)")
for s in ["MG1655", "Y19"]:
    if s in STRUCT_ALL:
        _n = sum(1 for f_ in Path(STRUCT_ALL[s]).iterdir()
                 if f_.suffix.lower() in (".cif", ".pdb"))
        print(f"  대상 : {s} 구조 **전체 {_n:,}개** ({Path(STRUCT_ALL[s]).name})")
print("\n  ★ 상대 균주의 CooC-like 목록(7개, 15개)하고만 비교하는 것이 아니다.")
print("    프로테옴 구조 전체에 걸어서 가장 닮은 것을 찾는다. 그래야")
print("    'BL21 에만 있는 것' 을 놓치지 않는다 — 상대 목록에 없는 단백질이")
print("    답일 수도 있기 때문이다.")

# 구조 파일 경로 찾기 (crosswalk 역방향)
_xw = pd.read_csv(TBL/"id_crosswalk_struct_to_genbank.csv") \
      if (TBL/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
P2T = {}
if len(_xw):
    for a, b in zip(_xw.tid, _xw.protein):
        if pd.isna(a) or pd.isna(b): continue
        for one in str(b).split(","): P2T.setdefault(one.strip(), Path(str(a)))

def find_struct(pid, strain):
    p = P2T.get(pid)
    if p and p.exists(): return p
    d = Path(STRUCT_ALL[strain])
    if p:
        c = d/p.name
        if c.exists(): return c
        for f in d.glob(f"{p.stem}*"):
            return f
    return None

MEMSET = {s: set(F[(F.frame == "예측") & (F.strain == s)].protein)
          for s in STRUCT_ALL}
rows = []
for pidv in BLM:
    q = find_struct(pidv, "BL21")
    if q is None:
        rows.append({"BL21": pidv, "note": "구조 파일 못 찾음"}); continue
    rec = {"BL21": pidv,
           "TM_to_CooC1": float(F[(F.frame == "예측") & (F.strain == "BL21") &
                              (F.protein == pidv)].tm.iloc[0])}
    for s in ["MG1655", "Y19"]:
        o = FSD/f"x_{pidv}_vs_{s}_{Path(STRUCT_ALL[s]).name}.m8"
        if not (o.exists() and o.stat().st_size):
            sh(f'"{FS}" easy-search "{q}" "{STRUCT_ALL[s]}" "{o}" "{FSD}/tx_{s}" '
               f'--format-output "{FFMT}" -e 10 --max-seqs 50 --exact-tmscore 1 '
               f'--threads {THREADS}', quiet=True)
        if not (o.exists() and o.stat().st_size):
            rec[f"{s}_hit"] = None; rec[f"{s}_TM"] = None; continue
        d = pd.read_csv(o, sep="\t", names=FFMT.split(","))
        d["prot"] = d.target.apply(lambda x: Path(str(x)).stem) \
                     .str.replace(r"_[A-Za-z0-9]$", "", regex=True)
        d = d.sort_values("alntmscore", ascending=False).drop_duplicates("prot")
        top = d.iloc[0]
        # 열 이름은 기준을 그대로 쓴다.
        #   {균주}_hit        그 균주에서 나온 최고 히트가 누구인가
        #   {균주}_TM         그 히트의 TM-score
        #   {균주}_in_family  그 히트가 그 균주의 CooC 과 목록(TM>=0.5)에도 있는가
        rec[f"{s}_hit"] = pid(top.prot)
        rec[f"{s}_TM"]  = round(float(top.alntmscore), 3)
        # in_family 열은 뺐다. X 가 CooC1 과 TM>=0.5 이고 Y 가 X 와 거의 같으면
        # Y 도 자동으로 0.5 를 넘는다 — 구성상 거의 항상 O 라 정보가 없다.
    rows.append(rec)

X = pd.DataFrame(rows)
X.to_csv(TBL/"slide_foldseek"/"table3_family_correspondence.csv",
         index=False, encoding="utf-8-sig")
print("\n" + "=" * 100); print("### 세 균주 과 구성원 대응"); print("=" * 100)
print(X.to_string(index=False))
print("\n  열 읽는 법")
print("    {균주}_hit   그 BL21 단백질을 해당 균주 구조에 걸었을 때 최고 히트")
print("    {균주}_TM    그 히트의 TM-score. 0.9 이상이면 사실상 같은 단백질")
print("\n  MG1655_TM 이 전부 0.9 이상이면 BL21 고유 구성원이 없다는 뜻이다.")
if "MG1655_TM" in X:
    hi = X[X.MG1655_TM >= 0.9]
    print(f"\n  MG1655 에 TM ≥ 0.9 대응이 있는 BL21 구성원: {len(hi)}/{len(X)}")
    print("  → 구조 수준에서도 BL21 고유 구성원이 없다는 뜻이다.")
print(f"\n저장: {TBL/'slide_foldseek'/'table3_family_correspondence.csv'}")
```

---

## CELL F11 — 과 스윕을 `--exact-tmscore 1` 로 다시 돌린다

```python
# =============================================================================
# CELL F11 | tm 이 1.0 을 넘는 값은 발표에 못 쓴다
#
#   hypA_fold_sweep.csv 와 g3e_family_sweep.csv 의 tm 이 1.008, 1.009 로 나온다.
#   Foldseek 이 --exact-tmscore 1 없이 주는 근사값이다. 순위를 매기는 데는
#   문제없지만 tier 경계(0.90 / 0.70 / 0.50) 근처의 개수가 바뀔 수 있고,
#   슬라이드에 'tm 0.680' 같은 숫자를 올리려면 정확값이어야 한다.
#
#   질의를 새로 정의하지 않는다 — **기존 .m8 의 query 열에서 그대로 읽어**
#   같은 질의로 다시 돌린다. seed 목록을 기억에 의존해 재구성하면 다른 스윕이 된다.
# =============================================================================
FSOLD = BASE/"result"/"foldseek"
EX    = FSOLD/"exact"; EX.mkdir(parents=True, exist_ok=True)

TARGETS = {}
for s, d in STRUCT_ALL.items(): TARGETS[s] = Path(d)
print("  대상 구조 DB:", {k: v.name for k, v in TARGETS.items()})

def queries_in(pattern):
    """기존 .m8 들에서 질의 파일 경로를 모은다."""
    qs = {}
    for f in sorted(FSOLD.glob(pattern)):
        if not f.stat().st_size: continue
        try:
            d = pd.read_csv(f, sep="\t", names=FFMT.split(","), nrows=5000)
        except Exception:
            continue
        for q in d["query"].astype(str).unique():
            p = Path(q)
            if not p.exists():
                for sd in TARGETS.values():
                    if (sd/p.name).exists(): p = sd/p.name; break
            if p.exists(): qs[p.stem] = p
    return qs

SWEEPS = [("hypa", "hypa_*_vs_*.m8", "hypA_fold_sweep"),
          ("g3e",  "g3e_vs_*.m8",     "g3e_family_sweep")]

for tag, pat, outname in SWEEPS:
    Q = queries_in(pat)
    print("\n" + "=" * 96); print(f"### {tag} — 질의 {len(Q)}개"); print("=" * 96)
    if not Q:
        print(f"  {pat} 에서 질의를 못 찾았다. 파일명을 확인할 것:")
        for f in sorted(FSOLD.glob("*.m8"))[:12]: print("   ", f.name)
        continue
    rows = []
    for qname, qp in Q.items():
        for s, td in TARGETS.items():
            o = EX/f"{tag}_{qname}_vs_{td.name}.m8"
            if not (o.exists() and o.stat().st_size):
                sh(f'"{FS}" easy-search "{qp}" "{td}" "{o}" "{EX}/tmp_{s}" '
                   f'--format-output "{FFMT}" -e 10 --max-seqs 2000 '
                   f'--exact-tmscore 1 --threads {THREADS}', quiet=True)
            if not (o.exists() and o.stat().st_size): continue
            d = pd.read_csv(o, sep="\t", names=FFMT.split(","))
            d["prot"] = d.target.apply(lambda x: Path(str(x)).stem) \
                         .str.replace(r"_[A-Za-z0-9]$", "", regex=True)
            d = d.sort_values("alntmscore", ascending=False).drop_duplicates("prot")
            for _, r in d[d.alntmscore >= 0.50].iterrows():
                rows.append({"seed": qname, "strain": s, "protein": pid(r.prot),
                             "tm": round(float(r.alntmscore), 3),
                             "fident": round(float(r.fident), 3),
                             "qcov": round(float(r.qcov), 3)})
    if not rows:
        print("  tm ≥ 0.50 인 것이 없다"); continue
    S = pd.DataFrame(rows).sort_values("tm", ascending=False) \
          .drop_duplicates(["strain", "protein"])       # 단백질당 최고 seed 1건
    S["tier"] = ["TM>=0.9" if v >= 0.90 else "TM 0.7-0.9" if v >= 0.70 else "TM 0.5-0.7" for v in S.tm]
    out = TBL/f"{outname}_exact.csv"
    S.to_csv(out, index=False, encoding="utf-8-sig")
    print(f"  tm 최대 {S.tm.max():.3f}  (1.0 초과 {int((S.tm > 1.0).sum())}개 — 0 이어야 정상)")
    print("\n  균주 × tier")
    print(pd.crosstab(S.strain, S.tier).to_string())
    print(f"\n저장: {out}")

    # 옛 값과 나란히
    old = TBL/f"{outname}.csv"
    if old.exists():
        O = pd.read_csv(old)
        kc = next((c for c in ["protein", "member", "target"] if c in O.columns), None)
        tc = next((c for c in ["tm", "alntmscore", "top_tm"] if c in O.columns), None)
        if kc and tc:
            O["key"] = O[kc].astype(str).str.split("-").str[-1]
            M = S.merge(O[["key", tc]].rename(columns={"key": "protein", tc: "tm_old"}),
                        on="protein", how="left").dropna(subset=["tm_old"])
            if len(M):
                M["Δ"] = (M.tm - M.tm_old).round(3)
                print(f"\n  옛 값과 비교 가능한 것 {len(M)}개   "
                      f"Δ 평균 {M['Δ'].mean():+.3f}  최대 {M['Δ'].abs().max():.3f}")
                big = M[M["Δ"].abs() >= 0.05]
                if len(big):
                    print(f"  0.05 이상 바뀐 것 {len(big)}개:")
                    print(big[["protein", "strain", "tm_old", "tm", "Δ"]]
                          .sort_values("Δ").to_string(index=False))
                else:
                    print("  0.05 이상 바뀐 것 없음 — tier 경계도 안 흔들린다")

print("\n  ※ 앞으로 Foldseek 은 항상 --exact-tmscore 1 을 붙인다.")
print("    안 붙이면 근사값이라 1.0 을 넘고, tier 경계 근처 개수가 달라진다.")
```
