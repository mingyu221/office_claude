# 서열 모티프 → 폴드 → Boltz — 새 파이프라인

## 왜 순서를 바꾸나

지금까지의 입구는 **"MG1655 에 없는 것"** 이었다. 그런데 이 연구의 관찰은
*lysate 차이* 이고, 보고서 자체의 결론은 *"차이는 유전자 구성이 아니라 발현량"* 이다.
둘이 모순이고, 그 결과 **identical + high-similarity 3,746개가 평가를 못 받았다.**

그리고 구조 모티프(Folddisco) 축은 쓸 수 없다는 것이 확인됐다 —
질의 잔기 Cys112/114 가 **루프**에 있어 예측 좌표를 믿을 수 없다.
같은 자리를 AF 가 Ch CooC1 에서 5.30 Å, Y19 CooC 에서 8.16 Å 으로 예측한다.

**서열 모티프는 좌표를 쓰지 않는다.** pLDDT 와 무관하고, GPU 도 필요 없고,
세 균주 프로테옴 전체(13,735개)를 몇 분에 훑는다.

## 새 순서

```
1. 서열 모티프  세 균주 프로테옴 전체에서 Cys 모티프 스캔        Q2   (CPU, 수 분)
2. 균주 비교    BL21 모티프 보유 → MG1655 에 모티프째 보존됐나   Q3   (mmseqs)
3. 폴드 교차    Foldseek 과(科) + pLDDT 로 좌표 신뢰도           Q4
4. 우선순위     점수 매겨 상위 N 선정 + Boltz 입력 생성           Q5
5. Boltz        후보 단독 + Ni. ChCODH2 를 빼서 5~10배 빠르다     Q6   (GPU)
6. 채점         기존과 같은 2.6 Å / Cys 개수 등급                 Q6
```

**입구가 "Ni 을 잡을 수 있는 자리를 가졌나" 로 바뀐다.** 균주 비교는 입구가 아니라
순위 요소가 된다. 그러면 3,746개도 심사를 받는다.

### Boltz 에서 ChCODH2 를 빼는 근거

기존 채점 코드가 `if comp == "NI" or ch != "B": continue` 로 **후보 단백질 사슬만**
셌다. ChCODH2(사슬 A)는 등급에 한 번도 안 쓰였다. 런타임은 ChCODH2 ~630 잔기가
거의 다 먹는다. 빼면 같은 숫자가 나오면서 훨씬 빨라진다.

---

## CELL Q1 — 설정

```python
# =============================================================================
# CELL Q1 | 독립 실행용 설정
# =============================================================================
import re, math, json, subprocess, itertools, tempfile
from pathlib import Path
from collections import Counter
import pandas as pd

BASE  = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS = Path("/mnt/af2results/mingyu")
TBL   = BASE/"result"/"table"
SRCH  = BASE/"result"/"search"
OUT   = BASE/"result"/"seqmotif"; OUT.mkdir(parents=True, exist_ok=True)
THREADS = 8

FS  = TOOLS/"foldseek"/"bin"/"foldseek"
FAA = {"BL21":   TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"bl21_db_match_qjz.faa",
       "MG1655": TOOLS/"database"/"protein_list"/"mg1655_protein.faa",
       "Y19":    TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"y19_db_match.faa"}
GEN_MG = BASE/"input"/"external"/"MG1655_U00096.3.fna"
STRUCT = {"BL21": TOOLS/"database"/"bacteriaDB"/"structures_UP000503272",
          "Y19":  TOOLS/"database"/"bacteriaDB"/"structures_UP000034085"}

def sh(cmd, quiet=False):
    r = subprocess.run(cmd, shell=True, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if not quiet: print(r.stdout.rstrip())
    return r.stdout

def read_faa(p):
    seqs, hdrs, k, buf = {}, {}, None, []
    for l in open(p, errors="ignore"):
        if l.startswith(">"):
            if k: seqs[k] = "".join(buf)
            t = l[1:].rstrip(); k = t.split()[0]; hdrs[k] = t[len(k):].strip(); buf = []
        else: buf.append(l.strip())
    if k: seqs[k] = "".join(buf)
    return seqs, hdrs

print("=" * 96); print("### 입력"); print("=" * 96)
SEQ, HDR = {}, {}
for st, p in FAA.items():
    if p.exists():
        SEQ[st], HDR[st] = read_faa(p)
        print(f"  O  {st:8s} {len(SEQ[st]):6d}개  {p.name}")
    else:
        print(f"  X  {st:8s} {p}")
print(f"  {'O' if GEN_MG.exists() else 'X'}  MG1655 게놈  {GEN_MG.name}")
```

---

## CELL Q2 — 세 균주 프로테옴 전체에서 Cys 모티프를 찾는다

```python
# =============================================================================
# CELL Q2 | 좌표를 쓰지 않는 검색
#
#   무엇을 찾나 — Ni(II) 같은 연한 금속을 잡는 자리는 티올레이트를 쓴다.
#   서열에서는 Cys 가 짧은 간격으로 붙어 나오는 패턴으로 나타난다.
#
#     CXC    C.C     CooC1 의 Cys112-x-Cys114 가 이것
#     CXXC   C..C    HypA · 티오레독신 · DsbA 형
#     CX3C   C...C
#     CX4C   C....C
#     CX5C   C.....C
#
#   그리고 Cys4 자리는 보통 이런 쌍이 **둘** 모여 만들어진다 (아연 손가락,
#   zinc ribbon, HypA 의 N말단 Ni 자리 전부 그렇다). 그래서 쌍의 개수와
#   창(window) 안의 Cys 밀도를 따로 센다.
#
#   ※ 이 검색은 '금속 자리다' 를 판정하지 않는다. **후보를 좁히는 그물**이다.
#     실제 판정은 Q6 의 Boltz 배위 등급이 한다. 그물이 넓은 것은 의도된 것이다.
# =============================================================================
PATTERNS = {"CXC":  r"(?=(C.C))",
            "CXXC": r"(?=(C..C))",
            "CX3C": r"(?=(C...C))",
            "CX4C": r"(?=(C....C))",
            "CX5C": r"(?=(C.....C))"}
WIN      = 40      # Cys4 자리가 들어갈 만한 서열 창
MIN_CYS4 = 4       # 창 안에 이만큼 있으면 Cys4 후보

def scan(s):
    pos = {n: [m.start() + 1 for m in re.finditer(p, s)] for n, p in PATTERNS.items()}
    cys = [i + 1 for i, a in enumerate(s) if a == "C"]
    best_win, best_at = 0, None
    for i, c in enumerate(cys):
        j = i
        while j + 1 < len(cys) and cys[j + 1] - c <= WIN: j += 1
        if j - i + 1 > best_win: best_win, best_at = j - i + 1, c
    npair = sum(len(v) for v in pos.values())
    return {"n_cys": len(cys), "n_pair": npair,
            **{f"n_{n}": len(v) for n, v in pos.items()},
            "cys4_win": best_win, "cys4_at": best_at,
            "pairs": ";".join(f"{n}@{p}" for n, v in pos.items() for p in v)[:120],
            "cys_pos": ",".join(map(str, cys))[:160]}

rows = []
for st in SEQ:
    for k, s in SEQ[st].items():
        r = scan(s)
        if r["n_pair"] == 0 and r["cys4_win"] < MIN_CYS4: continue
        r.update({"strain": st, "protein": k, "len": len(s),
                  "desc": HDR[st].get(k, "")[:70]})
        rows.append(r)
M = pd.DataFrame(rows)
M["cys4_cand"] = (M.cys4_win >= MIN_CYS4).astype(int)
M.to_csv(TBL/"seqmotif_all_strains.csv", index=False, encoding="utf-8-sig")

print("=" * 96); print("### 균주별 모티프 보유"); print("=" * 96)
tot = {st: len(SEQ[st]) for st in SEQ}
hdr = f"  {'':8s} {'프로테옴':>8s} {'모티프보유':>10s} {'비율':>7s} " \
      f"{'CXC':>6s} {'CXXC':>6s} {'Cys4후보':>8s} {'비율':>7s}"
print(hdr)
for st in SEQ:
    sub = M[M.strain == st]
    print(f"  {st:8s} {tot[st]:8d} {len(sub):10d} {len(sub)/tot[st]*100:6.1f}% "
          f"{int((sub.n_CXC>0).sum()):6d} {int((sub.n_CXXC>0).sum()):6d} "
          f"{int(sub.cys4_cand.sum()):8d} {sub.cys4_cand.sum()/tot[st]*100:6.2f}%")

print("\n  ※ 세 균주 비율이 비슷하게 나오는 것이 정상이다. 이 단계는 균주를")
print("    가르는 것이 아니라 심사 대상을 좁히는 것이다. 가르는 일은 Q3 가 한다.")
print(f"\n저장: {TBL/'seqmotif_all_strains.csv'}")

print("\n" + "=" * 96); print("### Cys4 후보 중 짧은 것 (샤페론 프로파일)"); print("=" * 96)
small = M[(M.strain == "BL21") & (M.cys4_cand == 1) & (M.len <= 250)]
print(f"  BL21 {len(small)}개")
print(small.sort_values("len").head(20)[
    ["protein", "len", "n_cys", "n_CXC", "n_CXXC", "cys4_win", "desc"]].to_string(index=False))
```

---

## CELL Q3 — MG1655 에 **모티프째** 보존돼 있는가

```python
# =============================================================================
# CELL Q3 | '있다/없다' 가 아니라 '그 Cys 가 상대에도 Cys 인가'
#
#   앞서 두 번 속았다.
#     - 구조 DB 버전 차이로 'BL21 에만 검출' 14개가 나왔는데 진짜 부재는 0개
#     - 프로테옴 부재 7개 중 5개가 게놈에는 온전 (어노테이션 누락)
#   그래서 여기서는 세 단계로 묻는다.
#     (a) MG1655 프로테옴에 상동체가 있는가        mmseqs, qcov >= 0.70
#     (b) 있다면 그 상동체에 **모티프 Cys 가 보존**됐는가   정렬을 잔기 단위로 따라감
#     (c) 상동체가 없다면 게놈에는 있는가          --search-type 2, fident >= 0.80
# =============================================================================
MFMT = "query,target,fident,alnlen,evalue,bits,qlen,tlen,qcov,tcov,qstart,qend,qaln,taln"
CAND = M[(M.strain == "BL21") & ((M.cys4_cand == 1) | (M.n_pair > 0))].copy()
print(f"  BL21 모티프 보유 {len(CAND)}개에게 묻는다")

QF = OUT/"bl21_motif.faa"
QF.write_text("".join(f">{k}\n{SEQ['BL21'][k]}\n" for k in CAND.protein))

def msearch(q, tgt, out, extra=""):
    if not (out.exists() and out.stat().st_size):
        with tempfile.TemporaryDirectory() as td:
            sh(f'mmseqs easy-search "{q}" "{tgt}" "{out}" "{td}" '
               f'--format-output "{MFMT}" -s 7.5 -e 1e-3 --max-seqs 5 '
               f'--threads {THREADS} {extra}', quiet=True)
    if not (out.exists() and out.stat().st_size): return pd.DataFrame()
    return pd.read_csv(out, sep="\t", names=MFMT.split(","))

d = msearch(QF, FAA["MG1655"], SRCH/"seqmotif_bl21_vs_MG1655.m8")
BEST = {}
if len(d):
    d = d[d.qcov >= 0.70].sort_values("bits", ascending=False).drop_duplicates("query")
    BEST = {str(r["query"]): r for _, r in d.iterrows()}
print(f"  (a) 프로테옴 상동체(qcov>=0.70) {len(BEST)} / {len(CAND)}")

def cys_kept(row, positions):
    """질의의 Cys 위치들이 정렬 상대에서도 C 인가. (보존수, 검사수)"""
    qa, ta = str(row["qaln"]), str(row["taln"])
    qi = int(row["qstart"]); ok = tot_ = 0
    want = set(positions)
    for a, b in zip(qa, ta):
        if a != "-":
            if qi in want:
                tot_ += 1
                if b == a: ok += 1
            qi += 1
    return ok, tot_

rows = []
for _, c in CAND.iterrows():
    pid = c.protein
    pos = [int(x) for x in str(c.cys_pos).split(",") if x]
    r = BEST.get(pid)
    if r is None:
        rows.append({"protein": pid, "mg_hit": None, "fident": None,
                     "cys_kept": None, "cys_total": len(pos), "verdict": "프로테옴 없음"})
        continue
    ok, tt = cys_kept(r, pos)
    v = ("동일" if float(r.fident) >= 0.999 else
         "모티프 보존" if tt and ok == tt else
         "모티프 달라짐" if tt else "정렬 구간 밖")
    rows.append({"protein": pid, "mg_hit": str(r.target), "fident": round(float(r.fident), 3),
                 "cys_kept": ok, "cys_total": tt, "verdict": v})
V = pd.DataFrame(rows)
print("\n  (b) 모티프 Cys 보존 판정")
print(V.verdict.value_counts().to_string())

NOHIT = list(V[V.verdict == "프로테옴 없음"].protein)
GH = {}
if NOHIT and GEN_MG.exists():
    gf = OUT/"bl21_motif_nohit.faa"
    gf.write_text("".join(f">{k}\n{SEQ['BL21'][k]}\n" for k in NOHIT))
    g = msearch(gf, GEN_MG, SRCH/"seqmotif_nohit_vs_MG1655genome.m8", "--search-type 2")
    if len(g):
        # --search-type 2 의 alnlen 은 염기 단위다. 커버리지는 qstart/qend 로 낸다
        g["qcov2"] = (g.qend - g.qstart + 1) / g.qlen
        g = g[(g.fident >= 0.80) & (g.qcov2 >= 0.70)] \
              .sort_values("bits", ascending=False).drop_duplicates("query")
        GH = {str(r["query"]): float(r.fident) for _, r in g.iterrows()}
print(f"\n  (c) 프로테옴 미검출 {len(NOHIT)} 중 게놈에 있는 것 {len(GH)}")

V["mg_status"] = [("게놈에만 있음" if p in GH else v) if v == "프로테옴 없음" else v
                  for p, v in zip(V.protein, V.verdict)]
V["bl21_only"] = (V.mg_status == "프로테옴 없음").astype(int)
V.to_csv(TBL/"seqmotif_strain_verdict.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 96); print("### 최종 균주 판정"); print("=" * 96)
print(V.mg_status.value_counts().to_string())
print(f"\n  ★ MG1655 에 프로테옴·게놈 둘 다 없는 것: {int(V.bl21_only.sum())}개")
print(V[V.bl21_only == 1].head(20).to_string(index=False))
print(f"\n저장: {TBL/'seqmotif_strain_verdict.csv'}")
```

---

## CELL Q3b — "MG1655 에 없다" 가 프로파지인지 가른다

```python
# =============================================================================
# CELL Q3b | 유전체 위치를 본다
#
#   BL21(DE3) 는 T7 발현을 위해 **λDE3 프로파지를 삽입해 만든 균주**다.
#   그 구간의 유전자는 MG1655 에 없는 것이 당연하다 — 균주 제작의 산물이지
#   생물학적 신호가 아니다. 다른 프로파지·이동성 구간도 마찬가지다.
#
#   구분하는 법: 진짜 단독 유전자는 유전체에 흩어져 있고, 삽입 구간은
#   **번호가 연속으로 뭉쳐 있다.** QJZ 번호와 HO396 로커스 태그가 유전체
#   순서를 따르므로 그대로 쓸 수 있다.
#
#   ※ 이 셀은 후보를 **버리지 않는다.** 표시만 한다. lysate 차이라는 관찰
#     앞에서 프로파지를 배제할 근거는 없지만, 사전확률이 낮다는 것은
#     표에 남아 있어야 한다.
# =============================================================================
GAP        = 20     # QJZ 번호 단위. 프로파지 한 덩어리가 이 정도로 벌어져 있다
MIN_BLOCK  = 3      # 이만큼 모이면 삽입 구간으로 의심
MOBILE_KW  = ["phage", "terminase", "tail", "capsid", "portal", "integrase",
              "transposase", "recombinase", "antitermination", "excisionase",
              "holin", "lysozyme", "baseplate", "host specificity", "Bet",
              "IS[0-9]", "insertion sequence", "prophage", "tRNA-", "repressor"]

def locus_num(pid, desc):
    """QJZ 번호만 쓴다.
       HO396 로커스 태그는 어노테이션이 있는 단백질 설명에만 붙어 있어서,
       HO396 우선으로 두면 어떤 행은 HO396(예: 3665), 어떤 행은 QJZ(예: 11292)가
       되어 **두 좌표계가 한 축에 섞인다.** 실제로 그렇게 나왔다.
       QJZ 는 모든 행에 있고 HO396 과 순서가 일치하므로 이쪽으로 통일한다."""
    m = re.search(r"QJZ(\d+)", str(pid))
    return int(m.group(1)) if m else None

ONLY = V[V.bl21_only == 1].merge(CAND[["protein", "desc", "len"]], on="protein", how="left")
# 열 이름을 'loc' 으로 두면 안 된다 — Series.loc 은 pandas 의 인덱서라
# r.loc 이 열 값이 아니라 _LocIndexer 를 준다. 이름을 바꿔 충돌을 없앤다.
ONLY["locnum"] = [locus_num(p, d) for p, d in zip(ONLY.protein, ONLY.desc)]
ONLY = ONLY.dropna(subset=["locnum"]).sort_values("locnum")

blocks, cur = [], []
for _, r in ONLY.iterrows():
    if cur and r["locnum"] - cur[-1]["locnum"] <= GAP: cur.append(r)
    else:
        if cur: blocks.append(cur)
        cur = [r]
if cur: blocks.append(cur)

print("=" * 100); print("### MG1655 에 없는 것들의 유전체 배치"); print("=" * 100)
MOBILE = set()
for b in blocks:
    ids = [x.protein for x in b]
    lo, hi = int(b[0]["locnum"]), int(b[-1]["locnum"])
    lo, hi = f"QJZ{lo}", f"QJZ{hi}"
    kw = sum(1 for x in b if re.search("|".join(MOBILE_KW), str(x.desc), re.I))
    tag = ""
    if len(b) >= MIN_BLOCK:
        tag = "  ★ 삽입 구간 의심"; MOBILE |= set(ids)
    elif kw:
        tag = "  (이동성 어노테이션)"; MOBILE |= set(ids)
    print(f"\n  {lo}–{hi}  {len(b)}개  이동성 어노테이션 {kw}/{len(b)}{tag}")
    for x in b:
        mark = "·" if re.search("|".join(MOBILE_KW), str(x.desc), re.I) else " "
        print(f"    {mark} {x.protein:12s} {int(x.len):5d}aa  {str(x.desc)[:66]}")

print("\n" + "=" * 100); print("### 판정"); print("=" * 100)
print(f"  MG1655 에 없는 것 {len(ONLY)}개")
print(f"    삽입 구간/이동성   {len(MOBILE)}개  ← 균주 제작·프로파지의 산물일 가능성")
print(f"    흩어진 단독 유전자 {len(ONLY) - len(MOBILE)}개  ← 이쪽이 관심 대상")
solo = ONLY[~ONLY.protein.isin(MOBILE)]
if len(solo):
    print("\n  단독:")
    print(solo[["protein", "len", "desc"]].to_string(index=False))

V["mobile"] = V.protein.isin(MOBILE).astype(int)
V.to_csv(TBL/"seqmotif_strain_verdict.csv", index=False, encoding="utf-8-sig")
print(f"\n  V 에 mobile 열 추가 후 재저장: {TBL/'seqmotif_strain_verdict.csv'}")
print("\n  ※ BL21(DE3) 는 λDE3 를 삽입해 만든 균주다. 그 구간이 MG1655 에 없는 것은")
print("    설계상 그런 것이지 발견이 아니다. 발표에서는 두 줄을 나눠 쓸 것 —")
print("    '균주 특이 33개 중 프로파지 N개를 빼면 M개' 가 정직한 문장이다.")
```

---

## CELL Q4 — 폴드와 좌표 신뢰도를 붙인다

```python
# =============================================================================
# CELL Q4 | 서열 모티프만으로는 그물이 너무 넓다. 두 가지를 붙여 좁힌다.
#   (1) 폴드 — Foldseek 으로 이미 뽑아 둔 과(科) 목록에 드는가
#       (CooC/G3E, HypA). 없으면 이 셀이 직접 돌린다
#   (2) pLDDT — 모티프 Cys 구간의 예측 신뢰도. 낮으면 Boltz 에 넣어도
#       좌표를 못 믿는다. **Folddisco 가 실패한 이유가 이것이었다**
# =============================================================================
FOLD = {}
for f, tag in [("hypA_fold_sweep.csv", "HypA"), ("g3e_family_sweep.csv", "G3E"),
               ("cooc_fold_afframe.csv", "CooC")]:
    p = TBL/f
    if not p.exists(): print(f"  - {f} 없음"); continue
    d = pd.read_csv(p)
    kc = next((c for c in ["protein", "member", "target"] if c in d.columns), d.columns[0])
    tc = next((c for c in ["tm", "alntmscore", "top_tm", "tm_CooC"] if c in d.columns), None)
    for _, r in d.iterrows():
        k = str(r[kc]).split("-")[-1]
        v = float(r[tc]) if tc and pd.notna(r.get(tc)) else None
        cur = FOLD.get(k)
        if cur is None or (v or 0) > (cur[1] or 0): FOLD[k] = (tag, v)
    print(f"  + {f}: {len(d)}행 → 과 정보 누적 {len(FOLD)}")
    if tc and pd.to_numeric(d[tc], errors="coerce").max() > 1.0:
        print(f"    ⚠ {f} 의 {tc} 가 1.0 을 넘는다 — --exact-tmscore 1 없이 돌린 근사값이다.")
        print("      순위에는 써도 되지만 발표 숫자로는 재실행한 값을 쓸 것.")

def _build_struct_index():
    """단백질 ID → 구조 파일 경로. 두 단계로 만든다.
       (1) 구조 디렉터리를 훑어 {파일 stem: 경로}
       (2) crosswalk 으로 {단백질 ID: stem}
       앞서 tid 를 그대로 경로로 믿었다가 926개 전부 실패했다 —
       tid 가 절대경로가 아니거나 서버 경로가 바뀌면 그대로 0 이 된다."""
    idx = {}
    for s in STRUCT.values():
        if not Path(s).is_dir(): continue
        for f in Path(s).iterdir():
            if f.suffix.lower() in (".cif", ".pdb"):
                idx[f.stem] = f
                # AF-XXXX-F1-model_v6 → XXXX 로도 걸 수 있게
                m = re.match(r"AF-([A-Z0-9]+)-F1", f.stem)
                if m: idx.setdefault(m.group(1), f)
    p2s = {}
    xw = TBL/"id_crosswalk_struct_to_genbank.csv"
    if xw.exists():
        d = pd.read_csv(xw)
        tc = next((c for c in ["tid", "structure", "path", "file"] if c in d.columns), None)
        pc = next((c for c in ["protein", "genbank", "gb"] if c in d.columns), None)
        if tc and pc:
            for a, b in zip(d[tc], d[pc]):
                if pd.isna(a) or pd.isna(b): continue
                for one in str(b).split(","):
                    p2s.setdefault(one.strip(), Path(str(a)).stem)
    print(f"  구조 색인 {len(idx)}개 / crosswalk 매핑 {len(p2s)}개")
    return idx, p2s

SIDX, P2S = _build_struct_index()

def struct_of(pid):
    stem = P2S.get(pid)
    if stem and stem in SIDX: return SIDX[stem]
    if stem:                                   # 확장자·버전이 다를 때
        base = re.sub(r"-model_v\d+$", "", stem)
        for k, v in SIDX.items():
            if k.startswith(base): return v
    return SIDX.get(pid)

def plddt_window(path, lo, hi):
    """구간 [lo,hi] 잔기의 평균 pLDDT 와 단백질 전체 평균."""
    vals, allv = [], []
    if Path(path).suffix.lower() in (".cif", ".mmcif"):
        hdr, started = [], False
        for l in open(path, errors="ignore"):
            if l.startswith("_atom_site."): hdr.append(l.strip().split(".")[1]); started = True; continue
            if started and l[:4] in ("ATOM", "HETA"):
                c = {n: i for i, n in enumerate(hdr)}; f = l.split()
                try:
                    rs = int(f[c.get("auth_seq_id", c.get("label_seq_id"))])
                    b = float(f[c["B_iso_or_equiv"]])
                except Exception: continue
                allv.append(b)
                if lo <= rs <= hi: vals.append(b)
            elif started and l.startswith("#") and allv: break
    else:
        for l in open(path, errors="ignore"):
            if not l.startswith("ATOM"): continue
            try: rs, b = int(l[22:26]), float(l[60:66])
            except ValueError: continue
            allv.append(b)
            if lo <= rs <= hi: vals.append(b)
    return (sum(vals)/len(vals) if vals else None,
            sum(allv)/len(allv) if allv else None)

rows = []
for _, c in CAND.iterrows():
    pid = c.protein
    fam, tm = FOLD.get(pid, (None, None))
    pw = pm = None
    sp = struct_of(pid)
    if sp and c.cys4_at:
        pw, pm = plddt_window(sp, int(c.cys4_at), int(c.cys4_at) + WIN)
    rows.append({"protein": pid, "fold": fam, "fold_tm": tm,
                 "plddt_motif": round(pw, 1) if pw else None,
                 "plddt_mean": round(pm, 1) if pm else None,
                 "struct": sp.name if sp else None})
Q = pd.DataFrame(rows)
Q.to_csv(TBL/"seqmotif_fold_plddt.csv", index=False, encoding="utf-8-sig")
print(f"\n  폴드 붙은 것 {int(Q.fold.notna().sum())} / {len(Q)}")
print(f"  구조 찾은 것 {int(Q.struct.notna().sum())} / {len(Q)}")
g = Q.dropna(subset=["plddt_motif"])
if len(g):
    print(f"  모티프 구간 pLDDT 평균 {g.plddt_motif.mean():.1f} "
          f"vs 단백질 평균 {g.plddt_mean.mean():.1f}")
    print(f"  모티프 구간 70 미만 {int((g.plddt_motif < 70).sum())} "
          f"({(g.plddt_motif < 70).mean()*100:.0f}%)")
print(f"\n저장: {TBL/'seqmotif_fold_plddt.csv'}")
```

---

## CELL Q5 — 우선순위 + Boltz 입력 생성

```python
# =============================================================================
# CELL Q5 | 점수를 매겨 상위 N 만 GPU 로 보낸다
#
#   가중치의 근거 — 재현성이 확인된 축에 높게 준다.
#     +3  Cys4 후보 (창 안에 Cys 4개 이상)   Boltz grade A 가 곧 Cys4 다
#     +3  MG1655 프로테옴·게놈 둘 다 없음     유일하게 0 이 아니었던 균주 축
#     +2  CooC/G3E/HypA 폴드                  Foldseek 은 대조군을 통과했다
#     +2  CXC 또는 CXXC 보유
#     +1  250 잔기 이하                        샤페론은 대체로 작다
#     +1  모티프 구간 pLDDT >= 70              좌표를 믿을 수 있다
#     -2  모티프 구간 pLDDT < 50               믿을 수 없다
#
#   ※ 균주 축이 **입구가 아니라 가중치**다. 이것이 이 파이프라인의 핵심 변경점이다.
#     MG1655 도 가진 단백질이 Cys4 + 좋은 폴드면 여전히 심사를 받는다.
# =============================================================================
TOPN = 400

_vc = ["protein", "mg_status", "bl21_only"] + (["mobile"] if "mobile" in V.columns else [])
D = CAND.merge(V[_vc], on="protein", how="left").merge(Q, on="protein", how="left")
if "mobile" not in D.columns:
    D["mobile"] = 0
    print("  ⚠ mobile 열이 없다 — CELL Q3b 를 먼저 돌리면 프로파지가 구분된다")

def _pos(v):
    """NaN 안전. float('nan') 은 파이썬에서 참이라 `x or 0` 로 쓰면
       값이 없을 때도 점수가 붙는다."""
    return pd.notna(v) and float(v) > 0

def score(r):
    p, why = 0, []
    if r.cys4_cand == 1:
        p += 3; why.append("+3 Cys4후보")
    elif _pos(r.n_CXC) or _pos(r.n_CXXC):
        # Cys4 후보와 배타로 둔다. 둘 다 주면 같은 증거를 두 번 세는 것이고,
        # 그 결과 CooC 형(단량체당 CXC 하나)이 3점 뒤처져 상위에서 밀린다.
        p += 2; why.append("+2 CXC/CXXC")
    if r.bl21_only == 1:
        if r.get("mobile", 0) == 1:
            # 프로파지/삽입 구간은 MG1655 에 없는 것이 설계상 당연하다.
            # 버리지는 않되 '균주 특이' 가산은 주지 않는다.
            p += 0; why.append("(MG없음-프로파지)")
        else:
            p += 3; why.append("+3 MG없음")
    elif r.mg_status == "모티프 달라짐":        p += 1; why.append("+1 모티프차이")
    if pd.notna(r.fold):                       p += 2; why.append(f"+2 {r.fold}폴드")
    if pd.notna(r.len) and r.len <= 250:       p += 1; why.append("+1 소형")
    if pd.notna(r.plddt_motif):
        if r.plddt_motif >= 70:                p += 1; why.append("+1 pLDDT")
        elif r.plddt_motif < 50:               p -= 2; why.append("-2 pLDDT낮음")
    return p, " ".join(why)

D[["score", "why"]] = D.apply(lambda r: pd.Series(score(r)), axis=1)
D = D.sort_values(["score", "cys4_win", "len"], ascending=[False, False, True])
D.to_csv(TBL/"seqmotif_ranked.csv", index=False, encoding="utf-8-sig")

print("=" * 96); print("### 점수 분포"); print("=" * 96)
print(D.score.value_counts().sort_index(ascending=False).to_string())
print("\n" + "=" * 96); print(f"### 상위 25"); print("=" * 96)
print(D.head(25)[["protein", "len", "score", "cys4_win", "n_CXC", "n_CXXC",
                  "fold", "fold_tm", "plddt_motif", "mg_status", "mobile",
                  "desc"]].to_string(index=False))
print(f"\n저장: {TBL/'seqmotif_ranked.csv'}")

# ---- 눈금 맞추기: 기존 채점의 grade A/B 를 반드시 포함시킨다 ----
#   이번 예측은 ChCODH2 를 빼고 후보 단독으로 접는다. 조건이 바뀌었으므로
#   등급이 그대로 재현되는지 확인할 기준점이 필요하다. 기존 grade A/B 를
#   점수와 무관하게 넣어 두고, 같은 등급이 다시 나오는지 본다.
CAL = []
_old = TBL/"ni_site_grade.csv"
if _old.exists():
    _o = pd.read_csv(_old)
    _o["key"] = _o.prey.astype(str).str.split("-").str[-1]
    CAL = sorted({k for k, g in zip(_o.key, _o.grade)
                  if g in ("A", "B") and k in SEQ["BL21"]})
    print(f"\n  눈금용(기존 grade A/B) {len(CAL)}개를 강제 포함: {', '.join(CAL[:8])}"
          + (" …" if len(CAL) > 8 else ""))
else:
    print(f"\n  ⚠ {_old.name} 없음 — 눈금 비교 없이 간다")

# ---- Boltz 입력: 후보 단독 + Ni ----
#   ★ 폴더를 먼저 비운다. 안 비우면 이전 실행에서 상위 400 에 들었다가
#     이번에 빠진 단백질의 YAML 이 그대로 남아 같이 돌아간다.
#     (실제로 그래서 어느 버전이 만든 입력인지 알 수 없게 됐다)
import yaml as _y, json, shutil, datetime
BIN = BASE/"result"/"boltz"/"inputs_seqmotif"
if BIN.exists() and any(BIN.glob("*.yaml")):
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    old_dir = BIN.parent/f"inputs_seqmotif.old_{stamp}"
    shutil.move(str(BIN), str(old_dir))
    print(f"  기존 입력 {len(list(old_dir.glob('*.yaml')))}개를 {old_dir.name} 로 옮겼다")
BIN.mkdir(parents=True, exist_ok=True)
MAXLEN = 1200
sel_ids = list(dict.fromkeys(list(D.head(TOPN).protein) + CAL))
LEN = dict(zip(D.protein, D.len))
n, skip, wrote = 0, 0, []
for pidv in sel_ids:
    s = SEQ["BL21"].get(pidv)
    if not s: continue
    if len(s) > MAXLEN: skip += 1; continue
    # 사슬 A = 후보 단백질, 리간드 B = Ni. ChCODH2 는 넣지 않는다
    doc = {"version": 1, "sequences": [
        {"protein": {"id": "A", "sequence": s}},
        {"ligand":  {"id": "B", "ccd": "NI"}}]}
    (BIN/f"{pidv.replace('|','_')}.yaml").write_text(_y.safe_dump(doc, sort_keys=False))
    wrote.append(pidv); n += 1
TOTRES = sum(len(SEQ["BL21"][k]) for k in wrote)
# 어떤 설정으로 만든 입력인지 남긴다. 나중에 '이거 어느 버전이지' 를 없앤다
# ★ 입력 폴더 안에 두면 안 된다. boltz 는 폴더의 모든 파일을 입력으로 읽고
#   .json 을 만나면 "Unable to parse filetype .json" 으로 죽는다. 옆에 둔다.
MANIFEST = BIN.parent/f"{BIN.name}_manifest.json"
MANIFEST.write_text(json.dumps({
    "made_at": datetime.datetime.now().isoformat(timespec="seconds"),
    "TOPN": TOPN, "MAXLEN": MAXLEN,
    "n_written": n, "n_skipped_len": skip,
    "calibration": CAL,
    "mobile_aware": bool("mobile" in D.columns and D.mobile.notna().any()),
    "score_rule": "Cys4(+3) XOR CXC/CXXC(+2) / MG없음 비프로파지(+3) / 폴드(+2) / 소형(+1) / pLDDT(±)",
    "proteins": sorted(wrote),
}, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"\n  Boltz 입력 {n}개 생성 (상위 {TOPN} + 눈금 {len(CAL)}, 길이 초과 제외 {skip})")
print(f"  설정 기록: {MANIFEST}")
print(f"  → {BIN}")
print(f"  총 잔기 {TOTRES:,} — 기존 방식이면 여기에 630×{n} = {630*n:,} 이 더 붙는다")
print(f"\n  ※ --use_msa_server 로 {n}건의 MSA 를 원격에서 받는다. 속도 제한이 걸리면")
print("    나눠 돌리거나 로컬 MSA 로 바꿀 것. 예측 자체보다 여기서 막히는 일이 흔하다.")
```

---

## CELL Q5b — CooC 형 전용 트랙 (CXC 하나짜리를 곧바로 이량체로)

```python
# =============================================================================
# CELL Q5b | Cys 2개짜리를 버리지 않는다
#
#   CooC1 은 단량체당 CXC **하나**뿐이다 (Cys112-x-Cys114).
#   Cys4 는 이량체 계면에서 두 CXC 가 만나야 생긴다. 그러니
#     - 진짜 CooC 형은 Cys 2개로 충분하고,
#     - Q2 의 `+3 Cys4후보` 를 못 받아 상위 400 밖으로 밀리고,
#     - 단량체로 접으면 어차피 grade C 로 나온다.
#   세 번 불리하다. 그래서 별도 트랙으로 빼서 **처음부터 이량체로** 접는다.
#
#   선정 기준 — CooC1 의 서열 특징을 그대로 쓴다
#     · CXC 또는 CXXC 를 가짐
#     · Cys4 창이 4 미만 (혼자서는 Cys4 를 못 만든다) ← 여기가 핵심
#     · 2 × 길이가 토큰 상한 안
#   여기에 폴드·균주·길이를 가중치로 얹어 정렬한다.
# =============================================================================
import yaml as _y, json, shutil, datetime
DIN = BASE/"result"/"boltz"/"inputs_cooclike"
if DIN.exists() and any(DIN.glob("*.yaml")):
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    old_dir = DIN.parent/f"inputs_cooclike.old_{stamp}"
    shutil.move(str(DIN), str(old_dir))
    print(f"  기존 입력 {len(list(old_dir.glob('*.yaml')))}개를 {old_dir.name} 로 옮겼다")
DIN.mkdir(parents=True, exist_ok=True)
MAXTOK, TOPD = 1200, 200

CL = D[( (D.n_CXC.fillna(0) > 0) | (D.n_CXXC.fillna(0) > 0) ) &
       (D.cys4_win < 4)].copy()
print("=" * 96); print("### CooC 형 후보 (혼자서는 Cys4 를 못 만드는 것)"); print("=" * 96)
print(f"  CXC/CXXC 보유 + Cys4 창 < 4 : {len(CL)}개")
print(f"    그중 폴드 정보 있음 {int(CL.fold.notna().sum())}, "
      f"MG1655 에 없음 {int((CL.mg_status == '프로테옴 없음').sum())}")

def dscore(r):
    p, why = 0, []
    if pd.notna(r.fold):                        p += 3; why.append(f"+3 {r.fold}폴드")
    if r.mg_status == "프로테옴 없음":
        if r.get("mobile", 0) == 1:
            why.append("(MG없음-프로파지)")     # Q5 와 같은 규칙. 가산 없음
        else:
            p += 3; why.append("+3 MG없음")
    elif r.mg_status == "모티프 달라짐":          p += 1; why.append("+1 모티프차이")
    if (r.n_CXC or 0) > 0:                      p += 2; why.append("+2 CXC")   # CooC1 과 동일 배치
    elif (r.n_CXXC or 0) > 0:                   p += 1; why.append("+1 CXXC")
    if r.len <= 350:                            p += 1; why.append("+1 소형")  # CooC1 은 ~256
    if pd.notna(r.plddt_motif):
        if r.plddt_motif >= 70:                 p += 1; why.append("+1 pLDDT")
        elif r.plddt_motif < 50:                p -= 2; why.append("-2 pLDDT낮음")
    return p, " ".join(why)

CL[["dscore", "dwhy"]] = CL.apply(lambda r: pd.Series(dscore(r)), axis=1)
CL = CL.sort_values(["dscore", "len"], ascending=[False, True])
CL.to_csv(TBL/"cooclike_ranked.csv", index=False, encoding="utf-8-sig")
print("\n  상위 20")
print(CL.head(20)[["protein", "len", "dscore", "n_CXC", "n_CXXC", "fold",
                   "mg_status", "desc"]].to_string(index=False))

n, skip = 0, 0
for _, r in CL.head(TOPD).iterrows():
    s = SEQ["BL21"].get(r.protein)
    if not s: continue
    if 2 * len(s) > MAXTOK: skip += 1; continue
    doc = {"version": 1, "sequences": [
        {"protein": {"id": "A", "sequence": s}},
        {"protein": {"id": "B", "sequence": s}},      # 동형이량체
        {"ligand":  {"id": "L", "ccd": "NI"}}]}
    (DIN/f"{r.protein.replace('|','_')}_dimer.yaml").write_text(_y.safe_dump(doc, sort_keys=False))
    n += 1

# 양성대조군 — Ch CooC1. 이게 A-계면으로 안 나오면 아래 전부 못 믿는다
CRY = TOOLS/"workspace"/"seek_ni_insertase"/"input"/"3kji.pdb"
AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}
if CRY.exists():
    seq, seen = [], set()
    for l in open(CRY, errors="ignore"):
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "A":
            rs = l[22:26].strip()
            if rs not in seen: seen.add(rs); seq.append(AA3.get(l[17:20].strip(), "X"))
    cs = "".join(seq)
    (DIN/"CTRL-CooC1_dimer.yaml").write_text(_y.safe_dump(
        {"version": 1, "sequences": [
            {"protein": {"id": "A", "sequence": cs}},
            {"protein": {"id": "B", "sequence": cs}},
            {"ligand":  {"id": "L", "ccd": "NI"}}]}, sort_keys=False))
    print(f"\n  ★ 양성대조군 CTRL-CooC1_dimer.yaml ({len(cs)} 잔기 ×2)")
    # 음성대조군 — CXC 가 아예 없는 비슷한 길이의 단백질
    # Cys 쌍이 없을 뿐 아니라 Cys 자체가 거의 없어야 음성대조군이다.
    # 앞서 고른 것은 Cys4 창이 5 여서 Ni 을 잡아도 이상하지 않았다.
    neg = D[(D.n_CXC.fillna(0) == 0) & (D.n_CXXC.fillna(0) == 0) &
            (D.n_cys.fillna(9) <= 1) &
            (D.len.between(len(cs) - 60, len(cs) + 60))].head(1)
    if not len(neg):
        neg = D[(D.n_cys.fillna(9) <= 1)].sort_values(
            "len", key=lambda s: (s - len(cs)).abs()).head(1)
    if len(neg):
        ns = SEQ["BL21"][neg.iloc[0].protein]
        (DIN/f"NEG-{neg.iloc[0].protein}_dimer.yaml").write_text(_y.safe_dump(
            {"version": 1, "sequences": [
                {"protein": {"id": "A", "sequence": ns}},
                {"protein": {"id": "B", "sequence": ns}},
                {"ligand":  {"id": "L", "ccd": "NI"}}]}, sort_keys=False))
        print(f"  ★ 음성대조군 NEG-{neg.iloc[0].protein}_dimer.yaml "
              f"(CXC 없음, {len(ns)} 잔기 ×2)")
        print("    이게 A 로 나오면 Boltz 가 Ni 을 아무 데나 넣는다는 뜻이다.")

DMANIFEST = DIN.parent/f"{DIN.name}_manifest.json"   # 입력 폴더 밖에 둔다
DMANIFEST.write_text(json.dumps({
    "made_at": datetime.datetime.now().isoformat(timespec="seconds"),
    "TOPD": TOPD, "MAXTOK": MAXTOK, "n_written": n, "n_skipped_len": skip,
    "선정": "CXC/CXXC 보유 + cys4_win < 4 (혼자서는 Cys4 를 못 만드는 것)",
    "대조군": [f.name for f in DIN.glob("CTRL-*.yaml")] + [f.name for f in DIN.glob("NEG-*.yaml")],
}, ensure_ascii=False, indent=1), encoding="utf-8")
print(f"\n  이량체 입력 {n}개 (2×길이 > {MAXTOK} 제외 {skip}) + 대조군  → {DIN}")
print(f"  설정 기록: {DMANIFEST}")
print(f"""
cd "{BASE/'result'/'boltz'}"
CUDA_VISIBLE_DEVICES=0 boltz predict inputs_cooclike \\
  --out_dir out_cooclike \\
  --use_msa_server --max_msa_seqs 2048 \\
  --recycling_steps 3 --diffusion_samples 1 \\
  --output_format mmcif --num_workers 2
echo DONE_cooclike
""")
print("  채점은 Q8 의 DOUT 을 out_cooclike 로 바꿔 돌리면 된다 (A-계면 등급).")
```

---

## CELL Q6 — Boltz 실행 + 채점

```python
# =============================================================================
# CELL Q6 | 예측하고, 기존과 같은 기준으로 등급을 매긴다
#   채점 기준은 CELL 28a-8b 와 동일하게 유지한다. 비교 가능해야 하기 때문이다.
#     TIGHT 2.6 Å 안의 도너를 잔기 단위로 세고, Cys 개수로 A/B/C/D
#   다른 점 하나 — 기존 코드는 사슬 B 만 셌다(사슬 A 가 ChCODH2 였으므로).
#   여기서는 후보가 사슬 A 이므로 **리간드가 아닌 모든 사슬**을 센다.
# =============================================================================
BOUT = BASE/"result"/"boltz"/"out_seqmotif"
GPU  = 0
RUN  = f"""
cd "{BASE/'result'/'boltz'}"
CUDA_VISIBLE_DEVICES={GPU} boltz predict inputs_seqmotif \\
  --out_dir out_seqmotif \\
  --use_msa_server \\
  --max_msa_seqs 2048 \\
  --recycling_steps 3 \\
  --diffusion_samples 1 \\
  --output_format mmcif \\
  --num_workers 2
echo DONE_seqmotif
"""
print("=" * 96); print("### 실행 명령 (백그라운드로 돌릴 것)"); print("=" * 96)
print(RUN)
print(f"  결과 예정: {BOUT}")
print("\n  ※ 이 셀은 명령만 찍는다. 기존 노트북의 sh_bg 로 띄우거나 터미널에서 돌린다.")

# ---- 채점 (예측이 끝난 뒤 다시 실행) ----
import glob
TIGHT = 2.6
BACKBONE = ("N", "O", "C", "CA", "CB")

def cif_atoms(p):
    hdr, started, out = [], False, []
    for l in open(p, errors="ignore"):
        if l.startswith("_atom_site."): hdr.append(l.strip().split(".")[1]); started = True; continue
        if started and l[:4] in ("ATOM", "HETA"):
            c = {n: i for i, n in enumerate(hdr)}; f = l.split()
            try:
                out.append((f[c["label_comp_id"]],
                            f[c.get("auth_asym_id", c.get("label_asym_id"))],
                            f[c.get("auth_seq_id", c.get("label_seq_id"))],
                            f[c["label_atom_id"]].strip('"'),
                            float(f[c["Cartn_x"]]), float(f[c["Cartn_y"]]), float(f[c["Cartn_z"]])))
            except Exception: continue
        elif started and l.startswith("#") and out: break
    return out

cifs = sorted(glob.glob(str(BOUT/"**"/"*_model_0.cif"), recursive=True))
print(f"\n  예측 결과 {len(cifs)}개")
rows = []
for cif in cifs:
    A = cif_atoms(cif)
    ni = [a for a in A if a[0] == "NI"]
    if not ni: continue
    _, nch, _, _, x, y, z = ni[0]
    donors = []
    for comp, ch, seq, at, ax, ay, az in A:
        if comp == "NI" or ch == nch: continue        # 리간드 사슬은 제외
        d = math.dist((x, y, z), (ax, ay, az))
        if d <= TIGHT and (at not in BACKBONE or comp in ("ASP", "GLU", "HIS", "CYS")):
            donors.append((round(d, 2), comp, seq, at))
    donors = sorted({(c, s): (d, c, s, a) for d, c, s, a in donors}.values())
    cnt = Counter(c for _, c, _, _ in donors)
    ncys, nhis = cnt.get("CYS", 0), cnt.get("HIS", 0)
    if   ncys >= 4:          g, why = "A", "Cys4 — 연한 금속 전달·저장형 자리"
    elif ncys == 3:          g, why = "B", "Cys3 + 보조"
    elif ncys == 2 and nhis: g, why = "B", "Cys2His"
    elif ncys or nhis >= 3:  g, why = "C", "His 중심"
    else:                    g, why = "D", "산소 공여체 위주 — 표면 부착 의심"
    rows.append({"grade": g, "protein": Path(cif).stem.replace("_model_0", ""),
                 "n_donor": len(donors), "Cys": ncys, "His": nhis,
                 "min_dist": donors[0][0] if donors else None,
                 "donors": " ".join(f"{c}{s}" for _, c, s, _ in donors), "why": why})
if rows:
    G = pd.DataFrame(rows).sort_values(["grade", "Cys", "min_dist"],
                                       ascending=[True, False, True])
    G = G.merge(D[["protein", "score", "mg_status", "fold", "len", "desc"]],
                on="protein", how="left")
    G.to_csv(TBL/"seqmotif_ni_grade.csv", index=False, encoding="utf-8-sig")
    print("\n" + "=" * 96); print("### 등급 분포"); print("=" * 96)
    print(G.groupby("grade").size().to_string())
    print("\n" + "=" * 96); print("### grade A"); print("=" * 96)
    print(G[G.grade == "A"][["protein", "len", "Cys", "min_dist", "donors",
                             "mg_status", "fold", "desc"]].to_string(index=False))
    print(f"\n저장: {TBL/'seqmotif_ni_grade.csv'}")
    print("\n  기존 grade A 4개(HslO·그 Y19 직교체·QJZ12568.1·GspE)가 여기에도")
    print("  나오는지 확인할 것. 나오면 두 파이프라인이 서로를 검증한 것이고,")
    print("  새 A 가 나오면 그것이 기존 입구 필터가 놓친 것이다.")
else:
    print("  아직 결과가 없다. 위 명령을 먼저 돌릴 것.")
```

---

## CELL Q6b — 백그라운드 실행 + 모니터

```python
# =============================================================================
# CELL Q6b | boltz 를 백그라운드로 띄운다
#   노트북 셀에서 직접 돌리면 커널이 묶이고, 커널이 죽으면 예측도 죽는다.
#   setsid + nohup 으로 떼어 놓고 로그만 본다. 이 셀을 다시 돌리면
#   띄우지 않고 **상태만** 보여준다 (중복 실행 방지).
# =============================================================================
import os, subprocess, shlex, time, glob
from pathlib import Path

BOLTZ_ROOT = BASE/"result"/"boltz"
IN_DIR     = "inputs_seqmotif"          # 이량체는 "inputs_cooclike"
OUT_DIR    = "out_seqmotif"             #           "out_cooclike"
CONDA_ENV  = "boltz"                    # boltz 가 설치된 env 이름. 다르면 고칠 것
FORCE_GPU  = None                       # 숫자를 주면 그 GPU 로 고정한다.
                                        # 다른 작업(CELL 54 등)이 쓸 카드를 피할 때 쓴다.
                                        # 그쪽이 아직 예측을 시작 안 했으면 여유 메모리로는
                                        # 구분이 안 되므로 여기서 직접 정하는 편이 안전하다.
LOG        = BOLTZ_ROOT/f"{OUT_DIR}.log"
PIDF       = BOLTZ_ROOT/f"{OUT_DIR}.pid"

MIN_FREE = 12000    # MiB. 이만큼 안 비어 있으면 안 띄운다

# ---- 남의 boltz 가 이미 돌고 있는지부터 본다 ----
def other_boltz():
    """이 셀이 띄운 것이 아닌 boltz 프로세스 목록."""
    mine = None
    if PIDF.exists():
        try: mine = int(PIDF.read_text().strip())
        except Exception: pass
    out = subprocess.run("ps -eo pid,etime,args", shell=True, text=True,
                         capture_output=True).stdout.splitlines()
    hits = []
    for l in out:
        if "boltz" not in l or "ps -eo" in l: continue
        try: pid = int(l.split()[0])
        except Exception: continue
        if pid == mine or pid == os.getpid(): continue
        if "predict" in l or "envs/boltz" in l: hits.append((pid, l.strip()[:150]))
    return hits

# ---- GPU 선택: 여유 메모리 ----
def gpu_rows():
    try:
        o = subprocess.run("nvidia-smi --query-gpu=index,memory.free,memory.total,"
                           "utilization.gpu --format=csv,noheader,nounits",
                           shell=True, text=True, capture_output=True).stdout
        return [tuple(int(x) for x in l.split(",")) for l in o.strip().splitlines() if l.strip()]
    except Exception as e:
        print(f"  nvidia-smi 실패 ({e})"); return []

def pick_gpu():
    rows = gpu_rows()
    for i, f, tt, u in rows:
        mark = "" if f >= MIN_FREE else f"   ← {MIN_FREE} MiB 미만"
        print(f"  GPU {i}: 여유 {f:6d} / {tt} MiB   사용률 {u:3d}%{mark}")
    ok = [r for r in rows if r[1] >= MIN_FREE]
    if not ok:
        return None
    return max(ok, key=lambda r: r[1])[0]

def running():
    """기록된 PID 가 살아 있는 것만으로는 부족하다.
       setsid 로 띄우면 PIDF 에 들어가는 것은 **셸의 PID** 이고 실제 예측은
       그 자식이다. 자식이 죽어도 셸이 남아 있으면 '돌고 있다' 로 잘못 읽힌다
       (실제로 그래서 죽은 작업을 살아 있다고 보고했다).
       예측 프로세스가 실재하는지 따로 확인한다."""
    if not PIDF.exists(): return None
    try: pid = int(PIDF.read_text().strip())
    except Exception: return None
    if not Path(f"/proc/{pid}").exists():
        PIDF.unlink(missing_ok=True); return None
    alive = subprocess.run("pgrep -f '[b]oltz predict'", shell=True,
                           capture_output=True, text=True).stdout.split()
    if not alive:
        print(f"  기록된 PID {pid} 의 셸은 살아 있으나 boltz predict 가 없다 "
              f"— 죽은 작업이다. PID 파일을 지우고 다시 띄울 수 있게 한다.")
        print(f"    남은 셸도 정리할 것:  kill {pid}")
        PIDF.unlink(missing_ok=True)
        return None
    return pid

n_in  = len([f for f in (BOLTZ_ROOT/IN_DIR).glob("*.yaml")])
_man = BOLTZ_ROOT/f"{IN_DIR}_manifest.json"
if not _man.exists() and (BOLTZ_ROOT/IN_DIR/"_manifest.json").exists():
    _man = BOLTZ_ROOT/IN_DIR/"_manifest.json"      # 옛 위치
_stray = [f.name for f in (BOLTZ_ROOT/IN_DIR).iterdir()
          if f.is_file() and f.suffix.lower() not in (".yaml", ".yml", ".fasta", ".fa")]
if _stray:
    print(f"  ★ 입력 폴더에 yaml/fasta 가 아닌 파일이 있다: {_stray}")
    print("    boltz 는 폴더의 모든 파일을 읽으려 해서 여기서 죽는다. 옮기고 다시 돌릴 것:")
    for _s in _stray:
        print(f"      mv {BOLTZ_ROOT/IN_DIR/_s} {BOLTZ_ROOT/(IN_DIR + '_' + _s.lstrip('_'))}")
if _man.exists():
    import json as _j
    _m = _j.loads(_man.read_text())
    print(f"  입력 생성 시각 {_m.get('made_at')}  (TOPN {_m.get('TOPN')}, "
          f"mobile 반영 {_m.get('mobile_aware')})")
else:
    print("  ⚠ _manifest.json 없음 — 어느 설정으로 만든 입력인지 알 수 없다. Q5 를 다시 돌릴 것")
n_out = len(glob.glob(str(BOLTZ_ROOT/OUT_DIR/"**"/"*_model_0.cif"), recursive=True))
print("=" * 96); print("### 상태"); print("=" * 96)
print(f"  입력  {BOLTZ_ROOT/IN_DIR}  {n_in}개")
print(f"  출력  {BOLTZ_ROOT/OUT_DIR}  {n_out}개  ({n_out/max(n_in,1)*100:.1f}%)")

pid = running()
if pid:
    print(f"\n  이미 돌고 있다 — PID {pid}")
    print(f"  로그: {LOG}")
    print("  중단하려면:  import os, signal; os.kill(%d, signal.SIGTERM)" % pid)
elif n_out >= n_in and n_in:
    print("\n  이미 다 끝났다. Q6 채점 부분을 돌릴 것.")
else:
    busy = other_boltz()
    if busy:
        print("\n  ★ 이미 다른 boltz 가 돌고 있다. 띄우지 않는다.")
        for pid_, line in busy: print(f"    PID {pid_}  {line}")
        print("\n    같은 작업이면 그대로 두고 이 셀로 진행률만 본다.")
        print("    다른 작업이면 끝나기를 기다리거나, 확인 후 아래로 종료:")
        print("      import os, signal; os.kill(<PID>, signal.SIGTERM)")
        GPU = None
    else:
        if FORCE_GPU is not None:
            rows = {i: (f, tt) for i, f, tt, _u in gpu_rows()}
            f, tt = rows.get(FORCE_GPU, (None, None))
            print(f"  GPU {FORCE_GPU} 로 고정 (여유 {f} / {tt} MiB)" if f is not None
                  else f"  GPU {FORCE_GPU} 로 고정 (상태 못 읽음)")
            if f is not None and f < MIN_FREE:
                print(f"  ★ 그런데 여유가 {MIN_FREE} MiB 미만이다. 그래도 진행한다 —")
                print("    직접 지정한 것이므로 판단은 사용자 몫이다.")
            GPU = FORCE_GPU
        else:
            GPU = pick_gpu()
        if GPU is None:
            print(f"\n  ★ 여유 {MIN_FREE} MiB 이상인 GPU 가 없다. 띄우지 않는다.")
            print("    다른 작업이 끝난 뒤 이 셀을 다시 돌릴 것.")
            print("    무리해서 띄우면 이쪽도 상대쪽도 OOM 으로 같이 죽는다.")
if 'GPU' in dir() and GPU is not None and not pid and not (n_out >= n_in and n_in):
    print(f"\n  GPU {GPU} 로 띄운다")
    cmd = (f'boltz predict {IN_DIR} --out_dir {OUT_DIR} --use_msa_server '
           f'--max_msa_seqs 2048 --recycling_steps 3 --diffusion_samples 1 '
           f'--output_format mmcif --num_workers 2')
    # conda run 으로 감싼다. env 가 없으면 현재 인터프리터의 boltz 를 쓴다
    # `conda env list | grep -q` 는 grep 이 먼저 닫혀 conda 가 BrokenPipeError
    # 리포트를 통째로 뱉는다. 출력을 받아서 파이썬에서 본다.
    _envs = subprocess.run("conda env list", shell=True, text=True,
                           capture_output=True).stdout
    has_env = any(l.split()[:1] == [CONDA_ENV] for l in _envs.splitlines()
                  if l.strip() and not l.startswith("#"))
    inner = f"conda run -n {CONDA_ENV} --no-capture-output {cmd}" if has_env else cmd
    print(f"  conda env '{CONDA_ENV}': {'있음' if has_env else '없음 → 현재 환경으로'}")
    script = (f'cd {shlex.quote(str(BOLTZ_ROOT))} && '
              f'CUDA_VISIBLE_DEVICES={GPU} '
              f'PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:512 '   # 1.x/2.x 둘 다 먹는 옵션
              f'{inner}')
    LOG.parent.mkdir(parents=True, exist_ok=True)
    if LOG.exists() and LOG.stat().st_size:
        # 이어 쓰면 지난 실행의 실패 메시지가 계속 잡힌다. 돌릴 때마다 새로 연다.
        LOG.rename(LOG.with_suffix(f".log.{time.strftime('%Y%m%d_%H%M%S')}"))
    with open(LOG, "wb") as lf:
        proc = subprocess.Popen(["setsid", "bash", "-lc", script],
                                stdout=lf, stderr=lf, start_new_session=True)
    PIDF.write_text(str(proc.pid))
    time.sleep(6)
    print(f"  PID {proc.pid}  로그 {LOG}")
    if not Path(f"/proc/{proc.pid}").exists():
        print("\n  ★ 6초 만에 죽었다. 로그 끝부분:")
        print("\n".join(LOG.read_text(errors='ignore').splitlines()[-25:]))

if LOG.exists():
    _txt = LOG.read_text(errors="ignore")
    for pat, msg in [("Missing folder", "출력 폴더가 실행 중에 사라졌다 (rm -rf 를 돌면서 했을 때 난다)"),
                     ("out of memory", "OOM — 같은 카드에 다른 작업이 있었는지 볼 것"),
                     ("429", "MSA 서버 속도 제한"),
                     ("CUDA error", "CUDA 초기화/실행 실패")]:
        if pat in _txt:
            print(f"\n  ★ 로그에 '{pat}' 가 있다 — {msg}")

print("\n" + "=" * 96); print("### 로그 끝 20줄"); print("=" * 96)
print("\n".join(LOG.read_text(errors="ignore").splitlines()[-20:])
      if LOG.exists() else "  아직 로그 없음")
print("\n  ※ 이 셀을 다시 돌리면 진행률과 로그만 갱신된다 (중복 실행 안 함).")
print("    MSA 서버에서 막히면 로그에 429 나 timeout 이 뜬다. 그때는 입력을")
print("    100개씩 하위 폴더로 나눠 IN_DIR 을 바꿔가며 돌린다.")
```

---

## CELL Q7 — 단량체로 접었다는 것을 보정한다 (이량체 계면 자리)

```python
# =============================================================================
# CELL Q7 | grade C 를 버리면 CooC 형을 버린다
#
#   CooC1 의 Ni 자리는 **두 단량체의 Cys112/114 가 계면에서 만나** Cys4 를 이룬다.
#   3kji 실측이 그렇다 — ZN 하나를 A112 2.23 · B114 2.24 · A112 2.27 · B112 2.28 로
#   두 사슬이 나눠 문다.
#
#   그런데 Q5 는 후보를 **단량체 한 사슬**로 접었다. 그러면 진짜 CooC 형 단백질도
#   자기 Cys 2개밖에 못 내놓는다. 나머지 둘은 파트너 사슬에 있기 때문이다.
#     → 단량체 예측의 **grade C (Cys 1~2)** 가 곧 이량체 계면 자리의 신호다.
#     → 특히 그 2개가 CXC / CXXC 쌍이면 CooC1 과 같은 배치다.
#
#   ★ 이 문제는 기존 Track B 결과(`ni_site_grade.csv`)에도 그대로 있다.
#     거기서도 prey 는 사슬 하나였다. 아래 (3) 에서 같이 다시 본다.
#
#   할 일 — 해당하는 것들을 **동형이량체 + Ni** 로 다시 접는다.
#     사슬 A = 후보, 사슬 A' = 같은 후보, 리간드 = NI
#     계면에서 Cys4 가 만들어지면 그때 비로소 grade A 다.
# =============================================================================

# ---- (1) 단량체 결과에서 이량체 후보를 고른다 ----
GP = TBL/"seqmotif_ni_grade.csv"
if not GP.exists():
    print(f"  {GP.name} 없음 — Q6 채점을 먼저 할 것"); G = None
else:
    G = pd.read_csv(GP)
    MOT = D.set_index("protein") if "D" in dir() else pd.read_csv(TBL/"seqmotif_ranked.csv").set_index("protein")

    def pair_type(pid):
        """배위한 Cys 2개가 CXC/CXXC 쌍인가 — 서열에서 확인한다."""
        s = SEQ["BL21"].get(pid)
        if not s: return None
        for n, rx in [("CXC", r"C.C"), ("CXXC", r"C..C")]:
            if re.search(rx, s): return n
        return None

    cand = G[(G.grade.isin(["C", "B"])) & (G.Cys.between(1, 3))].copy()
    cand["pair"] = [pair_type(p) for p in cand.protein]
    cand["dimer_flag"] = cand.pair.notna().astype(int)
    DIM = cand[cand.dimer_flag == 1].sort_values(["Cys", "min_dist"],
                                                 ascending=[False, True])
    print("=" * 96); print("### (1) 이량체 재예측 대상"); print("=" * 96)
    print(f"  단량체 grade B/C 중 Cys 1~3개: {len(cand)}")
    print(f"  그중 CXC/CXXC 쌍 보유: {len(DIM)}  ← CooC1 과 같은 배치")
    print(DIM.head(20)[["protein", "grade", "Cys", "His", "min_dist", "pair",
                        "mg_status", "fold"]].to_string(index=False))

# ---- (2) 동형이량체 + Ni 입력 생성 ----
import yaml as _y
DIN = BASE/"result"/"boltz"/"inputs_seqmotif_dimer"; DIN.mkdir(parents=True, exist_ok=True)
MAXTOK = 1200          # 두 사슬 합계. 단량체 상한의 절반이 후보 길이 상한이 된다
if G is not None:
    n, skip = 0, 0
    for _, r in DIM.iterrows():
        s = SEQ["BL21"].get(r.protein)
        if not s: continue
        if 2 * len(s) > MAXTOK: skip += 1; continue
        doc = {"version": 1, "sequences": [
            {"protein": {"id": "A", "sequence": s}},
            {"protein": {"id": "B", "sequence": s}},     # 같은 서열 = 동형이량체
            {"ligand":  {"id": "L", "ccd": "NI"}}]}
        (DIN/f"{r.protein.replace('|','_')}_dimer.yaml").write_text(
            _y.safe_dump(doc, sort_keys=False))
        n += 1
    print(f"\n  이량체 입력 {n}개 생성 (2×길이 > {MAXTOK} 제외 {skip})  → {DIN}")

    # 양성대조군 — Ch CooC1 자신을 같은 방식으로. 눈금을 여기서 잡는다
    CRY = TOOLS/"workspace"/"seek_ni_insertase"/"input"/"3kji.pdb"
    AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
           "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
           "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}
    if CRY.exists():
        seq, seen = [], set()
        for l in open(CRY, errors="ignore"):
            if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "A":
                rs = l[22:26].strip()
                if rs not in seen:
                    seen.add(rs); seq.append(AA3.get(l[17:20].strip(), "X"))
        cs = "".join(seq)
        (DIN/"CTRL-CooC1_dimer.yaml").write_text(_y.safe_dump(
            {"version": 1, "sequences": [
                {"protein": {"id": "A", "sequence": cs}},
                {"protein": {"id": "B", "sequence": cs}},
                {"ligand":  {"id": "L", "ccd": "NI"}}]}, sort_keys=False))
        print(f"  ★ 양성대조군 CTRL-CooC1_dimer.yaml 생성 ({len(cs)} 잔기 ×2)")
        print("    이게 grade A 로 안 나오면 이량체 예측 자체를 못 믿는다.")

print("\n" + "=" * 96); print("### 실행 명령"); print("=" * 96)
print(f"""
cd "{BASE/'result'/'boltz'}"
CUDA_VISIBLE_DEVICES=0 boltz predict inputs_seqmotif_dimer \\
  --out_dir out_seqmotif_dimer \\
  --use_msa_server --max_msa_seqs 2048 \\
  --recycling_steps 3 --diffusion_samples 1 \\
  --output_format mmcif --num_workers 2
echo DONE_seqmotif_dimer
""")

# ---- (3) 기존 Track B 결과도 같은 눈으로 다시 본다 ----
print("=" * 96); print("### (3) 기존 ni_site_grade.csv 재검토"); print("=" * 96)
OLD = TBL/"ni_site_grade.csv"
if not OLD.exists():
    print(f"  {OLD.name} 없음")
else:
    O = pd.read_csv(OLD)
    O["key"] = O.prey.astype(str).str.split("-").str[-1]
    O["pair"] = [pair_type(k) if k in SEQ.get("BL21", {}) else None for k in O.key]
    reconsider = O[(O.grade.isin(["B", "C", "D"])) & (O.Cys.between(1, 3)) & O.pair.notna()]
    print(f"  기존 채점에서도 prey 는 사슬 하나였다. Cys 1~3개 + CXC/CXXC 보유:")
    print(f"  {len(reconsider)}개가 이량체 계면 자리일 수 있다.")
    if len(reconsider):
        print(reconsider[["prey", "grade", "Cys", "His", "min_dist", "pair",
                          "donors", "desc"]].to_string(index=False))
        reconsider.to_csv(TBL/"trackB_dimer_reconsider.csv", index=False,
                          encoding="utf-8-sig")
        print(f"\n저장: {TBL/'trackB_dimer_reconsider.csv'}")
        print("\n  ※ YeiR 가 grade D 로 내려간 것도 접촉 1개(HIS209) 때문이었다.")
        print("    YeiR·HypB·UreG·CooC 는 전부 이량체다. 단량체 예측이 이들에게")
        print("    불리하게 작동했을 수 있다. 위 목록은 이량체로 다시 물어야 한다.")
```

---

## CELL Q8 — 이량체 예측 채점 (계면을 가로지르는 Cys4)

```python
# =============================================================================
# CELL Q8 | 두 사슬이 나눠 무는 것을 세야 한다
#   단량체 채점과 다른 점 — 배위 Cys 를 **사슬별로 나눠** 보고합니다.
#   A 2개 + B 2개 = Cys4 면 그것이 CooC1 과 같은 계면 자리다.
# =============================================================================
import glob
DOUT = BASE/"result"/"boltz"/"out_seqmotif_dimer"
cifs = sorted(glob.glob(str(DOUT/"**"/"*_model_0.cif"), recursive=True))
print(f"  이량체 예측 결과 {len(cifs)}개")

rows = []
for cif in cifs:
    A = cif_atoms(cif)
    ni = [a for a in A if a[0] == "NI"]
    if not ni: continue
    _, nch, _, _, x, y, z = ni[0]
    donors = []
    for comp, ch, seq, at, ax, ay, az in A:
        if comp == "NI" or ch == nch: continue
        d = math.dist((x, y, z), (ax, ay, az))
        if d <= TIGHT and (at not in BACKBONE or comp in ("ASP", "GLU", "HIS", "CYS")):
            donors.append((round(d, 2), comp, ch, seq, at))
    donors = sorted({(c, ch, s): (d, c, ch, s, a) for d, c, ch, s, a in donors}.values())
    cys = [(ch, s, d) for d, c, ch, s, _a in donors if c == "CYS"]
    by_ch = Counter(ch for ch, _s, _d in cys)
    nhis = sum(1 for d, c, ch, s, _a in donors if c == "HIS")
    ncys = len(cys)
    bridged = len(by_ch) >= 2           # 두 사슬이 같이 물었는가
    if   ncys >= 4 and bridged: g, why = "A_bridged", "두 사슬이 나눠 문 Cys4 — CooC1 형"
    elif ncys >= 4:             g, why = "A", "한 사슬 안의 Cys4"
    elif ncys == 3:             g, why = "B", "Cys3"
    elif ncys == 2 and nhis:    g, why = "B", "Cys2His"
    elif ncys or nhis >= 3:     g, why = "C", "His 중심 / Cys 부족"
    else:                       g, why = "D", "산소 공여체 위주"
    rows.append({"grade": g, "protein": Path(cif).stem.replace("_model_0", "")
                                  .replace("_dimer", ""),
                 "Cys": ncys, "His": nhis, "Cys_by_chain": dict(by_ch), "bridged": bridged,
                 "min_dist": donors[0][0] if donors else None,
                 "donors": " ".join(f"{c}{ch}{s}" for _d, c, ch, s, _a in donors)[:70],
                 "why": why})

if rows:
    E = pd.DataFrame(rows).sort_values(["grade", "Cys", "min_dist"],
                                       ascending=[True, False, True])
    E.to_csv(TBL/"seqmotif_dimer_grade.csv", index=False, encoding="utf-8-sig")
    print("\n" + "=" * 96); print("### 이량체 등급"); print("=" * 96)
    print(E.groupby("grade").size().to_string())

    ctrl = E[E.protein.str.contains("CTRL-CooC1", na=False)]
    print("\n  ★ 양성대조군 CooC1:")
    print(ctrl.to_string(index=False) if len(ctrl) else "    결과 없음 — 대조군을 먼저 볼 것")
    if len(ctrl) and not str(ctrl.iloc[0].grade).startswith("A"):
        print("\n    ⚠ 대조군이 A 가 아니다. 이량체 예측 자체가 이 자리를 못 만든다는")
        print("      뜻이므로, 아래 결과를 근거로 쓰면 안 된다.")

    print("\n" + "=" * 96); print("### A-계면 (두 사슬이 나눠 문 Cys4)"); print("=" * 96)
    hot = E[E.grade == "A_bridged"]
    print(hot.to_string(index=False) if len(hot) else "  없음")

    # 단량체 결과와 나란히
    if GP.exists():
        M1 = pd.read_csv(GP)[["protein", "grade", "Cys"]] \
               .rename(columns={"grade": "grade_monomer", "Cys": "Cys_monomer"})
        CMP = E[["protein", "grade", "Cys", "bridged"]] \
                .rename(columns={"grade": "grade_dimer", "Cys": "Cys_dimer"}) \
                .merge(M1, on="protein", how="left")
        CMP.to_csv(TBL/"seqmotif_mono_vs_dimer.csv", index=False, encoding="utf-8-sig")
        print("\n" + "=" * 96); print("### 단량체 → 이량체 변화"); print("=" * 96)
        print(CMP.sort_values("grade_dimer").to_string(index=False))
        up = CMP[(CMP["grade_monomer"].isin(["C", "D"])) & (CMP["grade_dimer"].str.startswith("A"))]
        print(f"\n  ★ 단량체에서 C/D 였다가 이량체에서 A 로 올라온 것: {len(up)}개")
        print("    이것이 단량체 예측이 놓치고 있던 것이다.")
else:
    print("  아직 결과가 없다. 위 명령을 먼저 돌릴 것.")
```
