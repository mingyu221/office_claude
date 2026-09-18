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

def struct_of(pid):
    xw = TBL/"id_crosswalk_struct_to_genbank.csv"
    if not hasattr(struct_of, "_m"):
        struct_of._m = {}
        if xw.exists():
            d = pd.read_csv(xw).dropna(subset=["protein"])
            for t, p in zip(d.tid, d.protein):
                struct_of._m[str(p).split(",")[0]] = Path(str(t))
    p = struct_of._m.get(pid)
    if p and p.exists(): return p
    if p:
        for s in STRUCT.values():
            if (s/p.name).exists(): return s/p.name
    return None

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

D = CAND.merge(V[["protein", "mg_status", "bl21_only"]], on="protein", how="left") \
        .merge(Q, on="protein", how="left")

def score(r):
    p, why = 0, []
    if r.cys4_cand == 1:                       p += 3; why.append("+3 Cys4후보")
    if r.bl21_only == 1:                       p += 3; why.append("+3 MG없음")
    elif r.mg_status == "모티프 달라짐":        p += 1; why.append("+1 모티프차이")
    if pd.notna(r.fold):                       p += 2; why.append(f"+2 {r.fold}폴드")
    if (r.n_CXC or 0) or (r.n_CXXC or 0):      p += 2; why.append("+2 CXC/CXXC")
    if r.len <= 250:                           p += 1; why.append("+1 소형")
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
                  "fold", "fold_tm", "plddt_motif", "mg_status", "desc"]].to_string(index=False))
print(f"\n저장: {TBL/'seqmotif_ranked.csv'}")

# ---- Boltz 입력: 후보 단독 + Ni ----
import yaml as _y
BIN = BASE/"result"/"boltz"/"inputs_seqmotif"; BIN.mkdir(parents=True, exist_ok=True)
MAXLEN = 1200
SEL = D.head(TOPN)
n, skip = 0, 0
for _, r in SEL.iterrows():
    s = SEQ["BL21"].get(r.protein)
    if not s: continue
    if len(s) > MAXLEN: skip += 1; continue
    # 사슬 A = 후보 단백질, 리간드 B = Ni. ChCODH2 는 넣지 않는다
    doc = {"version": 1, "sequences": [
        {"protein": {"id": "A", "sequence": s}},
        {"ligand":  {"id": "B", "ccd": "NI"}}]}
    (BIN/f"{r.protein.replace('|','_')}.yaml").write_text(_y.safe_dump(doc, sort_keys=False))
    n += 1
print(f"\n  Boltz 입력 {n}개 생성 (길이 초과 제외 {skip})  → {BIN}")
print(f"  총 잔기 {int(SEL.head(n).len.sum()):,} — 기존 방식이면 여기에 630×{n} 이 더 붙는다")
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
    if   ncys >= 4 and bridged: g, why = "A-계면", "두 사슬이 나눠 문 Cys4 — CooC1 형"
    elif ncys >= 4:             g, why = "A", "한 사슬 안의 Cys4"
    elif ncys == 3:             g, why = "B", "Cys3"
    elif ncys == 2 and nhis:    g, why = "B", "Cys2His"
    elif ncys or nhis >= 3:     g, why = "C", "His 중심 / Cys 부족"
    else:                       g, why = "D", "산소 공여체 위주"
    rows.append({"grade": g, "protein": Path(cif).stem.replace("_model_0", "")
                                  .replace("_dimer", ""),
                 "Cys": ncys, "His": nhis, "사슬별Cys": dict(by_ch), "계면": bridged,
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
    hot = E[E.grade == "A-계면"]
    print(hot.to_string(index=False) if len(hot) else "  없음")

    # 단량체 결과와 나란히
    if GP.exists():
        M1 = pd.read_csv(GP)[["protein", "grade", "Cys"]] \
               .rename(columns={"grade": "단량체등급", "Cys": "단량체Cys"})
        CMP = E[["protein", "grade", "Cys", "계면"]] \
                .rename(columns={"grade": "이량체등급", "Cys": "이량체Cys"}) \
                .merge(M1, on="protein", how="left")
        CMP.to_csv(TBL/"seqmotif_mono_vs_dimer.csv", index=False, encoding="utf-8-sig")
        print("\n" + "=" * 96); print("### 단량체 → 이량체 변화"); print("=" * 96)
        print(CMP.sort_values("이량체등급").to_string(index=False))
        up = CMP[(CMP.단량체등급.isin(["C", "D"])) & (CMP.이량체등급.str.startswith("A"))]
        print(f"\n  ★ 단량체에서 C/D 였다가 이량체에서 A 로 올라온 것: {len(up)}개")
        print("    이것이 단량체 예측이 놓치고 있던 것이다.")
else:
    print("  아직 결과가 없다. 위 명령을 먼저 돌릴 것.")
```
