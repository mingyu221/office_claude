# 모티프 매칭 표 — 독립 노트북

`ChCODH2_screen_antigravity.ipynb` 와 별개로 돌아간다. 앞 셀에 의존하지 않고
경로를 스스로 잡으므로, 커널을 새로 띄워 M1 → M2 → M3 순으로만 돌리면 된다.

산출물 세 개
- `result/table/motif_counts_by_strain.csv`
- `result/table/motif_metal_and_atp.csv`
- `result/table/motif_metal_ntp.csv`

---

## CELL M1 — 설정 (경로 · 원본 확인)

```python
# =============================================================================
# CELL M1 | 독립 실행용 설정
#   본 노트북(ChCODH2_screen_antigravity)의 CELL 01 에 의존하지 않는다.
#   경로가 서버에서 바뀌면 여기만 고치면 된다.
# =============================================================================
import re, glob, subprocess
from pathlib import Path
import pandas as pd

BASE    = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS   = Path("/mnt/af2results/mingyu")
INHOUSE = TOOLS/"database"/"bacteriaDB"/"inhouseDB"
PLIST   = TOOLS/"database"/"protein_list"

TBL = BASE/"result"/"table"
FDD = BASE/"result"/"folddisco"

def _pick(*c):
    for x in c:
        if x.exists(): return x
    return c[0]

FAA = {"BL21":   _pick(INHOUSE/"bl21_db_match_qjz.faa", PLIST/"bl21_de3_protein.faa"),
       "Y19":    _pick(INHOUSE/"y19_db_match.faa",      PLIST/"y19_protein.faa"),
       "MG1655": PLIST/"mg1655_protein.faa"}
STRUCT = {"BL21": TOOLS/"database"/"bacteriaDB"/"structures_UP000503272",
          "Y19":  TOOLS/"database"/"bacteriaDB"/"structures_UP000034085"}

# 구조 DB 를 만든 쪽이 준 정본 매핑. 시트당 한 균주이고, 각 구조에
#   UniParc accession / UniProtKB accessions / Genbank / Selected structure
# 가 들어 있다. Selected structure 가 "ColabFold predicted" 면 자체 예측이고
# 파일명이 cf_{UniParc}.pdb, 아니면 그 값이 AFDB accession 이라 AF-{그것}-F1 이다.
# 이 파일이 있으면 stem -> Genbank 가 100% 풀린다 (실측: 히트 2,232건 실패 0).
XLSX = _pick(BASE/"input"/"external"/"structure_accessions.xlsx",
             TOOLS/"database"/"bacteriaDB"/"structure_accessions.xlsx",
             PLIST/"structure_accessions.xlsx")

SETS    = {"metal": "metal_run01_{s}.tsv", "atp": "atp_run02_{s}.tsv"}
STRAINS = ["BL21", "MG1655", "Y19"]
# 구조 DB 에 이름 규칙이 두 가지 섞여 있다.
#   AF-{UniProt}-F1-model_v6.cif   AlphaFold DB
#   cf_{UniParc}.pdb               AlphaFold DB 에 없어 자체 예측한 것
# accession 만 뽑으면 후자가 통째로 빠진다 (BL21 66개, Y19 15개).
# 그러므로 키는 파일명 stem 으로 잡고, accession 은 UniProt 조회용으로만 쓴다.
ACC_RE  = re.compile(r"AF-([A-Za-z0-9]+)-F1")
def stem_of(tid):  return Path(str(tid)).stem
def acc_of(tid):
    m = ACC_RE.search(str(tid));  return m.group(1) if m else None

print("=" * 96); print("### 경로 점검"); print("=" * 96)
ok = True
for s in STRAINS:
    for k in SETS:
        p = FDD/SETS[k].format(s=s)
        print(f"  {k:6s} {s:8s} {'O' if p.exists() else 'X'}  {p}")
        ok &= p.exists()
for s in STRAINS:
    print(f"  faa    {s:8s} {'O' if FAA[s].exists() else 'X'}  {FAA[s]}")
for s, d in STRUCT.items():
    print(f"  struct {s:8s} {'O' if d.is_dir() else 'X'}  {d}")
print(f"  table  {'O' if TBL.is_dir() else 'X'}  {TBL}")
print(f"  xlsx   {'O' if XLSX.exists() else 'X'}  {XLSX}")
if not XLSX.exists():
    print("     ※ 없으면 id_crosswalk_struct_to_genbank.csv 로 넘어간다.")
    print("       그 경우 cf_* 구조의 이름이 안 붙을 수 있다.")
if not ok:
    print("\n⚠ X 가 있으면 그 경로부터 고칠 것. 아래 셀은 원본 TSV 가 있어야 돈다.")
```

---

## CELL M2 — 원본 읽기 + 이름 붙이기

```python
# =============================================================================
# CELL M2 | TSV 를 구조 단위로 정리하고 이름을 붙인다
#   한 행이 한 구조이고 한 단백질이다 — 중복은 없다. 다만 구조 파일 이름이
#   AF-*(AlphaFold DB) 와 cf_*(자체 예측) 두 가지라, accession 으로 키를 잡으면
#   후자가 조용히 빠진다. 키는 stem 으로 잡고 둘을 각각 세어 표에 남긴다.
#   이름은 crosswalk -> 프로테옴 FASTA -> UniProt 순으로 찾고, 조회 결과는
#   캐시에 저장해 두 번 묻지 않는다.
# =============================================================================
def read_set(kind, strain):
    p = FDD/SETS[kind].format(s=strain)
    if not p.exists(): return None, 0
    d = pd.read_csv(p, sep="\t")
    d.columns = [c.strip().lstrip("#") for c in d.columns]
    n_rows = len(d)
    d["stem"] = d.tid.map(stem_of)
    d["acc"]  = d.tid.map(acc_of)
    d["src"]  = d.acc.notna().map({True: "AFDB", False: "ColabFold"})
    if "idf" in d: d = d.sort_values("idf", ascending=False)
    return d.drop_duplicates("stem").set_index("stem"), n_rows

DATA, NROW = {}, {}
print("=" * 96); print("### 원본"); print("=" * 96)
for k in SETS:
    for s in STRAINS:
        DATA[(k, s)], NROW[(k, s)] = read_set(k, s)
        d = DATA[(k, s)]
        if d is None: print(f"  {k:6s} {s:8s} 없음"); continue
        c = d.src.value_counts().to_dict()
        print(f"  {k:6s} {s:8s} {NROW[(k,s)]:5d}행 / {len(d):5d}구조  "
              f"(AFDB {c.get('AFDB',0)}, ColabFold {c.get('ColabFold',0)})")
        if NROW[(k, s)] != len(d):
            print(f"      ⚠ 행 {NROW[(k,s)]} ≠ 구조 {len(d)} — 같은 구조가 두 번 나왔다")
assert all(DATA[(k, s)] is not None for k in SETS for s in STRAINS), "TSV 누락 — M1 점검"

# ---- 프로테옴 크기 ----
NPROT = {}
for s in STRAINS:
    d = STRUCT.get(s)
    # AF-* 만 세면 cf_* (ColabFold) 가 빠진다 — BL21 250, Y19 21 개.
    NPROT[s] = (len([f for f in d.iterdir() if f.name.startswith(("AF-", "cf_"))])
                if d and d.is_dir()
                else sum(1 for l in open(FAA[s], errors="ignore") if l.startswith(">")))
print("\n  프로테옴:", {s: NPROT[s] for s in STRAINS})

# ---- 이름 ----
HDR = {}
for s in STRAINS:
    if not FAA[s].exists(): continue
    for l in open(FAA[s], errors="ignore"):
        if l.startswith(">"):
            h = l[1:].rstrip(); k = h.split()[0]
            HDR[k] = h[len(k):].strip()

XW = {}       # stem -> GenBank protein ID
SRC_INFO = {}  # stem -> "AFDB" / "ColabFold"
if XLSX.exists():
    xl = pd.ExcelFile(XLSX)
    for sh in xl.sheet_names:
        d = xl.parse(sh)
        n_cf = 0
        for _, r in d.iterrows():
            gb = str(r["Genbank"]).split(",")[0].strip()
            upi = str(r["UniParc accession"]).strip()
            sel = str(r["Selected structure"]).strip()
            XW[f"cf_{upi}"] = gb; SRC_INFO[f"cf_{upi}"] = "ColabFold"
            if sel and sel != "ColabFold predicted":
                for v in ("v6", "v5", "v4", "v3"):
                    XW[f"AF-{sel}-F1-model_{v}"] = gb
                    SRC_INFO[f"AF-{sel}-F1-model_{v}"] = "AFDB"
            else:
                n_cf += 1
        print(f"  xlsx [{sh[:26]}] 구조 {len(d)} (ColabFold {n_cf}, AFDB {len(d)-n_cf})")
_x = TBL/"id_crosswalk_struct_to_genbank.csv"
if _x.exists():
    d = pd.read_csv(_x).dropna(subset=["protein"])
    for k, v in d.set_index("tid")["protein"].to_dict().items():
        XW.setdefault(stem_of(k), str(v).split(",")[0])   # 엑셀이 우선
print(f"  매핑 {len(XW)}건, FASTA 헤더 {len(HDR)}건")

# 히트가 이 매핑으로 다 풀리는지 먼저 본다 — 안 풀리면 이름 없는 행이 생긴다
for k in SETS:
    for s in STRAINS:
        idx = list(DATA[(k, s)].index)
        miss = [i for i in idx if i not in XW]
        if miss:
            print(f"  ⚠ {k} {s}: {len(miss)}/{len(idx)} 미매핑 예) {miss[:2]}")

UNI = {}
CACHE = TBL/"uniprot_name_cache.tsv"
if CACHE.exists():
    d = pd.read_csv(CACHE, sep="\t")
    UNI.update(dict(zip(d.acc.astype(str), d.name.astype(str))))
    print(f"  이름 캐시 {len(UNI)}건 재사용")

# crosswalk 로 이름이 안 붙은 것만 UniProt 에 묻는다 (AF 형식만 accession 이 있다)
NEED = sorted({DATA[("metal", s)].loc[k, "acc"] for s in STRAINS
               for k in DATA[("metal", s)].index
               if k not in XW and pd.notna(DATA[("metal", s)].loc[k, "acc"])
               and DATA[("metal", s)].loc[k, "acc"] not in UNI})
if NEED:
    import urllib.request, urllib.parse, time
    print(f"\n  UniProt 조회 {len(NEED)}개")
    got = {}
    for i in range(0, len(NEED), 40):
        q = " OR ".join(f"accession:{a}" for a in NEED[i:i+40])
        url = ("https://rest.uniprot.org/uniprotkb/search?format=tsv&size=500"
               "&fields=accession,protein_name,gene_primary,keyword"
               f"&query={urllib.parse.quote(q)}")
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                ls = r.read().decode().splitlines()
            h = ls[0].split("\t")
            for l in ls[1:]:
                rr = dict(zip(h, l.split("\t")))
                got[rr.get("Entry", "")] = " | ".join(
                    x for x in [rr.get("Protein names", ""),
                                rr.get("Gene Names (primary)", ""),
                                rr.get("Keywords", "")] if x)
            print(f"    batch {i//40+1}: 누적 {len(got)}")
        except Exception as e:
            print(f"    ⚠ 실패 {e} — 이름 없이 진행"); break
        time.sleep(0.4)
    if got:
        UNI.update(got)
        TBL.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"acc": list(UNI), "name": [UNI[a] for a in UNI]}
                     ).to_csv(CACHE, sep="\t", index=False)
        print(f"  캐시 저장: {CACHE}")

def name_of(stem, acc=None):
    """stem 으로 crosswalk -> FASTA, 없으면 accession 으로 UniProt."""
    pid = XW.get(stem)
    if pid and pid in HDR: return pid, HDR[pid]
    if acc and acc in UNI: return pid or acc, UNI[acc]
    if pid:                return pid, ""
    return acc or stem, ""

_n = sum(1 for s in STRAINS for k in DATA[("metal", s)].index
         if name_of(k, DATA[("metal", s)].loc[k, "acc"])[1])
print(f"\n  metal 히트 중 이름이 붙은 것 {_n} / "
      f"{sum(len(DATA[('metal', s)]) for s in STRAINS)}")
```

---

## CELL M3 — 표 1 · 2 · 3

```python
# =============================================================================
# CELL M3 | 표 세 개
#   표1  균주별 매칭 수 + 프로테옴 대비 비율
#   표2  metal ∩ ATP  (+ 세 균주 공통 직교체 묶음 표시)
#   표3  metal 중 NTP 소모 — 주석 기반이므로 근거 단어를 evidence 열에 남긴다
# =============================================================================
pd.set_option("display.max_rows", None); pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 46)
TBL.mkdir(parents=True, exist_ok=True)

# ---------------------------- 표 1 ----------------------------
print("=" * 96); print("### 표 1 — 균주별 모티프 매칭"); print("=" * 96)
def _cnt(k, s):
    d = DATA[(k, s)]; c = d.src.value_counts().to_dict()
    return len(d), c.get("AFDB", 0), c.get("ColabFold", 0)
rows = []
for s in STRAINS:
    m, m_af, m_cf = _cnt("metal", s); a, a_af, a_cf = _cnt("atp", s)
    rows.append({"균주": s, "metal": m, "metal(AFDB/cf)": f"{m_af}/{m_cf}",
                 "ATP": a, "ATP(AFDB/cf)": f"{a_af}/{a_cf}", "프로테옴": NPROT[s],
                 "metal %": round(100*m/NPROT[s], 2), "ATP %": round(100*a/NPROT[s], 2)})
T1 = pd.DataFrame(rows)
print(T1.to_string(index=False))
T1.to_csv(TBL/"motif_counts_by_strain.csv", index=False, encoding="utf-8-sig")
print("\n  한 행이 한 구조이고 한 단백질이다. AFDB/cf 는 구조의 출처다 —")
print("  AlphaFold DB 에 없어 ColabFold 로 접은 것이 BL21 250, Y19 21 개 있다.")
print("  둘 다 structure_accessions.xlsx 로 Genbank 에 붙으므로 이름은 다 나온다.")
print("  ATP 보유율이 20%대면 그 축은 선별력이 없다.")

# ---------------------------- 표 2 ----------------------------
print("\n" + "=" * 96); print("### 표 2 — metal ∩ ATP"); print("=" * 96)
rows = []
for s in STRAINS:
    M_, A_ = DATA[("metal", s)], DATA[("atp", s)]
    for a in sorted(set(M_.index) & set(A_.index),
                    key=lambda x: -float(M_.loc[x, "idf"])):
        pid, desc = name_of(a, M_.loc[a, "acc"])
        rows.append({"strain": s, "struct": a, "protein": pid,
                     "nres": M_.loc[a, "nres"], "plddt": M_.loc[a, "plddt"],
                     "metal_rmsd": M_.loc[a, "min_rmsd"],
                     "metal_res": M_.loc[a, "matching_residues"],
                     "atp_rmsd": A_.loc[a, "min_rmsd"],
                     "atp_res": A_.loc[a, "matching_residues"],
                     "desc": desc[:70]})
T2 = pd.DataFrame(rows)
print(T2.drop(columns=["atp_res"]).to_string(index=False))
print("\n  균주별: " + " / ".join(f"{s} {int((T2.strain == s).sum())}" for s in STRAINS))

# 길이와 잔기 번호가 세 균주에서 같으면 직교체다 — 균주를 가르지 못한다
# 길이를 키에 넣으면 직교체가 한두 잔기 차이로 안 묶인다 (190 vs 188).
# 잔기 시그니처만 쓴다 — 같은 자리를 같은 번호로 맞춘 것이 직교체다.
T2["_key"] = T2.metal_res.astype(str).str.split(":").str[0]
g = T2.groupby("_key").strain.nunique()
shared = set(g[g >= 3].index)
print(f"\n  세 균주 공통 묶음 {len(shared)}개 (길이·잔기 동일) — 균주 구분 불가")
for k in sorted(shared):
    sub = T2[T2._key == k]
    print(f"    {k:22s} " + " · ".join(f"{r.strain}:{r.protein}" for r in sub.itertuples()))
only = T2[(~T2._key.isin(shared)) & (T2.strain == "BL21")]
print(f"\n  ★ BL21 에만 남는 것 {len(only)}개 — 여기가 볼 자리다")
if len(only):
    print(only[["struct", "protein", "nres", "plddt", "metal_res", "desc"]].to_string(index=False))
T2.drop(columns=["_key"]).to_csv(TBL/"motif_metal_and_atp.csv",
                                 index=False, encoding="utf-8-sig")

# ---------------------------- 표 3 ----------------------------
print("\n" + "=" * 96); print("### 표 3 — metal 매칭 중 NTP 소모"); print("=" * 96)
NTP_RULES = [
    ("가수분해 확실", r"\bATPase\b|\bGTPase\b|ATP-dependent|GTP-dependent|"
                      r"\bhelicase\b|\btopoisomerase\b|\bgyrase\b|\bAAA\b"),
    ("전이·리가제",   r"\bkinase\b|\bligase\b|\bsynthetase\b|\bpolymerase\b|"
                      r"adenylyltransferase|nucleotidyltransferase"),
    ("결합만 표기",   r"ATP-binding|GTP-binding|Nucleotide-binding|P-loop|Walker"),
]
rows = []
for s in STRAINS:
    M_ = DATA[("metal", s)]; atp_hit = set(DATA[("atp", s)].index)
    for a in M_.index:
        pid, desc = name_of(a, M_.loc[a, "acc"])
        cat, ev = "해당 없음", ""
        for label, pat in NTP_RULES:
            m = re.search(pat, desc, re.I)
            if m: cat, ev = label, m.group(0); break
        rows.append({"strain": s, "struct": a, "protein": pid, "ntp": cat,
                     "evidence": ev, "ATP모티프": "O" if a in atp_hit else ".",
                     "nres": M_.loc[a, "nres"], "plddt": M_.loc[a, "plddt"],
                     "metal_rmsd": M_.loc[a, "min_rmsd"], "desc": desc[:70]})
T3 = pd.DataFrame(rows)
_named = int((T3.desc.astype(str).str.len() > 0).sum())
print(f"  이름이 붙은 것 {_named} / {len(T3)}")
if _named < len(T3):
    print("  ⚠ 이름이 빈 것은 어디에서도 못 찾은 것이다. 판정에서 빠진다.")
ORD = {"가수분해 확실": 0, "전이·리가제": 1, "결합만 표기": 2, "해당 없음": 3}
T3["_o"] = T3.ntp.map(ORD)
T3 = T3.sort_values(["_o", "strain", "metal_rmsd"]).drop(columns=["_o"])
print("\n  분류 x 균주")
print(pd.crosstab(T3.ntp, T3.strain).reindex(list(ORD)).fillna(0).astype(int).to_string())
hot = T3[T3.ntp != "해당 없음"]
print(f"\n  NTP 소모로 분류된 {len(hot)}개")
print(hot[["strain", "struct", "protein", "ntp", "evidence", "ATP모티프",
           "nres", "plddt", "desc"]].to_string(index=False) if len(hot) else "  없음")
T3.to_csv(TBL/"motif_metal_ntp.csv", index=False, encoding="utf-8-sig")

# ---------------------------- 요약 한 장 ----------------------------
# 네 열이 한 줄에 있어야 비교가 된다. 다만 앞의 셋과 마지막은 성격이 다르다 —
# metal·ATP·metal+ATP 는 구조 기하의 측정이고, metal+NTP 는 주석 문자열이다.
print("\n" + "=" * 96); print("### 요약"); print("=" * 96)
S = []
for s in STRAINS:
    m = set(DATA[("metal", s)].index); a = set(DATA[("atp", s)].index)
    ntp = set(T3[(T3.strain == s) & (T3.ntp != "해당 없음")].struct)
    S.append({"균주": s, "metal motif": len(m), "ATP motif": len(a),
              "metal + ATP motif": len(m & a), "metal + NTP motif": len(ntp),
              "셋 다": len(m & a & ntp)})
SUM = pd.DataFrame(S)
print(SUM.to_string(index=False))
SUM.to_csv(TBL/"motif_summary.csv", index=False, encoding="utf-8-sig")
print("\n  metal + ATP  : 구조 기하 둘이 한 단백질에서 겹친 것 — 측정이다.")
print("  metal + NTP  : 이름에 ATPase/GTPase/kinase/ligase 가 있는 것 — 주석이다.")
print("                 표3 의 evidence 열로 어느 단어에 걸렸는지 확인할 것.")
print("  셋 다        : 위 둘을 모두 만족. 여기가 가장 좁은 교집합이다.")
_ntp_break = pd.crosstab(T3[T3.ntp != "해당 없음"].ntp, T3[T3.ntp != "해당 없음"].strain)
print("\n  NTP 내역")
print(_ntp_break.to_string())

# ------------------------ 요약표의 구성원 ------------------------
# 요약표의 각 칸이 어느 단백질들인지 한 파일에 담는다.
# metal 히트 한 줄마다 어느 칸에 드는지를 O/. 로 표시한다 — ATP 열만
# 예외다 (그쪽은 metal 에 안 걸린 3,000여 개가 따로 있다).
print("\n" + "=" * 96); print("### 요약표 구성원"); print("=" * 96)
MEM = T3.copy()
MEM["metal"] = "O"
MEM["metal+ATP"] = MEM.ATP모티프.map(lambda v: "O" if v == "O" else ".")
MEM["metal+NTP"] = MEM.ntp.map(lambda v: "." if v == "해당 없음" else "O")
MEM["셋다"] = [("O" if a == "O" and n == "O" else ".")
               for a, n in zip(MEM["metal+ATP"], MEM["metal+NTP"])]
COLS = ["strain", "struct", "protein", "metal", "metal+ATP", "metal+NTP", "셋다",
        "ntp", "evidence", "nres", "plddt", "metal_rmsd", "desc"]
MEM = MEM[[c for c in COLS if c in MEM.columns]]
MEM = MEM.sort_values(["셋다", "metal+NTP", "metal+ATP", "strain"],
                      ascending=[False, False, False, True])
MEM.to_csv(TBL/"motif_summary_members.csv", index=False, encoding="utf-8-sig")
print(f"  {len(MEM)}행 -> {TBL/'motif_summary_members.csv'}")
print("\n  칸별 확인 (요약표와 맞아야 한다)")
chk = pd.DataFrame([{"균주": s,
                     "metal": int((MEM.strain == s).sum()),
                     "metal+ATP": int(((MEM.strain == s) & (MEM["metal+ATP"] == "O")).sum()),
                     "metal+NTP": int(((MEM.strain == s) & (MEM["metal+NTP"] == "O")).sum()),
                     "셋다": int(((MEM.strain == s) & (MEM["셋다"] == "O")).sum())}
                    for s in STRAINS])
print(chk.to_string(index=False))
print("\n  ★ '셋 다' 에 드는 것")
top = MEM[MEM["셋다"] == "O"]
print(top[["strain", "protein", "ntp", "evidence", "nres", "desc"]].to_string(index=False)
      if len(top) else "    없음")
print("\n  ATP 열의 전체 목록은 이 파일에 없다 — metal 에 안 걸린 3,000여 개가")
print("  따로 있기 때문이다. 그쪽은 folddisco/atp_run02_{균주}.tsv 원본을 볼 것.")

# ------------------------ 균주 특이성 ------------------------
# 균주 간에 단백질 ID 를 직접 겹칠 수 없다 (QJZ* / P* / AKE* 는 다른 체계다).
# 같은 자리를 같은 잔기 번호로 맞춘 히트를 직교체로 본다 — 길이는 키에서 뺀다
# (직교체는 한두 잔기 길이가 다를 수 있다).
#
# ★ 이 표의 '없음' 은 '그 균주에 그 단백질이 없다' 가 아니라
#   '그 균주의 히트 목록에 같은 시그니처가 없다' 이다. 셋이 섞여 있다.
#     (1) 진짜 없음
#     (2) 있는데 구조 버전 차이로 검출 안 됨 (BL21 v6 / MG1655 v4)
#     (3) 있는데 질의 민감도 밖 (Ch CooC1 질의가 Y19 CooC 를 놓치는 것과 같은 이유)
#   확정하려면 mmseqs 로 직접 묻고 게놈까지 봐야 한다.
print("\n" + "=" * 96); print("### 균주 특이성 (잔기 시그니처 기준)"); print("=" * 96)
sig = []
for s in STRAINS:
    d = DATA[("metal", s)]
    for k in d.index:
        pid, _ = name_of(k, d.loc[k, "acc"])
        sig.append({"strain": s, "pid": pid,
                    "sig": str(d.loc[k, "matching_residues"]).split(":")[0]})
G = pd.DataFrame(sig)
piv = G.pivot_table(index="sig", columns="strain", values="pid",
                    aggfunc=lambda v: ",".join(sorted(v))).reindex(columns=STRAINS)
piv["pattern"] = piv.notna().apply(
    lambda r: "".join("BMY"[i] if r.iloc[i] else "." for i in range(3)), axis=1)

def _n(pat, s):
    sub = piv[piv.pattern.isin(pat)][s].dropna()
    return int(sum(len(str(v).split(",")) for v in sub))
rows = []
for i, s in enumerate(STRAINS):
    tag = "BMY"[i]
    tot  = len(DATA[("metal", s)])
    allp = _n(["BMY"], s)
    uniq = _n(["".join(c if c == tag else "." for c in "BMY")], s)
    rows.append({"균주": s, "metal": tot, "세 균주 공통": allp,
                 "BL21+Y19 (MG없음)": _n(["B.Y"], s) if s in ("BL21", "Y19") else None,
                 "BL21+MG (Y19없음)": _n(["BM."], s) if s in ("BL21", "MG1655") else None,
                 "MG+Y19 (BL21없음)": _n([".MY"], s) if s in ("MG1655", "Y19") else None,
                 "단독": uniq})
SP = pd.DataFrame(rows)
print(SP.to_string(index=False))
SP.to_csv(TBL/"motif_strain_specificity.csv", index=False, encoding="utf-8-sig")

print("\n  ★ BL21+Y19 에 있고 MG1655 에 없는 자리 — 가설이 예측하는 칸")
hot = piv[piv.pattern == "B.Y"]
print(hot[STRAINS].to_string() if len(hot) else "    없음")
print("\n  BL21 단독")
b_only = piv[piv.pattern == "B.."]
print(b_only[["BL21"]].to_string() if len(b_only) else "    없음")
print("\n  ※ 위 두 목록의 '없음' 은 미검출이지 부재가 아니다. 확정하려면")
print("    mmseqs 로 MG1655 에 직접 묻고(프로테옴), 그다음 게놈까지 봐야 한다.")
print("    cys4_genome_check.csv 가 하던 일을 이 목록에 적용하면 된다.")

print("\n### 저장")
for f in ["motif_counts_by_strain.csv", "motif_metal_and_atp.csv", "motif_metal_ntp.csv",
          "motif_summary.csv", "motif_summary_members.csv", "motif_strain_specificity.csv"]:
    print(f"  {TBL/f}")
print("\n### 읽는 법")
print("  표3 의 분류는 주석 문자열에서 나온 것이지 측정이 아니다.")
print("  evidence 열의 단어를 보고 직접 걸러낼 것. ATP모티프 O 는 구조 쪽 근거가")
print("  하나 더 있다는 뜻이고, 이름이 빈 행은 아무 판정도 받지 않았다는 뜻이다.")
```

---

## CELL M4 — "MG1655 에 없음" 을 검증한다

```python
# =============================================================================
# CELL M4 | 미검출과 부재를 가른다
#   M3 의 'BL21+Y19 (MG없음)' 과 'BL21 단독' 은 MG1655 히트 목록에 같은
#   시그니처가 없다는 뜻일 뿐이다. 세 가지가 섞여 있다 —
#     진짜 없음 / 구조 버전 차이로 미검출 / 질의 민감도 밖.
#
#   가르는 법: MG1655 프로테옴에 직교체가 있는지 묻고, 있으면 그 직교체에
#   같은 Cys 쌍이 보존돼 있는지 정렬로 확인한다.
#     Cys 보존됨  -> 미검출이다. 이 자리는 균주 변별에 못 쓴다
#     직교체 없음 -> 게놈까지 확인한다 (어노테이션 누락·pseudogene 배제)
#
#   ★ HslO(QJZ13803.1) 가 목록에 있다. K-12 에 hslO 는 분명히 있으므로
#     적어도 하나는 미검출이다. 그 하나가 나머지의 성격을 시사한다.
# =============================================================================
import subprocess, urllib.request
_tmp = BASE/"tmp"/"m4"; _tmp.mkdir(parents=True, exist_ok=True)
def sh(c, quiet=False):
    r = subprocess.run(c, shell=True, text=True, stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT)
    if not quiet: print(r.stdout.rstrip())
    return r.stdout

# 검증 대상: BL21 쪽에서 MG1655 에 짝이 없는 것 전부
TGT = {}
for p_ in ("B.Y", "B.."):
    sub = piv[piv.pattern == p_]
    for s_, v in sub["BL21"].dropna().items():
        for one in str(v).split(","):
            TGT[one] = (p_, s_)     # protein -> (패턴, 시그니처)
print(f"검증 대상 {len(TGT)}개")
for k, (p_, s_) in sorted(TGT.items(), key=lambda x: x[1][0]):
    print(f"  {p_}  {k:14s} {s_}")

# ---------- 서열 ----------
SEQ = {}
for s in STRAINS:
    nm, buf = None, []
    for l in open(FAA[s], errors="ignore"):
        if l.startswith(">"):
            if nm: SEQ[nm] = "".join(buf)
            nm, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if nm: SEQ[nm] = "".join(buf)

q = BASE/"input"/"seq"/"m4_bl21_nohit.faa"; q.parent.mkdir(parents=True, exist_ok=True)
q.write_text("".join(f">{k}\n{SEQ[k]}\n" for k in TGT if k in SEQ))
print(f"\n질의 {sum(1 for k in TGT if k in SEQ)} / {len(TGT)}개 (FASTA 에 있는 것만)")

# ---------- 1. MG1655 프로테옴 ----------
print("\n" + "=" * 96); print("### 1. MG1655 프로테옴에 직교체가 있는가"); print("=" * 96)
o = BASE/"result"/"search"/"m4_vs_MG1655.m8"; o.parent.mkdir(parents=True, exist_ok=True)
FMT = "query,target,fident,alnlen,qcov,tcov,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
if not (o.exists() and o.stat().st_size):
    sh(f'mmseqs easy-search "{q}" "{FAA["MG1655"]}" "{o}" "{_tmp}" '
       f'--format-output "{FMT}" -e 1e-3 -s 7.5 -a 1 --threads 16 -v 1')
BEST = {}
if o.exists() and o.stat().st_size:
    d = pd.read_csv(o, sep="\t", names=FMT.split(","))
    d = d.sort_values("bits", ascending=False).drop_duplicates("query")
    BEST = {r["query"]: r for _, r in d.iterrows()}
print(f"  직교체 있음 {len(BEST)} / {len(TGT)}")

# ---------- 2. Cys 쌍이 보존됐는가 ----------
def cys_kept(row, cys_pos):
    """질의의 Cys 위치가 정렬 상대에서도 C 인가. (보존수, 확인가능수, 상대위치)"""
    qa, ta = str(row["qaln"]), str(row["taln"])
    qi = int(row["qstart"]); ti = int(row["tstart"])
    keep, seen, at = 0, 0, []
    for cq, ct in zip(qa, ta):
        if cq != "-":
            if qi in cys_pos:
                seen += 1
                if ct == "C": keep += 1; at.append(f"{ti}C")
                else:         at.append(f"{ti}{ct}")
            qi += 1
        if ct != "-": ti += 1
    return keep, seen, ",".join(at)

rows = []
for k, (pat, sig_) in TGT.items():
    cys = [int(x[1:]) for x in sig_.split(",")]
    r = BEST.get(k)
    if r is None:
        rows.append({"protein": k, "패턴": pat, "sig": sig_, "MG직교체": "없음",
                     "fident": None, "qcov": None, "Cys보존": "", "판정": "게놈 확인 필요"})
        continue
    keep, seen, at = cys_kept(r, set(cys))
    ver = ("미검출 (직교체·Cys 둘 다 있음)" if seen and keep == seen
           else "Cys 달라짐" if seen
           else "정렬 구간 밖")
    rows.append({"protein": k, "패턴": pat, "sig": sig_,
                 "MG직교체": str(r["target"]), "fident": round(float(r["fident"]), 3),
                 "qcov": round(float(r["qcov"]), 3),
                 "Cys보존": f"{keep}/{seen} ({at})", "판정": ver})
V = pd.DataFrame(rows)
pd.set_option("display.max_colwidth", 44)
print("\n" + V.to_string(index=False))

# ---------- 3. 직교체가 없는 것만 게놈 확인 ----------
NOHIT = [r["protein"] for r in rows if r["MG직교체"] == "없음"]
print("\n" + "=" * 96); print(f"### 2. 게놈 확인 ({len(NOHIT)}개)"); print("=" * 96)
GEN = {}
if NOHIT:
    g = BASE/"input"/"external"/"MG1655_U00096.3.fna"
    if not (g.exists() and g.stat().st_size > 1_000_000):
        url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
               "?db=nuccore&id=U00096.3&rettype=fasta&retmode=text")
        try:
            with urllib.request.urlopen(url, timeout=300) as r_: g.write_bytes(r_.read())
            print(f"  게놈 받음 {g.stat().st_size/1e6:.1f} MB")
        except Exception as e:
            print(f"  ⚠ 게놈 내려받기 실패 {e}")
    if g.exists() and g.stat().st_size > 1_000_000:
        q2 = BASE/"input"/"seq"/"m4_nohit.faa"
        q2.write_text("".join(f">{k}\n{SEQ[k]}\n" for k in NOHIT if k in SEQ))
        o2 = BASE/"result"/"search"/"m4_vs_MG1655genome.m8"
        GF = "query,target,fident,alnlen,qlen,qstart,qend,evalue,bits"
        if not (o2.exists() and o2.stat().st_size):
            sh(f'mmseqs easy-search "{q2}" "{g}" "{o2}" "{_tmp}" --search-type 2 '
               f'--format-output "{GF}" -e 1e-3 -s 7.5 --threads 16 -v 1')
        if o2.exists() and o2.stat().st_size:
            d2 = pd.read_csv(o2, sep="\t", names=GF.split(","))
            for k_, gg in d2.groupby("query"):
                top = gg.sort_values("bits", ascending=False).iloc[0]
                cov = (float(top.qend) - float(top.qstart) + 1)/float(top.qlen)
                GEN[k_] = (f"게놈에 온전 (fid {top.fident:.2f} cov {cov:.2f})"
                           if top.fident >= 0.80 and cov >= 0.70 else
                           f"게놈에 조각 (fid {top.fident:.2f} cov {cov:.2f})"
                           if top.fident >= 0.80 else
                           f"먼 유사체만 (fid {top.fident:.2f})")
    for k in NOHIT:
        v = GEN.get(k, "게놈에도 없음")
        print(f"  {k:14s} {v}")
        V.loc[V.protein == k, "판정"] = ("진짜 부재" if v == "게놈에도 없음"
                                         else f"미검출 — {v}")
else:
    print("  전부 프로테옴에 직교체가 있다.")

V.to_csv(TBL/"motif_strain_specific_verified.csv", index=False, encoding="utf-8-sig")
print("\n" + "=" * 96); print("### 판정 요약"); print("=" * 96)
print(V.판정.value_counts().to_string())
real = V[V.판정 == "진짜 부재"]
print(f"\n★ 검증을 통과한 '진짜 MG1655 부재' {len(real)}개")
print(real[["protein", "패턴", "sig", "Cys보존"]].to_string(index=False) if len(real) else "  없음")
print(f"\n저장: {TBL/'motif_strain_specific_verified.csv'}")
print("\n  '미검출' 로 분류된 것은 균주 변별에 쓸 수 없다. 구조 DB 의 버전 차이나")
print("  질의 민감도가 만든 차이지 생물학이 아니다.")
```
