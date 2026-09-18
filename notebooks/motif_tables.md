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

print("\n### 저장")
for f in ["motif_counts_by_strain.csv", "motif_metal_and_atp.csv", "motif_metal_ntp.csv"]:
    print(f"  {TBL/f}")
print("\n### 읽는 법")
print("  표3 의 분류는 주석 문자열에서 나온 것이지 측정이 아니다.")
print("  evidence 열의 단어를 보고 직접 걸러낼 것. ATP모티프 O 는 구조 쪽 근거가")
print("  하나 더 있다는 뜻이고, 이름이 빈 행은 아무 판정도 받지 않았다는 뜻이다.")
```
