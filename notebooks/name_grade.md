# Ni 배위 등급표에 이름 붙이기

`ni_site_grade.csv` 의 `desc` 열이 비어 있어 `prey` 가 접근번호뿐이다.
`MG1655-P06987` 이 무엇인지 표만 봐서는 알 수 없다.

이 노트북은 **셀 하나**로 이름을 붙인다. 네트워크를 쓰지 않는다.

- GenBank (`QJZ*` `AKE*` `AGE*` `AHZ*` `AAC*` …) → 프로테옴 FASTA 헤더
- UniProt (`P06987` …) → 구조 파일 제목, 없으면 **서열로 프로테옴에서 찾음**

원본은 건드리지 않고 `ni_site_grade_named.csv` 로 따로 저장한다.

---

## CELL N1 — 이름 붙이고 다시 정리

```python
# =============================================================================
# CELL N1 | ni_site_grade.csv 에 단백질 이름을 붙인다 (독립 실행)
# =============================================================================
import re, itertools
from pathlib import Path
import pandas as pd

BASE  = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS = Path("/mnt/af2results/mingyu")
TBL   = BASE/"result"/"table"
SRC   = TBL/"ni_site_grade.csv"
OUT   = TBL/"ni_site_grade_named.csv"

FAA = {"BL21":   TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"bl21_db_match_qjz.faa",
       "MG1655": TOOLS/"database"/"protein_list"/"mg1655_protein.faa",
       "Y19":    TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"y19_db_match.faa"}
STRUCT = {"BL21":   TOOLS/"database"/"bacteriaDB"/"structures_UP000503272",
          "Y19":    TOOLS/"database"/"bacteriaDB"/"structures_UP000034085",
          "MG1655": TOOLS/"database"/"mg1655_uni"}

# ---------------------------------------------------------------- 이름 사전
NAME, SEQ2ACC = {}, {}
print("=" * 100); print("### 이름 사전"); print("=" * 100)
for st, fp in FAA.items():
    if not fp.exists():
        print(f"  X  {st:8s} {fp}"); continue
    D_, acc, buf, n0 = {}, None, [], len(NAME)
    def flush():
        if acc and buf:
            s = "".join(buf); D_.setdefault(s, acc); D_.setdefault(s[:60], acc)
    for l in open(fp, errors="ignore"):
        if l.startswith(">"):
            flush()
            h = l[1:].rstrip(); acc = h.split()[0]; buf = []
            txt = re.sub(r"^\S+\s*", "", h)
            txt = re.split(r"\s*\[(?:Escherichia|Citrobacter|Carboxydothermus)", txt)[0]
            if txt.strip(): NAME[acc] = txt.strip()[:70]
        else:
            buf.append(l.strip())
    flush(); SEQ2ACC[st] = D_
    print(f"  O  {st:8s} 이름 +{len(NAME)-n0:5d}  서열 {len(D_):5d}")

AA3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E",
       "GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F",
       "PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M"}

def struct_title(p):
    try: txt = open(p, errors="ignore").read(20000)
    except Exception: return None
    for pat in (r"_struct\.title\s+(?:'([^']*)'|\"([^\"]*)\"|(\S.*))",
                r"_entity\.pdbx_description\s+(?:'([^']*)'|\"([^\"]*)\"|(\S.*))"):
        m = re.search(pat, txt)
        if m:
            s = (m.group(1) or m.group(2) or m.group(3) or "").strip().strip("'\"")
            if s and s not in ("?", ".", "None"): return s
    ttl = " ".join(l[10:].strip() for l in txt.splitlines() if l.startswith("TITLE"))
    if ttl:
        m = re.search(r"PREDICTION FOR\s+(.*)", ttl, re.I)
        return re.sub(r"\s*\([A-Z0-9]+\)\s*$", "", (m.group(1) if m else ttl).strip()) or None
    return None

def struct_seq(p):
    seq, seen = [], set()
    try:
        if Path(p).suffix.lower() in (".cif", ".mmcif"):
            hdr, on = [], False
            for l in open(p, errors="ignore"):
                if l.startswith("_atom_site."):
                    hdr.append(l.strip().split(".")[1]); on = True; continue
                if on and l[:4] == "ATOM":
                    c = {n: i for i, n in enumerate(hdr)}; f = l.split()
                    try:
                        if f[c["label_atom_id"]].strip('"') != "CA": continue
                        rs = f[c.get("auth_seq_id", c.get("label_seq_id"))]
                        rn = f[c["label_comp_id"]]
                    except Exception: continue
                    if rs in seen: continue
                    seen.add(rs); seq.append(AA3.get(rn, "X"))
                elif on and l.startswith("#") and seq: break
        else:
            for l in open(p, errors="ignore"):
                if l.startswith("ATOM") and l[12:16].strip() == "CA":
                    rs = l[22:27]
                    if rs in seen: continue
                    seen.add(rs); seq.append(AA3.get(l[17:20].strip(), "X"))
    except Exception: return ""
    return "".join(seq)

def seq_lookup(seq, prefer=None):
    order = ([prefer] if prefer in SEQ2ACC else []) + [k for k in SEQ2ACC if k != prefer]
    for st in order:
        d = SEQ2ACC[st]
        a = d.get(seq) or d.get(seq[:60])
        if a: return a, st
    return None, None

def find_struct(acc, prefer=None):
    order = ([prefer] if prefer in STRUCT else []) + [k for k in STRUCT if k != prefer]
    for st in order:
        d = STRUCT.get(st)
        if not d or not Path(d).is_dir(): continue
        for pat in (f"AF-{acc}-F1-model_v*.cif", f"AF-{acc}-F1-model_v*.pdb",
                    f"AF-{acc}-*.cif", f"AF-{acc}-*.pdb",
                    f"{acc}.cif", f"{acc}.pdb", f"*{acc}*.cif", f"*{acc}*.pdb"):
            for f in Path(d).glob(pat): return f, st
    return None, None

CACHE = {}
def label(acc, strain=None):
    key = (acc, strain)
    if key in CACHE: return CACHE[key]
    n = NAME.get(acc)                                   # GenBank 는 바로
    if not n:
        f, st = find_struct(acc, prefer=strain)         # UniProt 는 구조에서
        if f is not None:
            n = struct_title(f)
            if not n:
                gb, _ = seq_lookup(struct_seq(f), prefer=st or strain)
                if gb: n = f"{NAME.get(gb, '')} [{gb}]".strip()
    CACHE[key] = n or ""
    return CACHE[key]

# ---------------------------------------------------------------- 표
if not SRC.exists():
    raise SystemExit(f"{SRC} 없음 — CELL 28a-8b 를 먼저 돌릴 것")
G = pd.read_csv(SRC)
print(f"\n  {SRC.name}: {len(G)}행, 열 {list(G.columns)}")

def split_prey(v):
    s = str(v).strip()
    for st in ("Y19", "MG1655", "BL21"):
        if s.startswith(st + "-"): return st, s[len(st)+1:]
    return "BL21", s                                    # 접두사 없으면 BL21

G["strain"] = [split_prey(v)[0] for v in G.prey]
G["acc"]    = [split_prey(v)[1] for v in G.prey]
G["name"]   = [label(a, s) for a, s in zip(G.acc, G.strain)]

COLS = ["grade", "strain", "acc", "name", "Cys", "His", "Asp/Glu", "n_donor",
        "min_dist", "donors", "folddisco_metal", "set", "why"]
G = G[[c for c in COLS if c in G.columns] + [c for c in G.columns if c not in COLS]]
G = G.sort_values(["grade", "min_dist"])
G.to_csv(OUT, index=False, encoding="utf-8-sig")

pd.set_option("display.max_rows", None); pd.set_option("display.width", 220)
pd.set_option("display.max_colwidth", 46)
print("\n" + "=" * 100); print("### 등급표"); print("=" * 100)
print(G[["grade", "strain", "acc", "name", "Cys", "His", "n_donor",
         "min_dist", "folddisco_metal"]].to_string(index=False))

print("\n" + "=" * 100); print("### 등급별 요약"); print("=" * 100)
for g_, sub in G.groupby("grade"):
    print(f"\n  {g_}  {len(sub)}행   거리 {sub.min_dist.min():.2f}–{sub.min_dist.max():.2f} Å"
          f"   folddisco 검출 {int(sub.folddisco_metal.astype(str).eq('True').sum())}/{len(sub)}")
    for _, r in sub.iterrows():
        print(f"      {r.strain:7s} {r.acc:16s} {str(r['name'])[:44]:44s} "
              f"{r.min_dist:.2f} Å  {r.donors}")

print("\n" + "=" * 100); print("### 배위 잔기 서명으로 묶으면"); print("=" * 100)
G["sig"] = [" ".join(sorted(str(d).split())) for d in G.donors]
for sig, sub in G.groupby("sig"):
    if len(sub) < 2: continue
    print(f"\n  {sig}   ({len(sub)}개 — 같은 자리)")
    for _, r in sub.iterrows():
        print(f"      {r.grade} {r.strain:7s} {r.acc:16s} {str(r['name'])[:40]:40s} {r.min_dist:.2f} Å")
solo = G.groupby("sig").filter(lambda x: len(x) == 1)
print(f"\n  단독 서명 {len(solo)}개:")
for _, r in solo.iterrows():
    print(f"      {r.grade} {r.strain:7s} {r.acc:16s} {str(r['name'])[:40]:40s} {r.min_dist:.2f} Å")

miss = sorted({a for a, n in zip(G.acc, G["name"]) if not n})
print("\n" + "=" * 100); print("### 감사"); print("=" * 100)
print(f"  이름 못 찾은 것 {len(miss)}개" + (f": {miss}" if miss else ""))
for a in miss[:5]:
    base = ("https://www.ncbi.nlm.nih.gov/protein/" if re.match(r"^[A-Z]{3}\d+\.\d+$", a)
            else "https://www.uniprot.org/uniprotkb/")
    print(f"    {base}{a}")
print(f"\n저장: {OUT}   (원본 {SRC.name} 은 건드리지 않았다)")
```
