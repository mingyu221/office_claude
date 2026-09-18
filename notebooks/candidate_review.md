# 후보 검토용 통합 엑셀 — 독립 노트북

지금까지의 분석은 **엑셀이 아니라 `result/table/*.csv` 에 흩어져 있다.**
(업로드하신 `structure_accessions.xlsx` 는 구조↔Genbank 매핑 파일이고, 분석 결과가 아니다.)

이 노트북은 그 CSV 들을 **하나의 워크북**으로 모은다.
상위 5개만이 아니라 **가능성 있는 전체**를 시트별로 놓고 직접 보실 수 있게 하는 것이 목적이다.

산출물: `result/table/ChCODH2_candidates_review.xlsx`

---

## CELL R1 — 어떤 파일이 어디에 있고 무엇을 담는가

```python
# =============================================================================
# CELL R1 | 후보 판단에 쓰인 원본을 한눈에 놓는다
#   각 줄이 '무엇을 근거로 등급을 매겼는가' 에 대응한다.
# =============================================================================
from pathlib import Path
import pandas as pd

BASE = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TBL  = BASE/"result"/"table"

# (파일, 시트이름, 이 표가 답하는 질문)
SRC = [
    ("FINAL_candidates.csv",              "00_등급표",      "통합 등급 — 제외 없이 전부"),
    ("ALL_motif_candidates_full.csv",     "01_모티프전체",  "금속 모티프 80개 + 모든 증거 열"),
    ("motif_metal_3strain.csv",           "02_metal_3균주", "metal_rmsd · pLDDT · 균주 판정"),
    ("ni_site_grade.csv",                 "03_Ni배위등급",  "Ni 도너 등급 A~D 와 도너 잔기"),
    ("boltz_ni_placement.csv",            "04_Ni좌표",      "Boltz 가 놓은 Ni 의 배위 거리"),
    ("folddisco_x_foldseek.csv",          "05_모티프x폴드", "자리(Folddisco) × 과(Foldseek) 교차"),
    ("g3e_family_sweep.csv",              "06_G3E과",       "금속 샤페론 과 구성원 + tier"),
    ("hypA_fold_sweep.csv",               "07_HypA폴드",    "HypA 폴드를 가진 것 전수"),
    ("cys4_bl21_strain_specificity.csv",  "08_균주특이",    "프로테옴 수준 BL21 특이성"),
    ("cys4_genome_check.csv",             "09_게놈확인",    "게놈까지 봤을 때의 판정"),
    ("motif_counts_by_strain.csv",        "10_모티프계수",  "균주별 metal/ATP 매칭 수"),
    ("motif_metal_and_atp.csv",           "11_metal_ATP",  "두 모티프를 동시에 가진 것"),
    ("motif_metal_ntp.csv",               "12_NTP소모",     "metal 중 NTP 소모 주석"),
    ("trackA_RF2PPI_ranked.csv",          "13_TrackA",      "RF2-PPI 단량체 순위"),
    ("boltz_replicate_summary.csv",       "14_TrackB",      "Boltz ipTM 5표본 + 대조군"),
    ("ortholog_reproducibility.csv",      "15_잡음바닥",    "직교체 재현성 — 차이의 유의성 기준"),
    ("coo_operon_homologs.csv",           "16_coo오페론",   "Y19 오페론 상동체"),
    ("folddisco_sensitivity_sweep.csv",   "17_민감도",      "-d/-a 스윕 (민감도 노트북 산출)"),
    ("trackA_dimer_mode_assign.csv",      "18_이량체위상",  "이량체 재실행의 쌍별 위상"),
    ("trackA_dimer_oversize.csv",         "19_이량체제외",  "L 상한 초과로 점수 없는 쌍"),
]
print("=" * 100); print("### 원본 (있는 것만 엑셀에 들어간다)"); print("=" * 100)
HAVE = []
for f, sheet, what in SRC:
    p = TBL/f
    if p.exists():
        try: n = len(pd.read_csv(p))
        except Exception: n = -1
        HAVE.append((p, sheet, what, n))
        print(f"  O  {sheet:16s} {n:6d}행  {f:36s} {what}")
    else:
        print(f"  .  {sheet:16s}    —    {f:36s} {what}")
print(f"\n  {len(HAVE)} / {len(SRC)} 개가 있다.")
print("\n  ※ '.' 인 것은 아직 안 돌린 셀의 산출물이다. 나중에 이 셀을 다시 돌리면")
print("    자동으로 합쳐진다.")
```

---

## CELL R2 — 하나의 엑셀로 합친다

```python
# =============================================================================
# CELL R2 | 워크북 작성
#   첫 시트에 '무엇을 보고 있는지' 를 적어 둔다. 시트를 열었을 때 그 표가
#   어떤 질문에 답하는지 모르면 검토가 안 된다.
# =============================================================================
import re
XLSX = TBL/"ChCODH2_candidates_review.xlsx"
try:
    import xlsxwriter; ENGINE = "xlsxwriter"
except ImportError:
    ENGINE = "openpyxl"
print(f"엔진 {ENGINE} -> {XLSX}")

idx = pd.DataFrame([{"시트": s, "행": n, "원본": p.name, "이 표가 답하는 것": w}
                    for p, s, w, n in HAVE])
with pd.ExcelWriter(XLSX, engine=ENGINE) as xw:
    idx.to_excel(xw, sheet_name="목차", index=False)
    for p, sheet, what, n in HAVE:
        try:
            d = pd.read_csv(p)
        except Exception as e:
            print(f"  ⚠ {p.name} 읽기 실패 {e}"); continue
        # 엑셀 시트 이름 제한: 31자, \ / * ? [ ] : 불가
        nm = re.sub(r"[\\/\*\?\[\]:]", "_", sheet)[:31]
        d.to_excel(xw, sheet_name=nm, index=False)
        print(f"  {nm:18s} {len(d):6d}행 x {d.shape[1]}열")
print(f"\n저장: {XLSX}")

# ---------- 후보를 한 장에 모은다 ----------
# 등급표가 있으면 그걸 쓰고, 없으면 모티프 전체 + Ni 배위 등급을 붙여 만든다.
print("\n" + "=" * 100); print("### 후보 한 장 보기"); print("=" * 100)
cand = None
if (TBL/"FINAL_candidates.csv").exists():
    cand = pd.read_csv(TBL/"FINAL_candidates.csv")
    print(f"  FINAL_candidates.csv {len(cand)}행 — 제외 없이 전부 들어 있는 등급표")
elif (TBL/"ALL_motif_candidates_full.csv").exists():
    cand = pd.read_csv(TBL/"ALL_motif_candidates_full.csv")
    print(f"  ALL_motif_candidates_full.csv {len(cand)}행")
if cand is not None:
    pd.set_option("display.max_rows", None); pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 40)
    sortcol = next((c for c in ["score", "tier", "grade", "metal_rmsd"] if c in cand.columns), None)
    if sortcol: cand = cand.sort_values(sortcol, ascending=(sortcol == "metal_rmsd"))
    show = [c for c in ["protein", "strain", "tier", "score", "grade", "donors",
                        "metal_rmsd", "plddt", "vs_MG1655", "verdict", "reason", "desc"]
            if c in cand.columns]
    print(cand[show].to_string(index=False))
    print(f"\n  등급 분포: "
          + (cand[sortcol].value_counts().to_dict().__str__() if sortcol else "—"))
else:
    print("  등급표가 아직 없다 — CELL 39 를 돌리면 생긴다.")
    print("  그전까지는 01_모티프전체 · 03_Ni배위등급 시트를 직접 보시면 된다.")

print("\n### 읽는 순서 제안")
print("  1) 목차 시트에서 각 표가 무엇에 답하는지 본다")
print("  2) 00_등급표 — 제외 없이 전부. 점수가 낮다고 버린 게 아니라 낮게 둔 것이다")
print("  3) 03_Ni배위등급 — A/B 가 자리의 질. 이게 '금속을 다룰 개연성'의 1차 근거")
print("  4) 09_게놈확인 — '균주 특이' 를 액면 그대로 믿지 않기 위한 표")
print("  5) 17_민감도 — 이 결과들이 어느 민감도에서 나온 것인지. 0 의 해석이 여기 달렸다")
```

---

## CELL R3 — 후보 전체 × 4축 통과 여부

```python
# =============================================================================
# CELL R3 | "네 가지를 모두 만족" 을 검증 가능한 표로 만든다
#
#   네 축 (각각 어느 파일에서 오는지 같이 남긴다)
#     1. Cys4 재질의 히트   folddisco/cys4/cys4*_{템플릿}_{균주}.tsv
#          실제로 Ni 을 잡은 Cys4 기하를 질의로 썼을 때 잡혔는가
#     2. HypA 폴드          hypA_fold_sweep.csv  (Foldseek tm >= 0.5)
#     3. Ni 배위 grade A    ni_site_grade.csv    (Cys4 가 2.6 A 안)
#     4. MG1655 부재        cys4_bl21_strain_specificity.csv + cys4_genome_check.csv
#          프로테옴에도 게놈에도 없는가
#
#   ★ 축마다 신뢰도가 다르다. 표에 O 가 찍혔다고 같은 무게가 아니다.
#     1 번: 질의의 직교체 민감도가 0 으로 확인됐다 (Ch CooC1 질의가 Y19 CooC 를
#           못 찾는다). 즉 '안 잡힘' 이 '없음' 을 뜻하지 않는다. O 는 근거가 되지만
#           . 은 근거가 안 된다.
#     3 번: Boltz 예측 구조에서 잰 것이지 실측이 아니다.
#     4 번: 구조 DB 버전 차이(BL21 v6 / MG1655 v4)로 가짜 '특이' 가 나온 전례가
#           있다. 그래서 게놈까지 봐야 하고, 이 표는 둘을 따로 보여준다.
#     2 번이 셋 중 가장 단단하다 — Foldseek 은 대조군을 통과한 유일한 도구다.
# =============================================================================
import re
FDD = BASE/"result"/"folddisco"

def crosswalk():
    x = TBL/"id_crosswalk_struct_to_genbank.csv"
    if not x.exists(): return {}
    d = pd.read_csv(x).dropna(subset=["protein"])
    return {Path(str(k)).stem: str(v).split(",")[0]
            for k, v in d.set_index("tid")["protein"].to_dict().items()}
XW = crosswalk()
def pid_of(v):
    s = str(v)
    return XW.get(Path(s).stem, s.split("-")[-1])

# ---------- 축 1 : Cys4 재질의 ----------
AX1, ax1_src = set(), []
cy = FDD/"cys4"
if cy.is_dir():
    cov4 = sorted(cy.glob("cys4cov4_*.tsv"))
    files = cov4 if cov4 else sorted(cy.glob("cys4_*.tsv"))
    ax1_src = [f.name for f in files]
    for f in files:
        if not f.stat().st_size: continue
        d = pd.read_csv(f, sep="\t"); d.columns = [c.strip().lstrip("#") for c in d.columns]
        AX1 |= {pid_of(t) for t in d.tid}
    print(f"축1 Cys4 재질의: {len(AX1)}개  "
          f"({'covered-node 4' if cov4 else '부분일치 포함(구버전)'}, 파일 {len(files)})")
else:
    print("축1: cys4 디렉터리 없음 — CELL 29c 미실행")

# ---------- 축 2 : HypA 폴드 ----------
AX2, AX2_TM = set(), {}
p = TBL/"hypA_fold_sweep.csv"
if p.exists():
    d = pd.read_csv(p)
    kc = next((c for c in ["protein", "member", "target", "prey"] if c in d.columns), d.columns[0])
    tc = next((c for c in ["tm", "alntmscore", "top_tm"] if c in d.columns), None)
    for _, r in d.iterrows():
        k = pid_of(r[kc]); AX2.add(k)
        if tc: AX2_TM[k] = float(r[tc])
    print(f"축2 HypA 폴드: {len(AX2)}개  (원본 {p.name}, tm 열 {tc})")
else:
    print("축2: hypA_fold_sweep.csv 없음 — CELL 50 STEP5 미실행")

# ---------- 축 3 : Ni 배위 등급 ----------
GRADE, DONOR = {}, {}
p = TBL/"ni_site_grade.csv"
if p.exists():
    d = pd.read_csv(p)
    for _, r in d.iterrows():
        k = str(r["prey"]).split("-")[-1]
        GRADE[k] = str(r.get("grade", "")); DONOR[k] = str(r.get("donors", ""))[:44]
    print(f"축3 Ni 배위: {len(GRADE)}개 채점  "
          f"(A {sum(v=='A' for v in GRADE.values())}, B {sum(v=='B' for v in GRADE.values())})")
else:
    print("축3: ni_site_grade.csv 없음 — CELL 28a-8b 미실행")

# ---------- 축 4 : MG1655 부재 (프로테옴 / 게놈) ----------
PROT_ABS, GENOME = {}, {}
p = TBL/"cys4_bl21_strain_specificity.csv"
if p.exists():
    d = pd.read_csv(p)
    for _, r in d.iterrows():
        PROT_ABS[str(r["protein"])] = str(r.get("vs_MG1655", ""))
    print(f"축4a 프로테옴 판정: {len(PROT_ABS)}개")
p = TBL/"cys4_genome_check.csv"
if p.exists():
    d = pd.read_csv(p)
    for _, r in d.iterrows():
        GENOME[str(r["protein"])] = str(r.get("verdict", ""))
    print(f"축4b 게놈 판정: {len(GENOME)}개")
p = TBL/"g3e_family_sweep.csv"
if p.exists():
    d = pd.read_csv(p)
    if "genome_MG1655" in d.columns:
        for _, r in d.iterrows():
            v = str(r.get("genome_MG1655", ""))
            if v and v != "nan": GENOME.setdefault(str(r["member"]), v)

# ---------- 대상 모으기 ----------
UNIV, META = set(), {}
def add(keys, **kw):
    for k in keys:
        k = str(k)
        if not k or k == "nan": continue
        UNIV.add(k); META.setdefault(k, {}).update({a: b for a, b in kw.items() if b == b})
p = TBL/"motif_metal_3strain.csv"
if p.exists():
    d = pd.read_csv(p)
    for _, r in d.iterrows():
        k = str(r["key"]).split("-")[-1]
        add([k], strain=r.get("strain"), metal_rmsd=r.get("rmsd"), plddt=r.get("plddt"))
add(AX1); add(AX2); add(GRADE); add(PROT_ABS); add(GENOME)
p = TBL/"g3e_family_sweep.csv"
if p.exists():
    d = pd.read_csv(p)
    for _, r in d.iterrows():
        add([r["member"]], strain=r.get("strain"), g3e_tier=r.get("tier"))
print(f"\n대상 {len(UNIV)}개 (어느 축에든 한 번이라도 등장한 단백질 전부)")

# ---------- 표 ----------
rows = []
for k in sorted(UNIV):
    m = META.get(k, {})
    g = GRADE.get(k, "")
    pa = PROT_ABS.get(k, ""); gn = GENOME.get(k, "")
    a1 = "O" if k in AX1 else "."
    a2 = "O" if k in AX2 else "."
    a3 = "O" if g == "A" else ("~" if g == "B" else ".")
    # 4 번은 프로테옴·게놈 둘 다 '없음' 이어야 O. 하나만이면 ~, 판정이 없으면 ?
    prot_no = ("specific" in pa.lower()) or ("없" in pa)
    gen_no  = ("게놈에도 없음" in gn) or ("사실상 없음" in gn) or ("없음" == gn)
    if not pa and not gn: a4 = "?"
    elif prot_no and gen_no: a4 = "O"
    elif prot_no or gen_no:  a4 = "~"
    else:                    a4 = "."
    rows.append({"protein": k, "strain": m.get("strain", ""),
                 "1_Cys4": a1, "2_HypA": a2, "3_Ni등급": a3, "4_MG없음": a4,
                 "통과": sum(x == "O" for x in (a1, a2, a3, a4)),
                 "grade": g, "HypA_tm": round(AX2_TM[k], 3) if k in AX2_TM else None,
                 "metal_rmsd": m.get("metal_rmsd"), "plddt": m.get("plddt"),
                 "g3e_tier": m.get("g3e_tier", ""),
                 "프로테옴": pa[:22], "게놈": gn[:34], "donors": DONOR.get(k, "")})
T = pd.DataFrame(rows).sort_values(["통과", "3_Ni등급", "metal_rmsd"],
                                   ascending=[False, True, True])
pd.set_option("display.max_rows", None); pd.set_option("display.width", 260)
pd.set_option("display.max_colwidth", 40)
print("\n" + "=" * 110); print("### 4축 통과 표 (O 통과 / ~ 부분 / . 아님 / ? 판정없음)"); print("=" * 110)
print(T.to_string(index=False))
T.to_csv(TBL/"candidate_4axis_matrix.csv", index=False, encoding="utf-8-sig")
print(f"\n저장: {TBL/'candidate_4axis_matrix.csv'}")

print("\n" + "=" * 110); print("### 요약"); print("=" * 110)
print(T.통과.value_counts().sort_index(ascending=False).to_string())
for n in (4, 3):
    sub = T[T.통과 == n]
    print(f"\n  {n}축 통과 {len(sub)}개")
    if len(sub):
        print(sub[["protein", "strain", "1_Cys4", "2_HypA", "3_Ni등급", "4_MG없음",
                   "grade", "HypA_tm", "metal_rmsd", "게놈"]].to_string(index=False))
print("\n  ※ '.' 을 '아니다' 로 읽지 말 것. 특히 1번 축은 질의의 직교체 민감도가")
print("    0 으로 확인된 상태라, 안 잡힌 것이 없다는 뜻이 아니다 (민감도 노트북 참조).")
print("    4번 축의 '?' 는 그 단백질이 균주특이 판정 대상에 아예 안 들어갔다는 뜻이다.")
```
