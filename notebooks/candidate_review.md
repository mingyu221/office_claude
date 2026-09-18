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
