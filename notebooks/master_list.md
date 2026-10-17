# 마스터 리스트 — 후보 전체를 한 장으로

"전체 리스트 CSV 가 뭐냐" 에 대한 답. **지금까지 그런 파일이 없었다.**
단계별 CSV 가 20여 개 있었고, 각각 다른 질문의 답이라 한눈에 못 봤다.

이 노트북은 그것들을 **단백질 하나당 한 행**으로 합친다.
출력은 둘이고, 둘이 같은 내용이다.

| 파일 | 용도 |
|---|---|
| `result/table/MASTER_candidate_list.csv` | ★ **이것이 전체 리스트다** |
| `result/table/MASTER_candidate_list.xlsx` | 같은 내용 + 열 설명 시트 + 필터 걸린 시트 |

셀 하나만 돌리면 된다. 없는 원본은 알아서 건너뛰고 무엇이 빠졌는지 찍는다.

---

## CELL M1 — 전부 합친다

```python
# =============================================================================
# CELL M1 | 모든 단계의 결과를 단백질 단위로 조인한다
#   원칙 셋.
#     1. 있는 것만 쓴다. 없는 원본은 건너뛰고 목록에 남긴다
#     2. 한 열은 한 출처에서만 온다. 값이 겹치면 열을 따로 둔다
#     3. 판정을 덮어쓰지 않는다. 옛 파이프라인과 새 파이프라인을 나란히 둔다
# =============================================================================
import re
from pathlib import Path
import pandas as pd

BASE  = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS = Path("/mnt/af2results/mingyu")
TBL   = BASE/"result"/"table"
FAA_BL21 = TOOLS/"database"/"bacteriaDB"/"inhouseDB"/"bl21_db_match_qjz.faa"

def rd(name, **kw):
    p = TBL/name
    if not (p.exists() and p.stat().st_size):
        MISSING.append(name); return None
    try:
        return pd.read_csv(p, **kw)
    except Exception as e:
        MISSING.append(f"{name} (읽기 실패: {e})"); return None

def key_of(v):
    """'Y19-AKE60400.1' → 'AKE60400.1', 경로 → stem"""
    s = str(v).strip()
    if "/" in s: s = Path(s).stem
    return s.split("-")[-1] if s.split("-")[0] in ("BL21", "Y19", "MG1655") else s

MISSING, SRC = [], []
M = {}          # protein -> {열: 값}

def put(pid, **kw):
    if not pid or str(pid) == "nan": return
    M.setdefault(str(pid), {}).update({k: v for k, v in kw.items() if pd.notna(v)})

# ---------------------------------------------------------------- 1. 서열 모티프
d = rd("seqmotif_all_strains.csv")
if d is not None:
    SRC.append(("seqmotif_all_strains.csv", "서열 Cys 모티프 스캔", len(d)))
    for _, r in d[d.strain == "BL21"].iterrows():
        put(r.protein, len=r.get("len"), desc=r.get("desc"), n_cys=r.get("n_cys"),
            n_CXC=r.get("n_CXC"), n_CXXC=r.get("n_CXXC"),
            cys4_win=r.get("cys4_win"), cys4_cand=r.get("cys4_cand"))

# ---------------------------------------------------------------- 2. 균주 판정
d = rd("seqmotif_strain_verdict.csv")
if d is not None:
    SRC.append(("seqmotif_strain_verdict.csv", "MG1655 대비 (모티프 Cys 보존까지)", len(d)))
    for _, r in d.iterrows():
        put(r.protein, mg_status=r.get("mg_status"), mg_hit=r.get("mg_hit"),
            mg_fident=r.get("fident"), cys_kept=r.get("cys_kept"),
            cys_total=r.get("cys_total"), bl21_only=r.get("bl21_only"),
            mobile=r.get("mobile"))

# ---------------------------------------------------------------- 3. 폴드
for fn, col in [("hypA_fold_sweep_exact.csv", "HypA"), ("hypA_fold_sweep.csv", "HypA"),
                ("g3e_family_sweep_exact.csv", "G3E"), ("g3e_family_sweep.csv", "G3E")]:
    d = rd(fn)
    if d is None: continue
    SRC.append((fn, f"{col} 과 스윕", len(d)))
    kc = next((c for c in ["protein", "member", "target"] if c in d.columns), d.columns[0])
    tc = next((c for c in ["tm", "alntmscore", "top_tm"] if c in d.columns), None)
    exact = fn.endswith("_exact.csv")
    for _, r in d.iterrows():
        k = key_of(r[kc])
        v = float(r[tc]) if tc and pd.notna(r.get(tc)) else None
        cur = M.get(k, {}).get(f"tm_{col}")
        if v is not None and (cur is None or exact):
            put(k, **{f"tm_{col}": round(v, 3),
                      f"tm_{col}_exact": int(exact),
                      f"tier_{col}": r.get("tier")})
    # _exact 를 읽었으면 구버전은 건너뛴다
    if exact: 
        for nxt in [fn.replace("_exact", "")]:
            if (TBL/nxt).exists(): MISSING.append(f"{nxt} (exact 판 사용으로 생략)")

d = rd("cooc_fold_afframe.csv")
if d is not None:
    SRC.append(("cooc_fold_afframe.csv", "CooC 폴드 (F6)", len(d)))
    for _, r in d[(d.strain == "BL21") & (d.frame == "예측")].iterrows():
        put(key_of(r.protein), tm_CooC=r.get("tm"), tier_CooC=r.get("tier"))

d = rd("cooc_fold_mg1655.csv")
if d is not None:
    SRC.append(("cooc_fold_mg1655.csv", "CooC 과 × MG1655 (F7)", len(d)))
    for _, r in d.iterrows():
        put(key_of(r.protein), cooc_vs_MG=r.get("MG1655"))

# ---------------------------------------------------------------- 4. Ni 배위 등급
d = rd("ni_site_grade.csv")
if d is not None:
    SRC.append(("ni_site_grade.csv", "Ni 배위 등급 (구 Track B, ChCODH2 와 함께)", len(d)))
    for _, r in d.iterrows():
        put(key_of(r.prey), ni_grade_old=r.get("grade"), ni_Cys_old=r.get("Cys"),
            ni_His_old=r.get("His"), ni_dist_old=r.get("min_dist"),
            ni_donors_old=r.get("donors"))

d = rd("seqmotif_ni_grade.csv")
if d is not None:
    SRC.append(("seqmotif_ni_grade.csv", "Ni 배위 등급 (신 파이프라인, 단량체+Ni)", len(d)))
    for _, r in d.iterrows():
        put(key_of(r.protein), ni_grade_mono=r.get("grade"), ni_Cys_mono=r.get("Cys"),
            ni_dist_mono=r.get("min_dist"), ni_donors_mono=r.get("donors"))

for fn, pre in [("seqmotif_dimer_grade.csv", "dimer"), ("cooclike_dimer_grade.csv", "cooclike")]:
    d = rd(fn)
    if d is None: continue
    SRC.append((fn, f"Ni 배위 등급 (동형이량체+Ni, {pre})", len(d)))
    for _, r in d.iterrows():
        put(key_of(r.protein), **{f"ni_grade_{pre}": r.get("grade"),
                                  f"ni_Cys_{pre}": r.get("Cys"),
                                  f"ni_계면_{pre}": r.get("계면")})

# ---------------------------------------------------------------- 5. Folddisco (참고용)
d = rd("motif_metal_3strain.csv")
if d is not None:
    SRC.append(("motif_metal_3strain.csv", "Folddisco 금속 모티프 (※ 비검출은 근거 아님)", len(d)))
    for _, r in d[d.get("strain", "BL21") == "BL21"].iterrows():
        put(key_of(r.get("key")), fd_metal=1, fd_rmsd=r.get("rmsd"), plddt_struct=r.get("plddt"))

# ---------------------------------------------------------------- 6. 점수
d = rd("seqmotif_ranked.csv")
if d is not None:
    SRC.append(("seqmotif_ranked.csv", "서열 파이프라인 점수", len(d)))
    for _, r in d.iterrows():
        put(r.protein, score=r.get("score"), why=r.get("why"),
            plddt_motif=r.get("plddt_motif"), plddt_mean=r.get("plddt_mean"))

d = rd("cooclike_ranked.csv")
if d is not None:
    SRC.append(("cooclike_ranked.csv", "CooC 형 트랙 점수", len(d)))
    for _, r in d.iterrows():
        put(r.protein, cooclike_score=r.get("dscore"), cooclike_why=r.get("dwhy"))

# ---------------------------------------------------------------- 7. Track A / B 원점수
d = rd("trackA_RF2PPI_ranked.csv")
if d is not None:
    SRC.append(("trackA_RF2PPI_ranked.csv", "RF2-PPI (※ 대조군 0.266 실패)", len(d)))
    kc = next((c for c in ["prey", "protein"] if c in d.columns), None)
    vc = next((c for c in ["score", "mean", "max", "prob"] if c in d.columns), None)
    if kc and vc:
        for _, r in d.iterrows(): put(key_of(r[kc]), rf2ppi=r.get(vc))

d = rd("trackB_boltz2_ranked.csv")
if d is not None:
    SRC.append(("trackB_boltz2_ranked.csv", "Boltz ipTM (※ 눈금 미검증)", len(d)))
    kc = next((c for c in ["prey", "protein"] if c in d.columns), None)
    vc = next((c for c in ["iptm", "ipTM", "score"] if c in d.columns), None)
    if kc and vc:
        for _, r in d.iterrows(): put(key_of(r[kc]), boltz_iptm=r.get(vc))

# ---------------------------------------------------------------- 8. 표 만들기
T = pd.DataFrame.from_dict(M, orient="index").rename_axis("protein").reset_index()
T = T[T.protein.astype(str).str.startswith("QJZ")]          # BL21 만

# 설명이 비었으면 FASTA 헤더에서 채운다
if FAA_BL21.exists() and ("desc" not in T or T.desc.isna().any()):
    H, k = {}, None
    for l in open(FAA_BL21, errors="ignore"):
        if l.startswith(">"):
            t2 = l[1:].rstrip(); k = t2.split()[0]; H[k] = t2[len(k):].strip()[:70]
    T["desc"] = [d if pd.notna(d) else H.get(p) for p, d in zip(T.protein, T.get("desc"))]

def n_evi(r):
    """독립 근거 개수 — 서로 다른 축에서 온 것만 센다."""
    n = 0
    if r.get("cys4_cand") == 1 or (r.get("n_CXC") or 0) or (r.get("n_CXXC") or 0): n += 1
    if str(r.get("ni_grade_old")) in ("A", "B") or str(r.get("ni_grade_mono")) in ("A", "B") \
       or str(r.get("ni_grade_dimer", "")).startswith("A"): n += 1
    if pd.notna(r.get("tm_HypA")) or pd.notna(r.get("tm_G3E")) or pd.notna(r.get("tm_CooC")): n += 1
    if r.get("bl21_only") == 1 and r.get("mobile") != 1: n += 1
    return n

T["근거축수"] = T.apply(n_evi, axis=1)
ORDER = ["protein", "desc", "len", "근거축수", "score",
         "cys4_win", "cys4_cand", "n_cys", "n_CXC", "n_CXXC",
         "mg_status", "bl21_only", "mobile", "mg_hit", "mg_fident", "cys_kept", "cys_total",
         "tm_HypA", "tm_G3E", "tm_CooC", "tier_CooC", "cooc_vs_MG",
         "ni_grade_old", "ni_Cys_old", "ni_dist_old",
         "ni_grade_mono", "ni_Cys_mono", "ni_dist_mono",
         "ni_grade_dimer", "ni_계면_dimer",
         "plddt_motif", "plddt_mean", "fd_metal", "fd_rmsd",
         "rf2ppi", "boltz_iptm", "cooclike_score", "why"]
T = T[[c for c in ORDER if c in T.columns] +
      [c for c in T.columns if c not in ORDER]]
T = T.sort_values(["근거축수", "score"], ascending=False, na_position="last")

OUT = TBL/"MASTER_candidate_list.csv"
T.to_csv(OUT, index=False, encoding="utf-8-sig")

print("=" * 100); print("### 합친 원본"); print("=" * 100)
for f, what, n in SRC: print(f"  O  {f:36s} {n:6d}행   {what}")
if MISSING:
    print("\n  못 읽은 것 (해당 셀을 안 돌렸거나 이름이 다름):")
    for f in MISSING: print(f"  X  {f}")

print("\n" + "=" * 100); print("### 마스터 리스트"); print("=" * 100)
print(f"  단백질 {len(T)}개 × 열 {len(T.columns)}개")
print(f"  ★ {OUT}")
print("\n  근거축수 분포 (서로 다른 축 몇 개에서 걸렸나):")
print(T["근거축수"].value_counts().sort_index(ascending=False).to_string())

SHOW = [c for c in ["protein", "근거축수", "score", "cys4_win", "mg_status", "mobile",
                    "tm_HypA", "tm_CooC", "ni_grade_old", "ni_grade_mono",
                    "plddt_motif", "desc"] if c in T.columns]
print("\n  상위 25")
print(T.head(25)[SHOW].to_string(index=False))

# ---------------------------------------------------------------- 9. 엑셀
try:
    XL = TBL/"MASTER_candidate_list.xlsx"
    DICT = pd.DataFrame([
        ("protein",      "BL21 GenBank ID"),
        ("근거축수",       "서로 다른 축(모티프/Ni배위/폴드/균주) 중 몇 개에서 걸렸나 0~4"),
        ("score",        "서열 파이프라인 점수 — why 열에 내역"),
        ("cys4_win",     "40잔기 창 안의 Cys 최대 개수"),
        ("mg_status",    "동일 / 모티프 보존 / 모티프 달라짐 / 게놈에만 있음 / 프로테옴 없음"),
        ("bl21_only",    "MG1655 프로테옴·게놈 둘 다 없으면 1"),
        ("mobile",       "프로파지·삽입 구간이면 1 — λDE3 는 균주 제작의 산물이다"),
        ("cys_kept/total","정렬을 따라가 모티프 Cys 가 상대에도 Cys 인 개수 / 검사 수"),
        ("tm_HypA/G3E/CooC","해당 과 스윕에서의 TM-score. 0.5 같은 폴드, 0.9 사실상 동일"),
        ("ni_grade_old", "구 Track B (ChCODH2 와 함께 접음) 배위 등급 A~D"),
        ("ni_grade_mono","신 파이프라인 (단량체 + Ni) 배위 등급"),
        ("ni_grade_dimer","동형이량체 + Ni. 'A-계면' = 두 사슬이 나눠 배위한 Cys4 (CooC1 형)"),
        ("plddt_motif",  "모티프 구간 평균 pLDDT — plddt_mean 과 비교해야 의미가 있다"),
        ("fd_metal",     "Folddisco 금속 모티프 히트. ※ 비검출은 부재의 근거가 아니다"),
        ("rf2ppi",       "※ 대조군 0.266 으로 실패한 축. 참고만"),
        ("boltz_iptm",   "※ 대조군 0.300 인데 후보가 더 높았다. 순위에 쓰지 않는다"),
    ], columns=["열", "뜻"])
    sheets = [("전체", T)]
    sheets.append(("3축이상", T[T["근거축수"] >= 3]))
    if "bl21_only" in T.columns:
        _m = T.mobile if "mobile" in T.columns else 0
        sheets.append(("균주특이_비프로파지", T[(T.bl21_only == 1) & (_m != 1)]))
    for c, nm in [("ni_grade_old", "Ni배위_AB_구"), ("ni_grade_mono", "Ni배위_AB_신"),
                  ("ni_grade_dimer", "Ni배위_이량체")]:
        if c in T.columns:
            sub = T[T[c].astype(str).str.startswith(("A", "B"))]
            if len(sub): sheets.append((nm, sub))
    sheets += [("열설명", DICT), ("출처", pd.DataFrame(SRC, columns=["파일", "내용", "행수"]))]
    with pd.ExcelWriter(XL, engine="openpyxl") as w:
        for nm, df in sheets:
            if df is None or not len(df): continue          # 빈 시트는 안 쓴다
            df.to_excel(w, sheet_name=nm[:31], index=False)
    print(f"\n  ★ {XL}")
    print("    시트: " + " / ".join(nm for nm, df in sheets if df is not None and len(df)))
except Exception as e:
    print(f"\n  엑셀 저장 실패 ({e}) — CSV 는 만들어졌다")
```
