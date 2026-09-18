# 발표용 확정 표 — 붙여넣기 전용

이 문서에는 **슬라이드에 그대로 올릴 표와 문장만** 담는다.
근거·유도 과정은 `ChCODH2_slides_boltz_ni.md` 와 `ChCODH2_Ni_insertase_report.md` 에 있다.

---

## 슬라이드 8 — Foldseek: ChCooC1 질의 결과

**제목**: Structural homologs of ChCooC1 in BL21(DE3) and K-12 MG1655

### 표 8-1 (메인)

| BL21(DE3) | TM | qcov | MG1655 | TM | qcov | Name |
|---|---|---|---|---|---|---|
| QJZ13931.1 | 0.708 | 0.98 | P37655 | 0.731 | 0.96 | BcsQ |
| QJZ11902.1 | 0.638 | 0.98 | P0AEZ3 | 0.631 | 0.98 | MinD |
| QJZ12661.1 | 0.638 | 0.99 | P0AF08 | 0.563 | 0.99 | ApbC / Mrp |
| QJZ11745.1 | 0.605 | 0.74 | P38134 | 0.573 | 0.74 | Etk |
| QJZ12620.1 | 0.562 | 0.74 | P76387 | 0.560 | 0.74 | Wzc |
| QJZ14620.1 | 0.554 | *0.16* ⚠ | P39337 | 0.550 | *0.16* ⚠ | YjgM |
| QJZ12684.1 | 0.518 | *0.21* ⚠ | P33368 | 0.545 | *0.21* ⚠ | YohF |
| QJZ12719.1 | 0.502 | 0.85 | P33030 | *0.414* | — | YeiR |
| **8** | | | **7** | | | |

**양성대조군** — Y19 `AHZ96930.1` CooC: TM **0.837**, qcov **0.99**, tcov **1.00**, fident **0.385**

### ★ qcov 를 넣으면 두 줄이 빠진다

`YjgM` 과 `YohF` 는 qcov **0.16 / 0.21** — ChCooC1 의 16~21% 조각만 정렬됐다.
TM 0.55 는 그 조각 안에서의 값이다. **실질은 BL21 6 / MG1655 5**, BcsQ 를 슈도진으로
빼면 **5 / 4** 다.

`Etk` 와 `Wzc` 의 tcov 0.23 은 정상이다 — 두 단백질이 720 잔기가 넘어 ChCooC1 이
ATPase 도메인에만 맞는다. **qcov 0.74 로 질의 쪽은 대부분 설명된다.**

### 진짜 차이는 TM 이 아니라 fident 에 있다

| | fident |
|---|---|
| Y19 CooC (양성대조군) | **0.385** |
| BL21 · MG1655 히트 전부 | **0.13 ~ 0.20** |

대조군만 서열 동일성이 두 배다. 나머지는 **P-loop ATPase 초과(ParA/MinD 계열)를
공유할 뿐**이고, CooC 도 같은 초과에 속하니 당연한 결과다.

> **BL21 의 히트들은 CooC 와 같은 초과에 속할 뿐, CooC 직교체가 아니다.**

- TM = Foldseek `alntmscore`, `--exact-tmscore 1`. **ipTM 이 아니다** — 두 구조의 겹침을 재는 측정값이다 (0–1).
- TM ≥ 0.5 = 같은 폴드 (Xu & Zhang 2010). 위 표는 **TM ≥ 0.5** 컷.
- 8 vs 7 의 차이는 **YeiR 하나**다. 컷 0.5 를 사이에 두고 갈렸을 뿐 **같은 단백질**이다 (fident 0.996).
- **BL21 에서 TM ≥ 0.9 에 닿는 단백질은 없다.**

### 표 8-2 (참고 — Y19 양성대조군)

| 균주 | TM ≥ 0.5 히트 | 최고 히트 | TM |
|---|---|---|---|
| *C. amalonaticus* Y19 | 15 | AHZ96930.1 **CooC** | **0.837** |
| BL21(DE3) | 8 | QJZ13931.1 BcsQ | 0.708 |
| K-12 MG1655 | 7 | P37655 BcsQ | 0.731 |

> Y19 는 실제 CooC 를 TM 0.837 로 되찾는다 (결정구조 프레임에서는 0.913).
> **질의가 작동한다는 증거**이고, 동시에 **BL21 의 0.708 이 CooC 수준이 아니라는 기준선**이다.

### 슬라이드에 넣을 문장

> All top-ranking matches are present in **both** strains.

> → Fold-level search does not discriminate BL21 from MG1655.
> **We are therefore attempting a PPI-based search.**

---

## 슬라이드 9 — Foldseek 한 건만 다르다: BcsQ

| | BL21(DE3) | K-12 MG1655 |
|---|---|---|
| ID | QJZ13931.1 | P37655 |
| 길이 | 250 aa | **242 aa** (AFDB 모델) |
| 상태 | 온전 | **pseudogene** — codon 6 에 조기 종결, N 말단 8 잔기 결손 |
| 표현형 | — | 셀룰로스 합성 소실, reading frame 복구 시 회복 |
| 출처 | — | UniProt `PUTATIVE PSEUDOGENE`; 온전형은 P0DP92 (PMID 24097954) |

**말할 것**

> 폴드 수준에서 8 : 7 로 보였지만, 한 건은 **K-12 에서 기능이 없다.**
> UniProt 은 슈도진을 플래그와 함께 레퍼런스 프로테옴에 남기기 때문에 AFDB 가 모델을 만들었고,
> Foldseek 은 그 모델을 정상 히트로 돌려준 것이다. 기능 기준으로 세면 **BL21 8 / MG1655 6.**

**단서 (반드시 같이 말할 것)**

- BcsQ 의 알려진 역할은 **셀룰로스 합성효소 국재화**다. 금속 자리 근거는 없다.
- TM 0.708 은 폴드 수준이지 활성부위 일치가 아니다.
- 조기 종결은 bcs 오페론 하류에 **polar effect** 를 준다 — BcsQ 단독 효과로 볼 수 없다.

---

## 슬라이드 10 — Boltz co-folding: 왜 했는가

### 논리 (한 문단)

> MG1655 대비 **BL21 특이적 + 저유사도 단백질 326개**를 대상으로 Boltz 를 사용해
> ChCODH2 와 **단량체끼리 도킹하면서 Ni 을 함께 도킹**했다.
> 참조 자료로 쓸 수 있으며, 분석 결과 **기존에 어노테이션되어 있지 않던 HypA 폴드**를 추가로 찾았다.

### English (academic)

> Using Boltz, we co-folded each of the 326 BL21-specific or low-similarity proteins with
> ChCODH2 as a binary complex in the presence of Ni²⁺. The resulting models provide a
> reference set of candidate Ni-coordination sites, and the analysis additionally recovered a
> **HypA-like fold** that had not been annotated in this genome.

**단서**: 폴드를 식별한 것은 Foldseek 이다. BL21 에는 이미 어노테이션된 HypA 가 있다
(`QJZ13171.1`, `QJZ13424.1`). 주장은 "BL21 에 HypA 가 없다"가 아니라
**"이 특정 단백질(QJZ12568.1)이 HypA 폴드를 가지며 MG1655 에 부재한다"** 이다.

---

## 슬라이드 11 — Boltz Ni 배위 등급 결과

**제목**: Ni coordination grades from co-folded models

| Grade | 단백질 | BL21 | MG1655 | Y19 | 거리 (Å) |
|---|---|---|---|---|---|
| A | Hsp33 / HslO | O | − | O | 2.18 / − / 2.17 |
| **A** | **QJZ12568.1** | **O** | **−** | **−** | **2.23** |
| A | GspE | O | − | − | 2.29 |
| B | HisB | O | O | O | 1.91 / 1.92 / 1.86 |
| B | GmhB | O | − | − | 1.93 |
| B | YehQ | − | O | − | 2.09 |
| B | Tgt | O | O | O | 2.14 / 2.12 / 2.14 |
| B | SWIM zinc finger | − | − | O | 2.19 |
| C | HpaA | O | − | − | 1.96 |
| D | phage tail | O | − | − | 1.90 |
| D | YeiR | O | − | − | 1.91 |
| D | Ser/Thr phosphatase | O | − | − | 1.93 |
| D | NagA | O | − | − | 1.94 |
| D | PhnP | O | O | O | 2.05 / 2.08 / 2.00 |
| D | QJZ11612.1 | O | − | − | 2.09 |

**15 단백질 / 22 행** (한 단백질이 여러 균주에서 잡히면 행이 늘어난다)

### 등급 정의

| 등급 | 조건 |
|---|---|
| A | Cys ≥ 4 |
| B | Cys 3, 또는 Cys 2 + His |
| C | Cys 1–2, 또는 His ≥ 3 |
| D | 그 외 (산소 원자가 Ni 과 도킹된 경우) |

배위 판정선 `TIGHT = 2.6 Å`, 후보 단백질 사슬만 센다 (ChCODH2 쪽은 제외).

### 반드시 같이 말할 것

- **이 축은 균주를 못 가른다.** HisB · Tgt · PhnP 는 **세 균주 모두**에서 같은 등급으로 나온다.
- Grade A 는 **서로 다른 단백질 3종**이다. HslO 와 GspE 는 이미 알려진 아연 자리이므로,
  **기능 미상은 QJZ12568.1 하나뿐이다.**
- Folddisco 금속 모티프 목록과 겹치는 것은 22 행 중 **6 행**, grade A 4 행 중 **1 행**뿐이다.

---

## 슬라이드 12 — 최상위 후보

### QJZ12568.1

| 항목 | 값 |
|---|---|
| 어노테이션 | hypothetical protein HO396_09830 |
| 길이 | 133 aa |
| Boltz 등급 | **A** (Cys4), Ni–S 2.23–2.40 Å |
| 배위 잔기 | Cys99 · Cys102 · Cys112 · Cys115 (**전부 후보 단백질 자신의 사슬**) |
| MG1655 | 프로테옴·게놈 **둘 다 부재** |
| ~~폴드~~ | ~~HypA-like~~ → **철회**. 아래 참조 |

### 근거는 셋이다

1. Cys4 Ni 배위 자리 (Boltz, 2.23 Å)
2. MG1655 부재 — 프로테옴·게놈 둘 다. 실험 전제(BL21 lysate 만 활성 회복)와 방향이 같다
3. 기능 미상 + 133 aa 소형 — grade A 중 유일하게 기능이 안 알려진 것

> **Cys4 재질의(folddisco) 히트는 근거가 아니다.** 그 질의의 템플릿이 이 단백질 자신이라
> 자기참조다.

### ★ HypA 폴드 주장은 철회한다

`.m8` 을 열어 확인한 결과다.

| 비교 | TM | **qcov** |
|---|---|---|
| 진짜 HypA 파랄로그끼리 (`QJZ13171.1` ↔ `QJZ13424.1`) | **0.900** | **0.98 / 1.00** |
| QJZ12568.1 → `QJZ13424.1` (HypA #2) | 0.664 | **0.26** |
| QJZ12568.1 → `QJZ13171.1` (HypA #1) | **0.499** | 0.38 |

진짜 파랄로그는 전체가 정렬되는데(qcov ≈ 1.0) 이 후보는 **26% 만** 정렬되고,
두 HypA 중 **하나에는 0.5 컷 아래**다. **직교체도 파랄로그도 아니다.**

> **TM 은 qcov 없이 읽으면 안 된다.** qcov 가 낮으면 TM 이 높아도
> "일부 도메인만 공유"라는 뜻이다.

### 그림 캡션

> **Figure.** Boltz co-folding model of ChCODH2 (grey, chain A) with candidate QJZ12568.1
> (blue, chain B) and Ni²⁺. The Ni ion (green sphere) is coordinated by four cysteine
> thiolates (Cys99, Cys102, Cys112, Cys115; sticks, S in yellow) contributed entirely by
> QJZ12568.1, at 2.23 Å. Surfaces are shown semi-transparent.

**해석 한 줄**: 후보가 **자체 Cys4 자리로 Ni 을 쥐고 있다.** metallochaperone 에 기대하는 모습이고,
ChCODH2 로 넘기기 전의 적재 상태로 읽힌다. 다만 이것은 **도킹 모델의 기하일 뿐 전달의 증거는 아니다.**

---

## 슬라이드 13 — 다음 단계 (PPI based search)

### 가설 / 방법 (English, academic)

**Hypothesis**

> Among the 326 proteins that are either unique to BL21(DE3) or show low similarity to their
> K-12 MG1655 counterparts, at least one engages in a protein–protein interaction with ChCODH2.

**Rationale**

> Fold-level search (Foldseek, ChCooC1 as query) returned matches that are present in both
> strains, and therefore cannot account for the strain-dependent restoration of ChCODH2
> activity. A discriminating candidate must instead be identified by its capacity to form a
> complex with ChCODH2.

**Approach**

> Each candidate is co-folded with ChCODH2 and scored for interface quality, and the subset
> predicted to form a stable interface is examined for Ni-coordinating residues at the
> candidate chain. Homodimeric re-prediction is applied to candidates whose coordination site
> is expected to be shared across a dimer interface, with ChCooC1 itself as the positive control.

---

## 쓰는 파일

| 슬라이드 | 원본 |
|---|---|
| 8 | `result/table/foldseek_cooc1_bl21.m8`, `foldseek_cooc1_mg1655.m8`, `foldseek_cooc1_y19.m8` |
| 9 | UniProt P37655 / P0DP92, PMID 24097954 |
| 11 | `result/table/ni_site_grade_named.csv` |
| 12 | `boltz_results_inputs/predictions/ChCODH2_WT__QJZ12568.1/ChCODH2_WT__QJZ12568.1_model_0.cif` |

---

## 슬라이드 14 (마지막) — 종합: 후보와 근거

**제목**: Candidates for expression and in vitro Ni incorporation

### 선정 원칙

> Candidates were selected for **plausibility of metal handling**, not for strain
> specificity. Gene presence/absence did not account for the phenotype, and the
> in vitro reconstitution assay asks whether a protein inserts Ni — not whether
> MG1655 also carries the gene.

**근거 축은 둘만 썼다** — 재현성이 확인된 측정 두 가지.
Folddisco `metal_rmsd` · Boltz Ni 배위 도너 등급.

### 후보 표

| # | Protein | Ni-site grade | Evidence axis | Status |
|---|---|---|---|---|
| **1** | **QJZ12568.1** hypothetical HO396_09830, 133 aa | **A** — Cys₄, 2.23 Å | Ni coordination **+ absent from MG1655** (proteome and genome) **+ uncharacterised** | **Primary** |
| 2 | QJZ13396.1 **GspE**, 497 aa | **A** — Cys₄, 2.29 Å | Ni coordination; site quality highest | Established Zn enzyme |
| 3 | QJZ13803.1 **HslO** (Hsp33) | **A** — Cys₄, 2.18 Å | Ni coordination | **No strain discrimination** — Y19 ortholog identical |
| 4 | QJZ10977.1 **ErpA** | — | Folddisco `metal_rmsd` **0.133**; Fe-S insertion, C-cluster is [NiFe₄S₄] | **No Ni captured** |
| 5 | QJZ11953.1 **YchJ** / QJZ12547.1 **YecA** | — | Folddisco `metal_rmsd` 0.155 / **0.123** (lowest); SEC-C / zinc-ribbon, uncharacterised | **No Ni captured** |
| — | QJZ12719.1 **YeiR** | D | **Only protein recovered by both axes**; COG0523 / G3E family (HypB, UreG, CooC) | **Reserve** — MG1655 carries a near-identical copy (fident 0.996) |

**1–3 은 Boltz 축, 4–5 는 Folddisco 축이다.** 4·5 는 Ni 이 후보 쪽 사슬로 오지 않아
등급이 없다. 근거의 종류가 다르다는 것을 표에서 읽을 수 있어야 한다.

### 제외한 것

| | 이유 |
|---|---|
| Tgt · HisB · GmhB | grade B 로 자리는 좋으나 **세 균주가 모두 동일하게 잡는다**. 기능이 확립된 효소라 사전확률이 낮다 |
| QJZ14223.1 LysR | MG1655 최고 히트가 파랄로그 (fident 0.345) |
| QJZ13177.1 | pLDDT 56.7 |

### 무엇이 음성이었는지 — 먼저 말할 것

| 축 | 결과 |
|---|---|
| Foldseek (ChCooC1 질의) | 상위 히트가 **양쪽 균주에 모두 존재** → 균주 변별 불가 |
| Boltz Ni 배위 등급 | HisB · Tgt · PhnP 가 **세 균주 모두** 동일 등급 → 균주 변별 불가 |
| `ipTM` | 대조군 CooC1 = 0.300 인데 후보들이 그 위 → **폐기** |
| `ligand_ipTM` | 계면 배치 0/348, 사슬 길이와 ρ = −0.807 → **폐기** |
| RF2-PPI | 대조군이 떨어짐 → 순위를 근거로 못 씀 |

> 여섯 축 중 다섯이 0 이었다. **후보는 통과한 축에서만 나왔다.**

### 다음 단계

```
1.  PPI-based search
    326 BL21-specific / low-similarity proteins × ChCODH2
    → interface quality, then Ni coordination on the candidate chain

2.  Homodimeric re-prediction
    G3E chaperones (YeiR, HypB, UreG, CooC) share the metal site across
    the dimer interface; monomeric folding penalised them
    → positive control: ChCooC1 itself. If it does not grade at the
      interface, the track is not used

3.  Wet-lab
    Express candidates 1-5, in vitro Ni incorporation into apo-ChCODH2
```

### 발표에서 반드시 같이 말할 것

1. **예측 구조다.** 결정구조가 아니다. Boltz 가 Ni 위치를 정했다
2. **Cys₄ 는 Zn 과 Ni 이 공유하는 자리다.** 자리의 종류만으로 두 금속을 못 가른다 —
   grade A 셋 중 둘이 알려진 아연 단백질인 것이 그 증거다
3. **입구가 좁았다.** Track B 326개만 도킹했다. identical + high-similarity **3,746개(91%)는 평가받지 못했다**
4. **단량체로 접었다.** 이량체 계면 자리를 갖는 단백질에 불리했다
5. `O / −` 는 **도킹·채점 여부**이지 균주 내 유무가 아니다
