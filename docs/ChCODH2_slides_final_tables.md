# 발표용 확정 표 — 붙여넣기 전용

이 문서에는 **슬라이드에 그대로 올릴 표와 문장만** 담는다.
근거·유도 과정은 `ChCODH2_slides_boltz_ni.md` 와 `ChCODH2_Ni_insertase_report.md` 에 있다.

---

## 슬라이드 8 — Foldseek: ChCooC1 질의 결과

**제목**: Structural homologs of ChCooC1 in BL21(DE3) and K-12 MG1655

### 표 8-1 (메인)

| BL21(DE3) | TM | MG1655 | TM | Name |
|---|---|---|---|---|
| QJZ13931.1 | 0.708 | P37655 | 0.731 | BcsQ |
| QJZ11902.1 | 0.638 | P0AEZ3 | 0.631 | MinD |
| QJZ12661.1 | 0.638 | P0AF08 | 0.563 | ApbC / Mrp |
| QJZ11745.1 | 0.605 | P38134 | 0.573 | Etk |
| QJZ12620.1 | 0.562 | P76387 | 0.560 | Wzc |
| QJZ14620.1 | 0.554 | P39337 | 0.550 | YjgM |
| QJZ12684.1 | 0.518 | P33368 | 0.545 | YohF |
| QJZ12719.1 | 0.502 | P33030 | *0.414* | YeiR |
| **8** | | **7** | | |

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
| 폴드 | **HypA-like** (Foldseek) |

### 근거는 셋이다

1. Cys4 Ni 배위 자리 (Boltz, 2.23 Å)
2. HypA 폴드 — HypA 는 [NiFe] 하이드로게나제 성숙의 **Ni 삽입 인자**다
3. MG1655 부재 — 실험 전제(BL21 lysate 만 활성 회복)와 방향이 같다

> **넷이 아니다.** Cys4 재질의(folddisco) 히트는 그 질의의 템플릿이 이 단백질 자신이라
> 자기참조다. 근거에서 뺐다.

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
