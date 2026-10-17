# In silico search for a Ni insertase acting on ChCODH2

**Progress handover and request for methodological advice**
Prepared 2026-10-06

---

## 0. Status at a glance

| | |
|---|---|
| **Settled** | Fold-level search (Foldseek) is complete and reproducible in two reference frames. Ni-coordination grading over 411 co-folded structures is complete. Strain assignment of the lead candidate is confirmed at both the proteome and genome level |
| **Withdrawn** | The HypA-fold assignment for the lead candidate (see §4). The RF2-PPI ranking, the Boltz `ipTM` and `ligand_ipTM` axes — all three failed their controls |
| **Running** | Homodimeric re-prediction of the G3E-family candidates, with ChCooC1 as the positive control; a sequence-motif (CXC / CXXC / Cys₄) screen across all three proteomes, feeding the top 400 into Boltz |
| **Open** | Whether to pursue interface prediction at all (§5); how to make three heterogeneous structure sets commensurate (§3.3); whether the Folddisco query can be rebuilt usefully (§3.1) |

**Reproducibility note.** Every number in this document was read back from the
output files rather than from notes. Two earlier claims did not survive that
check and have been removed; they are listed in §4 so the record is honest.

---

## 1. The biological problem

*Carboxydothermus hydrogenoformans* CODH-II (ChCODH2) carries a [NiFe₄S₄]
C-cluster. Nickel insertion is normally attributed to the maturation factor
**CooC1**, a G3E/SIMIBI-family GTPase.

Our experimental observation does not fit that assignment:

- ChCODH2 is reconstituted into an **active, Ni-loaded** form in **BL21(DE3)**
  lysate **without CooC1 expression**
- The same reconstitution **fails in K-12 MG1655** lysate
- *Citrobacter amalonaticus* **Y19**, an enterobacterium that carries a genuine
  Ni-CODH, is used as a third reference strain

**Question.** Which BL21(DE3) protein substitutes for CooC1?

We treat this as a search problem over three proteomes:
BL21(DE3) **4,110** · MG1655 **4,300** · Y19 **5,325**.

---

## 2. What we ran, and what the controls said

| Track | Method | Control | Outcome |
|---|---|---|---|
| A | **RF2-PPI** coevolution | CooC1–ChCODH2 (a documented interacting pair) | **Failed.** Control scored **0.266** (SD 0.048, n=3); our exclusion threshold is < 0.30. **531 / 3,625** screening pairs scored above the control |
| B | **Boltz-2** co-folding with Ni²⁺ | CooC1 | **ipTM discarded.** Control ipTM 0.300 while candidates scored above it; replicate spread for identical input reached 0.309. **Ni coordination geometry retained** |
| C | **Folddisco** structural-motif search | ChCooC1 Cys112/Cys114 | **Query failed** (see §3.1) |
| D | **Foldseek** whole-structure search | ChCooC1 | **Works** — recovers Y19 CooC at TM 0.837 (0.913 against the crystal frame) |

Track A's control failure is not an MSA-depth problem: the control pair's paired
MSA depth is **1,451 genomes**, above the 75th percentile of the screening set
(996 median; 1,392 at p75). *cooC* and *cooS* are operonic, so this should be
the best case for coevolution, not the worst.

---

## 3. Three problems we would value your opinion on

### 3.1 Folddisco: a motif query taken from a holo crystal structure

We built the query from the **3kji** CooC1 crystal structure (Cys112/Cys114,
metal-bound) and searched **predicted** structures. It returned 80 hits in
BL21, none of them CooC-like, and **it does not retrieve the Y19 CooC** that
Foldseek finds at TM 0.837.

Diagnosing it, the query geometry is simply not present in the predictions:

| Structure | Cys112 SG – Cys114 SG |
|---|---|
| 3kji crystal (holo) | **3.58 Å** |
| ChCooC1, AlphaFold DB | **5.30 Å** |
| Y19 CooC, AlphaFold DB | **8.16 Å** |

Two AlphaFold predictions of the same site differ from each other by **2.86 Å**.
In the ChCooC1 model the query residues sit at **pLDDT 41.0 / 45.0** against a
chain mean of 79.8 (residues 106–115 are all in the 40s). The 80 hits that *were*
returned have matched residues averaging **pLDDT 85.3**.

So the query was drawn from a disordered loop and the hits are well-ordered sites
that happen to fit it. Relaxing `-d` (0.5 → 2.0) and `--rmsd` (1.0 → 2.0) changes
nothing — the maximum `min_rmsd` over all hits is 0.745, with nothing between
0.75 and 2.0, so the bottleneck is the index lookup, not the RMSD filter.

> **Q1.** Is structural-motif search over predicted models viable at all for
> metal sites, given that the apo prediction does not reproduce the holo
> geometry? Is there a recommended practice — building the query from a
> predicted structure rather than a crystal, filtering the index by pLDDT,
> or scaling tolerance by local confidence?

> **Q2.** CooC1's Ni site is built from **both protomers** (Cys112/Cys114 of each
> chain). We indexed monomer predictions. Does Folddisco have a sanctioned way to
> search interface motifs, or must we build dimer models first?

### 3.2 Foldseek: the fold is not the discriminator

ChCooC1 as query, TM ≥ 0.5, `--exact-tmscore 1`:

| BL21(DE3) | TM | qcov | MG1655 | TM | qcov | Protein |
|---|---|---|---|---|---|---|
| QJZ13931.1 | 0.708 | 0.98 | P37655 | 0.731 | 0.96 | BcsQ |
| QJZ11902.1 | 0.638 | 0.98 | P0AEZ3 | 0.631 | 0.98 | MinD |
| QJZ12661.1 | 0.638 | 0.99 | P0AF08 | 0.563 | 0.99 | ApbC / Mrp |
| QJZ11745.1 | 0.605 | 0.74 | P38134 | 0.573 | 0.74 | Etk |
| QJZ12620.1 | 0.562 | 0.74 | P76387 | 0.560 | 0.74 | Wzc |
| QJZ14620.1 | 0.554 | *0.16* | P39337 | 0.550 | *0.16* | YjgM |
| QJZ12684.1 | 0.518 | *0.21* | P33368 | 0.545 | *0.21* | YohF |
| QJZ12719.1 | 0.502 | 0.85 | P33030 | *0.414* | — | YeiR |

**Positive control** — Y19 `AHZ96930.1` CooC: TM **0.837**, qcov 0.99, tcov 1.00,
**fident 0.385**.

Two observations:

1. **Every hit pairs one-to-one across the two strains.** The fold axis does not
   separate BL21 from MG1655. The single functional difference we could establish
   is BcsQ, which is a **documented pseudogene in K-12** (premature stop at codon
   6; UniProt P37655 carries the flag, intact protein is P0DP92).
2. **What separates the control is not TM but sequence identity** — fident 0.385
   versus 0.13–0.20 for every *E. coli* hit. The *E. coli* matches share the
   P-loop ATPase (ParA/MinD) superfamily that CooC also belongs to, and nothing
   more.

> **Q3.** For a query whose fold is as promiscuous as the SIMIBI/P-loop NTPase
> class, is there a principled Foldseek-side way to separate "same superfamily"
> from "same function"? Would you restrict the alignment to the metal-binding
> insertion, use a different normalisation, or treat this as a question Foldseek
> should not be asked?

> **Q4.** Is "BL21 contains no CooC ortholog" a defensible negative result from
> this evidence, and if not, what would you require before stating it?

### 3.3 Mixed model provenance across three proteomes

Our structure sets are heterogeneous: AlphaFold DB `AF-{acc}-F1-model_v6` and
`model_v4`, plus in-house ColabFold `cf_{UniParc}.pdb`. We found four cases where
a motif was detected in one strain and not the other **at fident 1.000** — the
difference was the model version (v6 vs v4), not the biology. We removed those
four from the analysis.

> **Q5.** For a three-proteome comparison, is there a defensible way to make the
> structure sets commensurate — re-predicting everything under one pipeline, or
> restricting comparisons to matched model versions?

---

## 4. What survived, and the candidate

The only axis that reproduced across replicates was the **Ni coordination
geometry** in the Boltz models. Across orthologs in different strains, folded
separately on different GPUs, the same residues are held at the same distance:

| Protein | BL21 | MG1655 | Y19 | Coordinating residues |
|---|---|---|---|---|
| HslO | 2.18 Å | — | 2.17 Å | Cys230/232/263/266 (identical) |
| HisB | 1.91 | 1.92 | 1.86 | Cys93/101/103, His95 (identical) |
| Tgt | 2.14 | 2.12 | 2.14 | Cys302/304/307, His333 (identical) |

Δ 0.01–0.08 Å, against an ipTM spread of 0.28 for the same inputs. We therefore
graded sites by donor type within 2.6 Å of the docked Ni, counting the candidate
chain only (A: Cys ≥ 4; B: Cys 3 or Cys2+His; C: Cys 1–2 or His ≥ 3; D: otherwise).

The grading recovers published metal sites correctly — Tgt's Cys302/304/307/His333
and PhnP's His78/143/Asp164 match the literature — which is reassuring, but it
also means **the scheme finds zinc sites**, since Cys₄ is shared by Zn and Ni.
Two of the three grade-A proteins are established zinc enzymes (HslO, GspE).

**The remaining candidate is `QJZ12568.1`** — hypothetical protein HO396_09830,
133 aa, Cys₄ site at 2.23 Å, absent from the MG1655 proteome *and* genome, and
the only grade-A protein with no assigned function.

We had also claimed it carried a HypA-like fold. On checking the alignments this
week we withdrew that: the match covers **26% of the query** (TM 0.664), and
against BL21's second HypA paralog it scores **0.499**, below the fold cutoff —
whereas the two genuine HypA paralogs align over 98–100% of their length at
TM 0.900. A match over a quarter of the query is a shared subdomain, not a fold
assignment.

---

## 5. Where we are heading, and the question behind it

Because fold-level search returns proteins present in both strains, we are moving
to a **PPI-based search**: co-folding each of the 326 BL21-specific or
low-similarity proteins with ChCODH2 and scoring the interface, then examining
Ni coordination on the candidate chain. In parallel we are re-predicting the
G3E-family candidates as **homodimers**, since their metal sites are shared
across the dimer interface and monomeric folding penalised them; ChCooC1 itself
is the positive control, and if it does not grade at the interface we will not
use that track.

> **Q6.** Given that the coevolution control failed and ipTM proved unusable,
> is interface prediction a reasonable axis here at all — or is a transient
> metallochaperone–client interaction simply outside what these methods resolve?
> If so, what would you do instead?

---

## Appendix — tools and settings

| Tool | Settings |
|---|---|
| Foldseek | `easy-search --format-output "query,target,fident,alnlen,qcov,tcov,evalue,bits,prob,alntmscore,lddt" -e 10 --max-seqs 2000 --exact-tmscore 1` |
| Folddisco | `-d 0.5 / 1.0 / 2.0`, `-a 5.0`, `--rmsd 1.0 / 1.5 / 2.0`, `--covered-node 0 / 4` |
| MMseqs2 | `easy-search -s 7.5`; genome check with `--search-type 2` |
| Boltz-2 | `--use_msa_server --max_msa_seqs 2048 --recycling_steps 3 --diffusion_samples 1`; Ni as CCD `NI` |
| Structures | AlphaFold DB `model_v4` / `model_v6`; in-house ColabFold for the remainder |

**Counts.** 411 Boltz structures (321 Track B + 27 focus + 50 cross-strain + 13
replicates). Ni is placed within 3.5 Å of protein in 411/411, so placement alone
carries no information; 22 structures place it on the candidate chain. Folddisco
metal motif: 29 BL21 / 21 MG1655 / 30 Y19 hits, invariant across all tolerance
settings tested.

---

## 부록 B — 발표 구성 방향 (내부용)

자문 미팅용 슬라이드를 짤 때의 뼈대다. **이 문서를 그대로 읽는 발표가 아니라,
질문 여섯 개에 답을 받아오는 것이 목적**이므로 구성이 다르다.

### 원칙 셋

1. **실패를 먼저 말한다.** 후보 자랑이 아니라 "왜 이렇게 됐는지 봐 달라"가 목적이다.
   대조군 표(§2)를 3번 슬라이드에 두고 시작한다
2. **질문 슬라이드에 가장 많은 시간을 남긴다.** 앞부분은 맥락일 뿐이다.
   12장이면 7장까지를 15분 안에 끝낸다
3. **숫자는 파일에서 읽은 것만 올린다.** 근사값은 `≈` 로 쓰고, 미확인은 비워 둔다

### 슬라이드 뼈대 (12–13장)

| # | 내용 | 핵심 한 줄 |
|---|---|---|
| 1 | 표지 | 제목 + "세 축이 대조군을 통과하지 못했다" |
| 2 | 생물학적 문제 | CooC1 없이 BL21 에서만 Ni 조립이 된다 |
| 3 | 검색 설계 | 세 프로테옴 · 네 방법 · 각각의 양성대조군 |
| **4** | **대조군 결과표 (§2)** | **네 축 중 셋이 자기 대조군에 떨어졌다** ← 여기가 시작점 |
| 5 | Folddisco 실패 (§3.1) | SG–SG 3.58 / 5.30 / 8.16 Å + pLDDT 41·45 |
| 6 | Foldseek 표 (§3.2) | 8 : 7 이 1:1 로 짝지어진다 — 폴드는 균주를 못 가른다 |
| 7 | fident 가 가른다 | 대조군 0.385 vs 나머지 0.13–0.20 |
| 8 | 재현된 유일한 축 | 직교체 간 Ni 거리 Δ 0.01–0.08 Å (ipTM 은 0.28 벌어짐) |
| 9 | 등급 결과 | grade A 3종, 그중 둘이 알려진 Zn 자리 |
| 10 | 후보 + 철회한 것 | QJZ12568.1 / HypA 주장 철회 (qcov 0.26) |
| **11–12** | **질문 Q1–Q6** | **한 장에 셋씩. 각 질문은 바로 앞 슬라이드를 가리킨다** |
| 13 | 다음 단계 | PPI 기반 탐색 + 이량체 재예측, 대조군 설계까지 |

### 슬라이드별 주의

- **4번** — 네 행 다 "대조군 점수 / 판정 기준 / 결과" 세 칸으로 통일한다.
  RF2-PPI 0.266 (기준 0.30), Boltz ipTM 0.300 (후보가 그 위) 가 나란히 보여야 한다
- **5번** — 거리 세 개와 pLDDT 두 개가 전부다. 그림은 ChCooC1 예측 구조의
  해당 루프를 B-factor(pLDDT)로 칠한 것 하나면 충분하다
- **6번** — qcov 열을 반드시 같이 올린다. YjgM·YohF 가 0.16·0.21 로 빠지는 게
  이 슬라이드의 결론이다
- **9번** — Tgt 의 Cys302/304/307·His333 이 문헌 금속 자리와 일치한다는 걸
  한 줄 넣는다. 등급 체계가 작동한다는 유일한 외부 검증이다
- **10번** — 철회를 숨기지 않는다. 진짜 HypA 파랄로그끼리 qcov 0.98–1.00 에
  TM 0.900 인데 우리 후보는 qcov 0.26 이라는 비교를 그대로 보여준다
- **11–12번** — 질문마다 "우리가 생각한 답"을 한 줄 달아 두면 논의가 빨라진다.
  빈 질문만 던지면 "글쎄요"로 끝난다

### 올리지 말 것

- `ipTM`·`ligand_ipTM` 의 세부 분포 — 폐기한 축이다. 한 줄로 끝낸다
- 22행 전체 등급표 — grade A 와 "세 균주 모두 잡힌 셋"만 보이면 된다
- 아직 확정 안 된 수치: HypA 스윕의 정확 TM, Y19 `AKE61643.1` 의 정체

### 질문을 던질 때의 순서

> Q1·Q2 (Folddisco) → Q3·Q4 (Foldseek) → Q5 (모델 출처) → Q6 (계면 예측)

앞의 다섯은 **도구를 어떻게 쓰느냐**이고 Q6 만 **이 접근이 성립하느냐**이다.
시간이 모자라면 Q1·Q3·Q6 만 묻는다 — 각 축에서 하나씩이다.
