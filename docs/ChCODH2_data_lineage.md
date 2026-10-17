# 데이터 계보 — 원본 → 컷오프 → 결과

`BASE = /mnt/af2results/mingyu/workspace/ppi_discovery`
아래 경로는 전부 `BASE` 기준 상대경로.

---

## A. 공통 입력 (컷오프 이전)

| 항목 | 경로 | 건수 | 비고 |
|---|---|---|---|
| BL21(DE3) 프로테옴 | `TOOLS/database/bacteriaDB/inhouseDB/bl21_db_match_qjz.faa` | **4,110** | QJZ* — 구조 DB 와 ID 계열 일치 |
| Y19 프로테옴 | `TOOLS/database/bacteriaDB/inhouseDB/y19_db_match.faa` | **5,325** | AKE*/AGE*/AHZ* |
| MG1655 프로테옴 | `TOOLS/database/protein_list/mg1655_protein.faa` | **4,300** | 공식 K-12 어노테이션 |
| MG1655 게놈 | `input/external/MG1655_U00096.3.fna` | 4.7 MB | NCBI efetch |
| BL21 구조 | `TOOLS/database/bacteriaDB/structures_UP000503272/` | 4,096 | AlphaFold **model_v6** |
| Y19 구조 | `TOOLS/database/bacteriaDB/structures_UP000034085/` | 5,283 | AlphaFold |
| MG1655 구조 | (folddisco 인덱스만, 공식 배포본) | — | **model_v4** — 버전 불일치 |
| CooC1 구조 | `WS/seek_ni_insertase/input/3kji.pdb` | 2 사슬 | A·B 각 112·114 = CYS |
| bait 서열 | `input/seq/baits.fasta` | 8 | ChCODH2_WT + CooC1 + 절편 6 |
| 구조 tid ↔ GenBank | `result/table/id_crosswalk_struct_to_genbank.csv` | — | xlsx 기반, BL21/Y19 100% |
| 균주 분류 | `result/table/classify_MG1655_vs_BL21DE3.csv` | 4,300 | **MG1655 기준** — BL21 ID 와 조인 불가 |

---

## B. Track A — RF2-PPI (공진화)

| # | 단계 | 원본 | 건수 | 적용한 컷오프 | 남은 것 | 결과 파일 |
|---|---|---|---|---|---|---|
| 1 | bait 상동체 검색 | `result/search/bait_hits.m8` | 43,036 행 | taxid 당 bits 최고 1건 | 2,929 유전체 | (메모리) |
| 2 | **feasibility gate** | 위 | 2,929 | **깊이 ≥ 500** | 통과 | — |
| 3 | prey 상동체 검색 | `result/search/prey_hits.m8` | **62,642,570 행** | bait 보유 유전체만 + taxid 당 최고 | — | (메모리) |
| 4 | paired MSA | `result/paired/*.a3m` | 3,730쌍 시도 | **공유 유전체 ≥ 50** + hhfilter `-id 90` | 3,626 | `result/table/paired_msa_stats.csv` |
| 5 | 추론 | `result/rf2ppi/input_rep{1,2,3}.log` | 3,627 × 3 | — | 10,881 | `result/table/trackA_RF2PPI_all_pairs.csv` |
| 6 | prey 단위 집계 | 위 | 3,625 prey | bait 합집합 최대값 | 3,625 | `result/table/trackA_RF2PPI_ranked.csv` |
| 7 | **판정 컷오프** | 위 | 3,625 | **mean ≥ 0.74** (논문 95% precision) | **0개** | — |
| 8 | 양성대조군 | `CooC1__ChCODH2_WT.a3m` 깊이 1,451 | 1쌍 ×3 | **≥ 0.70** | **0.266 — 실패** | — |
| 9 | 구제 시도 | `result/paired/ctrlvar/` | 6 변형 ×3 | ≥ 0.70 | **최고 0.251 — 실패** | `result/table/trackA_control_variants.csv` |
| 10 | 이량체 재시도 | `result/paired/dimer/tandem.a3m` | 1쌍 ×3 | ≥ 0.70 | **진행 중** | `result/paired/dimer/in_t*.log` |

**현 상태:** 7번에서 0개, 8·9번 대조군 실패 → **순위를 후보 지명에 쓰지 않음.** 10번 결과 대기.

---

## C. Track B — Boltz-2 (co-folding)

| # | 단계 | 원본 | 건수 | 컷오프 | 남은 것 | 결과 파일 |
|---|---|---|---|---|---|---|
| 1 | prey 분류 | `result/table/classify_*.csv` | 4,110 | strain-specific + low-similarity | **326** → Track B | — |
| 2 | 입력 YAML | `result/boltz/inputs/` | 326 | **bait+prey ≤ 1,830 잔기** | 323 | — |
| 3 | 추론 | `result/boltz/out/` | 323 | `--max_msa_seqs 2048` | **321** (OOM 2건 손실) | `result/table/trackB_boltz2_ranked.csv` |
| 4 | focus (모티프 보유) | `result/boltz/inputs_focus/` → `out_focus/` | 27 | folddisco metal 히트 중 **pLDDT ≥ 70** | 27 | — |
| 5 | cmp (타 균주) | `result/boltz/inputs_cmp/` → `out_cmp/` | 51 | 길이 초과 1 제외 | **50** | `result/table/metal_motif_3strain_boltz.csv` |
| 6 | rep (반복+대조군) | `result/boltz/inputs_rep/` → `out_rep/` | 13 | `--diffusion_samples 5` | 13 (65 구조) | `result/table/boltz_replicate_summary.csv` |
| 7 | dimer (올리고머) | `result/boltz/inputs_dimer/` → `out_dimer/` | 5 | ≤ 1,830 | **진행 중** | — |
| 8 | Ni 좌표 판독 | 위 전부의 `*_model_0.cif` | **411 구조** | **배위 ≤ 3.5 Å** | Ni 있음 411 | `result/table/boltz_ni_placement.csv` |
| 9 | 배위 등급 | 위 | 22 (후보 쪽) | **≤ 2.6 Å 만 배위로 셈** | A:4 B:9 C:1 D:8 | `result/table/ni_site_grade.csv` |

**컷오프 판정:**
- `ipTM` — 통상 ≥0.8/0.6/0.4. **대조군 CooC1 = 0.300 (5표본)**, 후보 8개가 전부 위 → **순위에 안 씀**
- `ligand_ipTM` — 계면 배치 **0/348** , 길이와 ρ=−0.807 → **폐기**
- Ni 배위 좌표 — 직교체 간 Δ 0.01~0.08 Å → **[확정] 등급으로 사용**

---

## D. Track C — Folddisco (구조 모티프)

| # | 단계 | 원본 | 건수 | 컷오프 | 남은 것 | 결과 파일 |
|---|---|---|---|---|---|---|
| 1 | metal 질의 | `result/folddisco/metal_run01_{BL21,Y19,MG1655}.tsv` | — | `-q A112,A114` `-d 0.5` `-a 5.0` `--rmsd 1.0` `--top 2000` | **29 / 30 / 21** | — |
| 2 | ATP 질의 | `result/folddisco/atp_run02_*.tsv` | — | `-q A12,A13,A14` 동일 | **936 / 1,234 / 1,129** | — |
| 3 | crosswalk 매핑 | 위 | 3,382 행 | tid → GenBank ID, **자기 균주 프로테옴에 실재** | metal 80 / atp 3,299 | `folddisco_{metal,atp}_motif_detail.csv` |
| 4 | 단백질 단위 | 위 | — | **RMSD 최저 매치 1건**만 | 80 | `result/table/motif_metal_3strain.csv` |
| 5 | Tier 분류 | 위 | 80 | 금속+ATP=1 / 금속+어노테이션=2 / 금속만=3 | 1:24 2:49 3:7 | `result/table/FOLDDISCO_candidates.csv` |
| 6 | rmsd 민감도 | `result/folddisco/sweep/` | 3 값 ×3 균주 | 1.0 / 1.5 / 2.0 | **전부 동일** — 필터가 아니었음 | `result/table/folddisco_rmsd_sweep.csv` |
| 7 | Cys4 재질의 | `result/folddisco/cys4/` | 3 템플릿 ×3 | `--rmsd 1.0` | 103·56·47 (BL21) | `result/table/folddisco_cys4_requery.csv` |
| 8 | 허용폭 확장 | `result/folddisco/widen/` | 3 조합 ×3 | `-d 1.5 -a 15` | 52 / 50 / 64 | — |
| 9 | 이량체 질의 | `result/folddisco/{dimer,covnode}/` | — | `--covered-node 4` | **0** (node_count 전부 2) | — |

**대조군:** Y19 의 CooC(`AHZ96930.1`) 를 **기본·확장 설정 모두에서 못 찾음** → **Track C 대조군 실패.**
→ metal 29/21/30 은 "CooC 같은 단백질"이 아니라 **"A112/A114 간격을 가진 단백질"** 목록.

---

## E. Foldseek (구조 전체 비교) — 유일하게 대조군 통과

| # | 질의 | 대상 | 컷오프 | 결과 | 결과 파일 |
|---|---|---|---|---|---|
| 1 | CooC1 (3kji) | Y19 구조 5,283 | `-e 10 --max-seqs 2000` | **AHZ96930.1 tm 0.938** ✅ | `result/foldseek/cooc1_vs_structures_UP000034085.m8` |
| 2 | CooC1 | BL21 구조 4,096 | 동일 | 최고 **ApbC tm 0.690** (CooC 없음) | `result/foldseek/cooc1_vs_structures_UP000503272.m8` |
| 3 | **QJZ12568.1** | 양쪽 | 동일 | **HypA tm 0.664 / 0.680** | `result/foldseek/QJZ12568.1_vs_*.m8` |
| 4 | QJZ11471.1 | 양쪽 | 동일 | 산발 (tm 0.33~0.57) — 일관성 없음 | `result/foldseek/QJZ11471.1_vs_*.m8` |
| 5 | QJZ13396.1 | 양쪽 | 동일 | GspE tm 0.880 — 진짜 GspE | `result/foldseek/QJZ13396.1_vs_*.m8` |
| 6 | coo 오페론 세트 | BL21·MG1655 | 동일 | **진행 중** | `result/table/coo_operon_homologs.csv` |
| 7 | HypA 폴드 | 3 균주 | **tm ≥ 0.50** | **진행 중** | `result/table/hypA_fold_sweep.csv` |

**판정 기준:** `alntmscore ≥ 0.5` = 같은 폴드.

---

## F. 균주 특이성 판정

| # | 단계 | 원본 | 건수 | 컷오프 | 남은 것 | 결과 파일 |
|---|---|---|---|---|---|---|
| 1 | 모티프 ↔ 모티프 | `result/search/fdm_motif_vs_motif.m8` | 222 | `fident ≥ 0.50` `qcov ≥ 0.70` | — | → `motif_pattern` |
| 2 | 모티프 → 프로테옴 | `result/search/fdm_motif_vs_{균주}.m8` | 187/190/187 | 동일 | — | → `gene_pattern` |
| 3 | BL21 판정 | 위 | 29 | 두 축 교차 | 공통 15 / 검출차이 12 / **특이 2** | `motif_metal_3strain.csv` |
| 4 | 검출차이 검증 | 직교체 서열 매핑 | 4 | 모티프 잔기가 Cys 인가 | **4개 전부 보존 → 인공물** | — |
| 5 | Cys4 × 균주 | `result/search/cys4_bl21_vs_{균주}.m8` | 153 | `qcov ≥ 0.70`, fident 구간 | specific 7 / low 3 = **10** | `result/table/cys4_bl21_strain_specificity.csv` |
| 6 | 게놈 확인 | `result/search/cys4_specific_vs_MG1655genome.m8` | 7 | `--search-type 2`, **fident ≥ 0.80** 이어야 같은 유전자 | **게놈에도 없음 2** | `result/table/cys4_genome_check.csv` |

**최종 통과 2개:** `QJZ12568.1` (133 aa, HypA 폴드, Cys4) · `QJZ11471.1` (93 aa, 폴드 불명)

---

## G. 폐기·미사용 (기록만)

| 지표 | 파일 | 폐기 사유 |
|---|---|---|
| `ligand_iptm` | `trackB_boltz2_ranked.csv` 열 | 계면 배치 0/348, 길이와 ρ=−0.807 |
| `idf` (folddisco) | detail CSV 열 | idf 1위가 44잔기·pLDDT 56.7 — 희귀 ≠ 기능 |
| ATP 모티프 단독 | `motif_atp_3strain.csv` | 프로테옴의 23~26% 보유 |
| RF2-PPI 순위 | `trackA_RF2PPI_ranked.csv` | 대조군 0.266, 531/3,625 아래 |
| Boltz ipTM 순위 | `metal_motif_3strain_boltz.csv` | 대조군 0.300, 직교체 잡음 0.102(최대 0.309) |
| 이량체 folddisco 질의 | `folddisco/covnode/` | AlphaFold DB 가 단량체라 계면 질의 불가 |
| CELL 32 통합 랭킹 | `INTEGRATED_candidate_ranking.csv` | 가중치의 35%가 위 두 실패 지표 |

---

## H. 최종 산출물

| 파일 | 내용 |
|---|---|
| `result/table/FINAL_candidates.csv` | 등급표. **제외 없음**, `reason` 열에 가감점 내역 |
| `result/table/ALL_motif_candidates_full.csv` | 세 균주 80개 전체 + 모든 증거 열 |
| `result/table/ni_site_grade.csv` | Ni 배위 도너 등급 A~D |
| `result/table/cys4_genome_check.csv` | 게놈 수준 균주 특이성 |
| `result/ni_view/` | PyMOL 세션 (`view.pml`, `coordination.txt`) |
| `result/log/*.log` | 전 백그라운드 작업 로그 |
| `script/*.sh` | 생성된 실행 스크립트 |
