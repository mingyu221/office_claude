# 2026-09-17 작업 내역 — 인계문서 이후

`BASE = /mnt/af2results/mingyu/workspace/ppi_discovery` · `TBL = BASE/result/table`

---

## 1부. Track A 판정 (RF2-PPI)

| # | 셀 | 무엇을 / 왜 바꿨나 | 결과 | 산출 파일 |
|---|---|---|---|---|
| 1 | `15b` | CELL 15 는 로그 꼬리만 봐서 tqdm 진행률이 안 잡혔다. `\r`→`\n` 변환 후 파싱 | s/pair 로 GPU 오배치를 분 단위에 감지 | — |
| 2 | `23`·`24` | 3,626쌍 ×3 집계 + 양성대조군 | **대조군 0.266** (기준 0.74), `mean≥0.74` 0개 | `TBL/trackA_RF2PPI_{all_pairs,ranked}.csv` |
| 3 | `24b` | CELL 24 가 "깊이 부족"을 **고정 문자열**로 출력했다. 실제 깊이를 읽게 고침 | 깊이 **1,451** = 스크리닝 75분위(1,392) 위 → 깊이 탓 아님. 대조군 위에 **531/3,625** | — |
| 4 | `24c` | 대조군만 **체인 순서가 반대**였고 절편은 미사용(`SEGMENT_BAITS=[]`). 기존 a3m 을 잘라 8변형 생성 | 최고 **0.251**. 절편은 오히려 **0.039~0.114** 로 악화 | `TBL/trackA_control_variants.csv` |
| 5 | `24d`·`50` | ★ **로직 변경**: 지금까지 전부 단량체×단량체였다. 문헌상 CooC1 Ni 자리는 **이량체 계면**(각 monomer Cys112/114) | **진행 중** | `BASE/result/paired/dimer/in_t*.log` |

---

## 2부. Boltz 눈금 검증

| # | 셀 | 무엇을 / 왜 | 결과 | 파일 |
|---|---|---|---|---|
| 6 | `28a-11` | ipTM 잡음을 몰랐다. 세 균주 **직교체**는 입력이 같으니 차이가 곧 잡음 | `metal_rmsd` Δ **0.015** / `ipTM` Δ **0.102**(최대 0.309) | `TBL/ortholog_reproducibility.csv` |
| 7 | `28a-12`·`13` | Boltz 쪽엔 대조군이 **없었다**. CooC1–ChCODH2 를 5표본으로 | **ipTM 0.300** — 후보 8개가 전부 위 → **2번째 대조군 실패** | `TBL/boltz_replicate_summary.csv` |
| 8 | `28a-7` 수정 | 폴더를 `["out","out_focus"]` 로 박아놔 `out_cmp` 를 안 읽었다 → `out*` 전부 훑기 | 348 → **411 구조**. MG1655·Y19 Ni 배치 확보 | `TBL/boltz_ni_placement.csv` |
| 9 | `28a-8b` | "붙었다/안 붙었다"만 봤다. **도너 종류로 등급** | **A(Cys4) 4** · B 9 · C 1 · D 8. **YeiR 가 D 로 강등**(접촉 1개) | `TBL/ni_site_grade.csv` |

---

## 3부. Folddisco 노선 재구축

| # | 셀 | 무엇을 / 왜 | 결과 | 파일 |
|---|---|---|---|---|
| 10 | `33` | ID 계열이 균주마다 달랐다. **대응표를 둘로 분리** — 모티프↔모티프 / 모티프→프로테옴 | '검출 차이' 와 '유전자 부재' 를 구분 가능 | `search/fdm_motif_vs_*.m8` |
| 11 | `34` | 두 축 교차 판정 | BL21 29 중 **공통 15 / 검출차이 12 / 특이 2** | `TBL/motif_metal_3strain.csv` |
| 12 | `35` | ATP 는 필터가 아니라 표시로 | **23~26%** 보유 → 선별력 없음 확정 | `TBL/motif_atp_3strain.csv` |
| 13 | `36` | Tier 1/2/3, 삭제 없이 우선도만 | 1:24 2:49 3:7 | `TBL/FOLDDISCO_candidates.csv` |
| 14 | `36b` | BL21 만 찍고 PPI 점수를 뺐었다 → **80개 전부 + 모든 열**, 열마다 신뢰 등급 | 거르지 않음 | `TBL/ALL_motif_candidates_full.csv` |
| 15 | `36c` | 검출 차이 4개의 **Cys 보존 검사** | ErpA·PgaC·YeeO **fident 1.000**, Cys 보존 → **전부 구조예측 인공물**(v6 vs v4) | — |

---

## 4부. Folddisco 질의 자체를 의심

| # | 셀 | 무엇을 / 왜 | 결과 | 파일 |
|---|---|---|---|---|
| 16 | `29b` | 컷오프 1.0 Å 가 맞는지 검증한 적 없음 → 1.0/1.5/2.0 | **히트 불변**, 최대 min_rmsd 0.745 → `--rmsd` 는 병목이 아니었다 | `TBL/folddisco_rmsd_sweep.csv` |
| 17 | `29c` | ★ **로직 변경**: 질의를 바꿔야 한다. 실제 Ni 잡은 **Cys4 기하**를 템플릿으로 | 103·56·47 (2~4배 증가). 단 **HypA·HypB·Nik 전부 놓침** 확인 | `TBL/folddisco_cys4_requery.csv` |
| 18 | `29d` | Cys4 × 균주 특이 교차 | **버그** — `classify_MG1655_vs_BL21DE3.csv`(4,300행=MG1655 기준)를 집어 교집합 0. 내용 기반 선택으로 수정 | — |
| 19 | `29e` | 분류표 대신 **mmseqs 로 직접** 물음 | 153 중 specific 7 + low 3 = **10개** | `TBL/cys4_bl21_strain_specificity.csv` |
| 20 | `29f`·`29h`·`29i` | ★ 당신 지적 — CooC1 은 호모다이머. `A112,A114,B112,B114` | `--covered-node 4` 에서 **0**. `node_count` 전부 2 → **AlphaFold DB 가 단량체라 계면 질의 불가**. 노선 폐기 | `folddisco/{dimer,covnode,widen}/` |
| 21 | `29g` | 프로테옴 "없음"은 3가지를 뭉침 → **게놈**(`--search-type 2`)으로 재확인. 이후 pseudogene 기준을 `fident≥0.80` 으로 좁힘 | 게놈에도 없음 **2개** | `TBL/cys4_genome_check.csv` |

---

## 5부. Track C 대조군 → Foldseek 전환

| # | 셀 | 무엇을 / 왜 | 결과 | 파일 |
|---|---|---|---|---|
| 22 | `29j` | ★ 당신 지적 — **Y19 의 CooC 를 못 찾으면 대조군 실패다** | Y19 CooC(`AHZ96930.1`)를 기본·확장 **모두에서 못 찾음** → **3번째 대조군 실패** | — |
| 23 | `29k` | ★ **도구 전환**: 모티프(Folddisco) → **구조 전체(Foldseek)** | **CooC1 → Y19 CooC tm 0.938** ✅ 대조군 통과. BL21 최고는 ApbC 0.690 (CooC 없음) | `result/foldseek/cooc1_vs_*.m8` |
| 24 | `29k`(B) | 1순위 후보의 정체 규명 | **`QJZ12568.1` = HypA 폴드** (tm 0.664/0.680) | `result/foldseek/QJZ12568.1_vs_*.m8` |

---

## 6부. 정리 · 다음

| # | 셀 | 내용 | 파일 |
|---|---|---|---|
| 25 | `39` | ★ 당신 지시 — **필터링 대신 등급제**. 제외 0, `reason` 열에 가감점 | `TBL/FINAL_candidates.csv` |
| 26 | `40`·`41` | 논문 Fig.2 — Y19 는 CooC·CooT·CooJ·HypB 보유. **진짜 대조군 세트** | `TBL/coo_operon_homologs.csv` |
| 27 | `42`(50 STEP5) | HypA 폴드 전수 조사 | `TBL/hypA_fold_sweep.csv` |
| 28 | `50` | 위 전부를 한 셀로. GPU 둘은 `part8_gpu` 로 순차 백그라운드 | `script/part8_gpu.sh` · `log/part8_gpu.log` |

---

## 오늘 로직이 바뀐 지점 (되짚을 순서)

1. **깊이 탓 → 모델 탓** (#3) — 대조군 깊이가 오히려 상위 25% 였다
2. **`--rmsd` → `-d`/`-a`** (#16) — 필터가 아니라 검색 단계가 병목이었다
3. **CooC1 2잔기 → Cys4 4잔기** (#17) — 넓어졌지만 **잔기가 달라 부분집합이 아니었다**
4. **단량체 → 이량체 질의** (#20) — 논리는 맞았으나 **DB 가 단량체**라 불가
5. **프로테옴 → 게놈** (#21) — pseudogene 과 어노테이션 누락을 갈라야 했다
6. **Folddisco → Foldseek** (#23) — 대조군을 통과하는 유일한 도구
7. **필터 → 등급** (#25) — 거른 게 답이면 되돌릴 수 없다
8. **단량체 → 이량체 입력** (#5) — Track A 를 닫지 않고 전략 변경
