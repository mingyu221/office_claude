# 결과 파일 위치 — 발표자료 참조용

`BASE = /mnt/af2results/mingyu/workspace/ppi_discovery`
`DB   = /mnt/af2results/mingyu/database/bacteriaDB`

---

## 1. Folddisco — 원본 TSV

```
BASE/result/folddisco/
├── metal_run01_BL21.tsv              원래 질의 A112,A114        29
├── metal_run01_Y19.tsv                                          30
├── metal_run01_MG1655.tsv                                       21
├── atp_run02_{BL21,Y19,MG1655}.tsv   Walker A A12,A13,A14   936/1234/1129
│
├── sweep/        metal_r{10,15,20}_{균주}.tsv          rmsd 1.0/1.5/2.0
├── cys4/         cys4_{템플릿}_{균주}.tsv              템플릿 = QJZ13803.1·QJZ12568.1·QJZ13396.1
├── dimer/        {mono2,dimer3,dimer4}_{균주}.tsv      이량체 질의
├── widen/        mono2_d{d}_a{a}_{균주}.tsv            -d/-a 확장
└── covnode/      dimer4_cov{0,2,3,4}_{균주}.tsv        covered-node 기울기
```

**열:** `tid, idf, min_rmsd, nres, plddt, node_count, max_node_cov, matching_residues, …`

**`matching_residues` 형식:** `A27,A30:0.5054;A50,A52:0.81`
→ 콜론 뒤가 **그 매치의 RMSD**, 세미콜론이 **매치 구분**.

**질의 파라미터 (기본):** `-d 0.5 -a 5.0 --rmsd 1.0 --top 2000 --per-structure --sort-by idf`
※ `--rmsd` 는 기본이 *no limit* 이고, 실측 최대가 0.745 라 **실제로는 아무것도 자르지 않았다.**
  병목은 `-d`(거리 허용폭) / `-a`(각도 허용폭) 였다.

---

## 2. Folddisco — 가공된 CSV

```
BASE/result/table/
├── folddisco_metal_motif.csv          단백질당 1행 (crosswalk 매핑된 것만)
├── folddisco_metal_motif_detail.csv  ★ 전체 행 (MG1655 미매핑 포함)
├── folddisco_atp_motif{,_detail}.csv
├── motif_metal_3strain.csv          ★ 3균주 80개 + 두 축 균주 판정
├── motif_atp_3strain.csv
├── FOLDDISCO_candidates.csv           Tier 1/2/3
├── folddisco_rmsd_sweep.csv
├── folddisco_cys4_requery.csv
└── folddisco_dimer_query.csv
```

---

## 3. Foldseek — 결과

```
BASE/result/foldseek/
├── cooc1_vs_structures_UP000503272.m8      ★ Ch CooC1 → BL21
├── cooc1_vs_structures_UP000034085.m8      ★ Ch CooC1 → Y19   (AHZ96930.1 tm 0.938)
├── QJZ12568.1_vs_structures_UP000{503272,034085}.m8   ★ HypA 폴드 발견
├── QJZ11471.1_vs_*.m8
├── QJZ13396.1_vs_*.m8
├── coo_{단백질}_vs_{구조DB}.m8             CELL 50 STEP4 — Y19 coo 오페론 질의
├── hypa_{단백질}_vs_{구조DB}.m8            CELL 50 STEP5 — HypA 폴드 스윕
├── fdall_vs_{구조DB}.m8                    CELL 51 — Folddisco 히트 전체의 폴드
├── g3e_vs_{구조DB}.m8                      CELL 52 — G3E/SIMIBI 과 스윕
└── db_structures_UP000{503272,034085}*     구조 DB 캐시 (createdb, 재사용)
```

**질의 구조:** `/mnt/af2results/mingyu/workspace/seek_ni_insertase/input/3kji.pdb`
※ `Ignore 4 out of 6` — 3kji 에 짧은 체인 4개(리간드/펩타이드)가 있어 **A·B 두 사슬만** 질의로 쓰였다.

**열 (헤더 없음, `--format-output` 순서대로):**

| CELL 29k | `query target fident alnlen qcov tcov evalue bits prob alntmscore lddt` |
|---|---|
| CELL 51·52 | `query target fident alnlen qcov tcov evalue bits alntmscore lddt` (prob 없음) |

**`target` 은 `AF-XXXXXX-F1-model_v6.cif` 형식이라 단백질 ID 가 아니다.**
→ `BASE/result/table/id_crosswalk_struct_to_genbank.csv` 로 변환해야 QJZ*/AKE* 가 나온다.

**TM-score 기준:** `< 0.17` 무관 · **`≥ 0.5` 같은 폴드(표준 컷오프)** · `≥ 0.9` 사실상 같은 구조
※ 근사 계산(`--tmalign-fast`)에서 1.00~1.01 이 나올 수 있다. `--exact-tmscore 1` 로 해소.

---

## 4. 구조 디렉터리 (검색 대상)

```
DB/structures_UP000503272/    BL21   4,096 구조   AlphaFold model_v6
DB/structures_UP000034085/    Y19    5,283 구조
MG1655                        folddisco 인덱스만 (공식 배포본, model_v4)
```

**★ 버전 불일치 주의:** BL21 **v6** / MG1655 **v4**.
검출 차이 4개(ErpA·PgaC·YeeO·YeiR)가 **서열 100% 동일한데 한쪽만 검출**된 원인이 이것이다.

---

## 5. 교차·최종 산출물

```
BASE/result/table/
├── folddisco_x_foldseek.csv       ★ 모티프 × 폴드 교차 (둘 다 / Foldseek만 / Folddisco만)
├── g3e_family_sweep.csv           ★ G3E·SIMIBI 과 구성원
├── ni_site_grade.csv              ★ Ni 배위 도너 등급 A~D
├── boltz_ni_placement.csv           411 구조의 Ni 좌표
├── cys4_bl21_strain_specificity.csv 프로테옴 수준 균주 특이성
├── cys4_genome_check.csv          ★ 게놈 수준 확인
├── coo_operon_homologs.csv          Y19 오페론 상동체
├── hypA_fold_sweep.csv              HypA 폴드 구성원
├── ortholog_reproducibility.csv     직교체 재현성 (잡음 바닥)
├── boltz_replicate_summary.csv      5표본 + 양성대조군
├── FINAL_candidates.csv           ★ 등급표 (제외 없음)
└── ALL_motif_candidates_full.csv    80개 전체 + 모든 증거 열
```

```
BASE/result/ni_view/       PyMOL 세션 (view.pml, coordination.txt)
BASE/result/log/*.log      백그라운드 작업 로그
BASE/script/*.sh           생성된 실행 스크립트
```

---

## 6. 발표에 쓸 핵심 숫자

| 항목 | 값 | 출처 파일 |
|---|---|---|
| Track A 양성대조군 | **0.266** (기준 0.74), 531/3,625 아래 | `trackA_RF2PPI_ranked.csv` |
| Track B 양성대조군 | **ipTM 0.300** (5표본), 후보 8개가 전부 위 | `boltz_replicate_summary.csv` |
| Track C 양성대조군 | Y19 CooC 를 **못 찾음** | `folddisco/metal_run01_Y19.tsv` |
| **Foldseek 대조군** | **Y19 CooC tm 0.938** ✅ | `foldseek/cooc1_vs_structures_UP000034085.m8` |
| ipTM 잡음 (직교체 23쌍) | Δ 중앙값 **0.102**, 최대 0.309 (전체 폭 0.485) | `ortholog_reproducibility.csv` |
| metal_rmsd 재현성 | Δ 중앙값 **0.015 Å** | 동상 |
| Ni 배위 거리 | 평균 **2.24 Å** (1.97~3.27) | `boltz_ni_placement.csv` |
| ATP 모티프 보유율 | **23~26%** — 선별력 없음 | `motif_atp_3strain.csv` |
| BL21 금속모티프 중 MG1655 공유 | **27 / 29** | `motif_metal_3strain.csv` |
| 검출 차이 4개 | **fident 1.000 인데 한쪽만** → 구조예측 인공물 | CELL 36c 출력 |
| Cys4 × 균주특이 | 153 → **10** → 게놈 확인 후 **2** | `cys4_genome_check.csv` |

**1순위 후보 `QJZ12568.1`** (133 aa, hypothetical HO396_09830)
표 전체에서 **유일하게 네 가지를 모두 만족**:
Cys4 재질의 히트 · HypA 폴드(tm 0.680) · Ni 배위 **grade A**(2.23–2.40 Å) · MG1655 프로테옴·게놈 **둘 다 없음**
