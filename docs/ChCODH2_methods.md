# ChCODH2 Ni insertase 탐색 — 수행한 분석 (방법)

무엇을 어떤 도구로 어떤 설정으로 돌렸고, 그래서 몇 개가 남았는가.
해석과 결론은 `ChCODH2_Ni_insertase_report.md`, 파일 위치는 `ChCODH2_file_locations.md`.

작업 루트 `BASE = /mnt/af2results/mingyu/workspace/ppi_discovery`

---

## 0. 입력

| 항목 | 경로 | 건수 |
|---|---|---|
| BL21(DE3) 프로테옴 | `bl21_db_match_qjz.faa` | 4,110 |
| MG1655 프로테옴 | `mg1655_protein.faa` | 4,300 |
| Y19 프로테옴 | `y19_db_match.faa` | 5,325 |
| MG1655 게놈 | `MG1655_U00096.3.fna` | 4.7 MB |
| BL21 구조 (AlphaFold) | `structures_UP000503272/` | 4,096 |
| Y19 구조 (AlphaFold) | `structures_UP000034085/` | 5,283 |
| MG1655 구조 | folddisco 공식 인덱스만 | — |
| CooC1 결정구조 | `3kji.pdb` | 2 사슬, 112·114 = CYS |
| bait 서열 | `baits.fasta` | 8 (ChCODH2_WT + CooC1 + 절편 6) |

구조 DB 버전이 갈린다 — BL21·Y19 는 **model_v6**, MG1655 는 공식 배포 **v4**.
뒤에서 "BL21 에서만 검출" 로 보인 것들의 원인이 이것이다.

---

## 1. 균주 간 서열 분류 → Track 배정 (CELL 09–12)

BL21 단백질 하나하나가 MG1655 에 있는지 없는지를 먼저 가른다.

```bash
mmseqs easy-search bl21.faa mg1655.faa BL21_vs_MG1655.m8 tmp \
  --format-output "query,target,fident,alnlen,evalue,bits,qlen,tlen,qcov,tcov" \
  -s 7.5 -e 1e-5 --max-seqs 5
# 역방향(MG1655 → BL21)도 같은 설정으로 한 번 더
```

`-s 7.5` 는 mmseqs 최고 민감도. 부재 판정을 하려는 것이므로 놓치는 쪽이 치명적이다.
양방향을 다 돌리는 이유는 "BL21 에 있고 MG1655 에 없는 것" 과
"MG1655 에 있고 BL21 에 없는 것(억제자 후보)" 이 다른 질문이기 때문이다.

최고 bits 히트 1건만 남기고 4분류:

| 분류 | 기준 | 뜻 |
|---|---|---|
| identical | pident ≥ 99.9 **and** cov ≥ 0.99 | 같은 단백질 |
| high-similarity | pident ≥ 90 | 직교체 |
| low-similarity | pident < 90 | 발산한 직교체 |
| strain-specific | 히트 없음 | BL21 에만 |

`cov = min(qcov, tcov)` — 한쪽만 긴 부분정렬을 동일로 올리지 않기 위해서다.

**결과 — Track 배정**

| Track | 구성 | 수 | 왜 이 방법인가 |
|---|---|---|---|
| **A** | identical + high-similarity | **3,746** | 직교체가 많으니 **공진화 신호**를 쓸 수 있다 → RF2-PPI |
| **B** | low-similarity + strain-specific | **326** | 직교체가 없어 MSA 가 안 쌓인다 → 구조 예측으로 직접 co-fold |

→ `result/table/classify_*.csv`

---

## 2. Track A — RF2-PPI (공진화 기반 상호작용 예측)

### 2-1. 상동체 수집

```bash
mmseqs search baitDB bactDB bait_res tmp \
  -s 7.5 --num-iterations 3 -e 1e-3 --max-seqs 20000
mmseqs search preyDB bactDB prey_res tmp \
  -s 7.5 --num-iterations 3 -e 1e-3 --max-seqs 20000 --db-load-mode 2
```

`--num-iterations 3` = 반복검색(profile search). 먼 상동체까지 긁어야 MSA 깊이가 나온다.
bait 히트 43,036 행 → taxid 당 최고 1건 → **2,929 유전체**.
prey 히트 **62,642,570 행**.

### 2-2. paired MSA

같은 유전체에서 온 bait 서열과 prey 서열을 **이어 붙여** 한 줄로 만든다.
같이 진화했으면 두 줄의 변이가 같이 움직인다 — 그것이 RF2-PPI 가 보는 신호다.

- 공유 유전체 **≥ 50** 인 쌍만 (그 아래는 신호가 안 선다)
- `hhfilter -id 90 -M first` — 90% 이상 같은 서열을 접어 중복 가중 제거

3,730쌍 시도 → **3,626쌍** 생성.

### 2-3. 추론

```bash
CUDA_VISIBLE_DEVICES=N python predict_list_PPI.py \
  -list_fn input_repN -model_file RF2-PPI.pt
```

비결정적이라 **replicate 3회**. 3,627 × 3 = 10,881 예측 → prey 단위로 bait 합집합 최대값.

### 2-4. 대조군 — 여기서 멈췄다

CooC1–ChCODH2_WT (MSA 깊이 1,451) 를 넣었다. **알려진 상호작용이므로 높게 나와야 한다.**

| | 값 | 기준 |
|---|---|---|
| CooC1–ChCODH2 양성대조군 | **0.266** | ≥ 0.70 |
| 변형 6종 (절편·순서·이량체) | 최고 **0.251** | ≥ 0.70 |
| 전체 3,625 중 판정 컷 통과 | **0개** | mean ≥ 0.74 |

대조군이 실패했으므로 **Track A 는 후보를 내지 않는다.** 3,746개가 평가 없이 남았다.

이량체 재시도(ChCODH2 를 2개 이어 붙인 tandem)를 진행 중이다 — CODH 가 생리적으로
이량체이고 CooC1 도 이량체라, 단량체 쌍으로 물은 것 자체가 틀린 질문일 수 있다.

---

## 3. Track B — Boltz-2 co-folding + Ni 배위 판독

### 3-1. 예측

ChCODH2 와 후보 단백질, 그리고 **Ni 이온(CCD: NI)** 을 한 번에 넣고 접는다.

```bash
CUDA_VISIBLE_DEVICES=N boltz predict inputs --out_dir out \
  --use_msa_server --max_msa_seqs 2048 \
  --recycling_steps 3 --diffusion_samples 1 \
  --output_format mmcif --num_workers 2
```

326 → 길이 제한(bait+prey ≤ 1,830 잔기)으로 323 → OOM 2건 손실 → **321 완료**.
추가로 focus 27 · 타 균주 비교 50 · 반복 13(65 구조) = 총 **411 구조**.

### 3-2. 무엇을 읽었나 — 두 번 갈아탔다

| 지표 | 결과 | 판정 |
|---|---|---|
| `ipTM` | 대조군 CooC1 = **0.300** 인데 후보 8개가 전부 그 위 | 눈금이 안 맞음. **순위에 안 씀** |
| `ligand_ipTM` | 계면 배치 **0/348**, 길이와 ρ = **−0.807** | **폐기** |
| **Ni 배위 좌표** | 직교체 간 Δ **0.01–0.08 Å** | **재현됨 → 이것만 씀** |

ipTM 은 모델이 자기 계면을 자평한 값이다. 대조군이 0.300 인데 후보가 그보다 높다는 건
그 값이 "진짜 상호작용" 이 아니라 다른 것을 재고 있다는 뜻이다.
반면 **Ni 이 실제로 어디에 놓였는지**는 직교체끼리 0.01 Å 수준으로 일치했다.

### 3-3. 배위 등급

Ni 반경 **2.6 Å** 안의 원자만 배위로 세고, **도너 종류**로 등급을 매겼다.
거리만 보면 표면에 스친 것과 자리에 박힌 것이 구분되지 않는다.

| 등급 | 기준 | 수 |
|---|---|---|
| **A** | Cys 티올레이트 4개 — Ni(II) 같은 연한 금속의 전달·저장형 자리 | **4** |
| B | 자리는 좋으나 도너 구성이 섞임 | 9 |
| C | — | 1 |
| D | 산소 공여체 위주 / 접촉 1–2개 — 표면 부착 의심 | 8 |

**grade A 4개:** `QJZ13803.1` HslO 2.18 Å · `QJZ12568.1` 2.23 Å ·
`QJZ13396.1` GspE 2.29 Å · `Y19-AKE60400.1` (HslO 직교체) 2.17 Å

→ `result/table/ni_site_grade.csv`

---

## 4. Track C — Folddisco 구조 모티프 검색

CooC1 의 Ni 배위 Cys 쌍 기하를 질의로, 세 균주 구조 DB 전체를 훑는다.

```bash
folddisco query -p 3kji.pdb -q A112,A114 -i {strain_index} \
  -t 8 --per-structure --header --sort-by idf --rmsd 1.0 --top 2000
```

`-d` 0.5 Å / `-a` 5.0° 는 기본값. 히트 **BL21 29 / Y19 30 / MG1655 21**.
ATP(Walker A) 질의 `A12,A13,A14` 도 같은 설정으로 → 936 / 1,234 / 1,129.

crosswalk 로 구조 ID → GenBank ID 매핑 후 자기 균주 프로테옴에 실재하는 것만 →
금속 **80** / ATP 3,299. 단백질 단위로 RMSD 최저 매치 1건만 남김.

### 4-1. 확인한 것

- **`--rmsd` 는 병목이 아니었다** — 1.0 / 1.5 / 2.0 스윕 결과 **전부 동일**.
  folddisco 의 `--rmsd` 기본값은 제한 없음이라 애초에 걸린 적이 없다
- **`--covered-node` 를 안 주고 있었다** — 기본 0 이라 4잔기 질의에 2잔기만 맞아도
  히트로 잡힌다. `--covered-node 4` 로 재실행

### 4-2. 대조군 — 실패

Y19 에는 **CooC 직교체가 실재한다**(`AHZ96930.1`, Foldseek tm 0.938).
CooC1 질의는 자기 직교체를 찾아야 한다.

| 설정 | BL21 | MG1655 | Y19 | Y19 CooC |
|---|---|---|---|---|
| d 0.5 / a 5.0 | 29 | 21 | 30 | **X** |
| d 1.0 / a 10 | 49 | 48 | 61 | **X** |
| d 1.5 / a 15 | 52 | 50 | 64 | **X** |
| d 2.0 / a 20 | 48 | 46 | 60 | **X** |
| d 3.0 / a 30 | 48 | 46 | 60 | **X** |

**허용폭을 6배로 넓혀도 안 잡힌다.** 원인은 별도 노트북에서 실측했다 —
Y19 CooC 의 CXC 는 제자리에 있고(A115/A117) 자기 질의로 rmsd 0.000 에 잡히는데,
**SG–SG 가 8.16 Å** 이다. Ch CooC1 은 **3.58 Å**.
AlphaFold 가 그 자리를 **열린(apo) 상태**로 예측했다.

> **질의는 닫힌 결정구조, 대상은 전부 열린 예측구조.** 허용폭으로 풀 문제가 아니다.

**따라서 folddisco 의 "0" 을 부재로 읽을 수 없다.** 히트 80개 목록 자체는
"AF 모델에서 그 Cys 쌍이 닫혀 예측된 단백질" 목록으로는 유효하다.

---

## 5. Foldseek — 대조군을 통과한 유일한 축

모티프가 아니라 **구조 전체**를 정렬한다.

```bash
foldseek easy-search query.pdb {strain_structures} out.m8 tmp \
  -e 10 --max-seqs 2000 --exact-tmscore 1
```

판정 기준 `alntmscore ≥ 0.5` = 같은 폴드 (Xu & Zhang 2010).

| 질의 | 결과 |
|---|---|
| CooC1 → Y19 구조 5,283 | **`AHZ96930.1` tm 0.938** ✅ 대조군 통과 |
| CooC1 → BL21 구조 4,096 | 최고 ApbC **tm 0.690** — **BL21 에 CooC 가 없다** |
| `QJZ12568.1` → 양쪽 | **HypA 폴드 tm 0.664 / 0.680** |
| `QJZ13396.1` → 양쪽 | GspE tm 0.880 (진짜 GspE) |

---

## 6. 균주 특이성 확정 — 프로테옴 다음에 게놈

"프로테옴에 없다" 는 **진짜 없음 / 어노테이션 누락 / pseudogene** 을 뭉친 말이다.
게놈을 직접 봐야 갈린다.

```bash
mmseqs easy-search candidates.faa MG1655_U00096.3.fna out.m8 tmp \
  --search-type 2 -s 7.5 -e 1e-3
```

`--search-type 2` = 단백질 질의 vs 번역된 뉴클레오타이드. 어노테이션을 우회한다.
**주의:** 이 모드의 `alnlen` 은 염기 단위다. 커버리지는 `qstart/qend` 로 계산해야 한다
(처음에 3배 부풀려 읽었다).

| 단계 | 기준 | 남은 것 |
|---|---|---|
| Cys4 보유 BL21 단백질 | — | 153 |
| MG1655 대비 분류 | `qcov ≥ 0.70`, fident 구간 | specific 7 + low 3 = **10** |
| 게놈 확인 | `--search-type 2`, **fident ≥ 0.80** 이어야 같은 유전자 | **2** |

**통과 2개:** `QJZ12568.1` (133 aa) · `QJZ11471.1` (93 aa)

### 검증 — "BL21 에서만 검출" 은 대부분 인공물이었다

folddisco 히트 목록에서 BL21 에만 잡힌 14개를 MG1655 에 직접 물었다.
정렬을 잔기 단위로 따라가 **모티프 Cys 가 상대에서도 Cys 인지**까지 확인했다.

**진짜 부재 0개.** 12개는 두 Cys 가 모두 보존된 직교체가 있었고(6개는 fident **1.000**),
1개는 게놈에 온전했다. 구조 DB 버전 차이(v6 vs v4)가 만든 검출 차이였다.

---

## 7. 후보가 나온 지점

**두 축이 겹치는 곳.** 둘 다 대조군·재현성이 확인된 측정이다.

| | Ni 배위 (Boltz, 재현 Δ 0.01–0.08 Å) | 균주 특이 (mmseqs 프로테옴+게놈) |
|---|---|---|
| **QJZ12568.1** hypothetical HO396_09830, 133 aa | **grade A** 2.23–2.40 Å | **MG1655 프로테옴·게놈 둘 다 없음** |
| **QJZ13396.1** GspE, 497 aa | **grade A** 2.29–2.47 Å | fident 0.600 (저유사도) |
| QJZ13803.1 HslO | grade A 2.18 Å | 세 균주 동일 — 변별력 없음 |
| QJZ11471.1 | — | 게놈에도 없음 (폴드 불명) |

1순위 `QJZ12568.1` 의 근거는 **서로 독립인 세 측정**이다.

| 축 | 측정 | 값 |
|---|---|---|
| 1 | MG1655 프로테옴 · 게놈 부재 | mmseqs + `--search-type 2` |
| 2 | Ni 배위 등급 | Cys4, 2.23–2.40 Å |
| 3 | 폴드 정체 | **HypA** tm 0.664 / 0.680 |

HypA 는 **[NiFe] 하이드로게나제의 Ni 삽입 인자**다. 기능 미상 + 133 aa 소형 +
Cys4 + HypA 폴드 = 메탈로샤페론 프로파일.

> **Cys4 재질의(folddisco) 히트는 근거로 세지 않는다.** 그 질의의 템플릿이
> 이 단백질 자신이라 자기참조다.

---

## 8. 지금 알고 있는 한계

1. **Track A 3,746개가 평가되지 않았다.** 대조군이 떨어져 RF2-PPI 가 멈췄고 대체 평가가 없다.
   BL21·MG1655 가 **둘 다 가졌는데 발현량이 다른** 인서테이즈라면 — 보고서 자체가
   가장 유력하다고 본 시나리오인데 — 지금 파이프라인이 구조적으로 못 잡는다
2. **folddisco 축의 비검출은 음성이 아니다** (5장). 질의/대상 상태 불일치
3. **`QJZ14223.1` LysR** — 최고 MG1655 히트 fident 0.345. 직교체가 아니라 먼 파라로그를
   잡은 것이라 따로 봐야 한다
4. 진행 중 — Track A 이량체 재시도, Cys4 `--covered-node 4` 재집계

---

## 부록. 대조군을 어디에 걸었나

방법마다 "답을 아는 입력" 을 하나씩 넣어 눈금을 쟀다. **이것이 이 분석의 척추다.**

| 방법 | 대조군 | 기대 | 실측 | 판정 |
|---|---|---|---|---|
| RF2-PPI | CooC1–ChCODH2 | ≥ 0.70 | **0.266** | ✗ |
| Boltz ipTM | CooC1 | 후보보다 높게 | **0.300** (후보가 더 높음) | ✗ |
| Boltz Ni 좌표 | 직교체 간 일치 | 작은 Δ | **0.01–0.08 Å** | ✅ |
| Folddisco | CooC1 질의 → Y19 CooC | 잡혀야 함 | **못 잡음** (d 3.0 까지) | ✗ |
| Foldseek | CooC1 → Y19 CooC | tm ≥ 0.5 | **tm 0.938** | ✅ |
| mmseqs 부재 판정 | fident 1.000 직교체 14건 | 부재로 안 잡혀야 함 | **진짜 부재 0** | ✅ |

다섯 중 셋이 떨어졌고, 후보는 통과한 셋에서만 나왔다.
