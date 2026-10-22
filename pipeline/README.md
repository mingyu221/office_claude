# ChCODH2 상동체 분석 파이프라인 — IPR010047

ChCODH2 (*Carboxydothermus hydrogenoformans* CO dehydrogenase II, Q9F8A8, 636 aa, 호모다이머,
금속 클러스터 B/B′/C/C′/D) 를 inverse folding 으로 재설계하기 위해, IPR010047 상동체의
클러스터 구조와 잔기 보존도를 정량화한다.

최종 목표는 **636행짜리 위치 테이블 하나**다. 행은 ChCODH2 잔기 1~636, 열은 보존도·구조 역할·
설계 판정(Fix / Bias / Free). 이 테이블이 그대로 LigandMPNN 입력 스펙으로 바뀐다.

---

## 두 개의 축

혼동하기 쉬운 두 가지를 분리해 둔다. 단계마다 어느 축에 열을 더하는지가 정해져 있다.

| | **축 A — 서열** | **축 B — 위치** |
|---|---|---|
| 단위 | 1,827개 상동체 (accession) | 636개 잔기 (ChCODH2 번호) |
| 기준표 | `master_1827.tsv` (기존) | 단계 6에서 생성 |
| 담는 것 | 클러스터, fident, 분류군, QC 플래그 | 보존도, 구조 역할, Fix/Bias |

두 축을 잇는 것이 **단계 3의 `equiv_map`** 이다. 파이프라인 전체의 척추이고,
여기가 틀리면 뒤의 모든 숫자가 조용히 틀린다.

---

## 단계

| | 단계 | 입력 | 산출 | 축 | 상태 |
|---|---|---|---|---|---|
| 0 | 재현성 고정 | — | `env_probe.{json,md}` | — | 본 레포 |
| 1 | 서열 세트 정의 | IPR FASTA | `sequence_set.tsv` | A | 본 레포 |
| 2 | **클러스터 스윕** | IPR FASTA | `sweep_metrics.tsv` 외 | A | 본 레포 |
| 3 | **좌표계 구축** | `chcodh2_alignment.tsv` | `residue_matrix.tsv` 외 | A→B | 본 레포 |
| 4 | 클러스터 특성화 | 1·2·3 산출 | 클러스터 주석표 | A | 미착수 |
| 5 | 근연성 정의 | 2·4 산출 | 티어 T1/T2/T3 | A | 미착수 |
| 6 | 계층형 보존도 | 3·5 산출 | 위치 테이블 생성 | B | 미착수 |
| 7 | 구조 주석 결합 | PDB + CAVER | 위치 테이블에 추가 | B | 미착수 |
| 8 | 설계 판정 | 6·7 | Fix/Bias/Free | B | 미착수 |
| 9 | 규칙 검증 | 논문 기지 위치 | 합격/회귀 | B | 미착수 |

**이번 작업 범위는 0~3** — mmseqs 를 직접 다루는 부분이다. 4 이후는 2·3의 결과를 보고 설계한다.

```
IPR010047 FASTA (1,827)
   │
   ├─[1]→ sequence_set.tsv ──────┐  QC 는 삭제가 아니라 플래그
   │                              │
   └─[2]→ sweep_metrics.tsv       ├─→ master 테이블 (축 A)
          sweep_assignments.tsv ──┘
                                        │
chcodh2_alignment.tsv ─[3]→ equiv_map ──┴─→ residue_matrix.tsv (1827 × 636)
                              │                    │
                              └→ known_position_equivalents.tsv   ← 클러스터 특성화에 바로 쓰임
                                                   │
                                                  [6] 보존도 → 위치 테이블 (축 B)
```

---

## 설계 원칙

1. **삭제하지 않는다.** 1,827행은 끝까지 1,827행이다. QC 는 불리언 플래그 컬럼으로만.
   기준을 바꿀 때 처음부터 다시 돌릴 필요가 없어야 한다.
2. **임계값은 목적별로 따로.** 클러스터 특성 비교와 보존도 계산용 근연 세트는 서로 다른
   해상도를 쓸 수 있다. 하나로 통일하려 하지 않는다. 단계 2가 그 근거를 만든다.
3. **좌표계는 하나.** ChCODH2 서열 번호, Met1 = 1, 636 aa. PDB 번호는 단계 7에서 한 번만 변환한다.
4. **조용히 틀린 숫자가 최악이다.** 가정이 깨지면 경고가 아니라 예외. 행 수·조인율·결측률을
   항상 세어서 stdout 과 `params.json` 양쪽에 남긴다.
5. **명령과 버전을 산출물 옆에.** 모든 산출 디렉터리에 `params.json` (명령·버전·입력 sha256·행 수).
   2026-04 의 `search` 명령이 유실된 것이 이 원칙이 생긴 이유다.
6. **의존성 최소.** pandas + numpy + 표준 라이브러리뿐. scipy / sklearn / Biopython 불필요 —
   ARI·NMI 는 `lib/common.py` 에 numpy 로 직접 구현되어 있다. conda base 에서 설치 없이 돈다.

---

## 번호 체계 검증

`lib/chcodh2.fasta` 의 636 aa 참조서열에 대해 논문 보고 잔기 16개가 전부 일치함을 확인했다
(`common.verify_known_positions()`).

| 위치 | 역할 | 출처 |
|---|---|---|
| F41 | D-cluster 방향족, F41A 에서 활성 소실 | Nat Commun 2024 |
| E43 | viologen 정전기 결합 / 터널 #1 입구 | Nat Commun 2024, ChemBio 2025 |
| R57, N59 | R57G/N59L 에서 EV 친화도 10배 개선 | Nat Commun 2024 |
| A559, V565, I580 | 터널 병목 triad. A559W/H 로 O2 내성 61–148배 | Nat Catal 2022 |
| V610 | 터널 #2 병목. V610H 추가로 O2 내성 개선 | ChemBio 2025 |
| Q206, K450, L583, I586, T593, T597, S599, I603 | 비선택 터널 #1–#3 입구 | ChemBio 2025 |

이 16개는 **단계 9의 정답지**다. 설계 판정 규칙이 이들을 올바로 분류하지 못하면 규칙이 틀린 것이다.
단계 3의 `known_position_equivalents.tsv` 가 각 상동체의 해당 자리 잔기를 표로 만든다.

---

## epel-af2 에서 실행

```bash
cd <이 레포>
source pipeline/config.sh          # 경로를 한 곳에서 고친다

# 단계 0 — 먼저 이것부터. 인계문서 8절 빈칸을 채우고 219 를 재현한다.
python3 pipeline/stages/00_env_probe.py --reproduce --recover-search

#   → env_probe.md 의 mmseqs 절대경로를 config.sh 의 CHCODH2_MMSEQS 에 박고 다시 source
#   → 복원된 search 명령을 selection_params.txt 에 기록

python3 pipeline/stages/01_sequence_set.py
python3 pipeline/stages/02_cluster_sweep.py --markers pipeline/lib/markers.txt
python3 pipeline/stages/03_equiv_map.py
```

각 스크립트는 `--help` 로 옵션을 보여주고, 경로는 전부 `--옵션` 으로 덮어쓸 수 있다.
데이터가 없는 환경에서는 `02` 의 `--dry-run` 과 `--from-existing` 으로 경로를 확인할 수 있다.

---

## 먼저 정해야 할 것

파이프라인은 이 셋을 정하지 않아도 0~3 까지 돈다. 단계 7 이후에 필요하다.

1. **템플릿 구조** — 1SU7 (WT) / 9IYM (A559W/V610H) / 8X9D (R57G/N59L).
   O2 내성을 유지하려면 559·610 변이를 Fix 해야 하는데, 그러면 "야생형 보존도" 와 "설계 목표" 가
   충돌하는 위치가 생긴다. 이 충돌을 테이블에서 어떻게 표기할지 단계 8 규칙에 미리 넣어야 한다.
2. **모노머 설계인가 다이머 설계인가** — D-cluster 가 다이머 계면에 걸쳐 있어 tied positions 처리가 달라진다.
3. **유전자 주변 맥락 주석** (CooC / CooF / CooT / CooJ) — 넣으면 maturase PPI 를 클러스터 단위로
   다룰 수 있지만 1,827개의 지놈 맥락을 받아와야 해서 단계 4에서 가장 비싸다. 2순위로 미뤄도 나머지는 돈다.

## 알려진 미해결 항목

- `chcodh2_alignment.tsv` 를 만든 **mmseqs search 명령이 유실**되었다 (인계문서 2절).
  `fident` 가 `cf` → `cband` → 대표 선정 → 도킹 박스 중심을 전부 결정하므로 논문 Methods 에 필요하다.
  단계 0의 `--recover-search` 가 노트북 셀 이력과 셸 히스토리를 훑는다.
- 같은 파일이 1,827행보다 적을 가능성이 있다. mmseqs search 기본값은 `--max-seqs 300`,
  e-value 0.001 이라 상동체가 잘렸을 수 있다. 단계 3의 감사 단계가 이것을 먼저 센다.
- 인계문서는 `--cov-mode 1` 을 "더 짧은 쪽 길이 기준" 으로 적었다. 실제 mmseqs 정의와 대조가 필요하다
  (219 라는 숫자 자체는 실제 옵션으로 얻은 것이므로 재현에는 영향 없음 — 논문 Methods 서술의 문제다).
