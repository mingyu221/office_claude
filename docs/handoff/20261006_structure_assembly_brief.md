# 구조는 어떻게 만들어졌나 — 예측부터 조립까지

| | |
|---|---|
| 작성일 | 2026-10-06 |
| 대상 | 도킹 세션 |
| 묻는 것 | 「이 321개 `*_dimer.pdb` 는 어디서 온 파일인가」 |
| 코드 | `structure_pipeline.ipynb` (셀 번호를 본문에 적어둠) |

```
BASE  = /mnt/af2results/mingyu/workspace/docking_simulation
INPUT = $BASE/input
SPLIT = $BASE/result/homolog_split
```

> **먼저 결론** — 구조를 **다시 만들 필요도, 만지실 필요도 없습니다.**
> `$INPUT/assembled/{accession}_dimer.pdb` 321개가 완성품입니다.
> 이 문서는 그 파일의 출처와, 공변량 컬럼이 무슨 뜻인지를 설명합니다.

---

## 0. 전체 흐름 한 장

```
IPR010047 서열 1,827
      │
      │  선정 (별도 문서)
      ▼
  대상 321 서열          core 212 + expansion 109
      │
      ├─ ① 단량체 확보 ───────────────────────────────┐
      │     AFDB 에 있으면 내려받기          298개    │
      │     없으면 ColabFold 로 예측           23개    │  셀 6~8
      │                                               │
      ▼                                               │
  단량체 321개  (사슬 하나뿐. 이대로는 도킹 못 함)  ◀─┘
      │
      ├─ ② 이량체 조립 ───────────────────────────────┐
      │     결정 이량체(1SU6 / 6OND) 위에 중첩         │  셀 9
      │     아키텍처에 따라 참조를 고름                │
      │     전장·국소 두 방식 중 하나 채택             │
      │                                               │
      ▼                                               │
  이량체 321개  $INPUT/assembled/*_dimer.pdb        ◀─┘
      │
      └─ ③ 품질 검증 → dock_targets.tsv            셀 9~10
```

---

## 1. 왜 구조를 만들어야 했나

상동체 321개 중 **AFDB 에 모델이 없는 것이 23개**였다.

| 출처 | 개수 | 비율 |
|---|---|---|
| AlphaFold DB | 298 | 93% |
| **ColabFold (직접 예측)** | **23** | **7%** |

### 1.1 왜 23개가 AFDB 에 없었나

**전부 `A0AB*` 또는 `A0AC*` accession 이다.** 이 대역은 UniProt 에 최근
등록된 서열이고, AFDB v4 (2022) 이후에 들어왔다. 데이터베이스가 아직
따라오지 못한 것뿐이다.

### 1.2 그럼 23개를 버리면 안 됐나 — 안 된다

버리면 **선정이 "AFDB 가 무엇을 담았나" 에 좌우된다.**
AFDB 등재는 단백질의 성질이 아니라 데이터베이스의 수록 상태다.

그리고 실제로 편향이 확인됐다 — 누락분이 **전부 최근 등록 대역**이므로,
그것으로 걸렀다면 **최근 발견된 계통(대개 환경·메타게놈 유래)을 통째로
잃었을 것**이다. 이 과의 다양성 중 가장 새로운 부분이다.

그래서 23개를 직접 예측해 넣었다.

---

## 2. ColabFold 예측 조건 (셀 7)

```
colabfold_batch  todo.fasta  $INPUT/colabfold
    --num-models 1  --num-recycle 3
    --model-type alphafold2_ptm  --rank plddt
```

| | |
|---|---|
| 실행 환경 | `/home/unsit/localcolabfold/.pixi/envs/default` |
| GPU | **0번 고정** (1번은 RoseTTAFold2-PPI 작업 중이었음) |
| 산출 | `$INPUT/colabfold/{accession}.pdb` (rank_001 복사본) |
| 요약 | `$INPUT/colabfold/prediction_summary.tsv` |
| pLDDT | 평균 **94.3** (AFDB 298개는 95.6) |

### 2.1 디렉터리에 무엇이 들어 있나

ColabFold 원본 출력은 파일명이 길고 `{accession}.pdb` 가 아니다.
**조립 단계에서 쓰기 쉽게 rank_001 을 `{accession}.pdb` 로 복사해 둔다** (셀 8).

```
$INPUT/colabfold/
    A0ABU0B233.pdb                                  ← 쓰는 것. rank_001 복사본
    A0ABU0B233_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
    A0ABU0B233_scores_rank_001_....json             ← pLDDT/PAE 원본
    A0ABU0B233.a3m                                  ← MMseqs2 MSA
    A0ABU0B233_coverage.png  A0ABU0B233_plddt.png
    ...  (23개 × 위 세트)
    prediction_summary.tsv       accession · n_res · plddt · tier · cband
    todo_HHMMSS.fasta            그 실행에 넣은 입력
    colabfold.log                실행 로그. 'no GPU detected' 확인용
    colabfold.pid                실행 중 PID (끝나면 삭제됨)
```

**AFDB 쪽은 이름 규칙이 다르다.** 내려받은 그대로 두었다.

```
$INPUT/afdb/
    AF-A0A2N7PLW0-F1-model_v4.pdb      ← accession 이 중간에 박혀 있다
    ...  (298개 + 과거에 받아둔 것들)
```

그래서 단량체를 찾는 코드는 두 규칙을 모두 본다.

```python
def mono_path(acc):
    for root, pat in ((AFDIR, f"*{acc}*.pdb"),        # AF-{acc}-F1-model_v4.pdb
                      (PRED,  f"{acc}.pdb"),          # 정규화한 복사본
                      (PRED,  f"{acc}*rank_001*.pdb")):  # 복사 전이면 원본
        hits = sorted(root.glob(pat))
        if hits:
            return hits[0]
    return None
```

**조립본은 규칙이 하나다** — `$INPUT/assembled/{accession}_dimer.pdb`.
도킹에서 실제로 여는 파일은 이것뿐이므로 위 혼선은 상관없다.

**AFDB 와 같은 AlphaFold2 지만 MSA 파이프라인이 다르다.**
AFDB 는 JackHMMER/HHblits, ColabFold 는 MMseqs2 를 쓴다. 같은 서열을 줘도
구조가 미세하게 다를 수 있으므로 **`struct_source` 컬럼으로 기록해 공변량으로
둔다.** 조립 지표상 차이는 아래와 같다.

| `struct_source` | RMSD | 충돌 | 임계충돌 | `center_dev` | pLDDT |
|---|---|---|---|---|---|
| AFDB (298) | 1.38 | 13.45 | 1.17 | 0.54 | 95.59 |
| ColabFold (23) | 1.46 | **33.05** | **1.18** | 0.93 | 94.31 |

`충돌` 만 2.5배이고 **`임계충돌`(1.18 vs 1.17)은 사실상 같다.**
겹침이 결합 부위가 아니라 **주변부**에 있다는 뜻이다. 그래서 ColabFold
구조를 그대로 쓴다.

### 2.2 겪은 함정 두 개 (참고)

주피터에서 돌리느라 환경변수 문제를 두 번 맞았다. 재예측할 일이 생기면
`build_env()` 를 그대로 쓸 것.

| 증상 | 원인 | 조치 |
|---|---|---|
| 임포트 단계에서 즉사 | 주피터의 `MPLBACKEND=module://matplotlib_inline...` 가 자식 프로세스로 새어 들어감 | 주피터·conda 변수 제거 후 `MPLBACKEND=Agg` |
| `no GPU detected` → CPU 폴백 | pixi 환경에 `libcusparse.so.12` 가 있는데 로더 경로에 없었음 | `LD_LIBRARY_PATH` 에 `nvidia/*/lib` 전부 추가 |

> 두 번째는 처음에 "환경변수 누수" 로 오진했다. 터미널에서 `jax.devices()` 를
> 직접 불러도 같은 증상이 나와서 알았다. **CPU 폴백은 조용히 일어나고
> 결과물은 정상으로 보인다** — 로그에서 `no GPU detected` 를 확인해야 한다.

---

## 3. 왜 단량체로는 안 되나

AlphaFold 가 주는 것은 **사슬 하나**다. 그런데 EV 결합 부위는 **이량체 계면에
있다.**

1SU6 (ChCODH2 결정구조) 의 F41 주변 12 Å 을 실측하면 이렇다 (셀 2).

```
잔기 36개 중  방향족 1개 (PHE61, 10.77 Å)  ·  CYS 7개

CYS   A/47(5.56)   A/39(5.58)   A/48(9.34)   A/56(9.55)
      A-2/39(10.26)   A-2/56(11.62)   A-2/47(11.75)
        ↑ 반대편 protomer
```

**반대편 사슬의 Cys 세 개가 10~12 Å 안에 들어온다.**
[4Fe-4S] D클러스터가 두 subunit 에 걸쳐 배위되므로 **단량체로는 결합 자리가
완성되지 않는다.** 단량체에 도킹하면 한쪽 벽이 없는 포켓에 리간드를 넣는 셈이다.

### 3.1 덤으로 확인된 것 — 서열 분류가 실제 배위 자리를 짚는다

같은 출력에서 Cys39·Cys47 이 F41 을 5.5 Å 거리에서 양쪽으로 끼고 있다.
잔기 번호를 세면

```
C39 - x40 - F41 - x42 x43 x44 x45 x46 - C47
 C  .  [FY]  .  .  .  .  .   C              ←  C.[FY].{5}C
```

아키텍처 분류(SF4/SF2)에 쓴 **정규식 그대로**다. 1,827개를 서열로 긁은 분류가
휴리스틱이 아니라 실제 금속 배위 기하를 짚고 있었다는 구조적 증거다.

---

## 4. 왜 결정학적 대칭을 못 쓰나

결정구조라면 PyMOL `symexp` 로 대칭 짝을 꺼내면 된다.
**AFDB·ColabFold 모델에는 그게 없다** — 공간군·단위격자 정보가 애초에 없다.
예측 모델은 결정이 아니기 때문이다.

그래서 **결정 이량체의 기하를 예측 단량체에 옮겨 심는 방식**을 썼다.

---

## 5. 조립 방법 (셀 3 · 셀 9)

### 5.1 원리

```
참조 이량체        사슬 A ─── 사슬 B        (1SU6 또는 6OND, 실제 결정구조)
                     │         │
                     │ 중첩    │ 중첩        ← 등가 잔기 CA 끼리 맞춤
                     ▼         ▼
상동체 단량체      복사본1   복사본2
                     └─ 결과 ─┘
                   사슬 A      사슬 B        → {accession}_dimer.pdb
```

같은 단량체를 **두 번 복사**해서, 하나는 참조의 A 사슬 위에, 다른 하나는
B 사슬 위에 각각 중첩한다. 그러면 참조의 이량체 기하가 상동체에 전사된다.

### 5.2 어떤 잔기끼리 맞추나 — `equiv_map`

mmseqs 정렬 (`chcodh2_alignment.tsv`) 에서 **ChCODH2 잔기번호 → 상동체
잔기번호** 대응표를 만든다.

```python
def equiv_map(r):
    qi, ti, out = int(r.qstart) - 1, int(r.tstart) - 1, {}
    for qc, tc in zip(r.qaln, r.taln):   # q = ChCODH2(query), t = 상동체(target)
        if qc != "-": qi += 1
        if tc != "-": ti += 1
        if qc != "-" and tc != "-": out[qi] = ti
    return out
```

상동체마다 **잔기번호가 전부 다르다.** `resSeq == 41` 로 찾으면 엉뚱한 잔기가
잡힌다. 이 대응표가 번호 차이를 흡수한다.

### 5.3 아키텍처마다 참조가 다르다 — 이게 중요하다

이 과에는 Fe-S 클러스터 배위가 다른 두 계통이 있고, **이량체 기하도 다르다.**

| 아키텍처 | Cys 간격 모티프 | 참조 이량체 |
|---|---|---|
| SF4 · SF4v · SFx | `C.[FY].{5}C` | **1SU6** |
| SF2 | `C.[FY]C` | **6OND** |

**2×2 로 교차해 확인했다** (별도 노트북 `veri_dimer_generation`).

| | 1SU6 에 붙임 | 6OND 에 붙임 |
|---|---|---|
| SF4 구조 | **임계충돌 0.50** | 3.00 |
| SF2 구조 | 4.57 | **임계충돌 1.14** |

자기 참조를 쓰면 임계충돌이 3~4배 낮아진다. **참조를 잘못 고르면 반대편
사슬이 결합 부위로 밀려 들어온다.** 예외 없이 아키텍처대로 고른다.

321개 분포 — SF4 254 · SF2 55 · SF4v 7 · SFx 2 · both 2 · none 1.
`both`·`none`(모티프로 못 가린 3개)은 기본값 1SU6 로 간다.

### 5.4 6OND 를 먼저 1SU6 좌표계로 옮긴다 (셀 2)

두 참조가 각자의 좌표계에 있으면 조립 결과도 두 좌표계로 갈린다.
그러면 박스 중심을 비교할 수 없다. 그래서 **6OND 를 1SU6 위에 먼저 중첩**해서
좌표계를 통일했다.

```
6OND → 1SU6 중첩    RMSD 1.66 Å   등가쌍 619
```

**결과적으로 321개 전부가 1SU6 좌표계에 놓여 있다.** 그래서 `center_dev`
(각 구조 41등가 CA 와 1SU6 PHE41 CA 의 거리) 라는 값이 의미를 갖는다.

### 5.5 전장이냐 국소냐 — 두 방식을 다 시도한다

| 방식 | 맞추는 범위 | 최소 쌍 수 |
|---|---|---|
| `global` | 사슬 전체 | 250 |
| `local` | 41등가 ± 60 잔기 창 | 40 |

전장 중첩은 전체 접힘을 잘 맞추지만 결합 부위가 조금 틀어질 수 있다.
국소 중첩은 반대다. **둘 다 만들어 보고 하나를 고른다.**

채택 순서는 이렇다.

```
1순위   center_dev 2.0 Å 통과 여부
2순위   임계충돌이 적은 쪽
3순위   RMSD 가 작은 쪽
```

`RMSD > 5.0 Å` 후보는 아예 버린다 (중첩 자체가 실패한 것).

구간별 채택 비율 — 유사도가 낮을수록 `local` 이 많이 뽑힌다.

| `cband` | global | local | local % |
|---|---|---|---|
| ≥60 | 7 | 3 | 30.0 |
| 50-60 | 9 | 14 | 60.9 |
| 40-50 | 64 | 97 | 60.2 |
| 30-40 | 42 | 71 | 62.8 |
| <30 | 8 | 6 | 42.9 |

**`mode` 컬럼이 어느 쪽을 썼는지 기록한다.** 공변량으로 쓸 것.

---

## 6. 품질 지표 — 공변량 컬럼의 뜻

`$SPLIT/assembly_manifest_all.tsv` 에 321행이 있다.

| 컬럼 | 무엇을 재나 | 어떻게 쓰나 |
|---|---|---|
| **`RMSD`** | 중첩이 제대로 됐나 (Å) | **배제 기준.** > 5.0 이면 실패 |
| `쌍` | 중첩에 쓴 등가 CA 쌍 수 | 참고 |
| `충돌` | A사슬 중원자 중 B사슬 중원자가 2.0 Å 안에 있는 것 | 참고 (주변부 포함) |
| **`임계충돌`** | 그중 **박스 중심 12 Å 안**에 있는 것 | 공변량. ≥10 이면 태그 |
| **`center_dev`** | 41등가 CA 와 1SU6 원점의 거리 | 공변량. **배제 기준 아님** |
| `d_arom` | 가장 가까운 방향족 CA 까지 거리 | 기록만 |
| `n_B_near` | 중심 최근접 6 CA 중 B사슬 개수 | 321개 **전부 0** — 정보 없음 |
| `mode` | `global` / `local` | 공변량 |
| `struct_source` | `AFDB` / `ColabFold` | 공변량 |
| `pLDDT` | 단량체 예측 품질 (평균 B-factor) | 공변량 |

### 6.1 `충돌` 과 `임계충돌` 은 다른 것을 잰다

```
충돌      사슬 전체에서 센다.          원자 수. A→B 단방향
임계충돌  박스 중심 12 Å 안에서만 센다
```

ColabFold 구조가 `충돌` 33 vs AFDB 13 으로 2.5배인데 `임계충돌` 은 같다 —
**겹침이 결합 부위가 아닌 데 있다**는 뜻이고, 그래서 문제가 아니다.
**점수에 영향을 줄 수 있는 쪽은 `임계충돌` 이다.**

짝수/홀수가 번갈아 나오는 것은 C2 대칭 때문이며 정상이다.

### 6.2 `center_dev` 는 오차가 아니다

```
박스 중심 = 각 구조 자신의 41등가 잔기 CA      ← 1SU6 고정 좌표가 아님
```

규칙이 **단백질 자기 기준**이므로, 그 상동체의 축이 5~10 Å 이동해 있으면
박스도 같이 옮겨간다. **박스는 언제나 그 단백질 자신의 부위를 반변 15 Å 로
온전히 덮는다.**

따라서 `center_dev` 는 **「이 상동체의 축이 ChCODH2 로부터 얼마나 이동했나」**
를 재는 기술자(descriptor)다. 큰 값은 이상이 아니다.

실측 분포 — 295/321 이 1.0 Å 이내, **최대 9.69 Å (4건)**.
그 4건 목록과 "문제로 보고하지 않아도 된다" 는 설명은
`docking_brief.md` 4절에 있다.

참고로 **한 잔기 안에서도 CA 와 실제 스태킹 고리가 3.77 Å 떨어져 있다**
(1SU6 PHE41 실측). `center_dev 5.25` 는 그 1.4배 수준이다.

---

## 7. 한 번 크게 틀렸던 것 (왜 믿을 만한가)

**ColabFold 23개가 전부 SF4 로 기록되어 1SU6 에 붙었던 적이 있다.**

선정 단계에서 AFDB 게이트를 뺐을 때(1,480 → 1,722) 새로 들어온 242개는
주석 테이블에 없던 서열이라 `f41_pos` 와 `arch` 가 **전부 NaN** 이 되었다.
그중 선정된 23개가 **ColabFold 23개와 정확히 일치**했다.

| 결측 | 결과 |
|---|---|
| `f41_pos` | 박스 중심을 계산할 수 없음 |
| `arch` | **전부 SF4 기본값 → 전부 1SU6 참조** |

두 번째가 컸다. IPR010047 의 SF2 비율이 23% 인데 ColabFold 23개가 전부 SF4 라는
것 자체가 이상 신호였다. 서열 모티프로 다시 스캔하니 **5개가 실제로 SF2**였다.

```
A0ABU0B233 · A0ABM8AV39 · A0ABT2EXU1 · A0ABY1NSA6 · A0ABR5SJT4
   →  6OND 로 다시 조립
```

**지금 쓰는 321개는 이 복구 이후의 것이다** (셀 4 가 `f41_pos` 를 정렬에서,
`arch` 를 서열 모티프에서 다시 계산한다). 현재 `arch` 결측 0, 참조 배정은
전부 모티프 근거다.

> 이 이력을 적어두는 이유 — 도킹 세션이 예전에 `center_dev` 큰 구조를
> 문제로 보고해 주셨고, 추적하니 원인이 이 하나였다. **그 보고가 버그를
> 잡았다.** 지금은 고쳐진 상태다.

---

## 8. 조립이 점수를 왜곡하지 않는다는 검증

조립이라는 인공 절차를 거쳤으니, 그것 자체가 점수를 바꾸는지 먼저 확인했다.

| | med | 차이 |
|---|---|---|
| 288 run Wildtype 결정구조 재도킹 | **−4.380** | 정본 −4.38 과 **0.000** |
| 조립한 ChCODH2 이량체 | **−4.30** | **+0.080** (0.9 × IQR) |

첫 줄은 **파이프라인이 288 run 과 일치**함을 보인다.
둘째 줄은 **조립 효과가 0.080** 이라는 것인데, 이는 결정구조 자신의 50런
범위 **0.11 보다 작다.**

**조립 절차는 도킹 점수를 왜곡하지 않는다.**

| | |
|---|---|
| 기록 | `$BASE/result/validation/validation_stage1.tsv` |
| 자기대조 구조 | `$INPUT/assembled/Q9F8A8_ChCODH2_assembled.pdb` |
| 정본 | `$BASE/result/Veri_Ca_diameter/engines_288_CA_r15.csv` |

ChCODH2 자기대조 구조는 **잔기번호가 실제로 41** 이므로, 중심 계산 코드를
검산할 때 쓰면 편하다.

---

## 9. 최종 상태

```
도킹 대상 321        core 212 + expansion 109      클러스터 212
배제         1        A0AC61SBG3 (RMSD 14.91, 중첩 실패)
대체         1        A0AC61SBG3 → A0A7J2S3E9 (같은 클러스터 다음 medoid)
태그         5        유지. 공변량으로 기록
실패         0
```

`A0AC61SBG3` 은 두 방식 모두 RMSD 5.0 을 넘겨 좌표계가 달라졌다.
다른 구조들과 비교가 불가능하므로 제외하고, 같은 클러스터의 다음 후보로
바꿨다. **클러스터 212개가 전부 대표를 갖는다.**

제외된 파일은 `$INPUT/assembled/_excluded/` 로 옮겨 두었다 —
`$INPUT/assembled/*.pdb` 글롭에 섞이지 않는다.

### 9.1 결합부위 실측 분포 (321개)

| 항목 | 값 |
|---|---|
| 41등가 잔기 | PHE 256 · TYR 38 · ILE 18 · HIS 3 · MET 3 · ASP 2 · LEU 1 |
| 방향족(PHE+TYR) | **294/321 = 91.6%** |
| Gly-x-Gly 보존 | 134/321 = 41.7% |
| `n_B_near` | 전부 0 |
| `center_dev` | 295개가 ≤1.0, 최대 9.69 |

**비방향족 27건은 정상이다.** 41 자리가 방향족이냐 아니냐가 이 연구의
종속변수이므로, 비방향족이라고 해서 잘못 잡은 것이 아니다.

---

## 10. 파일

| 무엇 | 경로 |
|---|---|
| **조립 구조 321개** | `$INPUT/assembled/{accession}_dimer.pdb` |
| 제외분 | `$INPUT/assembled/_excluded/` |
| 자기대조 | `$INPUT/assembled/Q9F8A8_ChCODH2_assembled.pdb` |
| 단량체 (AFDB) | `$INPUT/afdb/*.pdb` |
| 단량체 (ColabFold) | `$INPUT/colabfold/{accession}.pdb` |
| 예측 요약 | `$INPUT/colabfold/prediction_summary.tsv` |
| 참조 이량체 | `$INPUT/1su6_dimer.cif` · `$INPUT/6ond_dimer.cif` |
| **도킹 목록** | `$SPLIT/dock_targets.tsv` |
| **공변량** | `$SPLIT/assembly_manifest_all.tsv` |
| 결합부위 원자료 | `$SPLIT/site_check.tsv` |
| 태그 5건 | `$SPLIT/tagged_highclash.tsv` |
| 제외 기록 | `$SPLIT/excluded_assembly.tsv` |
| 대체 기록 | `$SPLIT/replacements.tsv` |
| 정렬 | `$BASE/result/homolog_selection/chcodh2_alignment.tsv` |
| 전체 기준표 | `$SPLIT/master_1827.tsv` |

### 10.1 노트북 셀 대응

```
셀 2    참조 이량체 로드 · 6OND → 1SU6 좌표계 통일 · 1SU6 F41 주변 실측
셀 3    equiv_map · assemble_on · clash_stats · assemble_best
셀 4    master 정합성 복구 (7절의 버그를 고치는 곳)
셀 6    AFDB 다운로드
셀 7    ColabFold 실행 + 완료 대기 (build_env 가 환경 함정 처리)
셀 8    예측 정리 + pLDDT
셀 9    조립 고리 — 대상 산출 · 조립 · 검증 · 대체를 반복
셀 10   dock_targets.tsv 확정
```

---

## 11. 도킹 세션이 알아야 할 것 다섯 줄

1. `$INPUT/assembled/*_dimer.pdb` **321개가 완성품이다.** 다시 만들지 않는다.
2. **박스 중심은 `dock_targets.tsv` 의 `f41_pos`** 다. `resSeq == 41` 로 찾으면 안 된다.
3. `center_dev` 가 커도 정상이다. 그 상동체의 축이 이동한 양이다.
4. 공변량 여덟 개를 점수와 함께 돌려주면 이쪽에서 회귀로 가른다.
5. 조립이 점수를 왜곡하지 않는다는 것은 이미 검증됐다 (8절). 다시 할 필요 없다.
