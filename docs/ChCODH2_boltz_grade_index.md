# Boltz Ni 배위 등급 — 자료 위치와 내용

`BASE = /mnt/af2results/mingyu/workspace/ppi_discovery`

---

## 1. 파일 위치

### 예측 구조 (원본)

| 경로 | 무엇 | 구조 수 |
|---|---|---|
| `result/boltz/out/` | Track B 본 예측 (326 입력) | 321 |
| `result/boltz/out_focus/` | folddisco 금속 히트 중 pLDDT ≥ 70 | 27 |
| `result/boltz/out_cmp/` | 타 균주 직교체 비교 | 50 |
| `result/boltz/out_rep/` | 반복 + 대조군 (13건 × 5샘플) | 65 |
| | **합계** | **411** |

파일명 규칙: `{입력이름}_model_0.cif`

### 표 (가공)

| 파일 | 무엇이 들어 있나 | 만든 셀 |
|---|---|---|
| **`result/table/ni_site_grade.csv`** | **★ 등급표.** 단백질별 grade / Cys / His / Asp·Glu / n_donor / min_dist / donors / why | `28a-8b` |
| `result/table/boltz_ni_placement.csv` | 411 구조의 Ni 좌표와 배위 거리 (`CUT_NI 3.5` 기준) | `28a-8b` |
| `result/table/ortholog_reproducibility.csv` | 같은 단백질을 다른 배치·다른 균주에서 잰 값 — **재현성 근거** | `28a` |
| `result/table/boltz_replicate_summary.csv` | 반복 13건 × 5샘플 요약 + 대조군 ipTM | `part7_boltz_rep` |
| `result/table/metal_motif_3strain_boltz.csv` | 타 균주 직교체 50건 | `out_cmp` 집계 |
| `result/table/trackB_boltz2_ranked.csv` | 321건 전체 (ipTM 포함 — **순위에 쓰지 않음**) | Track B |
| `result/table/trackB_focus_metal_motif.csv` | focus 27건 | |
| `result/table/trackB_interface_match.csv` | Ni 이 계면에 놓였는지 | |
| `result/table/trackB_dimer_reconsider.csv` | **이량체 재검토 대상** (Q7 이 만든다 — 아직 없음) | `Q7` |

---

## 2. 어떤 기준으로 매겼나

```
1. 예측 구조에서 NI 좌표를 찾는다
2. 2.6 Å 안의 원자를 모은다        ← 후보 단백질 사슬만. ChCODH2 는 제외
3. 잔기 단위로 중복 제거
4. 도너 종류로 등급
```

| 상수 | 값 | 무엇 |
|---|---|---|
| `CUT_NI` | **3.5 Å** | Ni 이 어디엔가 놓였나 → **411 / 411 전부 통과.** 변별력 0 |
| `TIGHT` | **2.6 Å** | 배위로 세는 선 → **여기서 등급이 갈린다** |

| 등급 | 조건 | 뜻 |
|---|---|---|
| **A** | Cys ≥ 4 | 티올레이트 4개 — 금속을 잡았다 놓는 운반형 자리 |
| **B** | Cys 3, 또는 Cys 2 + His | 구조적 금속 자리 |
| **C** | Cys 1–2, 또는 His ≥ 3 | 촉매 금속 자리 가능성 |
| **D** | 그 외 (산소 위주) | 표면 부착 의심 |

---

## 3. 결과

**411 구조 중 Ni 을 후보 쪽 사슬로 가져간 것 = 22개.** 그 22개를 채점했다.

| 등급 | 수 | 단백질 | 최단 배위 |
|---|---|---|---|
| **A** | **4** | `QJZ13803.1` **HslO** | 2.18 Å |
| | | `Y19-AKE60400.1` HslO 직교체 | 2.17 Å |
| | | **`QJZ12568.1`** hypothetical HO396_09830 | **2.23 Å** |
| | | `QJZ13396.1` **GspE** | 2.29 Å |
| **B** | 9 | HisB · GmhB · Tgt (세 균주 직교체 포함) | 1.86–2.14 Å |
| **C** | 1 | HpaA | — |
| **D** | 8 | **YeiR** · PhnP · NagA · Ser/Thr 포스파타제 · phage tail | — |

> grade A 4개 중 **둘은 같은 단백질의 직교체**(BL21 HslO · Y19 HslO)다.
> 서로 다른 단백질로 세면 **3종**이다.

### grade A 의 배위 잔기

```
QJZ12568.1      CYS99  CYS102 CYS112 CYS115
QJZ13396.1      CYS391 CYS394 CYS424 CYS427     GspE
Y19-AKE60400.1  CYS230 CYS232 CYS263 CYS266     HslO 직교체
QJZ13803.1      CYS230 CYS232 CYS263 CYS266     HslO (BL21)
```

### YeiR 가 D 로 내려간 이유

```
n_donor = 1        2.6 Å 안에 든 것이 HIS209 하나뿐
CYS65  = 2.69 Å    잘림 (Ni–S 는 2.2 Å 근처여야 한다)
```

앞선 후보 목록의 1번이었다. **계통으로 고른 후보를 기하가 내렸다.**
단, `TIGHT` 를 3.0 으로 올리면 D 가 C 로 올라오기 시작한다 — 경계는 우리가 정한 것이다.
그리고 YeiR 는 이량체라 **단량체로 접은 것이 불리했을 수 있다**(Q7 에서 재검토 중).

---

## 4. 재현성 — 이 축을 믿는 근거

`ortholog_reproducibility.csv`

| 단백질 | BL21 | MG1655 | Y19 | 배위 잔기 |
|---|---|---|---|---|
| HslO | 2.18 | — | 2.17 | CYS232 CYS230 CYS263 CYS266 (동일) |
| HisB | 1.91 | 1.92 | 1.86 | CYS103 HIS95 CYS101 CYS93 (동일) |
| Tgt | 2.14 | 2.12 | 2.14 | HIS333 CYS302 CYS304 CYS307 (동일) |
| PhnP | 2.05 | 2.08 | 2.00 | ASP164 HIS78 … |

**Δ 0.01–0.08 Å** — 다른 배치, 다른 GPU, 다른 균주의 직교체인데 같은 잔기를 같은 거리로 잡는다.
ipTM 이 같은 입력끼리 0.28 씩 벌어지는 것과 정반대다.

**동시에 이 축은 균주를 못 가른다** — HslO·HisB·Tgt·PhnP 가 세 균주 모두 똑같이 잡는다.

---

## 5. 버린 지표 두 개

| 지표 | 대조군 / 진단 | 판정 |
|---|---|---|
| `ipTM` | CooC1 = **0.300** 인데 후보 8개가 전부 그 위 | **폐기** |
| `ligand_ipTM` | 계면 배치 **0**, 단백질 길이와 **ρ = −0.807** | **폐기** |

**Ni 이 놓인 위치** (`trackB_interface_match.csv`, 348 구조 기준)

| 위치 | 수 |
|---|---|
| 계면 (양쪽 사슬 배위) | **0** |
| ChCODH2 쪽만 | 335 |
| 후보 쪽만 | 13 |

전달 복합체는 한 번도 모델링되지 않았다.

> **숫자 두 개를 구분할 것** — 이 표는 **348 구조** 기준이고(후보 쪽 13),
> 3장의 "후보 쪽으로 가져간 22개"는 **411 구조** 기준이다. 서로 다른 패스에서
> 나온 값이라 분모가 다르다. 발표에 쓸 때 분모를 같이 적을 것.

**단, 배위 거리 자체는 타당하다** — 평균 **2.24 Å** (1.97–3.27). Ni–S 2.2 / Ni–N 2.1 근처다.

---

## 6. 열 이름

`ni_site_grade.csv`

| 열 | 뜻 |
|---|---|
| `grade` | A / B / C / D |
| `prey` | 후보 단백질 (`Y19-` `MG1655-` 접두사가 붙기도 한다) |
| `set` | 어느 폴더에서 나온 구조인가 (`out` / `out_focus` / `out_cmp` / `out_rep`) |
| `n_donor` | 2.6 Å 안의 배위 **잔기** 수 (원자 아님) |
| `Cys` / `His` / `Asp/Glu` | 종류별 개수 |
| `min_dist` | 가장 가까운 배위 거리 |
| `donors` | 배위 잔기 목록 (`CYS112 CYS114` 형태) |
| `strain` | BL21 / MG1655 / Y19 |
| `folddisco_metal` | folddisco 금속 모티프 목록에도 있나 |
| `why` | 등급 사유 |
| `desc` | 어노테이션 |

---

## 7. 이 자료로 만든 발표 구성

`docs/ChCODH2_slides_boltz_ni.md` — S1~S9 + 부록.
슬라이드마다 **올릴 것 / 숫자 / 예상 질문**을 적어 뒀다.
