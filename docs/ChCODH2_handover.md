# 인계 자료 — 2026-09-19 04:00 기준

## 0. 지금 돌고 있는 것

| GPU | 작업 | 시작 | 예상 | 결과 파일 |
|---|---|---|---|---|
| 1 | **CELL 54** Track A 이량체 스크린 (`in_both1_b00` ~) | 09-19 03:40 | **63시간** | `result/rf2ppi/in_{both,fall}{rep}_b{nn}.log` |
| 0 | **Q6b** Boltz `inputs_seqmotif` 400건 | 09-19 03:56 | 미상 (MSA 대기 중) | `result/boltz/out_seqmotif/` |

**확인 방법**

```bash
# RF2-PPI: 끝난 버킷 수 / 54
ls /mnt/af2results/mingyu/workspace/ppi_discovery/result/rf2ppi/in_{both,fall}*_b*.log 2>/dev/null | wc -l

# Boltz: MSA 진행 / 예측 완료
grep -c "Generating MSA" .../result/boltz/out_seqmotif.log
ls .../result/boltz/out_seqmotif/**/*_model_0.cif 2>/dev/null | wc -l
```

또는 노트북에서 **Q6b 를 다시 실행**하면 상태만 갱신된다 (중복 실행 안 함).

> `RATELIMIT` 로그는 정상이다. ColabFold MSA 서버가 무인증 요청에 속도를 걸고
> Boltz 가 5–9초씩 쉬었다 재시도한다. MSA 를 다 받은 뒤에야 예측이 시작되므로
> 초반에 `out_seqmotif` 가 비어 있는 것도 정상이다.

---

## 1. 오늘 무엇이 바뀌었나 (한 문단)

**Track C(Folddisco)를 쓸 수 없다는 것이 확정됐고, 그 자리를 두 축이 대신한다.**
질의 잔기 Cys112/114 의 pLDDT 가 **41 / 45** 다(단백질 평균 79.8). AlphaFold 가
모른다고 표시해 둔 좌표로 기하 검색을 하고 있었다. 대신 **서열 모티프**(좌표를
안 씀)와 **Foldseek**(대조군 통과한 유일한 축)로 간다. 그리고 파이프라인의 입구를
"MG1655 에 없는 것"에서 **"Ni 을 잡을 수 있는 자리를 가진 것"** 으로 바꿨다 —
기존 입구가 3,746개를 평가 없이 남겨두고 있었기 때문이다.

---

## 2. 노트북 지도

| 파일 | 셀 | 무엇 |
|---|---|---|
| `ChCODH2_screen_antigravity.ipynb` | 1–56 | 본체. **CELL 54** 가 지금 돌고 있는 이량체 스크린 |
| `ChCODH2_seq_motif_pipeline.ipynb` | Q1–Q8 | **새 주 파이프라인.** 서열 모티프 → 균주 → 폴드/pLDDT → Boltz |
| `ChCODH2_folddisco_afframe.ipynb` | F1–F11 | Folddisco 실패 규명 + **Foldseek 축** + 발표 표 |
| `ChCODH2_folddisco_sensitivity.ipynb` | S0–S4 | 3kji 조립체 검증, `-d`/`-a` 스윕 |
| `ChCODH2_master_list.ipynb` | M1 | **전체 후보 한 장** (`MASTER_candidate_list.csv/xlsx`) |
| `ChCODH2_candidate_review.ipynb` | R1–R4 | 4축 행렬, 검토용 엑셀 |
| `ChCODH2_motif_tables.ipynb` | M1–M4 | 모티프 계수 표 |

### 실행 순서 (서열 파이프라인)

```
Q1 설정 → Q2 모티프 스캔 → Q3 균주 판정 → Q3b 프로파지 구분
   → Q4 폴드·pLDDT → Q5 점수+입력생성 → Q5b CooC형 이량체 입력
   → Q6b 백그라운드 실행 → (끝나면) Q6 채점 → Q7 이량체 선별 → Q8 이량체 채점
```

---

## 3. 내일 아침에 할 것

1. **Q6b 재실행** — 진행률 확인. 400/400 이면 **Q6 채점**
   - ★ 눈금용 6개(기존 grade A/B: `QJZ11021.1` `QJZ11161.1` `QJZ12568.1`
     `QJZ12582.1` `QJZ13396.1` `QJZ13803.1`)가 **같은 등급으로 재현되는지 먼저** 볼 것.
     재현 안 되면 ChCODH2 를 뺀 것이 결과를 바꾼 것이므로 옛 결과와 나란히 못 놓는다
2. **Q7 → Q5b → Q8** (이량체). GPU 가 비면 `IN_DIR="inputs_cooclike"` 로 Q6b
   - ★ 양성대조군 `CTRL-CooC1_dimer` 가 **`A-계면`** 으로 나와야 한다.
     안 나오면 이량체 예측 자체를 못 믿으므로 아래 결과를 쓰면 안 된다
   - ★ 음성대조군 `NEG-*` 가 A 로 나오면 Boltz 가 Ni 을 아무 데나 넣는 것이다
3. **F11** — HypA·G3E 스윕을 `--exact-tmscore 1` 로 재실행.
   발표에 `tm 0.680` 을 쓰려면 필요하다 (현재 값은 1.0 을 넘는 근사값)
4. **M1** 마스터 리스트 — 위가 쌓인 뒤에 돌리면 전체가 한 장으로 합쳐진다
5. CELL 54 는 63시간짜리다. 끝나면 `mode` 열로 갈라 **위상별로** 집계

---

## 4. 밟으면 안 되는 지뢰 (전부 오늘 실제로 밟았다)

| | 내용 |
|---|---|
| **캐시 키** | 파일명에 **입력을 전부** 넣는다. 질의 잔기·구조 DB 이름. 안 넣으면 코드를 고쳐도 옛 결과를 재사용한다. 오늘 두 번 당했다 |
| **입력 폴더** | Boltz 는 폴더의 **모든 파일**을 입력으로 읽는다. `.json` 하나가 있으면 죽는다. 매니페스트는 폴더 **밖**에 |
| **폴더 정리** | 입력 폴더를 안 비우면 지난 실행의 YAML 이 남아 같이 돌아간다. Q5/Q5b 가 이제 `*.old_<시각>` 으로 옮긴다 |
| **GPU 선택** | 여유 메모리로만 고르면 안 된다. 상대가 아직 메모리를 안 잡았으면 구분이 안 된다. `FORCE_GPU` 로 직접 정한다 |
| **L 상한 프로브** | **일부러 OOM 을 내는** 코드다. 같은 카드에 다른 작업이 있으면 상한이 낮게 조작된다 |
| **allocator** | `expandable_segments` 는 이 환경 PyTorch 1.x 에서 **CUDA 초기화가 통째로 실패**한다. `max_split_size_mb:512` 만 쓴다 |
| **NaN** | 파이썬에서 `float('nan')` 은 **참**이다. `x or 0` 로 쓰면 값이 없을 때도 점수가 붙는다 |
| **`.loc` 열 이름** | `Series.loc` 은 인덱서다. `loc` 이라는 열을 만들면 `r.loc` 이 값이 아니라 인덱서를 준다 |
| **좌표계** | HO396 로커스 태그는 일부 단백질에만 있다. QJZ 번호로 통일한다 |
| **`--search-type 2`** | `alnlen` 이 **염기 단위**다. 커버리지는 `qstart/qend` 로 낸다 |
| **Foldseek** | `--exact-tmscore 1` 을 안 붙이면 TM 이 1.0 을 넘는 근사값이 나온다 |
| **프레임 혼용** | 결정 프레임 값과 AFDB 프레임 값을 한 문장에 섞지 않는다 |
| **λDE3** | BL21(DE3) 는 λDE3 를 삽입해 만든 균주다. 그 구간이 MG1655 에 없는 것은 **설계**이지 발견이 아니다 |

---

## 5. 현재 결론 (발표에 쓸 수 있는 것)

### 대조군 6개 중 셋이 떨어졌다 — 후보는 통과한 셋에서만 나왔다

| 방법 | 대조군 | 실측 | |
|---|---|---|---|
| RF2-PPI | CooC1–ChCODH2 | **0.266** (기준 0.70) | ✗ |
| Boltz ipTM | CooC1 | **0.300** (후보가 더 높음) | ✗ |
| Folddisco | CooC1 질의 → Y19 CooC | **못 잡음** (`-d` 4.0 까지) | ✗ |
| **Boltz Ni 좌표** | 직교체 간 일치 | **Δ 0.01–0.08 Å** | ✅ |
| **Foldseek** | CooC1 → Y19 CooC | **tm 0.913 / 0.837** | ✅ |
| **mmseqs 부재 판정** | fident 1.000 직교체 14건 | **진짜 부재 0** | ✅ |

### Foldseek — 명확한 음성 결과

| Strain | Proteome | Total aln | ≥0.4 | **≥0.5** | per proteome | Top TM |
|---|---|---|---|---|---|---|
| BL21(DE3) | 4,110 | 899 | 46 | **8** | 0.195% | 0.708 |
| K-12 MG1655 | 4,300 | 926 | 44 | **7** | 0.163% | 0.731 |
| Y19 | 5,325 | 1,035 | **64** | **15** | **0.282%** | **0.837** |

- **BL21 8개 중 7개는 MG1655 에도 온전히 있다** (6개는 fident ≥ 0.99)
- **남은 하나 `QJZ13931.1` BcsQ 는 MG1655 에서 pseudogene 이다** — 6번 코돈에
  조기 종결(UniProt P37655, PMID 24097954). AFDB 모델도 **242 aa 단편**으로 확인했다
  (온전한 단백질 P0DP92 는 250 aa). 기능형으로 세면 **BL21 8 / MG1655 6**
- **BL21 에 진짜 CooC 은 없다** (최고 0.708, Y19 는 0.837)
- **Y19 만 1.4–1.9배 많다.** `≥0.4` 와 `≥0.5` 양쪽에서 같은 방향이라 컷오프가
  결론을 만든 것이 아니다
- → CooC 상동체로는 BL21 lysate 활성이 설명되지 않는다. **HypA 계열 쪽이 남는다**

### 1순위 `QJZ12568.1` (133 aa, hypothetical HO396_09830)

독립적인 **세** 측정이 겹친다.

| 축 | 값 |
|---|---|
| MG1655 프로테옴·게놈 둘 다 부재 | mmseqs + `--search-type 2` |
| Ni 배위 grade **A** | Cys4, 2.23–2.40 Å |
| 폴드 | **HypA** tm 0.664 / 0.680 (※ exact 재실행 필요) |

그리고 입구를 완전히 바꾼 새 파이프라인에서도 **score 10 으로 1위**다 (2위 8점).
`mobile`(프로파지) 구간 밖이다.

> **Cys4 재질의(folddisco) 히트는 근거로 세지 않는다.** 그 질의의 템플릿이
> 이 단백질 자신이라 자기참조다.

---

## 6. 아직 안 끝난 것

- CELL 54 이량체 스크린 (63시간)
- Boltz `inputs_seqmotif` 400건 → 채점
- Boltz `inputs_cooclike` 이량체 (아직 안 띄움)
- F11 exact 재실행
- `QJZ14223.1` LysR — 최고 MG1655 히트 fident 0.345 (직교체 아닌 먼 파라로그)
- `AHZ96932` (CooJ 추정) — Y19 FASTA 와 xlsx 양쪽에 없음. NCBI 에서 받아야 함
