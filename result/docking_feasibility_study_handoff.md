# 도킹 실행 방법 인계 — CODH/EV feasibility study

| | |
|---|---|
| **작성일** | 2026-09-02 |
| **문서 버전** | 1.0 |
| **선행 문서** | `HANDOFF_docking_simulation.md` (2026-08-31) — 172개 시점. 이 문서가 방법 부분을 대체한다 |
| **도킹 실행일** | 2026-09-02 (288 구조 × 2 엔진 × 2 리간드, 신규 57,600 런) |
| **기존 172개 배치** | 2026-08-26 (`prep/*.dlg`, `RX_*_out.sd` — 보존 중) |
| **실측 Km 기준** | 연구실 측정 21건 + Kim *et al.* Nat Commun 2024;15:2732 Table 1 6건 |
| **작업 환경** | epel-af2, conda env `dock`, RTX 3090 × 2 / 64 core |

**이 문서의 범위는 "어떻게 돌리는가"이다.** 결과 해석과 엔진 판정은
`RESULTS_288.md` 및 보고서 페이지 참조.

문서와 실물이 어긋난 사례가 이 프로젝트에서 이미 여러 번 나왔다. 아래 값은 전부
**산출물 파일에서 직접 회수해 재현 검증까지 마친 것**이며, 추측으로 채운 항목은
없다. 새 세션은 이 값들을 재조사하지 말고 그대로 쓸 것.

---

## 0. 경로 규약

```
BASE = /mnt/af2results/mingyu/workspace/docking_simulation
  input/1su6_varient/{위치}/        수용체 구조 288개  ← 기준 트리
  input/veri_engines/               RxDock 작업공간 (평면, bare filename 필수)
  input/veri_engines/prep/          수용체 전처리 + 구 172개 AutoDock 산출물
  input/veri_engines/prep/ds/       신규 AutoDock 결과 (DS 리간드)
  input/veri_engines/prep/lig/      신규 AutoDock 결과 (참조 리간드)
  input/ev_source/                  EV 리간드 원본
  input/variants/                   구 172개 트리 (좌표 대조용, 실행에는 미사용)
  result/veri_engines/              집계 CSV + 문서 + figs/
```

`input/veri_engines/` 가 평면이어야 하는 이유: RxDock의 `GetDataFileName` 이 cwd를
먼저 본다. 수용체 mol2, 리간드 sd, 점수함수 JSON 23개가 모두 이 디렉터리에 있어야
하고 커맨드에는 **경로 없이 파일명만** 넘긴다.

---

## 1. 도구 (절대경로 · 검증됨)

| 도구 | 경로 / 버전 |
|---|---|
| `autodock_gpu_128wi` | `/mnt/af2results/mingyu/soft/docking/AutoDock-GPU/bin/` |
| `autogrid4` | `/mnt/af2results/mingyu/soft/docking/ADFRsuite/bin/` |
| `prepare_receptor` | `/mnt/af2results/mingyu/soft/docking/ADFRsuite/bin/` |
| `rxcmd` | conda env `dock`. `RXDOCK_ROOT=/mnt/.../soft/docking/rxdock` |
| `pdb2pqr` | 3.6.1 (propka 3.5.1 동반) |
| `obabel` | conda env `dock` |
| `mk_prepare_ligand.py` | Meeko 0.7.1, conda env `dock` |
| GPU | RTX 3090 × 2. `--devnum` 은 **1-based** |

**"AutoDock4"는 점수함수 이름이고 실행 바이너리는 AutoDock-GPU다.** 둘은 같은 AD4
scoring function을 쓰고 `autogrid4` 산출 `.map`/`.fld` 를 공유하며 출력도 동일한
`.dlg` 이다. 엔진 정체를 재검증하지 말 것.

---

## 2. 입력 자산

### 수용체 288개
15개 위치(G40 G42 E43 T44 L46 R57 I58 N59 P60 F61 G62 D63 E64 P65 K66) × 19 치환
= 285, + `N59N_A559W`, + `A559W`(배경), + `Wildtype` = 288.

- 라벨 무결성 전수 검증 통과 — 파일명이 주장하는 잔기와 좌표가 288/288 일치,
  559번은 `Wildtype`(ALA)을 제외하고 전량 TRP.
- 그리드 중심 = F41 CE2. **288개 전부 중심과의 거리 0.000 Å** (구조들이 완전히
  중첩되어 있음). 15개 위치 모두 30 Å 큐브 안 (최대 |dz| ≈ 10.6 Å).
- 구 트리 `input/variants/` 와 겹치는 172개 중 **171개는 좌표 동일, `A559W` 만
  상이**. 그래서 A559W는 새 트리 기준으로 재전처리했다.
- 과거 `K66Y.pdb` / `K66W.pdb` 가 `_A559W` 접미사 없이 내보내져 있었다. 559번이
  TRP임을 확인한 뒤 `K66Y_A559W.pdb` / `K66W_A559W.pdb` 로 개명 완료.

### 리간드 2종
둘 다 원자 16 / 형식전하 +2 / 회전결합 3 이며 **시작 컨포머만 다르다.**

| 파일 | 출처 | H | 고리간 이면각 |
|---|---|---|---|
| `ev_ds.sd` / `ev_ds.pdbqt` | `ev_source/oxEV_prepared.sd` — **C-Docker가 쓴 것** | 18 | −30.0° |
| `ev_lig.sd` / `ev_lig.pdbqt` | `ev_source/ev_ref.sd` — 결정구조 참조 컨포머 | 0 | +2.02° (평면) |

두 PDBQT 모두 `(16 atoms, charge 1.998, types A C N, TORSDOF 3)`.

EV = diethyl viologen C14H18N2(2+). 18개 H는 전량 C–H 비극성이라 AutoDock의
united-atom 관례상 PDBQT에서 병합되어 사라진다. **H 0개는 결함이 아니라 정답.**

---

## 3. 수용체 전처리

구조당 약 30초. 리간드와 무관하므로 한 번만 만들면 두 리간드가 공유한다.

```bash
# 1) HET 제거 — ATOM 레코드만 보존 (SF4 / FES / NFS / HOH 전량 제거)
grep '^ATOM' {src}.pdb > {tag}_strip.pdb ; echo END >> {tag}_strip.pdb

# 2) 프로톤화
pdb2pqr --ff=AMBER --with-ph=8.0 --titration-state-method=propka \
        --pdb-output {tag}_ph8.pdb  {tag}_strip.pdb  {tag}.pqr

# 3) RxDock 용 (cwd = input/veri_engines)
obabel {tag}.pqr -O {tag}.mol2 --partialcharge gasteiger

# 4) AutoDock 용 (cwd = input/veri_engines/prep)
prepare_receptor -r {tag}_ph8.pdb -o {tag}.pdbqt -U nphs -A None
```

`prepare_receptor` 의 플래그는 **결과에 영향이 없다** — `-U nphs -A None` / `-A None`
/ 기본값 세 후보가 전부 동일한 11,388 원자 PDBQT를 냈다. 인계문서 표기를 따라
`-U nphs -A None` 을 쓴다.

### 낡은 전처리 탐지 (중요)
"산출물이 있으면 skip" 만으로 판단하면 **소스 구조가 바뀐 경우를 놓친다.** 실제로
A559W가 이 방식으로 조용히 건너뛰어져 구 좌표 기준 전처리가 새 배치에 섞일 뻔했다.
반드시 내용을 대조할 것:

```python
md5(ATOM 레코드[12:54] of {tag}_strip.pdb) == md5(ATOM 레코드[12:54] of {src}.pdb)
```

---

## 4. 리간드 준비 — **Meeko 외 대체 불가**

```bash
cp oxEV_prepared.sd  ev_ds.sdf          # .sd 확장자를 Meeko가 거부한다
mk_prepare_ligand.py -i ev_ds.sdf -o ev_ds.pdbqt
```

- Meeko는 **명시적 H를 요구한다.** H가 없는 소스는 `obabel -h` 로 먼저 붙일 것.
- **obabel 이나 prepare_ligand4 로 리간드를 재생성하면 PDBQT 전하합이 +2가 아니라
  0.000 이 된다.** 점수의 대부분을 전하가 설명하므로 신호가 통째로 사라지는데,
  출력 숫자는 정상 범위라 알아채기 어렵다. 350 homolog 확장 시에도 기존
  `ev_*.pdbqt` / `ev_*.sd` 를 **그대로 재사용**할 것.
- RxDock 용 `.sd` 는 MDL 원자블록의 **legacy 전하열**(37–39 컬럼)을 읽는다. `M CHG`
  가 아니다. 코드맵 `{+1:3, +2:2, +3:1, −1:5, −2:6, −3:7}`.
  `oxEV_prepared.sd` 는 이미 N 두 개에 코드 3 을 갖고 있어 그대로 복사하면 된다.

---

## 5. AutoDock-GPU

### GPF (수기 작성 — 생성기 사용 금지)
ADFRsuite의 GPF 생성기는 `parameter_file None` 이라는 줄을 뱉어 `autogrid4` 를
SIGSEGV 로 죽인다.

```
npts 80 80 80
gridfld {tag}.maps.fld
spacing 0.375
receptor_types A C HD N NA OA SA
ligand_types A C N
receptor {tag}.pdbqt
gridcenter 3.266 -25.350 -3.469
smooth 0.5
map {tag}.A.map
map {tag}.C.map
map {tag}.N.map
elecmap {tag}.e.map
dsolvmap {tag}.d.map
dielectric -0.1465
```

`map` 줄의 순서는 `ligand_types` 순서와 일치해야 한다.

### 실행
```bash
autogrid4 -p {tag}.gpf -l {tag}.glg          # 약 2초

autodock_gpu_128wi --ffile {tag}.maps.fld --lfile {ligand}.pdbqt \
                   --resnam {tag} --nrun 50 --seed 1 \
                   --heuristics 0 --autostop 0 --devnum {1|2}    # 약 3초
```

**`--heuristics 0 --autostop 0` 은 빼면 안 된다.** 빠지면 AutoDock-GPU가 런 수를
스스로 줄이거나 수렴 시 조기 종료해서, 50런 고정과 seed 재현성이 동시에 깨진다.
변이체마다 표본 수가 달라지면 best/median 비교가 성립하지 않는다.

### 그리드 공유와 정리
`ligand_types` 가 두 리간드에서 동일하므로 **그리드는 구조당 한 번만 만들어 두
리간드에 쓴다.** 맵은 구조당 약 27 MB이므로 도킹 직후 삭제해 디스크를 평평하게
유지한다 (288개를 남기면 약 8 GB).

---

## 6. RxDock

### 수용체 파라미터 JSON — `RX_{tag}.json`
```json
{
  "media-type": "application/vnd.rxdock.parameters",
  "title": "RX_{tag}",
  "version": "0.1.0",
  "sections": ["receptor", "mapper", "cavity"],
  "receptor": {"file": "{tag}.mol2"},
  "mapper": {
    "site-mapper": "SphereSiteMapper",
    "center": "(3.266,-25.350,-3.469)",
    "radius": 15.0,
    "small-sphere-radius": 1.5,
    "large-sphere-radius": 6.0,
    "excluded-volume-radius-increment": 0,
    "grid-step": 0.5,
    "minimum-cavity-volume": 50,
    "maximum-cavities": 1
  },
  "cavity": {"scoring-function": "CavityGridSF", "weight": 1}
}
```

- `sections` 키가 없으면 동작하지 않는다.
- `center` 는 **소괄호 문자열** `"(x,y,z)"` 여야 한다. 배열이 아니다.

### 실행
```bash
# cwd = input/veri_engines  (평면 · bare filename)
rxcmd cavity-search -r RX_{tag}.json -W                        # 캐비티 검출

rxcmd dock -i {ligand}.sd -o {out}_out.sd \
           -r RX_{tag}.json -p dock.json -n 50 -s 1            # 약 12초
```

- 서브커맨드는 `cavity-search` / `dock` 이다. `rbcavity` / `rbdock` 이 아니고
  `-was` 도 아니다.
- **캐비티는 리간드와 무관**하므로 한 번 검출하면 두 리간드가 공유한다.
- **`-H` 는 쓰지 않는다.** 기본값이 "극성 H만 읽음"이라, H 18개를 가진 DS 리간드와
  H 0개인 참조 리간드가 같은 중원자 표현으로 들어간다. 두 리간드 비교에서 차이를
  컨포머 하나로 좁히려면 이 기본값을 유지해야 한다.

### 15 Å의 의미 (혼동 주의)
AutoDock의 30 Å는 **정육면체 반변**이고 288개 전부 동일 고정이다. RxDock의 15 Å는
**캐비티를 찾을 구의 반경**이며, 실제 탐색공간은 검출된 캐비티다. 구조마다 달라져
288개에서 371–1,723 Å³ (4.6배) 변동했다. **`cav_vol` 을 반드시 함께 기록할 것.**

```python
cav_vol = len(json["docking-site"]["cavities"][0]["coordinates"]) * 0.5**3
```

---

## 7. 파싱 규약

```python
# AutoDock .dlg → 50개 값
re.compile(r"Estimated Free Energy of Binding\s*=\s*([-\d.]+)")

# AutoDock 포즈 좌표
"DOCKED: ATOM" 으로 시작하는 줄에서 앞 8자를 떼면 PDBQT ATOM 줄이 된다.
16원자씩 끊으면 런 순서대로 50개 포즈, 위 정규식의 값 순서와 일치한다.

# RxDock _out.sd → 50개 값
SD 데이터 필드명은  "> <rxdock.score>"  이다.  "SCORE" 가 아니다.
(rxdock.score.inter 등 하위 필드가 많으므로 정확히 매칭할 것)
```

---

## 8. 배치 전 재현 검증 (반드시 수행)

커맨드를 문서에서 읽지 말고 **기존 산출물을 재현해서 확정한다.** `D63A_A559W` 로
검증했고 두 엔진 모두 완전 일치했다.

| 단계 | 검증 방법 | 결과 |
|---|---|---|
| strip | 기존 `_strip.pdb` 와 md5 대조 | 일치 |
| pdb2pqr | 기존 `.pqr` / `_ph8.pdb` 와 좌표 md5 대조 | 일치 |
| prepare_receptor | 후보 플래그별 생성 후 기존 `.pdbqt` 와 대조 | 3후보 전부 일치 |
| AutoDock | 구 리간드로 재실행 → 기존 `.dlg` 의 50개 에너지 정렬 대조 | **전량 일치** |
| RxDock | 구 리간드로 재실행 → 기존 `_out.sd` 의 50개 점수 정렬 대조 | **전량 일치** |

새 환경·새 노드에서 시작할 때 이 절차를 먼저 돌릴 것. 재현되지 않으면 배치를
시작하지 말 것.

---

## 9. 함정 목록

1. **리간드를 obabel/prepare_ligand4 로 재생성** → 전하합 0, 신호 소실, 발견 어려움
2. **`--heuristics 0 --autostop 0` 누락** → 런 수 가변, 시드 재현성 붕괴
3. **Meeko에 `.sd` 확장자 전달** → `Format [sd] not in supported formats` (`.sdf` 로)
4. **Meeko에 H 없는 분자 전달** → `RDKit molecule has implicit Hs`
5. **SD 점수 필드를 `SCORE` 로 파싱** → 빈 리스트, 조용히 실패
6. **산출물 존재만으로 skip** → 소스가 바뀐 구조를 건너뜀 (A559W 사례)
7. **GPF 생성기 사용** → `parameter_file None` 로 autogrid4 SIGSEGV
8. **RxDock JSON에 `sections` 키 누락 / `center` 를 배열로** → 동작 안 함
9. **RxDock을 하위 디렉터리에서 실행** → 데이터 파일 탐색 실패 (cwd 평면 유지)
10. **`.xlsx` 읽기 전 `openpyxl` 미설치** → ImportError

---

## 10. 산출물과 보존 정책

| 산출물 | 개수 | 보존 |
|---|---|---|
| `prep/ds/*.dlg`, `prep/lig/*.dlg` | 288 × 2 | **재실행 없이 복원 불가 — 보존** |
| `RXDS_*_out.sd`, `RXLG_*_out.sd` | 288 × 2 | **보존** |
| `prep/*.dlg`, `RX_*_out.sd` (구 172/171) | — | **보존** (신규분과 이름이 다르므로 충돌 없음) |
| `RX_*-docking-site.json` | 288, 약 6 MB/개 | 재생성 가능 — 공간 부족 시 1순위 삭제 |
| `*.map`, `*.fld` | — | 도킹 직후 삭제 (GPF로 재생성) |
| `{tag}.pdbqt`, `.mol2`, `.pqr`, `_ph8.pdb`, `_strip.pdb` | 288씩 | 재생성 가능 (구조당 30초) |

집계: `result/veri_engines/` 의 `manifest.csv`, `engines_288.csv`,
`analysis_288.csv`, `km_merge_288.csv`, `km_merge_paper.csv`, `summary_288.csv`,
그림은 `figs/`.

---

## 11. 소요 시간 (실측)

| 단계 | 구조당 | 288개 배치 |
|---|---|---|
| 수용체 전처리 | 약 30초 | 3.3분 (12 워커) |
| autogrid4 | 약 2초 | — |
| AutoDock-GPU 도킹 | 약 3초 | 약 15분 (16 워커, 리간드 2종 = 576런) |
| RxDock cavity+dock | 약 12초 | 6.9분 (24 워커, 리간드 2종) |

전체 재실행은 30분 이내다. **이 비용이 낮다는 사실이 의사결정에 중요하다** —
"돌려야 하나"를 오래 논의하는 것보다 그냥 돌리는 편이 싸다.

---

## 12. 미확정으로 남은 것

- **`ev_lig` 가 결정구조 참조 컨포머라는 가설.** centroid가 결합부위 실좌표이고
  파일명이 `ev_ref.sd` 인 것이 근거이나 **결정구조와 직접 대조하지 않았다.**
  사실이면 재도킹 RMSD 평가가 부풀려졌을 수 있다. 8X9F의 EV(잔기명 `S8I`)와
  centroid·이면각을 대조하면 즉시 판정된다.
- **산출 포즈의 실제 이면각 수렴값.** 리간드 교체 실험에서 AutoDock 순위상관이
  0.94로 나와 "시작 컨포머는 잊혀진다"가 간접 확인됐으나, 포즈의 이면각 분포를
  직접 잰 적은 없다. `prep/*/*.dlg` 의 DOCKED 좌표에서 **재실행 없이** 측정 가능.
  고리 대칭 때문에 0~90°로 접어서(fold) 봐야 한다.
- **`best` 와 `med` 의 우열.** n=27에서 리간드에 따라 순서가 뒤집힌다. 이 표본으로
  주 지표를 고르는 것은 과적합이다.

---

## 13. 다음 단계 (350 homolog)

동일 절차를 그대로 적용하되 세 가지가 달라진다.

1. **그리드 중심을 구조마다 정해야 한다.** 288개는 F41 CE2가 좌표까지 동일했지만
   이종 homolog는 중첩되어 있지 않다. 정렬 후 대응 잔기의 CE2를 쓰거나, 구조별
   중심을 계산해야 한다. 이 결정이 전체 결과를 좌우한다.
2. **RxDock 캐비티 실패율을 기록할 것.** 288개에서 이미 4.6배 변동했고 이종
   구조에서는 악화될 가능성이 크다. `cav_vol` 과 검출 실패 건수가 3차 관문의
   판정 자료다.
3. **리간드는 손대지 말 것.** 기존 `ev_*.pdbqt` / `ev_*.sd` 를 그대로 복사해 쓴다.

350개에는 정답지(실측 Km도, C-Docker 점수도)가 없다. **3차 관문은 적용가능성
판정이지 정확도 판정이 아니다.** 정확도는 288개 + 실측 27건에서 이미 결정됐다.
