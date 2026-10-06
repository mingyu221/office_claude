# 도킹 의뢰 — ChCODH2 상동체 321개

| | |
|---|---|
| 작성일 | 2026-09-24 (2차 개정) |
| 워크스테이션 | `epel-af2` (같은 장비) |
| 상태 | 구조 확보·조립 완료. **도킹만 하면 됨** |
| 전체 맥락 | `handoff_homolog_selection.md` (필요할 때만) |

```
BASE = /mnt/af2results/mingyu/workspace/docking_simulation
```

> ## ⚠ 1차 판에서 바뀐 것 — 먼저 읽을 것
>
> 앞서 `center_dev`(41등가 CA 가 REF 에서 떨어진 거리)가 큰 구조 4건을
> **문제로 보고해 주셨는데, 그 판단 기준을 철회했습니다.**
>
> `center_dev` 가 큰 이유는 결합부위가 ChCODH2 와 많이 다르기 때문이고
> **그게 바로 이 연구가 예측하려는 것**입니다. 그걸로 표본을 자르면 분포의
> 꼬리만 골라서 빠집니다. 넷 다 전역 중첩은 RMSD 2.03 이하로 멀쩡했습니다.
>
> **이제 `center_dev` 는 이상 신호가 아니라 기록할 값입니다.**
> 목록과 예상값은 4절에 있습니다. 그 4건도 똑같이 도킹해 주시면 됩니다.
>
> 검산 문턱도 바뀌었습니다 — 예전 `REF 에서 1 Å 안쪽` 은 더 이상 기준이
> 아닙니다. 5절 참조.

---

## 1. 무엇을

**조립된 이량체 321개**에 ethyl viologen 을 도킹하고 점수를 돌려준다.

| | 경로 | 개수 |
|---|---|---|
| 대상 목록 | `$BASE/result/homolog_split/dock_targets.tsv` | 321 |
| 구조 | `$BASE/input/assembled/{accession}_dimer.pdb` | 321 |
| 리간드 | `$BASE/input/Veri_Ca_diameter/ev_ds.pdbqt` | 1 |

`dock_targets.tsv` 의 `tier` 컬럼이 core 212 / expansion 109 를 나눈다.
**도킹은 둘을 구분하지 않는다. 전부 같은 조건으로 돌린다.**

---

## 2. 조건 — 하나도 바꾸지 말 것

```
중심   각 구조의 41번 등가 잔기 CA      ← 고정 좌표 아님. 3절 참조
박스   npts 80 × spacing 0.375 = 반변 15.000 Å
리간드 ev_ds.pdbqt 를 그대로 복사
실행   --nrun 50 --seed 1 --heuristics 0 --autostop 0
지표   med (50런 중앙값)
```

**실행 파일**

```
/mnt/af2results/mingyu/soft/docking/ADFRsuite/bin/prepare_receptor
/mnt/af2results/mingyu/soft/docking/ADFRsuite/bin/autogrid4
/mnt/af2results/mingyu/soft/docking/AutoDock-GPU/bin/autodock_gpu_128wi
```

**GPF 는 수기 작성** (생성기 금지 — 6절 함정)

```
npts 80 80 80
gridfld {tag}.maps.fld
spacing 0.375
receptor_types A C HD N NA OA SA
ligand_types A C N
receptor {tag}.pdbqt
gridcenter {x} {y} {z}
smooth 0.5
map {tag}.A.map
map {tag}.C.map
map {tag}.N.map
elecmap {tag}.e.map
dsolvmap {tag}.d.map
dielectric -0.1465
```

`map` 줄 순서는 `ligand_types` 순서와 같아야 한다.
`parameter_file` 줄은 **없어야** 한다.

---

## 3. 박스 중심 — `dock_targets.tsv` 에 이미 들어 있다

상동체마다 41번 등가 잔기의 **잔기번호가 다르다.** PDB 에서 `resSeq == 41` 을
찾으면 엉뚱한 잔기가 잡힌다.

**이번 판부터 `dock_targets.tsv` 가 직접 `f41_pos` 를 싣는다.**
`master_1827.tsv` 를 따로 열 필요가 없다.

| 컬럼 | 내용 |
|---|---|
| `f41_pos` | **박스 중심으로 쓸 chain A 잔기번호** |
| `resname` | 그 자리의 잔기 이름 (PHE 256 · TYR 38 · 그 외 27) |
| `center_dev` | 그 CA 와 1SU6 원점의 거리 Å — **검산용이자 공변량** |
| `gxg` | 40·42 등가가 둘 다 GLY 인가 |

```python
import pandas as pd, numpy as np
from pathlib import Path

BASE = Path("/mnt/af2results/mingyu/workspace/docking_simulation")
D    = pd.read_csv(BASE/"result/homolog_split/dock_targets.tsv", sep="\t")
POS  = dict(zip(D.accession, D.f41_pos.astype(int)))
REF  = np.array([48.884, -6.486, 31.065])        # 1SU6 chain A PHE41 CA

def center(acc):
    p = next((BASE/"input/assembled").glob(f"{acc}*_dimer.pdb"))
    want = POS[acc]
    for l in p.read_text().splitlines():
        if (l.startswith("ATOM") and l[21] == "A"
                and l[12:16].strip() == "CA" and int(l[22:26]) == want):
            c = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
            return c, l[17:20].strip(), float(np.linalg.norm(c - REF))
    return None, None, None
```

**ChCODH2 자기대조 구조** (검증에 쓴 것, 잔기번호가 실제로 41)
`$BASE/input/assembled/Q9F8A8_ChCODH2_assembled.pdb`

---

## 4. `center_dev` 가 큰 4건 — 정상이다

```
accession      tier       cband  f41_pos resname  RMSD  center_dev
A0A7J2S3E9     expansion  30-40       43     ASP  2.03        9.69
A0AC61KYD2     core       30-40       40     ASP  1.90        9.55
A0A484HNT9     core       40-50       99     PHE  1.94        8.01
A0A2N7PLW0     core       30-40      130     PHE  1.85        5.25
```

**넷 다 그대로 도킹한다. 문제로 보고하지 않아도 된다.**

박스 중심 규칙은 **각 구조 자신의 41등가 CA** 다. 1SU6 고정 좌표가 아니다.
그 상동체의 축이 실제로 5~10 Å 이동해 있으면 박스도 같이 옮겨가므로,
박스는 언제나 그 단백질 자신의 부위를 반변 15 Å 로 온전히 덮는다.
**포켓이 잘리거나 무관한 표면이 끌려 들어오지 않는다.**

참고로 한 잔기 안에서도 CA 와 실제 스태킹 고리가 3.77 Å 떨어져 있다
(1SU6 PHE41 실측). `center_dev 5.25` 는 그 1.4배 수준이다.

다만 이 4건은 **점수가 이상치로 나오면 가장 먼저 볼 대상**이므로,
결과에 `center_dev` 를 그대로 실어 주면 이쪽에서 회귀로 확인한다.

**따로 알려줄 것 한 건** — `A0A1V4SNR3` 은 임계충돌이 10 으로 다른 구조
(대부분 0~3)보다 크다. RMSD 1.29 로 중첩 자체는 정상이라 **제외하지 않고
도킹한다.** 태그 목록: `$BASE/result/homolog_split/tagged_highclash.tsv` (5건)

---

## 5. 검산 — 문턱이 바뀌었다

예전 판의 `REF 에서 1 Å 안쪽` 은 **더 이상 기준이 아니다.**
실측 분포는 이렇다.

```
center_dev (Å)      개수
  0.0 ~ 1.0          295
  1.0 ~ 1.5           20
  1.5 ~ 2.0            2
  2.0 ~ 5.0            0
  5.0 초과             4      ← 4절의 네 건
```

**확인할 것은 하나다** — 계산한 중심의 `center_dev` 가
`dock_targets.tsv` 의 `center_dev` 값과 **일치하는가**.
일치하면 올바른 잔기를 잡은 것이다. 값 자체의 크기는 판단 기준이 아니다.

잔기 이름은 PHE 256 · TYR 38 이고 **그 외 27건은 비방향족**이다
(ILE 18 · HIS 3 · MET 3 · ASP 2 · LEU 1). **정상이다.**
41 자리가 방향족이냐 아니냐가 이 연구의 종속변수이므로,
비방향족이라고 해서 잘못 잡은 것이 아니다.

---

## 6. 함정 — 전부 실제로 겪은 것

| 하면 안 되는 것 | 결과 |
|---|---|
| 리간드를 obabel / prepare_ligand4 로 **재생성** | 전하합 +2 → **0.000**. 출력 숫자는 정상 범위라 알아채기 어렵고 신호가 통째로 사라진다 |
| GPF **생성기** 사용 | `parameter_file None` 줄을 뱉어 autogrid4 가 SIGSEGV |
| `--heuristics 0 --autostop 0` 누락 | 런 수가 가변이 되어 50런 고정과 시드 재현성이 동시에 깨진다 |
| **npts 108** 사용 | 반변 18 이상에서 위치내 상관이 0.93 → 0.78 로 붕괴. 80 이어야 한다 |
| PDB 에서 `resSeq == 41` 로 중심 찾기 | 상동체는 잔기번호가 다르다. `f41_pos` 를 쓸 것 |
| 산출물 존재만으로 skip | 소스 구조가 바뀐 경우를 놓친다 |

**리간드 검산** (돌리기 전에 한 번)

```
원자 16개   전하합 1.998   타입 A C N   TORSDOF 3
```

**맵 관리** — 구조당 약 27 MB. 321개를 남기면 8 GB 가 넘는다.
스크래치에서 돌리고 **도킹 직후 맵을 지울 것.** `.dlg` 만 보관.

---

## 7. 이미 검증된 것 — 다시 하지 않아도 된다

| | med | 차이 |
|---|---|---|
| 288 run Wildtype 재도킹 | **−4.380** | 정본 −4.38 과 **0.000** |
| 조립한 ChCODH2 이량체 | **−4.30** | **+0.080** (0.9×IQR) |

파이프라인이 288 run 과 일치하고, 조립 효과는 결정구조 자신의 50런 범위(0.11)보다
작다. **조립 절차는 도킹 점수를 왜곡하지 않는다.**

기록: `$BASE/result/validation/validation_stage1.tsv`
정본: `$BASE/result/Veri_Ca_diameter/engines_288_CA_r15.csv` (Wildtype `ad_med_ds` = −4.38, IQR 0.09)

**조립이 왜 필요했는지** — 1SU6 실측 결과 D클러스터를 배위하는 Cys 가
반대편 protomer 에서도 10~12 Å 안으로 들어온다. 단량체로는 결합 자리가
완성되지 않는다.

---

## 8. 돌려줄 것

```
$BASE/result/homolog_split/dock_scores.tsv
```

| 컬럼 | 내용 |
|---|---|
| `accession` | |
| `ad_n_ds` | 런 수 — **전부 50 이어야 함** |
| `ad_best_ds` | 50런 최저 |
| `ad_med_ds` | **50런 중앙값 — 주 지표** |
| `ad_iqr_ds` | 사분위 범위 |
| `gridcenter` | 실제로 쓴 좌표 (검산용) |
| `center_dev` | 그 좌표와 REF 의 거리 Å |

`ad_n_ds` 가 50 이 아닌 행이 있으면 그것만 따로 알려줄 것. `--heuristics 0
--autostop 0` 이 빠졌다는 신호다.

---

## 9. 그 다음에 할 일 (참고)

도킹이 끝나면 이쪽에서 두 가지를 한다. **도킹 세션이 할 일은 아니지만,
왜 공변량이 필요한지 알아두면 좋다.**

**노이즈 바닥** — expansion 109개는 학습 클러스터의 2번째 멤버다. 대표와 짝지어
`|score(대표) − score(2번째)|` 를 내면 그게 상동체 도킹의 잡음 크기다.
288 세트는 골격을 공유해 이 잡음원이 없었으므로 그쪽 SNR 을 그대로 쓸 수 없다.

**공변량 회귀**

```
score ~ cband + RMSD + center_dev + struct_source + mode + 임계충돌 + pLDDT
```

| 공변량 | 왜 |
|---|---|
| **`RMSD`** | **필수.** cband 와 교락 — 아래 |
| **`center_dev`** | **필수.** 〃 |
| `struct_source` | 321개 중 23개(7%)가 ColabFold. 나머지는 AFDB |
| `mode` | 조립이 global 인지 local 인지 |
| `임계충돌` | 결합 부위 근처 원자 겹침 수 |
| `pLDDT` | 구조 예측 품질 |

**교락이 있다.** 유사도가 낮을수록 조립이 나빠진다.

```
cband   RMSD   충돌   center_dev
≥60     0.61    1.4      0.14
<30     2.51   22.9      0.65
```

앞의 둘을 보정하지 않으면 `cband` 계수를 읽을 수 없다.

전부 여기 있다 — `$BASE/result/homolog_split/assembly_manifest_all.tsv`
결합부위 원자료 — `$BASE/result/homolog_split/site_check.tsv`
