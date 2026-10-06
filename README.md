# codhqsar — ChCODH2 상동체 EV 결합의 3D-field QSAR

CODH/ethyl-viologen 도킹 파이프라인에 **CoMFA/CoMSIA 방식의 3차원 장 기술자**를
붙이고, 그 결과를 상동체 순위 → Km 눈금으로 잇는 패키지.

설계 근거와 해석 규칙은 **[`docs/design_qsar_comfa.md`](docs/design_qsar_comfa.md)**
에 있다. 먼저 읽을 것.

---

## 왜 이것이 되는가

CoMFA는 보통 정렬에서 죽는다. 이 프로젝트는 그 조건이 도킹 때문에 **이미 전부
충족돼 있다.**

- 321개 조립 이량체가 **전부 1SU6 좌표계**에 중첩돼 있다 (6OND도 먼저 옮겨짐)
- 결합부위가 **F41 등가 CA, 반변 15 Å**으로 확정돼 있다 (`f41_pos`)
- 격자 파라미터가 검증돼 있다 (반변 18 Å에서 위치내 상관 0.93 → 0.78 붕괴)
- 라벨이 균질하다 — 같은 엔진·같은 조건·50런 중앙값

그래서 **추가 도킹 없이** 장 기술자를 만들 수 있고, 그 장은 조성 24차원이 담지
못하는 것을 담는다: 전하·부피·방향족이 **어디에** 있는가.

---

## 설치

```bash
pip install -r requirements.txt
export CODHQSAR_BASE=/mnt/af2results/mingyu/workspace/docking_simulation
```

경로는 전부 `CODHQSAR_BASE` 기준이다. `src/codhqsar/config.py`의 상수는
인계문서가 고정한 값이므로 근거 없이 바꾸지 말 것.

## 실행

```bash
python scripts/00_preflight.py --with-labels      # 무결성 · 노이즈 바닥 · 공변량
python scripts/01_build_fields.py                 # 321 × 59,058, 10~15분
python scripts/02_fit_comfa.py                    # 적합 + 기준선 + y-rand + AD
python scripts/03_contours.py                     # .dx + PyMOL
python scripts/04_score_candidates.py --candidates pool.tsv
```

민감도 (주 분석 후, 설계는 건드리지 않고):

```bash
python scripts/01_build_fields.py --origin REF        # 좌표계 선택
python scripts/01_build_fields.py --source autogrid   # 고전 CoMFA 장
python scripts/02_fit_comfa.py --residualized         # 조립 인공물 제거
```

## 테스트

```bash
python -m pytest tests/ -q      # 33 passed
```

실제 데이터는 `epel-af2`에 있어 이 저장소에서 접근할 수 없다. 그래서
`tests/test_end_to_end.py`가 **가짜 구조에 기하학적 신호를 심어** 전 과정을
돌린다 — 장 계산 → 클러스터 분할 → PLS → 홀드아웃 평가 → y-randomization →
등고선 → AD 게이트. 신호를 기술자가 아니라 기하에 심었으므로, 모델이 그것을
**표현을 통해** 찾아야 통과한다.

---

## 구조

```
src/codhqsar/
  config.py              인계문서가 고정한 상수 전부 (격자 · 노이즈 바닥 · 시드)
  labels.py              score = -ad_med_ds. 부호 누락/이중반전 가드
  tables.py              정본 테이블 리더. 결측은 조용히 넘기지 않고 실패
  splits.py              클러스터 단위 분할 · GroupKFold · 중첩 검사
  domain.py              적용 범위 — leverage · kNN (DModX는 모델 쪽)
  evaluate.py            OLS 공변량 회귀 · 노이즈 바닥 · 구간별 보고 · claimable
  contours.py            StDev*Coeff → .dx + PyMOL
  km_bridge.py           프록시 검증 (집단 분리 · dchg 편상관) + Km 보정
  pdbio.py               고정컬럼 PDB 리더 + 원자 성질표
  features/
    grid.py              격자 기하. SELF/REF 원점
    fields_comsia.py     가우시안 유사도 지수 6종 (주력)
    fields_autogrid.py   autogrid 맵 리더 + GPF 작성기 (고전 CoMFA)
    fieldset.py          행렬 조립 · min-sigma 필터 · 블록 스케일링
    composition.py       24차원 조성 — 측정된 기준선
    site.py              결합부위 잔기 기술자 (S3) + 정렬 불필요 shell 버전
  models/pls_comfa.py    PLS + 클러스터 CV 성분 선택 + y-rand + DModX
```

### 장 6종

| 장 | 프로브 | 왜 |
|---|---|---|
| `steric` | vdW 반경 | 포켓 수용 |
| `electrostatic` | 부분전하 | EV가 2가 양이온 |
| `hydrophobic` | Kyte-Doolittle | 소수성 접촉 |
| `hbd` / `hba` | 공여체/수용체 지시자 | 극성 상호작용 |
| `aromatic` | 고리 중심 | **바이피리디늄 스태킹 — 이 계의 종속변수** |

고전 CoMFA 장(LJ `A`/`C`/`N` + Coulomb + 탈용매화)은 프로젝트 자신의 autogrid
맵에서 읽는다. 맵은 도킹 직후 삭제되는 설정이므로 필요하면 `write_gpf`로
재생성한다 (**GPF 생성기는 쓰지 말 것** — `parameter_file None` 줄로 autogrid4가
SIGSEGV).

---

## 코드가 강제하는 규칙

주석이 아니라 실행 경로로 막아둔 것들.

| 규칙 | 위치 |
|---|---|
| 라벨 부호 (반전 누락·이중반전 둘 다) | `labels.assert_score_sign` |
| 노이즈 바닥 0.280 미만은 주장 불가 | `labels.resolvable` |
| 분할은 항상 클러스터 단위 (내부 선택 루프까지) | `splits.*` |
| 확장분은 평가셋 금지 | `splits.split_from_targets` → AssertionError |
| 필터·스케일러는 fold마다 재적합 | `CoMFAModel.cross_validate` |
| `cband`는 잔차화 대상이 아니다 | `evaluate.TECH_COVARIATES`에서 제외 |
| `cband` 계수는 RMSD·center_dev 없이 못 읽는다 | `covariate_regression` → KeyError |
| 양 끝 구간(2개·3개)은 개별 보고 금지 | `report_by_band`가 합침 |
| `ad_n_ds != 50` 행은 거부 | `tables.load_scores` |
| `X` 잔기 삭제 금지 | `composition` |
| FASTA 헤더 숫자는 라벨이 아니다 | 라벨은 항상 `dock_scores.tsv`에서 조인 |

---

## 상태

| | |
|---|---|
| 코드 | 작성 완료. 단위 33개 테스트 통과 |
| 실제 데이터 실행 | **미실행** — 데이터가 `epel-af2`에 있고 이 컨테이너에서 접근 불가 |
| 전제 | 도킹 321개 완료 (`dock_scores.tsv`) |

미정 파라미터(`min_sigma` · stride · radius · CoMSIA `alpha`)는
`docs/design_qsar_comfa.md` §10에 모아 뒀다. **전부 라벨을 보기 전에**
학습셋 교차검증으로 정해야 한다.
