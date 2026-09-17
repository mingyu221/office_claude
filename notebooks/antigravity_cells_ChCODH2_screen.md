# ChCODH2 Ni-insertase 스크리닝 — Antigravity 노트북 셀 모음

인계문서(`HANDOVER_ChCODH2_screen.md`, `HANDOVER_bacteriaDB.md`, `ChCODH2_Ni_insertase_in_silico_design.md`)와
원본 노트북(`RF2PPI_ChCODH2_screen.ipynb`)을 합쳐, **epel-af2 서버의 실제 경로와 이미 구축된 자산**을 반영한 셀 모음.

사용법: 아래 `CELL 00` ~ `CELL 37`을 **위에서부터 하나씩 복사해 Antigravity 노트북 셀에 붙여넣고 실행**한다.
전부 파이썬 셀이다(`%%bash` 매직 없음 — Antigravity/커널 환경에 상관없이 돌아가도록 `subprocess`로 통일).
셀 사이에 변수가 이어지므로 순서대로 실행해야 한다.

원본 노트북과 달라진 핵심 4가지:
1. `BASE`를 `/data/chcodh2_ppi` → `/mnt/af2results/mingyu/workspace/ppi_discovery` 로 변경.
   서버 workspace 관례(`seek_ni_insertase`, `ppi_discovery`)대로 **`input/` `result/` `script/`** 3층 구조
2. BL21 proteome ID를 **UP000503272**(우리 구조 DB 기준)로 통일, Y19(UP000034085) 추가
3. 프로테옴 서열·4-category 분류는 **UniProt 재다운로드 대신 기존 GenBank 자산 재사용** (ID 체계를 GenBank로 통일)
4. Folddisco 결과(tid = `AF-<UniProt>` / `cf_<UniParc>`)를 **서열 매칭 crosswalk**로 GenBank ID에 붙여 Part 8 통합

---

## CELL 00 — 작업 순서 개요 (제일 먼저 실행)

```python
# =============================================================================
# CELL 00 | 작업 순서 개요 · 블로커 체크리스트
#   - 실행해도 아무것도 바꾸지 않는다. 전체 흐름과 "지금 없는 것"을 확인만 한다.
# =============================================================================
PLAN = """
[ChCODH2 Ni-insertase in silico 스크리닝 — 실행 절차서]

■ 규칙 3가지
  1. 위에서부터 순서대로. 건너뛰면 각 셀 맨 앞의 need() 가 막고 무엇을 먼저 할지 알려준다.
  2. (BG) 표시된 셀은 "시작만" 한다. 끝나기 전에 다음 단계로 가면 안 된다.
     완료 확인 = CELL 15 를 돌려 로그에 DONE_ 이 보이는지 본다.
  3. 커널을 재시작했으면 CELL 01, 02 부터 다시. stage() 로 어디까지 됐는지 확인한다.

──────────────────────────────────────────────────────────────────────
단계 A. 준비                                                  (한 번만)
  CELL 01 → 02                설정, 자산·커널·GPU 점검
  CELL 03 → [터미널에서 bash] → CELL 04
                              RF2-PPI 설치          ★체크포인트 1
                              예제 판정 구간이 5/5 일치해야 통과

단계 B. bait 서열                                               (5분)
  CELL 05 → 06 → 07           ChCODH2 서열 입력, CooC1 추출, segment 생성
                              ★서열을 안 넣으면 여기서 멈춘다

단계 C. 프로테옴·차등 분류                             (5분 또는 1시간)
  CELL 08                     프로테옴 로드 + 구조 DB 와 ID 계열 대조
  CELL 09                     분류 TSV 가 있으면 여기서 끝
    └ [없음] 이면 → CELL 10 (BG, 30~60분) → DONE_easy_search 대기 → CELL 11
  CELL 12                     Track A / Track B 확정

단계 D. 세균 레퍼런스 DB                              (2~3시간, 최장)
  CELL 13 (BG, ~1시간)        17,992개 세균 프로테옴 다운로드 (~16GB)
    ↓ DONE_download 확인 (CELL 15)          ※ 확인 전에 14 로 가면 반쪽 DB 가 된다
  CELL 14 (BG, 1~2시간)       병합 + mmseqs DB + 인덱스
    ↓ DONE_bactdb 확인
  ※ 이 두 시간 동안 단계 F(Track B)를 병렬로 돌릴 수 있다. GPU 와 CPU 라 안 겹친다.

단계 E. Track A — RF2-PPI                        (gate 통과 시에만)
  CELL 16 (1분) → CELL 17     ★체크포인트 2 — 여기서 진행 여부가 갈린다
                              깊이 >= 500 진행 / 200~500 진행하되 replicate 5회
                              200 미만 → 단계 E 를 건너뛰고 단계 F 에 집중
  CELL 18 → CELL 19 (BG, 4~10시간) → DONE_prey 대기
  CELL 20 → 21 → CELL 22 (BG·GPU, 6~15시간) → DONE_rf2ppi 대기
  CELL 23 → 24                랭킹 + 양성대조군 판정

단계 F. Track B — Boltz-2                  (레퍼런스 DB 불필요, 병렬 가능)
  CELL 25 → [터미널에서 bash] Boltz-2 설치
  CELL 26 → CELL 27 (BG·GPU, 10~20시간) → DONE_boltz 대기 → CELL 28
  CELL 28a → 28a-2 → 28a-3   ★랭킹 해석. ipTM 은 짧은 prey 를 부풀리므로
                              길이 보정 없이 상위권을 믿으면 안 된다
  CELL 28a-4 (BG·GPU, 30분)   metal 모티프 보유 단백질에 co-folding 증거 추가
  CELL 28a-5                  그 결과를 Track B 분포 대비 백분위로 판정
  CELL 28a-6                  막단백질 교란 점검 (상위권이 수송체로 채워졌는가)
  CELL 28a-7 → 28a-8          Ni 이 실제로 어디 놓였는지 좌표로 확인
  CELL 28a-9 (BG·GPU, 1시간)   MG1655·Y19 금속 모티프 단백질 co-folding
  CELL 28a-10                 세 균주 비교 판정 (대조군이 있어야 해석된다)
                              (계면 0/348 — ligand_iptm 폐기 근거)
  CELL 28b → 28c              (선택) Foldseek-Interface 계면 대조

단계 G. Track C + 통합                                       (수십 분)
  CELL 29                     Folddisco metal/ATP 모티프 검색
  CELL 30                     구조 tid → GenBank ID crosswalk (수 분)
  CELL 31 → 32 → 33           모티프 CSV 변환 → 통합 랭킹 → 최종 리포트
──────────────────────────────────────────────────────────────────────

■ 상황별 대처
  "어디까지 했더라"          → stage()
  "백그라운드 다 됐나"        → CELL 15
  "선행조건 미충족 에러"      → 메시지가 지목한 셀을 먼저 돌린다
  "결과가 이상하다"           → 산출물 시각을 확인한다. 입력보다 오래된 결과는 무효다
  "GPU 두 장 쓰고 싶다"       → 단계 E 는 GPU_ID=1, 단계 F 는 CELL 27 의 BOLTZ_GPU=0

[디렉터리 구조 — 서버 workspace 관례(input/result/script)]
  /mnt/af2results/mingyu/workspace/ppi_discovery/
  |-- input/
  |   |-- seq/         baits.fasta, prey_trackA.fasta
  |   |-- db/          프로젝트 mmseqs DB (baitDB, preyDB)
  |   +-- external/    Mac 에서 가져온 BLASTp 분류 TSV
  |-- result/
  |   |-- search/      *.m8 (homolog 검색 결과)
  |   |-- paired/      paired MSA (.a3m)
  |   |-- rf2ppi/      input_file, replicate 로그
  |   |-- boltz/       inputs/, out/
  |   |-- folddisco/   run01(metal) / run02(atp) TSV
  |   |-- table/       ★최종 CSV (분류·랭킹·crosswalk)
  |   +-- log/         백그라운드 실행 로그
  |-- script/          생성된 .sh (설치·검색·추론)
  +-- tmp/             mmseqs 스크래치 — Part 2 끝나면 지워도 되는 곳

  ※ 세균 레퍼런스 DB는 공용이라 워크스페이스 밖에 둔다:
    /mnt/af2results/mingyu/database/reference_proteomes/Bacteria/
      fasta/               UP*_taxid.fasta.gz  17,992개 (~20GB)
      bacteria_ref.fasta   taxid 를 헤더에 박아 병합한 것
      bactDB*              mmseqs DB + 인덱스
  ※ RF2-PPI 코드는 워크스페이스가 아니라 /mnt/af2results/mingyu/RoseTTAFold2-PPI 에 설치한다
    (folddisco 바이너리와 같은 층 — 툴은 밖, 데이터는 워크스페이스 안)

[지금 없는 것 = 블로커]
  1. ChCODH2(A559W) native 서열      -> CELL 05 (필수, 실험에 쓰는 그 서열로)
  2. CooC1/2/CooT/CooJ 양성대조군 서열 -> CELL 05-06
  3. RF2-PPI 설치 + 가중치            -> CELL 03
  4. Boltz-2 설치                     -> CELL 25
  5. Bacteria reference proteome DB   -> CELL 13-14
  6. Feasibility gate 결과            -> CELL 16-17
  7. (선택) 계면 지원 foldseek + PDB 계면 대표 DB(Zenodo) -> CELL 28b

[이미 있는 것 = 재사용]
  - BL21/Y19/MG1655 프로테옴 서열 (GenBank ID)
  - BL21/Y19 vs MG1655 4-category 분류 (2171/1575/125/217)
  - BL21/Y19 구조 DB + Folddisco/Foldseek 인덱스 (협업팀 구축)
  - Folddisco 실행 바이너리 (Rust, 파이썬 환경 무관)
  - GPU RTX 3090 x2 (driver 580.173.02)

[설계문서 경고 — 결과 해석 시 반드시 반영]
  - Track A는 BL21/MG1655 점수가 사실상 같다. 균주 특이성은 gene content 분류로 잡는다.
  - CooS 계통 분포가 좁아 Track B가 주력, Track A가 보조가 될 확률이 높다.
  - insertase는 apo 중간체에 결합할 수 있다 -> segment bait 병행, ipTM 컷오프 완화 후 수동 검토.
"""
print(PLAN)


# 진행 상태 확인과 선행조건 검사 함수는 CELL 01 에 있다 (stage(), need(), already()).
# 커널 재시작 후 어디까지 됐는지 보려면 CELL 01 을 돌린 뒤 stage() 를 부른다.
```

---

## CELL 01 — CONFIG · 공용 헬퍼

```python
# =============================================================================
# CELL 01 | CONFIG — 경로/스레드/GPU + 백그라운드 실행 헬퍼
#   - 이 셀만 환경에 맞게 고치면 나머지는 그대로 실행 가능
#   - 원본 BASE=/data/chcodh2_ppi 는 서버에 없어 workspace/ppi_discovery 로 변경
# =============================================================================
import os, re, sys, json, gzip, shutil, subprocess, textwrap, glob, time
import numpy as np
from pathlib import Path

# ---------------- 작업 루트 ----------------
# 서버 workspace 관례(seek_ni_insertase, ppi_discovery)를 따라 input / result / script 3층으로 잡는다.
BASE      = Path("/mnt/af2results/mingyu/workspace/ppi_discovery")
TOOLS     = Path("/mnt/af2results/mingyu")      # folddisco 바이너리와 같은 층. 툴은 워크스페이스 밖에 둔다.
THREADS   = 24
GPU_ID    = 1          # ★GPU 0 은 디스플레이가 물려 있어 13배 느리다 (예제 47~84s vs 3.5~5s).
                       #  Track A/B 를 두 GPU 에 나눠 돌릴 때만 0 을 쓴다.
CONDA_ENV_RF2 = "rf2ppi"
CONDA_ENV_BOLTZ = "boltz"

DIR = {
    # --- input/ : 넣는 것 ---
    "seq":       BASE/"input"/"seq",         # baits.fasta, prey_trackA.fasta
    "db":        BASE/"input"/"db",          # mmseqs DB (bactDB/baitDB/preyDB) — 대용량
    "external":  BASE/"input"/"external",    # Mac 에서 가져온 BLASTp 분류 TSV
    # --- result/ : 나오는 것 ---
    "search":    BASE/"result"/"search",     # *.m8
    "paired":    BASE/"result"/"paired",     # paired MSA (.a3m)
    "rf2ppi":    BASE/"result"/"rf2ppi",     # input_file, *.log
    "boltz":     BASE/"result"/"boltz",      # inputs/, out/
    "folddisco": BASE/"result"/"folddisco",  # run01/run02 TSV
    "table":     BASE/"result"/"table",      # 최종 CSV (분류·랭킹·crosswalk)
    "log":       BASE/"result"/"log",        # 백그라운드 실행 로그
    # --- script/ : 돌리는 것 ---
    "script":    BASE/"script",              # 생성된 .sh
    # --- tmp/ : mmseqs 스크래치. Part 2 가 끝나면 언제든 지워도 되는 곳 ---
    "tmp":       BASE/"tmp",
}
for _p in DIR.values():
    _p.mkdir(parents=True, exist_ok=True)
RF2PPI_DIR = TOOLS / "RoseTTAFold2-PPI"       # 워크스페이스가 아니라 툴 디렉터리에 설치

# 공용 레퍼런스 DB — 이 프로젝트 전용이 아니므로 database/ 아래에 두고 재사용한다.
# (colabfold_db, bacteriaDB, folddisco 인덱스와 같은 층)
REFDB = TOOLS / "database" / "reference_proteomes"
RB    = REFDB / "Bacteria"      # 세균 세트. CELL 16/19 가 여기 bactDB 를 참조한다
(RB / "fasta").mkdir(parents=True, exist_ok=True)

# ---------------- 이미 구축된 자산 (인계문서 §1) ----------------
# ★ 프로테옴 FASTA 는 구조 DB 와 같은 GenBank 제출 기록에서 뽑은 것을 쓴다.
#   inhouseDB/ 가 그 목적으로 다시 받은 것(bl21 = QJZ… 계열)이고, protein_list/ 는
#   그보다 앞서 쓰던 다른 기록(QZI… 계열)이다. 후자로 돌리면 Part 8 에서 구조 DB 와
#   교집합이 0 이 된다 — CELL 02 의 "ID 공간 대조" 가 이걸 잡는다.
INHOUSE = TOOLS / "database" / "bacteriaDB" / "inhouseDB"
PLIST   = TOOLS / "database" / "protein_list"

def _pick(*cands):
    """앞의 것을 우선하되 없으면 뒤로 넘어간다."""
    for c in cands:
        if c.exists(): return c
    return cands[0]

ASSET = {
    # 프로테옴 서열 (GenBank ID — 구조 DB 와 같은 기록)
    "faa_bl21":   _pick(INHOUSE/"bl21_db_match_qjz.faa", PLIST/"bl21_de3_protein.faa"),
    "faa_y19":    _pick(INHOUSE/"y19_db_match.faa",      PLIST/"y19_protein.faa"),
    "faa_mg1655": PLIST/"mg1655_protein.faa",   # 분류의 비교 대상일 뿐 조인 키가 아니다
    "tsv_bl21":   PLIST/"bl21_de3_protein_extracted.tsv",   # ⚠ QZI 계열 — ID 로 조인하지 말 것
    # 구조 DB + 인덱스 (협업팀 Cameron Gilchrist 구축)
    "struct_bl21": Path("/mnt/af2results/mingyu/database/bacteriaDB/structures_UP000503272"),
    "struct_y19":  Path("/mnt/af2results/mingyu/database/bacteriaDB/structures_UP000034085"),
    "fd_idx_bl21": Path("/mnt/af2results/mingyu/database/bacteriaDB/folddisco/UP000503272"),
    "fd_idx_y19":  Path("/mnt/af2results/mingyu/database/bacteriaDB/folddisco/UP000034085"),
    "fd_idx_mg":   Path("/mnt/af2results/mingyu/database/folddisco/ecoli_folddisco/e_coli_folddisco"),
    # Folddisco 실행 바이너리 (Rust — 파이썬 환경 무관)
    "folddisco":   Path("/mnt/af2results/mingyu/folddisco/bin/folddisco"),
    # Track C 워크스페이스
    "fd_ws":       Path("/mnt/af2results/mingyu/workspace/seek_ni_insertase"),
    "cooc1_pdb":   Path("/mnt/af2results/mingyu/workspace/seek_ni_insertase/input/3kji.pdb"),
}

# BLASTp 4-category 분류 결과(원본은 Mac 로컬). 서버로 복사한 경로를 지정.
# 없으면 CELL 10-11 에서 mmseqs 로 재산출한다.
CATEGORY_TSV = {
    "BL21": DIR["external"] / "bl21_mg1655_category_detail.tsv",
    "Y19":  DIR["external"] / "y19_mg1655_category_detail.tsv",
}

# ---------------- 대상 프로테옴 ----------------
# ★ 조인 키는 upid 가 아니라 faa 안의 단백질 ID 다. upid 는 표시용 라벨일 뿐이어서,
#    여기를 구조 DB 의 UP 번호로 맞춰도 faa 가 다른 제출 기록에서 왔으면 소용이 없다.
#    실측: upid 를 UP000503272 로 적어두고 faa 는 협업팀 파일(QZI…)을 그대로 썼더니
#    구조 DB(QJZ…)와 교집합이 0 이었고, Part 8 조인이 전멸할 뻔했다.
#    -> faa 는 반드시 구조 DB 와 같은 기록에서 뽑은 것을 쓰고, CELL 02 의
#       "ID 공간 대조" 로 교집합을 눈으로 확인할 것.
PROTEOMES = {
    "BL21DE3": {"upid": "UP000503272", "taxid": 469008,
                "name": "Escherichia coli BL21(DE3)", "faa": ASSET["faa_bl21"]},
    "MG1655":  {"upid": "UP000000625", "taxid": 511145,
                "name": "Escherichia coli K-12 MG1655", "faa": ASSET["faa_mg1655"]},
    "Y19":     {"upid": "UP000034085", "taxid": 1261127,
                "name": "Citrobacter amalonaticus Y19", "faa": ASSET["faa_y19"]},
    "CHY":     {"upid": "UP000001320", "taxid": 246194,
                "name": "Carboxydothermus hydrogenoformans Z-2901", "faa": None},
}
# ---------------- 스크리닝 대상 스위치 ----------------
# PREY_STRAIN 만 바꾸면 CELL 08~12(서열 로드·분류·Track 배정)가 전부 따라간다.
#   "BL21DE3" = 이 프로젝트의 기본값. in vivo 에서 BL21 lysate 만 활성을 회복시켰으므로 정답에 가장 가깝다.
#   "Y19"     = Citrobacter amalonaticus Y19 로 같은 스크리닝을 돌릴 때만 변경.
PREY_STRAIN  = "BL21DE3"
REF_STRAIN   = "MG1655"        # 차등 분류의 기준 균주 (BL21/Y19 를 이것과 비교)
CATEGORY_KEY = {"BL21DE3": "BL21", "Y19": "Y19"}[PREY_STRAIN]   # CATEGORY_TSV 선택용

# 아래 둘은 고칠 일이 없다.
ID_SPACE = "GenBank"                            # 표시용 상수. 코드가 분기에 쓰지 않는다.
os.environ["MMSEQS_NUM_THREADS"] = str(THREADS)  # 위 THREADS 에서 자동 반영

# ---------------- 공용 헬퍼 ----------------
def sh(cmd, cwd=None, check=True, quiet=False):
    """포그라운드 셸 실행. 짧은 명령용."""
    r = subprocess.run(cmd, shell=True, cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if not quiet:
        print(r.stdout.rstrip())
    if check and r.returncode != 0:
        raise RuntimeError(f"명령 실패 (rc={r.returncode}): {cmd}")
    return r.stdout

def sh_bg(name, script, env=None):
    """긴 작업용 백그라운드 실행. scripts/<name>.sh 로 저장하고 nohup 실행."""
    sp = DIR["script"] / f"{name}.sh"
    lg = DIR["log"] / f"{name}.log"
    header = "#!/usr/bin/env bash\nset -euo pipefail\n"
    if env:
        header += f'source "$(conda info --base)/etc/profile.d/conda.sh"\nconda activate {env}\n'
    sp.write_text(header + textwrap.dedent(script))
    sp.chmod(0o755)
    subprocess.Popen(f'nohup bash "{sp}" > "{lg}" 2>&1 & echo $! > "{lg}.pid"',
                     shell=True, cwd=str(BASE))
    time.sleep(1)
    print(f"[백그라운드 시작] {name}\n  script: {sp}\n  log   : {lg}\n"
          f"  진행확인: CELL 15 에서 bg_tail('{name}')")
    return lg

def bg_tail(name, n=40):
    lg = DIR["log"] / f"{name}.log"
    pidf = Path(str(lg) + ".pid")
    alive = False
    if pidf.exists():
        pid = pidf.read_text().strip()
        alive = Path(f"/proc/{pid}").exists()
    print(f"=== {name} === {'실행 중' if alive else '종료됨'}")
    if lg.exists():
        print("".join(lg.read_text(errors="ignore").splitlines(keepends=True)[-n:]))
    else:
        print("(로그 없음)")
    return alive

def stage():
    """커널 재시작 후 어디까지 돼 있는지 확인하고 다음에 돌릴 셀을 알려준다."""
    g = globals()
    steps = [
        ("CELL 05-07  bait 서열",        "BAIT_ALL"),
        ("CELL 08     프로테옴 로드",     "hdr2seq"),
        ("CELL 09/11  4-category 분류",  "cls_prey"),
        ("CELL 12     Track A/B 배정",   "TRACK_A"),
        ("CELL 17     feasibility gate", "BAIT_DEPTH"),
        ("CELL 18     prey FASTA",       "sel"),
        ("CELL 20     best-hit 정리",     "PREY_MSA"),
        ("CELL 21     paired MSA",       "sdf"),
        ("CELL 23     RF2-PPI 집계",      "R_best"),
        ("CELL 28     Boltz 파싱",        "B"),
        ("CELL 30     ID crosswalk",     "XW"),
    ]
    nxt = None
    for label, var in steps:
        ok = var in g and g[var] is not None
        print(f"  [{'O' if ok else ' '}] {label}")
        if not ok and nxt is None:
            nxt = label
    print(f"\n  -> 다음에 돌릴 것: {nxt}" if nxt else "\n  -> 전부 완료")

# 커널을 재시작했거나 어디까지 돌았는지 헷갈릴 때 아무 셀에서나 stage() 를 부르면 된다.

def done(job, *artifacts):
    """작업이 끝났는가. 로그의 DONE_ 표시가 1순위,
    없으면 산출물 존재로 판정한다. 로그는 재실행 때 덮어써지므로
    산출물이 더 믿을 만한 증거다 (실측: 다운로드를 다시 눌렀다 끊어
    DONE_download 가 사라졌는데 파일 17,992개는 멀쩡했다)."""
    lg = DIR["log"] / f"{job}.log"
    if lg.exists() and "DONE_" in lg.read_text(errors="ignore")[-4000:]:
        return True
    return bool(artifacts) and all(Path(a).exists() for a in artifacts)

def have(*names):
    """앞 셀이 만든 변수가 살아 있는가 (커널 재시작 후 확인용)."""
    g = globals()
    return all(n in g and g[n] is not None for n in names)

def newer(a, b):
    """a 가 b 보다 나중에 만들어졌는가. 오래된 DB 로 검색한 결과를 걸러낸다."""
    a, b = Path(a), Path(b)
    return a.exists() and b.exists() and a.stat().st_mtime >= b.stat().st_mtime

def already(cond, what, redo=""):
    """이미 만들어 둔 산출물이 있으면 True 를 돌려준다.
    노트북을 새로 받아 CELL 01 부터 다시 돌릴 때 무거운 재계산을 건너뛰기 위한 것."""
    if cond:
        print(f"[건너뜀] {what} 이(가) 이미 있다.")
        if redo:
            print(f"  다시 만들려면: {redo}")
    return cond

def need(*conds):
    """셀 맨 앞에서 선행조건을 검사한다. 어긋나면 무엇을 먼저 해야 하는지 알리고 멈춘다."""
    bad = [m for ok, m in conds if not ok]
    if bad:
        raise RuntimeError("선행조건 미충족 — 아래를 먼저 해결할 것:\n  - " + "\n  - ".join(bad))
    print("[선행조건 OK]")

def conda_run(env, cmd, cwd=None, check=True):
    """지정 conda 환경에서 포그라운드 실행."""
    full = (f'source "$(conda info --base)/etc/profile.d/conda.sh" && '
            f'conda activate {env} && {cmd}')
    return sh(full, cwd=cwd, check=check)

print("BASE  =", BASE)
print("TOOLS =", TOOLS, "(RF2-PPI 설치 위치, folddisco 와 같은 층)")
print("REFDB =", RB, "(공용 세균 레퍼런스 DB)")
for grp in ["input", "result", "script", "tmp"]:
    print(f"  {grp}/")
    for k, v in DIR.items():
        if v.parent.name == grp or (grp == "tmp" and k == "tmp") or (grp == "script" and k == "script"):
            print(f"    {k:10s} {v}")
print("\nprey strain =", PREY_STRAIN, "| ID space =", ID_SPACE)
print("BL21 proteome ID =", PROTEOMES['BL21DE3']['upid'], "(구조 DB와 통일)")
```

---

## CELL 02 — 기존 자산 · 환경 점검

```python
# =============================================================================
# CELL 02 | 기존 자산 존재 확인 + GPU/디스크 상태
#   - 여기서 빨간 [X]가 뜨는 항목은 그 항목을 쓰는 Part에서 반드시 막힌다
# =============================================================================
# ---- 커널 확인 (제일 먼저) ----
# 노트북이 rf2ppi 환경이 아닌 다른 파이썬으로 돌면 뒤에서 ImportError 가 난다.
# 예전 환경의 ~/.local/lib/pythonX.Y/site-packages 를 물고 도는 경우가 흔하다.
import sys                     # CELL 01 을 건너뛰고 이 셀만 돌려도 되게 여기서도 import
print("### 커널 ###")
print("  python  :", sys.executable)
print("  version :", sys.version.split()[0])
_env_ok = CONDA_ENV_RF2 in sys.executable
print(f"  판정    : {'OK' if _env_ok else '⚠ rf2ppi 환경이 아니다'}")
if not _env_ok:
    print(f"    -> 커널을 'Python ({CONDA_ENV_RF2})' 로 다시 선택하거나 Antigravity 를 재시작할 것.")
    print(f"    -> 커널 등록 확인: cat ~/.local/share/jupyter/kernels/{CONDA_ENV_RF2}/kernel.json")

print("\n### 기존 자산 ###")
missing = []
for k, p in ASSET.items():
    ok = p.exists()
    if not ok:
        missing.append(k)
    size = ""
    if ok and p.is_file():
        size = f"{p.stat().st_size/1e6:.1f} MB"
    elif ok and p.is_dir():
        size = f"{len(list(p.iterdir()))} files"
    print(f"  [{'O' if ok else 'X'}] {k:14s} {p}  {size}")

# ---- ★ID 공간 대조 (이 검사가 없어서 한 번 크게 당했다) ----
# 프로테옴 faa 의 단백질 ID 와 구조 DB 파일이 같은 GenBank 제출 기록에서 왔는지 본다.
# "GenBank 로 통일" 만으로는 부족하다. 같은 균주라도 제출 기록이 다르면 accession
# 계열이 통째로 달라진다 (실측: 프로테옴 QZI… vs 구조 DB QJZ… -> 교집합 0).
# PROTEOMES 의 upid 는 라벨일 뿐 faa 를 고르지 않으므로, 실제 ID 를 맞대봐야 한다.
print("\n### ID 공간 대조 (프로테옴 faa vs 구조 DB) ###")
_xw = DIR["table"]/"id_crosswalk_struct_to_genbank.csv"
if not _xw.exists():
    print(f"  (crosswalk 아직 없음 — CELL 30 뒤에 다시 확인할 것)")
else:
    import pandas as pd
    _x = pd.read_csv(_xw).dropna(subset=["protein"])
    for _strain, _db in [("BL21DE3", "BL21"), ("Y19", "Y19")]:
        _f = PROTEOMES.get(_strain, {}).get("faa")
        if not _f or not Path(_f).exists():
            print(f"  {_strain:8s} faa 없음"); continue
        _ids = {l[1:].split()[0] for l in open(_f, errors="ignore") if l.startswith(">")}
        _st  = set(_x[_x.db == _db].protein)
        _ov  = len(_ids & _st)
        _pct = _ov / max(1, len(_ids))
        _mark = "OK" if _pct >= 0.95 else "⚠ ID 계열 불일치 — Part 8 조인이 전멸한다"
        print(f"  {_strain:8s} faa {len(_ids):5d} / 구조 {len(_st):5d} / 교집합 {_ov:5d} ({_pct:.1%})  {_mark}")
        if _pct < 0.95:
            print(f"    faa 예시   : {sorted(_ids)[:3]}")
            print(f"    구조 DB 예시: {sorted(_st)[:3]}")
            print("    -> 구조 DB 쪽 기록으로 faa 를 다시 준비할 것. upid 만 바꾸는 것으로는 해결되지 않는다.")

# folddisco 인덱스는 확장자 없는 본체 + .lookup/.offset/.type 4종 세트
print("\n### Folddisco 인덱스 세트 ###")
for k in ["fd_idx_bl21", "fd_idx_y19", "fd_idx_mg"]:
    base = ASSET[k]
    parts = {ext: (Path(str(base) + ext)).exists() for ext in ["", ".lookup", ".offset", ".type"]}
    print(f"  {k:12s} {base.name:16s} " + " ".join(f"{e or 'body'}={'O' if v else 'X'}"
                                                  for e, v in parts.items()))

print("\n### 공용 레퍼런스 DB (Part 2) ###")
for n in ["Bacteria/fasta", "Bacteria/bacteria_ref.fasta", "Bacteria/bactDB"]:
    q = REFDB / n
    extra = f"  ({len(list(q.iterdir()))} files)" if q.is_dir() else ""
    print(f"  [{'O' if q.exists() else 'X'}] {n:28s} {q}{extra}")

print("\n### 분류 TSV (없으면 CELL 10-11 fallback) ###")
for k, p in CATEGORY_TSV.items():
    print(f"  [{'O' if p.exists() else 'X'}] {k}: {p}")

print("\n### GPU ###")
sh("nvidia-smi --query-gpu=index,name,driver_version,memory.total,memory.used "
   "--format=csv,noheader || echo 'nvidia-smi 실패 — 드라이버 확인 필요'", check=False)

print("\n### 디스크 (워크스페이스 / 홈=conda 환경) ###")
sh(f'df -h "{BASE}" "$HOME"', check=False)

print("\n### 툴 확인 ###")
for t in ["mmseqs", "hhfilter", "aria2c", "conda", "curl"]:
    sh(f"command -v {t} >/dev/null && echo '  [O] {t}' || echo '  [X] {t} (미설치)'", check=False)
print("  ※ hhfilter 는 hhsuite 패키지에 들어 있다 (CELL 21 paired MSA 에서 필수).")
print("     없으면: conda install -y -c conda-forge -c bioconda hhsuite")

if missing:
    print("\n⚠ 없는 자산:", ", ".join(missing))
print("\n### 워크스페이스 구조 (CELL 01 에서 생성됨) ###")
sh(f"find {BASE} -maxdepth 2 -type d | sort | sed 's|{BASE}|  .|'", check=False)
print(f"  RF2-PPI 설치 위치: {RF2PPI_DIR}  [{'O' if RF2PPI_DIR.exists() else 'X — CELL 03 필요'}]")

print("\n### NVIDIA 드라이버 고정 (apt-mark hold) ###")
_hold = subprocess.run("apt-mark showhold 2>/dev/null | grep -c '^\\(lib\\)\\?nvidia'",
                       shell=True, capture_output=True, text=True).stdout.strip()
_n = int(_hold) if _hold.isdigit() else 0
if _n:
    print(f"  [O] nvidia 관련 {_n} 개 패키지가 hold 상태 — 커널 업데이트로 드라이버가 깨지지 않는다")
else:
    print("  [X] 미적용 — 커널 업데이트 시 NVML mismatch 재발 가능")
    print("      sudo apt-mark hold nvidia-* libnvidia-*")
```

---

## CELL 03 — Part 0. RF2-PPI 설치 스크립트 생성 (터미널 실행 권장)

```python
# =============================================================================
# CELL 03 | Part 0-1. RF2-PPI 설치 스크립트 생성
#   - conda create 는 노트북 커널보다 터미널이 안정적이라 "스크립트만 만들고" 출력한다
#   - 실행: 터미널에서 bash <출력된 경로>   (약 20~40분)
# =============================================================================
# False(기본) = 기존 rf2ppi 환경을 그대로 쓰고, RF2-PPI 가 직접 쓰는 경량 패키지와
#               코드/가중치만 받는다. epel-af2 는 이미 환경이 잡혀 있으므로 이쪽이다.
# True        = torch/numpy 를 RF2-PPI 공식 스펙(1.12.1+cu113 / 1.21.2)으로 덮어쓴다.
#               환경을 처음부터 만들 때만 쓸 것. 잘 붙어 있는 CUDA 조합이 되돌아간다.
INSTALL_PY_DEPS = False

env_block = f"""
source "$(conda info --base)/etc/profile.d/conda.sh"
# 같은 이름의 환경이 이미 있으면 create 가 실패한다(set -e 로 스크립트 중단).
# 처음부터 다시 만들려면 먼저:  conda env remove -n {CONDA_ENV_RF2}
if conda env list | awk '{{print $1}}' | grep -qx "{CONDA_ENV_RF2}"; then
  echo "[info] 기존 {CONDA_ENV_RF2} 환경을 그대로 쓴다"
else
  conda create -y -n {CONDA_ENV_RF2} python=3.9
fi
conda activate {CONDA_ENV_RF2}
"""

if INSTALL_PY_DEPS:
    deps_block = f"""
# ---- 파이썬 의존성 (RF2-PPI 공식 스펙: python 3.9 / torch 1.12.1+cu113) ----
conda install -y -c conda-forge -c bioconda hhsuite mmseqs2 aria2
pip install numpy==1.21.2 pandas==1.5.3 biopython==1.79 scipy==1.7.1 einops
pip install torch==1.12.1+cu113 -f https://download.pytorch.org/whl/torch_stable.html
# 드라이버가 최신(CUDA 13)이라 cu113 이 안 붙으면 위 줄 대신:
# pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cu118
# 노트북용 추가 패키지. ★numpy 를 다시 못박아야 한다 —
# matplotlib 3.7+ 는 numpy>=1.23 을 요구해서 위에서 깐 1.21.2 를 밀어내고,
# 그러면 scipy 1.7.1 과 torch 가 동시에 깨진다. 3.6.3 이 numpy 1.21 과 맞는 마지막 계열.
conda install -y -c conda-forge -c bioconda mmseqs2 aria2
pip install "numpy==1.21.2" "matplotlib==3.6.3" seaborn tqdm pyyaml jupyter ipykernel
python -m ipykernel install --user --name {CONDA_ENV_RF2} --display-name "Python ({CONDA_ENV_RF2})"
python -c "import numpy, scipy, torch, matplotlib; print('numpy', numpy.__version__, '| scipy', scipy.__version__, '| torch', torch.__version__, '| cuda', torch.cuda.is_available())"
"""
else:
    deps_block = """
# ---- INSTALL_PY_DEPS=False : 기존 환경의 torch/numpy 를 건드리지 않는다 ----
command -v hhfilter >/dev/null || conda install -y -c conda-forge -c bioconda hhsuite
# RF2-PPI 가 직접 import 하는 경량 의존성만 채운다.
# "numpy<2" 를 함께 적어야 pip 이 numpy 2.x 를 끌어올려 torch 연결을 깨뜨리지 않는다.
pip install "numpy<2" scipy einops biopython pandas
python -c "import torch, numpy; print('numpy', numpy.__version__, '| torch', torch.__version__, \
'| cuda', torch.cuda.is_available()); torch.zeros(3).numpy(); print('numpy<->torch OK')"
"""

code_block = f"""
# ---- RF2-PPI 코드 + 가중치 (워크스페이스가 아니라 툴 디렉터리에) ----
cd "{TOOLS}"
[ -d RoseTTAFold2-PPI ] || git clone https://github.com/CongLabCode/RoseTTAFold2-PPI.git
cd RoseTTAFold2-PPI/src/models
[ -f RF2-PPI.pt ] || wget --no-check-certificate https://conglab.swmed.edu/humanPPI/downloads/RF2-PPI.pt
ls -la RF2-PPI.pt

# ---- 예제 실행 (판정은 CELL 04 에서) ----
cd "{RF2PPI_DIR}/examples"
# 저장소 버그 보정: segment_pairs_input 이 존재하지 않는 segment_paired_msas_new/ 를
# 가리킨다. 실제 폴더는 segment_paired_msas/ 이고 expected_output 도 그 경로를 쓴다.
sed -i 's|segment_paired_msas_new/|segment_paired_msas/|' segment_pairs_input
python ../src/predict_list_PPI.py -list_fn segment_pairs_input \\
       -model_file ../src/models/RF2-PPI.pt
echo "=== 설치 스크립트 완료 ==="
"""

p = DIR["script"] / "install_rf2ppi.sh"
p.write_text("#!/usr/bin/env bash\nset -euo pipefail\n"
             + textwrap.dedent(env_block + deps_block + code_block))
p.chmod(0o755)
print(f"INSTALL_PY_DEPS = {INSTALL_PY_DEPS}")
print("\n※ 이 셀은 스크립트를 만들기만 한다. 설치는 아래를 터미널에서 실행할 때 일어난다:\n")
print(f"    bash {p}\n")
print("-" * 70)
print(p.read_text())
```

---

## CELL 04 — Part 0. 설치 검증 (★체크포인트 1)

```python
# =============================================================================
# CELL 04 | Part 0-2. 설치 검증 — 예제 결과가 expected_output 과 일치하는가
#   ★체크포인트 1: 판정 기준은 "값이 똑같은가"가 아니라 "같은 결론이 나오는가"다.
#     RF2-PPI 는 비결정적이고(README: 약 5% 쌍에서 SD>0.1, 중간 점수대일수록 심함),
#     expected_output 자체도 평균이 아니라 한 번의 추출값이다.
#     그래서 (1) 판정 구간(high/gray/low)이 일치하는지, (2) 극단값 쌍이 잘 맞는지로 본다.
# =============================================================================
ex = RF2PPI_DIR / "examples"
assert ex.exists(), f"RF2-PPI 미설치. CELL 03 스크립트를 터미널에서 먼저 실행하세요: {ex}"

log = ex / "segment_pairs_input.log"
if not log.exists():
    conda_run(CONDA_ENV_RF2,
              f'cd "{ex}" && python ../src/predict_list_PPI.py '
              f'-list_fn segment_pairs_input -model_file ../src/models/RF2-PPI.pt')

import pandas as pd
def read_rf2_log(path):
    """RF2-PPI 로그 읽기. 헤더 줄이 있을 수도 없을 수도 있어 숫자 행만 남긴다."""
    d = pd.read_csv(path, sep=r"\s+", names=["msa", "prob", "sec"])
    d = d[pd.to_numeric(d.prob, errors="coerce").notna()].copy()
    d["prob"] = d.prob.astype(float)
    d["key"] = d.msa.apply(lambda x: Path(x).name)      # 디렉터리 차이는 무시
    return d

got = read_rf2_log(log)
exp_files = list((ex / "expected_output").glob("segment*.log"))
print("=== 예제 실행 결과 ===")
print(got[["key", "prob"]].to_string(index=False))

if exp_files:
    expd = read_rf2_log(exp_files[0])
    m = got.merge(expd, on="key", suffixes=("_got", "_exp"))
    m["diff"] = (m.prob_got - m.prob_exp).abs()
    print("\n=== expected_output 대조 ===")
    print(m[["key", "prob_got", "prob_exp", "diff"]].to_string(index=False))
    BIN = lambda x: "high" if x >= 0.74 else ("gray" if x >= 0.30 else "low")
    m["bin_got"], m["bin_exp"] = m.prob_got.apply(BIN), m.prob_exp.apply(BIN)
    bins_ok = bool((m.bin_got == m.bin_exp).all())

    # 극단값(확실히 붙거나 확실히 안 붙는 쌍)은 재현성이 높아야 한다
    ext = m[(m.prob_exp >= 0.74) | (m.prob_exp < 0.15)]
    ext_ok = bool((ext["diff"] < 0.10).all()) if len(ext) else True

    print(f"\n  판정 구간 일치 : {int((m.bin_got == m.bin_exp).sum())}/{len(m)}")
    print(f"  극단값 최대 편차: {ext['diff'].max():.4f}" if len(ext) else "  극단값 쌍 없음")
    print(f"  전체 최대 편차  : {m['diff'].max():.4f}"
          f"  (중간 점수대의 큰 편차는 문서화된 비결정성이다)")
    if bins_ok and ext_ok:
        print("\n★체크포인트 1: 통과 — 설치 정상. CELL 05 로 진행한다.")
    else:
        print("\n★체크포인트 1: 실패 — CUDA/torch 버전을 점검한다.")
        if not bins_ok:
            print("   판정 구간이 어긋난 쌍:")
            print(m.loc[m.bin_got != m.bin_exp, ["key","prob_got","prob_exp"]].to_string(index=False))
        if not ext_ok:
            print("   극단값 쌍이 0.10 이상 벗어났다 — 이건 비결정성으로 설명되지 않는다.")
else:
    print("\nexpected_output 을 못 찾음 — 수동 대조 필요")
```

---

## CELL 05 — Part 1-1. Bait 서열 입력 (★필수 블로커)

```python
# =============================================================================
# CELL 05 | Part 1-1. Bait 서열 입력
#   ★이 셀이 비어 있으면 뒤 전부 안 돌아간다.
#   - 실제 발현·정제하는 그 서열을 쓴다 (A559W 변이체, His-tag 등 제거한 native)
#   - C. hydrogenoformans Z-2901 은 CODH 가 5개(I~V)라 gene name 만으로 II번 특정 불가
# =============================================================================
BAIT_SEQ = {}

# ---- ChCODH2 (C. hydrogenoformans CODH-II), 636 aa ----
# 559번이 A 인 것을 확인했다. A559W 변이의 번호 체계와 일치한다.
BAIT_WT = (
    "MAKQNLKSTDRAVQQMLDKAKREGIQTVWDRYEAMKPQCGFGETGLCCRHCLQGPCRINPFGDEPKVGIC"
    "GATAEVIVARGLDRSIAAGAAGHSGHAKHLAHTLKKAVQGKAASYMIKDRTKLHSIAKRLGIPTEGQKDE"
    "DIALEVAKAALADFHEKDTPVLWVTTVLPPSRVKVLSAHGLIPAGIDHEIAEIMHRTSMGCDADAQNLLL"
    "GGLRCSLADLAGCYMGTDLADILFGTPAPVVTESNLGVLKADAVNVAVHGHNPVLSDIIVSVSKEMENEA"
    "RAAGATGINVVGICCTGNEVLMRHGIPACTHSVSQEMAMITGALDAMILDYQCIQPSVATIAECTGTTVI"
    "TTMEMSKITGATHVNFAEEAAVENAKQILRLAIDTFKRRKGKPVEIPNIKTKVVAGFSTEAIINALSKLN"
    "ANDPLKPLIDNVVNGNIRGVCLFAGCNNVKVPQDQNFTTIARKLLKQNVLVVATGCGAGALMRHGFMDPA"
    "NVDELCGDGLKAVLTAIGEANGLGGPLPPVLHMGSCVDNSRAVALVAALANRLGVDLDRLPVVASAAEAM"
    "HEKAVAIGTWAVTIGLPTHIGVLPPITGSLPVTQILTSSVKDITGGYFIVELDPETAADKLLAAINERRA"
    "GLGLPW"
)
assert len(BAIT_WT) == 636 and BAIT_WT[558] == "A", "서열이 다릅니다"

# 기본은 WT. RF2-PPI 는 공진화(paired MSA)로 작동하므로 점 돌연변이 1개는
# 점수에 사실상 영향이 없고, Boltz-2 구조 예측에서도 마찬가지다.
BAIT_SEQ["ChCODH2_WT"] = BAIT_WT

# 실험 서열과 정확히 맞추려면 위 줄 대신 아래를 쓴다:
# BAIT_SEQ["ChCODH2_A559W"] = BAIT_WT[:558] + "W" + BAIT_WT[559:]

# 양성대조군: CooC1 은 CELL 06 이 3kji.pdb 에서 자동으로 뽑는다.
# CooC2 / CooT / CooJ 는 아래 UniProt 조회로 찾아 여기에 직접 넣는다.
# BAIT_SEQ["CooT"] = "..."

# --- 아직 없으면 UniProt 에서 후보 목록만 조회해 본다 (인터넷 필요) ---
import urllib.request, urllib.parse
import pandas as pd

def uniprot_tsv(query, fields="accession,id,gene_names,protein_name,length,sequence", size=100):
    url = ("https://rest.uniprot.org/uniprotkb/search?"
           + urllib.parse.urlencode({"query": query, "format": "tsv",
                                     "fields": fields, "size": size}))
    with urllib.request.urlopen(url, timeout=120) as r:
        txt = r.read().decode()
    rows = [l.split("\t") for l in txt.strip().split("\n")]
    return pd.DataFrame(rows[1:], columns=rows[0])

LOOKUP_CONTROLS = False   # True 면 CooC/CooT/CooJ 후보를 UniProt 에서 조회한다
print("bait 등록:", {k: len(v) for k, v in BAIT_SEQ.items()})

if LOOKUP_CONTROLS:
    try:
        df_coo = uniprot_tsv("organism_id:246194 AND "
                             "(gene:cooS OR gene:cooC OR gene:cooT OR gene:cooJ OR gene:cooF)")
        pd.set_option("display.max_colwidth", 60)
        print(df_coo[["Entry", "Gene Names", "Protein names", "Length"]].to_string(index=False))
        UNIPROT_COO = df_coo.set_index("Entry")["Sequence"].to_dict()
        print("\n>> 위 표에서 CODH-II(CooS-II) / CooC / CooT / CooJ 를 골라 BAIT_SEQ 에 넣으세요.")
        print(">> 서열은 UNIPROT_COO['<Entry>'] 로 꺼낼 수 있습니다.")
        print(">> ⚠ gene name 만으로 II번 특정 불가 — 실험에 쓰는 서열로 직접 대조 필수.")
    except Exception as e:
        print("UniProt 조회 실패:", e)
```

---

## CELL 06 — Part 1-1b. CooC1 서열을 3kji.pdb에서 추출

```python
# =============================================================================
# CELL 06 | Part 1-1b. 양성대조군 CooC1 서열 추출 (Folddisco 쿼리 구조 3kji.pdb 재사용)
#   - 양성대조군(CooC-CooS)에서 높은 점수가 안 나오면 E. coli 스크리닝 결과도 해석 불가
# =============================================================================
AA3to1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
          "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
          "THR":"T","TRP":"W","TYR":"Y","VAL":"V","MSE":"M","SEC":"U","PYL":"O"}

def pdb_chain_seqs(pdb_path):
    """PDB 파일에서 체인별 서열 추출. SEQRES 우선, 없으면 CA 원자 기준."""
    seqres, ca, seen = {}, {}, set()
    for line in Path(pdb_path).read_text(errors="ignore").splitlines():
        if line.startswith("SEQRES"):
            ch = line[11]
            seqres.setdefault(ch, []).extend(AA3to1.get(r, "X") for r in line[19:].split())
        elif line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            ch, resi = line[21], line[22:27]
            if (ch, resi) in seen:
                continue
            seen.add((ch, resi))
            ca.setdefault(ch, []).append(AA3to1.get(line[17:20].strip(), "X"))
    src = seqres if seqres else ca
    return {ch: "".join(v) for ch, v in src.items()}, ("SEQRES" if seqres else "CA")

if ASSET["cooc1_pdb"].exists():
    chains, src = pdb_chain_seqs(ASSET["cooc1_pdb"])
    print(f"3kji.pdb 체인 ({src} 기준):")
    for ch, s in chains.items():
        print(f"  chain {ch}  L={len(s)}  {s[:60]}...")
    # Folddisco 쿼리 residue 112/114 확인용 (인계문서 TODO: 체인 ID 확정)
    for ch, s in chains.items():
        if len(s) >= 114:
            print(f"  -> chain {ch}: res112={s[111]} res114={s[113]}  "
                  f"(metal-binding motif, Cys 계열이면 체인 확정)")
    # 가장 긴 체인을 CooC1 으로 채택 (필요하면 수동 지정)
    ch_best = max(chains, key=lambda c: len(chains[c]))
    BAIT_SEQ.setdefault("CooC1", chains[ch_best])
    print(f"\nBAIT_SEQ['CooC1'] <- chain {ch_best} (L={len(chains[ch_best])})")
else:
    print("3kji.pdb 없음:", ASSET["cooc1_pdb"], "\n-> CooC1 은 UniProt(CELL 05)에서 가져오세요.")
```

---

## CELL 07 — Part 1-1c. Bait 검증 + segment bait + baits.fasta

```python
# =============================================================================
# CELL 07 | Part 1-1c. Bait 검증 · segment bait 생성 · baits.fasta 저장
#   - 설계문서 §2: insertase는 apo 중간체에 붙을 가능성이 높다
#     -> full-length 뿐 아니라 C-cluster 주변 segment bait 를 병행하고 합집합을 취한다
# =============================================================================
assert BAIT_SEQ, "BAIT_SEQ 가 비었습니다. CELL 05 에서 서열을 채우세요."
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

for name, seq in list(BAIT_SEQ.items()):
    s = "".join(seq.split()).upper()
    BAIT_SEQ[name] = s
    bad = set(s) - VALID_AA
    print(f"{name:20s} L={len(s):5d}  비표준문자={bad if bad else '없음'}")

VARIANT_POS = 559        # A559W 변이 위치(1-based) = C-cluster 근방. 번호 체계가 다르면 여기만 고친다.

# primary bait 지정 (ChCODH2)
cands = [k for k in BAIT_SEQ if "CODH" in k.upper() or "A559W" in k.upper()]
assert cands, "ChCODH2 bait 를 찾을 수 없습니다. 이름에 'CODH' 또는 'A559W' 를 넣으세요."
BAIT_KEY = cands[0]
print("\nprimary bait =", BAIT_KEY)

# A559W 변이 확인
s = BAIT_SEQ[BAIT_KEY]
if len(s) >= VARIANT_POS:
    aa = s[VARIANT_POS - 1]
    tag = {"W": "A559W 변이체",
           "A": "야생형(WT) — 번호 체계 정상"}.get(aa, f"⚠ A 도 W 도 아님 — 번호 체계 재확인 필요")
    print(f"  {VARIANT_POS}번 잔기 = {aa}   ({tag})")
else:
    print(f"  ⚠ 길이 {len(s)} < {VARIANT_POS} — 서열/번호 체계 재확인 필요")

# ---- segment bait: C-cluster 주변(변이 위치 중심) + 전체 타일링 ----
SEG_WIN, SEG_STEP = 200, 100
SEGMENTS = {}
c = VARIANT_POS
lo, hi = max(0, c - 120), min(len(s), c + 80)          # C-cluster 주변
SEGMENTS[f"{BAIT_KEY}_seg{lo+1}-{hi}"] = s[lo:hi]
for st in range(0, max(1, len(s) - SEG_WIN + 1), SEG_STEP):
    SEGMENTS[f"{BAIT_KEY}_seg{st+1}-{min(st+SEG_WIN, len(s))}"] = s[st:st + SEG_WIN]

USE_SEGMENT_BAIT = True     # False 로 두면 full-length 만 사용
BAIT_ALL = dict(BAIT_SEQ)
if USE_SEGMENT_BAIT:
    BAIT_ALL.update(SEGMENTS)

(DIR["seq"] / "baits.fasta").write_text(
    "".join(f">{n}\n{v}\n" for n, v in BAIT_ALL.items()))
print(f"\nbait {len(BAIT_SEQ)}개 + segment {len(SEGMENTS)}개 = {len(BAIT_ALL)}개")
print("저장:", DIR["seq"] / "baits.fasta")
```

---

## CELL 08 — Part 1-2. 기존 프로테옴 서열 로드 (재다운로드 대신 재사용)

```python
# =============================================================================
# CELL 08 | Part 1-2. 프로테옴 서열 로드 — 기존 GenBank FASTA 재사용
#   - 조인 키는 GenBank ID 다. 다만 "GenBank 면 된다" 가 아니라 "구조 DB 와 같은
#     제출 기록의 GenBank" 여야 한다. 같은 균주라도 기록이 다르면 accession 계열이
#     통째로 달라진다 (QZI… vs QJZ…). CELL 02 의 ID 공간 대조로 반드시 확인할 것.
#   - 기대값: BL21 4,088 / Y19 5,325
# =============================================================================
import pandas as pd

def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path, errors="ignore"):
        if line.startswith(">"):
            if name: seqs[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.strip())
    if name: seqs[name] = "".join(buf)
    return seqs

for key, meta in PROTEOMES.items():
    f = meta.get("faa")
    if not f or not Path(f).exists():
        print(f"{key:9s} (서열 파일 없음 — 필요 시 UniProt 다운로드)")
        continue
    raw = read_fasta(f)
    meta["seqs"] = {h.split()[0]: sq for h, sq in raw.items()}   # ID = 헤더 첫 토큰
    meta["desc"] = {h.split()[0]: " ".join(h.split()[1:]) for h in raw}
    print(f"{key:9s} {len(meta['seqs']):5d} proteins   {f.name}")

hdr2seq = PROTEOMES[PREY_STRAIN]["seqs"]          # prey 서열 사전 (GenBank ID -> seq)
hdr2desc = PROTEOMES[PREY_STRAIN]["desc"]
ex_id = next(iter(hdr2seq))
print(f"\nprey = {PREY_STRAIN}, {len(hdr2seq)}개, ID 예시: {ex_id}  ({hdr2desc[ex_id][:60]})")

# 개수를 고정값과 맞대는 건 의미가 없다 — 어느 제출 기록에서 뽑았느냐에 따라 달라진다.
# 정말 확인해야 할 것은 "이 faa 의 ID 가 구조 DB 에 있느냐" 다.
# (옛 파일 기준 개수는 BL21 4,088 / Y19 5,325 였고, 그건 QZI 계열이라 구조 DB 와 안 붙었다.)
_xw = DIR["table"]/"id_crosswalk_struct_to_genbank.csv"
print()
for k, db in [("BL21DE3", "BL21"), ("Y19", "Y19")]:
    if "seqs" not in PROTEOMES.get(k, {}): continue
    ids = set(PROTEOMES[k]["seqs"])
    line = f"  {k:8s} {len(ids):5d}개  {Path(PROTEOMES[k]['faa']).name}"
    if _xw.exists():
        st = set(pd.read_csv(_xw).query("db == @db").protein.dropna())
        ov = len(ids & st)
        line += f"  | 구조 DB 교집합 {ov}/{len(ids)} ({ov/max(1,len(ids)):.1%})"
        line += "  OK" if ov/max(1,len(ids)) >= 0.95 else "  ⚠ ID 계열 불일치"
    else:
        line += "  | (crosswalk 없음 — CELL 30 뒤 CELL 02 로 확인)"
    print(line)
```

---

## CELL 09 — Part 1-3a. 기존 4-category 분류 로드 + 값 검증

```python
# =============================================================================
# CELL 09 | Part 1-3a. 기존 BLASTp 4-category 분류 로드 (재계산 스킵)
#   - 기대값(BL21->MG1655): 2171 / 1575 / 125 / 217 (합 4088)
#   - 기대값(Y19->MG1655) :   39 / 1856 / 2067 / 1363 (합 5325)
#   - PREY_STRAIN 에 맞는 TSV 를 자동 선택한다 (CELL 01 의 CATEGORY_KEY)
#   - 파일이 없으면 cls_prey = None 이 되고, CELL 10-11 fallback 으로 넘어간다
# =============================================================================
import pandas as pd, numpy as np

CANON = {"identical": "identical",
         "high similarity": "high-similarity", "high-similarity": "high-similarity",
         "high_similarity": "high-similarity", "high": "high-similarity",
         "low similarity": "low-similarity", "low-similarity": "low-similarity",
         "low_similarity": "low-similarity", "low": "low-similarity",
         "strain specific": "strain-specific", "strain-specific": "strain-specific",
         "strain_specific": "strain-specific", "specific": "strain-specific",
         "no hit": "strain-specific", "no_hit": "strain-specific"}

def load_category_tsv(path):
    """컬럼명이 파일마다 달라 자동 탐지."""
    df = pd.read_csv(path, sep=None, engine="python")
    cols = {c.lower(): c for c in df.columns}
    idc = next((cols[c] for c in cols if c in
                ("protein", "protein_id", "query", "qseqid", "id", "accession")), df.columns[0])
    catc = next((cols[c] for c in cols if "categ" in c or c == "class"), None)
    assert catc, f"category 컬럼을 못 찾음: {list(df.columns)}"
    out = pd.DataFrame({"protein": df[idc].astype(str),
                        "category": df[catc].astype(str).str.strip().str.lower().map(
                            lambda x: CANON.get(x, x))})
    for extra, name in [("pident", "pident"), ("best_hit", "best_hit"), ("cov", "cov")]:
        c = next((cols[k] for k in cols if k == extra or k.startswith(extra)), None)
        if c: out[name] = df[c]
    return out

# ⚠ 아래 값은 옛 faa(QZI 계열, BL21 4,088개)로 산출된 것이다. 구조 DB 와 같은
#   기록의 faa 로 바꾸면 단백질 세트가 달라져 숫자도 달라진다. 정확히 맞기를 기대하지
#   말고, 자릿수와 비율이 비슷한지만 본다 (identical 이 절반쯤, strain-specific 이 수백).
EXPECT_CAT = {   # 인계문서 §1-2 의 확정값 (옛 ID 공간 기준)
    "BL21DE3": {"identical": 2171, "high-similarity": 1575,
                "low-similarity": 125, "strain-specific": 217},
    "Y19":     {"identical": 39,   "high-similarity": 1856,
                "low-similarity": 2067, "strain-specific": 1363},
}

cls_prey = None
p = CATEGORY_TSV[CATEGORY_KEY]
if p.exists():
    cls_prey = load_category_tsv(p)
    vc = cls_prey.category.value_counts()
    print(f"=== {PREY_STRAIN} -> {REF_STRAIN} (로드) ===")
    print(vc.to_string())
    print("총:", len(cls_prey))
    expected = EXPECT_CAT[PREY_STRAIN]
    print("\n기대값 대조:")
    for k, v in expected.items():
        g = int(vc.get(k, 0))
        print(f"  {k:16s} {g:5d} / {v:5d}  {'OK' if g == v else '⚠ 불일치'}")
else:
    print(f"[없음] {p}")
    print("-> Mac 로컬 /Users/mingyukim/Desktop/EPEL/bacteriaDB_construction/ 에서 서버로 복사하거나,")
    print("   CELL 10-11 로 mmseqs 재산출 (값이 2171/1575/125/217 과 맞는지 대조).")
```

---

## CELL 10 — Part 1-3b. (fallback) mmseqs easy-search 재실행

```python
# =============================================================================
# CELL 10 | Part 1-3b. [fallback] 분류 TSV 가 없을 때만 — 양방향 mmseqs easy-search
#   - 역방향(MG1655->BL21)도 산출한다: MG1655 에만 있는 단백질 = 억제자 후보(rcnAB 등)
#   - CELL 09 에서 cls_prey 이 로드됐으면 이 셀은 건너뛴다
# =============================================================================
_FWD = f"{PREY_STRAIN}_vs_{REF_STRAIN}"
_REV = f"{REF_STRAIN}_vs_{PREY_STRAIN}"
if cls_prey is not None:
    print("CELL 09 에서 이미 로드됨 — 이 셀은 건너뜁니다.")
elif already((DIR["search"]/f"{_FWD}.m8").exists() and (DIR["search"]/f"{_REV}.m8").exists(),
             "mmseqs 검색 결과 (.m8)", f"rm {DIR['search']}/*_vs_*.m8"):
    print("  -> CELL 11 로 바로 가면 된다. 30~60분을 아낀다.")
else:
    fa_prey = PROTEOMES[PREY_STRAIN]["faa"]
    fa_ref  = PROTEOMES[REF_STRAIN]["faa"]
    FWD, REV = f"{PREY_STRAIN}_vs_{REF_STRAIN}", f"{REF_STRAIN}_vs_{PREY_STRAIN}"
    FMT_EASY = "query,target,fident,alnlen,evalue,bits,qlen,tlen,qcov,tcov"  # 10컬럼 (CELL 11 의 COLS 와 짝)
    S, T = DIR["search"], DIR["tmp"]
    script = f"""
    mmseqs easy-search "{fa_prey}" "{fa_ref}" "{S}/{FWD}.m8" "{T}/fwd" \\
      --format-output "{FMT_EASY}" -s 7.5 -e 1e-5 --max-seqs 5 --threads {THREADS}
    mmseqs easy-search "{fa_ref}" "{fa_prey}" "{S}/{REV}.m8" "{T}/rev" \\
      --format-output "{FMT_EASY}" -s 7.5 -e 1e-5 --max-seqs 5 --threads {THREADS}
    wc -l "{S}/{FWD}.m8" "{S}/{REV}.m8"
    echo DONE_easy_search
    """
    sh_bg("easy_search", script)
    print("\n30분~1시간 소요. CELL 15 로 진행 확인 후 CELL 11 실행.")
```

---

## CELL 11 — Part 1-3c. (fallback) m8 → 4-category 분류

```python
# =============================================================================
# CELL 11 | Part 1-3c. [fallback] m8 -> 4-category 분류 + 기대값 대조
#   - identical: pident>=99.9 & cov>=0.99 / high: >=90 / low: <90 / no hit: strain-specific
# =============================================================================
COLS = ["query","target","fident","alnlen","evalue","bits","qlen","tlen","qcov","tcov"]

def classify(m8_path, all_ids):
    df = pd.read_csv(m8_path, sep="\t", names=COLS)
    best = df.sort_values("bits", ascending=False).drop_duplicates("query").set_index("query")
    rows = []
    for pid in all_ids:
        if pid not in best.index:
            rows.append((pid, None, 0.0, 0.0, "strain-specific")); continue
        r = best.loc[pid]
        pident, cov = r.fident * 100, min(r.qcov, r.tcov)
        cat = ("identical" if (pident >= 99.9 and cov >= 0.99)
               else "high-similarity" if pident >= 90 else "low-similarity")
        rows.append((pid, r.target, pident, cov, cat))
    return pd.DataFrame(rows, columns=["protein","best_hit","pident","cov","category"])

if cls_prey is not None:
    print("이미 분류가 로드되어 있어 건너뜁니다.")
else:
    FWD, REV = f"{PREY_STRAIN}_vs_{REF_STRAIN}", f"{REF_STRAIN}_vs_{PREY_STRAIN}"
    prey_ids = list(PROTEOMES[PREY_STRAIN]["seqs"])
    ref_ids  = list(PROTEOMES[REF_STRAIN]["seqs"])
    cls_prey = classify(DIR["search"]/f"{FWD}.m8", prey_ids)
    cls_rev  = classify(DIR["search"]/f"{REV}.m8", ref_ids)
    print(f"=== {PREY_STRAIN} 기준 ===")
    print(cls_prey.category.value_counts().to_string())
    print(f"\n=== {REF_STRAIN} 기준 ({PREY_STRAIN} 결손 = 억제자 후보) ===")
    print(cls_rev.category.value_counts().to_string())
    cls_rev.to_csv(DIR["table"]/f"classify_{REV}.csv", index=False)
    print(f"\n기대값({PREY_STRAIN}):", EXPECT_CAT[PREY_STRAIN])
    print("역방향 MG1655->BL21 기대값: 2199 / 1569 / 152 / 380")
```

---

## CELL 12 — Part 1-3d. Track A / Track B 확정

```python
# =============================================================================
# CELL 12 | Part 1-3d. Track 배정 확정 + 저장
#   Track A (RF2-PPI, 공진화 가능)  = identical + high-similarity  (~3,746)
#   Track B (Boltz-2, 공진화 불가)  = low-similarity + strain-specific (~342)
#   ※ Track A 는 균주 특이성을 못 잡는다. 균주 특이성은 이 분류 자체가 잡는 것이다.
# =============================================================================
assert cls_prey is not None, "분류 결과가 없습니다. CELL 09 또는 CELL 10-11 을 먼저 실행하세요."
cls_prey.to_csv(DIR["table"]/f"classify_{PREY_STRAIN}_vs_{REF_STRAIN}.csv", index=False)

known = set(hdr2seq)
cls_prey["in_proteome"] = cls_prey.protein.isin(known)
if (~cls_prey.in_proteome).any():
    n = int((~cls_prey.in_proteome).sum())
    print(f"⚠ 분류에는 있으나 프로테옴 FASTA 에 없는 ID {n}개 — ID 체계(GenBank 버전 .1 등) 확인")
    print("  예시:", cls_prey.loc[~cls_prey.in_proteome, "protein"].head(5).tolist())

TRACK_A = cls_prey.query("category in ['identical','high-similarity'] and in_proteome").protein.tolist()
TRACK_B = cls_prey.query("category in ['low-similarity','strain-specific'] and in_proteome").protein.tolist()

pd.Series(TRACK_A).to_csv(DIR["table"]/"track_A_ids.txt", index=False, header=False)
pd.Series(TRACK_B).to_csv(DIR["table"]/"track_B_ids.txt", index=False, header=False)
# 기대 규모는 옛 ID 공간 기준 Track A ~3,746 / Track B ~342 였다.
# faa 를 바꿨으면 단백질 세트가 달라져 숫자도 달라진다 — 자릿수만 본다.
print(f"Track A (RF2-PPI): {len(TRACK_A):5d}")
print(f"Track B (Boltz-2): {len(TRACK_B):5d}")
print("저장:", DIR["table"]/"track_A_ids.txt", "/", DIR["table"]/"track_B_ids.txt")
```

---

## CELL 13 — Part 2. Bacteria reference proteome 다운로드 (백그라운드)

```python
# =============================================================================
# CELL 13 | Part 2-1. 세균 레퍼런스 프로테옴 다운로드 (세균만 선별)
#   ⚠ 전체 tarball(Reference_Proteomes_*.tar.gz)은 328GB 다. 고세균·진핵·바이러스와
#     DNA 서열까지 들어 있어 우리한테는 대부분 낭비다 (실측 ETA 9시간).
#   - README 표에서 superregnum=bacteria 인 프로테옴만 골라 단백질 fasta 만 받는다.
#     17,992개, 약 20GB, aria2 병렬로 30~60분.
# =============================================================================
REL_URL = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/"
           "knowledgebase/reference_proteomes")
# README 표 형식: UP_ID  TaxID  OSCODE  superregnum  #(1) #(2) #(3)  species
script = f"""
cd "{RB}"
REL={REL_URL}
curl -s "$REL/README" | awk -v rel="$REL" '$4=="bacteria" {{printf "%s/Bacteria/%s/%s_%s.fasta.gz\\n", rel, $1, $1, $2}}' > urls.txt
echo "URL 개수: $(wc -l < urls.txt)"
aria2c -i urls.txt -d fasta -j 16 -x 4 -c --auto-file-renaming=false --console-log-level=warn --summary-interval=120
echo "받은 파일: $(ls fasta | wc -l)"
du -sh fasta
echo DONE_download
"""
_n_fa = len(list((RB/"fasta").glob("*.fasta.gz")))
if not already(_n_fa >= 17000,
               f"세균 프로테옴 {_n_fa}개", f"rm -rf {RB}/fasta"):
    sh_bg("part2_download", script)
print("\n디스크 여유:")
sh(f'df -h "{REFDB}"', check=False)
```

---

## CELL 14 — Part 2b. 병합 + mmseqs createdb (백그라운드)

```python
# =============================================================================
# CELL 14 | Part 2-2. taxid 를 헤더에 박아 하나로 병합 -> mmseqs DB 생성
#   - 헤더를 ">taxid|accession" 으로 만들어 두면 pairing 시 organism key 가 공짜
#   - CELL 13 의 DONE_download 확인 후 실행. 병합 10~20분 + createdb/index 1~2시간.
# =============================================================================
need((len(list((RB/"fasta").glob("*.fasta.gz"))) >= 17000,
      f"fasta 파일이 부족하다 ({len(list((RB/'fasta').glob('*.fasta.gz')))}개) — CELL 13 을 끝낼 것"))
script = f"""
cd "{RB}"
# 이전 실행이 남긴 DB/인덱스를 먼저 치운다.
# 반쪽 DB 위에 덮어쓰면 조각이 섞여 검색이 조용히 틀어진다 (실측: 28M vs 60M 서열).
rm -f bactDB* bacteria_ref.fasta
OUT=bacteria_ref.fasta
: > "$OUT"
# 진행률은 stderr 로 보낸다 (stdout 은 $OUT 으로 리다이렉트되므로).
# 이게 없으면 병합이 끝날 때까지 로그가 비어 있어 멈춘 것처럼 보인다.
TOT=$(ls fasta/*.fasta.gz | wc -l)
N=0
for f in fasta/*.fasta.gz; do
  tax=$(basename "$f" .fasta.gz | cut -d_ -f2)
  zcat "$f" | awk -v t="$tax" '/^>/{{split($1,a,"|"); acc=(length(a)>=2?a[2]:substr($1,2)); print ">"t"|"acc; next}}{{print}}'
  N=$((N+1))
  [ $((N % 1000)) -eq 0 ] && echo "  병합 $N/$TOT" >&2
done >> "$OUT"
echo "병합 완료: $TOT 개" >&2
echo "sequences: $(grep -c '^>' "$OUT")"
echo "proteomes: $(grep '^>' "$OUT" | cut -d'|' -f1 | tr -d '>' | sort -u | wc -l)"
mmseqs createdb "$OUT" bactDB
mmseqs createindex bactDB "{DIR['tmp']}/idx" --threads {THREADS} --search-type 1 || true
ls -la
echo DONE_bactdb
"""
if not already(done("part2_bactdb", RB/"bactDB.index"),
               "bactDB", f"rm -f {RB}/bactDB* {RB}/bacteria_ref.fasta"):
    sh_bg("part2_bactdb", script, env=CONDA_ENV_RF2)
```

---

## CELL 15 — 백그라운드 작업 모니터 (수시 실행)

```python
# =============================================================================
# CELL 15 | 백그라운드 작업 상태 확인 — 긴 단계 돌릴 때마다 이 셀로 확인
# =============================================================================
# 이름을 고정 목록으로 들면 나중에 추가한 작업(part7_boltz_focus 등)이 안 보인다.
# 로그 폴더를 훑되, 파이프라인 순서를 아는 것부터 먼저 보여준다.
ORDER = ["easy_search", "part2_download", "part2_bactdb", "part3_gate",
         "part4_prey", "part6_rf2ppi", "part7_boltz"]
found = sorted(f.stem for f in DIR["log"].glob("*.log"))
jobs  = [j for j in ORDER if j in found] + [j for j in found if j not in ORDER]
for job in jobs:
    bg_tail(job, n=8)
    print("-" * 70)

print("\n### 디스크 ###"); sh(f"df -h {BASE}", check=False)
print("### GPU ###");   sh("nvidia-smi --query-gpu=index,utilization.gpu,memory.used,memory.total "
                           "--format=csv,noheader", check=False)
```

---

## CELL 15b — 장기 작업 진행률 (숫자로 확인)

```python
# =============================================================================
# CELL 15b | boltz / RF2-PPI 진행률을 숫자로 본다
#   CELL 15 는 로그 꼬리만 보여준다. tqdm 은 \r 로 같은 줄을 덮어쓰므로 tail 에는
#   진행률이 안 잡힌다. 여기서는 \r 을 줄바꿈으로 바꾼 뒤 마지막 상태를 파싱한다.
#   GPU 0 은 디스플레이가 물려 있어 RF2-PPI 에서 13배 느렸다. Boltz 도 같은 함정에
#   빠질 수 있으니 s/pair 를 반드시 확인하고, 예상의 3배를 넘으면 죽이고 GPU 1 로
#   옮기는 편이 낫다.
# =============================================================================
def _tqdm_state(log):
    """tqdm 마지막 상태 -> (done, total, sec_per_item, 경과, 남은)"""
    txt = log.read_text(errors="ignore").replace("\r", "\n")
    hits = re.findall(r"(\d+)/(\d+)\s*\[([\d:]+)<([\d:?]+),\s*([\d.]+)(s/it|it/s)", txt)
    if not hits:
        return None
    d, t, el, eta, rate, unit = hits[-1]
    sec = float(rate) if unit == "s/it" else 1.0 / max(float(rate), 1e-9)
    return int(d), int(t), sec, el, eta

_alive = subprocess.run(["pgrep", "-af", "[b]oltz predict"],
                        capture_output=True, text=True).stdout.strip().splitlines()
print("=" * 74)
print(f"살아있는 boltz 프로세스: {len(_alive)}개")
for _l in _alive:
    print("   ", _l[:110])
print("-" * 74)

EXPECT_SEC = 60          # GPU 1 · max_msa_seqs 2048 · 1000aa 기준 실측치
for _job in ["part7_boltz", "part7_boltz_focus", "part7_boltz_cmp"]:
    _log = DIR["log"] / f"{_job}.log"
    if not _log.exists():
        continue
    _txt = _log.read_text(errors="ignore")
    _st, _oom = _tqdm_state(_log), _txt.count("ran out of memory")
    if _st is None:
        print(f"[{_job}] tqdm 줄 없음 — 전처리 중이거나 이미 끝남")
    else:
        _d, _t, _sec, _el, _eta = _st
        _flag = ("" if _sec <= EXPECT_SEC * 3 else
                 f"   ⚠ {EXPECT_SEC}초 예상인데 {_sec:.0f}초 — GPU 배치를 의심할 것")
        print(f"[{_job}] {_d}/{_t}  {_sec:.1f}s/pair  경과 {_el}  남은 {_eta}{_flag}")
    if _oom:
        print(f"   ⚠ OOM skip {_oom}건 — 그만큼 결과가 조용히 비어 있다 (CELL 28 로 확인)")
    if _txt.strip() and _txt.strip().splitlines()[-1].startswith("DONE_"):
        print("   완료 마커 확인 — 다음 셀로 넘어가도 된다")
print("-" * 74)

# --- RF2-PPI: replicate 별 완료 줄 수 + 실측 속도로 ETA ---
_NREP = globals().get("N_REPLICATE", 3)
_inp  = DIR["rf2ppi"] / "input_file"
_tot  = sum(1 for _ in open(_inp)) if _inp.exists() else 0
_sum, _secs = 0, []
for _rep in range(1, _NREP + 1):
    _rl = DIR["rf2ppi"] / f"input_rep{_rep}.log"
    if not _rl.exists():
        print(f"[part6_rf2ppi] rep{_rep}: 대기")
        continue
    _n = 0
    for _line in open(_rl, errors="ignore"):
        _f = _line.split()
        if len(_f) >= 3:
            _n += 1
            try:    _secs.append(float(_f[2]))
            except ValueError: pass
    _sum += _n
    print(f"[part6_rf2ppi] rep{_rep}: {_n}/{_tot}")
if _tot and _secs:
    _left = _tot * _NREP - _sum
    _rate = sum(_secs[-200:]) / len(_secs[-200:])
    print(f"  전체 {_sum}/{_tot*_NREP} ({100*_sum/(_tot*_NREP):.1f}%)  "
          f"{_rate:.2f}s/pair  남은 시간 약 {_left*_rate/3600:.1f}시간")
print("-" * 74)
sh("nvidia-smi --query-gpu=index,utilization.gpu,memory.used,memory.total "
   "--format=csv,noheader", check=False)
sh("nvidia-smi --query-compute-apps=pid,used_memory --format=csv,noheader", check=False)
```

---

## CELL 16 — Part 3. Feasibility gate 검색

```python
# =============================================================================
# CELL 16 | Part 3-1. Feasibility gate — bait orthologue 보유 genome 수 세기
#   - paired MSA 깊이 상한 = min(bait orthologue genome 수, prey orthologue genome 수)
#   - CooS 는 CO-oxidizing anaerobe 에만 분포 -> 이 값이 낮을 수 있다. 절대 건너뛰지 말 것.
# =============================================================================
need((done("part2_bactdb", RB/"bactDB.index"),
      "CELL 14 미완료 — 미완성 DB 로 검색하면 결과가 조용히 틀어진다"),
     ((RB/"bactDB").exists(), "bactDB 가 없다 — CELL 14 를 먼저 돌릴 것"),
     ((DIR["seq"]/"baits.fasta").exists(), "baits.fasta 가 없다 — CELL 07 을 먼저 돌릴 것"))
# 13컬럼. CELL 17/20 의 AC 컬럼명과 반드시 같아야 한다 (어긋나면 조용히 잘못 파싱된다)
FMT_ALN = "query,target,fident,evalue,bits,qstart,qend,qlen,tstart,tend,tlen,qaln,taln"
script = f"""
cd "{DIR['db']}"
# 재실행 시 mmseqs 가 "exists already!" 로 거부하므로 중간 산출물도 함께 치운다
rm -f baitDB* "{DIR['search']}"/bait_res*
mmseqs createdb "{DIR['seq']}/baits.fasta" baitDB
mmseqs search baitDB "{RB}/bactDB" "{DIR['search']}/bait_res" "{DIR['tmp']}/bs" \\
  -s 7.5 --num-iterations 3 -e 1e-3 --max-seqs 20000 --threads {THREADS}
mmseqs convertalis baitDB "{RB}/bactDB" "{DIR['search']}/bait_res" "{DIR['search']}/bait_hits.m8" \\
  --format-output "{FMT_ALN}" --threads {THREADS}
wc -l "{DIR['search']}/bait_hits.m8"
echo DONE_gate
"""
if not already(newer(DIR["search"]/"bait_hits.m8", RB/"bactDB"),
               "bait_hits.m8 (현재 bactDB 기준)", f"rm {DIR['search']}/bait_hits.m8"):
    sh_bg("part3_gate", script, env=CONDA_ENV_RF2)
    print("\n★ 이 셀은 시작만 한다. 1~2분 뒤 CELL 15 에서 DONE_gate 를 확인하고 CELL 17 로 갈 것.")
```

---

## CELL 17 — Part 3b. 깊이 판정 (★체크포인트 2)

```python
# =============================================================================
# CELL 17 | Part 3-2. Gate 판정 — Track A 를 진행할지 여기서 결정한다
#   >=500 : Track A 진행 / 200-500 : 진행하되 replicate 5회, 경계 점수 불신
#   <200  : Track A 포기, Track B(CELL 25~)에 자원 집중
# =============================================================================
need(((DIR["search"]/"bait_hits.m8").exists(), "CELL 16 을 먼저 돌릴 것"),
     # 로그의 DONE_ 표시만 보면, 로그를 치운 뒤에는 멀쩡한 결과를 두고도 막힌다.
     # 산출물이 더 믿을 만한 증거다 (CELL 14 에서 같은 이유로 한 번 막혔다).
     (done("part3_gate", DIR["search"]/"bait_hits.m8"),
      "CELL 16 검색이 아직 돌고 있다 — CELL 15 에서 DONE_gate 를 확인한 뒤 다시 실행할 것 (1~2분)"),
     (newer(DIR["search"]/"bait_hits.m8", RB/"bactDB"),
      "bait_hits.m8 이 bactDB 보다 오래됐다 — 옛 DB 로 검색한 결과다. CELL 16 을 다시 돌릴 것"))
AC = ["query","target","fident","evalue","bits","qstart","qend","qlen",
      "tstart","tend","tlen","qaln","taln"]
bh = pd.read_csv(DIR["search"]/"bait_hits.m8", sep="\t", names=AC)
bh["taxid"] = bh.target.str.split("|").str[0]

gate = bh.groupby("query").taxid.nunique().sort_values(ascending=False)
print("Bait 별 orthologue 보유 genome 수 (= paired MSA 깊이 상한)\n")
for q, n in gate.items():
    v = "OK" if n >= 500 else ("주의(replicate 5회)" if n >= 200 else "Track A 부적합")
    print(f"  {q:34s} {n:6d}  -> {v}")

BAIT_DEPTH = gate.to_dict()
d = BAIT_DEPTH.get(BAIT_KEY, 0)
TRACK_A_GO = d >= 200
N_REPLICATE = 3 if d >= 500 else 5
print(f"\n★체크포인트 2: {BAIT_KEY} depth={d}")
print(f"  Track A 진행 = {TRACK_A_GO}, replicate = {N_REPLICATE}회")
if not TRACK_A_GO:
    print("  -> CELL 18~24 를 건너뛰고 CELL 25(Track B, Boltz-2)로 바로 간다.")
    print("     이건 실패가 아니라 이 문제의 생물학적 성질이다 "
          "(E. coli host factor 가 CODH 와 공진화했을 이유가 없다).")
```

---

## CELL 18 — Part 4. Track A prey FASTA

```python
# =============================================================================
# CELL 18 | Part 4-1. Track A prey FASTA (길이 필터)
#   - bait(~630) + prey 합이 너무 길면 GPU 메모리 초과. 긴 건 제외 후 별도 segment 처리.
# =============================================================================
need((have("TRACK_A", "hdr2seq"), "CELL 08 과 CELL 12 를 먼저 돌릴 것 (stage() 로 확인)"))
MAX_PREY_LEN = 1000

# ---- 파일럿 모드 ----
# 전체 3,750개를 돌리기 전에 알려진 후보 몇 개로 파이프라인을 먼저 검증할 때 쓴다.
# 비워두면 Track A 전체. 유전자/제품명 키워드를 넣으면 그것만 추린다.
PILOT_GENES = []   # 예: ["hypB","slyD","mrp","yeiR","yjiA","hypA","nikA","iscU","erpA","nfuA"]

if PILOT_GENES:
    import re
    cat_of = cls_prey.set_index("protein")["category"].to_dict()
    hit, miss = {}, []
    for g in PILOT_GENES:
        pat = re.compile(rf"(?<![A-Za-z]){re.escape(g)}(?![A-Za-z])", re.I)
        found = [pid for pid in (TRACK_A + TRACK_B) if pat.search(hdr2desc.get(pid, ""))]
        (hit.setdefault(g, found) if found else miss.append(g))
    for g, v in hit.items():
        for pid in v[:3]:
            print(f"  {g:8s} {pid:16s} [{cat_of.get(pid,'?'):15s}] {hdr2desc[pid][:55]}")
    if miss:
        print(f"\n  ⚠ 못 찾은 유전자: {miss}")
        print("     GenBank .faa 헤더는 gene symbol 이 아니라 product 설명이다.")
        print("     'nickel', 'metallochaperone', 'ATP-binding', 'Fe-S' 같은 product 키워드로 다시 찾을 것:")
        print("     [d for d in hdr2desc.values() if 'nickel' in d.lower()][:10]")
    TRACK_A_FULL = TRACK_A
    TRACK_A = sorted({pid for v in hit.values() for pid in v})
    print(f"\n파일럿 모드: Track A {len(TRACK_A_FULL)} -> {len(TRACK_A)}개")

sel, skipped = {}, []
for pid in TRACK_A:
    s = hdr2seq[pid]
    (skipped.append((pid, len(s))) if len(s) > MAX_PREY_LEN else sel.update({pid: s}))

(DIR["seq"]/"prey_trackA.fasta").write_text(
    "".join(f">{k}\n{v}\n" for k, v in sel.items()))
print(f"Track A prey: {len(sel)}개  (길이>{MAX_PREY_LEN} 제외: {len(skipped)}개)")
if skipped:
    print("제외된 긴 단백질 상위 10:", sorted(skipped, key=lambda x: -x[1])[:10])
pd.DataFrame(skipped, columns=["protein","length"]).to_csv(
    DIR["table"]/"trackA_skipped_long.csv", index=False)
```

---

## CELL 19 — Part 4b. prey homolog 검색 (최장 단계)

```python
# =============================================================================
# CELL 19 | Part 4-2. prey homolog 검색 — 이 파이프라인에서 가장 오래 걸린다 (4~10h)
#   - 반드시 백그라운드. 이 시간 동안 다른 GPU 에서 Track B(CELL 25~)를 먼저 돌려도 된다.
# =============================================================================
need(((DIR["seq"]/"prey_trackA.fasta").exists(), "CELL 18 을 먼저 돌릴 것"),
     (done("part2_bactdb", RB/"bactDB.index"), "CELL 14 미완료 — bactDB 가 완성되지 않았다"))
FMT_ALN = "query,target,fident,evalue,bits,qstart,qend,qlen,tstart,tend,tlen,qaln,taln"  # CELL 16 과 동일해야 함

# ---- 메모리 상한 --------------------------------------------------------
# mmseqs 는 "시작 시점의" 여유 RAM 으로 prefilter 분할 수를 정한다. 그래서 같은
# 머신에서 Boltz(CELL 27)가 병렬로 돌면 나중에 밀려 커널 OOM killer 에 죽는다
# (실측: "[====Killed / Error: Prefilter died"). 상한을 명시해 막는다.
# --db-load-mode 2 = mmap. DB 를 RAM 에 통째로 읽지 않고 페이지 캐시로 두므로
# 메모리 압박 시 커널이 프로세스를 죽이는 대신 캐시를 회수한다.
_avail_gb = int(re.search(r"MemAvailable:\s+(\d+)", open("/proc/meminfo").read()).group(1)) >> 20
SPLIT_MEM = f"{max(8, int(_avail_gb * 0.5))}G"
print(f"MemAvailable = {_avail_gb}G  ->  --split-memory-limit {SPLIT_MEM}")

script = f"""
cd "{DIR['db']}"
rm -f preyDB* "{DIR['search']}"/prey_res*
rm -rf "{DIR['tmp']}/ps"          # 죽은 검색이 남긴 중간 상태를 지운다
mmseqs createdb "{DIR['seq']}/prey_trackA.fasta" preyDB
mmseqs search preyDB "{RB}/bactDB" "{DIR['search']}/prey_res" "{DIR['tmp']}/ps" \\
  -s 7.5 --num-iterations 3 -e 1e-3 --max-seqs 20000 --threads {THREADS} \\
  --split-memory-limit {SPLIT_MEM} --db-load-mode 2
mmseqs convertalis preyDB "{RB}/bactDB" "{DIR['search']}/prey_res" "{DIR['search']}/prey_hits.m8" \\
  --format-output "{FMT_ALN}" --threads {THREADS}
wc -l "{DIR['search']}/prey_hits.m8"
echo DONE_prey
"""
if not already(done("part4_prey", DIR["search"]/"prey_hits.m8"),
               "prey_hits.m8", f"rm {DIR['search']}/prey_hits.m8"):
    sh_bg("part4_prey", script, env=CONDA_ENV_RF2)
print("\n※ Track B 는 이 DB 가 필요 없다 — CELL 25~28 을 지금 병렬로 시작할 수 있다.")
```

---

## CELL 20 — Part 5. organism별 best-hit 정리

```python
# =============================================================================
# CELL 20 | Part 5-1. RF2-PPI README 절차 (2): organism 별 최고 유사도 hit 1개만 남긴다
#   - 정렬을 query 좌표계로 투영해 고정폭 행으로 만든다 (query gap 컬럼은 버림)
#   - prey_hits.m8 이 수 GB 라 chunk 로 읽는다. 수 분~십수 분.
# =============================================================================
need((done("part4_prey", DIR["search"]/"prey_hits.m8"), "CELL 19 미완료"),
     ((DIR["search"]/"prey_hits.m8").exists(), "prey_hits.m8 이 없다"))
from collections import defaultdict

# CELL 16/19 의 --format-output 과 같은 순서여야 한다. CELL 17 에도 같은 정의가
# 있지만, 17 을 건너뛰고 20 만 돌려도 되도록 여기서도 세운다.
AC = ["query","target","fident","evalue","bits","qstart","qend","qlen",
      "tstart","tend","tlen","qaln","taln"]

def to_query_frame(qaln, taln, qstart, qlen):
    row = ["-"] * qlen
    qi = int(qstart) - 1
    for qc, tc in zip(qaln, taln):
        if qc == "-":
            continue
        if qi < qlen:
            row[qi] = tc if tc != "-" else "-"
        qi += 1
    return "".join(row)

def best_hit_per_taxid(m8_path, chunksize=2_000_000, keep_taxids=None):
    """taxid 별 최고 bits hit 하나만 남기고, 그것만 query 좌표계로 투영한다.

    to_query_frame 은 정렬을 한 글자씩 도는 파이썬 루프라 행당 비용이 크다.
    prey_hits.m8 은 실측 36GB / 6,262만 행이어서 전부 투영하면 수 시간이 걸린다.
    그래서 (1) 페어링 불가능한 유전체를 먼저 버리고 (2) bits 비교는 pandas 에서
    끝낸 뒤 (3) 이긴 행에만 투영을 건다.
    """
    best = defaultdict(dict)          # query -> taxid -> (bits, qaln, taln, qstart, qlen)
    n_read = n_kept = 0
    for chunk in pd.read_csv(m8_path, sep="\t", names=AC, chunksize=chunksize):
        n_read += len(chunk)
        chunk["taxid"] = chunk.target.str.split("|", n=1).str[0]
        if keep_taxids is not None:
            chunk = chunk[chunk.taxid.isin(keep_taxids)]
        chunk = (chunk.sort_values("bits", ascending=False)
                      .drop_duplicates(["query", "taxid"]))
        n_kept += len(chunk)
        for r in chunk.itertuples(index=False):
            cur = best[r.query].get(r.taxid)
            if cur is None or r.bits > cur[0]:
                best[r.query][r.taxid] = (r.bits, r.qaln, r.taln, r.qstart, r.qlen)
        print(f"  읽음 {n_read:,} / 후보 {n_kept:,}", end="\r", file=sys.stderr)
    print(f"  읽음 {n_read:,} / 후보 {n_kept:,}", file=sys.stderr)
    return {q: {t: (b, to_query_frame(qa, ta, qs, ql))
                for t, (b, qa, ta, qs, ql) in d.items()}
            for q, d in best.items()}

print("bait hits 정리 중...")
BAIT_MSA = best_hit_per_taxid(DIR["search"]/"bait_hits.m8")
print("bait depth:", {k: len(v) for k, v in BAIT_MSA.items()})

# paired MSA 는 bait 와 prey 가 같은 유전체에 둘 다 있어야 성립한다. bait 오솔로그가
# 없는 유전체의 prey hit 은 CELL 21 에서 어차피 버려지므로 여기서 미리 거른다.
# 실측: bait 는 2,929개 유전체, prey 검색은 17,505개 전체를 훑었다.
PAIRABLE = set().union(*BAIT_MSA.values()) if BAIT_MSA else None
print(f"\n페어링 가능한 유전체 {len(PAIRABLE):,}개로 prey hit 을 제한한다")

print("prey hits 정리 중... (36GB / 6,262만 행 — 십수 분)")
PREY_MSA = best_hit_per_taxid(DIR["search"]/"prey_hits.m8", keep_taxids=PAIRABLE)
print("prey 개수:", len(PREY_MSA))
```

---

## CELL 21 — Part 5b. paired MSA 생성 (백그라운드)

```python
# =============================================================================
# CELL 21 | Part 5-2. paired MSA 생성 — 같은 organism 끼리 이어 붙이고 hhfilter 90%
#   3,730쌍에 hhfilter 를 한 번씩 부르는 1~2시간짜리라 백그라운드로 뺀다.
#   BAIT_MSA/PREY_MSA 는 메모리에만 있어 별도 프로세스가 못 쓰므로, 스크립트가
#   m8 를 다시 읽는다 (15분 추가). 대신 커널이 풀려 다른 셀을 돌릴 수 있다.
#   이미 만들어둔 a3m 은 건너뛰므로 중단 후 재실행해도 이어서 한다.
# =============================================================================
need(((DIR["seq"]/"prey_trackA.fasta").exists(), "CELL 18 을 먼저 돌릴 것"),
     ((DIR["search"]/"prey_hits.m8").exists(), "CELL 19 미완료"),
     (have("BAIT_KEY", "SEGMENTS"), "CELL 07 을 먼저 돌릴 것"))

MIN_PAIRED    = 50
SEGMENT_BAITS = []      # 예: ["ChCODH2_WT_seg440-636"]. 하나 늘 때마다 3,730쌍이 는다.
BAITS_TO_RUN  = [BAIT_KEY] + [b for b in SEGMENT_BAITS if b in SEGMENTS]

_n = len(BAITS_TO_RUN) * 3730
print(f"bait {len(BAITS_TO_RUN)}개 -> 약 {_n:,}쌍")
print(f"  디스크 약 {_n*5/1024:.0f} GB / RF2-PPI 추론 약 {_n*3*4/3600:.0f} 시간(3 replicate)")
need((_n <= 8000, f"{_n:,}쌍은 너무 많다 — SEGMENT_BAITS 를 줄일 것"))

SRC = r"""
import sys, re, subprocess
from pathlib import Path
from collections import defaultdict
import pandas as pd

BASE   = Path("__BASE__")
SEARCH = BASE/"result"/"search"
PAIRED = BASE/"result"/"paired"; PAIRED.mkdir(parents=True, exist_ok=True)
TABLE  = BASE/"result"/"table"
RF2    = BASE/"result"/"rf2ppi";  RF2.mkdir(parents=True, exist_ok=True)
SEQ    = BASE/"input"/"seq"
MIN_PAIRED   = __MINP__
BAIT_KEY     = "__BAITKEY__"
BAITS_TO_RUN = __BAITS__

AC = ["query","target","fident","evalue","bits","qstart","qend","qlen",
      "tstart","tend","tlen","qaln","taln"]

def read_fasta(path):
    out, name, buf = {}, None, []
    for line in open(path, errors="ignore"):
        if line.startswith(">"):
            if name: out[name] = "".join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name: out[name] = "".join(buf)
    return out

def to_query_frame(qaln, taln, qstart, qlen):
    row = ["-"] * qlen
    qi = int(qstart) - 1
    for qc, tc in zip(qaln, taln):
        if qc == "-":
            continue
        if qi < qlen:
            row[qi] = tc if tc != "-" else "-"
        qi += 1
    return "".join(row)

def best_hit_per_taxid(m8, keep=None, chunksize=2_000_000):
    best, n_read, n_kept = defaultdict(dict), 0, 0
    for ch in pd.read_csv(m8, sep="\t", names=AC, chunksize=chunksize):
        n_read += len(ch)
        ch["taxid"] = ch.target.str.split("|", n=1).str[0]
        if keep is not None:
            ch = ch[ch.taxid.isin(keep)]
        ch = ch.sort_values("bits", ascending=False).drop_duplicates(["query","taxid"])
        n_kept += len(ch)
        for r in ch.itertuples(index=False):
            cur = best[r.query].get(r.taxid)
            if cur is None or r.bits > cur[0]:
                best[r.query][r.taxid] = (r.bits, r.qaln, r.taln, r.qstart, r.qlen)
        print(f"  read {n_read:,} / keep {n_kept:,}", flush=True)
    return {q: {t: (b, to_query_frame(qa, ta, qs, ql))
                for t, (b, qa, ta, qs, ql) in d.items()}
            for q, d in best.items()}

baits = read_fasta(SEQ/"baits.fasta")
preys = read_fasta(SEQ/"prey_trackA.fasta")
print(f"bait {len(baits)} / prey {len(preys)}", flush=True)

print("bait hits...", flush=True)
BAIT_MSA = best_hit_per_taxid(SEARCH/"bait_hits.m8")
print("bait depth:", {k: len(v) for k, v in BAIT_MSA.items()}, flush=True)

PAIRABLE = set().union(*BAIT_MSA.values())
print(f"pairable genomes {len(PAIRABLE):,}", flush=True)

print("prey hits...", flush=True)
PREY_MSA = best_hit_per_taxid(SEARCH/"prey_hits.m8", keep=PAIRABLE)
print(f"prey {len(PREY_MSA)}", flush=True)

def make_pair(bkey, bseq, brows, pname, pseq, prows):
    # 이미 만들어둔 것은 건너뛴다 (중단 후 이어서 하기)
    flt = PAIRED / (bkey + "__" + pname.replace("|", "_") + ".a3m")
    if flt.exists():
        return str(flt), sum(1 for l in open(flt) if l.startswith(">"))
    shared = set(brows) & set(prows)
    if len(shared) < MIN_PAIRED:
        return None, len(shared)
    raw = Path(str(flt).replace(".a3m", ".raw.a3m"))
    with open(raw, "w") as fh:
        fh.write(">query\n" + bseq + pseq + "\n")
        for tx in sorted(shared):
            fh.write(">" + str(tx) + "\n" + brows[tx][1] + prows[tx][1] + "\n")
    try:
        subprocess.run(["hhfilter","-i",str(raw),"-o",str(flt),"-id","90","-M","first"],
                       check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print("  [hhfilter fail]", bkey, pname, e.stderr.decode(errors="ignore")[:100], flush=True)
        raw.unlink(missing_ok=True)
        return None, len(shared)
    depth = sum(1 for l in open(flt) if l.startswith(">"))
    raw.unlink()
    return str(flt), depth

CTRL = [k for k in baits
        if k != BAIT_KEY and k.upper().startswith(("COOC","COOT","COOJ","COOF"))
        and k in BAIT_MSA]
print(f"baits_to_run {BAITS_TO_RUN} / controls {CTRL}", flush=True)

lines, stats = [], []
for bkey in BAITS_TO_RUN:
    if bkey not in BAIT_MSA:
        print("  [no bait MSA]", bkey, flush=True); continue
    bseq, brows = baits[bkey], BAIT_MSA[bkey]
    for i, (pid, prows) in enumerate(PREY_MSA.items(), 1):
        if pid not in preys:
            continue
        flt, depth = make_pair(bkey, bseq, brows, pid, preys[pid], prows)
        stats.append((bkey, pid, depth, "skip" if flt is None else "ok"))
        if flt:
            lines.append(flt + " " + str(len(bseq)))
        if i % 200 == 0:
            print(f"  {bkey} {i}/{len(PREY_MSA)}", flush=True)
    print(f"  {bkey} done", flush=True)

# 양성대조군: CooC 를 bait, ChCODH2 를 prey 로. 알려진 상호작용이라 여기서
# 점수가 안 나오면 스크리닝 전체를 신뢰할 수 없다 (CELL 24 가 판정).
for ck in CTRL:
    flt, depth = make_pair(ck, baits[ck], BAIT_MSA[ck],
                           BAIT_KEY, baits[BAIT_KEY], BAIT_MSA[BAIT_KEY])
    stats.append((ck, BAIT_KEY, depth, "skip" if flt is None else "ok"))
    if flt:
        lines.append(flt + " " + str(len(baits[ck])))
    print(f"  control {ck}-{BAIT_KEY}: depth={depth} {'ok' if flt else 'SKIPPED'}", flush=True)
if not CTRL:
    print("  [WARN] no positive control bait", flush=True)

(RF2/"input_file").write_text("\n".join(lines) + "\n")
sdf = pd.DataFrame(stats, columns=["bait","prey","paired_depth","status"])
sdf.to_csv(TABLE/"paired_msa_stats.csv", index=False)
print(f"made {len(lines)} / skipped {(sdf.status=='skip').sum()}", flush=True)
print(sdf.query("status=='ok'").paired_depth.describe().to_string(), flush=True)
"""

SRC = (SRC.replace("__BASE__", str(BASE))
          .replace("__MINP__", str(MIN_PAIRED))
          .replace("__BAITKEY__", BAIT_KEY)
          .replace("__BAITS__", repr(BAITS_TO_RUN)))
sp = DIR["script"]/"part5_paired.py"
sp.write_text(SRC)
print("스크립트:", sp)

if not already(done("part5_paired", DIR["rf2ppi"]/"input_file"),
               "paired MSA", f"rm {DIR['rf2ppi']}/input_file"):
    sh_bg("part5_paired", f'python "{sp}"\necho DONE_paired', env=CONDA_ENV_RF2)
```

---

## CELL 22 — Part 6. RF2-PPI ×N replicate 실행

```python
# =============================================================================
# CELL 22 | Part 6-1. RF2-PPI 추론 — replicate 필수
#   - 비결정적 모델이다. README 기준 ~5% 쌍에서 SD>0.1, 중간 점수대 변동이 크다.
#   - N_REPLICATE 는 CELL 17 gate 결과에 따라 3 또는 5 (깊이 200-500 이면 5회)
# =============================================================================
need(((DIR["rf2ppi"]/"input_file").exists(), "CELL 21 을 먼저 돌릴 것"),
     (have("N_REPLICATE"), "CELL 17 을 먼저 돌릴 것 (replicate 횟수가 gate 결과로 정해진다)"))
reps = " ".join(str(i) for i in range(1, N_REPLICATE + 1))
script = f"""
cd "{DIR['rf2ppi']}"
for rep in {reps}; do
  cp input_file input_rep${{rep}}
  CUDA_VISIBLE_DEVICES={GPU_ID} python "{RF2PPI_DIR}/src/predict_list_PPI.py" \\
      -list_fn input_rep${{rep}} \\
      -model_file "{RF2PPI_DIR}/src/models/RF2-PPI.pt"
  echo "replicate ${{rep}} done"
done
ls -la *.log
echo DONE_rf2ppi
"""
if not already(done("part6_rf2ppi", DIR["rf2ppi"]/"input_rep1.log"), "RF2-PPI replicate 결과",
               f"rm {DIR['rf2ppi']}/input_rep*.log"):
    sh_bg("part6_rf2ppi", script, env=CONDA_ENV_RF2)
print(f"\nreplicate {N_REPLICATE}회. GPU {GPU_ID} 사용. 6~15시간 예상.")
```

---

## CELL 23 — Part 6b. RF2-PPI 결과 집계

```python
# =============================================================================
# CELL 23 | Part 6-2. replicate 집계 -> mean/sd 랭킹
#   컷오프: mean>=0.74 strict(95% precision) / 0.30-0.74 회색지대 / <0.30 배제
# =============================================================================
need((done("part6_rf2ppi", DIR["rf2ppi"]/"input_rep1.log"), "CELL 22 미완료"))
reps_df = []
for rep in range(1, N_REPLICATE + 1):
    log = DIR["rf2ppi"]/f"input_rep{rep}.log"
    if not log.exists():
        print(f"[없음] {log}"); continue
    d = pd.read_csv(log, sep=r"\s+", names=["msa","prob","sec"])
    d = d[pd.to_numeric(d.prob, errors="coerce").notna()].copy()   # 헤더 줄 제거
    d["prob"] = d.prob.astype(float)
    stem = d.msa.apply(lambda p: Path(p).stem)
    d["bait"] = stem.apply(lambda s: s.split("__")[0])
    d["prey"] = stem.apply(lambda s: s.split("__")[-1])
    reps_df.append(d.set_index(["bait","prey"])["prob"].rename(f"rep{rep}"))

assert reps_df, "replicate 로그가 없습니다. CELL 22 가 끝났는지 CELL 15 로 확인하세요."
R = pd.concat(reps_df, axis=1)
R["mean"], R["sd"] = R.mean(axis=1), R.std(axis=1)
R = R.sort_values("mean", ascending=False)

# 양성대조군(prey = ChCODH2)은 스크리닝 랭킹에서 분리한다 — CELL 24 에서 따로 본다
RP     = R.reset_index()
R_ctrl = RP[RP.prey == BAIT_KEY].copy()
R_scr  = RP[RP.prey != BAIT_KEY]

# full-length / segment 를 prey 단위로 합집합(최대값)으로 묶는다
R_best = R_scr.sort_values("mean", ascending=False).groupby("prey").agg(
    best_bait=("bait", "first"), mean=("mean", "max"), sd=("sd", "first")
).sort_values("mean", ascending=False)

meta = cls_prey.set_index("protein")
R_best = R_best.join(meta[["category","pident"]], how="left")
# paired_depth 를 함께 둔다. hhfilter 는 90% 동일성으로 중복을 걷어내므로 필터 후
# 깊이가 크게 떨어질 수 있다 (실측: 공유 유전체는 50개 넘었는데 필터 후 4줄인 쌍이
# 있었다). 얕은 MSA 의 공진화 점수는 못 믿으니 해석 단계에서 이 값을 함께 본다.
# CELL 21 이 백그라운드로 도므로 변수가 아니라 파일에서 읽는다.
_sdf = pd.read_csv(DIR["table"]/"paired_msa_stats.csv")
R_best = R_best.join(_sdf.query("status=='ok'").groupby("prey").paired_depth.max(), how="left")
R_best["desc"] = [hdr2desc.get(i, "")[:70] for i in R_best.index]

R.to_csv(DIR["table"]/"trackA_RF2PPI_all_pairs.csv")
R_best.to_csv(DIR["table"]/"trackA_RF2PPI_ranked.csv")
print("=== 상위 40 후보 (bait 합집합) ===")
print(R_best.head(40).to_string())
print("\nmean>=0.74 (strict):", int((R_best["mean"] >= 0.74).sum()))
print("mean>=0.05 (loose) :", int((R_best["mean"] >= 0.05).sum()))
print("sd>0.1 (재현성 낮음, 재실행 권장):", int((R_best["sd"] > 0.1).sum()))
print("\n※ paired_depth 가 100 미만인 후보는 점수를 그대로 믿지 말 것.")
print(f"   깊이 100 미만 후보: {int((R_best.paired_depth < 100).sum())}개")
```

---

## CELL 24 — Part 6c. 양성대조군 검증 + 진단 플롯

```python
# =============================================================================
# CELL 24 | Part 6-3. 양성대조군(CooC/CooT/CooJ - CooS) 확인 + depth-score 상관 진단
#   - 양성대조군이 0.7 미만이면 paired MSA 깊이 부족 -> 결과 전체 신뢰도 하향
#   - 깊이와 점수가 상관되면 그건 생물학이 아니라 artifact 다
# =============================================================================
print("[해석 가이드]")
print("  mean >= 0.74 : strict. 논문 기준 95% precision. 최우선 검증 대상.")
print("  0.30 ~ 0.74  : 회색지대. Boltz-2 로 재검증 후 판단.")
print("  mean <  0.30 : 배제 (단 paired_depth<200 이면 위음성 가능 -> Track B 로 넘김)\n")

print("[양성대조군 판정] CooC/CooT/CooJ - ChCODH2 (CELL 21 에서 생성한 쌍)")
if len(R_ctrl):
    print(R_ctrl[["bait", "prey", "mean", "sd"]].to_string(index=False))
    top = float(R_ctrl["mean"].max())
    if top >= 0.7:
        print(f"\n  최고 {top:.3f} >= 0.70 -> 파이프라인 검증 통과. 아래 결과를 그대로 해석해도 된다.")
    else:
        print(f"\n  ⚠ 최고 {top:.3f} < 0.70 -> 검증 실패. 원인은 아직 모른다.")
        print("     깊이 부족일 수도, 모델이 이 계를 못 맞히는 것일 수도 있다.")
        print("     CELL 24b 로 둘을 가른 뒤에 해석할 것. 그 전까지 아래 순위표는 인용 금지.")
else:
    print("  ⚠ 대조군 결과가 없다. CELL 05 에 CooC/CooT/CooJ 서열을 넣고")
    print("     CELL 07 -> 16 -> 21 -> 22 를 다시 돌려야 검증된 결과가 된다.")
    print("     이 상태의 순위표는 '검증되지 않은 값'으로 취급할 것.")

d = R_best.dropna(subset=["paired_depth"])
print("\nSpearman r (depth vs score) =",
      round(d[["paired_depth","mean"]].corr(method="spearman").iloc[0,1], 3),
      "  ← 0 에 가까워야 정상. 크면 얕은 MSA 가 점수를 만든 artifact 다.")

# 그림은 있으면 좋고 없어도 그만이다. matplotlib 이 numpy 핀과 충돌해 빠질 수 있어 감싼다.
try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    ax[0].scatter(d.paired_depth, d["mean"], s=6, alpha=.4)
    ax[0].set_xlabel("paired MSA depth"); ax[0].set_ylabel("RF2-PPI prob (mean)")
    ax[0].set_title("depth vs score")   # 라벨은 영문 — 한글 폰트 없으면 네모로 깨진다
    ax[1].hist(R_best["mean"], bins=60); ax[1].set_yscale("log")
    ax[1].set_xlabel("RF2-PPI prob"); ax[1].set_title("score distribution")
    plt.tight_layout(); plt.show()
except Exception as e:
    print(f"[그림 생략] matplotlib 사용 불가: {type(e).__name__}: {e}")
```

---

## CELL 24b — Part 6d. 양성대조군 실패의 원인 분해

```python
# =============================================================================
# CELL 24b | Part 6-4. 대조군이 왜 떨어졌나 — 깊이 탓인가, 모델 탓인가
#   CELL 24 는 0.7 미만이면 "paired MSA 깊이 부족"이라고 출력한다. 그건 가설이지
#   진단이 아니다. 깊이가 충분한데도 낮다면 원인은 다른 데 있고, 그때는 Track A
#   순위표를 후보 지명 근거로 쓸 수 없다. 여기서 그 둘을 가른다.
#     (1) 대조군 쌍의 paired MSA 깊이가 실제로 얕은가
#     (2) 대조군이 스크리닝 전체에서 몇 등인가          <- 결정적인 숫자
#     (3) 대조군보다 높은 쌍과 낮은 쌍의 깊이가 다른가  <- 랭킹이 깊이의 그림자인가
# =============================================================================
need(((DIR["table"]/"trackA_RF2PPI_all_pairs.csv").exists(), "CELL 23 을 먼저 돌릴 것"))
_A    = pd.read_csv(DIR["table"]/"trackA_RF2PPI_all_pairs.csv")
_best = pd.read_csv(DIR["table"]/"trackA_RF2PPI_ranked.csv", index_col=0)
_sd   = pd.read_csv(DIR["table"]/"paired_msa_stats.csv").query("status=='ok'")
_dep  = {(r.bait, r.prey): r.paired_depth for r in _sd.itertuples(index=False)}

# ---- (1) 대조군 쌍의 깊이 ----
_ctrl = _A[_A.prey == BAIT_KEY].copy()
_ctrl["paired_depth"] = [_dep.get((b, p)) for b, p in zip(_ctrl.bait, _ctrl.prey)]
print("=== (1) 양성대조군 쌍 ===")
print(_ctrl[["bait", "prey", "mean", "sd", "paired_depth"]].round(3).to_string(index=False))
_q = _sd.paired_depth
print(f"\n스크리닝 쌍 깊이 분포: 중앙값 {_q.median():.0f} "
      f"(25% {_q.quantile(.25):.0f} / 75% {_q.quantile(.75):.0f} / 최대 {_q.max():.0f})")

# ---- (2) 대조군의 순위 ----
_scr   = _best.dropna(subset=["mean"])
_top   = float(_ctrl["mean"].max())
_above = int((_scr["mean"] > _top).sum())
print(f"\n=== (2) 대조군 순위 ===")
print(f"대조군 최고점 {_top:.3f}")
print(f"이보다 높은 점수를 받은 스크리닝 후보: {_above} / {len(_scr)} 개 "
      f"({100*_above/len(_scr):.1f}%)")
print("알려진 참(true positive)이 무작위 대사효소들 아래에 깔리면, 그 위의 순위는")
print("'붙을 가능성' 순서가 아니다.")

# ---- (3) 랭킹이 깊이의 그림자인가 ----
_hi = _scr[_scr["mean"] > _top].paired_depth.dropna()
_lo = _scr[_scr["mean"] <= _top].paired_depth.dropna()
print(f"\n=== (3) 깊이 비교 ===")
print(f"대조군보다 높은 쌍 {len(_hi)}개: 깊이 중앙값 {_hi.median():.0f}")
print(f"대조군보다 낮은 쌍 {len(_lo)}개: 깊이 중앙값 {_lo.median():.0f}")
print("두 값이 크게 벌어지면 순위는 생물학이 아니라 MSA 깊이를 재고 있는 것이다.")

# ---- 판정 ----
_cd = _ctrl["paired_depth"].max()
print("\n=== 판정 ===")
if pd.isna(_cd):
    print("대조군 깊이를 못 읽었다. paired_msa_stats.csv 의 bait/prey 이름을 확인할 것.")
elif _cd < 200:
    print(f"대조군 깊이 {_cd:.0f} — 얕다. 낮은 점수가 깊이 탓일 수 있다.")
    print("Track A 를 접기 전에 유전체 DB 를 넓혀 깊이를 올리는 쪽을 먼저 시도한다.")
else:
    print(f"대조군 깊이 {_cd:.0f} — 얕지 않다. 깊이로는 설명되지 않는다.")
    print("즉 RF2-PPI 가 이 계에서 '알려진 참'을 위로 올리지 못한다는 뜻이다.")
    print("→ Track A 순위는 후보 지명 근거로 쓸 수 없다. 이 점수를 증거로 인용하지 말 것.")
    print("  상위 40 표는 '계산은 끝났으나 검증 실패' 상태로 남긴다.")
    print("  (되살리려면 대조군을 맞히는 설정 — 다른 bait 절편, 다른 컷오프, 다른")
    print("   MSA 페어링 — 을 먼저 찾아야 한다. 찾기 전 결과는 해석 대상이 아니다.)")
```

---

## CELL 24c — Part 6e. Track A 구제 가능성 검사 (대조군 한 쌍)

```python
# =============================================================================
# CELL 24c | 버리기 전에 '설정 탓'인지 한 번만 확인한다
#   CELL 24b 결론: RF2-PPI 가 알려진 참(CooC1-ChCODH2, 깊이 1451)을 531등 아래에
#   둔다. 깊이 탓이 아니다. 그런데 모델 탓이라고 단정하기 전에, 우리가 준 입력
#   형태 탓일 가능성 두 가지가 남아 있다.
#     (1) 방향  스크리닝은 ChCODH2 가 chain A 인데 대조군은 CooC1 이 chain A 였다.
#               같은 쌍이라도 순서가 바뀌면 점수가 달라질 수 있다.
#     (2) 절편  full-length 630잔기 중 CooC 가 닿는 곳은 일부다. 나머지 500잔기의
#               공진화는 잡음으로 들어간다. C-cluster 주변만 남기면 신호 대 잡음이
#               올라간다. (CELL 17 gate 에서도 C-말단 절편이 깊이를 다 갖고 있었다)
#   둘 다 이미 만들어둔 a3m 을 자르고 붙이는 것으로 끝난다. 새 검색이 필요 없고
#   GPU 도 몇 분이면 된다.
#   판정: 하나라도 0.7 을 넘으면 그 설정으로 스크리닝을 다시 돌릴 가치가 있다
#         (약 12시간). 전부 낮으면 설정 문제가 아니라 모델이 이 계를 못 푸는 것이다.
# =============================================================================
need((have("BAIT_KEY", "BAIT_ALL", "SEGMENTS", "VARIANT_POS"), "CELL 07 을 먼저 돌릴 것"),
     ((DIR["paired"]).exists(), "CELL 21 을 먼저 돌릴 것"))

CTRL_KEY   = next((k for k in BAIT_ALL if k.upper().startswith("COOC")), None)
need((CTRL_KEY is not None, "BAIT_ALL 에 CooC 서열이 없다"))
src = DIR["paired"]/f"{CTRL_KEY}__{BAIT_KEY}.a3m"
need((src.exists(), f"{src} 가 없다 — CELL 21 의 control 생성 로그를 확인할 것"))

cs, bs = BAIT_ALL[CTRL_KEY], BAIT_ALL[BAIT_KEY]
LC, LB = len(cs), len(bs)
print(f"대조군 a3m: {src.name}   {CTRL_KEY} {LC}aa + {BAIT_KEY} {LB}aa = {LC+LB}")

# ---- a3m 읽기 (query 가 첫 줄) ----
names, seqs, nm, buf = [], [], None, []
for l in open(src, errors="ignore"):
    if l.startswith(">"):
        if nm is not None: names.append(nm); seqs.append("".join(buf))
        nm, buf = l[1:].strip(), []
    else: buf.append(l.strip())
if nm is not None: names.append(nm); seqs.append("".join(buf))
_bad = [len(s) for s in seqs if len(s) != LC + LB]
need((not _bad, f"a3m 줄 길이가 {LC+LB} 가 아니다 (예: {_bad[:3]}). 소문자 삽입열이 "
                f"있으면 고정 위치로 못 자른다 — hhfilter -M first 설정을 확인할 것"))
print(f"행 {len(seqs)}개 — 길이 일치 확인")

vdir = DIR["paired"]/"ctrlvar"; vdir.mkdir(exist_ok=True)
for f in vdir.glob("*.a3m"): f.unlink()

def _write(tag, rows, first_len):
    """all-gap 이 된 행을 버리고 hhfilter 로 90% 중복 제거한 뒤 저장한다."""
    raw = vdir/f"{tag}.raw.a3m"
    keep = [(n, s) for n, s in rows
            if n == names[0] or s.count("-") < 0.7 * len(s)]
    raw.write_text("".join(f">{n}\n{s}\n" for n, s in keep))
    out = vdir/f"{tag}.a3m"
    try:
        subprocess.run(["hhfilter", "-i", str(raw), "-o", str(out),
                        "-id", "90", "-M", "first"], check=True, capture_output=True)
    except Exception as e:
        shutil.copy(raw, out); print(f"  [hhfilter 실패 -> 원본 사용] {tag}: {e}")
    raw.unlink(missing_ok=True)
    d = sum(1 for l in open(out) if l.startswith(">"))
    return str(out), first_len, d

# ---- 변형 만들기 ----
# 절편은 변이 위치(C-cluster)를 포함하는 것만 쓴다. 이름에서 좌표를 읽는다.
segs = {}
for sname in SEGMENTS:
    m = re.search(r"_seg(\d+)-(\d+)$", sname)
    if not m: continue
    lo, hi = int(m.group(1)) - 1, int(m.group(2))
    if lo < VARIANT_POS <= hi:
        segs[sname] = (lo, hi)
print(f"C-cluster({VARIANT_POS}) 를 포함하는 절편 {len(segs)}개: {list(segs)}")

jobs = []
jobs.append(_write("v0_base", [(n, s) for n, s in zip(names, seqs)], LC))
jobs.append(_write("v1_swap", [(n, s[LC:] + s[:LC]) for n, s in zip(names, seqs)], LB))
for sname, (lo, hi) in segs.items():
    tag = re.sub(r"[^A-Za-z0-9]+", "_", sname.replace(BAIT_KEY, ""))
    jobs.append(_write(f"v2_seg{tag}",
                       [(n, s[:LC] + s[LC+lo:LC+hi]) for n, s in zip(names, seqs)], LC))
    jobs.append(_write(f"v3_swapseg{tag}",
                       [(n, s[LC+lo:LC+hi] + s[:LC]) for n, s in zip(names, seqs)], hi-lo))

print("\n만든 변형:")
for p, fl, d in jobs:
    print(f"  {Path(p).stem:28s} chainA={fl:4d}  깊이 {d}")

lst = vdir/"input_ctrlvar"
lst.write_text("\n".join(f"{p} {fl}" for p, fl, _ in jobs) + "\n")

# ---- 실행: 변형 수가 적어 전경에서 돌린다 (replicate 3회, 수 분) ----
# GPU 를 이미 누가 쓰고 있으면 붙지 않는다. Boltz 가 14GB 를 쥔 채로 여기에
# RF2-PPI 를 얹으면 둘 다 OOM 으로 죽는다 — 기다렸다 돌리는 편이 빠르다.
CV_GPU = globals().get("GPU_ID", 1)
_mem = subprocess.run(
    "nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv,noheader,nounits",
    shell=True, capture_output=True, text=True).stdout.strip().splitlines()
print("\nGPU 상태:")
_free = {}
for _l in _mem:
    _i, _u, _t = [x.strip() for x in _l.split(",")]
    _free[int(_i)] = int(_t) - int(_u)
    print(f"  GPU {_i}: {_u} / {_t} MiB 사용  (여유 {_free[int(_i)]} MiB)")
if _free.get(CV_GPU, 0) < 8000:
    print(f"\n⚠ GPU {CV_GPU} 의 여유가 {_free.get(CV_GPU, 0)} MiB 뿐이다. 실행하지 않았다.")
    print("   돌고 있는 작업이 끝나기를 기다리거나(CELL 15b), 여유 있는 GPU 로")
    print("   CV_GPU 를 바꿔 이 셀을 다시 실행할 것.")
    print("   ※ GPU 0 은 디스플레이가 물려 있어 RF2-PPI 가 13배 느렸다. 8쌍x3회라")
    print("     20~30분이면 끝나므로 급하면 GPU 0 도 쓸 만하다.")
    raise SystemExit
sh(f'''cd "{vdir}"
for rep in 1 2 3; do
  cp input_ctrlvar in_rep${{rep}}
  CUDA_VISIBLE_DEVICES={CV_GPU} python "{RF2PPI_DIR}/src/predict_list_PPI.py" \\
      -list_fn in_rep${{rep}} -model_file "{RF2PPI_DIR}/src/models/RF2-PPI.pt"
done''', check=False)

# ---- 집계 ----
_r = []
for rep in (1, 2, 3):
    lg = vdir/f"in_rep{rep}.log"
    if not lg.exists(): continue
    for line in open(lg, errors="ignore"):
        f = line.split()
        if len(f) >= 2:
            try: _r.append({"variant": Path(f[0]).stem, "rep": rep, "prob": float(f[1])})
            except ValueError: pass
if _r:
    V = (pd.DataFrame(_r).groupby("variant").prob
           .agg(n="count", mean="mean", sd="std", max="max")
           .sort_values("mean", ascending=False).round(3))
    V["depth"] = [dict((Path(p).stem, d) for p, _, d in jobs).get(i) for i in V.index]
    print("\n=== 변형별 대조군 점수 ===")
    print(V.to_string())
    V.to_csv(DIR["table"]/"trackA_control_variants.csv", encoding="utf-8-sig")
    _best = float(V["mean"].max())
    print("\n=== 판정 ===")
    if _best >= 0.7:
        print(f"{V['mean'].idxmax()} 가 {_best:.3f} — 설정 문제였다.")
        print("→ 이 설정으로 스크리닝을 다시 돌릴 가치가 있다 (CELL 21 의 SEGMENT_BAITS")
        print("   또는 체인 순서를 맞춰 CELL 21 -> 22 재실행, 약 12시간).")
    elif _best >= 0.45:
        _b0 = float(V.loc["v0_base", "mean"]) if "v0_base" in V.index else float("nan")
        print(f"최고 {_best:.3f} — 기준 0.7 에는 못 미치지만 base({_b0:.3f}) 보다 올랐다.")
        print("방향은 맞다. 절편을 더 좁혀 한 번 더 시도할 만하다.")
    else:
        print(f"최고 {_best:.3f} — 방향·절편 어느 쪽으로도 안 오른다.")
        print("  절편이 오히려 더 낮은 것은, C-cluster 주변만 남기면 신호가 농축되는")
        print("  게 아니라 사라진다는 뜻이다.")
        print("\n★ 다만 여기서 Track A 를 닫지 않는다. 이 검사는 전부 '단량체 x 단량체'")
        print("  안에서만 바꿔본 것이다. 올리고머 상태는 한 번도 안 건드렸다.")
        print("  문헌(Merrouch 2018 / Jeoung et al.): CooC1 의 Ni 자리는 각 monomer 의")
        print("  Cys112/Cys114 가 **이량체 계면**에서 만드는 4배위 자리다. 단량체")
        print("  CooC1 에는 그 자리가 존재하지 않는다 — 실재하지 않는 상태를 물어온 것이다.")
        print("  -> CELL 24d 에서 ChCODH2 x2 + CooC1 x2 로 다시 묻는다.")
else:
    print("\n로그를 못 읽었다. vdir 의 in_rep*.log 를 직접 확인할 것:", vdir)
```

---

## CELL 24d — 올리고머 상태로 대조군 재시도 (Boltz 4사슬 + RF2-PPI tandem)

```python
# =============================================================================
# CELL 24d | 단량체가 아니라 실제 올리고머로 묻는다
#   지금까지 대조군은 전부 '단량체 CooC1 × 단량체 ChCODH2' 였다. 그런데
#   Merrouch et al. 2018 / Jeoung et al. 이 보인 것은:
#     CooC1 의 Ni 자리는 각 monomer 의 Cys112/Cys114 가 **이량체 계면**에서
#     만나 만드는 4배위 자리다. 단량체 CooC1 에는 그 자리가 존재하지 않는다.
#     ADP 결합 시 Zn-S(Cys114) 거리가 2.5 -> 2.2 A 로 줄고 CAP loop 가 열린다.
#   즉 우리는 실재하지 않는 상태를 놓고 '붙느냐'고 물어 왔다.
#   Track A/B 대조군 실패(0.266 / ipTM 0.300)의 설명으로 아직 안 지운 가설이다.
#
#   두 가지를 돌린다.
#   (A) Boltz — 4사슬을 네이티브로 지원한다. 이게 진짜 검사다.
#         chain A,B = ChCODH2 x2 (636 x2)   chain C,D = CooC1 x2 (254 x2)
#         ligand    = NI                     합계 1,780 (상한 1,830 안)
#       비교군으로 단량체판(A,C + Ni)도 같이 돌려 같은 축에서 본다.
#   (B) RF2-PPI — 2사슬만 받는다. 각 사슬을 tandem 반복으로 만들어 우회한다.
#         chain1 = [CooC1|CooC1] 508    chain2 = [ChCODH2|ChCODH2] 1272
#       ※ 솔직히: MSA 열을 복제하는 것이라 공진화 정보량은 늘지 않는다.
#         점수가 공진화 통계에서 나오므로 큰 변화는 기대하기 어렵다. 다만
#         기존 a3m 을 자르고 붙이는 것뿐이라 비용이 분 단위다. 안 해볼 이유가 없다.
# =============================================================================
need((have("BAIT_ALL", "BAIT_KEY"), "CELL 07 을 먼저 돌릴 것"))
import yaml as _y
CTRL_KEY = next((k for k in BAIT_ALL if k.upper().startswith("COOC")), None)
need((CTRL_KEY is not None, "BAIT_ALL 에 CooC 서열이 없다"))
bseq, cseq = BAIT_ALL[BAIT_KEY], BAIT_ALL[CTRL_KEY]
print(f"{BAIT_KEY} {len(bseq)}aa · {CTRL_KEY} {len(cseq)}aa")
print(f"이량체 조합 총 토큰 {2*len(bseq) + 2*len(cseq)}  (상한 1830)")

# ---------- (A) Boltz ----------
DIMER_GPU  = 1
REP_N      = 5
BOLTZ_MAXMSA = globals().get("BOLTZ_MAXMSA", 2048)
din = DIR["boltz"]/"inputs_dimer"; din.mkdir(parents=True, exist_ok=True)
for f in din.glob("*.yaml"): f.unlink()

def wr(name, chains, ni=True):
    seqs = [{"protein": {"id": cid, "sequence": s}} for cid, s in chains]
    if ni: seqs.append({"ligand": {"id": "Z", "ccd": "NI"}})
    tot = sum(len(s) for _, s in chains)
    if tot > 1830:
        print(f"  길이 초과 제외 {name} ({tot})"); return 0
    (din/f"{name}.yaml").write_text(_y.safe_dump({"version": 1, "sequences": seqs},
                                                 sort_keys=False))
    print(f"  {name}: {len(chains)}사슬 {tot} 토큰")
    return 1

n = 0
n += wr("CTRL-dimer2x2", [("A", bseq), ("B", bseq), ("C", cseq), ("D", cseq)])
n += wr("CTRL-mono1x1",  [("A", bseq), ("C", cseq)])
n += wr("CTRL-codhdimer", [("A", bseq), ("B", bseq), ("C", cseq)])
n += wr("APO-codhdimer",  [("A", bseq), ("B", bseq)])       # 파트너 없이 Ni 만
# 1순위 후보도 같은 축에서
for pid in ["QJZ12568.1"]:
    _s = None
    for l in open(ASSET["faa_bl21"], errors="ignore"):
        pass
    seqs_bl = {}
    nm, buf = None, []
    for l in open(ASSET["faa_bl21"], errors="ignore"):
        if l.startswith(">"):
            if nm: seqs_bl[nm] = "".join(buf)
            nm, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if nm: seqs_bl[nm] = "".join(buf)
    if pid in seqs_bl:
        n += wr(f"CAND-{pid}-codhdimer", [("A", bseq), ("B", bseq), ("C", seqs_bl[pid])])

print(f"\n입력 {n}개 -> {din}")
script = f"""
cd "{DIR['boltz']}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
CUDA_VISIBLE_DEVICES={DIMER_GPU} boltz predict inputs_dimer \\
  --out_dir out_dimer --use_msa_server --max_msa_seqs {BOLTZ_MAXMSA} \\
  --recycling_steps 3 --diffusion_samples {REP_N} \\
  --output_format mmcif --num_workers 2
echo DONE_boltz_dimer
"""
_busy = subprocess.run("pgrep -f '[b]oltz predict'", shell=True,
                       capture_output=True, text=True).stdout.split()
if _busy:
    print(f"⚠ boltz 가 이미 돌고 있다 (pid {' '.join(_busy)}). 안 띄웠다.")
elif n and not already(done("part8_boltz_dimer", DIR["boltz"]/"out_dimer"),
                       "이량체 대조군", f"rm -rf {DIR['boltz']}/out_dimer"):
    sh_bg("part8_boltz_dimer", script, env=CONDA_ENV_BOLTZ)

# ---------- (B) RF2-PPI tandem ----------
src = DIR["paired"]/f"{CTRL_KEY}__{BAIT_KEY}.a3m"
if not src.exists():
    print(f"\n[건너뜀] {src} 없음")
else:
    LC, LB = len(cseq), len(bseq)
    names, seqs, nm, buf = [], [], None, []
    for l in open(src, errors="ignore"):
        if l.startswith(">"):
            if nm is not None: names.append(nm); seqs.append("".join(buf))
            nm, buf = l[1:].strip(), []
        else: buf.append(l.strip())
    if nm is not None: names.append(nm); seqs.append("".join(buf))
    ok = all(len(s) == LC + LB for s in seqs)
    print(f"\nRF2-PPI tandem: a3m {len(seqs)}행, 길이 일치 {ok}")
    if ok:
        vdir = DIR["paired"]/"dimer"; vdir.mkdir(exist_ok=True)
        raw = vdir/"tandem.raw.a3m"
        raw.write_text("".join(
            f">{n_}\n{s[:LC]}{s[:LC]}{s[LC:]}{s[LC:]}\n" for n_, s in zip(names, seqs)))
        out = vdir/"tandem.a3m"
        try:
            subprocess.run(["hhfilter", "-i", str(raw), "-o", str(out),
                            "-id", "90", "-M", "first"], check=True, capture_output=True)
        except Exception as e:
            shutil.copy(raw, out); print(f"  [hhfilter 실패 -> 원본] {e}")
        d = sum(1 for l in open(out) if l.startswith(">"))
        print(f"  tandem a3m: chainA={2*LC} chainB={2*LB} 깊이 {d}")
        (vdir/"input_tandem").write_text(f"{out} {2*LC}\n")
        CV = globals().get("GPU_ID", 1)
        sh(f'''cd "{vdir}"
for rep in 1 2 3; do
  cp input_tandem in_t${{rep}}
  CUDA_VISIBLE_DEVICES={CV} python "{RF2PPI_DIR}/src/predict_list_PPI.py" \\
      -list_fn in_t${{rep}} -model_file "{RF2PPI_DIR}/src/models/RF2-PPI.pt"
done''', check=False)
        vals = []
        for rep in (1, 2, 3):
            lg = vdir/f"in_t{rep}.log"
            if lg.exists():
                for line in open(lg, errors="ignore"):
                    f = line.split()
                    if len(f) >= 2:
                        try: vals.append(float(f[1]))
                        except ValueError: pass
        if vals:
            import statistics as _st
            print(f"\n=== RF2-PPI tandem 이량체 대조군 ===")
            print(f"  mean {_st.mean(vals):.3f}  sd {_st.pstdev(vals):.3f}  n {len(vals)}")
            print(f"  단량체판 기준: 0.251 (CELL 24c v0_base)")
            if _st.mean(vals) >= 0.7:
                print("  → 올리고머 상태가 문제였다. Track A 를 이 설정으로 되살린다.")
            elif _st.mean(vals) > 0.35:
                print("  → 올랐다. 방향이 맞다. 더 밀어볼 가치가 있다.")
            else:
                print("  → 안 오른다. 예상대로 MSA 열 복제로는 공진화 정보가 늘지 않는다.")
                print("     Track A 의 판단은 Boltz 4사슬 결과(A)를 보고 내린다.")
print("\n다음: CELL 15b 로 part8_boltz_dimer 확인 -> 끝나면 CELL 28a-13 방식으로 파싱")
```

---

## CELL 25 — Part 7. Boltz-2 설치

```python
# =============================================================================
# CELL 25 | Part 7-1. Boltz-2 설치 스크립트 (Track B — 여기가 주력일 가능성이 높다)
#   - 공진화 불필요. strain-specific 217 + low-sim 125 = 342개가 대상.
#   - Ni 을 ligand(CCD: NI)로 명시할 수 있는 게 이 문제에서 결정적 장점.
# =============================================================================
# pip install boltz 는 cuEquivariance 가속 커널을 안 딸려온다. 그게 없으면
# triangular_mult 에서 ModuleNotFoundError 로 첫 쌍부터 죽는다. 커널을 끄면
# (--no_kernels) 돌긴 하지만 중간 텐서를 전부 물어서 쌍당 VRAM 이 23.8GB 로
# 뛰고, 3090(24.5GB)에서 15% 가 "ran out of memory, skipping batch" 로
# 조용히 버려진다. 그래서 커널은 선택이 아니라 필수다.
#   ops 휠은 torch 의 CUDA 메이저 버전과 맞춰야 한다 (cu13 빌드면 -cu13).
boltz_sh = f"""
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -n {CONDA_ENV_BOLTZ} python=3.11
conda activate {CONDA_ENV_BOLTZ}
pip install boltz -U
pip install pyyaml pandas

# torch 가 어느 CUDA 로 빌드됐는지 보고 맞는 ops 휠을 고른다
CU=$(python -c "import torch;print(torch.version.cuda.split('.')[0])")
echo "torch CUDA major = $CU"
pip install cuequivariance-torch "cuequivariance-ops-torch-cu${{CU}}"

python -c "import boltz; print('boltz ok')"
python -c "from cuequivariance_torch.primitives.triangle import triangle_multiplicative_update; print('kernel ok')"
"""
p = DIR["script"] / "install_boltz.sh"
p.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + textwrap.dedent(boltz_sh))
p.chmod(0o755)
print(f"터미널에서 실행:\n\n    bash {p}\n")
print(p.read_text())
```

---

## CELL 26 — Part 7b. Boltz 입력 YAML 생성 (Ni ligand 포함)

```python
# =============================================================================
# CELL 26 | Part 7-2. Track B 입력 YAML 생성 — Ni 을 ligand 로 명시
#   - bait 630 + prey 1200 = 1830 잔기가 24GB 상한 근처. 더 긴 건 제외.
# =============================================================================
need((have("TRACK_B", "hdr2seq", "BAIT_ALL"), "CELL 08/12 와 CELL 07 을 먼저 돌릴 것"))
import yaml
boltz_in = DIR["boltz"]/"inputs"; boltz_in.mkdir(parents=True, exist_ok=True)
MAX_B_LEN = 1200

def write_boltz_yaml(name, seqA, seqB, msaA=None, msaB=None, with_ni=True):
    seqs = [
        {"protein": {"id": "A", "sequence": seqA, **({"msa": str(msaA)} if msaA else {})}},
        {"protein": {"id": "B", "sequence": seqB, **({"msa": str(msaB)} if msaB else {})}},
    ]
    if with_ni:
        seqs.append({"ligand": {"id": "C", "ccd": "NI"}})
    p = boltz_in/f"{name}.yaml"
    p.write_text(yaml.safe_dump({"version": 1, "sequences": seqs}, sort_keys=False))
    return p

# Track A 는 segment bait 를 쓰는데 Track B 는 기본적으로 full-length 만 쓴다.
# 설계문서상 segment 가 더 민감할 수 있지만 GPU 시간이 segment 개수만큼 배로 든다
# (342쌍 x 7 segment = 2,400쌍, 수일). Feasibility 를 보고 True 로 켠다.
BOLTZ_USE_SEGMENT = False
BOLTZ_BAITS = [BAIT_KEY] + (list(SEGMENTS) if BOLTZ_USE_SEGMENT else [])

n, long_skip = 0, []
for bk in BOLTZ_BAITS:
    bseq = BAIT_ALL[bk]
    for pid in TRACK_B:
        s = hdr2seq[pid]
        if len(bseq) + len(s) > MAX_B_LEN + 630:
            if bk == BAIT_KEY:
                long_skip.append((pid, len(s)))
            continue
        write_boltz_yaml(f"{bk}__{pid.replace('|','_')}", bseq, s, with_ni=True)
        n += 1
print(f"Boltz-2 입력 {n}개 생성 (bait {len(BOLTZ_BAITS)}종) -> {boltz_in}")
print(f"길이 초과 제외: {len(long_skip)}개  {long_skip[:5]}")
```

---

## CELL 27 — Part 7c. Boltz-2 실행

```python
# =============================================================================
# CELL 27 | Part 7-3. Boltz-2 추론 (GPU). 342쌍 기준 10~20시간.
#   - MSA 는 ColabFold 서버 사용(--use_msa_server). 수백 건이면 서버 부하가 크니
#     가능하면 Part 4 의 mmseqs MSA 를 a3m 으로 변환해 yaml 의 msa: 필드로 넣는 게 낫다.
#   - Track A 가 다른 GPU 를 쓰고 있으면 GPU_ID 를 1 로 바꿔 병렬 실행.
# =============================================================================
need((len(list((DIR["boltz"]/"inputs").glob("*.yaml"))) > 0, "CELL 26 을 먼저 돌릴 것"))
BOLTZ_GPU   = GPU_ID
BOLTZ_MAXMSA = 2048     # 기본 8192. VRAM 을 가장 크게 좌우한다 (아래 설명)

# ---- VRAM 튜닝 근거 (RTX 3090 24.5GB, bait 636 + prey <=1200 = 최대 1836 토큰) ----
# 실측 3단계:
#   1) --no_kernels (커널 없이)      쌍당 23.8GB -> 15% 가 OOM 으로 버려짐
#   2) cuEquivariance 커널 켬         쌍당 20.2GB -> 그래도 일부 OOM
#   3) + --max_msa_seqs 2048          MSA 텐서가 1/4 로
# (2)에서 실패한 할당이 4.24GB 였는데 MSA 표현 크기와 일치한다:
#   8192(깊이) x 1836(토큰) x 64(채널) x 4바이트 = 3.85GB  -> 2048 이면 0.96GB
# Boltz 는 OOM 이 나면 예외를 던지지 않고 그 쌍을 버린다
# ("WARNING: ran out of memory, skipping batch"). 결과만 조용히 줄어드니
# 끝난 뒤 반드시 입력 개수와 산출 개수를 대조할 것 (CELL 28 이 한다).

script = f"""
cd "{DIR['boltz']}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True   # 할당자 단편화 완화
CUDA_VISIBLE_DEVICES={BOLTZ_GPU} boltz predict inputs \\
  --out_dir out \\
  --use_msa_server \\
  --max_msa_seqs {BOLTZ_MAXMSA} \\
  --recycling_steps 3 \\
  --diffusion_samples 1 \\
  --output_format mmcif \\
  --num_workers 2
echo DONE_boltz
"""

# out/ 디렉터리는 boltz 가 시작하자마자 만든다. 중간에 죽어도 남기 때문에
# 디렉터리 존재만으로 "완료"로 보면 재실행이 영영 건너뛰어진다.
# 실제로 나온 구조 개수를 세서 판정한다.
_n_in  = len(list((DIR["boltz"]/"inputs").glob("*.yaml")))
_n_out = len(list((DIR["boltz"]/"out").rglob("*_model_0.cif")))
print(f"입력 {_n_in}개 / 예측 완료 {_n_out}개")

# ---- BOLTZ_MAXMSA 가 실제로 반영되는지 검사 -------------------------------
# 함정: --max_msa_seqs 는 추론이 아니라 "전처리" 단계에서 적용된다. a3m 을 읽어
# processed/ 밑에 npz 로 구울 때 깊이를 자른다. 그래서 이전 전처리 결과가 남아
# 있으면 ("All inputs are already processed.") 값을 바꿔도 조용히 무시되고
# VRAM 도 안 줄어든다. 실측으로 한 번 당했다.
_proc = DIR["boltz"]/"out"/"boltz_results_inputs"/"processed"
_mark = DIR["boltz"]/".maxmsa"
_prev = _mark.read_text().strip() if _mark.exists() else None
if _proc.exists() and _prev != str(BOLTZ_MAXMSA):
    _was = f"max_msa_seqs={_prev}" if _prev else "알 수 없는 max_msa_seqs"
    print(f"\n⚠ 기존 전처리({_was})가 남아 있어 BOLTZ_MAXMSA={BOLTZ_MAXMSA} 이 무시된다.")
    print(f"  적용하려면 전처리만 지울 것:  rm -rf {_proc}")
    print(f"  (형제 디렉터리 msa/ 는 ColabFold 에서 받은 a3m 이다. 남겨두면 서버에 다시 안 물어본다)")
    print(f"  깊이를 통일하려면 predictions/ 도 함께 지운다 — 섞이면 ipTM 비교가 불공정해진다.")

# 이미 돌고 있는 프로세스가 있는데 또 띄우면 같은 GPU 에 두 개가 붙어 즉시 OOM 이다.
# 패턴을 "[b]oltz" 로 쓰는 이유: shell=True 가 띄우는 sh 의 cmdline 에도 검색어가
# 그대로 들어가서, 평범하게 쓰면 pgrep 이 자기 자신을 잡아 항상 "돌고 있다"가 된다.
# 대괄호는 정규식으로는 boltz 와 매칭되지만 문자열로는 달라 자기 매칭을 피한다.
_busy = subprocess.run("pgrep -f '[b]oltz predict'", shell=True,
                       capture_output=True, text=True).stdout.split()
if _busy:
    print(f"⚠ 이미 boltz 가 돌고 있다 (pid {' '.join(_busy)}). 새로 띄우지 않았다.")
    print("  설정을 바꿔 다시 돌리려면 먼저:  pkill -9 -f 'boltz predict'")
    print("  (Boltz 는 이미 끝낸 쌍을 건너뛰므로 재시작 비용은 작다)")
elif not already(done("part7_boltz") or (_n_in > 0 and _n_out >= _n_in),
                 "Boltz-2 예측 결과", f"rm -rf {DIR['boltz']}/out"):
    _mark.write_text(str(BOLTZ_MAXMSA))   # 다음 실행 때 깊이가 바뀌었는지 비교하려고
    sh_bg("part7_boltz", script, env=CONDA_ENV_BOLTZ)
```

---

## CELL 28 — Part 7d. Boltz-2 결과 파싱

```python
# =============================================================================
# CELL 28 | Part 7-4. confidence json 파싱 -> ipTM / ligand_ipTM 랭킹
#   ipTM  >=0.80 강한 복합체 / 0.60-0.80 유의 / 0.40-0.60 회색 / <0.40 배제
#   ligand_ipTM: Ni 이 두 사슬 경계면에 놓였는가 -> 'Ni 전달 복합체' 시나리오와 부합
#   ⚠ 설계문서: insertase 는 apo 중간체에 붙을 수 있으니 컷오프를 낮춰 수동 검토할 것
# =============================================================================
need((done("part7_boltz", DIR["boltz"]/"out"), "CELL 27 미완료"))
rows = []
for f in glob.glob(str(DIR["boltz"]/"out"/"**"/"confidence_*.json"), recursive=True):
    c = json.load(open(f))
    name = re.sub(r"^confidence_|_model_\d+$", "", Path(f).stem)
    rows.append({"pair": name, "prey": name.split("__")[-1],
                 "iptm": c.get("iptm"), "ptm": c.get("ptm"),
                 "complex_plddt": c.get("complex_plddt"),
                 "ligand_iptm": c.get("ligand_iptm"),
                 "confidence_score": c.get("confidence_score")})
B = pd.DataFrame(rows).sort_values("iptm", ascending=False)
if len(B):
    B["category"] = B.prey.map(cls_prey.set_index("protein")["category"])
    B["desc"] = B.prey.map(lambda i: hdr2desc.get(i, "")[:70])
B.to_csv(DIR["table"]/"trackB_boltz2_ranked.csv", index=False)
print(f"파싱 {len(B)}건")

# Boltz 는 VRAM 이 모자라면 예외 없이 그 쌍을 버린다 ("skipping batch").
# 결과만 조용히 줄어드므로 입력과 대조해 누락을 드러낸다.
_want = {f.stem for f in (DIR["boltz"]/"inputs").glob("*.yaml")}
_miss = sorted(_want - set(B.pair)) if len(B) else sorted(_want)
if _miss:
    print(f"\n⚠ 누락 {len(_miss)}/{len(_want)}건 — OOM 으로 버려졌을 가능성이 높다.")
    print("  로그에서 확인:  grep -c 'skipping batch' result/log/part7_boltz.log")
    print(f"  회수: CELL 27 의 BOLTZ_MAXMSA 를 1024 로 낮춰 다시 돌린다")
    print("        (이미 끝난 쌍은 건너뛰므로 누락분만 계산한다)")
    print("  누락 예시:", _miss[:5])
else:
    print("누락 없음 — 입력 전량이 예측됐다.")

print(B.head(40).to_string(index=False))
```

---

## CELL 28a — Part 7d-1. Track B 랭킹의 길이 편향 점검

```python
# =============================================================================
# CELL 28a | ipTM 랭킹이 결합을 재는가, 길이를 재는가
#   ipTM 은 계면 잔기 기준이라 짧은 prey 일수록 부풀려진다. 20잔기 펩타이드가
#   636잔기 bait 표면에 닿기만 해도 계면 대부분이 "맞은" 것으로 잡히기 때문이다.
#   실측(BL21 321쌍): ipTM ρ=-0.33 / ligand_iptm ρ=-0.81 / confidence ρ=-0.57
#   -> ligand_iptm 은 랭킹 축으로 못 쓴다. "Ni 이 계면에 있나"가 아니라
#      "복합체가 작은가"를 재고 있다.
# =============================================================================
need((have("hdr2seq"), "CELL 08 을 먼저 돌릴 것"),
     ((DIR["table"]/"trackB_boltz2_ranked.csv").exists(), "CELL 28 을 먼저 돌릴 것"))
B = pd.read_csv(DIR["table"]/"trackB_boltz2_ranked.csv")
B["prey_len"] = B.prey.map(lambda p: len(hdr2seq.get(p, "")))

print("prey 길이:", B.prey_len.describe()[["min","25%","50%","75%","max"]].astype(int).to_dict())
for c in ["iptm", "ligand_iptm", "confidence_score"]:
    if c in B.columns:
        r = B[[c, "prey_len"]].corr(method="spearman").iloc[0, 1]
        print(f"  {c:18s} vs 길이  Spearman rho = {r:+.3f}")

B["len_bin"] = pd.cut(B.prey_len, [0, 50, 100, 200, 400, 10**5])
print("\n길이 구간별:")
print(B.groupby("len_bin", observed=True)
       .agg(n=("iptm","size"), iptm_평균=("iptm","mean"), iptm_최대=("iptm","max"),
            ligand_평균=("ligand_iptm","mean")).round(3).to_string())

print("\n=== 100잔기 이상만 추린 상위 20 ===")
cols = ["prey","prey_len","iptm","ligand_iptm","complex_plddt","category","desc"]
print(B[B.prey_len >= 100].nlargest(20, "iptm")[cols].to_string(index=False))

B.to_csv(DIR["table"]/"trackB_boltz2_ranked.csv", index=False)
print("\n길이 컬럼을 추가해 다시 저장했다.")
```

---

## CELL 28a-2 — Part 7d-2. 길이 보정 재랭킹 + Track C 교차

```python
# =============================================================================
# CELL 28a-2 | 구간 내 z-점수로 재랭킹하고 구조 모티프와 교차한다
#   ipTM 을 길이 구간 안에서 표준화하면 구간 간 비교가 된다.
#   metal/ATP 모티프(Track C)는 구조 기하 기반이라 길이 편향이 없는 독립 축이다.
#   ⚠ 짧고 pLDDT 낮은 단백질은 두 축을 동시에 속인다. 겹쳤다고 곧바로 후보로
#     올리지 말고 CELL 28a-3 에서 nres/plddt 를 확인할 것.
# =============================================================================
need((have("B"), "CELL 28a 를 먼저 돌릴 것"))
fd_metal = pd.read_csv(DIR["table"]/"folddisco_metal_motif.csv")
fd_atp   = pd.read_csv(DIR["table"]/"folddisco_atp_motif.csv")
Mset, Aset = set(fd_metal.protein), set(fd_atp.protein)

B["iptm_z"]      = B.groupby("len_bin", observed=True).iptm.transform(
                       lambda x: (x - x.mean()) / (x.std(ddof=0) or 1))
B["metal_motif"] = B.prey.isin(Mset)
B["atp_motif"]   = B.prey.isin(Aset)

print(f"Track B {len(B)}개 중")
print(f"  metal 모티프 보유 : {int(B.metal_motif.sum())}")
print(f"  ATP  모티프 보유 : {int(B.atp_motif.sum())}")
print(f"  둘 다            : {int((B.metal_motif & B.atp_motif).sum())}")
print("  ※ ATP 보유율이 프로테옴 전체 비율(약 23%)과 비슷하면 농축이 없는 것이다 —")
print("     Walker A 는 P-loop NTPase 전체가 공유해 단독 필터로는 무력하다.")

cols = ["prey","prey_len","iptm","iptm_z","metal_motif","atp_motif","category","desc"]
print("\n=== 길이 보정 상위 20 (구간 내 z-점수) ===")
print(B.nlargest(20, "iptm_z")[cols].to_string(index=False))

print("\n=== metal 모티프를 가진 Track B 후보 ===")
sub = B[B.metal_motif].sort_values("iptm", ascending=False)
print(sub[cols].to_string(index=False) if len(sub) else "  없음")

B.to_csv(DIR["table"]/"trackB_boltz2_ranked.csv", index=False)
```

---

## CELL 28a-3 — Part 7d-3. metal 모티프 보유 단백질이 어느 Track 에 있는가

```python
# =============================================================================
# CELL 28a-3 | 금속 배위 기하를 가진 단백질의 소속 확인
#   Track B(균주 특이)에 거의 없고 Track A(보존)에 몰려 있다면, 균주 차이는
#   유전자 유무가 아니라 발현량/조절 문제로 넘어간다.
#   실측(BL21 29개): A 27 / B 2 — 설계문서의 "Track B 가 주력" 예측과 어긋난다.
#   pLDDT 순으로 정렬한다. 70 미만은 AlphaFold 가 접힘을 확신하지 못한 영역이라
#   그 좌표 위에서 잰 모티프 기하를 근거로 쓰기 어렵다.
# =============================================================================
need((have("TRACK_A", "TRACK_B", "hdr2desc"), "CELL 08 과 CELL 12 를 먼저 돌릴 것"))
fdd = pd.read_csv(DIR["table"]/"folddisco_metal_motif_detail.csv")
bl  = fdd[fdd.strain == "BL21"].dropna(subset=["protein"])
A, Bs = set(TRACK_A), set(TRACK_B)

bl = bl.assign(track=bl.protein.map(lambda x: "A(보존)" if x in A else
                                              "B(균주특이)" if x in Bs else "어느 쪽도 아님"))
print(bl.track.value_counts().to_string())

cols = ["protein","track","idf","min_rmsd","nres","plddt","matching_residues"]
out = bl.sort_values("plddt", ascending=False)[cols].copy()
out["desc"] = out.protein.map(lambda x: hdr2desc.get(x, "")[:55])
print()
print(out.to_string(index=False))

print("\n※ pLDDT 70 미만은 근거로 쓰지 말 것. 짧고 못 접힌 단백질은 folddisco idf 와")
print("   ipTM 을 동시에 부풀리므로, 독립적인 두 증거처럼 보여도 같은 아티팩트다.")
out.to_csv(DIR["table"]/"metal_motif_by_track.csv", index=False)
print("저장:", DIR["table"]/"metal_motif_by_track.csv")
```

---

## CELL 28a-4 — Part 7d-4. metal 모티프 보유 단백질을 Boltz 로 직접 검증

```python
# =============================================================================
# CELL 28a-4 | 금속 배위 기하를 가진 단백질에 co-folding 증거를 붙인다
#   Boltz 는 Track B 에만 돌렸던 탓에, 정작 metal 모티프를 가진 쪽(대부분 Track A)이
#   빠졌다. 수십 개뿐이라 GPU 30분이면 끝난다 — RF2-PPI 를 기다릴 필요가 없다.
#   pLDDT 70 미만은 좌표를 못 믿으니 제외한다.
#   판정 기준: 같은 길이 구간의 Track B 분포(CELL 28a 표)보다 뚜렷이 위인가.
# =============================================================================
need((have("BAIT_KEY", "BAIT_ALL", "hdr2seq"), "CELL 07 과 CELL 08 을 먼저 돌릴 것"))
import yaml as _y
BOLTZ_MAXMSA = globals().get("BOLTZ_MAXMSA", 2048)   # CELL 27 을 안 거쳤을 때 대비
FOCUS_PLDDT  = 70

fdd = pd.read_csv(DIR["table"]/"folddisco_metal_motif_detail.csv")
focus = (fdd[(fdd.strain == "BL21") & (fdd.plddt >= FOCUS_PLDDT)]
           .dropna(subset=["protein"]).protein.unique().tolist())
print(f"대상 {len(focus)}개 (pLDDT >= {FOCUS_PLDDT})")

fin = DIR["boltz"]/"inputs_focus"; fin.mkdir(parents=True, exist_ok=True)
for f in fin.glob("*.yaml"): f.unlink()
n = 0
for pid in focus:
    sq = hdr2seq.get(pid)
    if not sq:
        print(f"  [서열 없음] {pid}"); continue
    (fin/f"{BAIT_KEY}__{pid}.yaml").write_text(_y.safe_dump({
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": BAIT_ALL[BAIT_KEY]}},
            {"protein": {"id": "B", "sequence": sq}},
            {"ligand":  {"id": "C", "ccd": "NI"}},
        ]}, sort_keys=False))
    n += 1
print(f"입력 {n}개 -> {fin}")

script = f"""
cd "{DIR['boltz']}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
CUDA_VISIBLE_DEVICES={GPU_ID} boltz predict inputs_focus \\
  --out_dir out_focus \\
  --use_msa_server \\
  --max_msa_seqs {BOLTZ_MAXMSA} \\
  --recycling_steps 3 \\
  --diffusion_samples 1 \\
  --output_format mmcif \\
  --num_workers 2
echo DONE_boltz_focus
"""
_busy = subprocess.run("pgrep -f '[b]oltz predict'", shell=True,
                       capture_output=True, text=True).stdout.split()
if _busy:
    print(f"⚠ boltz 가 이미 돌고 있다 (pid {' '.join(_busy)}). 새로 띄우지 않았다.")
elif not already(done("part7_boltz_focus", DIR["boltz"]/"out_focus"),
                 "focus 예측 결과", f"rm -rf {DIR['boltz']}/out_focus"):
    sh_bg("part7_boltz_focus", script, env=CONDA_ENV_BOLTZ)
```

---

## CELL 28a-5 — Part 7d-5. focus 예측 판정 (Track B 분포 대비)

```python
# =============================================================================
# CELL 28a-5 | metal 모티프 단백질의 ipTM 이 정말 높은가
#   ipTM 절대값은 길이에 따라 기준선이 달라서 그대로 비교하면 안 된다.
#   Track B 321쌍을 배경 분포로 삼아, 같은 길이 구간 안에서 몇 번째인지를 본다.
#   백분위 90 이상이면 "같은 길이의 무작위 균주특이 단백질보다 뚜렷이 위"다.
# =============================================================================
need((have("hdr2seq"), "CELL 08 을 먼저 돌릴 것"),
     ((DIR["boltz"]/"out_focus").exists(), "CELL 28a-4 를 먼저 돌릴 것"),
     ((DIR["table"]/"trackB_boltz2_ranked.csv").exists(), "CELL 28 을 먼저 돌릴 것"))

BINS = [0, 50, 100, 200, 400, 10**5]
rows = []
for f in glob.glob(str(DIR["boltz"]/"out_focus"/"**"/"confidence_*.json"), recursive=True):
    c = json.load(open(f))
    name = re.sub(r"^confidence_|_model_\d+$", "", Path(f).stem)
    rows.append({"pair": name, "prey": name.split("__")[-1],
                 "iptm": c.get("iptm"), "ptm": c.get("ptm"),
                 "complex_plddt": c.get("complex_plddt"),
                 "ligand_iptm": c.get("ligand_iptm")})
F = pd.DataFrame(rows)
want = len(list((DIR["boltz"]/"inputs_focus").glob("*.yaml")))
print(f"입력 {want} / 파싱 {len(F)}")
if len(F) < want:
    print(f"  ⚠ {want - len(F)}건 누락 — CELL 15 로 part7_boltz_focus 가 끝났는지 확인할 것")
if not len(F):
    raise SystemExit

# 배경 분포: Track B. CSV 를 거치면 len_bin 이 문자열이 되므로 길이에서 다시 만든다.
Bg = pd.read_csv(DIR["table"]/"trackB_boltz2_ranked.csv")
if "prey_len" not in Bg.columns:
    Bg["prey_len"] = Bg.prey.map(lambda x: len(hdr2seq.get(x, "")))
Bg["len_bin"] = pd.cut(Bg.prey_len, BINS)

F["prey_len"] = F.prey.map(lambda x: len(hdr2seq.get(x, "")))
F["len_bin"]  = pd.cut(F.prey_len, BINS)

def percentile_in_bin(r):
    ref = Bg.loc[Bg.len_bin == r.len_bin, "iptm"].dropna()
    return round(100.0 * (ref < r.iptm).mean(), 1) if len(ref) else float("nan")
F["pct_vs_trackB"] = F.apply(percentile_in_bin, axis=1)

# 구조 모티프 쪽 근거를 붙인다
fdd = pd.read_csv(DIR["table"]/"folddisco_metal_motif_detail.csv")
fdd = (fdd[fdd.strain == "BL21"].dropna(subset=["protein"])
         .sort_values("idf", ascending=False).drop_duplicates("protein"))
F = F.merge(fdd[["protein","idf","min_rmsd","nres","plddt"]],
            left_on="prey", right_on="protein", how="left").drop(columns=["protein"])
F["desc"] = F.prey.map(lambda x: hdr2desc.get(x, "")[:50])

cols = ["prey","prey_len","iptm","pct_vs_trackB","ligand_iptm","complex_plddt",
        "idf","min_rmsd","plddt","desc"]
F = F.sort_values("pct_vs_trackB", ascending=False)
print()
print(F[cols].to_string(index=False))

F.to_csv(DIR["table"]/"trackB_focus_metal_motif.csv", index=False)
print("\n저장:", DIR["table"]/"trackB_focus_metal_motif.csv")
print("\n※ pct_vs_trackB 는 같은 길이 구간의 Track B 단백질 중 몇 %보다 높은지다.")
print("   90 이상 = 구조 기하(모티프)와 co-folding(ipTM) 두 축이 함께 가리키는 것.")
print("   ipTM 절대값만 보고 순위를 매기면 길이에 속는다 (CELL 28a 참조).")
```

---

## CELL 28a-6 — Part 7d-6. 막단백질 교란 점검

```python
# =============================================================================
# CELL 28a-6 | 상위권이 막단백질로 채워졌는가
#   구조 예측기는 막단백질의 소수성 표면에 다른 사슬을 잘 붙인다. 실제로는 지질이
#   놓일 자리인데 파트너가 대신 들어가는 것이다. ChCODH2 는 세포질 효소라
#   막횡단 나선 표면에 도킹하는 결과는 생물학적으로 성립하지 않는다.
#   실측(27개): 상위 7개 중 넷이 수송체·막결합 효소였다.
#   ※ 주석 키워드로 거르는 것은 어림짐작이다. 최종 후보는 TMHMM/DeepTMHMM 같은
#     막횡단 예측으로 확인할 것. 여기서는 편향의 존재 여부만 본다.
# =============================================================================
need(((DIR["table"]/"trackB_focus_metal_motif.csv").exists(),
      "CELL 28a-5 를 먼저 돌릴 것"))
F = pd.read_csv(DIR["table"]/"trackB_focus_metal_motif.csv")

MEMBRANE_KW = ("transport|export|import|permease|channel|efflux|porin|"
               "synthase|phosphatase|membrane|secretion|pilin|flagell")
F["막추정"] = F.desc.fillna("").str.contains(MEMBRANE_KW, case=False, regex=True)

print("=== 막추정 여부별 백분위 ===")
print(F.groupby("막추정").pct_vs_trackB
       .agg(n="count", 평균="mean", 중앙값="median", 최대="max").round(1).to_string())
gap = (F[F.막추정].pct_vs_trackB.mean() - F[~F.막추정].pct_vs_trackB.mean()) if F.막추정.any() else 0
print(f"\n평균 차이: {gap:+.1f}%p", "-> 막 편향이 있다" if gap > 15 else "-> 뚜렷한 편향은 없다")

cols = ["prey","prey_len","iptm","pct_vs_trackB","idf","min_rmsd","plddt","desc"]
print("\n=== 막추정 제외 상위 10 ===")
print(F[~F.막추정].nlargest(10, "pct_vs_trackB")[cols].to_string(index=False))

F.to_csv(DIR["table"]/"trackB_focus_metal_motif.csv", index=False)
print("\n막추정 컬럼을 추가해 다시 저장했다.")
print("\n※ 해석 시 유의: 지금까지 확인된 교란 세 가지")
print("   1) 길이   — 짧은 prey 일수록 ipTM 이 부풀려진다 (CELL 28a)")
print("   2) 구조품질 — pLDDT 낮은 단편은 folddisco idf 와 ipTM 을 동시에 부풀린다 (CELL 28a-3)")
print("   3) 막소수성 — 지질이 놓일 자리에 파트너가 대신 들어간다 (이 셀)")
print("   RF2-PPI(공진화)는 셋 모두와 무관한 축이므로 Track A 결과가 판정의 중심이다.")
```

---

## CELL 28a-7 — Part 7d-7. Boltz 가 Ni 을 실제로 어디에 놓았나

```python
# =============================================================================
# CELL 28a-7 | ligand_iptm 이 무엇을 재고 있었는지 좌표로 확인한다
#   Ni 을 리간드로 넣은 이유는 "Ni 이 두 사슬 경계면에 놓이는가"를 묻기 위해서였다.
#   그게 실제로 일어났는지는 점수가 아니라 좌표를 봐야 안다.
#   실측(out+out_focus 348 구조): 계면 0 / ChCODH2 쪽 335 / 후보 쪽 13.
#   후보 쪽으로 간 것: Tgt · YeiR · GmhB · PhnP · HslO · HisB (focus 27개 중 6개)
#   -> 전달 복합체는 한 번도 모델링되지 않았다. ligand_iptm 은 "Ni 이 계면에 있나"가
#      아니라 "ChCODH2 안 어디에 박혔나"를 재고 있었고, 그래서 길이와 rho=-0.81 로
#      묶였다 (CELL 28a). 이 지표는 랭킹 축에서 제외한다.
#   배위 거리는 타당했다: A측 평균 2.24A, B측 1.90~2.29A (Ni-S 2.2 / Ni-O 2.0 근처).
#   모델이 아무 데나 얹은 게 아니라 진짜 금속 자리를 고른 것이다.
# =============================================================================
import math
from collections import Counter
CUT_NI = 3.5          # 배위 거리 상한 (A)

def cif_atoms(cif):
    # _atom_site 루프를 헤더 순서대로 읽는다 (컬럼 순서가 파일마다 다를 수 있다)
    cols, rows, inloop = [], [], False
    for ln in open(cif, errors="ignore"):
        if ln.startswith("_atom_site."):
            cols.append(ln.strip().split(".")[1]); inloop = True; continue
        if inloop:
            if ln.startswith(("#", "loop_", "_")) or not ln.strip():
                if cols and rows: break
                continue
            f = ln.split()
            if len(f) == len(cols): rows.append(f)
    i = {c: k for k, c in enumerate(cols)}
    out = []
    for f in rows:
        try:
            out.append((f[i["label_comp_id"]], f[i["label_asym_id"]], f[i["label_seq_id"]],
                        f[i["label_atom_id"]],
                        float(f[i["Cartn_x"]]), float(f[i["Cartn_y"]]), float(f[i["Cartn_z"]])))
        except (KeyError, ValueError):
            pass
    return out

rows = []
# out_cmp(MG1655·Y19) 와 out_rep(반복 측정) 이 나중에 생겼다. 고정 목록으로 들면
# 새로 만든 예측이 조용히 빠진다 — 실제로 MG1655·Y19 의 Ni 배치가 통째로 비어
# 있었던 원인이다. 있는 out* 를 전부 훑는다.
_SETS = sorted(d.name for d in DIR["boltz"].glob("out*") if d.is_dir())
print("읽을 예측 폴더:", _SETS)
for d in _SETS:
    # diffusion_samples>1 인 폴더(out_rep)도 model_0 한 개만 본다. 배치 재현성은
    # ipTM 이 아니라 좌표로 따로 봐야 하므로 여기서는 대표값으로 충분하다.
    for cif in glob.glob(str(DIR["boltz"]/d/"**"/"*_model_0.cif"), recursive=True):
        name = Path(cif).stem.replace("_model_0", "")
        A = cif_atoms(cif)
        ni = [a for a in A if a[0] == "NI"]
        if not ni:
            rows.append({"set": d, "pair": name, "Ni": "없음"}); continue
        _, _, _, _, x, y, z = ni[0]
        near, mind = Counter(), {}
        for comp, ch, seq, at, ax, ay, az in A:
            if comp == "NI": continue
            dist = math.dist((x, y, z), (ax, ay, az))
            mind[ch] = min(mind.get(ch, 9e9), dist)
            if dist <= CUT_NI: near[(ch, seq)] = 1
        chains = Counter(ch for ch, _ in near)
        rows.append({"set": d, "pair": name, "Ni": "있음",
                     "A_배위": chains.get("A", 0), "B_배위": chains.get("B", 0),
                     "A_최단": round(mind.get("A", float("nan")), 2),
                     "B_최단": round(mind.get("B", float("nan")), 2),
                     "계면": chains.get("A", 0) > 0 and chains.get("B", 0) > 0})

NIP = pd.DataFrame(rows)
def _where(r):
    if r.get("Ni") == "없음": return "Ni 없음"
    if r.get("계면"):          return "계면 (양쪽 사슬에 배위)"
    if r.get("A_배위", 0) > 0: return "ChCODH2 쪽만"
    if r.get("B_배위", 0) > 0: return "후보 쪽만"
    return f"어디에도 안 붙음 (>{CUT_NI}A)"
NIP["위치"] = NIP.apply(_where, axis=1)

print(f"구조 {len(NIP)}개\n")
print(NIP.위치.value_counts().to_string())
print("\n=== ChCODH2 쪽 Ni 최단거리 분포 ===")
print(NIP[NIP.A_배위 > 0].A_최단.describe().round(2).to_string())
NIP.to_csv(DIR["table"]/"boltz_ni_placement.csv", index=False)
print("\n저장:", DIR["table"]/"boltz_ni_placement.csv")
```

---

## CELL 28a-8 — Part 7d-8. Ni 을 후보 쪽에 뺏어온 단백질 + 구조 묶기

```python
# =============================================================================
# CELL 28a-8 | 636잔기 CODH 의 C-cluster 자리를 제치고 Ni 을 가져간 단백질들
#   insertase 라면 금속을 든 쪽이 공여자이므로, 이 구도가 의도했던 "전달"에 가깝다.
#   실측: Track B 321개 중 7개(2.2%) vs metal 모티프 27개 중 6개(22.2%) — 10배.
#   -> Folddisco 가 CooC1 의 Cys 쌍 기하로 고른 것들이 독립 모델(Boltz)에서도
#      강한 금속 자리로 확인됐다. 모티프 검색의 "정밀도"는 신뢰할 만하다
#      (MinD 를 놓친 데서 보듯 민감도는 별개 문제다 — CELL 28a-6 참조).
#
#   금속의 정체는 미해결로 남는다. 모델은 전자 구조를 계산하지 않아 Ni(2+)와
#   Zn(2+)를 구분하지 못한다. 다만 그것이 후보를 무효화하지는 않는다:
#     - DB 의 "zinc-binding" 주석은 대개 처음 구조를 푼 조건에서 붙은 이름이지
#       "Zn 전용"이라는 주장이 아니다.
#     - Irving-Williams 계열상 Ni(2+)는 대부분의 배위자 조합에서 Zn(2+)과 같거나
#       더 강하게 결합한다. 세포 안에서 Zn 이 그 자리를 차지하는 것은 자리의
#       선호가 아니라 Zn 완충·전달 체계의 배분 결과다.
#     - 알려진 Ni 관련 단백질도 그렇다. HypA 는 구조적 Zn 자리와 Ni 자리를 함께
#       갖고, SlyD 는 Zn 결합 도메인을 가진 채 HypB 에 Ni 을 넘긴다.
#   기능 정보가 없는 discovery 단계에서는 어노테이션이 아니라 "금속을 잡을 줄
#   아는가"로 좁히는 것이 맞다. 어노테이션으로 거르면 이미 알려진 것만 남는다.
#   따라서 아래 목록은 배위 수와 거리로 줄을 세운다. 금속 정체성은 실험
#   (정제 단백질 금속 재구성·ICP-MS)이 답할 몫이다.
# =============================================================================
need((have("NIP"), "CELL 28a-7 을 먼저 돌릴 것"))
NIP["prey"] = NIP.pair.str.split("__").str[-1]
_B = pd.read_csv(DIR["table"]/"trackB_boltz2_ranked.csv").set_index("prey")
_M = set(pd.read_csv(DIR["table"]/"folddisco_metal_motif.csv").protein)

sub = NIP[NIP.B_배위 > 0].copy()
sub["iptm"]  = sub.prey.map(lambda x: _B.loc[x, "iptm"] if x in _B.index else float("nan"))
sub["metal"] = sub.prey.isin(_M)
sub["desc"]  = sub.prey.map(lambda x: hdr2desc.get(x, "")[:55])
print("=== Ni 이 후보 단백질 쪽에 놓인 쌍 ===")
print(sub.sort_values("B_배위", ascending=False)
         [["set","prey","B_배위","B_최단","iptm","metal","desc"]].to_string(index=False))

for tag, n in [("Track B", (NIP.set == "out").sum()), ("focus", (NIP.set == "out_focus").sum())]:
    k = ((NIP.B_배위 > 0) & (NIP.set == ("out" if tag == "Track B" else "out_focus"))).sum()
    print(f"\n{tag}: {k}/{n} = {100*k/max(n,1):.1f}% 가 Ni 을 후보 쪽에 가져감")

# ---- 구조 파일 묶기 + 배위 잔기 목록 + PyMOL 스크립트 ----
NIV = DIR["boltz"].parent/"ni_view"; NIV.mkdir(exist_ok=True)
for f in NIV.iterdir(): f.unlink()
cifs = {Path(c).stem.replace("_model_0",""): c
        for d in ["out","out_focus"]
        for c in glob.glob(str(DIR["boltz"]/d/"**"/"*_model_0.cif"), recursive=True)}

pml, report = ["bg_color white", "set cartoon_transparency, 0.4", ""], []
for name in sub.pair:
    cif = cifs.get(name)
    if not cif: continue
    prey = name.split("__")[-1]
    tag  = f"{prey}_B측"
    dst  = NIV/f"{tag}.cif"
    shutil.copy(cif, dst)
    A  = cif_atoms(cif)
    _, _, _, _, x, y, z = [a for a in A if a[0] == "NI"][0]
    near = {}
    for comp, ch, seq, at, ax, ay, az in A:
        if comp == "NI": continue
        dist = math.dist((x, y, z), (ax, ay, az))
        k = (ch, seq, comp)
        if dist <= 4.0 and (k not in near or dist < near[k][0]):
            near[k] = (dist, at)
    report.append(f"\n### {tag}\n    {hdr2desc.get(prey,'')[:70]}")
    for dist, ch, seq, comp, at in sorted((d_, c_, s_, m_, a_)
                                          for (c_, s_, m_), (d_, a_) in near.items()):
        report.append(f"    {dist:4.2f} A  chain{ch}({'ChCODH2' if ch=='A' else '후보'})  {comp}{seq}  {at}")
    pml += [f"load {dst.name}, {tag}", f"hide everything, {tag}", f"show cartoon, {tag}",
            f"color grey70, {tag} and chain A", f"color skyblue, {tag} and chain B",
            f"show spheres, {tag} and resn NI", f"color green, {tag} and resn NI",
            f"set sphere_scale, 0.5, {tag} and resn NI",
            f"show sticks, byres ({tag} and polymer within 4 of ({tag} and resn NI))",
            f"color orange, byres ({tag} and polymer within 4 of ({tag} and resn NI))", ""]

(NIV/"view.pml").write_text("\n".join(pml + ["disable all", "orient"]))
(NIV/"coordination.txt").write_text("\n".join(report))
print("\n".join(report))
print(f"\n저장: {NIV}")
print("  로컬에서:  scp -r unsit@100.81.7.34:" + str(NIV) + " ~/Desktop/")
print("             cd ~/Desktop/ni_view && pymol view.pml")
print("  ※ 배위 잔기가 Cys 여러 개면 진짜 금속 자리, Asp/Glu 산소 한둘에 3.5A 언저리면")
print("     모델이 표면에 얹어놓은 것이다. coordination.txt 로 뷰어 없이 판단 가능하다.")
```

---

## CELL 28a-8b — Ni 배위 자리를 도너 종류로 분류 (Cys4 > Cys3His > His/Asp)

```python
# =============================================================================
# CELL 28a-8b | "Ni 가 후보 쪽에 붙었다"를 화학으로 등급 매긴다
#   앞 셀은 '붙었나 안 붙었나'만 봤다. 그런데 같은 '붙었다'도 자리의 성격이 다르다.
#     Cys4 (티올레이트 4개)  Ni(II) 같은 연한 금속에 가장 맞는 자리. 전달·저장형
#     Cys3His / Cys2His2     구조적 아연 자리의 전형. Ni 도 들어갈 수 있다
#     His/Asp/Glu 중심       가수분해효소의 촉매 금속 자리. Zn/Mg/Mn 용이고
#                            배위 거리도 길어지는 경향. 전달 단백질일 가능성 낮음
#   Irving-Williams 로는 Ni(II) 가 Zn(II) 자리에 들어갈 수 있으니 아연 어노테이션을
#   이유로 빼지는 않는다. 다만 '어떤 자리인가'는 우선순위의 근거가 된다.
#
#   ★ Track B(out) 에서 온 것도 같이 본다. 그쪽은 folddisco 를 안 거친 후보라
#     CooC1 의 A112/A114 간격과 다른 금속 자리를 가진 단백질이 들어올 수 있다.
#     실제로 Cys4 자리를 가진 hypothetical 이 거기서 나왔다.
# =============================================================================
need((have("NIP"), "CELL 28a-7 을 먼저 돌릴 것"))
import math
from collections import Counter

# folddisco 금속 목록과의 대조. prey 는 'Y19-AKE60400.1' 처럼 균주 접두사가 붙어
# 있고, folddisco_metal_motif.csv 는 crosswalk 이 붙은 것만 담아 MG1655 가 0건이다.
# 둘 다 무시하면 out_cmp 행이 전부 False 로 찍힌다 — 실제로 그랬다.
# 균주별 (strain, key) 집합으로 비교한다.
_mp = DIR["table"]/"motif_metal_3strain.csv"
_MSET = set()
if _mp.exists():
    _m3 = pd.read_csv(_mp)
    _MSET = set(zip(_m3.strain.astype(str), _m3.key.astype(str)))
else:
    print("⚠ motif_metal_3strain.csv 없음 (CELL 34) — folddisco 대조를 건너뛴다")
def _in_motif(prey):
    t = str(prey)
    st, key = (t.split("-", 1) if t.split("-", 1)[0] in ("Y19", "MG1655", "BL21")
               else ("BL21", t))
    return (st, key) in _MSET
CUT_NI = globals().get("CUT_NI", 3.5)
TIGHT  = 2.6      # Ni-S 2.2 / Ni-N 2.1 근처. 이보다 멀면 배위로 세지 않는다

rows = []
for _, r in NIP[NIP.get("B_배위", 0) > 0].iterrows():
    cif = next(iter(glob.glob(str(DIR["boltz"]/r["set"]/"**"/f"{r['pair']}_model_0.cif"),
                              recursive=True)), None)
    if not cif: continue
    A = cif_atoms(cif)
    ni = [a for a in A if a[0] == "NI"]
    if not ni: continue
    _, _, _, _, x, y, z = ni[0]
    donors = []
    for comp, ch, seq, at, ax, ay, az in A:
        if comp == "NI" or ch != "B": continue
        d = math.dist((x, y, z), (ax, ay, az))
        if d <= TIGHT and at not in ("N", "O", "C", "CA", "CB"):
            donors.append((round(d, 2), comp, seq, at))
        elif d <= TIGHT and comp in ("ASP", "GLU", "HIS", "CYS"):
            donors.append((round(d, 2), comp, seq, at))
    donors = sorted({(c, s): (d, c, s, a) for d, c, s, a in donors}.values())
    cnt = Counter(c for _, c, _, _ in donors)
    ncys, nhis = cnt.get("CYS", 0), cnt.get("HIS", 0)
    nacid = cnt.get("ASP", 0) + cnt.get("GLU", 0)
    if   ncys >= 4:            grade, why = "A", "Cys4 — 연한 금속 전달·저장형 자리"
    elif ncys == 3:            grade, why = "B", "Cys3 + 보조 — 구조적 금속 자리"
    elif ncys == 2 and nhis:   grade, why = "B", "Cys2His — 구조적 금속 자리"
    elif ncys or nhis >= 3:    grade, why = "C", "His 중심 — 촉매 금속 자리 가능성"
    else:                      grade, why = "D", "산소 공여체 위주 — 표면 부착 의심"
    rows.append({"grade": grade, "set": r["set"], "prey": r["prey"],
                 "n_donor": len(donors), "Cys": ncys, "His": nhis, "Asp/Glu": nacid,
                 "min_dist": donors[0][0] if donors else None,
                 "strain": (str(r["prey"]).split("-", 1)[0]
                            if str(r["prey"]).split("-", 1)[0] in ("Y19", "MG1655") else "BL21"),
                 "folddisco_metal": _in_motif(r["prey"]),
                 "donors": " ".join(f"{c}{s}" for _, c, s, _ in donors),
                 "why": why, "desc": str(r.get("desc", ""))[:55]})

G = pd.DataFrame(rows).sort_values(["grade", "Cys", "min_dist"],
                                   ascending=[True, False, True])
if not _MSET:
    G["folddisco_metal"] = None
G.to_csv(DIR["table"]/"ni_site_grade.csv", index=False, encoding="utf-8-sig")
pd.set_option("display.max_colwidth", 60); pd.set_option("display.width", 220)
print(G.to_string(index=False))
print("\n=== 등급별 ===")
print(G.groupby("grade").size().to_string())
print("\n=== folddisco 금속 모티프에 없던 것 (Track B 에서만 나온 후보) ===")
_new = G[~G.folddisco_metal]
print(_new[["grade", "prey", "Cys", "His", "donors", "desc"]].to_string(index=False)
      if len(_new) else "  없음")
print("\n※ folddisco 질의는 CooC1 의 A112/A114 간격 하나만 찾는다. 간격이 다른 금속")
print("  자리는 원리상 안 잡힌다. 위 목록이 그 사각지대에 있던 후보들이다.")
print(f"\n저장: {DIR['table']/'ni_site_grade.csv'}")
```

---

## CELL 28a-9 — Part 7d-9. MG1655 · Y19 metal 모티프 단백질 co-folding (대조군)

```python
# =============================================================================
# CELL 28a-9 | 다른 두 균주의 금속 모티프 단백질도 같은 방식으로 co-folding 한다
#   지금까지 계산은 BL21 만 했다. 그런데 실험의 전제가 "BL21 은 되고 MG1655 는
#   안 된다" 이므로, 대조군 없이 BL21 후보만 보는 것은 비교가 성립하지 않는다.
#     BL21   : folddisco 29 / Boltz 321+27 / RF2-PPI 진행 중
#     Y19    : folddisco 30 / Boltz 없음
#     MG1655 : folddisco 21 / Boltz 없음
#   -> 여기서 Y19 30 + MG1655 21 = 51쌍을 채워 세 균주를 같은 축에 올린다.
#
#   서열 출처가 균주마다 다르다:
#     Y19    crosswalk 로 GenBank ID 가 붙어 있으므로 y19_db_match.faa 에서 가져온다
#     MG1655 공식 배포 인덱스라 crosswalk 이 없다. tid 의 UniProt accession 으로
#            UniProt REST 에서 직접 받는다 (실측: 21건 모두 조회됨)
#   ※ AlphaFold 버전 차이(BL21 v6 / MG1655 v4)는 folddisco 선택 단계에만 영향을
#     준다. Boltz 는 구조가 아니라 서열에서 접으므로 이 비교 자체는 공정하다.
# =============================================================================
need((have("BAIT_KEY", "BAIT_ALL"), "CELL 07 을 먼저 돌릴 것"))
import yaml as _y, urllib.request

CMP_GPU      = 0          # RF2-PPI 가 GPU 1 을 쓰는 동안 여기는 0
BOLTZ_MAXMSA = globals().get("BOLTZ_MAXMSA", 2048)
MAX_TOTAL    = 1830       # bait+prey 토큰 상한 (24GB 에서 안전한 선)

def _read_fasta(path):
    out, nm, buf = {}, None, []
    for l in open(path, errors="ignore"):
        if l.startswith(">"):
            if nm: out[nm] = "".join(buf)
            nm, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if nm: out[nm] = "".join(buf)
    return out

fdd = pd.read_csv(DIR["table"]/"folddisco_metal_motif_detail.csv")
targets = {}      # (strain, id) -> seq

# ---- Y19: crosswalk 으로 GenBank ID 가 붙어 있다 ----
y19_seq = _read_fasta(ASSET["faa_y19"])
for pid in fdd[fdd.strain == "Y19"].protein.dropna().unique():
    if pid in y19_seq:
        targets[("Y19", pid)] = y19_seq[pid]
print(f"Y19  {len([k for k in targets if k[0]=='Y19'])}개")

# ---- MG1655: tid 의 UniProt accession 으로 REST 조회 ----
mg_acc = sorted({m.group(1) for t in fdd[fdd.strain == "MG1655"].tid_stem.dropna()
                 for m in [re.search(r"AF-([A-Z0-9]+)-F", str(t))] if m})
print(f"MG1655 accession {len(mg_acc)}개 조회 중...")
try:
    url = ("https://rest.uniprot.org/uniprotkb/stream?query="
           + "%20OR%20".join(f"accession:{a}" for a in mg_acc) + "&format=fasta")
    with urllib.request.urlopen(url, timeout=120) as r:
        fa = r.read().decode()
    nm, buf = None, []
    for l in fa.splitlines():
        if l.startswith(">"):
            if nm: targets[("MG1655", nm)] = "".join(buf)
            nm, buf = l.split("|")[1], []
        else: buf.append(l.strip())
    if nm: targets[("MG1655", nm)] = "".join(buf)
except Exception as e:
    print(f"  ⚠ UniProt 조회 실패: {e}")
print(f"MG1655 {len([k for k in targets if k[0]=='MG1655'])}개")

# ---- 입력 YAML ----
bseq = BAIT_ALL[BAIT_KEY]
cin = DIR["boltz"]/"inputs_cmp"; cin.mkdir(parents=True, exist_ok=True)
for f in cin.glob("*.yaml"): f.unlink()
n, skip = 0, []
for (strain, pid), sq in sorted(targets.items()):
    if len(bseq) + len(sq) > MAX_TOTAL:
        skip.append((strain, pid, len(sq))); continue
    (cin/f"{BAIT_KEY}__{strain}-{pid}.yaml").write_text(_y.safe_dump({
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": bseq}},
            {"protein": {"id": "B", "sequence": sq}},
            {"ligand":  {"id": "C", "ccd": "NI"}},
        ]}, sort_keys=False))
    n += 1
print(f"\n입력 {n}개 -> {cin}")
if skip: print(f"길이 초과 제외 {len(skip)}개: {skip}")
print(f"예상 소요: {n} x 약 60초 = 약 {n/60:.1f}시간 (GPU {CMP_GPU})")

script = f"""
cd "{DIR['boltz']}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
CUDA_VISIBLE_DEVICES={CMP_GPU} boltz predict inputs_cmp \\
  --out_dir out_cmp \\
  --use_msa_server \\
  --max_msa_seqs {BOLTZ_MAXMSA} \\
  --recycling_steps 3 \\
  --diffusion_samples 1 \\
  --output_format mmcif \\
  --num_workers 2
echo DONE_boltz_cmp
"""
_busy = subprocess.run("pgrep -f '[b]oltz predict'", shell=True,
                       capture_output=True, text=True).stdout.split()
if _busy:
    print(f"⚠ boltz 가 이미 돌고 있다 (pid {' '.join(_busy)}). 새로 띄우지 않았다.")
elif not already(done("part7_boltz_cmp", DIR["boltz"]/"out_cmp"),
                 "MG1655·Y19 비교 예측", f"rm -rf {DIR['boltz']}/out_cmp"):
    sh_bg("part7_boltz_cmp", script, env=CONDA_ENV_BOLTZ)
```

---

## CELL 28a-10 — Part 7d-10. 세 균주 비교 판정

```python
# =============================================================================
# CELL 28a-10 | BL21 · Y19 · MG1655 의 금속 모티프 단백질을 같은 축에서 비교
#   묻는 것: 금속 자리를 가진 단백질 중 BL21 것이 ChCODH2 와 더 잘 붙는가.
#   실험 관찰("BL21 lysate 만 활성을 회복시킨다")을 계산으로 옮긴 형태다.
#   길이 구간별로 나눠 본다 — ipTM 은 짧은 prey 를 부풀리므로 (CELL 28a) 전체를
#   한 줄로 세우면 길이 분포 차이가 균주 차이로 오독된다.
# =============================================================================
need(((DIR["boltz"]/"out_cmp").exists(), "CELL 28a-9 를 먼저 돌릴 것"))
BINS = [0, 50, 100, 200, 400, 10**5]

def collect(dirname, strain_from):
    rows = []
    for f in glob.glob(str(DIR["boltz"]/dirname/"**"/"confidence_*.json"), recursive=True):
        c = json.load(open(f))
        nm = re.sub(r"^confidence_|_model_\d+$", "", Path(f).stem).split("__")[-1]
        st, pid = strain_from(nm)
        rows.append({"strain": st, "protein": pid, "iptm": c.get("iptm"),
                     "ligand_iptm": c.get("ligand_iptm"),
                     "complex_plddt": c.get("complex_plddt")})
    return pd.DataFrame(rows)

C = collect("out_cmp", lambda n: tuple(n.split("-", 1)) if "-" in n else ("?", n))
F = collect("out_focus", lambda n: ("BL21", n))
ALL = pd.concat([F, C], ignore_index=True).drop_duplicates(["strain","protein"])
print("예측 건수:"); print(ALL.strain.value_counts().to_string())

# 길이: BL21/Y19 는 faa, MG1655 는 folddisco 의 nres 로 대신한다
fdd = pd.read_csv(DIR["table"]/"folddisco_metal_motif_detail.csv")
nres = {}
for r in fdd.itertuples(index=False):
    m = re.search(r"AF-([A-Z0-9]+)-F", str(r.tid_stem))
    if pd.notna(r.protein): nres[str(r.protein)] = r.nres
    if m: nres[m.group(1)] = r.nres
ALL["prey_len"] = ALL.protein.map(nres)
ALL["len_bin"]  = pd.cut(ALL.prey_len, BINS)

print("\n=== 균주별 ipTM ===")
print(ALL.groupby("strain").iptm.agg(n="count", mean="mean", median="median",
                                     max="max").round(3).to_string())
print("\n=== 길이 구간 x 균주 (평균 ipTM) ===")
print(ALL.pivot_table(index="len_bin", columns="strain", values="iptm",
                      aggfunc="mean", observed=True).round(3).to_string())
print("\n=== 균주별 상위 8 ===")
for st in ["BL21", "Y19", "MG1655"]:
    sub = ALL[ALL.strain == st].nlargest(8, "iptm")
    if len(sub):
        print(f"\n[{st}]")
        print(sub[["protein","prey_len","iptm","ligand_iptm","complex_plddt"]].to_string(index=False))

ALL.to_csv(DIR["table"]/"metal_motif_3strain_boltz.csv", index=False, encoding="utf-8-sig")
print("\n저장:", DIR["table"]/"metal_motif_3strain_boltz.csv")
print("\n※ 해석 주의")
print("   - 균주별 ipTM 평균이 비슷하면 '금속 자리 보유'만으로는 균주를 못 가른다는 뜻이다.")
print("   - BL21 이 뚜렷이 높아야 실험 관찰과 방향이 맞는다.")
print("   - 표본이 작다 (BL21 27 / Y19 30 / MG1655 21). 차이를 단정하지 말 것.")
print("   - diffusion_samples=1 이라 각 값은 단일 표본이다. 재현성 정보가 없다.")
```

---

## CELL 28a-11 — Part 7d-11. 같은 단백질끼리 짝지어 비교 (직교체 페어링)

```python
# =============================================================================
# CELL 28a-11 | 균주 평균이 아니라 '같은 단백질의 균주 간 차이'를 본다
#   CELL 28a-10 은 균주별 평균을 냈다. 그런데 세 균주의 금속모티프 단백질 목록은
#   구성이 다르다 (BL21 27 / Y19 30 / MG1655 20). 구성이 다른 집합의 평균 차이는
#   균주 차이인지 목록 차이인지 구분되지 않는다.
#   여기서는 서열로 직교체를 짝지어 두 가지를 얻는다:
#     (a) 서열이 사실상 동일한 짝의 |Δ ipTM| = 이 측정의 잡음 바닥
#         같은 서열이면 Boltz 입력이 같다. 남는 차이는 전부 diffusion 표본 오차다.
#         diffusion_samples=1 로 돌렸으므로 반복이 없다. 동일 서열 짝이 현재
#         가진 유일한 잡음 추정치다.
#     (b) 그 잡음 바닥을 넘는 균주 차이가 실제로 있는가
#   (a) 를 모르면 (b) 를 말할 수 없다. 0.68 대 0.40 이 의미 있는 차이인지는
#   잡음이 0.05 인지 0.30 인지에 달려 있다.
# =============================================================================
need(((DIR["table"]/"metal_motif_3strain_boltz.csv").exists(), "CELL 28a-10 을 먼저 돌릴 것"))
import difflib, urllib.request

A3 = pd.read_csv(DIR["table"]/"metal_motif_3strain_boltz.csv")

def _read_fasta(path):
    out, nm, buf = {}, None, []
    for l in open(path, errors="ignore"):
        if l.startswith(">"):
            if nm: out[nm] = "".join(buf)
            nm, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if nm: out[nm] = "".join(buf)
    return out

# MG1655 서열은 CELL 28a-9 와 같은 출처(UniProt)에서 받아 캐시해 둔다.
_cache = DIR["table"]/"mg1655_metal_motif.faa"
if not _cache.exists():
    _acc = sorted(A3[A3.strain == "MG1655"].protein.astype(str).unique())
    _url = ("https://rest.uniprot.org/uniprotkb/stream?query="
            + "%20OR%20".join(f"accession:{a}" for a in _acc) + "&format=fasta")
    with urllib.request.urlopen(_url, timeout=120) as r:
        _fa = r.read().decode()
    _cache.write_text("\n".join((">" + l.split("|")[1]) if l.startswith(">") else l
                                for l in _fa.splitlines()))
    print(f"MG1655 서열 {_fa.count('>')}개 캐시 -> {_cache}")

POOL = {"BL21":   _read_fasta(ASSET["faa_bl21"]),
        "Y19":    _read_fasta(ASSET["faa_y19"]),
        "MG1655": _read_fasta(_cache)}
have_seq = {st: {p: POOL[st][p] for p in A3[A3.strain == st].protein.astype(str)
                 if p in POOL[st]} for st in POOL}
for st in ["BL21", "Y19", "MG1655"]:
    print(f"{st}: 서열 확보 {len(have_seq[st])} / 예측 {int((A3.strain == st).sum())}")

def best_match(seq, pool, cutoff=0.80):
    """길이로 1차 거르고 quick_ratio(상한)로 2차 거른 뒤 정확히 잰다."""
    best, bid = 0.0, None
    for pid, s in pool.items():
        if abs(len(s) - len(seq)) / max(len(s), len(seq)) > 0.25:
            continue
        sm = difflib.SequenceMatcher(None, seq, s)
        if sm.quick_ratio() < cutoff:      # quick_ratio 는 ratio 의 상한이다
            continue
        r = sm.ratio()
        if r > best: best, bid = r, pid
    return bid, best

IPTM = {(r.strain, str(r.protein)): r.iptm for r in A3.itertuples(index=False)}
DESC = {}
if "hdr2desc" in dir():
    DESC = hdr2desc

rows = []
for pid, sq in have_seq["BL21"].items():
    row = {"BL21": pid, "len": len(sq), "BL21_iptm": IPTM.get(("BL21", pid)),
           "desc": str(DESC.get(pid, ""))[:50]}
    for st in ["MG1655", "Y19"]:
        m, r = best_match(sq, have_seq[st])
        row[st] = m
        row[f"{st}_id"] = round(r, 3) if m else None
        row[f"{st}_iptm"] = IPTM.get((st, m)) if m else None
    rows.append(row)
P = pd.DataFrame(rows).sort_values("BL21_iptm", ascending=False)
P.to_csv(DIR["table"]/"metal_motif_ortholog_paired.csv", index=False, encoding="utf-8-sig")

# ---- (a) 잡음 바닥 ----
NOISE = None
_same = P[(P.MG1655_id >= 0.98) & P.MG1655_iptm.notna()].copy()
_same["dIPTM"] = (_same.BL21_iptm - _same.MG1655_iptm).abs()
print(f"\n=== (a) 서열 98% 이상 일치하는 BL21-MG1655 짝: {len(_same)}개 ===")
if len(_same):
    print(_same[["BL21", "MG1655", "MG1655_id", "BL21_iptm", "MG1655_iptm",
                 "dIPTM", "desc"]].round(3).to_string(index=False))
    NOISE = float(_same.dIPTM.median())
    print(f"\n|Δ ipTM| 중앙값 {NOISE:.3f} / 평균 {_same.dIPTM.mean():.3f} / "
          f"최대 {_same.dIPTM.max():.3f}")
    print("   ← 입력이 사실상 같은데도 이만큼 벌어진다. 이보다 작은 균주 차이는")
    print("     해석 대상이 아니다. (Boltz 는 diffusion 이라 비결정적이다)")
else:
    print("동일 서열 짝이 없다. 잡음 바닥을 못 구했으므로 (b) 는 판정 보류.")

# ---- (b) 짝이 없는 BL21 단백질 ----
_only = P[P.MG1655.isna()]
print(f"\n=== (b) MG1655 목록에 짝이 없는 BL21 단백질: {len(_only)}개 ===")
print("※ '유전자가 없다'가 아니라 'MG1655 folddisco 목록에 안 잡혔다'는 뜻이다.")
print("  (이전 검토에서 이들 다수가 대장균 공통 유전자로 확인됐다 — 검출 차이다)")
if len(_only):
    print(_only[["BL21", "len", "BL21_iptm", "Y19", "Y19_id", "Y19_iptm",
                 "desc"]].round(3).to_string(index=False))

# ---- 판정 ----
print("\n=== 판정 ===")
if NOISE is None:
    print("잡음 바닥 미측정 — CELL 28a-12 로 반복 측정을 돌린 뒤 다시 본다.")
else:
    _gap = A3.groupby("strain").iptm.mean()
    _spread = float(_gap.max() - _gap.min())
    print(f"균주 평균 ipTM 최대-최소 차이 {_spread:.3f}  vs  잡음 바닥 {NOISE:.3f}")
    if _spread <= NOISE:
        print("→ 균주 간 차이가 잡음보다 작다. 세 균주를 구분할 신호가 없다.")
        print("  '금속 모티프 보유 + ChCODH2 와의 ipTM' 만으로는 BL21 을 못 고른다.")
    else:
        print("→ 균주 차이가 잡음보다 크다. 다만 표본이 20~30개라 아직 단정하지 말 것.")
    print("  어느 쪽이든 다음 단계는 CELL 28a-12 (반복 측정 + 양성대조군) 다.")
print("\n저장:", DIR["table"]/"metal_motif_ortholog_paired.csv")
```

---

## CELL 28a-12 — Part 7d-12. Boltz 양성대조군 + 상위 후보 반복 측정 (백그라운드)

```python
# =============================================================================
# CELL 28a-12 | Boltz 쪽에 처음으로 대조군을 물린다 + 오차막대를 만든다
#   Track A 는 대조군이 있었고 그걸 틀렸다 (CELL 24b). 반면 Track B/C 의 Boltz 는
#   대조군 자체가 없었다. 즉 ipTM 0.68 이 이 계에서 무슨 뜻인지 아직 모른다.
#   CooC1-ChCODH2 는 문헌상 실제로 결합하는 쌍이다. 이것을 접어서 '맞는 답이
#   몇 점을 받는가'를 먼저 정한다.
#     대조군 ipTM 높음(>=0.6) -> 눈금을 믿고 그 선에서 후보를 자른다
#     대조군 ipTM 낮음         -> 0.68 짜리 후보도 근거가 못 된다. 눈금이 안 통한다
#   동시에 상위 후보를 diffusion_samples 5 로 다시 접어 표준편차를 얻는다.
#   지금 값은 전부 단일 표본이라 0.68 과 0.40 의 차이가 잡음인지 알 수 없다.
#
#   비용: 쌍당 약 1~2분 (표본 5개라도 trunk 는 한 번만 돈다). 15쌍 = 30분 안쪽.
# =============================================================================
need((have("BAIT_KEY", "BAIT_ALL"), "CELL 07 을 먼저 돌릴 것"),
     (have("have_seq"), "CELL 28a-11 을 먼저 돌릴 것 (서열 사전이 거기서 만들어진다)"))
import yaml as _y

REP_GPU      = 1      # RF2-PPI 가 끝났으면 1 이 비어 있다. nvidia-smi 로 확인하고 고칠 것
REP_SAMPLES  = 5      # diffusion 표본 수 = 오차막대의 근거
TOPN         = 4      # 균주별 상위 몇 개를 다시 접을지
BOLTZ_MAXMSA = globals().get("BOLTZ_MAXMSA", 2048)
MAX_TOTAL    = 1830

bseq = BAIT_ALL[BAIT_KEY]
rin  = DIR["boltz"]/"inputs_rep"; rin.mkdir(parents=True, exist_ok=True)
for f in rin.glob("*.yaml"): f.unlink()

def _yam(name, sq):
    if len(bseq) + len(sq) > MAX_TOTAL:
        print(f"  길이 초과 제외: {name} ({len(sq)}aa)"); return 0
    (rin/f"{name}.yaml").write_text(_y.safe_dump({
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": bseq}},
            {"protein": {"id": "B", "sequence": sq}},
            {"ligand":  {"id": "C", "ccd": "NI"}},
        ]}, sort_keys=False))
    return 1

# ---- 1. 양성대조군: CooC1 / CooT / CooJ 중 서열이 있는 것 전부 ----
n_ctrl = 0
for k, sq in BAIT_ALL.items():
    if k == BAIT_KEY or k.startswith("seg"):
        continue
    if not k.lower().startswith(("cooc", "coot", "cooj")):
        continue
    n_ctrl += _yam(f"CTRL-{k}", sq)
print(f"양성대조군 {n_ctrl}개")
if n_ctrl == 0:
    print("  ⚠ BAIT_ALL 에 Coo 계열 서열이 없다. CELL 05 에 넣고 CELL 07 을 다시 돌릴 것.")
    print("     대조군 없이 아래 후보 반복만 돌리면 눈금 문제는 그대로 남는다.")

# ---- 2. 상위 후보: 균주별 TOPN ----
A3 = pd.read_csv(DIR["table"]/"metal_motif_3strain_boltz.csv")
n_cand = 0
for st in ["BL21", "MG1655", "Y19"]:
    for r in A3[A3.strain == st].nlargest(TOPN, "iptm").itertuples(index=False):
        sq = have_seq.get(st, {}).get(str(r.protein))
        if sq is None:
            print(f"  서열 없음: {st}/{r.protein}"); continue
        n_cand += _yam(f"{st}-{r.protein}", sq)
print(f"상위 후보 {n_cand}개 (균주별 {TOPN})")

_tot = n_ctrl + n_cand
print(f"\n입력 {_tot}개 -> {rin}   표본 {REP_SAMPLES}개/쌍   GPU {REP_GPU}")

script = f"""
cd "{DIR['boltz']}"
export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
CUDA_VISIBLE_DEVICES={REP_GPU} boltz predict inputs_rep \\
  --out_dir out_rep \\
  --use_msa_server \\
  --max_msa_seqs {BOLTZ_MAXMSA} \\
  --recycling_steps 3 \\
  --diffusion_samples {REP_SAMPLES} \\
  --output_format mmcif \\
  --num_workers 2
echo DONE_boltz_rep
"""
_busy = subprocess.run("pgrep -f '[b]oltz predict'", shell=True,
                       capture_output=True, text=True).stdout.split()
if _busy:
    print(f"⚠ boltz 가 이미 돌고 있다 (pid {' '.join(_busy)}). 새로 띄우지 않았다.")
elif _tot == 0:
    print("⚠ 입력이 0개다. 띄우지 않았다.")
elif not already(done("part7_boltz_rep", DIR["boltz"]/"out_rep"),
                 "반복 측정 + 대조군", f"rm -rf {DIR['boltz']}/out_rep"):
    sh_bg("part7_boltz_rep", script, env=CONDA_ENV_BOLTZ)
```

---

## CELL 28a-13 — Part 7d-13. 반복 측정 집계 + 눈금 판정

```python
# =============================================================================
# CELL 28a-13 | 표본 5개로 평균과 표준편차를 내고, 대조군으로 눈금을 정한다
#   읽는 순서:
#     1) 대조군 ipTM  — 이 계에서 '참'이 받는 점수
#     2) 후보의 SD    — 단일 표본 값이 얼마나 흔들리는지
#     3) 후보가 대조군 선을 넘는가, 그 차이가 SD 보다 큰가
# =============================================================================
# 로그가 아니라 결과 파일로 판정한다. done() 은 백그라운드 로그의 완료 마커를
# 보는데, 로그가 지워지거나 다른 경로로 돌린 경우 결과가 있어도 막힌다.
_nrep = len(glob.glob(str(DIR["boltz"]/"out_rep"/"**"/"confidence_*.json"), recursive=True))
need((_nrep > 0, "out_rep 에 결과가 없다 — CELL 28a-12 를 먼저 돌릴 것"))
print(f"out_rep confidence 파일 {_nrep}개")
rows = []
for f in glob.glob(str(DIR["boltz"]/"out_rep"/"**"/"confidence_*.json"), recursive=True):
    c  = json.load(open(f))
    nm = re.sub(r"^confidence_", "", Path(f).stem)
    pair, _, _ = nm.rpartition("_model_")
    rows.append({"pair": (pair or nm).split("__")[-1],
                 "iptm": c.get("iptm"), "ptm": c.get("ptm"),
                 "complex_plddt": c.get("complex_plddt"),
                 "ligand_iptm": c.get("ligand_iptm")})
D = pd.DataFrame(rows)
G = (D.groupby("pair")
       .agg(n=("iptm", "size"), iptm_mean=("iptm", "mean"), iptm_sd=("iptm", "std"),
            iptm_min=("iptm", "min"), iptm_max=("iptm", "max"),
            plddt=("complex_plddt", "mean"))
       .sort_values("iptm_mean", ascending=False).round(3))
G["kind"] = ["대조군" if p.startswith("CTRL-") else p.split("-")[0] for p in G.index]
G.to_csv(DIR["table"]/"boltz_replicate_summary.csv", encoding="utf-8-sig")
print(G.to_string())

_ctrl = G[G.kind == "대조군"]
_cand = G[G.kind != "대조군"]
print("\n=== 1) 눈금 ===")
if len(_ctrl):
    _c = float(_ctrl.iptm_mean.max())
    print(f"양성대조군 최고 ipTM {_c:.3f} (SD {float(_ctrl.iptm_sd.max()):.3f})")
    if _c >= 0.6:
        print("→ 눈금이 통한다. 알려진 참이 높은 점수를 받았다.")
        print(f"   이 선({_c:.2f}) 위의 후보만 다음 단계로 넘긴다.")
    else:
        print("→ ⚠ 알려진 참조차 낮다. ipTM 으로 이 계의 결합을 못 가린다.")
        print("   Track B/C 의 ipTM 순위를 '결합 가능성'으로 읽으면 안 된다.")
        print("   Track A 에 이어 두 번째 대조군 실패다. 계산으로 후보를 좁히는")
        print("   현재 설계 자체를 다시 봐야 한다.")
else:
    print("대조군 결과 없음 — Coo 계열 서열을 넣고 CELL 28a-12 를 다시 돌릴 것.")

print("\n=== 2) 단일 표본의 흔들림 ===")
if len(_cand):
    print(f"후보 SD: 중앙값 {_cand.iptm_sd.median():.3f} / 최대 {_cand.iptm_sd.max():.3f}")
    print(f"후보 (max-min): 중앙값 {(_cand.iptm_max - _cand.iptm_min).median():.3f}")
    print("   ← CELL 28a-10 의 표는 전부 표본 1개다. 이 폭 안의 순위 차이는 무의미하다.")

print("\n=== 3) 대조군 대비 ===")
if len(_ctrl) and len(_cand):
    _c = float(_ctrl.iptm_mean.max())
    _win = _cand[_cand.iptm_mean - _cand.iptm_sd.fillna(0) > _c]
    print(f"대조군보다 (평균-SD) 기준으로도 높은 후보: {len(_win)}개")
    if len(_win):
        print(_win.to_string())
    else:
        print("없다. 지금 후보 중 알려진 참을 넘어서는 것은 없다.")
print("\n저장:", DIR["table"]/"boltz_replicate_summary.csv")
```

---

## CELL 28b — Part 7e. Foldseek-Interface 검색 (선택, Track B 뒤)

```python
# =============================================================================
# CELL 28b | Part 7e. Foldseek-Interface — 예측 복합체의 계면을 PDB 계면과 대조
#   출처: Strom, Cha, Kim et al. "Foldseek-Interface reveals a protein interface
#         universe far from complete" (bioRxiv 2026). 저자에 Cameron Gilchrist 포함.
#   - 이건 PPI 예측 도구가 아니다. 이미 있는 복합체에서 계면 잔기만 뽑아 3Di 로
#     인코딩하고, PDB 의 77,167개 계면 클러스터와 대조해 "이 계면이 알려진
#     결합 방식인가"를 묻는다. ipTM 과 독립적인 증거 축이 하나 생긴다.
#   - 논문의 humanPPI 워크플로(예측 복합체 21,048개 → PDB 계면 검색)와 같은 구조.
#     다만 해석은 반대다. 그쪽은 hit 없는 것을 "novel"로 봤지만, 우리는
#     금속 샤페론 계면(HypA-HypB, UreE-UreG, Atx1-Ccc2 등)에 걸리는 쪽이 신호다.
#   - 전제 2가지: ① 계면 지원 foldseek  ② PDB 계면 대표 DB (Zenodo 22040892)
# =============================================================================
FOLDSEEK_BIN = "foldseek"                       # 계면 명령을 지원하는 버전이어야 한다
PDB_INT_REP  = REFDB.parent / "foldseek_interface" / "PDB_int_rep"   # Zenodo 에서 받아 둘 곳
IFACE_TM_CUT = 0.4                              # 논문 기준: qTM 또는 tTM >= 0.4 이면 hit

ifc_dir = DIR["boltz"] / "interface"; ifc_dir.mkdir(parents=True, exist_ok=True)
ok_bin = "multimersearch" in sh(f"{FOLDSEEK_BIN} --help 2>&1 || true", check=False, quiet=True)
print(f"[{'O' if ok_bin else 'X'}] foldseek 계면 지원")
print(f"[{'O' if PDB_INT_REP.exists() else 'X'}] PDB 계면 대표 DB  {PDB_INT_REP}")

if not ok_bin:
    print("\n-> foldseek 를 계면 지원 버전으로 올린다:")
    print("   conda install -y -c conda-forge -c bioconda 'foldseek>=11'")
    print("   (또는 https://github.com/steineggerlab/foldseek 최신 바이너리)")
elif not PDB_INT_REP.exists():
    print("\n-> PDB 계면 클러스터 자원을 먼저 받는다 (공용이므로 database/ 아래):")
    print(f"   mkdir -p {PDB_INT_REP.parent} && cd {PDB_INT_REP.parent}")
    print("   # https://doi.org/10.5281/zenodo.22040892 에서 PDB 계면 대표 DB 내려받기")
else:
    cif = sorted(glob.glob(str(DIR["boltz"]/"out"/"**"/"*.cif"), recursive=True))
    print(f"\nBoltz-2 예측 구조 {len(cif)}개")
    assert cif, "Boltz-2 출력이 없다. CELL 27 을 먼저 끝낼 것."
    for f in cif:                                # foldseek 은 디렉터리 단위로 읽는다
        shutil.copy(f, ifc_dir / Path(f).name)
    script = f"""
cd "{ifc_dir}"
{FOLDSEEK_BIN} createinterfacedb . trackB_intdb
{FOLDSEEK_BIN} easy-multimersearch trackB_intdb "{PDB_INT_REP}" aln tmp --cov-mode 0
wc -l aln_report
echo DONE_interface
"""
    sh_bg("part7_interface", script, env=CONDA_ENV_RF2)
    print("\n검색이 끝나면 CELL 28c 로 결과를 표로 만든다.")
```

---

## CELL 28c — Part 7f. 계면 일치 표 만들기

```python
# =============================================================================
# CELL 28c | Part 7f. Foldseek-Interface 결과 -> prey 별 계면 일치 표
#   aln_report 컬럼(easy-multimersearch): query target ... qTM tTM ...
#   논문 기준대로 qTM 또는 tTM >= 0.4 를 hit 로 본다.
# =============================================================================
rep = ifc_dir / "aln_report"
if not rep.exists():
    print(f"[없음] {rep} — CELL 28b 가 끝났는지 CELL 15 로 확인할 것.")
    IF = pd.DataFrame(columns=["prey", "iface_tm", "iface_hit"])
else:
    A = pd.read_csv(rep, sep="\t", header=None)
    A.columns = [f"c{i}" for i in range(A.shape[1])]
    # 6, 7번째 컬럼이 qTM / tTM (humanppi 스크립트의 $6, $7)
    A["iface_tm"] = A[["c5", "c6"]].max(axis=1)
    A["prey"] = A.c0.astype(str).apply(lambda x: Path(x).stem.split("__")[-1])
    IF = (A.groupby("prey", as_index=False)
            .agg(iface_tm=("iface_tm", "max"), best_pdb=("c1", "first")))
    IF["iface_hit"] = (IF.iface_tm >= IFACE_TM_CUT).astype(int)
    IF = IF.sort_values("iface_tm", ascending=False)
    IF.to_csv(DIR["table"]/"trackB_interface_match.csv", index=False)
    print(f"계면 검색된 prey {len(IF)}개, PDB 계면과 일치(>= {IFACE_TM_CUT}) {int(IF.iface_hit.sum())}개")
    print(IF.head(25).to_string(index=False))
    print("\n[해석] 일치가 있다 = 이 결합 방식이 PDB 에 전례가 있다.")
    print("       금속 샤페론-표적 계면(HypA-HypB, UreE-UreG, Atx1-Ccc2 등)에 걸리면 강한 방증.")
    print("       일치가 없다 = 새 계면일 수도, Boltz-2 위양성일 수도 있다. 단독으로는 근거가 약하다.")
```

---

## CELL 29 — Part 8a. Folddisco 검색 (Track C)

```python
# =============================================================================
# CELL 29 | Part 8-1. Folddisco 구조 모티프 검색 (Rust 바이너리 — 파이썬 환경 무관)
#   - run01: CooC1 metal-binding motif (residue 112,114). 체인 ID 는 CELL 06 으로 확인.
#   - run02: ATP-binding motif — RESIDUES_ATP 를 지정해야 돌아간다 (인계문서 TODO)
#   ⚠ residue 2개면 관계쌍 1개뿐이라 위양성이 많다. --rmsd 로 조이고 idf 로 정렬.
# =============================================================================
FD = ASSET["folddisco"]
FD_OUT = DIR["folddisco"]; FD_OUT.mkdir(exist_ok=True)
# 3kji.pdb ATOM 레코드로 직접 확인한 값 (체인 A, 잔기번호 1-254, SEQRES 와 어긋남 없음):
#   A112=CYS, A114=CYS            -> metal-binding motif
#   A12=GLY, A13=LYS, A14=THR     -> Walker A (GKGGVGKT 의 G-K-T). CooC 가 ATPase 인 것과 일치
RESIDUES_METAL = "A112,A114"
RESIDUES_ATP   = "A12,A13,A14"   # 잔기 2개면 관계쌍이 1개뿐이라 위양성이 많다 -> 3개로 특이성 확보

# ⚠ MG1655 인덱스는 공식 배포본이라 우리 구조 DB(structures_UP*)에 대응 파일이 없다.
#    -> CELL 30 crosswalk 로 GenBank ID 를 못 붙이므로 CELL 31 에서 제외된다.
#    MG1655 결과는 "BL21 에만 있는 모티프인가"를 눈으로 비교하는 용도로만 쓴다.
IDX = {"BL21": ASSET["fd_idx_bl21"], "Y19": ASSET["fd_idx_y19"], "MG1655": ASSET["fd_idx_mg"]}

def folddisco_query(tag, residues, top=2000, rmsd=1.0):
    if not residues:
        print(f"[스킵] {tag}: residue 미지정"); return {}
    outs = {}
    for strain, idx in IDX.items():
        out = FD_OUT / f"{tag}_{strain}.tsv"
        cmd = (f'"{FD}" query -p "{ASSET["cooc1_pdb"]}" -q {residues} -i "{idx}" '
               f'-t {min(THREADS,8)} --per-structure --header --sort-by idf '
               f'--rmsd {rmsd} --top {top} -o "{out}"')
        print(f"\n>>> {tag} / {strain}")
        sh(cmd, check=False)
        outs[strain] = out
    return outs

FD_METAL = folddisco_query("metal_run01", RESIDUES_METAL)
FD_ATP   = folddisco_query("atp_run02",   RESIDUES_ATP)

# 이미 돌려둔 run01 결과가 있으면 그것도 후보에 넣는다
prev = sorted((ASSET["fd_ws"]/"result"/"run01").glob("*.tsv")) if (ASSET["fd_ws"]/"result"/"run01").exists() else []
print("\n기존 run01 결과:", [p.name for p in prev])
```

---

## CELL 29b — Folddisco --rmsd 민감도 (1.0 / 1.5 / 2.0)

```python
# =============================================================================
# CELL 29b | 컷오프 1.0 A 가 맞는 선인지 재 본다
#   지금 후보 목록은 '--rmsd 1.0' 이라는 선을 정해서 만든 표다. 그런데 그 선이
#   충분한지 검증한 적이 없다. 너무 빡빡하면 진짜 후보가 밖에 있고, 너무 느슨하면
#   아무 Cys 쌍이나 들어온다.
#   기준점으로 MinD 를 쓴다 — CooC1 의 최근접 BL21 구조 이웃(MinD/ParA 폴드)이다.
#   1.0 에서 잡히면 현재 선이 충분한 것이고, 1.5~2.0 에서야 들어오면 지금 목록이
#   너무 좁은 것이다.
#   또 하나: 임계값을 올릴 때 히트가 몇 배로 늘어나는지가 이 모티프의 특이도를 말한다.
#   잔기 2개짜리 질의는 느슨하게 풀면 급격히 무너진다.
# =============================================================================
need((have("ASSET", "sh"), "CELL 01 을 먼저 돌릴 것"))
FD  = ASSET["folddisco"]
IDX = {"BL21": ASSET["fd_idx_bl21"], "Y19": ASSET["fd_idx_y19"], "MG1655": ASSET["fd_idx_mg"]}
RES_METAL = globals().get("RESIDUES_METAL", "A112,A114")
SWEEP = [1.0, 1.5, 2.0]

out_dir = DIR["folddisco"]/"sweep"; out_dir.mkdir(exist_ok=True)
for rm in SWEEP:
    tag = f"metal_r{int(rm*10):02d}"
    for strain, idx in IDX.items():
        out = out_dir/f"{tag}_{strain}.tsv"
        if out.exists() and out.stat().st_size:
            print(f"  [재사용] {out.name}"); continue
        sh(f'"{FD}" query -p "{ASSET["cooc1_pdb"]}" -q {RES_METAL} -i "{idx}" '
           f'-t {min(THREADS,8)} --per-structure --header --sort-by idf '
           f'--rmsd {rm} --top 20000 -o "{out}"', check=False)

# ---------- 집계 ----------
xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(xw) else {}

def load(tag, strain):
    p = out_dir/f"{tag}_{strain}.tsv"
    if not p.exists() or not p.stat().st_size: return pd.DataFrame()
    d = pd.read_csv(p, sep="\t")
    d.columns = [c.strip().lstrip("#") for c in d.columns]
    d["tid_stem"] = d.tid.apply(lambda x: Path(str(x)).stem)
    d["protein"]  = d.tid_stem.map(xwm)
    d["uniprot"]  = d.tid_stem.astype(str).str.extract(r"AF-([A-Z0-9]+)-F", expand=False)
    d["key"] = d.protein.where(d.protein.notna(), d.uniprot).astype(str).str.split(",").str[0]
    return d

TAB, KEYS = [], {}
for rm in SWEEP:
    tag = f"metal_r{int(rm*10):02d}"
    for strain in IDX:
        d = load(tag, strain)
        if not len(d): continue
        KEYS[(rm, strain)] = set(d.key)
        TAB.append({"rmsd": rm, "strain": strain, "hits": len(d),
                    "max_min_rmsd": round(float(d.min_rmsd.max()), 3) if "min_rmsd" in d else None})
S = pd.DataFrame(TAB)
print("\n=== 임계값별 히트 수 ===")
print(S.pivot_table(index="rmsd", columns="strain", values="hits").to_string())
print("\n=== 각 파일의 실제 최대 min_rmsd (컷오프가 하드 필터인지 확인) ===")
print(S.pivot_table(index="rmsd", columns="strain", values="max_min_rmsd").to_string())
print("  임계값과 최대값이 거의 같으면 하드 필터다. 훨씬 작으면 애초에 그 위가 없던 것이다.")

print("\n=== 임계값을 올릴 때 새로 들어오는 것 (BL21) ===")
prev = set()
for rm in SWEEP:
    cur = KEYS.get((rm, "BL21"), set())
    new = cur - prev
    print(f"  rmsd {rm}: 총 {len(cur)}개 (신규 {len(new)})")
    prev = cur

# ---------- MinD 기준점 ----------
_bl = {}
for l in open(ASSET["faa_bl21"], errors="ignore"):
    if l.startswith(">"):
        h = l[1:].rstrip(); _bl[h.split()[0]] = h[len(h.split()[0]):]
mind = [k for k, v in _bl.items() if re.search(r"\bMinD\b|septum site-determining", str(v), re.I)]
print(f"\n=== 기준점: MinD (CooC1 의 최근접 BL21 구조 이웃) ===")
print(f"BL21 프로테옴에서 찾은 MinD: {mind}")
for k in mind:
    hit_at = [rm for rm in SWEEP if k in KEYS.get((rm, "BL21"), set())]
    print(f"  {k} {_bl[k][:60]}")
    print(f"    잡히는 임계값: {hit_at if hit_at else '어느 값에서도 안 잡힘'}")
if mind and not any(k in KEYS.get((SWEEP[0], "BL21"), set()) for k in mind):
    print("\n  → 1.0 A 에서 MinD 가 안 잡힌다. 두 가지 중 하나다:")
    print("     (a) 현재 컷오프가 너무 빡빡하다 -> 더 큰 값에서 잡히는지 위 표로 확인")
    print("     (b) MinD 는 금속 자리가 없다 (Walker A 만 공유) -> 애초에 기준점으로")
    print("         부적합. 그렇다면 metal 모티프의 민감도는 다른 방법으로 재야 한다.")
    print("     CooC1 의 Ni 자리는 CXC 인데 MinD 에는 그 Cys 가 없을 수 있다.")
    print("     위 '신규' 열이 임계값마다 몇 배로 느는지를 보고 특이도 손실을 판단할 것.")
S.to_csv(DIR["table"]/"folddisco_rmsd_sweep.csv", index=False, encoding="utf-8-sig")
print(f"\n저장: {DIR['table']/'folddisco_rmsd_sweep.csv'}")
```

---

## CELL 29c — 금속 질의의 민감도 진단 + Cys4 기하로 재질의

```python
# =============================================================================
# CELL 29c | 왜 29개밖에 안 나오는가 — 병목이 어디인지 찾고, 질의를 바꿔 본다
#   CELL 29b 결과: --rmsd 를 1.0->2.0 으로 풀어도 히트가 하나도 안 는다.
#   최대 min_rmsd 가 0.745 인데 그 위로 아무것도 없다. 즉 병목은 RMSD 필터가
#   아니라 그 앞단(인덱스 검색)이다. folddisco 는 기하 특징을 해시로 찾아온 뒤
#   RMSD 를 재는데, 가져온 후보가 전부 0.75 이하라는 뜻이다.
#   -> 임계값을 아무리 풀어도 검색 범위는 안 넓어진다. 질의 자체를 바꿔야 한다.
#
#   이 셀이 하는 것:
#     (1) 민감도 추정 — Boltz 가 Cys4 자리를 찾아준 단백질 중 folddisco 가 몇 개나
#         잡았나. 정답을 아는 소수 집합에서 재는 recall 이다.
#     (2) 알려진 Ni 단백질(hyp 오페론·SlyD·Fe-S 운반)이 29개 안에 있는가
#     (3) 재질의 — CooC1 의 A112/A114 대신, 실제로 Ni 을 잡은 Cys4 기하를 질의로
#         써서 세 균주를 다시 훑는다. 간격이 다른 자리를 잡기 위해서다.
# =============================================================================
need(((DIR["table"]/"ni_site_grade.csv").exists(), "CELL 28a-8b 를 먼저 돌릴 것"))
G  = pd.read_csv(DIR["table"]/"ni_site_grade.csv")
MET_SET = set(pd.read_csv(DIR["table"]/"folddisco_metal_motif.csv").protein) \
          if (DIR["table"]/"folddisco_metal_motif.csv").exists() else set()

# ---------- (1) recall ----------
A_ = G[G.grade == "A"].copy()
print("=== (1) Cys4 자리를 가진 단백질을 folddisco 가 잡았는가 ===")
print(A_[["strain", "prey", "donors", "folddisco_metal"]].to_string(index=False))
_miss = A_[A_.folddisco_metal == False]
print(f"\nfolddisco 금속 목록에 없던 것: {len(_miss)} / {len(A_)}")
if len(_miss):
    print(_miss[["strain", "prey", "donors"]].to_string(index=False))
print("""
※ 이 숫자를 recall(재현율)로 부르지 말 것. 세 가지 이유로 비율이 성립하지 않는다.
   1. 표본이 한 자릿수다. 신뢰구간이 사실상 0~90% 다.
   2. 분모가 정답이 아니다. Boltz 의 Cys4 배치는 또 다른 모델의 예측이지 실측이
      아니다. 이건 재현율이 아니라 두 방법의 '일치도'다.
   3. 분모가 편향돼 있다. out_focus 27개와 out_cmp 51개는 전부 folddisco 히트로
      만든 입력이다. 분모가 folddisco 쪽으로 기울어 일치도가 실제보다 높게 나온다.
   말할 수 있는 것: 질의 간격(A112/A114)과 다른 Cys4 자리가 실재하고 그것들이
   목록에 없었다. 사각지대가 있다는 것은 확실하고, 크기는 모른다.
   크기를 알려면 아래 (3) 재질의로 늘어나는 수를 봐야 한다.""")

# ---------- (2) 알려진 금속 단백질이 29개 안에 있나 ----------
def read_hdr(path):
    out = {}
    for l in open(path, errors="ignore"):
        if l.startswith(">"):
            h = l[1:].rstrip(); out[h.split()[0]] = h[len(h.split()[0]):].strip()
    return out
BLD = read_hdr(ASSET["faa_bl21"])
KNOWN = {
    "hyp 오페론 (Ni 삽입, KO 음성)": r"hydrogenase (matura|nickel)|\bHyp[A-F]\b|HybF",
    "SlyD (Ni 샤페론)":              r"\bSlyD\b|peptidyl-prolyl.*SlyD",
    "Ni 수송/조절":                  r"nickel|\bNik[ABCDER]\b|\bRcn[AR]\b",
    "Fe-S 운반":                     r"iron-sulfur cluster (insertion|assembly|carrier)|\bIsc[AUS]\b|\bSufA\b|\bNfuA\b|\bErpA\b",
    "COG0523 금속 샤페론":           r"\bYeiR\b|\bYjiA\b|\bCobW\b|\bZigA\b|COG0523",
    "urease":                        r"urease|\bUre[EFGDH]\b",
}
print("\n=== (2) 알려진 금속 관련 단백질이 metal 모티프 29개 안에 있는가 ===")
for label, pat in KNOWN.items():
    hits = [k for k, v in BLD.items() if re.search(pat, v, re.I)]
    inm  = [k for k in hits if k in MET_SET]
    print(f"\n[{label}]  프로테옴 {len(hits)}개 / 모티프 목록에 {len(inm)}개")
    for k in hits[:8]:
        print(f"   {'O' if k in MET_SET else '.'} {k}  {BLD[k][:72]}")
    if len(hits) > 8: print(f"   ... 외 {len(hits)-8}개")
print("\n  'O' 가 모티프 목록에 든 것. 알려진 Ni 단백질이 대부분 '.' 이면")
print("  이 질의의 민감도가 낮다는 직접 증거다.")

# ---------- (3) Cys4 기하로 재질의 ----------
# 질의 템플릿은 실제로 Ni 을 잡은 구조에서 가져온다. 구조 파일은 folddisco 인덱스가
# 가리키는 디렉터리에 있으므로 tid 로 찾는다.
print("\n" + "=" * 90)
print("=== (3) Cys4 기하로 재질의 ===")
_roots = [TOOLS/"database"/"bacteriaDB", TOOLS/"database", DIR["external"]]
def find_struct(pid):
    xw = DIR["table"]/"id_crosswalk_struct_to_genbank.csv"
    if not xw.exists(): return None
    x = pd.read_csv(xw)
    row = x[x.protein.astype(str).str.contains(str(pid), na=False)]
    if not len(row): return None
    stem = Path(str(row.iloc[0].tid)).stem
    for r in _roots:
        if not Path(r).exists(): continue
        for ext in ("pdb", "cif", "pdb.gz", "cif.gz"):
            f = next(iter(glob.glob(str(Path(r)/"**"/f"{stem}*.{ext}"), recursive=True)), None)
            if f: return f
    return None

QUERIES = {}
for _, r in A_.iterrows():
    pid = str(r.prey).split("-")[-1]
    if not pid.startswith("QJZ"): continue        # BL21 구조만 (인덱스 대응이 있다)
    res = [t for t in str(r.donors).split() if t.startswith("CYS")]
    if len(res) < 4: continue
    nums = sorted(int(t[3:]) for t in res)
    QUERIES[pid] = ",".join(f"A{n}" for n in nums)

print("질의 후보:", QUERIES if QUERIES else "없음")
FD  = ASSET["folddisco"]
IDX = {"BL21": ASSET["fd_idx_bl21"], "Y19": ASSET["fd_idx_y19"], "MG1655": ASSET["fd_idx_mg"]}
outd = DIR["folddisco"]/"cys4"; outd.mkdir(exist_ok=True)
found = {}
for pid, qres in QUERIES.items():
    sp = find_struct(pid)
    if not sp:
        print(f"  [구조 못 찾음] {pid} — 인덱스 디렉터리를 확인할 것"); continue
    print(f"\n>>> 질의 {pid}  {qres}\n    템플릿 {sp}")
    for strain, idx in IDX.items():
        out = outd/f"cys4_{pid}_{strain}.tsv"
        if not (out.exists() and out.stat().st_size):
            sh(f'"{FD}" query -p "{sp}" -q {qres} -i "{idx}" -t {min(THREADS,8)} '
               f'--per-structure --header --sort-by idf --rmsd 1.0 --top 5000 -o "{out}"',
               check=False)
        if out.exists() and out.stat().st_size:
            d = pd.read_csv(out, sep="\t"); d.columns = [c.strip().lstrip("#") for c in d.columns]
            found[(pid, strain)] = len(d)
            print(f"    {strain:7s} {len(d)}개")

if found:
    F = (pd.Series(found).rename("hits").rename_axis(["query", "strain"])
           .reset_index().pivot(index="query", columns="strain", values="hits"))
    print("\n=== Cys4 질의 히트 수 ===")
    print(F.to_string())
    F.to_csv(DIR["table"]/"folddisco_cys4_requery.csv", encoding="utf-8-sig")
    print(f"\n저장: {DIR['table']/'folddisco_cys4_requery.csv'}")
    print("\n원래 질의(A112/A114)의 29/21/30 과 비교할 것.")
    print("  크게 늘면 → 원래 질의가 좁았던 것이고, 새 목록을 후보에 합쳐야 한다.")
    print("  비슷하면 → 잔기 4개 질의는 원래 특이도가 높다. 지금 목록이 맞다.")
else:
    print("\n재질의를 못 돌렸다. 구조 파일 경로(_roots)를 확인할 것:")
    for r in _roots: print(f"   {r}  {'있음' if Path(r).exists() else '없음'}")
```

---

## CELL 29d — Cys4 재질의 결과 판정 + 균주 특이성과 교차

```python
# =============================================================================
# CELL 29d | 넓힌 질의가 무엇을 얻었고 무엇을 못 얻었는지 정리한다
#   CELL 29c 결과: Cys4 질의가 균주당 47~128개를 잡는다 (원래 질의는 29/21/30).
#   질의를 넓힌 것 자체는 성공이다. 문제는 넓힌 뒤에도 균주를 가르느냐다.
#   여기서 세 가지를 본다.
#     (1) 프로테옴 크기로 정규화한 비율 — 균주마다 다른가
#     (2) 벤치마크 양성 — HypA·HypB·Nik·IscA 처럼 알려진 Ni/금속 단백질이
#         넓힌 질의에서는 잡히는가. 원래 질의는 10개 중 0개였다
#     (3) ★ Cys4 x BL21 균주 특이 교차 — 두 축을 모두 만족하는 후보
#         'Cys4 자리를 가졌고 + MG1655 에 상동체가 없거나 멀다'
#         QJZ12568.1 이 그 예다. 다른 것이 몇 개나 더 있는지 본다
# =============================================================================
need(((DIR["folddisco"]/"cys4").exists(), "CELL 29c 를 먼저 돌릴 것"))
FAA = {"BL21": ASSET["faa_bl21"], "Y19": ASSET["faa_y19"], "MG1655": ASSET["faa_mg1655"]}
NPROT = {}
HDR = {}
for st, fa in FAA.items():
    h = {}
    for l in open(fa, errors="ignore"):
        if l.startswith(">"):
            t = l[1:].rstrip(); h[t.split()[0]] = t[len(t.split()[0]):].strip()
    HDR[st] = h; NPROT[st] = len(h)

xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(xw) else {}

CY = {}      # (query, strain) -> set(key)
for f in sorted((DIR["folddisco"]/"cys4").glob("cys4_*.tsv")):
    if not f.stat().st_size: continue
    q, st = re.match(r"cys4_(.+)_(BL21|Y19|MG1655)\.tsv", f.name).groups()
    d = pd.read_csv(f, sep="\t"); d.columns = [c.strip().lstrip("#") for c in d.columns]
    d["stem"] = d.tid.apply(lambda x: Path(str(x)).stem)
    d["prot"] = d.stem.map(xwm)
    d["uni"]  = d.stem.astype(str).str.extract(r"AF-([A-Z0-9]+)-F", expand=False)
    d["key"]  = d.prot.where(d.prot.notna(), d.uni).astype(str).str.split(",").str[0]
    d = d[[k in HDR[st] for k in d.key]] if st != "MG1655" else d
    CY[(q, st)] = set(d.key)

# ---------- (1) 정규화 ----------
print("=== (1) Cys4 질의 히트 — 프로테옴 대비 ===")
rows = []
for (q, st), ks in CY.items():
    rows.append({"query": q, "strain": st, "hits": len(ks),
                 "proteome": NPROT[st], "pct": round(100*len(ks)/NPROT[st], 2)})
R = pd.DataFrame(rows)
print(R.pivot(index="query", columns="strain", values="hits").to_string())
print("\n비율(%)")
print(R.pivot(index="query", columns="strain", values="pct").to_string())
print("\n원래 질의(A112/A114): BL21 29 (0.71%) / MG1655 21 (0.49%) / Y19 30 (0.56%)")
_sp = R.groupby("query").pct.agg(lambda x: x.max() - x.min())
print(f"\n질의별 균주 간 비율 차이(최대-최소): " +
      ", ".join(f"{k} {v:.2f}%p" for k, v in _sp.items()))
print("  ← 0 에 가까우면 Cys4 모티프 보유율은 세 균주가 같다는 뜻이다.")
print("    질의를 넓혀도 균주를 가르지 못한다.")

# ---------- (2) 벤치마크 양성 ----------
UNION = {st: set().union(*[v for (q, s), v in CY.items() if s == st]) for st in FAA}
BENCH = {
  "HypA/HypB (Ni 메탈로샤페론)": r"nickel metallochaperone|nickel incorporation|\bHypA\b|\bHypB\b",
  "Nik 수송체":                  r"\bNik[ABCDE]\b|nickel ABC|nickel import",
  "Fe-S 운반":                   r"iron-sulfur cluster (insertion|assembly|carrier)|\bIscA\b|\bIscU\b|\bNfuA\b|\bApbC\b",
}
print("\n" + "=" * 90)
print("=== (2) 벤치마크: 알려진 금속 단백질이 넓힌 질의에서는 잡히는가 (BL21) ===")
_m3 = pd.read_csv(DIR["table"]/"motif_metal_3strain.csv") if (DIR["table"]/"motif_metal_3strain.csv").exists() else pd.DataFrame()
OLD = set(_m3[_m3.strain == "BL21"].key.astype(str)) if len(_m3) else set()
for label, pat in BENCH.items():
    hits = [k for k, v in HDR["BL21"].items() if re.search(pat, v, re.I)]
    print(f"\n[{label}] 프로테옴 {len(hits)}개")
    for k in hits:
        o = "O" if k in OLD else "."
        n = "O" if k in UNION["BL21"] else "."
        print(f"   원래 {o} / Cys4 {n}   {k}  {HDR['BL21'][k][:66]}")
print("\n  원래 질의에서는 HypA·HypB·Nik 이 전부 '.' 이었다. Cys4 쪽에서 'O' 로")
print("  바뀌면 넓힌 질의의 민감도가 실제로 올라간 것이다.")

# ---------- (3) Cys4 x BL21 균주 특이 ----------
# ★ 분류표는 방향이 두 가지로 저장돼 있다. classify_MG1655_vs_BL21DE3.csv 는
#   MG1655 단백질로 키가 잡혀 있어(4,300행) BL21 의 QJZ* 와 대조하면 교집합이
#   0 이 된다 — 실제로 그렇게 나왔었다. 파일 이름이 아니라 내용으로 고른다.
print("\n" + "=" * 90)
print("=== (3) Cys4 자리 + BL21 균주 특이 (두 축 동시) ===")
cls, chosen = None, None
for _f in sorted(DIR["table"].glob("classify_*.csv")):
    _d = pd.read_csv(_f)
    if "protein" not in _d.columns: continue
    _ov = sum(1 for k in _d.protein.astype(str) if k in HDR["BL21"])
    print(f"  {_f.name}: {len(_d)}행, BL21 ID 와 일치 {_ov}")
    if _ov > 0.5 * len(_d): cls, chosen = _d, _f
if cls is None:
    print("\n⚠ BL21 단백질로 키가 잡힌 분류표가 없다. 이 교차는 판정 불가.")
    print("   CELL 09~12 를 BL21 을 prey 로 두고 다시 돌려야 한다.")
    print("   (지금 있는 표는 MG1655 기준이라 BL21 의 균주 특이 여부를 말해주지 않는다)")
else:
    print(f"\n사용: {chosen.name}  ({len(cls)}행)  카테고리 {cls.category.value_counts().to_dict()}")
    SPEC = {"strain-specific", "low-similarity"}
    tgt = cls[cls.category.isin(SPEC)].set_index("protein")
    hit = sorted(UNION["BL21"] & set(tgt.index))
    print(f"\nBL21 Cys4 히트 {len(UNION['BL21'])}개 중 균주 특이·저유사도: **{len(hit)}개**")
    if hit:
        T = pd.DataFrame({"protein": hit,
                          "category": [tgt.loc[k, "category"] for k in hit],
                          "pident": [tgt.loc[k].get("pident") for k in hit],
                          "in_old_motif": [k in OLD for k in hit],
                          "queries": ["+".join(sorted(q for (q, s), v in CY.items()
                                                      if s == "BL21" and k in v)) for k in hit],
                          "desc": [HDR["BL21"].get(k, "")[:70] for k in hit]})
        T = T.sort_values(["category", "protein"])
        T.to_csv(DIR["table"]/"cys4_bl21_strain_specific.csv", index=False, encoding="utf-8-sig")
        pd.set_option("display.max_rows", None); pd.set_option("display.width", 220)
        print(T.to_string(index=False))
        print(f"\n저장: {DIR['table']/'cys4_bl21_strain_specific.csv'}")
        print("\n※ 이 목록이 '금속 자리 + BL21 특이' 두 축을 모두 만족하는 후보다.")
        print("  QJZ12568.1 이 여기 있어야 한다 (없으면 분류표 ID 계열을 확인할 것).")
    else:
        print("없다. Cys4 자리를 가진 BL21 단백질은 전부 MG1655 에도 상동체가 있다.")
        print("→ 금속 자리 보유로는 균주를 가를 수 없다는 결론이 한 번 더 확인된 것이다.")
```

---

## CELL 29e — Cys4 히트의 균주 특이성 (분류표 없이 mmseqs 로 직접)

```python
# =============================================================================
# CELL 29e | Cys4 자리를 가진 BL21 단백질 중 MG1655 에 없는 것을 찾는다
#   CELL 29d 의 분류표 교차는 ID 계열이 안 맞아 판정 불가였다. 분류표를 기다리지
#   않고 mmseqs 로 직접 묻는다 — "이 단백질이 MG1655 프로테옴에 있는가".
#   판정 기준은 기존 분류와 같은 뜻으로 맞춘다.
#     히트 없음            -> strain-specific   (MG1655 에 상동체 없음)
#     fident < 0.90        -> low-similarity    (있지만 멀다)
#     0.90 <= fident < 1.0 -> high-similarity
#     fident = 1.000       -> identical
#   앞의 둘이 '균주 특이' 축이다. 여섯 축 중 유일하게 아직 안 본 축이다.
# =============================================================================
need(((DIR["folddisco"]/"cys4").exists(), "CELL 29c 를 먼저 돌릴 것"))
FAA = {"BL21": ASSET["faa_bl21"], "Y19": ASSET["faa_y19"], "MG1655": ASSET["faa_mg1655"]}
SEQ, HDR = {}, {}
for st, fa in FAA.items():
    s_, h_, nm, buf = {}, {}, None, []
    for l in open(fa, errors="ignore"):
        if l.startswith(">"):
            if nm: s_[nm] = "".join(buf)
            t = l[1:].rstrip(); nm = t.split()[0]; buf = []
            h_[nm] = t[len(nm):].strip()
        else: buf.append(l.strip())
    if nm: s_[nm] = "".join(buf)
    SEQ[st], HDR[st] = s_, h_

xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(xw) else {}

# BL21 Cys4 히트 합집합
CY4, PERQ = set(), {}
for f in sorted((DIR["folddisco"]/"cys4").glob("cys4_*_BL21.tsv")):
    if not f.stat().st_size: continue
    q = re.match(r"cys4_(.+)_BL21\.tsv", f.name).group(1)
    d = pd.read_csv(f, sep="\t"); d.columns = [c.strip().lstrip("#") for c in d.columns]
    d["stem"] = d.tid.apply(lambda x: Path(str(x)).stem)
    d["key"]  = d.stem.map(xwm).astype(str).str.split(",").str[0]
    ks = {k for k in d.key if k in SEQ["BL21"]}
    PERQ[q] = ks; CY4 |= ks
    print(f"  {q}: {len(ks)}개")
print(f"\nBL21 Cys4 합집합 {len(CY4)}개")
need((len(CY4) > 0, "Cys4 히트가 비었다 — CELL 29c 출력을 확인할 것"))

QF = DIR["seq"]/"bl21_cys4.faa"
QF.write_text("".join(f">{k}\n{SEQ['BL21'][k]}\n" for k in sorted(CY4)))
_tmp = DIR["tmp"]/"cys4map"; _tmp.mkdir(parents=True, exist_ok=True)
FMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits"
BEST = {}
for st in ["MG1655", "Y19"]:
    out = DIR["search"]/f"cys4_bl21_vs_{st}.m8"
    if not (out.exists() and out.stat().st_size):
        sh(f'mmseqs easy-search "{QF}" "{FAA[st]}" "{out}" "{_tmp}" '
           f'--format-output "{FMT}" -e 1e-3 --threads {min(THREADS,16)} -v 1', check=False)
    if out.exists() and out.stat().st_size:
        d = pd.read_csv(out, sep="\t", names=FMT.split(","))
        d = d[d.qcov >= 0.70].sort_values("bits", ascending=False).drop_duplicates("query")
        BEST[st] = d.set_index("query")

def verdict(k, st):
    d = BEST.get(st)
    if d is None or k not in d.index: return "strain-specific", None
    f = float(d.loc[k, "fident"])
    return ("identical" if f >= 0.999 else
            "high-similarity" if f >= 0.90 else "low-similarity"), round(f, 3)

rows = []
for k in sorted(CY4):
    c_mg, f_mg = verdict(k, "MG1655")
    c_y, f_y   = verdict(k, "Y19")
    rows.append({"protein": k, "len": len(SEQ["BL21"][k]),
                 "vs_MG1655": c_mg, "fident_MG1655": f_mg,
                 "vs_Y19": c_y, "fident_Y19": f_y,
                 "queries": "+".join(q for q, s in PERQ.items() if k in s),
                 "desc": HDR["BL21"].get(k, "")[:75]})
V = pd.DataFrame(rows)
V["bl21_specific"] = V.vs_MG1655.isin(["strain-specific", "low-similarity"]).astype(int)
V = V.sort_values(["bl21_specific", "vs_MG1655", "protein"], ascending=[False, True, True])
V.to_csv(DIR["table"]/"cys4_bl21_strain_specificity.csv", index=False, encoding="utf-8-sig")

print("\n=== MG1655 대비 판정 ===")
print(V.vs_MG1655.value_counts().to_string())
_hot = V[V.bl21_specific == 1]
print(f"\n=== ★ Cys4 자리 + MG1655 에 없거나 먼 것: {len(_hot)}개 ===")
pd.set_option("display.max_rows", None); pd.set_option("display.width", 230)
if len(_hot):
    print(_hot[["protein", "len", "vs_MG1655", "fident_MG1655", "vs_Y19",
                "queries", "desc"]].to_string(index=False))
    print("\n두 축(금속 자리 + 균주 특이)을 모두 만족하는 후보다. 여섯 축 중")
    print("유일하게 균주를 가르는 축이 될 수 있다. 다만 이 목록 자체가")
    print("'Cys4 기하'라는 넓은 그물에서 나온 것이므로 개별 검증이 필요하다.")
else:
    print("없다. Cys4 자리를 가진 BL21 단백질 전부가 MG1655 에 가까운 상동체를 갖는다.")
    print("→ 금속 자리 보유로는 균주를 가를 수 없다는 결론이 여섯 번째로 확인된 것이다.")
print(f"\n저장: {DIR['table']/'cys4_bl21_strain_specificity.csv'}")
```

---

## CELL 29f — CooC1 이량체 계면 질의 (A112,A114 + B112,B114)

```python
# =============================================================================
# CELL 29f | CooC1 은 호모다이머다. Ni 자리를 계면으로 놓고 다시 묻는다
#   CooC1 의 CXC(112/114)가 이량체 계면에서 반대편 CXC 와 만나 Cys4 를 이룬다면,
#   진짜 질의는 A112,A114 가 아니라 A112,A114,B112,B114 다. 그러면 HslO·GspE 의
#   Cys4 기하를 빌려 쓸 필요 없이 CooC1 자신의 기하로 4잔기 질의가 된다.
#
#   그리고 원래 의도대로 된다 — 2잔기로 넓게 뽑고 4잔기로 좁히는 것. 4잔기 질의가
#   2잔기의 부분집합이므로(잔기를 포함하므로) 순서가 성립한다. HslO 템플릿은
#   잔기가 아예 달라서 부분집합이 아니었고, 그래서 오히려 히트가 늘었다.
#
#   전제: PDB 파일에 chain B 가 있어야 한다. 비대칭 단위가 단량체면 생물학적
#   조립체(.pdb1)를 받아야 한다. 아래에서 확인하고 없으면 받는다.
# =============================================================================
need((have("ASSET", "sh"), "CELL 01 을 먼저 돌릴 것"))
import urllib.request, gzip
PDB = Path(ASSET["cooc1_pdb"])
print("현재 질의 구조:", PDB)

def chains_of(p):
    ch = {}
    op = gzip.open if str(p).endswith(".gz") else open
    for l in op(p, "rt", errors="ignore"):
        if l.startswith(("ATOM", "HETATM")) and len(l) > 26:
            c, comp, seq = l[21], l[17:20].strip(), l[22:26].strip()
            ch.setdefault(c, {})[seq] = comp
    return ch

ch = chains_of(PDB)
print("체인:", {c: len(v) for c, v in ch.items()})
for c, v in ch.items():
    print(f"  chain {c}: 112={v.get('112','없음')}  114={v.get('114','없음')}")

TARGET = PDB
if len(ch) < 2:
    print("\n단량체다. RCSB 에서 생물학적 조립체(.pdb1)를 받는다.")
    _id = re.search(r"([0-9][a-zA-Z0-9]{3})", PDB.name)
    if not _id:
        print("  ⚠ 파일명에서 PDB ID 를 못 읽었다. 직접 지정할 것.")
    else:
        pid = _id.group(1).lower()
        dst = DIR["external"]/f"{pid}_bioassembly.pdb"
        if not dst.exists():
            try:
                with urllib.request.urlopen(
                        f"https://files.rcsb.org/download/{pid.upper()}.pdb1.gz",
                        timeout=120) as r:
                    dst.write_bytes(gzip.decompress(r.read()))
                print(f"  받음 -> {dst}")
            except Exception as e:
                print(f"  ⚠ 내려받기 실패: {e}")
        if dst.exists():
            ch2 = chains_of(dst)
            print("  조립체 체인:", {c: len(v) for c, v in ch2.items()})
            for c, v in ch2.items():
                print(f"    chain {c}: 112={v.get('112','없음')}  114={v.get('114','없음')}")
            if len(ch2) >= 2: TARGET, ch = dst, ch2

# ---- 질의 실행 ----
IDX = {"BL21": ASSET["fd_idx_bl21"], "Y19": ASSET["fd_idx_y19"], "MG1655": ASSET["fd_idx_mg"]}
FD  = ASSET["folddisco"]
outd = DIR["folddisco"]/"dimer"; outd.mkdir(exist_ok=True)

cs = [c for c in sorted(ch) if ch[c].get("112") == "CYS" and ch[c].get("114") == "CYS"]
print(f"\n112/114 가 CYS 인 체인: {cs}")
QSET = {}
if len(cs) >= 2:
    a, b = cs[0], cs[1]
    QSET["dimer4"] = f"{a}112,{a}114,{b}112,{b}114"
    QSET["dimer3"] = f"{a}112,{a}114,{b}112"      # 4잔기가 너무 빡빡할 때의 중간 단계
QSET["mono2"] = f"{cs[0]}112,{cs[0]}114" if cs else "A112,A114"
print("질의:", QSET)

res = {}
for tag, q in QSET.items():
    for strain, idx in IDX.items():
        out = outd/f"{tag}_{strain}.tsv"
        if not (out.exists() and out.stat().st_size):
            sh(f'"{FD}" query -p "{TARGET}" -q {q} -i "{idx}" -t {min(THREADS,8)} '
               f'--per-structure --header --sort-by idf --rmsd 1.0 --top 20000 -o "{out}"',
               check=False)
        n = 0
        if out.exists() and out.stat().st_size:
            d = pd.read_csv(out, sep="\t"); n = len(d)
        res[(tag, strain)] = n
if res:
    D = (pd.Series(res).rename("hits").rename_axis(["query", "strain"]).reset_index()
           .pivot(index="query", columns="strain", values="hits"))
    print("\n=== 질의별 히트 수 ===")
    print(D.to_string())
    D.to_csv(DIR["table"]/"folddisco_dimer_query.csv", encoding="utf-8-sig")
    print("\n읽는 법")
    print("  dimer4 < dimer3 < mono2 이면 의도대로 좁혀진 것이다.")
    print("  dimer4 가 0 이면 4잔기 제약이 너무 빡빡하다 — dimer3 를 쓴다.")
    print("  mono2 가 여전히 29개면, 문제는 잔기 수가 아니라 인덱스 검색 단계다.")

# ---- 2잔기 질의가 왜 29개뿐인지: 검색 단계 옵션 확인 ----
print("\n" + "=" * 80)
print("=== folddisco query 옵션 (검색 단계 조절 인자를 찾는다) ===")
print("--rmsd 는 후처리 필터라 풀어도 안 늘었다 (CELL 29b). 늘리려면 후보를")
print("가져오는 단계의 인자를 건드려야 한다. 아래에서 그런 옵션을 찾을 것:")
sh(f'"{FD}" query --help', check=False)
```

---

## CELL 29g — 게놈으로 확인 (어노테이션 누락 · pseudogene 구분)

```python
# =============================================================================
# CELL 29g | "MG1655 에 없다"를 게놈에서 다시 묻는다
#   CELL 29e 는 MG1655 의 '어노테이션된 프로테옴'을 봤다. 그런데 유전자가 게놈에
#   있어도 (a) 어노테이션이 놓쳤거나 (b) pseudogene 으로 분류돼 단백질 FASTA 에서
#   빠졌으면 똑같이 '없음'으로 나온다. 세 경우가 뭉쳐 있다.
#
#     단백질 없음 + 게놈 없음        -> 진짜 균주 특이           [강한 후보]
#     단백질 없음 + 게놈 있음(온전)  -> 어노테이션 누락          [후보 아님]
#     단백질 없음 + 게놈 있음(깨짐)  -> pseudogene               [★ 가장 흥미롭다]
#
#   세 번째는 "조상은 둘 다 가졌는데 K-12 에서 망가졌다"는 뜻이고, 그게 사실이면
#   BL21 lysate 만 rescue 하는 표현형을 그대로 설명한다. K-12 는 B 계열 대비
#   IS 삽입·프레임시프트로 깨진 유전자가 적지 않다.
#
#   tblastn 대신 mmseqs --search-type 2 (번역 검색) 를 쓴다. 새 도구가 필요 없다.
# =============================================================================
need(((DIR["table"]/"cys4_bl21_strain_specificity.csv").exists(),
      "CELL 29e 를 먼저 돌릴 것"))
import urllib.request

# ---- 게놈 확보 ----
GENOME = {"MG1655": ("U00096.3", DIR["external"]/"MG1655_U00096.3.fna")}
for name, (acc, dst) in GENOME.items():
    if dst.exists() and dst.stat().st_size > 1_000_000:
        print(f"{name}: 있음 {dst.name} ({dst.stat().st_size/1e6:.1f} MB)"); continue
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={acc}&rettype=fasta&retmode=text")
    try:
        with urllib.request.urlopen(url, timeout=300) as r:
            dst.write_bytes(r.read())
        print(f"{name}: 받음 {acc} -> {dst} ({dst.stat().st_size/1e6:.1f} MB)")
    except Exception as e:
        print(f"{name}: ⚠ 내려받기 실패 {e}")
        print(f"   수동: curl -o '{dst}' '{url}'")

V = pd.read_csv(DIR["table"]/"cys4_bl21_strain_specificity.csv")
TGT = V[V.vs_MG1655 == "strain-specific"].protein.astype(str).tolist()
print(f"\n프로테옴에서 'MG1655 에 없음' 으로 나온 것: {len(TGT)}개")
if not TGT:
    print("없다. 게놈 확인이 필요한 대상이 없다 — 29e 결론을 그대로 쓴다.")
else:
    print(TGT)

_g = GENOME["MG1655"][1]
if TGT and _g.exists():
    SEQB = {}
    nm, buf = None, []
    for l in open(ASSET["faa_bl21"], errors="ignore"):
        if l.startswith(">"):
            if nm: SEQB[nm] = "".join(buf)
            nm, buf = l[1:].split()[0], []
        else: buf.append(l.strip())
    if nm: SEQB[nm] = "".join(buf)

    q = DIR["seq"]/"cys4_specific_for_genome.faa"
    q.write_text("".join(f">{k}\n{SEQB[k]}\n" for k in TGT if k in SEQB))
    out = DIR["search"]/"cys4_specific_vs_MG1655genome.m8"
    _tmp = DIR["tmp"]/"gsearch"; _tmp.mkdir(parents=True, exist_ok=True)
    FMT = "query,target,fident,alnlen,qlen,qstart,qend,tstart,tend,evalue,bits"
    if not (out.exists() and out.stat().st_size):
        # --search-type 2 : 단백질 질의 vs 핵산 대상 (tblastn 계열)
        sh(f'mmseqs easy-search "{q}" "{_g}" "{out}" "{_tmp}" --search-type 2 '
           f'--format-output "{FMT}" -e 1e-3 -s 7.5 --threads {min(THREADS,16)} -v 1',
           check=False)

    print("\n" + "=" * 100)
    print("=== 게놈 대조 결과 ===")
    if not (out.exists() and out.stat().st_size):
        print("히트 없음 — 질의 전부가 MG1655 게놈에 없다.")
        print("→ 진짜 균주 특이. 프로테옴 판정이 그대로 유지된다.")
    else:
        d = pd.read_csv(out, sep="\t", names=FMT.split(","))
        d["qcov"] = (d.qend - d.qstart + 1) / d.qlen
        best = d.sort_values("bits", ascending=False).drop_duplicates("query")
        rows = []
        for k in TGT:
            if k not in set(best["query"]):
                rows.append({"protein": k, "genome": "없음",
                             "verdict": "진짜 균주 특이 [강한 후보]",
                             "fident": None, "qcov": None, "n_frag": 0}); continue
            r = best[best["query"] == k].iloc[0]
            nfrag = int((d["query"] == k).sum())
            # pseudogene 의 서명은 '같은 유전자인데 깨짐' 이다 -> fident 가 높아야 한다.
            # fident 0.4 대는 깨진 사본이 아니라 먼 상동체(파라로그)다. 예: K-12 의
            # dksA 는 TraR 의 파라로그이고, K-12 자신의 프로파지 terminase 도 있다.
            # fident 조건 없이 qcov/n_frag 만 보면 그것들이 전부 pseudogene 으로 찍힌다.
            if r.fident < 0.80:
                v = "먼 상동체 — 이 유전자가 아니라 파라로그일 가능성 [균주 특이 유지]"
            elif r.qcov >= 0.90 and nfrag <= 2:
                v = "어노테이션 누락 — 게놈에 온전히 있다 [후보 아님]"
            elif r.qcov < 0.70:
                # 잘린 것만 pseudogene 으로 본다. n_frag 는 번역 검색의 프레임
                # 경계로도 늘어나므로 단독 근거로 쓰지 않는다 (실측: fident 1.000,
                # qcov 0.933 인 온전한 유전자가 조각 4개로 나왔다).
                v = f"★ pseudogene 의심 — 질의의 {r.qcov:.0%} 만 게놈에 있다"
            elif nfrag >= 3:
                v = f"온전하나 조각 {nfrag}개 — 번역 프레임 경계일 가능성. 좌표 확인"
            else:
                v = "부분 일치 — 좌표를 직접 볼 것"
            rows.append({"protein": k, "genome": "있음", "verdict": v,
                         "fident": round(float(r.fident), 3),
                         "qcov": round(float(r.qcov), 3), "n_frag": nfrag})
        G = pd.DataFrame(rows)
        G["desc"] = G.protein.map(V.set_index("protein").desc)
        pd.set_option("display.max_rows", None); pd.set_option("display.width", 230)
        print(G.to_string(index=False))
        G.to_csv(DIR["table"]/"cys4_genome_check.csv", index=False, encoding="utf-8-sig")
        print("\n=== 판정별 ===")
        print(G.verdict.value_counts().to_string())
        print("\n★ 이 판정은 표시일 뿐 후보에서 빼지 않는다. 어떤 행도 삭제되지 않는다.")
        print("  '어노테이션 누락' 조차 제외 사유가 아니다 — 이 계의 차이는 유전자")
        print("  구성이 아니라 발현량 쪽이라는 것이 이미 확인됐으므로(4-12), MG1655 에")
        print("  유전자가 있다는 사실이 후보를 배제하지 않는다.")
        print("\n※ 판정 기준: fident 가 낮으면(<0.80) 깨진 사본이 아니라 다른 유전자다.")
        print("  fident 가 높은데 qcov 가 낮거나 조각이 많으면 프레임시프트·IS 삽입을")
        print("  의심한다. n_frag 는 번역 검색의 프레임 경계로도 늘 수 있으니 참고값이다.")
        print("  그 경우 해당 좌표를 뽑아 직접 확인할 것 — 이게 사실이면 '조상은 둘 다")
        print("  가졌는데 K-12 에서 깨졌다' 가 되어 표현형을 그대로 설명한다.")
        print(f"\n저장: {DIR['table']/'cys4_genome_check.csv'}")
elif TGT:
    print("\n게놈 파일이 없어 확인을 못 했다. 위 curl 명령으로 받은 뒤 다시 돌릴 것.")
```

---

## CELL 29h — 검색 허용폭을 풀고(넓게) 4잔기 전부 맞은 것만 남긴다(좁게)

```python
# =============================================================================
# CELL 29h | 원래 의도대로 — 2잔기로 넓게, 4잔기로 좁게
#   --help 로 드러난 사실 셋:
#     -d/--distance [0.5] · -a/--angle [5.0]  검색(해시) 단계 허용폭. 진짜 병목이다.
#     --covered-node [0]                       질의 잔기 중 몇 개가 맞아야 보고할지.
#                                              기본 0 이라 4잔기 질의도 2개만 맞으면
#                                              보고된다. dimer4 가 mono2 와 같았던 이유.
#     --rmsd [no limit]                        기본이 제한 없음. 우리가 1.0 을 준 것이
#                                              오히려 필터를 추가한 셈이었다.
#   그래서:
#     1단계  mono2 + 넉넉한 -d/-a       -> 후보를 크게 늘린다
#     2단계  dimer4 + --covered-node 4  -> 네 잔기가 모두 맞은 것만 남긴다
#   CooC1 은 호모다이머이고 A112/A114 가 반대편 것과 계면에서 Cys4 를 이룬다.
#   dimer4 는 그 자신의 기하이므로 HslO/GspE 템플릿을 빌릴 필요가 없다.
# =============================================================================
need((have("ASSET", "sh"), "CELL 01 을 먼저 돌릴 것"))
PDBQ = Path(globals().get("TARGET", ASSET["cooc1_pdb"]))
IDX  = {"BL21": ASSET["fd_idx_bl21"], "Y19": ASSET["fd_idx_y19"], "MG1655": ASSET["fd_idx_mg"]}
FD   = ASSET["folddisco"]
outd = DIR["folddisco"]/"widen"; outd.mkdir(exist_ok=True)
FMT  = "tid,idf,min_rmsd,nres,plddt,node_count,matching_residues"

# (d, angle) 조합. 기본은 (0.5, 5.0)
GRID  = [(0.5, 5.0), (1.0, 10.0), (1.5, 15.0)]
MONO  = "A112,A114"
DIMER = "A112,A114,B112,B114"

def run(tag, q, strain, d, a, extra=""):
    out = outd/f"{tag}_d{d}_a{a}_{strain}.tsv"
    if not (out.exists() and out.stat().st_size):
        sh(f'"{FD}" query -p "{PDBQ}" -q {q} -i "{IDX[strain]}" -t {min(THREADS,8)} '
           f'--per-structure --header --sort-by idf -d {d} -a {a} {extra} '
           f'--format-output "{FMT}" -o "{out}"', check=False)
    if not (out.exists() and out.stat().st_size): return None
    df = pd.read_csv(out, sep="\t")
    df.columns = [c.strip().lstrip("#") for c in df.columns]
    return df

# ---------- 1단계: 허용폭을 풀어 2잔기로 넓게 ----------
print("=== 1단계: mono2 (A112,A114) — 검색 허용폭을 푼다 ===")
W = []
for d, a in GRID:
    for st in IDX:
        df = run("mono2", MONO, st, d, a)
        W.append({"d": d, "angle": a, "strain": st, "hits": 0 if df is None else len(df)})
        print(f"  d={d} a={a} {st:7s} {W[-1]['hits']}")
Wd = pd.DataFrame(W).pivot_table(index=["d", "angle"], columns="strain", values="hits")
print("\n" + Wd.to_string())
print("\n기본 (0.5, 5.0) 대비 얼마나 늘었는지가 핵심이다. 안 늘면 병목이 또 다른 데 있다.")

# ---------- 2단계: 4잔기 전부 맞은 것만 ----------
print("\n" + "=" * 90)
print("=== 2단계: dimer4 (A112,A114,B112,B114) + --covered-node 4 ===")
print("질의 잔기 4개가 모두 맞은 구조만 남긴다. 부분 매치를 세지 않는다.\n")
KEEP = {}
D = []
for d, a in GRID:
    for st in IDX:
        df = run("dimer4cov", DIMER, st, d, a, extra="--covered-node 4")
        n = 0 if df is None else len(df)
        D.append({"d": d, "angle": a, "strain": st, "hits": n})
        if df is not None and len(df): KEEP[(d, a, st)] = df
        print(f"  d={d} a={a} {st:7s} {n}")
Dd = pd.DataFrame(D).pivot_table(index=["d", "angle"], columns="strain", values="hits")
print("\n" + Dd.to_string())

# ---------- 깔때기 ----------
print("\n=== 깔때기 (BL21) ===")
for d, a in GRID:
    m = Wd.loc[(d, a), "BL21"] if (d, a) in Wd.index else None
    k = Dd.loc[(d, a), "BL21"] if (d, a) in Dd.index else None
    if m is not None and k is not None:
        print(f"  d={d} a={a}:  2잔기 {int(m):5d}  ->  4잔기 전부 맞음 {int(k):4d}"
              f"  ({100*k/max(1,m):.1f}%)")
print("\n이게 '2개로 넓게 뽑고 4개로 줄인다' 의 실제 모습이다.")

# ---------- 가장 넓은 설정의 dimer4 결과를 이름과 함께 ----------
xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(xw) else {}
hdr = {}
for l in open(ASSET["faa_bl21"], errors="ignore"):
    if l.startswith(">"):
        t = l[1:].rstrip(); hdr[t.split()[0]] = t[len(t.split()[0]):].strip()

best = max((k for k in KEEP if k[2] == "BL21"), key=lambda k: (k[0], k[1]), default=None)
if best:
    df = KEEP[best].copy()
    df["stem"] = df.tid.apply(lambda x: Path(str(x)).stem)
    df["protein"] = df.stem.map(xwm).astype(str).str.split(",").str[0]
    df["desc"] = df.protein.map(lambda k: hdr.get(k, "")[:70])
    df = df.sort_values("min_rmsd")
    cols = [c for c in ["protein", "idf", "min_rmsd", "nres", "plddt",
                        "node_count", "matching_residues", "desc"] if c in df.columns]
    pd.set_option("display.max_rows", None); pd.set_option("display.width", 240)
    print(f"\n=== BL21 · d={best[0]} a={best[1]} · 4잔기 전부 맞음: {len(df)}개 ===")
    print(df[cols].to_string(index=False))
    df.to_csv(DIR["table"]/"folddisco_dimer4_widened.csv", index=False, encoding="utf-8-sig")
    print(f"\n저장: {DIR['table']/'folddisco_dimer4_widened.csv'}")
    print("\n※ 이 목록은 CooC1 자신의 이량체 Ni 자리 기하를 네 잔기 모두 만족한다.")
    print("  지금까지 중 가장 근거가 좁고 분명한 집합이다. CELL 29e 의 균주 특이")
    print("  판정과 교차해서 보면 된다.")
else:
    print("\n4잔기 전부 맞은 구조가 없다. --covered-node 3 으로 한 단계 낮춰 볼 것.")
```

---

## CELL 29i — dimer4 가 0 인 이유 가르기 (covered-node 기울기)

```python
# =============================================================================
# CELL 29i | 4잔기 전부 맞은 구조가 0 인 것이 사실인가, 질의가 성립을 안 한 것인가
#   AlphaFold 구조 DB 는 전부 단량체 예측이다. CooC1 의 A112,A114 + B112,B114 는
#   이량체 계면을 가로지르는 배치이고, 그런 배치는 단량체 모델에 존재할 수 없다.
#   0 이 나온 것이
#     (a) 기하학적 사실 — 그런 Cys4 배치를 가진 단백질이 없다
#     (b) 구조적 모순 — 단량체 인덱스에 이량체 질의를 던졌다
#   중 어느 쪽인지 지금 데이터로는 구분이 안 된다.
#
#   가르는 법: --covered-node 를 4 -> 3 -> 2 로 낮추며 기울기를 본다.
#     3 에서 상당수 나오면 (a) 쪽 — 4잔기 배치만 없는 것이다
#     3 에서도 0 이면 (b) 쪽 — 질의 자체가 인덱스와 안 맞는다
#   그리고 필터 없이 node_count 분포를 보면 실제로 몇 잔기가 맞고 있는지 나온다.
# =============================================================================
need((have("ASSET", "sh"), "CELL 01 을 먼저 돌릴 것"))
PDBQ  = Path(globals().get("TARGET", ASSET["cooc1_pdb"]))
IDX   = {"BL21": ASSET["fd_idx_bl21"], "Y19": ASSET["fd_idx_y19"], "MG1655": ASSET["fd_idx_mg"]}
FD    = ASSET["folddisco"]
outd  = DIR["folddisco"]/"covnode"; outd.mkdir(exist_ok=True)
FMT   = "tid,idf,min_rmsd,nres,plddt,node_count,max_node_cov,matching_residues"
DIMER = "A112,A114,B112,B114"
D, A  = 1.5, 15.0          # 가장 넓은 설정으로 고정

rows, KEEP = [], {}
for cov in [0, 2, 3, 4]:
    for st in IDX:
        out = outd/f"dimer4_cov{cov}_{st}.tsv"
        if not (out.exists() and out.stat().st_size):
            extra = f"--covered-node {cov}" if cov else ""
            sh(f'"{FD}" query -p "{PDBQ}" -q {DIMER} -i "{IDX[st]}" -t {min(THREADS,8)} '
               f'--per-structure --header --sort-by idf -d {D} -a {A} {extra} '
               f'--format-output "{FMT}" -o "{out}"', check=False)
        n = 0
        if out.exists() and out.stat().st_size:
            df = pd.read_csv(out, sep="\t")
            df.columns = [c.strip().lstrip("#") for c in df.columns]
            KEEP[(cov, st)] = df; n = len(df)
        rows.append({"covered_node": cov, "strain": st, "hits": n})
        print(f"  covered-node {cov}  {st:7s} {n}")

R = pd.DataFrame(rows).pivot(index="covered_node", columns="strain", values="hits")
print("\n=== covered-node 별 히트 수 (d=1.5, a=15) ===")
print(R.to_string())
print("  covered-node 0 = 필터 없음 (부분 매치 포함)")

# ---- 실제로 몇 잔기가 맞고 있나 ----
base = KEEP.get((0, "BL21"))
if base is not None and "node_count" in base.columns:
    print("\n=== BL21, 필터 없음 — node_count 분포 ===")
    print(base.node_count.value_counts().sort_index().to_string())
    print("  질의 잔기 4개 중 실제로 맞은 개수다. 전부 2 에 몰려 있으면")
    print("  이량체 배치가 단량체 모델에서 재현되지 않는다는 뜻이다.")

print("\n=== 판정 ===")
_b3 = int(R.loc[3, "BL21"]) if 3 in R.index else 0
_b4 = int(R.loc[4, "BL21"]) if 4 in R.index else 0
if _b4 > 0:
    print(f"4잔기 전부 맞은 구조가 {_b4}개 있다. 그 목록이 가장 좁은 후보 집합이다.")
elif _b3 > 0:
    print(f"4잔기는 0 이지만 3잔기는 {_b3}개다. 기하학적 사실 쪽이다 —")
    print("CooC1 의 이량체 Cys4 배치를 그대로 가진 단백질이 BL21 에 없다.")
    print("3잔기 목록을 후보로 쓰되, 네 번째 Cys 가 없다는 점을 명시할 것.")
else:
    print("3잔기에서도 0 이다. 질의가 인덱스와 안 맞는 쪽이다 —")
    print("단량체 예측 DB 에 이량체 계면 질의를 던진 것이 원인일 가능성이 높다.")
    print("→ 이량체 노선은 여기서 접고, 단일 사슬 Cys4 템플릿(CELL 29c 의")
    print("  QJZ12568.1 · GspE · HslO)으로 뽑은 목록을 쓴다. 그것이 지금 가진")
    print("  가장 넓고 재현 가능한 금속 자리 집합이다.")

# ---- 나온 것이 있으면 이름을 붙인다 ----
xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(xw) else {}
hdr = {}
for l in open(ASSET["faa_bl21"], errors="ignore"):
    if l.startswith(">"):
        t = l[1:].rstrip(); hdr[t.split()[0]] = t[len(t.split()[0]):].strip()
for cov in [4, 3]:
    for st in IDX:
        df = KEEP.get((cov, st))
        if df is None or not len(df): continue
        df = df.copy()
        df["stem"] = df.tid.apply(lambda x: Path(str(x)).stem)
        df["protein"] = df.stem.map(xwm).astype(str).str.split(",").str[0]
        df["desc"] = df.protein.map(lambda k: hdr.get(k, "")[:65])
        cols = [c for c in ["protein", "idf", "min_rmsd", "node_count", "plddt",
                            "matching_residues", "desc"] if c in df.columns]
        print(f"\n--- covered-node {cov} / {st} ({len(df)}개) ---")
        print(df.sort_values("min_rmsd")[cols].head(30).to_string(index=False))
        df.to_csv(DIR["table"]/f"dimer4_cov{cov}_{st}.csv", index=False, encoding="utf-8-sig")
```

---

## CELL 29j — ★ Track C 양성대조군: Y19 의 CooC 를 질의가 찾는가

```python
# =============================================================================
# CELL 29j | 빠져 있던 대조군 — Folddisco 질의가 자기 부류를 찾는가
#   Y19 는 Ni 클러스터를 가진 CODH 를 보유한다. 그러면 자기 CooC(또는 동등한 Ni
#   삽입 단백질)를 갖고 있어야 한다. 우리 질의는 CooC1 자신의 A112/A114 기하다.
#   그 질의가 Y19 의 CooC 를 못 찾는다면, 찾으려는 부류조차 못 찾는 것이다.
#
#   Track A 는 대조군 0.266 으로 실패했고 (CELL 24b),
#   Track B 는 대조군 ipTM 0.300 으로 실패했다 (CELL 28a-13).
#   Track C 만 대조군이 없었다. 여기서 세운다.
#
#   순서
#     (1) 세 균주 프로테옴에서 CooC/CooS 상동체를 서열로 찾는다 (mmseqs)
#     (2) 이름으로도 찾는다 (어노테이션에 CODH/Coo 가 있는지)
#     (3) 찾은 단백질이 folddisco 금속 목록에 있는가 — 기본/넓힌 설정 모두
#   판정
#     있다  -> 질의가 작동한다. BL21/MG1655 의 음성은 진짜 음성이다
#     없다  -> Track C 도 대조군 실패. 세 트랙 전부 자기 대조군을 못 맞힌다
# =============================================================================
need((have("ASSET", "sh"), "CELL 01 을 먼저 돌릴 것"))
FAA = {"BL21": ASSET["faa_bl21"], "Y19": ASSET["faa_y19"], "MG1655": ASSET["faa_mg1655"]}
SEQ, HDR = {}, {}
for st, fa in FAA.items():
    s_, h_, nm, buf = {}, {}, None, []
    for l in open(fa, errors="ignore"):
        if l.startswith(">"):
            if nm: s_[nm] = "".join(buf)
            t = l[1:].rstrip(); nm = t.split()[0]; buf = []
            h_[nm] = t[len(nm):].strip()
        else: buf.append(l.strip())
    if nm: s_[nm] = "".join(buf)
    SEQ[st], HDR[st] = s_, h_

# ---- (1) 서열로 찾기 ----
bf = DIR["seq"]/"baits.fasta"
need((bf.exists(), "baits.fasta 가 없다 — CELL 07 을 먼저 돌릴 것"))
baits, nm, buf = {}, None, []
for l in open(bf, errors="ignore"):
    if l.startswith(">"):
        if nm: baits[nm] = "".join(buf)
        nm, buf = l[1:].split()[0], []
    else: buf.append(l.strip())
if nm: baits[nm] = "".join(buf)
CTRLQ = {k: v for k, v in baits.items()
         if not k.startswith("seg") and (k.upper().startswith("COO") or "CODH" in k.upper())}
print("질의로 쓸 bait:", {k: len(v) for k, v in CTRLQ.items()})

qf = DIR["seq"]/"coo_controls.faa"
qf.write_text("".join(f">{k}\n{v}\n" for k, v in CTRLQ.items()))
_tmp = DIR["tmp"]/"coosearch"; _tmp.mkdir(parents=True, exist_ok=True)
FMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits"
FOUND = {}
print("\n=== (1) 서열 검색 (mmseqs, e<=1e-3) ===")
for st in FAA:
    out = DIR["search"]/f"coo_vs_{st}.m8"
    if not (out.exists() and out.stat().st_size):
        sh(f'mmseqs easy-search "{qf}" "{FAA[st]}" "{out}" "{_tmp}" '
           f'--format-output "{FMT}" -e 1e-3 -s 7.5 --threads {min(THREADS,16)} -v 1',
           check=False)
    if not (out.exists() and out.stat().st_size):
        print(f"\n[{st}] 히트 없음"); continue
    d = pd.read_csv(out, sep="\t", names=FMT.split(","))
    d = d.sort_values("bits", ascending=False).groupby("query").head(3)
    print(f"\n[{st}]")
    for r in d.itertuples(index=False):
        print(f"  {r.query:12s} -> {r.target:14s} fident {r.fident:.3f} "
              f"qcov {r.qcov:.2f} e {r.evalue:.1e}  {HDR[st].get(r.target,'')[:55]}")
        FOUND.setdefault(st, []).append(r.target)

# ---- (2) 이름으로 찾기 ----
PAT = r"carbon monoxide dehydrogenase|\bCooC\b|\bCooS\b|\bCooF\b|\bCooT\b|\bCooJ\b|CO dehydrogenase|acetyl-CoA synthase"
print("\n=== (2) 이름 검색 ===")
for st in FAA:
    hits = [(k, v) for k, v in HDR[st].items() if re.search(PAT, v, re.I)]
    print(f"\n[{st}] {len(hits)}개")
    for k, v in hits[:10]:
        print(f"   {k}  {v[:70]}")
        FOUND.setdefault(st, []).append(k)

# ---- (3) folddisco 가 그것을 잡았는가 ----
print("\n" + "=" * 100)
print("=== (3) folddisco 금속 목록에 있는가 ===")
MSET = {}
_p = DIR["table"]/"motif_metal_3strain.csv"
if _p.exists():
    m3 = pd.read_csv(_p)
    for st in FAA: MSET[st] = set(m3[m3.strain == st].key.astype(str))
    print("기본 설정 목록 크기:", {k: len(v) for k, v in MSET.items()})

# 넓힌 설정(d=1.5,a=15) 결과도 본다
xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(xw) else {}
WSET = {}
for f in (DIR["folddisco"]/"widen").glob("mono2_d1.5_a15.0_*.tsv"):
    st = f.stem.split("_")[-1]
    if not f.stat().st_size: continue
    d = pd.read_csv(f, sep="\t"); d.columns = [c.strip().lstrip("#") for c in d.columns]
    d["stem"] = d.tid.apply(lambda x: Path(str(x)).stem)
    WSET[st] = {str(k).split(",")[0] for k in d.stem.map(xwm) if pd.notna(k)}
if WSET: print("넓힌 설정(d=1.5,a=15) 목록 크기:", {k: len(v) for k, v in WSET.items()})

hit_any = False
for st in FAA:
    cands = sorted(set(FOUND.get(st, [])))
    if not cands:
        print(f"\n[{st}] CooC/CODH 상동체 자체가 없다 — 대조군으로 못 쓴다"); continue
    print(f"\n[{st}] 상동체 {len(cands)}개")
    for k in cands:
        a = "O" if k in MSET.get(st, set()) else "."
        b = "O" if k in WSET.get(st, set()) else "."
        if a == "O" or b == "O": hit_any = True
        print(f"   기본 {a} / 넓힘 {b}   {k}  {HDR[st].get(k,'')[:62]}")

print("\n" + "=" * 100)
print("=== 판정 ===")
if hit_any:
    print("질의가 CooC/CODH 계열을 잡는다. Track C 대조군 통과.")
    print("→ BL21/MG1655 에서 그런 단백질이 안 나온 것은 진짜 음성이다.")
else:
    print("⚠ 질의가 CooC/CODH 계열을 하나도 못 잡는다. Track C 도 대조군 실패다.")
    print("   Track A (0.266) · Track B (ipTM 0.300) 에 이어 세 번째다.")
    print("   즉 지금 금속 모티프 목록은 'CooC 같은 단백질'을 찾은 결과가 아니라")
    print("   'CooC1 의 두 Cys 간격과 우연히 같은 기하를 가진 단백질' 목록이다.")
    print("   후보를 이 목록에서 고르는 근거가 약해진다 — 미팅에서 이 점을 말할 것.")
    print("\n   다음 수: Y19 의 CooC 구조를 질의 템플릿으로 삼아 (CooC1 대신)")
    print("   다시 검색한다. 같은 기능을 하는 두 단백질의 기하가 얼마나 다른지가")
    print("   이 방법의 한계를 그대로 보여준다.")
```

---

## CELL 29k — Foldseek: 구조 전체로 대조군 확인 + 후보 정체 규명

```python
# =============================================================================
# CELL 29k | Folddisco 가 모티프를 본다면 Foldseek 은 구조 전체를 본다
#   두 질문에 답한다.
#
#   (A) Track C 대조군, 다른 방법으로
#       CooC1 구조로 세 균주를 훑어 CooC 상동체가 나오는가.
#       Foldseek 은 찾는데 Folddisco 는 못 찾는다 -> 문제는 생물학이 아니라
#       '두 잔기 기하' 라는 질의 설계다. 그게 지금 가장 알고 싶은 것이다.
#
#   (B) ★ QJZ12568.1 은 대체 무엇인가
#       어노테이션이 'hypothetical protein' 이라 정체가 없다. 133잔기에 Cys4
#       자리를 갖고 MG1655 프로테옴·게놈 어디에도 없는 단백질인데, 구조가
#       무엇을 닮았는지 알면 기능 가설이 선다. 미팅에서 "이게 뭐냐"에 답이 된다.
#       PDB 전체를 대상으로 하면 알려진 금속 단백질과 비교된다.
# =============================================================================
need((have("ASSET", "sh"), "CELL 01 을 먼저 돌릴 것"))
FS = shutil.which("foldseek") or globals().get("FOLDSEEK_BIN", "foldseek")
print("foldseek:", FS)
sh(f'"{FS}" version', check=False)

# 구조 디렉터리 찾기 (AlphaFold 예측 묶음)
ROOTS = sorted({Path(p).parent for r in [TOOLS/"database"/"bacteriaDB", TOOLS/"database"]
                if Path(r).exists()
                for p in glob.glob(str(Path(r)/"**"/"AF-*-F1-model_v*.cif"), recursive=True)[:1]})
STRUCT = {}
for r in [TOOLS/"database"/"bacteriaDB", TOOLS/"database"]:
    if not Path(r).exists(): continue
    for d in sorted(Path(r).glob("structures_*")):
        if d.is_dir(): STRUCT[d.name] = d
print("\n구조 디렉터리:")
for k, v in STRUCT.items():
    print(f"  {k}: {len(list(v.glob('AF-*'))[:1]) and sum(1 for _ in v.iterdir())}개  {v}")
need((bool(STRUCT), "구조 디렉터리를 못 찾았다 — 경로를 직접 지정할 것"))

fsd = DIR["folddisco"].parent/"foldseek"; fsd.mkdir(exist_ok=True)
_tmp = DIR["tmp"]/"fs"; _tmp.mkdir(parents=True, exist_ok=True)
FMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits,prob,alntmscore,lddt"

def fs_search(qpath, tgt, name, extra=""):
    out = fsd/f"{name}.m8"
    if out.exists() and out.stat().st_size:
        print(f"  [재사용] {out.name}")
    else:
        sh(f'"{FS}" easy-search "{qpath}" "{tgt}" "{out}" "{_tmp}" '
           f'--format-output "{FMT}" -e 10 --max-seqs 2000 {extra} '
           f'--threads {min(THREADS,16)}', check=False)
    if not (out.exists() and out.stat().st_size): return pd.DataFrame()
    d = pd.read_csv(out, sep="\t", names=FMT.split(","))
    d["stem"] = d.target.apply(lambda x: Path(str(x)).stem.replace(".cif", "").replace(".pdb", ""))
    return d

xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(xw) else {}
def name_of(stem):
    for k, v in xwm.items():
        if Path(str(k)).stem == stem: return str(v).split(",")[0]
    return None
HDRB = {}
for st, fa in [("BL21", ASSET["faa_bl21"]), ("Y19", ASSET["faa_y19"]),
               ("MG1655", ASSET["faa_mg1655"])]:
    for l in open(fa, errors="ignore"):
        if l.startswith(">"):
            t = l[1:].rstrip(); HDRB[t.split()[0]] = (st, t[len(t.split()[0]):].strip())

# ---------- (A) CooC1 -> 세 균주 ----------
print("\n" + "=" * 100)
print("=== (A) CooC1 구조로 세 균주 훑기 — Track C 대조군 ===")
MSET = {}
_p = DIR["table"]/"motif_metal_3strain.csv"
if _p.exists():
    m3 = pd.read_csv(_p)
    for st in ["BL21", "Y19", "MG1655"]: MSET[st] = set(m3[m3.strain == st].key.astype(str))

for sname, sdir in STRUCT.items():
    d = fs_search(ASSET["cooc1_pdb"], str(sdir), f"cooc1_vs_{sname}")
    if not len(d):
        print(f"\n[{sname}] 히트 없음"); continue
    d["protein"] = d.stem.map(name_of)
    d["desc"] = d.protein.map(lambda k: HDRB.get(str(k), ("", ""))[1][:58] if k else "")
    d = d.sort_values("bits", ascending=False).head(12)
    print(f"\n[{sname}] 상위 12")
    for r in d.itertuples(index=False):
        st = HDRB.get(str(r.protein), ("?",))[0]
        inm = "O" if (r.protein and r.protein in MSET.get(st, set())) else "."
        print(f"  모티프{inm}  {str(r.protein):14s} tm {r.alntmscore:.3f} "
              f"fident {r.fident:.3f} e {r.evalue:.1e}  {r.desc}")
print("\n모티프 열이 'O' 면 folddisco 도 잡은 것이다. 구조적으로 CooC 를 닮은 것이")
print("나오는데 전부 '.' 이면, 문제는 생물학이 아니라 '두 잔기 기하' 라는 질의 설계다.")

# ---------- (B) 후보의 정체 ----------
print("\n" + "=" * 100)
print("=== (B) 후보 단백질이 무엇을 닮았는가 ===")
TARGETS = ["QJZ12568.1", "QJZ11471.1", "QJZ13396.1"]
def find_struct(pid):
    for k, v in xwm.items():
        if str(v).split(",")[0] == pid:
            stem = Path(str(k)).stem
            for sdir in STRUCT.values():
                for ext in ("cif", "pdb", "cif.gz", "pdb.gz"):
                    f = next(iter(glob.glob(str(sdir/f"{stem}*.{ext}"))), None)
                    if f: return f
    return None

# PDB 데이터베이스 (알려진 구조와 비교). 몇 GB 라 기본은 끔.
USE_PDB = False          # True 로 바꾸면 foldseek 이 PDB DB 를 내려받는다 (~2GB)
PDBDB = fsd/"pdb"
if USE_PDB and not Path(str(PDBDB) + ".dbtype").exists():
    print("PDB DB 내려받는 중 (수 분 ~ 수십 분)...")
    sh(f'"{FS}" databases PDB "{PDBDB}" "{_tmp}" --threads {min(THREADS,16)}', check=False)

for pid in TARGETS:
    sp = find_struct(pid)
    print(f"\n--- {pid}  {HDRB.get(pid,('','?'))[1][:60]}")
    if not sp:
        print("    구조 파일을 못 찾았다"); continue
    print(f"    구조 {sp}")
    for sname, sdir in STRUCT.items():
        d = fs_search(sp, str(sdir), f"{pid}_vs_{sname}")
        if not len(d): continue
        d["protein"] = d.stem.map(name_of)
        d = d[d.protein.astype(str) != pid].sort_values("bits", ascending=False).head(5)
        if not len(d): continue
        print(f"    [{sname}]")
        for r in d.itertuples(index=False):
            print(f"      {str(r.protein):14s} tm {r.alntmscore:.3f} fident {r.fident:.3f}"
                  f"  {HDRB.get(str(r.protein), ('',''))[1][:52]}")
    if USE_PDB and Path(str(PDBDB) + ".dbtype").exists():
        d = fs_search(sp, str(PDBDB), f"{pid}_vs_PDB")
        if len(d):
            print("    [PDB]")
            for r in d.sort_values("bits", ascending=False).head(8).itertuples(index=False):
                print(f"      {r.target:20s} tm {r.alntmscore:.3f} prob {r.prob:.3f} "
                      f"e {r.evalue:.1e}")
print("\n※ alntmscore 0.5 이상이면 같은 폴드로 본다. QJZ12568.1 이 알려진 금속")
print("  샤페론·아연 손가락 계열과 닮으면 'hypothetical' 이라는 이름이 걷힌다.")
print("  PDB 비교를 하려면 위 USE_PDB 를 True 로 두고 다시 돌릴 것 (~2GB 다운로드).")
```

---

## CELL 30 — Part 8b. ID crosswalk (구조 tid → GenBank protein ID)

```python
# =============================================================================
# CELL 30 | Part 8-2. ★ID crosswalk — Folddisco tid 를 GenBank ID 에 붙인다
#   문제: 구조 파일명이 AF-<UniProtID>-F1-model_v6.cif / cf_<UniParcID>.pdb 인데
#         우리 분류·프로테옴은 GenBank ID (QJZ…, AKE…) 라 그냥은 조인이 안 된다.
#   1순위: 구조 DB 를 만든 쪽이 준 대응표 structure_accessions.xlsx
#          - 시트당 UniParc / UniProtKB / Genbank / Selected structure / AFDB index
#          - "Selected structure" 가 AF 파일명의 accession. AFDB 에 없어 ColabFold 로
#            직접 예측한 것은 그 칸이 "ColabFold predicted" 이고 UniParc 로 식별된다.
#          - 실측: BL21 3,846 AF + 250 CF = 4,096 / Y19 5,262 AF + 21 CF = 5,283.
#            키가 전부 1:1 이라 손실이 없다.
#   2순위(대응표가 없을 때): 구조 파일에서 CA 서열을 뽑아 프로테옴과 완전일치 매칭.
#          정확하지만 9,000여 파일을 읽어 수 분 걸리고, isoform/부분서열은 놓친다.
# =============================================================================
XWALK_CSV = DIR["table"]/"id_crosswalk_struct_to_genbank.csv"
import pandas as pd

# 대응표는 워크스페이스 안(external/)에 둘 수도, 구조 DB 옆에 둘 수도 있다.
# 파일명 접두사도 제각각이라 두 자리를 다 훑는다.
_XLSX_DIRS = [DIR["external"], TOOLS/"database"/"bacteriaDB", TOOLS/"database"]
XLSX = next((f for d in _XLSX_DIRS if d.exists()
             for f in sorted(d.glob("*structure_accessions*.xlsx"))), None)
print("대응표:", XLSX if XLSX else f"없음 (찾은 자리: {[str(d) for d in _XLSX_DIRS]})")

def xwalk_from_xlsx(path):
    """대응표에서 accession -> GenBank 사전을 만든다 (구조 파일을 안 읽는다)."""
    frames = []
    for sheet in pd.ExcelFile(path).sheet_names:
        d = pd.read_excel(path, sheet_name=sheet)
        tag = ("BL21" if "UP000503272" in sheet else
               "Y19"  if "UP000034085" in sheet else sheet.strip())
        sel   = d["Selected structure"].astype(str).str.strip()
        is_cf = sel.str.contains("ColabFold", case=False, na=False)
        frames.append(pd.DataFrame({
            "db":       tag,
            "acc_type": is_cf.map({True: "UniParc", False: "UniProt"}),
            # ColabFold 로 예측한 것은 Selected structure 가 accession 이 아니라
            # 안내 문구다. 그 행은 UniParc 로 파일명(cf_<UniParc>.pdb)과 맞춘다.
            "acc":      sel.where(~is_cf, d["UniParc accession"].astype(str).str.strip()),
            "protein":  d["Genbank"].astype(str).str.strip(),
            "matched_strain": {"BL21": "BL21DE3"}.get(tag, tag),
        }))
    return pd.concat(frames, ignore_index=True)

def tid_to_acc(tid):
    """구조 파일 stem 에서 accession 을 뽑는다."""
    m = re.match(r"AF-([A-Z0-9]+)-F\d+", tid)
    if m:              return "UniProt", m.group(1)
    if tid.startswith("cf_"): return "UniParc", tid[3:]
    return "unknown", tid

if XWALK_CSV.exists():
    XW = pd.read_csv(XWALK_CSV)
    print("캐시 로드:", XWALK_CSV, len(XW), "행")

elif XLSX is not None:
    acc2gb = {}
    for r in xwalk_from_xlsx(XLSX).itertuples(index=False):
        acc2gb[r.acc] = (r.protein, r.matched_strain)
    print(f"대응표 로드: {XLSX.name}  accession {len(acc2gb)}개")

    rows = []
    for tag, d in [("BL21", ASSET["struct_bl21"]), ("Y19", ASSET["struct_y19"])]:
        if not d.exists():
            print("[없음]", d); continue
        files = sorted(list(d.glob("*.cif")) + list(d.glob("*.pdb")))
        for f in files:                       # 파일명만 본다 — 내용을 안 읽으니 즉시 끝난다
            acc_type, acc = tid_to_acc(f.stem)
            gb, strain = acc2gb.get(acc, (None, None))
            rows.append({"tid": f.stem, "db": tag, "acc_type": acc_type, "acc": acc,
                         "protein": gb, "matched_strain": strain})
        print(f"  {tag}: 구조 {len(files)}개")
    XW = pd.DataFrame(rows)
    XW.to_csv(XWALK_CSV, index=False)
    print("저장:", XWALK_CSV)

else:
    print("[대응표 없음]")
    print(f"  협업팀이 준 *structure_accessions*.xlsx 를 {DIR['external']} 에 두면 즉시 끝난다.")
    print("  없으면 구조 서열을 직접 읽어 매칭한다 (9,000여 파일, 수 분).")
    seq2gb = {}
    for strain in ["BL21DE3", "Y19", "MG1655"]:
        for pid, sq in PROTEOMES.get(strain, {}).get("seqs", {}).items():
            seq2gb.setdefault(sq.rstrip("*"), (pid, strain))
    # 이 사전은 CELL 08 이 채운다. 비어 있는데 그냥 진행하면 9,000여 파일을 수 분간
    # 읽고도 매칭 0건이 나오고, 그 0건이 캐시로 저장돼 이후 셀이 전부 빈손이 된다.
    # 실측으로 한 번 당했다 (CELL 01 -> 29 -> 30 으로 건너뛴 경우).
    need((bool(seq2gb), "CELL 08 을 먼저 돌릴 것 — 비교할 프로테옴 서열이 메모리에 없다"))
    print(f"비교 대상 서열 {len(seq2gb)}개")

    def struct_seq(path):
        """cif/pdb 에서 첫 체인의 CA 기준 서열."""
        txt = Path(path).read_text(errors="ignore")
        if path.suffix == ".cif":
            seq, seen = [], set()
            for line in txt.splitlines():
                if not line.startswith("ATOM"): continue
                f = line.split()
                if len(f) < 9 or f[3] != "CA": continue
                key = (f[6], f[8]) if len(f) > 8 else (f[6],)
                if key in seen: continue
                seen.add(key); seq.append(AA3to1.get(f[5], "X"))
            return "".join(seq)
        chains, _ = pdb_chain_seqs(path)
        return max(chains.values(), key=len) if chains else ""

    rows = []
    for tag, d in [("BL21", ASSET["struct_bl21"]), ("Y19", ASSET["struct_y19"])]:
        if not d.exists():
            print("[없음]", d); continue
        files = sorted(list(d.glob("*.cif")) + list(d.glob("*.pdb")))
        print(f"{tag}: {len(files)} 구조 파싱 중...")
        for i, f in enumerate(files):
            if i and i % 1000 == 0: print(f"  {i}/{len(files)}")
            acc_type, acc = tid_to_acc(f.stem)
            sq = struct_seq(f)
            gb, strain = seq2gb.get(sq, (None, None))
            rows.append({"tid": f.stem, "db": tag, "acc_type": acc_type, "acc": acc,
                         "struct_len": len(sq), "protein": gb, "matched_strain": strain})
    XW = pd.DataFrame(rows)
    XW.to_csv(XWALK_CSV, index=False)
    print("저장:", XWALK_CSV)

if len(XW):
    hit = XW.protein.notna().mean()
    print(f"\n매칭률 {hit:.1%}  ({XW.protein.notna().sum()}/{len(XW)})")
    print(XW.groupby(["db","acc_type"]).agg(n=("tid","size"),
                                            matched=("protein", lambda x: x.notna().sum())).to_string())
    if hit < 0.99:
        print("\n미매칭이 남으면: 구조가 isoform/부분서열이거나 프로테옴 버전이 다른 것.")
        print("  -> UniProt idmapping API (UniProtKB_AC-ID -> EMBL-GenBank-DDBJ_CDS) 로 보완 가능")
    print("\n※ MG1655 는 공식 배포 인덱스라 구조 파일이 우리 DB 에 없다 -> crosswalk 대상이 아니다.")
    print("   CELL 31 에서 MG1655 매칭 0건은 정상이며, 비교용으로만 본다.")
```

---

## CELL 31 — Part 8c. Folddisco 결과 → motif CSV 변환

```python
# =============================================================================
# CELL 31 | Part 8-3. Folddisco TSV -> Part 8 이 기대하는 (protein, score) CSV 로 변환
#   - score 는 idf 를 기본으로 쓰고, min_rmsd 는 낮을수록 좋으므로 역수로 보정 결합
# =============================================================================
# FD_METAL / FD_ATP 는 CELL 29 가 메모리에 올린다. tsv 파일이 디스크에 남아 있어도
# 커널을 새로 띄웠으면 변수는 없다 — 그냥 두면 NameError 로 끝난다.
need((have("XW"), "CELL 30 을 먼저 돌릴 것"),
     (have("FD_METAL", "FD_ATP"), "CELL 29 를 먼저 돌릴 것 (tsv 가 있어도 커널에 변수가 없다)"))

xw_map = XW.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict() if len(XW) else {}

def fd_to_csv(tsv_paths, out_csv, score_col="idf"):
    frames = []
    for strain, p in (tsv_paths or {}).items():
        p = Path(p)
        if not p.exists() or p.stat().st_size == 0:
            print(f"[스킵] {p}"); continue
        d = pd.read_csv(p, sep="\t")
        d.columns = [c.strip().lstrip("#") for c in d.columns]
        if "tid" not in d.columns:
            print(f"[형식 확인 필요] {p}: {list(d.columns)[:6]}"); continue
        d["tid_stem"] = d.tid.apply(lambda x: Path(str(x)).stem)
        d["protein"] = d.tid_stem.map(xw_map)
        sc = d[score_col].astype(float) if score_col in d.columns else 0.0
        if "min_rmsd" in d.columns:
            sc = sc * (1.0 / (1.0 + d["min_rmsd"].astype(float)))
        d["score"] = sc
        d["strain"] = strain
        # nres/plddt 는 여기서 거르지 않고 그대로 들고 간다. Folddisco 는 --sort-by idf
        # 로 "희귀한" 기하를 위로 올리는데, AlphaFold DB 에서 희귀한 기하는 드문 기능이
        # 아니라 엉망으로 예측된 단편인 경우가 많다 (실측: idf 1위가 nres=44, plddt=56.7).
        # 판단은 CELL 32 통합 단계에서 하되, 판단 근거는 detail CSV 에 남겨둔다.
        frames.append(d[["protein","score","strain","tid","tid_stem"] +
                        [c for c in ["idf","min_rmsd","nres","plddt","matching_residues"]
                         if c in d.columns]])
    if not frames:
        print("변환할 결과 없음"); return pd.DataFrame(columns=["protein","score"])
    A = pd.concat(frames)
    rep = A.groupby("strain").agg(hits=("tid", "size"),
                                  mapped=("protein", lambda x: x.notna().sum()))
    print(rep.to_string())
    if "MG1655" in rep.index and rep.loc["MG1655", "mapped"] == 0:
        print("  (MG1655 0건 매칭은 정상 — 구조 DB 가 없어 crosswalk 불가. 비교용으로만 본다)")
    out = (A.dropna(subset=["protein"]).groupby("protein", as_index=False)
             .score.max().sort_values("score", ascending=False))
    out.to_csv(out_csv, index=False)
    A.to_csv(str(out_csv).replace(".csv", "_detail.csv"), index=False)
    print(f"{out_csv.name}: {len(out)}개 단백질 (미매칭 {int(A.protein.isna().sum())}행 제외)")
    return out

FOLDDISCO_NI  = DIR["table"]/"folddisco_metal_motif.csv"
FOLDDISCO_ATP = DIR["table"]/"folddisco_atp_motif.csv"
ni_df  = fd_to_csv(FD_METAL, FOLDDISCO_NI)
atp_df = fd_to_csv(FD_ATP,   FOLDDISCO_ATP)
print("\n상위 10 (metal motif):")
print(ni_df.head(10).to_string(index=False))
```

---

## CELL 32 — Part 8d. 통합 랭킹

```python
# =============================================================================
# CELL 32 | Part 8-4. Track A + B + C 통합 랭킹
#   가중치 근거:
#     strain_specific 0.30 — 실험(BL21 lysate 만 rescue)과 직결. 최대 가중치
#     Ni_motif 0.20 / ATP_motif 0.15 — Folddisco. nucleotide 의존성 실험으로 확정 가능
#     RF2PPI 0.20 / Boltz ipTM 0.15 — 결합 자체의 증거
#   ★한 방법에서만 높은 후보보다 독립적인 두 방법에서 동시에 높은 후보가 훨씬 신뢰도 높다
#
#   ⚠ 이 셀의 가중치는 검증 전에 정한 것이다. 이후 확인된 사실:
#       rf2 0.20  — Track A 는 양성대조군 0.266 으로 실패했다 (CELL 24b). 쓸 수 없다.
#       boltz 0.15 — ipTM 눈금이 이 계에서 검증되지 않았다 (CELL 28a-13 대기).
#     즉 priority 의 35% 가 못 믿을 값에서 나온다. 후보 지명에는 CELL 36 의
#     FOLDDISCO_candidates.csv 를 쓰고, 이 표는 '전체 계산 기록'으로만 남긴다.
#   선행: CELL 12(TRACK_A/B) · 23(R_best) · 28(B) · 31(ni_df/atp_df) 가 메모리에 있어야 한다.
#     커널을 새로 띄웠다면 변수가 전부 사라졌으므로 처음부터 다시 돌려야 한다.
# =============================================================================
need((have("TRACK_A", "TRACK_B", "cls_prey", "hdr2desc"),
      "CELL 12 까지 다시 돌릴 것 (커널을 새로 띄우면 TRACK_A/B 가 사라진다)"),
     (have("ni_df", "atp_df"), "CELL 31 을 먼저 돌릴 것"))
def safe_series(df, key, col):
    try:
        return df.set_index(key)[col]
    except Exception:
        return pd.Series(dtype=float)

M = pd.DataFrame(index=sorted(set(TRACK_A) | set(TRACK_B)))
M["category"]      = cls_prey.set_index("protein")["category"]
M["desc"]          = [hdr2desc.get(i, "")[:70] for i in M.index]
M["rf2ppi"]        = R_best["mean"]        if "R_best" in dir() else np.nan
M["rf2ppi_sd"]     = R_best["sd"]          if "R_best" in dir() else np.nan
M["paired_depth"]  = R_best["paired_depth"] if "R_best" in dir() else np.nan
M["boltz_iptm"]    = safe_series(B, "prey", "iptm")        if "B" in dir() and len(B) else np.nan
M["boltz_ligiptm"] = safe_series(B, "prey", "ligand_iptm") if "B" in dir() and len(B) else np.nan
M["ni_motif"]      = safe_series(ni_df, "protein", "score")  if len(ni_df)  else np.nan
M["atp_motif"]     = safe_series(atp_df, "protein", "score") if len(atp_df) else np.nan

def nz(s):
    s = s.astype(float)
    if s.notna().sum() == 0 or s.max() == s.min():
        return s.fillna(0) * 0
    return ((s - s.min()) / (s.max() - s.min())).fillna(0)

W = {"strain": 0.30, "ni": 0.20, "atp": 0.15, "rf2": 0.20, "boltz": 0.15}
M["strain_flag"] = M.category.isin(["strain-specific", "low-similarity"]).astype(float)
M["priority"] = (W["strain"]*M.strain_flag + W["ni"]*nz(M.ni_motif) + W["atp"]*nz(M.atp_motif)
                 + W["rf2"]*nz(M.rf2ppi) + W["boltz"]*nz(M.boltz_iptm))
# Foldseek-Interface (CELL 28b-c). 가중치에는 넣지 않고 참고 열로만 둔다 —
# 계면 전례 유무는 결합의 증거라기보다 해석의 근거다.
_ifp = DIR["table"]/"trackB_interface_match.csv"
M["iface_tm"] = (pd.read_csv(_ifp).set_index("prey")["iface_tm"]
                 if _ifp.exists() else np.nan)

M["n_evidence"] = ((M.rf2ppi > 0.3).fillna(False).astype(int)
                   + (M.boltz_iptm > 0.6).fillna(False).astype(int)
                   + M.ni_motif.notna().astype(int) + M.atp_motif.notna().astype(int))
M = M.sort_values(["priority", "n_evidence"], ascending=False)
M.to_csv(DIR["table"]/"INTEGRATED_candidate_ranking.csv")
print(M.head(50).to_string())
print("\n독립적 2개 이상 증거를 가진 후보:", int((M.n_evidence >= 2).sum()))
```

---

## CELL 33 — Folddisco 노선 (0) 준비: 정체 확인 + 균주 간 대응표

```python
# =============================================================================
# CELL 33 | Folddisco 우선 노선 — 0단계. 누가 누구인지부터 확정한다
#   ID 체계가 균주마다 다르다.
#     BL21 / Y19  GenBank protein ID (QJZ*, AKE*) — crosswalk 으로 붙어 있음
#     MG1655      UniProt accession — 구조 DB 가 공식 배포본이라 crosswalk 없음
#
#   ★ 균주 간 대응은 metal 모티프 단백질에만 건다.
#     ATP(Walker A) 는 균주당 1,000개 안팎이 걸린다 (프로테옴의 20% 이상). 이걸
#     전부 대응시키면 계산이 몇 배로 늘어나는데, ATP 는 3단계에서 '같은 균주 안에서
#     금속과 겹치는가'로만 쓰므로 균주 간 대응이 애초에 필요 없다.
#
#   대응표는 두 개를 따로 만든다. 섞으면 '검출 차이'와 '유전자 부재'가 구분되지 않는다.
#     (2) 모티프 집합끼리        -> 다른 균주에서도 folddisco 가 잡았나
#     (3) 모티프 집합 vs 프로테옴 -> 그 균주에 유전자가 있기는 한가
# =============================================================================
need(((DIR["table"]/"folddisco_metal_motif_detail.csv").exists(),
      "CELL 29 -> 30 -> 31 을 먼저 돌릴 것"))
import urllib.request
import numpy as np

# ---------- 프로테옴 먼저 (ID 검증에 쓴다) ----------
def read_fasta(path, with_desc=False):
    seq, desc, nm, buf = {}, {}, None, []
    for l in open(path, errors="ignore"):
        if l.startswith(">"):
            if nm: seq[nm] = "".join(buf)
            h = l[1:].rstrip(); nm = h.split()[0]; buf = []
            desc[nm] = h[len(nm):].strip()
        else: buf.append(l.strip())
    if nm: seq[nm] = "".join(buf)
    return (seq, desc) if with_desc else seq

FAA  = {"BL21": ASSET["faa_bl21"], "Y19": ASSET["faa_y19"], "MG1655": ASSET["faa_mg1655"]}
PROT, DESC, GENE = {}, {}, {}
for st, fa in FAA.items():
    s, d = read_fasta(fa, with_desc=True)
    PROT[st] = s
    for k, v in d.items(): DESC[(st, k)] = v
    print(f"{st:7s} 프로테옴 {len(s):5d}개  ({Path(fa).name})")

# ---------- folddisco 결과 읽기 ----------
FDD = {}
for tag, fn in [("metal", "folddisco_metal_motif_detail.csv"),
                ("atp",   "folddisco_atp_motif_detail.csv")]:
    p = DIR["table"]/fn
    if not p.exists():
        print(f"[없음] {fn}"); continue
    d = pd.read_csv(p); d["motif"] = tag
    FDD[tag] = d
    print(f"{tag:6s} {len(d):5d}행   {d.strain.value_counts().to_dict()}")
need((bool(FDD), "folddisco detail CSV 가 하나도 없다"))

# ---------- 정체(key) 확정 ----------
# crosswalk 이 tid 하나에 단백질 여러 개를 콤마로 붙여 준 행이 있다 -> 첫 번째만 쓴다.
# 그리고 균주가 뒤섞인 행이 있다 (BL21 행에 Y19 의 AKE* 가 붙어 있었다). 자기 균주
# 프로테옴에 그 ID 가 실제로 있는지로 검증해 걸러낸다.
UP_RE = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")

def resolve(d):
    d = d.copy()
    d["uniprot"] = d.tid_stem.astype(str).str.extract(r"AF-([A-Z0-9]+)-F", expand=False)
    k = d.protein.where(d.protein.notna(), d.uniprot).astype(str)
    d["key"] = k.str.split(",").str[0].str.strip()
    ok = []
    for st, key in zip(d.strain, d.key):
        if st == "MG1655": ok.append(bool(UP_RE.match(str(key))))
        else:              ok.append(str(key) in PROT.get(st, {}))
    d["key_ok"] = ok
    return d

RAW = {t: resolve(v) for t, v in FDD.items()}
for t, d in RAW.items():
    bad = d[~d.key_ok]
    print(f"\n{t}: 정체 확정 {int(d.key_ok.sum())}행 / 버림 {len(bad)}행")
    if len(bad):
        print("  버린 예:", list(zip(bad.strain, bad.key))[:5])
        print("  (crosswalk 이 균주를 넘나들거나 tid 하나에 단백질 여럿을 붙인 행이다)")
MET_RAW = RAW["metal"][RAW["metal"].key_ok].copy()
ATP_RAW = RAW["atp"][RAW["atp"].key_ok].copy() if "atp" in RAW else pd.DataFrame()
print(f"\nmetal 단백질 {MET_RAW.groupby(['strain','key']).ngroups}개 / "
      f"atp 단백질 {ATP_RAW.groupby(['strain','key']).ngroups if len(ATP_RAW) else 0}개")

# ---------- MG1655 이름·서열 (metal 만, 묶어서 조회) ----------
MGSEQ  = {}
_acc   = sorted(MET_RAW[MET_RAW.strain == "MG1655"].key.unique())
_cache = DIR["table"]/"mg1655_metal_uniprot.tsv"

def uniprot_fetch(accs, chunk=90):
    """URL 길이 제한 때문에 묶어서 받는다. 한 번에 1,000개를 넣으면 HTTP 400 이다."""
    hdr, body = None, []
    for i in range(0, len(accs), chunk):
        url = ("https://rest.uniprot.org/uniprotkb/accessions?accessions="
               + ",".join(accs[i:i+chunk])
               + "&format=tsv&fields=accession,gene_primary,protein_name,sequence")
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                ls = r.read().decode().splitlines()
        except Exception as e:
            print(f"  ⚠ {i}~{i+chunk} 실패: {e}"); continue
        if not ls: continue
        if hdr is None: hdr = ls[0]
        body += ls[1:]
    return hdr, body

if _acc and not _cache.exists():
    _h, _b = uniprot_fetch(_acc)
    if _h: _cache.write_text("\n".join([_h] + _b))
    print(f"UniProt: {len(_b)}/{len(_acc)}개 조회 -> {_cache.name}")
if _cache.exists():
    _u = pd.read_csv(_cache, sep="\t")
    for r in _u.itertuples(index=False):
        acc = str(r[0])
        GENE[("MG1655", acc)] = str(r[1]) if len(r) > 1 and pd.notna(r[1]) else ""
        DESC[("MG1655", acc)] = str(r[2]) if len(r) > 2 and pd.notna(r[2]) else ""
        if len(r) > 3 and isinstance(r[3], str): MGSEQ[acc] = r[3]
    print(f"MG1655 이름 {len(_u)}개 / 서열 {len(MGSEQ)}개")

# BL21/Y19 유전자명은 설명문의 관용 표기(예: FlgH, GmhB)에서 뽑는다
for (st, k), d in list(DESC.items()):
    if st != "MG1655":
        m = re.search(r"\b([a-z]{3}[A-Z])\b", str(d))
        GENE[(st, k)] = m.group(1) if m else ""

# ---------- metal 단백질 서열 -> mmseqs 대응표 ----------
QF, missing = DIR["seq"]/"folddisco_metal_proteins.faa", []
with open(QF, "w") as fh:
    for st, k in sorted(set(zip(MET_RAW.strain, MET_RAW.key))):
        sq = MGSEQ.get(k) if st == "MG1655" else PROT.get(st, {}).get(k)
        if not sq: missing.append((st, k)); continue
        fh.write(f">{st}|{k}\n{sq}\n")
_n = sum(1 for l in open(QF) if l.startswith(">"))
print(f"\nmetal 단백질 서열 {_n}개 -> {QF.name}  (못 찾음 {len(missing)}개 {missing[:5]})")
need((_n > 0, "서열을 하나도 못 모았다"))

_tmp = DIR["tmp"]/"fdmap"; _tmp.mkdir(parents=True, exist_ok=True)
FMT  = "query,target,fident,alnlen,qcov,tcov,evalue,bits"
def _search(target_fa, out_name):
    out = DIR["search"]/out_name
    if out.exists() and out.stat().st_size:
        print(f"  [재사용] {out.name}"); return out
    sh(f'mmseqs easy-search "{QF}" "{target_fa}" "{out}" "{_tmp}" '
       f'--format-output "{FMT}" -e 1e-5 --threads {min(THREADS,16)} -v 1', check=False)
    return out

print("\n(2) 모티프 집합끼리 — 다른 균주에서도 folddisco 가 잡았나")
_self = _search(QF, "fdm_motif_vs_motif.m8")
print("(3) 전체 프로테옴 — 그 균주에 유전자가 있기는 한가")
_full = {st: _search(fa, f"fdm_motif_vs_{st}.m8") for st, fa in FAA.items()}

COLS = FMT.split(",")
def _rd(m8):
    p = Path(m8)
    if not p.exists() or not p.stat().st_size: return pd.DataFrame(columns=COLS + ["q_strain","q_key"])
    d = pd.read_csv(p, sep="\t", names=COLS)
    d[["q_strain", "q_key"]] = d["query"].str.split("|", n=1, expand=True)
    return d
ORTH_SELF = _rd(_self)
ORTH_FULL = {st: _rd(v) for st, v in _full.items()}
print(f"\n(2) {len(ORTH_SELF):,}행   (3) " + " / ".join(f"{k} {len(v):,}" for k, v in ORTH_FULL.items()))
print("\n다음: CELL 34")
```

---

## CELL 34 — Folddisco 노선 (1) metal 모티프 3자 비교

```python
# =============================================================================
# CELL 34 | 1단계. 금속 모티프를 두 축으로 분류한다
#     축1  모티프 검출: folddisco 가 그 균주에서도 이 자리를 잡았는가
#     축2  유전자 존재: 애초에 그 균주 프로테옴에 이 단백질이 있는가
#   판정:
#     BL21 검출 + MG1655 유전자 없음 -> 진짜 균주 특이 (강한 후보)
#     BL21 검출 + MG1655 유전자 있음 -> 검출 차이일 뿐 (약한 후보, 버리지 않음)
#   지난번 'BL21 특이 14개'가 사실 대장균 공통 유전자였던 건 축2 를 안 봤기 때문이다.
# =============================================================================
need((have("MET_RAW", "ORTH_SELF", "ORTH_FULL", "DESC", "GENE"), "CELL 33 를 먼저 돌릴 것"))
STRAINS = ["BL21", "MG1655", "Y19"]
ID_THR, COV_THR = 0.50, 0.70

def parse_match(s):
    """'A27,A30:0.5054;A50,A52:0.81' -> (가장 잘 맞은 잔기, 그 RMSD, 매치 수)"""
    if not isinstance(s, str) or not s.strip(): return ("", float("nan"), 0)
    items = []
    for part in s.split(";"):
        part = part.strip()
        if not part: continue
        res, sep, rm = part.rpartition(":")
        try: items.append((res if sep else part, float(rm)))
        except ValueError: items.append((part, float("nan")))
    if not items: return ("", float("nan"), 0)
    ok = [i for i in items if i[1] == i[1]]
    b = min(ok, key=lambda x: x[1]) if ok else items[0]
    return (b[0], b[1], len(items))

_s = ORTH_SELF[(ORTH_SELF.fident >= ID_THR) & (ORTH_SELF.qcov >= COV_THR)].copy()
MOTIF_HIT = {}
if len(_s):
    _s[["t_strain", "t_key"]] = _s.target.str.split("|", n=1, expand=True)
    for r in _s.itertuples(index=False):
        d = MOTIF_HIT.setdefault((r.q_strain, r.q_key), {})
        if r.fident > d.get(r.t_strain, (0, ""))[0]: d[r.t_strain] = (r.fident, r.t_key)

GENE_HIT = {}
for st, df in ORTH_FULL.items():
    if not len(df): continue
    for r in df[(df.fident >= ID_THR) & (df.qcov >= COV_THR)].itertuples(index=False):
        d = GENE_HIT.setdefault((r.q_strain, r.q_key), {})
        if r.fident > d.get(st, (0, ""))[0]: d[st] = (r.fident, r.target)

def motif_summary(d):
    """(strain, key) 단위로 가장 잘 맞은 매치 한 줄만 남긴다."""
    if not len(d): return pd.DataFrame()
    d = d.copy()
    if "matching_residues" in d.columns:
        pm = d.matching_residues.apply(parse_match)
        d["residues"] = [x[0] for x in pm]
        d["rmsd"]     = [x[1] for x in pm]
        d["n_match"]  = [x[2] for x in pm]
    else:
        d["residues"], d["rmsd"], d["n_match"] = "", float("nan"), 0
    if "min_rmsd" in d.columns:
        d["rmsd"] = d["rmsd"].fillna(d["min_rmsd"])
    d = d.sort_values(["rmsd", "idf"], ascending=[True, False])
    g = d.groupby(["strain", "key"], as_index=False).first()
    keep = ["strain", "key", "idf", "rmsd", "residues", "n_match", "nres", "plddt", "tid_stem"]
    return g[[c for c in keep if c in g.columns]]

def annotate(g):
    g = g.copy()
    g["gene"] = [GENE.get((s, k), "") for s, k in zip(g.strain, g.key)]
    g["description"] = [str(DESC.get((s, k), ""))[:90] for s, k in zip(g.strain, g.key)]
    for st in STRAINS:
        g[f"motif_in_{st}"] = [int(st in MOTIF_HIT.get((s, k), {})) for s, k in zip(g.strain, g.key)]
        g[f"gene_in_{st}"]  = [int(st in GENE_HIT.get((s, k), {}))  for s, k in zip(g.strain, g.key)]
        g.loc[g.strain == st, [f"motif_in_{st}", f"gene_in_{st}"]] = 1
    g["motif_pattern"] = ["+".join(s for s in STRAINS if r[f"motif_in_{s}"]) for _, r in g.iterrows()]
    g["gene_pattern"]  = ["+".join(s for s in STRAINS if r[f"gene_in_{s}"])  for _, r in g.iterrows()]
    g["annotated"] = [0 if re.search(r"hypothetical|uncharacteri|DUF\d|putative protein",
                                     str(x), re.I) else 1 for x in g.description]
    def verdict(r):
        if r.strain != "BL21": return "-"
        if not r.motif_in_MG1655 and not r.gene_in_MG1655: return "BL21 특이 (유전자 없음)"
        if not r.motif_in_MG1655 and r.gene_in_MG1655:     return "검출 차이 (유전자는 있음)"
        return "공통"
    g["strain_verdict"] = [verdict(r) for _, r in g.iterrows()]
    return g

MET = annotate(motif_summary(MET_RAW))
need((len(MET) > 0, "metal 요약이 비었다"))
MET.to_csv(DIR["table"]/"motif_metal_3strain.csv", index=False, encoding="utf-8-sig")

print("=== 균주별 metal 모티프 단백질 수 ===")
print(MET.strain.value_counts().to_string())
print("\n=== 검출 패턴 (folddisco 가 잡은 균주) ===")
print(MET.motif_pattern.value_counts().to_string())
print("\n=== 유전자 존재 패턴 (프로테옴에 있는 균주) ===")
print(MET.gene_pattern.value_counts().to_string())
print("\n=== BL21 판정 ===")
print(MET[MET.strain == "BL21"].strain_verdict.value_counts().to_string())
_hot = MET[(MET.strain == "BL21") & MET.strain_verdict.str.startswith("BL21 특이")]
print(f"\n--- BL21 특이 (MG1655 에 유전자 자체가 없음): {len(_hot)}개 ---")
if len(_hot):
    print(_hot[["key", "gene", "rmsd", "idf", "residues", "plddt", "description"]]
          .sort_values("rmsd").to_string(index=False))
else:
    print("  없음 — 모든 BL21 금속 모티프 단백질이 MG1655 에도 유전자가 있다는 뜻이다.")
    print("  그러면 BL21/MG1655 차이는 유전자 유무가 아니라 발현·안정성·구조 예측")
    print("  품질 쪽이라는 것이고, '검출 차이' 목록을 후보로 봐야 한다.")
print(f"\n저장: {DIR['table']/'motif_metal_3strain.csv'}")
```

---

## CELL 35 — Folddisco 노선 (2) ATP 모티프

```python
# =============================================================================
# CELL 35 | 2단계. ATP(Walker A) 모티프 — 필터가 아니라 표시로 쓴다
#   실측: 균주당 1,000개 안팎이 걸린다. 프로테옴의 20% 이상이다.
#   이 해상도로는 후보를 좁힐 수 없다. 그래서 균주 간 대응은 만들지 않고
#   '같은 균주 안에서 금속 모티프와 겹치는가'(3단계)에만 쓴다.
# =============================================================================
need((have("ATP_RAW", "motif_summary"), "CELL 33 -> 34 를 먼저 돌릴 것"))
ATP = motif_summary(ATP_RAW)
if not len(ATP):
    print("ATP 결과 없음 — CELL 29 의 RESIDUES_ATP 와 CELL 31 변환을 확인할 것.")
else:
    ATP["gene"] = [GENE.get((s, k), "") for s, k in zip(ATP.strain, ATP.key)]
    ATP["description"] = [str(DESC.get((s, k), ""))[:90] for s, k in zip(ATP.strain, ATP.key)]
    ATP.to_csv(DIR["table"]/"motif_atp_3strain.csv", index=False, encoding="utf-8-sig")
    print("=== 균주별 ATP 모티프 단백질 수 / 프로테옴 대비 ===")
    for st in STRAINS:
        n, tot = int((ATP.strain == st).sum()), len(PROT.get(st, {}))
        if tot: print(f"  {st:7s} {n:5d} / {tot:5d} = {100*n/tot:5.1f}%")
    print("\n→ 이 비율이면 ATP 모티프 단독 선별력은 없다. 3단계 교집합에서만 쓴다.")
    print(f"저장: {DIR['table']/'motif_atp_3strain.csv'}")
```

---

## CELL 36 — Folddisco 노선 (3)(4) 교집합 · 우선순위 · 최종표

```python
# =============================================================================
# CELL 36 | 3단계+4단계. 두 모티프를 겹쳐 우선순위를 매긴다
#   Tier 1  금속 + ATP 둘 다             — CooC1 의 두 성질을 모두 가진 것
#   Tier 2  금속만 + 기능 어노테이션 있음 — ATP 는 안 잡혔지만 정체가 분명한 것
#   Tier 3  금속만 + 어노테이션 없음      — 삭제하지 않는다. 우선도만 낮춘다
#   같은 Tier 안에서는 (1) BL21 특이 (2) 검출 차이 (3) RMSD 낮은 순.
#   ※ PPI 점수는 순위에 넣지 않는다. RF2-PPI 는 대조군 실패, Boltz ipTM 은 눈금
#     미검증이다. 참고 열로만 붙인다 — 컷오프로 쓰지 말 것.
# =============================================================================
need((have("MET"), "CELL 34 을 먼저 돌릴 것"))
_atp = ATP if ("ATP" in dir() and len(ATP)) else pd.DataFrame(columns=["strain","key","idf","rmsd","residues","n_match"])

T = MET.rename(columns={"idf": "metal_idf", "rmsd": "metal_rmsd",
                        "residues": "metal_residues", "n_match": "metal_n_match"}).copy()
T = T.merge(_atp[["strain","key","idf","rmsd","residues","n_match"]].rename(
                columns={"idf": "atp_idf", "rmsd": "atp_rmsd",
                         "residues": "atp_residues", "n_match": "atp_n_match"}),
            on=["strain", "key"], how="left")
T["has_atp"] = T.atp_idf.notna().astype(int)
T["tier"] = [1 if r.has_atp else (2 if r.annotated else 3) for _, r in T.iterrows()]
T["bl21_specific"] = T.strain_verdict.str.startswith("BL21 특이").astype(int)
T["bl21_detect_diff"] = T.strain_verdict.str.startswith("검출 차이").astype(int)

# 참고 열: Boltz ipTM (컷오프 아님)
_ip = {}
for _f in ["metal_motif_3strain_boltz.csv", "trackB_boltz2_ranked.csv"]:
    _p = DIR["table"]/_f
    if not _p.exists(): continue
    _b = pd.read_csv(_p)
    _k = "protein" if "protein" in _b.columns else "prey"
    if "iptm" in _b.columns:
        _ip.update({str(x): y for x, y in zip(_b[_k], _b["iptm"]) if pd.notna(y)})
if _ip: T["boltz_iptm"] = T.key.map(_ip)

T = T.sort_values(["tier", "bl21_specific", "bl21_detect_diff", "metal_rmsd"],
                  ascending=[True, False, False, True])
OUT = ["tier","strain","key","gene","description","annotated",
       "metal_idf","metal_rmsd","metal_residues","metal_n_match",
       "has_atp","atp_idf","atp_rmsd","atp_residues",
       "nres","plddt","motif_pattern","gene_pattern","strain_verdict","bl21_specific","boltz_iptm"]
T = T[[c for c in OUT if c in T.columns]].rename(
        columns={"key": "protein_id", "nres": "length", "plddt": "struct_plddt"})
T.to_csv(DIR["table"]/"FOLDDISCO_candidates.csv", index=False, encoding="utf-8-sig")

print("=== Tier x 균주 ===")
print(pd.crosstab(T.tier, T.strain).to_string())
print("\n=== Tier 1 (금속 + ATP 둘 다) — BL21 ===")
_t1 = T[(T.tier == 1) & (T.strain == "BL21")]
print(_t1[["protein_id","gene","metal_rmsd","metal_residues","atp_residues",
           "strain_verdict","description"]].to_string(index=False) if len(_t1) else "  없음")
print(f"\n=== Tier 2 (금속만 + 어노테이션) — BL21 상위 15 ===")
_t2 = T[(T.tier == 2) & (T.strain == "BL21")].head(15)
print(_t2[["protein_id","gene","metal_rmsd","metal_residues","strain_verdict",
           "description"]].to_string(index=False) if len(_t2) else "  없음")
print(f"\n=== Tier 3 (어노테이션 없음): {int((T.tier == 3).sum())}개 — 보관 ===")
print(f"\n최종표: {DIR['table']/'FOLDDISCO_candidates.csv'}")
print("  metal_residues = 가장 잘 맞은 매치의 잔기 / metal_rmsd = 그 매치의 RMSD")
print("  motif_pattern  = folddisco 가 잡은 균주 / gene_pattern = 유전자가 있는 균주")
```

---

## CELL 36b — 전체 후보표 (세 균주 · 모든 증거 열 · 숨김 없음)

```python
# =============================================================================
# CELL 36b | 금속 모티프 단백질 전부를 한 표에 올린다 — 거르지 않는다
#   세 균주 80개 전부, 가진 증거 열을 전부 붙인다. 검증에 실패한 점수도 지우지 않고
#   '검증 상태'를 명시해 함께 싣는다. 무엇을 근거로 삼고 무엇을 참고로만 볼지는
#   표를 다 본 사람이 정한다.
#
#   열의 신뢰 등급 (표 아래에도 다시 출력한다)
#     [확정]  metal_rmsd / metal_residues / motif_pattern / gene_pattern
#             구조 모티프 검색과 서열 대응. 대조군이 필요 없는 직접 측정이다.
#     [주의]  atp_match   질의 3잔기 중 몇 개가 맞았나. 2/3 부분 매치가 많다.
#                         그리고 ATP 모티프는 프로테옴의 23~26% 가 갖고 있다.
#     [참고]  ni_captured Boltz 가 Ni 를 후보 쪽 사슬로 가져갔나. 대상 27개가 전부
#                         folddisco 히트라 집합 밖으로 일반화 불가.
#     [미검증] boltz_iptm  이 계의 눈금이 정해지지 않았다 (part7_boltz_rep 대기)
#     [검증실패] rf2ppi   양성대조군 CooC1-ChCODH2 가 0.266, 531등 아래 (CELL 24b)
# =============================================================================
need(((DIR["table"]/"FOLDDISCO_candidates.csv").exists(), "CELL 36 을 먼저 돌릴 것"))
C = pd.read_csv(DIR["table"]/"FOLDDISCO_candidates.csv")

# --- ATP 부분 매치 ---
def atp_frac(s):
    if not isinstance(s, str) or not s.strip(): return ""
    p = s.split(",")
    return f"{sum(1 for x in p if x.strip() != '_')}/{len(p)}"
C["atp_match"] = C.get("atp_residues", pd.Series([""]*len(C))).apply(atp_frac)

# --- Boltz Ni 배치 (있는 열은 전부 가져온다) ---
_nip = DIR["table"]/"boltz_ni_placement.csv"
if _nip.exists():
    _n = pd.read_csv(_nip)
    _n["prot"] = _n.pair.astype(str).str.split("__").str[-1].str.split("-").str[-1]
    _n = _n.drop_duplicates("prot").set_index("prot")
    for col, new in [("위치", "ni_where"), ("A_배위", "ni_coord_chainA"),
                     ("B_배위", "ni_coord_chainB"), ("B_최단", "ni_dist_chainB")]:
        if col in _n.columns: C[new] = C.protein_id.map(_n[col])
    C["ni_captured"] = (C.get("ni_where", pd.Series([""]*len(C)))
                         .astype(str).str.startswith("후보 쪽").astype(int))
else:
    C["ni_captured"] = 0
    print("boltz_ni_placement.csv 없음 (CELL 28a-7)")

# --- RF2-PPI 점수: 검증 실패했지만 지우지 않는다 ---
_rf = DIR["table"]/"trackA_RF2PPI_ranked.csv"
if _rf.exists():
    _r = pd.read_csv(_rf, index_col=0)
    C["rf2ppi"]       = C.protein_id.map(_r["mean"])
    C["rf2ppi_sd"]    = C.protein_id.map(_r["sd"])
    C["paired_depth"] = C.protein_id.map(_r["paired_depth"])
    _ctrl = 0.266
    C["rf2_vs_control"] = np.where(C.rf2ppi.notna(),
                                   np.where(C.rf2ppi > _ctrl, "대조군보다 위", "아래"), "")
    print(f"RF2-PPI 점수 붙은 단백질: {int(C.rf2ppi.notna().sum())} / {len(C)}")
    print(f"  (양성대조군 {_ctrl} 기준 — 이 수치는 순위 근거가 아니라 기록이다)")

# --- 추천 사유: 다는 것뿐, 거르지 않는다 ---
PICK = {
 "QJZ12719.1": "COG0523 금속 샤페론 GTPase (UreG·HypB·CobW 과). 금속+뉴클레오타이드 두 모티프",
 "QJZ10977.1": "Fe-S 클러스터 삽입 단백질. C-cluster 가 [NiFe4S4] 라 직접 관련. "
               "Hyp KO 음성은 hydrogenase 전용 경로를 배제한 것이지 Fe-S 운반 경로는 아니다",
 "QJZ13803.1": "Hsp33 — Cys4 아연 중심을 산화환원에 따라 잡았다 놓는 것이 메커니즘",
 "QJZ11021.1": "HAD 포스파타제의 2가 금속 자리. 금속+뉴클레오타이드 두 모티프",
 "QJZ11953.1": "SEC-C/zinc-ribbon 계열, 기능 미상. 강한 금속 기하 + 기능 정보 없음",
 "QJZ12547.1": "UPF0149 계열, 기능 불명확. metal_rmsd 최저",
 "QJZ11161.1": "Tgt. 기능 확립된 효소라 사전확률 낮음. Ni 포획에는 포함",
 "QJZ14223.1": "BL21 특이 2개 중 구조 품질이 좋은 쪽. 전사조절자라 lysate rescue 와 "
               "기전이 안 맞고 ATP 도 2/3 부분 매치",
 "QJZ13177.1": "BL21 특이지만 pLDDT 56.7 — 구조 예측 신뢰 구간 밖",
}
C["note"] = [PICK.get(p, "") for p in C.protein_id]

C.to_csv(DIR["table"]/"ALL_motif_candidates_full.csv", index=False, encoding="utf-8-sig")

# ---------- 출력: 전부 ----------
pd.set_option("display.max_rows", None)
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 55)

SHOW = [c for c in ["strain", "protein_id", "gene", "tier", "metal_rmsd", "metal_residues",
                    "atp_match", "ni_captured", "ni_where", "rf2ppi", "boltz_iptm",
                    "struct_plddt", "length", "motif_pattern", "gene_pattern",
                    "strain_verdict", "description"] if c in C.columns]
for st in ["BL21", "MG1655", "Y19"]:
    sub = C[C.strain == st].sort_values(["tier", "metal_rmsd"])
    print(f"\n{'='*120}\n=== {st} — {len(sub)}개 (전부) ===\n{'='*120}")
    print(sub[SHOW].to_string(index=False))

print(f"\n{'='*120}")
print("=== 메모가 붙은 항목 ===")
_p = C[C.note != ""]
for _, r in _p.iterrows():
    print(f"\n{r.protein_id} ({r.strain})  rmsd {r.metal_rmsd:.3f}  Tier {r.tier}"
          f"  ATP {r.atp_match or '-'}  Ni포획 {'O' if r.ni_captured else '-'}")
    print(f"  {r.description}")
    print(f"  {r.note}")

print(f"\n{'='*120}")
print("열 신뢰 등급")
print("  [확정]    metal_rmsd, metal_residues, motif_pattern, gene_pattern")
print("  [주의]    atp_match — 2/3 은 부분 매치. ATP 모티프는 프로테옴의 23~26% 가 보유")
print("  [참고]    ni_captured — 대상 27개가 전부 folddisco 히트. 집합 밖 일반화 불가")
print("  [미검증]  boltz_iptm — 이 계의 눈금 미정 (part7_boltz_rep 대기)")
print("  [검증실패] rf2ppi — 양성대조군 0.266, 스크리닝 531등 아래 (CELL 24b)")
print(f"\n전체표(모든 열): {DIR['table']/'ALL_motif_candidates_full.csv'}")
print(f"3균주 원본     : {DIR['table']/'motif_metal_3strain.csv'}")
print(f"Tier 표        : {DIR['table']/'FOLDDISCO_candidates.csv'}")
```

---

## CELL 36c — 직교체 재현성 + BL21 단독 검출의 원인 규명

```python
# =============================================================================
# CELL 36c | 표에서 바로 읽히는 세 가지를 수치로 확정한다
#   (A) 같은 단백질이 균주마다 같은 값을 받는가 = 각 열의 잡음 바닥
#       세 균주에 같은 단백질(직교체)이 여럿 있다. 입력이 사실상 같으므로 균주 간
#       차이는 측정 잡음이다. part7_boltz_rep 을 기다릴 것 없이 지금 계산된다.
#   (B) BL21 에서만 '검출'된 단백질 — 유전자는 세 균주에 다 있는데 folddisco 가
#       BL21 에서만 잡았다. 서열 자체가 다른가(진짜), 아니면 같은 서열인데 구조
#       예측만 다른가(인공물)? 모티프 잔기가 직교체에서도 Cys 인지 보면 갈린다.
#   (C) 앞 셀의 버그 둘: ni_captured 문자열 불일치, BL21/Y19 유전자명 미추출
# =============================================================================
need((have("GENE_HIT", "MOTIF_HIT", "PROT", "C"), "CELL 33 -> 34 -> 36 -> 36b 를 먼저 돌릴 것"))
import difflib

# ---------- (C1) ni_captured 수정 ----------
# CELL 28a-7 이 쓰는 문자열은 'B쪽만 (후보)' 인데 앞 셀은 '후보 쪽' 으로 찾았다.
if "ni_where" in C.columns:
    C["ni_captured"] = C.ni_where.astype(str).str.contains("후보").astype(int)
    print(f"[수정] Ni 를 후보 쪽 사슬로 가져간 단백질: {int(C.ni_captured.sum())}개")
    print(C[C.ni_captured == 1][["protein_id", "description"]].to_string(index=False))
    print("\n※ ni_where 가 NaN 인 행은 Boltz Ni 배치를 계산한 적이 없다는 뜻이다")
    print(f"   (계산된 것 {int(C.ni_where.notna().sum())} / {len(C)}. MG1655·Y19 는 out_cmp 를")
    print("    CELL 28a-7 이 안 읽어서 통째로 비어 있다 — 균주 비교에는 못 쓴다)")

# ---------- (C2) BL21/Y19 유전자명 ----------
# 설명문의 관용 표기(Tgt, YeiR, HslO ...)를 잡고, 없으면 MG1655 직교체 이름을 빌린다.
_mgname = {}
for (st, k), gname in GENE.items():
    if st == "MG1655" and gname: _mgname[k] = gname
def pick_gene(row):
    if isinstance(row.get("gene"), str) and row["gene"]: return row["gene"]
    d = str(row.get("description", "")).split("[")[0]
    cands = re.findall(r"\b([A-Z][a-z]{2}[A-Z0-9]{0,3})\b", d)
    if cands: return cands[-1][0].lower() + cands[-1][1:]
    m = MOTIF_HIT.get((row["strain"], row["protein_id"]), {}).get("MG1655")
    return _mgname.get(m[1], "") if m else ""
C["gene"] = [pick_gene(r) for _, r in C.iterrows()]

# ---------- (A) 직교체 재현성 ----------
def counterparts(strain, key):
    d = MOTIF_HIT.get((strain, key), {})
    return {s: v[1] for s, v in d.items() if s != strain}

idx = C.set_index(["strain", "protein_id"])
rows = []
for _, r in C[C.strain == "BL21"].iterrows():
    cp = counterparts("BL21", r.protein_id)
    if not cp: continue
    e = {"protein": r.protein_id, "gene": r.gene, "desc": str(r.description)[:45],
         "BL21_rmsd": r.metal_rmsd, "BL21_iptm": r.get("boltz_iptm")}
    for s, k in cp.items():
        if (s, k) in idx.index:
            o = idx.loc[(s, k)]
            e[f"{s}_rmsd"] = float(o.metal_rmsd)
            e[f"{s}_iptm"] = float(o.boltz_iptm) if pd.notna(o.get("boltz_iptm")) else np.nan
    rows.append(e)
O = pd.DataFrame(rows)
for col in ["rmsd", "iptm"]:
    cs = [c for c in O.columns if c.endswith(f"_{col}")]
    O[f"d_{col}"] = O[cs].max(axis=1) - O[cs].min(axis=1)

print("\n" + "=" * 100)
print("=== (A) 같은 단백질의 균주 간 차이 = 잡음 바닥 ===")
print(O[[c for c in O.columns if c != "desc"]].round(3).to_string(index=False))
for col, label in [("rmsd", "metal_rmsd"), ("iptm", "boltz_iptm")]:
    d = O[f"d_{col}"].dropna()
    if not len(d): continue
    print(f"\n{label}: 직교체 {len(d)}쌍  |Δ| 중앙값 {d.median():.3f}  "
          f"평균 {d.mean():.3f}  최대 {d.max():.3f}")

_rng = C[C.strain == "BL21"].boltz_iptm
if _rng.notna().any() and "d_iptm" in O:
    _n = O.d_iptm.median()
    print(f"\n비교: BL21 ipTM 전체 범위는 {_rng.min():.3f} ~ {_rng.max():.3f} "
          f"(폭 {_rng.max()-_rng.min():.3f})")
    print(f"      같은 단백질끼리도 {_n:.3f} 씩 벌어진다.")
    if _n > (_rng.max() - _rng.min()) / 3:
        print("      → 잡음이 전체 폭의 1/3 을 넘는다. ipTM 순위는 읽을 수 없다.")

# ---------- (B) BL21 에서만 검출된 것들: 서열이 다른가 ----------
def map_pos(s1, s2, i):
    """s1 의 0-based 위치 i 가 s2 의 어느 위치에 대응하는지"""
    for tag, i1, i2, j1, j2 in difflib.SequenceMatcher(None, s1, s2).get_opcodes():
        if i1 <= i < i2:
            return (j1 + (i - i1)) if tag == "equal" else (j1 if tag == "replace" else None)
    return None

only = C[(C.strain == "BL21") & (C.motif_pattern == "BL21") &
         (C.gene_pattern.astype(str).str.contains("MG1655"))]
print("\n" + "=" * 100)
print(f"=== (B) 유전자는 다른 균주에도 있는데 BL21 에서만 검출된 것: {len(only)}개 ===")
print("모티프 잔기가 직교체에서도 Cys 인지 본다. Cys 가 보존돼 있으면 서열이 아니라")
print("구조 예측이 달랐다는 뜻이다 (인공물). Cys 가 없으면 진짜 차이다.\n")
for _, r in only.iterrows():
    bs = PROT["BL21"].get(r.protein_id, "")
    pos = [int(x[1:]) for x in str(r.metal_residues).split(",") if x[1:].isdigit()]
    print(f"--- {r.protein_id} {r.gene}  {str(r.description)[:60]}")
    print(f"    BL21 {r.metal_residues} = " +
          ",".join(f"{bs[p-1]}{p}" if 0 < p <= len(bs) else "?" for p in pos))
    for st in ["MG1655", "Y19"]:
        hit = GENE_HIT.get(("BL21", r.protein_id), {}).get(st)
        if not hit:
            print(f"    {st:7s} 직교체 없음"); continue
        fid, tid = hit
        os_ = PROT[st].get(tid, "")
        if not os_:
            print(f"    {st:7s} {tid} 서열 없음"); continue
        aa = []
        for p in pos:
            j = map_pos(bs, os_, p - 1)
            aa.append(f"{os_[j]}{j+1}" if j is not None and j < len(os_) else "-")
        same = all(x and x[0] == "C" for x in aa)
        print(f"    {st:7s} {tid} (fident {fid:.3f}) -> " + ",".join(aa) +
              ("   Cys 보존 → 구조 예측 차이" if same else "   Cys 아님 → 서열 차이"))
    print()

C.to_csv(DIR["table"]/"ALL_motif_candidates_full.csv", index=False, encoding="utf-8-sig")
O.to_csv(DIR["table"]/"ortholog_reproducibility.csv", index=False, encoding="utf-8-sig")
print("=" * 100)
print(f"갱신: {DIR['table']/'ALL_motif_candidates_full.csv'}")
print(f"신규: {DIR['table']/'ortholog_reproducibility.csv'}")
```

---

## 부록 A — 검증 지표 · 단위 · 컷오프 · 근거 (참조용, 실행 셀 아님)

각 숫자가 무엇이고, 어디서 온 기준으로 자르며, 그 기준을 왜 쓰는지. 그리고
**이 계에서 검증됐는지**. 마지막 열이 핵심이다 — 벤치마크에서 통하는 컷오프가
여기서도 통한다는 보장은 없고, 실제로 두 개는 이미 깨졌다.

### A-1. Folddisco (구조 모티프 검색) — Track C

| 지표 | 단위 | 쓴 값 | 왜 | 이 계에서 |
|---|---|---|---|---|
| `--rmsd` | Å | **1.0** | 질의 잔기 2개면 관계쌍이 1개뿐이라 기하 제약이 약하다. 1.0 Å 이 넘어가면 아무 Cys 쌍이나 걸린다 | **확정** — 직교체끼리 \|Δrmsd\| 가 대개 0.02 미만. 재현된다 |
| `metal_rmsd` | Å | 낮을수록 좋음 | CooC1 의 A112/A114 배치와 얼마나 같은 기하인가 | 동상 |
| `idf` | 없음 | 정렬용 | 모티프 희귀도. **후보 선별에는 쓰지 않는다** — AlphaFold DB 에서 희귀한 기하는 드문 기능이 아니라 엉망으로 예측된 단편인 경우가 많다 (실측: idf 1위가 44잔기, pLDDT 56.7) | **쓰지 않음** |
| `--top` | 개 | 2000 | 균주당 상한. 실제 히트는 metal 21~30개라 걸리지 않는다 | — |
| `atp_match` | n/3 | **3/3 만 완전** | Walker A 질의 3잔기 중 몇 개가 맞았나. 2/3 은 부분 매치다 | **주의** — Tier 1 대부분이 2/3 |
| ATP 모티프 자체 | — | 단독 필터 금지 | 프로테옴의 **23~26%** 가 보유. 선별력이 없다. CooC1 의 최근접 이웃이 MinD 로 나온 것도 이것 때문 | **선별력 없음 확정** |
| `struct_plddt` | 0–100 | **< 70 은 기하를 믿지 말 것** | AlphaFold 신뢰도. 낮은 영역의 곁사슬 배치는 예측이 아니라 잡음이다 | RseC 3균주 rmsd 가 0.31/0.56/0.48 로 흩어졌는데 셋 다 pLDDT 55~73 |

### A-2. mmseqs2 (균주 간 대응)

| 지표 | 단위 | 쓴 값 | 왜 |
|---|---|---|---|
| `fident` | 0–1 | **≥ 0.50** | 직교체 판정 하한. 장내세균 3종 사이 직교체는 보통 0.8 이상이라 여유 있는 선 |
| `qcov` | 0–1 | **≥ 0.70** | 도메인 하나만 겹친 것을 직교체로 세지 않기 위해 |
| `evalue` | 기대값 | ≤ 1e-5 | 검색 단계 관용 컷오프 |

두 대응표를 **분리해서** 만든다. 섞으면 '검출 차이'와 '유전자 부재'가 구분되지 않는다.
- 모티프 집합 ↔ 모티프 집합 → `motif_pattern` (folddisco 가 잡은 균주)
- 모티프 집합 → 전체 프로테옴 → `gene_pattern` (유전자가 있는 균주)

### A-3. RF2-PPI (공진화 PPI 분류) — Track A

| 지표 | 단위 | 기준 | 왜 | 이 계에서 |
|---|---|---|---|---|
| `mean` (prob) | 확률 0–1 | ≥0.74 strict / 0.30–0.74 회색 / <0.30 배제 | 논문 벤치마크에서 0.74 가 precision 95% 선 | **검증 실패** |
| `sd` | 확률 | > 0.1 이면 재실행 | 비결정적 모델. README 기준 ~5% 쌍이 SD>0.1 | 유효 |
| `N_REPLICATE` | 회 | 깊이≥500 이면 3, 아니면 5 | 위와 같은 이유 | 유효 |
| `paired_depth` | **유전체 개수** | MIN_PAIRED 50 / <200 위음성 의심 | 공진화를 보려면 두 단백질이 같은 유전체에 있어야 한다. 그 짝의 수 | 유효 |
| gate 깊이 | 유전체 개수 | ≥ 500 | 이보다 얕으면 공진화 신호가 안 선다 | 통과 (2,929) |

**왜 검증 실패인가.** 양성대조군 CooC1–ChCODH2 가 **0.266**. 0.30 미만이니 모델이
"안 붙는다"고 답한 것이다. 그런데 이 쌍은 문헌상 실제로 붙고, cooC/cooS 는 1,451개
유전체에서 같은 오페론에 산다 — 공진화가 최대로 잡혀야 할 교과서적 케이스다.
게다가 깊이 1,451 은 스크리닝 쌍의 75분위(1,392)보다 위라 재료 부족도 아니다.
이보다 높은 점수를 받은 스크리닝 후보가 **531 / 3,625 (14.6%)**. 정답을 목록에
넣으려면 531개를 다 가져와야 한다 → 후보 지명에 쓸 수 없다.

### A-4. Boltz-2 (co-folding) — Track B

| 지표 | 단위 | 통상 기준 | 왜 | 이 계에서 |
|---|---|---|---|---|
| `iptm` | 0–1 | ≥0.8 강함 / 0.6–0.8 유의 / <0.4 배제 | 계면만 보는 predicted TM-score. PPI 판정의 표준 | **눈금 미검증 + 잡음 과다** |
| `ptm` | 0–1 | — | 복합체 전체. 큰 사슬이 지배해 계면 판정에 못 씀 | 쓰지 않음 |
| `ligand_iptm` | 0–1 | — | Ni 이 계면에 놓였나 | **폐기** — 길이와 ρ=−0.807, 348개 중 계면 배치 0건 |
| `complex_plddt` | 0–1 | 높을수록 | 구조 신뢰도 | 참고 |
| Ni 배위 거리 | **Å** | ≤ **3.5** | Ni(II)–S/N 배위 결합 길이(1.9–2.2 Å)에 여유를 둔 상한 | 유효 (실측 1.91–2.18 Å) |
| `--max_msa_seqs` | 줄 | 2048 | VRAM 23.8→14.6 GB. 8192 로는 15%가 OOM 으로 조용히 버려졌다 | 유효 |
| 입력 길이 | 잔기 | ≤ 1830 | bait+prey 합. 24 GB 에서 안전한 선 | 유효 |

**왜 눈금 미검증인가.** Track B/C 에는 양성대조군을 물린 적이 없다. ipTM 0.68 이
이 계에서 무슨 뜻인지 정하는 기준이 없다. `part7_boltz_rep` 이 CooC1–ChCODH2 를
5표본으로 접어 그 기준을 만드는 중이다.

**왜 잡음 과다인가.** 직교체(사실상 같은 입력)끼리도 ipTM 이 크게 벌어진다.
- yjjP/ThrE (BL21 256aa 0.681 vs MG1655 256aa 0.397) → **Δ0.284**
- fumA/fumB (같은 BL21 안의 유사 파라로그, 548aa 동일) 0.365 vs 0.213 → **Δ0.152**
- 반대로 yraI (0.627 vs 0.612) 는 Δ0.015

BL21 ipTM 전체 범위가 0.17–0.68 (폭 0.51) 인데 같은 단백질끼리 0.15~0.28 이 벌어진다.
순위 차이의 상당 부분이 잡음이다. CELL 36c 가 이걸 전체 직교체에서 계산한다.

### A-5. 지금 후보 선별에 실제로 쓰는 기준

검증된 것만 남기면 이렇게 된다.

| 순서 | 기준 | 근거 등급 |
|---|---|---|
| 1 | 금속 모티프 존재 (folddisco, rmsd 1.0 Å 이내) | 확정 — 직교체 재현성 확인됨 |
| 2 | `struct_plddt` ≥ 70 | 확정 — 이하는 기하 자체가 예측이 아님 |
| 3 | 금속 + 뉴클레오타이드 두 모티프 (Tier 1), 단 ATP 3/3 만 온전 | 주의 |
| 4 | 기능 어노테이션 유무 (Tier 2/3) | 보조 |
| 5 | 생화학적 개연성 — 금속 샤페론 과, Fe-S 운반, 금속 잡았다 놓는 기전 | **사람 판단** |

**PPI 점수(rf2ppi, boltz_iptm)는 1~5 어디에도 들어가지 않는다.** 표에는 싣되
순위 근거로 쓰지 않는다. 하나는 대조군에서 떨어졌고, 하나는 눈금이 없다.

**균주 특이성은 필터가 아니다.** BL21 금속 모티프 29개 중 27개가 MG1655 에도
유전자가 있다 — 유전자 유무로는 표현형이 설명되지 않는다. 그리고 *in vitro*
재구성 assay 는 "MG1655 에도 있나"를 묻지 않고 "이 단백질이 Ni 를 넣나"를 묻는다.

---

## CELL 39 — 기준 통일 + 최종 후보표 (미팅용)

```python
# =============================================================================
# CELL 39 | 흩어진 결과를 하나로 합치고 임계값을 통일한다
#   지금 판정이 셀마다 다른 기준을 쓰고 있다. 여기서 한 벌로 고정하고, 경계선에
#   있던 후보(GspE fident 0.600, GspH 0.285)를 그 기준으로 다시 판정한다.
#
#   통일 기준
#     직교체        fident >= 0.50 & qcov >= 0.70 & e <= 1e-5
#     균주 구간     identical >=0.999 / high >=0.90 / low <0.90 / 히트없음 = specific
#     게놈 동일유전자 fident >= 0.80         (이보다 낮으면 파라로그이지 깨진 사본이 아니다)
#     게놈 pseudogene fident >=0.80 & (qcov <0.70 또는 조각 >=3)
#     구조 신뢰      plddt >= 70            (이하는 기하 자체가 예측이 아니다)
#     Ni 배위        <= 2.6 A 만 배위로 센다 (Ni-S 2.2 / Ni-N 2.1)
#
#   합치는 표
#     FOLDDISCO_candidates.csv        금속/ATP 모티프 + 균주 판정
#     ni_site_grade.csv               Ni 배위 도너 등급 (A/B/C/D)
#     cys4_bl21_strain_specificity.csv Cys4 x 프로테옴 균주 특이
#     cys4_genome_check.csv           게놈 확인
#     ALL_motif_candidates_full.csv   PPI 점수 (참고 열)
# =============================================================================
TH = {"orth_fident": 0.50, "orth_qcov": 0.70, "orth_evalue": 1e-5,
      "identical": 0.999, "high": 0.90,
      "genome_same_gene": 0.80, "pseudo_qcov": 0.70, "pseudo_frag": 3,
      "plddt": 70.0, "ni_coord": 2.6}
print("통일 기준:", json.dumps(TH, ensure_ascii=False))

T = DIR["table"]
def rd(name):
    p = T/name
    if p.exists():
        d = pd.read_csv(p); print(f"  {name}: {len(d)}행"); return d
    print(f"  {name}: 없음"); return pd.DataFrame()

print("\n읽는 표:")
FC  = rd("FOLDDISCO_candidates.csv")
NG  = rd("ni_site_grade.csv")
CS  = rd("cys4_bl21_strain_specificity.csv")
GC  = rd("cys4_genome_check.csv")
AL  = rd("ALL_motif_candidates_full.csv")

# ---- 후보 풀: 금속 모티프 BL21 + Cys4 BL21 합집합 ----
pool = set()
if len(FC): pool |= set(FC[FC.strain == "BL21"].protein_id.astype(str))
if len(CS): pool |= set(CS.protein.astype(str))
if len(NG):
    pool |= {p for p in NG.prey.astype(str) if p.startswith("QJZ")}
print(f"\n후보 풀 (BL21): {len(pool)}개")

hdr = {}
for l in open(ASSET["faa_bl21"], errors="ignore"):
    if l.startswith(">"):
        t = l[1:].rstrip(); hdr[t.split()[0]] = t[len(t.split()[0]):].strip()

def get(df, key, col, k):
    if not len(df) or key not in df.columns or col not in df.columns: return None
    r = df[df[key].astype(str) == k]
    return r.iloc[0][col] if len(r) else None

# ---- 게놈 판정을 통일 기준으로 다시 ----
def genome_verdict(k):
    if not len(GC): return ""
    r = GC[GC.protein.astype(str) == k]
    if not len(r): return ""
    r = r.iloc[0]
    if str(r.get("genome")) == "없음": return "게놈에도 없음 [확정]"
    f, c = r.get("fident"), r.get("qcov")
    n = r.get("n_frag", 0)
    if pd.isna(f) or float(f) < TH["genome_same_gene"]:
        return "먼 상동체 — 다른 유전자 [균주 특이 유지]"
    if float(c) < TH["pseudo_qcov"] or (pd.notna(n) and int(n) >= TH["pseudo_frag"]):
        return "★ pseudogene 의심"
    return "게놈에 온전히 있음 [어노테이션 누락]"

rows = []
for k in sorted(pool):
    ni_g = get(NG, "prey", "grade", k)
    plddt = get(FC, "protein_id", "struct_plddt", k)
    rows.append({
        "protein": k,
        "desc": hdr.get(k, "")[:70],
        "length": get(FC, "protein_id", "length", k),
        "metal_rmsd": get(FC, "protein_id", "metal_rmsd", k),
        "metal_residues": get(FC, "protein_id", "metal_residues", k),
        "struct_plddt": plddt,
        "plddt_ok": (pd.notna(plddt) and float(plddt) >= TH["plddt"]) if plddt is not None else None,
        "ni_grade": ni_g,
        "ni_donors": get(NG, "prey", "donors", k),
        "cys4_query": get(CS, "protein", "queries", k),
        "vs_MG1655": get(CS, "protein", "vs_MG1655", k) or get(FC, "protein_id", "strain_verdict", k),
        "fident_MG1655": get(CS, "protein", "fident_MG1655", k),
        "genome": genome_verdict(k),
        "boltz_iptm": get(AL, "protein_id", "boltz_iptm", k),
        "rf2ppi": get(AL, "protein_id", "rf2ppi", k),
    })
F = pd.DataFrame(rows)

# ---- 등급제: 어떤 행도 제외하지 않는다. 우선순위만 낮춘다 ----
#   필터링은 위험하다. 걸러낸 것이 답이면 되돌릴 방법이 없다. 그래서 감점만 하고
#   전부 표에 남긴다. 감점 사유는 reason 열에 그대로 적는다.
#
#   가점 (검증된 축만)
#     +3  Ni 배위 grade A (Cys4)   — 독립 실행 간 Δ 0.01~0.08 A. 가장 재현성 높다
#     +2  Ni 배위 grade B (Cys3His)
#     +2  금속 모티프 rmsd <= 1.0   — 직교체 간 Δ 중앙값 0.015 A
#     +1  Cys4 재질의 히트
#     +2  균주 특이 (프로테옴)      — MG1655 에 상동체 없음/멂
#     +1  게놈에서도 확인           — 어노테이션 문제가 아님
#   감점 (제외가 아니라 강등)
#     -2  pLDDT < 70               — 기하 자체가 예측이 아니다
#     -1  MG1655 에 온전한 유전자   — 균주 축의 근거가 없다. 생화학 근거는 남는다
#     -1  Ni 배위 grade D / 없음
#   ※ PPI 점수(rf2ppi, boltz_iptm)는 가점에도 감점에도 쓰지 않는다. 둘 다 자기
#     양성대조군에서 떨어졌다 (0.266 / ipTM 0.300).
def score_row(r):
    pts, why = 0, []
    if r.ni_grade == "A":            pts += 3; why.append("+3 Ni배위 Cys4")
    elif r.ni_grade == "B":          pts += 2; why.append("+2 Ni배위 Cys3His")
    elif r.ni_grade in ("C", "D"):   pts -= 1; why.append("-1 Ni배위 약함")
    if pd.notna(r.metal_rmsd) and float(r.metal_rmsd) <= 1.0:
        pts += 2; why.append(f"+2 금속모티프 {float(r.metal_rmsd):.3f}A")
    if isinstance(r.cys4_query, str) and r.cys4_query:
        pts += 1; why.append("+1 Cys4 재질의")
    if str(r.vs_MG1655) in ("strain-specific", "low-similarity") or \
       str(r.vs_MG1655).startswith(("BL21 특이", "검출 차이")):
        pts += 2; why.append(f"+2 균주특이({r.vs_MG1655})")
    g = str(r.genome)
    if g.startswith(("게놈에도 없음", "★")):
        pts += 1; why.append("+1 게놈확인")
    elif g.startswith("게놈에 온전히"):
        pts -= 1; why.append("-1 MG1655 에 유전자 온전")
    if r.plddt_ok is False:
        pts -= 2; why.append(f"-2 pLDDT {r.struct_plddt:.0f} (<70)")
    return pts, " / ".join(why)

_sc = [score_row(r) for _, r in F.iterrows()]
F["score"]  = [x[0] for x in _sc]
F["reason"] = [x[1] for x in _sc]
F["tier"] = pd.cut(F.score, bins=[-99, 0, 2, 4, 6, 99],
                   labels=["E", "D", "C", "B", "A"]).astype(str)
F = F.sort_values(["score", "metal_rmsd"], ascending=[False, True], na_position="last")
F.to_csv(T/"FINAL_candidates.csv", index=False, encoding="utf-8-sig")

pd.set_option("display.max_rows", None); pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 42)
COLS = ["tier", "score", "protein", "ni_grade", "metal_rmsd", "vs_MG1655",
        "genome", "struct_plddt", "desc"]
print("\n" + "=" * 130)
print(f"=== 전체 {len(F)}개 — 제외된 행 없음 ===")
print(F[COLS].to_string(index=False))
print("\n=== 등급 분포 ===")
print(F.tier.value_counts().reindex(["A","B","C","D","E"]).fillna(0).astype(int).to_string())
print("\n=== 상위 등급의 감점·가점 내역 ===")
for _, r in F[F.tier.isin(["A", "B"])].iterrows():
    print(f"\n[{r.tier}] {r.protein}  score {r.score}")
    print(f"    {r.desc}")
    print(f"    {r.reason}")
print(f"\n저장: {T/'FINAL_candidates.csv'}")
print("\n※ 감점은 제외가 아니다. E 등급도 표에 그대로 있고, 실험에서 뒤집힐 수 있다.")
print("\n※ evidence 에 PPI 점수는 들어가지 않는다. rf2ppi 는 대조군 실패,")
print("  boltz_iptm 은 직교체 간 잡음 0.102(최대 0.309)로 전체 폭 0.485 의 21~64%.")
print("  둘 다 열로만 싣는다.")
```

---

## CELL 40 — Y19 coo 오페론 확보 (양성대조군 세트)

```python
# =============================================================================
# CELL 40 | Y19 의 coo 오페론 전체를 뽑는다
#   Merrouch et al. 2018 (JBIC 23:613) Fig.2: CaY19 의 cooS 뒤에 CooC · CooT ·
#   CooJ · HypB 가 붙어 있다. 그리고 그 오페론에서 CODH 활성에 필수인 것은
#   놀랍게도 CooF 하나뿐이었다.
#   이 다섯이 우리가 한 번도 가져본 적 없는 '진짜 양성대조군 세트' 다.
#   CooC1 하나로는 대조군이 하나뿐이었고 세 트랙이 전부 그것을 못 맞혔다.
#
#   문헌상 크기 (상동체 찾을 때 기준으로 쓴다)
#     CooS  ~630 aa   CODH 본체
#     CooC  ~270 aa   MRP/MinD 계열 ATPase, CXC (Cys112/114) 이량체 계면 Ni 자리
#     CooT  ~65 aa    7.1 kDa, 이량체, Cys2 로 Ni 결합 (Kd 9 nM)
#     CooJ  ~115 aa   12.6 kDa, C-말단 His-rich (His 16개), Ni 4개/monomer
#     HypB  ~290 aa   GTPase
#     CooF  ~200 aa   Fe-S 단백질. Y19 에서 유일하게 필수였던 것
# =============================================================================
need((have("ASSET",), "CELL 01 을 먼저 돌릴 것"))
Y19 = {}
nm, buf, desc = None, [], {}
for l in open(ASSET["faa_y19"], errors="ignore"):
    if l.startswith(">"):
        if nm: Y19[nm] = "".join(buf)
        t = l[1:].rstrip(); nm = t.split()[0]; buf = []
        desc[nm] = t[len(nm):].strip()
    else: buf.append(l.strip())
if nm: Y19[nm] = "".join(buf)
print(f"Y19 프로테옴 {len(Y19)}개")

# AHZ* 은 coo 오페론 별도 제출본이다 (Ainala et al. 2015). 전부 본다.
AHZ = sorted(k for k in Y19 if k.startswith("AHZ"))
print(f"\n=== AHZ 계열 {len(AHZ)}개 (coo 오페론 제출본으로 추정) ===")
for k in AHZ:
    print(f"  {k}  {len(Y19[k]):5d} aa  {desc[k][:70]}")

# 이름으로도 훑는다
PAT = r"\bCoo[CSTJF]\b|carbon monoxide dehydrogenase|hydrogenase nickel|\bHypB\b|nickel-binding"
byname = [k for k, v in desc.items() if re.search(PAT, v, re.I)]
print(f"\n=== 이름으로 잡힌 것 {len(byname)}개 ===")
for k in sorted(byname):
    print(f"  {k}  {len(Y19[k]):5d} aa  {desc[k][:70]}")

# 크기로 후보를 좁힌다 (AHZ 안에서)
SIZE = {"CooS": (550, 720), "CooC": (220, 320), "CooT": (50, 90),
        "CooJ": (90, 150), "HypB": (240, 340), "CooF": (150, 250)}
print("\n=== 크기로 본 배정 후보 (AHZ 계열 안에서) ===")
ASSIGN = {}
for role, (lo, hi) in SIZE.items():
    cands = [k for k in AHZ if lo <= len(Y19[k]) <= hi]
    print(f"  {role:6s} ({lo}-{hi} aa): {cands if cands else '없음'}")
    for c in cands: ASSIGN.setdefault(role, []).append(c)

# 질의용 fasta
SET = sorted(set(AHZ) | set(byname))
qf = DIR["seq"]/"y19_coo_operon.faa"
qf.write_text("".join(f">{k}\n{Y19[k]}\n" for k in SET))
print(f"\n대조군 세트 {len(SET)}개 -> {qf}")
print("이 세트를 CELL 41 이 BL21·MG1655 에 던진다.")
print("\n※ AHZ 계열이 2개뿐이면 오페론 나머지가 AKE/AGE 로 들어가 있다는 뜻이다.")
print("  그 경우 CELL 41 의 서열·구조 검색이 이름과 무관하게 찾아낸다.")
```

---

## CELL 41 — coo 오페론 상동체를 BL21·MG1655 에서 찾는다 (Foldseek + mmseqs)

```python
# =============================================================================
# CELL 41 | 검증된 도구로 다시 찾는다
#   Folddisco 의 두-Cys 질의는 Y19 CooC 를 못 찾았고 (CELL 29j), Foldseek 은
#   tm 0.938 로 찾았다 (CELL 29k). 그러니 이번엔 Foldseek 을 쓴다.
#   서열(mmseqs)도 같이 돌려 둘이 일치하는지 본다 — CooT 는 65잔기라 구조
#   검색이 짧은 서열에서 불안정할 수 있고, 반대로 CooJ 의 His-rich 꼬리는
#   서열은 달라도 구조가 남을 수 있다.
# =============================================================================
need(((DIR["seq"]/"y19_coo_operon.faa").exists(), "CELL 40 을 먼저 돌릴 것"))
FS  = shutil.which("foldseek") or globals().get("FOLDSEEK_BIN", "foldseek")
QF  = DIR["seq"]/"y19_coo_operon.faa"
FAA = {"BL21": ASSET["faa_bl21"], "MG1655": ASSET["faa_mg1655"], "Y19": ASSET["faa_y19"]}
STRUCT = {}
for r in [TOOLS/"database"/"bacteriaDB", TOOLS/"database"]:
    if Path(r).exists():
        for d in sorted(Path(r).glob("structures_*")):
            if d.is_dir(): STRUCT[d.name] = d
_tmp = DIR["tmp"]/"coo"; _tmp.mkdir(parents=True, exist_ok=True)
outd = DIR["folddisco"].parent/"foldseek"; outd.mkdir(exist_ok=True)

HDR, SEQ = {}, {}
for st, fa in FAA.items():
    nm, buf = None, []
    for l in open(fa, errors="ignore"):
        if l.startswith(">"):
            if nm: SEQ[nm] = "".join(buf)
            t = l[1:].rstrip(); nm = t.split()[0]; buf = []
            HDR[nm] = (st, t[len(nm):].strip())
        else: buf.append(l.strip())
    if nm: SEQ[nm] = "".join(buf)

# ---------- 서열 검색 ----------
SFMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits"
print("=== 서열 검색 (mmseqs -s 7.5) ===")
SEQHIT = {}
for st in ["BL21", "MG1655"]:
    out = DIR["search"]/f"cooop_vs_{st}.m8"
    if not (out.exists() and out.stat().st_size):
        sh(f'mmseqs easy-search "{QF}" "{FAA[st]}" "{out}" "{_tmp}" '
           f'--format-output "{SFMT}" -e 1e-2 -s 7.5 --threads {min(THREADS,16)} -v 1',
           check=False)
    if out.exists() and out.stat().st_size:
        d = pd.read_csv(out, sep="\t", names=SFMT.split(","))
        SEQHIT[st] = d.sort_values("bits", ascending=False)
    print(f"  {st}: {len(SEQHIT.get(st, []))}행")

# ---------- 구조 검색 ----------
xw = pd.read_csv(DIR["table"]/"id_crosswalk_struct_to_genbank.csv") \
     if (DIR["table"]/"id_crosswalk_struct_to_genbank.csv").exists() else pd.DataFrame()
xwm = {Path(str(k)).stem: str(v).split(",")[0]
       for k, v in (xw.dropna(subset=["protein"]).set_index("tid")["protein"].to_dict().items()
                    if len(xw) else [])}
rev = {v: k for k, v in xwm.items()}
FFMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits,prob,alntmscore,lddt"

def struct_of(pid):
    stem = rev.get(pid)
    if not stem: return None
    for sdir in STRUCT.values():
        for ext in ("cif", "pdb"):
            f = next(iter(glob.glob(str(sdir/f"{stem}*.{ext}"))), None)
            if f: return f
    return None

print("\n=== 구조 검색 (foldseek) ===")
STRHIT = {}
for pid in [l[1:].strip() for l in open(QF) if l.startswith(">")]:
    sp = struct_of(pid)
    if not sp:
        print(f"  {pid}: 구조 없음 (서열 검색 결과만 씀)"); continue
    for sname, sdir in STRUCT.items():
        out = outd/f"coo_{pid}_vs_{sname}.m8"
        if not (out.exists() and out.stat().st_size):
            sh(f'"{FS}" easy-search "{sp}" "{sdir}" "{out}" "{_tmp}" '
               f'--format-output "{FFMT}" -e 10 --max-seqs 500 '
               f'--threads {min(THREADS,16)} -v 1', check=False)
        if out.exists() and out.stat().st_size:
            d = pd.read_csv(out, sep="\t", names=FFMT.split(","))
            d["stem"] = d.target.apply(lambda x: Path(str(x)).stem)
            d["protein"] = d.stem.map(xwm)
            STRHIT[(pid, sname)] = d

# ---------- 종합 ----------
print("\n" + "=" * 110)
print("=== 질의별 최상위 상동체 ===")
rows = []
for pid in [l[1:].strip() for l in open(QF) if l.startswith(">")]:
    qd = HDR.get(pid, ("Y19", ""))[1][:48]
    print(f"\n--- {pid}  {len(SEQ.get(pid,'')):4d} aa  {qd}")
    for st in ["BL21", "MG1655"]:
        s = SEQHIT.get(st)
        best_s = None
        if s is not None and len(s):
            m = s[s["query"] == pid]
            if len(m): best_s = m.iloc[0]
        if best_s is not None:
            print(f"    [서열 {st:7s}] {best_s.target:14s} fident {best_s.fident:.3f} "
                  f"qcov {best_s.qcov:.2f} e {best_s.evalue:.1e}  "
                  f"{HDR.get(best_s.target,('',''))[1][:44]}")
        else:
            print(f"    [서열 {st:7s}] 없음")
        # 구조
        bt = None
        for (q, sname), d in STRHIT.items():
            if q != pid: continue
            dd = d[d.protein.map(lambda k: HDR.get(str(k), ("?",))[0]) == st]
            dd = dd[dd.protein.astype(str) != pid]
            if len(dd):
                c = dd.sort_values("bits", ascending=False).iloc[0]
                if bt is None or c.bits > bt.bits: bt = c
        if bt is not None:
            print(f"    [구조 {st:7s}] {str(bt.protein):14s} tm {bt.alntmscore:.3f} "
                  f"fident {bt.fident:.3f}  {HDR.get(str(bt.protein),('',''))[1][:44]}")
        else:
            print(f"    [구조 {st:7s}] 없음")
        rows.append({"query": pid, "query_desc": qd, "strain": st,
                     "seq_hit": best_s.target if best_s is not None else None,
                     "seq_fident": round(float(best_s.fident), 3) if best_s is not None else None,
                     "str_hit": str(bt.protein) if bt is not None else None,
                     "str_tm": round(float(bt.alntmscore), 3) if bt is not None else None,
                     "str_desc": HDR.get(str(bt.protein), ("", ""))[1][:60] if bt is not None else ""})
C41 = pd.DataFrame(rows)
C41.to_csv(DIR["table"]/"coo_operon_homologs.csv", index=False, encoding="utf-8-sig")
print(f"\n저장: {DIR['table']/'coo_operon_homologs.csv'}")
print("\n※ 읽는 법")
print("  BL21 에는 있고 MG1655 에는 없는 행이 있으면 그것이 이 프로젝트가")
print("  찾던 것이다. 둘 다 있으면 발현량 가설로 넘어간다.")
print("  tm >= 0.5 는 같은 폴드, fident 는 서열 유사도다. 둘이 어긋나면")
print("  '기능은 같고 서열은 먼' 경우이므로 오히려 흥미롭다.")
```

---

## CELL 50 — ★ 한 번에 실행 (이량체 대조군 + coo 오페론 + HypA 폴드 조사)

```python
# =============================================================================
# CELL 50 | CELL 24d · 40 · 41 · 42 를 순서대로 한 번에 돌린다
#   GPU 경합이 없도록 순서를 짰다.
#     STEP 1  RF2-PPI tandem  GPU 전경 ~3분   (먼저 끝내고 GPU 를 비운다)
#     STEP 2  Boltz 이량체     GPU 백그라운드   (뒤 단계가 도는 동안 돈다)
#     STEP 3  Y19 coo 오페론   CPU
#     STEP 4  오페론 상동체    CPU  mmseqs + foldseek
#     STEP 5  HypA 폴드 조사   CPU  foldseek
#     STEP 6  종합
#   각 단계는 실패해도 다음으로 넘어간다. 마지막에 무엇이 되고 무엇이 안 됐는지 찍는다.
#   선행: CELL 01 (경로·헬퍼) · CELL 07 (BAIT_ALL)
# =============================================================================
need((have("ASSET", "DIR", "BAIT_ALL", "BAIT_KEY"), "CELL 01 과 CELL 07 을 먼저 돌릴 것"))
import yaml as _y, statistics as _st
RESULT = {}
def step(tag):
    print("\n" + "=" * 110); print(f"### {tag}"); print("=" * 110)

# ---------- 준비 ----------
FS   = shutil.which("foldseek") or globals().get("FOLDSEEK_BIN", "foldseek")
FAA  = {"BL21": ASSET["faa_bl21"], "MG1655": ASSET["faa_mg1655"], "Y19": ASSET["faa_y19"]}
SEQ, HDR = {}, {}
for _st_, _fa in FAA.items():
    nm, buf = None, []
    for l in open(_fa, errors="ignore"):
        if l.startswith(">"):
            if nm: SEQ[nm] = "".join(buf)
            t = l[1:].rstrip(); nm = t.split()[0]; buf = []
            HDR[nm] = (_st_, t[len(nm):].strip())
        else: buf.append(l.strip())
    if nm: SEQ[nm] = "".join(buf)
STRUCT = {}
for r in [TOOLS/"database"/"bacteriaDB", TOOLS/"database"]:
    if Path(r).exists():
        for d in sorted(Path(r).glob("structures_*")):
            if d.is_dir(): STRUCT[d.name] = d
_xwp = DIR["table"]/"id_crosswalk_struct_to_genbank.csv"
_xw  = pd.read_csv(_xwp) if _xwp.exists() else pd.DataFrame()
XWM  = {Path(str(k)).stem: str(v).split(",")[0]
        for k, v in (_xw.dropna(subset=["protein"]).set_index("tid")["protein"]
                     .to_dict().items() if len(_xw) else [])}
REV  = {v: k for k, v in XWM.items()}
_tmp = DIR["tmp"]/"c50"; _tmp.mkdir(parents=True, exist_ok=True)
FSD  = DIR["folddisco"].parent/"foldseek"; FSD.mkdir(exist_ok=True)
def struct_of(pid):
    stem = REV.get(pid)
    if not stem: return None
    for sdir in STRUCT.values():
        for ext in ("cif", "pdb"):
            f = next(iter(glob.glob(str(sdir/f"{stem}*.{ext}"))), None)
            if f: return f
    return None
bseq  = BAIT_ALL[BAIT_KEY]
CTRLK = next((k for k in BAIT_ALL if k.upper().startswith("COOC")), None)
cseq  = BAIT_ALL.get(CTRLK, "")
print(f"준비 완료 — 구조 디렉터리 {list(STRUCT)}  crosswalk {len(XWM)}건")
print(f"{BAIT_KEY} {len(bseq)}aa · {CTRLK} {len(cseq)}aa")

# ================= STEP 1-2. GPU 작업 준비 =================
step("STEP 1-2  GPU 작업 준비 — 입력만 만들고 실행은 백그라운드로 보낸다")
print("RF2-PPI 와 Boltz 를 한 스크립트에 순차로 묶는다. 둘 다 GPU 를 쓰므로")
print("동시에 띄우면 서로 죽인다. 순차면 경합이 없고, 커널도 안 묶인다.")
GPU_SCRIPT, n_boltz = [], 0
_GID = globals().get("GPU_ID", 1)

# --- (a) RF2-PPI tandem a3m 준비 (CPU, 몇 초) ---
try:
    src = DIR["paired"] / f"{CTRLK}__{BAIT_KEY}.a3m"
    if not src.exists():
        raise FileNotFoundError(src)
    LC, LB = len(cseq), len(bseq)
    names, seqs, nm, buf = [], [], None, []
    for l in open(src, errors="ignore"):
        if l.startswith(">"):
            if nm is not None:
                names.append(nm); seqs.append("".join(buf))
            nm, buf = l[1:].strip(), []
        else:
            buf.append(l.strip())
    if nm is not None:
        names.append(nm); seqs.append("".join(buf))
    assert all(len(s) == LC + LB for s in seqs), "a3m 길이 불일치"
    vdir = DIR["paired"] / "dimer"; vdir.mkdir(exist_ok=True)
    raw = vdir / "tandem.raw.a3m"
    raw.write_text("".join(">%s\n%s%s%s%s\n" % (n_, s[:LC], s[:LC], s[LC:], s[LC:])
                           for n_, s in zip(names, seqs)))
    out = vdir / "tandem.a3m"
    try:
        subprocess.run(["hhfilter", "-i", str(raw), "-o", str(out), "-id", "90",
                        "-M", "first"], check=True, capture_output=True)
    except Exception:
        shutil.copy(raw, out)
    _d = sum(1 for l in open(out) if l.startswith(">"))
    (vdir / "input_tandem").write_text("%s %d\n" % (out, 2 * LC))
    print("  (a) tandem a3m: chainA=%d chainB=%d 깊이 %d" % (2 * LC, 2 * LB, _d))
    print("      ※ MSA 열 복제라 공진화 정보량은 안 는다. 기대는 낮지만 비용도 낮다.")
    GPU_SCRIPT.append(
        'echo "=== RF2-PPI tandem ==="\n'
        'cd "%s"\n'
        'for rep in 1 2 3; do\n'
        '  cp input_tandem in_t${rep}\n'
        '  CUDA_VISIBLE_DEVICES=%s python "%s/src/predict_list_PPI.py" '
        '-list_fn in_t${rep} -model_file "%s/src/models/RF2-PPI.pt" || true\n'
        'done\n'
        'echo DONE_rf2_tandem\n' % (vdir, _GID, RF2PPI_DIR, RF2PPI_DIR))
except Exception as e:
    print("  (a) RF2-PPI 준비 실패: %s: %s" % (type(e).__name__, e))

# --- (b) Boltz 이량체 YAML 준비 (CPU, 즉시) ---
try:
    MAXMSA = globals().get("BOLTZ_MAXMSA", 2048)
    din = DIR["boltz"] / "inputs_dimer"; din.mkdir(parents=True, exist_ok=True)
    for f in din.glob("*.yaml"):
        f.unlink()

    def wr(name, chains, ni=True):
        tot = sum(len(s) for _, s in chains)
        if tot > 1830:
            print("      제외 %s (%d 토큰 > 1830)" % (name, tot)); return 0
        sq = [{"protein": {"id": c, "sequence": s}} for c, s in chains]
        if ni:
            sq.append({"ligand": {"id": "Z", "ccd": "NI"}})
        (din / (name + ".yaml")).write_text(
            _y.safe_dump({"version": 1, "sequences": sq}, sort_keys=False))
        print("      %s: %d사슬 %d 토큰" % (name, len(chains), tot)); return 1

    print("  (b) Boltz 입력")
    n_boltz += wr("CTRL-dimer2x2",  [("A", bseq), ("B", bseq), ("C", cseq), ("D", cseq)])
    n_boltz += wr("CTRL-mono1x1",   [("A", bseq), ("C", cseq)])
    n_boltz += wr("CTRL-codhdimer", [("A", bseq), ("B", bseq), ("C", cseq)])
    n_boltz += wr("APO-codhdimer",  [("A", bseq), ("B", bseq)])
    for pid in ["QJZ12568.1"]:
        if pid in SEQ:
            n_boltz += wr("CAND-%s-codhdimer" % pid,
                          [("A", bseq), ("B", bseq), ("C", SEQ[pid])])
    if n_boltz:
        shutil.rmtree(DIR["boltz"] / "out_dimer", ignore_errors=True)
        GPU_SCRIPT.append(
            'echo "=== Boltz dimer ==="\n'
            'conda activate %s\n'
            'cd "%s"\n'
            'export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True\n'
            'CUDA_VISIBLE_DEVICES=%s boltz predict inputs_dimer '
            '--out_dir out_dimer --use_msa_server --max_msa_seqs %s '
            '--recycling_steps 3 --diffusion_samples 5 '
            '--output_format mmcif --num_workers 2\n'
            'echo DONE_boltz_dimer\n' % (CONDA_ENV_BOLTZ, DIR["boltz"], _GID, MAXMSA))
except Exception as e:
    print("  (b) Boltz 준비 실패: %s: %s" % (type(e).__name__, e))

# --- (c) 한 스크립트로 묶어 백그라운드 ---
_busy = subprocess.run("pgrep -f '[b]oltz predict|[p]redict_list_PPI'", shell=True,
                       capture_output=True, text=True).stdout.split()
if _busy:
    print("\n⚠ GPU 작업이 이미 돈다 (pid %s). 새로 안 띄웠다." % " ".join(_busy))
elif GPU_SCRIPT:
    sh_bg("part8_gpu", "\n".join(GPU_SCRIPT), env=CONDA_ENV_RF2)
    RESULT["gpu_bg"] = "RF2 tandem + Boltz %d개" % n_boltz
    print("\n>>> GPU 작업은 백그라운드. 아래 STEP 들은 CPU 라 그 위에서 같이 돈다.")
else:
    print("\n띄울 GPU 작업이 없다.")

# ================= STEP 3. Y19 coo 오페론 =================
step("STEP 3  Y19 coo 오페론 — 양성대조군 세트 확보")
print("문헌(Merrouch 2018 Fig.2): CaY19 는 cooS 뒤에 CooC·CooT·CooJ·HypB 를 갖는다.")
print("크기 기준  CooS ~630 / CooC ~270 / CooT ~65 / CooJ ~115 / HypB ~290 / CooF ~200 aa")
COO_SET = []
try:
    y19 = {k: v for k, v in SEQ.items() if HDR[k][0] == "Y19"}
    AHZ = sorted(k for k in y19 if k.startswith("AHZ"))
    print(f"\nAHZ 계열 {len(AHZ)}개 (coo 오페론 제출본)")
    for k in AHZ:
        print(f"  {k}  {len(y19[k]):5d} aa  {HDR[k][1][:66]}")
    PAT = r"\bCoo[CSTJF]\b|carbon monoxide dehydrogenase|hydrogenase nickel|\bHypB\b|nickel"
    byname = sorted(k for k in y19 if re.search(PAT, HDR[k][1], re.I))
    print(f"\n이름으로 잡힌 것 {len(byname)}개")
    for k in byname:
        print(f"  {k}  {len(y19[k]):5d} aa  {HDR[k][1][:66]}")
    SIZE = {"CooS": (550, 720), "CooC": (220, 320), "CooT": (45, 95),
            "CooJ": (90, 155), "HypB": (240, 340), "CooF": (150, 250)}
    print("\n크기 기반 배정 후보 (AHZ 안에서)")
    for role, (lo, hi) in SIZE.items():
        c = [k for k in AHZ if lo <= len(y19[k]) <= hi]
        print(f"  {role:6s} {lo:>3}-{hi:<3} aa : {c if c else '없음'}")
    COO_SET = sorted(set(AHZ) | set(byname))
    qf = DIR["seq"]/"y19_coo_operon.faa"
    qf.write_text("".join(f">{k}\n{y19[k]}\n" for k in COO_SET))
    RESULT["coo_set"] = len(COO_SET)
    print(f"\n대조군 세트 {len(COO_SET)}개 -> {qf.name}")
    if len(AHZ) <= 2:
        print("※ AHZ 가 2개뿐이다 — 오페론 나머지는 AKE/AGE 로 들어가 있다.")
        print("  STEP 4 의 서열·구조 검색이 이름과 무관하게 찾아낸다.")
except Exception as e:
    print(f"[STEP 3 실패] {type(e).__name__}: {e}")

# ================= STEP 4. 오페론 상동체 탐색 =================
step("STEP 4  coo 오페론 상동체를 BL21·MG1655 에서 — mmseqs + foldseek")
print("Folddisco 두-Cys 질의는 Y19 CooC 를 놓쳤고 Foldseek 은 tm 0.938 로 찾았다.")
print("그래서 이번엔 Foldseek 을 쓰고, 짧은 CooT(65aa) 때문에 서열 검색도 같이 돌린다.")
C41 = pd.DataFrame()
try:
    qf = DIR["seq"]/"y19_coo_operon.faa"
    SFMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits"
    SEQHIT = {}
    for st_ in ["BL21", "MG1655"]:
        o = DIR["search"]/f"cooop_vs_{st_}.m8"
        if not (o.exists() and o.stat().st_size):
            sh(f'mmseqs easy-search "{qf}" "{FAA[st_]}" "{o}" "{_tmp}" '
               f'--format-output "{SFMT}" -e 1e-2 -s 7.5 --threads {min(THREADS,16)} -v 1',
               check=False)
        if o.exists() and o.stat().st_size:
            SEQHIT[st_] = pd.read_csv(o, sep="\t", names=SFMT.split(",")) \
                            .sort_values("bits", ascending=False)
    FFMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits,prob,alntmscore,lddt"
    STRHIT = {}
    for pid in COO_SET:
        sp = struct_of(pid)
        if not sp: continue
        for sname, sdir in STRUCT.items():
            o = FSD/f"coo_{pid}_vs_{sname}.m8"
            if not (o.exists() and o.stat().st_size):
                sh(f'"{FS}" easy-search "{sp}" "{sdir}" "{o}" "{_tmp}" '
                   f'--format-output "{FFMT}" -e 10 --max-seqs 500 '
                   f'--threads {min(THREADS,16)} -v 1', check=False)
            if o.exists() and o.stat().st_size:
                d = pd.read_csv(o, sep="\t", names=FFMT.split(","))
                d["protein"] = d.target.apply(lambda x: XWM.get(Path(str(x)).stem))
                STRHIT[(pid, sname)] = d
    rows = []
    for pid in COO_SET:
        print(f"\n--- {pid}  {len(SEQ.get(pid,'')):4d} aa  {HDR.get(pid,('',''))[1][:52]}")
        for st_ in ["BL21", "MG1655"]:
            bs = None
            d = SEQHIT.get(st_)
            if d is not None and len(d):
                m = d[d["query"] == pid]
                if len(m): bs = m.iloc[0]
            print(f"   [서열 {st_:7s}] " + (f"{bs.target:14s} fident {bs.fident:.3f} "
                  f"qcov {bs.qcov:.2f}  {HDR.get(bs.target,('',''))[1][:40]}"
                  if bs is not None else "없음"))
            bt = None
            for (q, sn), dd in STRHIT.items():
                if q != pid: continue
                z = dd[dd.protein.map(lambda k: HDR.get(str(k), ("?",))[0]) == st_]
                z = z[z.protein.astype(str) != pid]
                if len(z):
                    c = z.sort_values("bits", ascending=False).iloc[0]
                    if bt is None or c.bits > bt.bits: bt = c
            print(f"   [구조 {st_:7s}] " + (f"{str(bt.protein):14s} tm {bt.alntmscore:.3f} "
                  f"fident {bt.fident:.3f}  {HDR.get(str(bt.protein),('',''))[1][:40]}"
                  if bt is not None else "없음"))
            rows.append({"query": pid, "query_desc": HDR.get(pid, ("", ""))[1][:50],
                         "strain": st_,
                         "seq_hit": bs.target if bs is not None else None,
                         "seq_fident": round(float(bs.fident), 3) if bs is not None else None,
                         "str_hit": str(bt.protein) if bt is not None else None,
                         "str_tm": round(float(bt.alntmscore), 3) if bt is not None else None,
                         "str_desc": HDR.get(str(bt.protein), ("", ""))[1][:55] if bt is not None else ""})
    C41 = pd.DataFrame(rows)
    C41.to_csv(DIR["table"]/"coo_operon_homologs.csv", index=False, encoding="utf-8-sig")
    RESULT["coo_homologs"] = len(C41)
    # BL21 에만 있는 것
    piv = C41.pivot_table(index="query", columns="strain",
                          values="str_tm", aggfunc="max")
    if "BL21" in piv and "MG1655" in piv:
        only = piv[(piv.BL21.notna()) & (piv.MG1655.isna() |
                   (piv.BL21 - piv.MG1655.fillna(0) > 0.15))]
        print(f"\n>>> BL21 쪽이 뚜렷이 가까운 질의: {len(only)}개")
        if len(only): print(only.round(3).to_string())
except Exception as e:
    print(f"[STEP 4 실패] {type(e).__name__}: {e}")

# ================= STEP 5. HypA 폴드 전수 조사 =================
step("STEP 5  HypA 폴드로 세 균주 전수 조사")
print("QJZ12568.1 이 HypA 폴드(tm 0.66~0.68)였다. 같은 폴드가 몇 개나 있고")
print("그중 BL21 특이가 정말 얘 하나뿐인지 본다.")
HYPA = pd.DataFrame()
try:
    seeds = [k for k, v in HDR.items()
             if v[0] == "BL21" and re.search(r"HypA|nickel metallochaperone", v[1], re.I)]
    seeds = sorted(set(seeds) | {"QJZ12568.1"})
    print(f"질의 seed: {seeds}")
    FFMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits,prob,alntmscore,lddt"
    rows = []
    for pid in seeds:
        sp = struct_of(pid)
        if not sp:
            print(f"  {pid}: 구조 없음"); continue
        for sname, sdir in STRUCT.items():
            o = FSD/f"hypa_{pid}_vs_{sname}.m8"
            if not (o.exists() and o.stat().st_size):
                sh(f'"{FS}" easy-search "{sp}" "{sdir}" "{o}" "{_tmp}" '
                   f'--format-output "{FFMT}" -e 10 --max-seqs 500 '
                   f'--threads {min(THREADS,16)} -v 1', check=False)
            if not (o.exists() and o.stat().st_size): continue
            d = pd.read_csv(o, sep="\t", names=FFMT.split(","))
            d["protein"] = d.target.apply(lambda x: XWM.get(Path(str(x)).stem))
            d = d[d.protein.notna() & (d.alntmscore >= 0.50)]
            for r in d.itertuples(index=False):
                rows.append({"seed": pid, "protein": str(r.protein),
                             "strain": HDR.get(str(r.protein), ("?",))[0],
                             "tm": round(float(r.alntmscore), 3),
                             "fident": round(float(r.fident), 3),
                             "desc": HDR.get(str(r.protein), ("", ""))[1][:58]})
    HYPA = (pd.DataFrame(rows).sort_values("tm", ascending=False)
              .drop_duplicates(["seed", "protein"]))
    if len(HYPA):
        HYPA.to_csv(DIR["table"]/"hypA_fold_sweep.csv", index=False, encoding="utf-8-sig")
        RESULT["hypA_fold"] = len(HYPA)
        print(f"\ntm >= 0.50 히트 {len(HYPA)}건")
        print(HYPA.groupby(["seed", "strain"]).size().to_string())
        pd.set_option("display.max_rows", None); pd.set_option("display.width", 220)
        print("\n" + HYPA.to_string(index=False))
        # BL21 히트의 MG1655 상동체 유무
        bl = sorted(set(HYPA[HYPA.strain == "BL21"].protein))
        if bl:
            q = DIR["seq"]/"hypa_fold_bl21.faa"
            q.write_text("".join(f">{k}\n{SEQ[k]}\n" for k in bl if k in SEQ))
            o = DIR["search"]/"hypafold_vs_MG1655.m8"
            SFMT = "query,target,fident,alnlen,qcov,tcov,evalue,bits"
            if not (o.exists() and o.stat().st_size):
                sh(f'mmseqs easy-search "{q}" "{FAA["MG1655"]}" "{o}" "{_tmp}" '
                   f'--format-output "{SFMT}" -e 1e-3 -s 7.5 '
                   f'--threads {min(THREADS,16)} -v 1', check=False)
            hit = set()
            if o.exists() and o.stat().st_size:
                dd = pd.read_csv(o, sep="\t", names=SFMT.split(","))
                hit = set(dd[dd.qcov >= 0.70]["query"])
            print(f"\n>>> HypA 폴드 BL21 단백질 {len(bl)}개 중 MG1655 에 상동체 없음:")
            for k in bl:
                mark = "없음 ★" if k not in hit else "있음"
                print(f"     {mark:6s}  {k}  {HDR.get(k,('',''))[1][:60]}")
            RESULT["hypA_bl21_only"] = [k for k in bl if k not in hit]
    else:
        print("tm >= 0.50 히트 없음")
except Exception as e:
    print(f"[STEP 5 실패] {type(e).__name__}: {e}")

# ================= STEP 6. 종합 =================
step("STEP 6  종합")
for k, v in RESULT.items():
    print(f"  {k:18s} : {v}")
print("\n저장된 표")
for f in ["coo_operon_homologs.csv", "hypA_fold_sweep.csv"]:
    p = DIR["table"]/f
    print(f"  {'O' if p.exists() else 'X'} {p}")
print("\n### GPU 백그라운드 (part8_gpu) — 지금도 돌고 있다")
print("  진행확인:  bg_tail('part8_gpu')  또는 CELL 15b")
print("  순서   :  RF2-PPI tandem 3회 -> Boltz dimer 5쌍 x 5표본")
print("  다 끝나면 로그에 DONE_rf2_tandem 과 DONE_boltz_dimer 가 찍힌다.")
print("\n### 끝난 뒤 읽는 법")
print("  (1) RF2-PPI tandem")
print("      result/paired/dimer/in_t*.log 의 2번째 열이 점수다.")
print("      단량체판 0.251 (CELL 24c v0_base) 과 비교한다.")
print("  (2) Boltz dimer — CELL 28a-13 의 파싱 방식을 out_dimer 에 적용")
print("      ★ CTRL-dimer2x2 의 ipTM 이 CTRL-mono1x1(0.300)보다 높은가")
print("        높으면 지금까지의 Track B 348개 구조 전체를 이량체로 다시 해야 한다.")
print("      ★ APO-codhdimer 에서 Ni 이 C-cluster 자리에 들어가는가")
print("        75 C 무단백 조립 결과와 직접 대응하는 항목이다.")
```

---

## CELL 37 — 최종 리포트 + 남은 TODO

```python
# =============================================================================
# CELL 37 | 최종 요약 + 남은 TODO
# =============================================================================
print("### 산출물 ###")
for f in sorted(DIR["table"].glob("*.csv")):
    print(f"  {f.name:44s} {f.stat().st_size/1e3:8.1f} KB")

print(f"""
### 요약 ###
  Track A (RF2-PPI)  : {len(TRACK_A) if 'TRACK_A' in dir() else 0}개 대상, """
      f"""depth={BAIT_DEPTH.get(BAIT_KEY,'?') if 'BAIT_DEPTH' in dir() else '?'}, """
      f"""진행={'예' if 'TRACK_A_GO' in dir() and TRACK_A_GO else '아니오'}
  Track B (Boltz-2)  : {len(TRACK_B) if 'TRACK_B' in dir() else 0}개 대상
  Track C (Folddisco): metal={len(ni_df) if 'ni_df' in dir() else 0}, """
      f"""atp={len(atp_df) if 'atp_df' in dir() else 0}
  통합 후보          : {len(M) if 'M' in dir() else 0}개, 상위 20개가 실험 검증 대상

### 다음 단계 (설계문서 §6) ###
  1. 통합 랭킹 상위 20개
  2. BL21(DE3) 단일 KO / CRISPRi -> ChCODH2 활성 소실 여부
  3. ★MG1655 complementation -> 활성 회복 여부 (Keio/ASKA 자원 사용 가능, 가장 빠르고 결정적)
  4. 정제 단백질 재구성: 후보 + apo-ChCODH2 + NiCl2 (+-ATP/GTP)
  5. His-apo-ChCODH2 pulldown + LC-MS/MS (BL21 vs MG1655 lysate 비교)

### 계산과 병행할 실험 (Priority 0) ###
  열처리 / Proteinase K / 분자량 분획 / nucleotide 의존성 / 활성 유도 분획 + LC-MS/MS
  ※ 5번(활성 유도 분획화)이 답에 가장 빨리 도달하는 경로다. 병렬로 돌릴 것.
  ※ mrp / yeiR / yjiA 단일 KO 미실시 — 즉시 시험 가치 높음
  ※ hypA·hybF, hypC·hybG 이중 KO 로 paralog 보완 배제

### 남은 TODO (인계문서 §7) ###
  [!] BL21 proteome ID: faa(QZI…)와 구조 DB(QJZ…)의 ID 계열이 달랐다 (교집합 0).
      -> 구조 DB 기록으로 faa 재구축 후 CELL 08~ 재실행 필요
  [x] Folddisco ATP-binding motif residue 지정 -> A12,A13,A14 (Walker A) 확정
  [x] 3kji.pdb 체인 ID 확정 -> 체인 A, A112/A114 = CYS 확인
  [ ] crosswalk 미매칭분 UniProt idmapping 으로 보완 (CELL 30)
  [x] apt-mark hold nvidia-* -> 적용 완료 (CELL 02 에서 상태 확인)
  [ ] 정량 프로테오믹스 병행 — 서열 동일·발현량만 다를 경우의 보험
""")
```

---

## 부록 — 실행 체크리스트

```
[ ] CELL 01-02  CONFIG / 자산 점검
[ ] CELL 03-04  RF2-PPI 설치 + 예제 검증           ★체크포인트 1
[ ] CELL 05-07  ChCODH2(A559W) 서열 입력           ← 없으면 아무것도 안 돌아감
[ ] CELL 08-12  프로테옴 로드 + 4-category 분류 + Track 배정
[ ] CELL 13-15  Bacteria reference DB              (~4h, 백그라운드)
[ ] CELL 16-17  Feasibility gate                   ★체크포인트 2
[ ] CELL 18-21  prey MSA + paired MSA              (~10h, 백그라운드)
[ ] CELL 22-24  RF2-PPI x3(또는 5) replicate       (GPU)
[ ] CELL 25-28  Boltz-2 Track B                    (GPU, 병렬 가능)
[ ] CELL 29-31  Folddisco + ID crosswalk
[ ] CELL 32-37  통합 랭킹 + 리포트
```

**병렬화**: CELL 13-21(CPU)이 도는 동안 CELL 25-28(Boltz-2 Track B)을 다른 GPU에서 먼저 시작할 수 있다.
Track B는 MSA 서버를 쓰므로 Part 2가 필요 없다.
