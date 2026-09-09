# ChCODH2 Ni-insertase 스크리닝 — Antigravity 노트북 셀 모음

인계문서(`HANDOVER_ChCODH2_screen.md`, `HANDOVER_bacteriaDB.md`, `ChCODH2_Ni_insertase_in_silico_design.md`)와
원본 노트북(`RF2PPI_ChCODH2_screen.ipynb`)을 합쳐, **epel-af2 서버의 실제 경로와 이미 구축된 자산**을 반영한 셀 모음.

사용법: 아래 `CELL 00` ~ `CELL 33`을 **위에서부터 하나씩 복사해 Antigravity 노트북 셀에 붙여넣고 실행**한다.
전부 파이썬 셀이다(`%%bash` 매직 없음 — Antigravity/커널 환경에 상관없이 돌아가도록 `subprocess`로 통일).
셀 사이에 변수가 이어지므로 순서대로 실행해야 한다.

원본 노트북과 달라진 핵심 4가지:
1. `BASE`를 `/data/chcodh2_ppi` → `/mnt/af2results/mingyu/chcodh2_ppi` 로 변경 (서버에 `/data` 없음)
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
[ChCODH2 Ni-insertase in silico 스크리닝 — 셀 실행 순서]

 CELL 01  CONFIG (경로/스레드/GPU/헬퍼 함수)            <- 여기만 본인 환경에 맞게 수정
 CELL 02  기존 자산 + 환경 점검 (GPU/디스크/인덱스/서열)

 [Part 0] 설치
 CELL 03  RF2-PPI 설치 스크립트 생성 (터미널에서 실행)
 CELL 04  설치 검증 — 예제 재현            ★체크포인트 1

 [Part 1] 서열 확보
 CELL 05  bait 서열 입력 (ChCODH2 A559W)   ★없으면 아무것도 안 돌아감
 CELL 06  CooC1 서열을 3kji.pdb에서 추출 (양성대조군)
 CELL 07  bait 검증 + segment bait 생성 + baits.fasta 저장

 [Part 1-2/1-3] 프로테옴 · 차등 분류  (기존 자산 재사용)
 CELL 08  기존 프로테옴 서열 로드 (BL21/Y19/MG1655, GenBank ID)
 CELL 09  기존 BLASTp 4-category 분류 로드 + 값 검증 (2171/1575/125/217)
 CELL 10  (fallback) mmseqs easy-search 재실행
 CELL 11  (fallback) m8 -> 4-category 분류
 CELL 12  Track A / Track B 확정 및 저장

 [Part 2] paired MSA용 세균 레퍼런스 DB
 CELL 13  Bacteria reference proteome 다운로드   (백그라운드, ~1h)
 CELL 14  병합 + mmseqs createdb                 (백그라운드, ~2h)
 CELL 15  백그라운드 작업 모니터 (수시 실행)

 [Part 3] Feasibility gate
 CELL 16  bait homolog 검색
 CELL 17  깊이 판정 — Track A 진행 여부 결정      ★체크포인트 2

 [Part 4-6] Track A = RF2-PPI (공진화)
 CELL 18  Track A prey FASTA
 CELL 19  prey homolog 검색                       (백그라운드, 4~10h, 최장 단계)
 CELL 20  organism별 best-hit 정리
 CELL 21  paired MSA 생성 + hhfilter
 CELL 22  RF2-PPI x3 replicate 실행               (GPU)
 CELL 23  결과 집계 (mean/sd 랭킹)
 CELL 24  양성대조군 검증 + depth-score 상관 진단

 [Part 7] Track B = Boltz-2 (구조 co-folding, 공진화 불필요) ★주력 가능성 높음
 CELL 25  Boltz-2 설치 스크립트
 CELL 26  입력 YAML 생성 (Ni을 ligand CCD:NI로 명시)
 CELL 27  boltz predict 실행                      (GPU)
 CELL 28  결과 파싱 (ipTM / ligand_ipTM)

 [Part 8] Track C 연계 + 통합
 CELL 29  Folddisco 검색 (metal motif / ATP motif)
 CELL 30  ID crosswalk: 구조 tid(UniProt/UniParc) -> GenBank protein ID
 CELL 31  Folddisco 결과 -> folddisco_*_motif.csv 변환
 CELL 32  Track A + B + C 통합 랭킹
 CELL 33  최종 리포트 + 남은 TODO

[지금 없는 것 = 블로커]
  1. ChCODH2(A559W) native 서열      -> CELL 05 (필수, 실험에 쓰는 그 서열로)
  2. CooC1/2/CooT/CooJ 양성대조군 서열 -> CELL 05-06
  3. RF2-PPI 설치 + 가중치            -> CELL 03
  4. Boltz-2 설치                     -> CELL 25
  5. Bacteria reference proteome DB   -> CELL 13-14
  6. Feasibility gate 결과            -> CELL 16-17

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
```

---

## CELL 01 — CONFIG · 공용 헬퍼

```python
# =============================================================================
# CELL 01 | CONFIG — 경로/스레드/GPU + 백그라운드 실행 헬퍼
#   - 이 셀만 환경에 맞게 고치면 나머지는 그대로 실행 가능
#   - 원본 노트북의 BASE=/data/chcodh2_ppi 는 이 서버에 없어 /mnt/af2results 로 변경
# =============================================================================
import os, re, json, gzip, shutil, subprocess, textwrap, glob, time
from pathlib import Path

# ---------------- 작업 루트 ----------------
BASE      = Path("/mnt/af2results/mingyu/chcodh2_ppi")   # 디스크 여유 있는 곳 (/data 없음)
THREADS   = 24
GPU_ID    = 0
CONDA_ENV_RF2 = "rf2ppi"
CONDA_ENV_BOLTZ = "boltz"

DIR = {k: BASE / k for k in
       ["seq", "db", "search", "msa", "paired", "rf2ppi", "boltz",
        "results", "tmp", "scripts", "logs", "folddisco"]}
for p in DIR.values():
    p.mkdir(parents=True, exist_ok=True)
RF2PPI_DIR = BASE / "RoseTTAFold2-PPI"

# ---------------- 이미 구축된 자산 (인계문서 §1) ----------------
ASSET = {
    # 프로테옴 서열 (GenBank ID 기준)
    "faa_bl21":   Path("/mnt/af2results/mingyu/database/protein_list/bl21_de3_protein.faa"),
    "faa_y19":    Path("/mnt/af2results/mingyu/database/protein_list/y19_protein.faa"),
    "faa_mg1655": Path("/mnt/af2results/mingyu/database/protein_list/mg1655_protein.faa"),
    "tsv_bl21":   Path("/mnt/af2results/mingyu/database/protein_list/bl21_de3_protein_extracted.tsv"),
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
    "BL21": BASE / "external" / "bl21_mg1655_category_detail.tsv",
    "Y19":  BASE / "external" / "y19_mg1655_category_detail.tsv",
}
(BASE / "external").mkdir(exist_ok=True)

# ---------------- 대상 프로테옴 ----------------
# ★ 원본 노트북의 BL21 = UP000002032 는 우리 구조 DB(UP000503272)와 단백질 세트가 어긋난다.
#    Part 8 조인 실패를 막기 위해 구조 DB 쪽 ID로 통일한다.
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
PREY_STRAIN = "BL21DE3"        # prey 프로테옴 (Y19로 바꾸려면 여기만 변경)
ID_SPACE    = "GenBank"        # 이 파이프라인의 표준 ID 체계

os.environ["MMSEQS_NUM_THREADS"] = str(THREADS)

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
    sp = DIR["scripts"] / f"{name}.sh"
    lg = DIR["logs"] / f"{name}.log"
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
    lg = DIR["logs"] / f"{name}.log"
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

def conda_run(env, cmd, cwd=None, check=True):
    """지정 conda 환경에서 포그라운드 실행."""
    full = (f'source "$(conda info --base)/etc/profile.d/conda.sh" && '
            f'conda activate {env} && {cmd}')
    return sh(full, cwd=cwd, check=check)

print("BASE =", BASE)
for k, v in DIR.items():
    print(f"  {k:10s} {v}")
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
print("### 기존 자산 ###")
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

# folddisco 인덱스는 확장자 없는 본체 + .lookup/.offset/.type 4종 세트
print("\n### Folddisco 인덱스 세트 ###")
for k in ["fd_idx_bl21", "fd_idx_y19", "fd_idx_mg"]:
    base = ASSET[k]
    parts = {ext: (Path(str(base) + ext)).exists() for ext in ["", ".lookup", ".offset", ".type"]}
    print(f"  {k:12s} {base.name:16s} " + " ".join(f"{e or 'body'}={'O' if v else 'X'}"
                                                  for e, v in parts.items()))

print("\n### 분류 TSV (없으면 CELL 10-11 fallback) ###")
for k, p in CATEGORY_TSV.items():
    print(f"  [{'O' if p.exists() else 'X'}] {k}: {p}")

print("\n### GPU ###")
sh("nvidia-smi --query-gpu=index,name,driver_version,memory.total,memory.used "
   "--format=csv,noheader || echo 'nvidia-smi 실패 — 드라이버 확인 필요'", check=False)

print("\n### 디스크 (여유 빠듯: /mnt/af2results 93% 사용 중이었음) ###")
sh(f"df -h {BASE} /mnt/af2results 2>/dev/null | sort -u", check=False)

print("\n### 툴 확인 ###")
for t in ["mmseqs", "hhfilter", "aria2c", "conda", "curl"]:
    sh(f"command -v {t} >/dev/null && echo '  [O] {t}' || echo '  [X] {t} (미설치)'", check=False)

if missing:
    print("\n⚠ 없는 자산:", ", ".join(missing))
print("\n⚠ 인계문서 경고: `sudo apt-mark hold nvidia-*` 미적용 — 커널 업데이트 시 드라이버 재발 가능")
```

---

## CELL 03 — Part 0. RF2-PPI 설치 스크립트 생성 (터미널 실행 권장)

```python
# =============================================================================
# CELL 03 | Part 0-1. RF2-PPI 설치 스크립트 생성
#   - conda create 는 노트북 커널보다 터미널이 안정적이라 "스크립트만 만들고" 출력한다
#   - 실행: 터미널에서 bash <출력된 경로>   (약 20~40분)
# =============================================================================
install_sh = f"""
# ---- 1) conda 환경 (RF2-PPI 공식 스펙: python 3.9 / torch 1.12.1+cu113) ----
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -n {CONDA_ENV_RF2} python=3.9
conda activate {CONDA_ENV_RF2}
conda install -y -c conda-forge -c bioconda hhsuite mmseqs2 aria2
pip install numpy==1.21.2 pandas==1.5.3 biopython==1.79 scipy==1.7.1 einops
pip install torch==1.12.1+cu113 -f https://download.pytorch.org/whl/torch_stable.html
# 드라이버가 CUDA 13 (580.173.02) 이라 cu113 wheel 이 안 뜨면 아래로 교체:
# pip install torch==2.1.2 --index-url https://download.pytorch.org/whl/cu118
pip install jupyter ipykernel matplotlib seaborn tqdm pyyaml
python -m ipykernel install --user --name {CONDA_ENV_RF2} --display-name "Python ({CONDA_ENV_RF2})"

# ---- 2) RF2-PPI 코드 + 가중치 ----
cd "{BASE}"
[ -d RoseTTAFold2-PPI ] || git clone https://github.com/CongLabCode/RoseTTAFold2-PPI.git
cd RoseTTAFold2-PPI/src/models
[ -f RF2-PPI.pt ] || wget --no-check-certificate https://conglab.swmed.edu/humanPPI/downloads/RF2-PPI.pt

# ---- 3) 예제 실행 (검증은 CELL 04 에서) ----
cd "{BASE}/RoseTTAFold2-PPI/examples"
python ../src/predict_list_PPI.py -list_fn segment_pairs_input \\
       -model_file ../src/models/RF2-PPI.pt
echo "설치 스크립트 완료"
"""
p = DIR["scripts"] / "install_rf2ppi.sh"
p.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + textwrap.dedent(install_sh))
p.chmod(0o755)
print("터미널에서 실행하세요:\n")
print(f"    bash {p}\n")
print("-" * 70)
print(p.read_text())
```

---

## CELL 04 — Part 0. 설치 검증 (★체크포인트 1)

```python
# =============================================================================
# CELL 04 | Part 0-2. 설치 검증 — 예제 결과가 expected_output 과 일치하는가
#   ★체크포인트 1: 소수 둘째 자리까지 비슷해야 설치 성공.
#     크게 다르면 CUDA/torch 버전 문제 -> 여기서 잡고 넘어간다 (뒤에서 못 잡는다)
# =============================================================================
ex = RF2PPI_DIR / "examples"
assert ex.exists(), f"RF2-PPI 미설치. CELL 03 스크립트를 터미널에서 먼저 실행하세요: {ex}"

log = ex / "segment_pairs_input.log"
if not log.exists():
    conda_run(CONDA_ENV_RF2,
              f'cd "{ex}" && python ../src/predict_list_PPI.py '
              f'-list_fn segment_pairs_input -model_file ../src/models/RF2-PPI.pt')

import pandas as pd
got = pd.read_csv(log, sep=r"\s+", names=["msa", "prob", "sec"])
exp_files = list((ex / "expected_output").glob("*.log"))
print("=== 예제 실행 결과 ===")
print(got.head(10).to_string(index=False))

if exp_files:
    expd = pd.read_csv(exp_files[0], sep=r"\s+", names=["msa", "prob", "sec"])
    m = got.merge(expd, on="msa", suffixes=("_got", "_exp"))
    m["diff"] = (m.prob_got - m.prob_exp).abs()
    print("\n=== expected_output 대조 ===")
    print(m[["msa", "prob_got", "prob_exp", "diff"]].to_string(index=False))
    ok = (m["diff"] < 0.05).all()
    print(f"\n★체크포인트 1: {'통과' if ok else '실패 — CUDA/torch 버전 점검 필요'}"
          f" (최대 편차 {m['diff'].max():.4f})")
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
BAIT_SEQ = {
    # "ChCODH2_A559W": "MSEK...",     # <<<< 여기에 붙여넣기 (필수)
    # "CooC1": "...", "CooC2": "...", "CooT": "...", "CooJ": "...",   # 양성대조군
}

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

if not BAIT_SEQ:
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
        print("UniProt 조회 실패:", e, "\n-> BAIT_SEQ 를 수동으로 채우세요.")
else:
    print("BAIT_SEQ 입력됨:", {k: len(v) for k, v in BAIT_SEQ.items()})
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

# primary bait 지정 (ChCODH2)
cands = [k for k in BAIT_SEQ if "CODH" in k.upper() or "A559W" in k.upper()]
assert cands, "ChCODH2 bait 를 찾을 수 없습니다. 이름에 'CODH' 또는 'A559W' 를 넣으세요."
BAIT_KEY = cands[0]
print("\nprimary bait =", BAIT_KEY)

# A559W 변이 확인
s = BAIT_SEQ[BAIT_KEY]
if len(s) >= 559:
    print(f"  559번 잔기 = {s[558]}   (A559W 변이체라면 'W' 여야 함)")
else:
    print(f"  ⚠ 길이 {len(s)} < 559 — 서열/번호 체계 재확인 필요")

# ---- segment bait: C-cluster 주변(559 중심) + 전체 타일링 ----
SEG_WIN, SEG_STEP = 200, 100
SEGMENTS = {}
c = 559
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
#   - 원본 노트북은 UniProt 에서 새로 받지만, 우리 구조 DB/분류가 GenBank ID 기준이라
#     여기서 UniProt 을 새로 받으면 Part 8 에서 조인이 깨진다. ID 체계를 GenBank 로 통일.
#   - 기대값: BL21 4,088 / Y19 5,325
# =============================================================================
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

exp = {"BL21DE3": 4088, "Y19": 5325}
for k, n in exp.items():
    if "seqs" in PROTEOMES.get(k, {}):
        got = len(PROTEOMES[k]["seqs"])
        print(f"  {k}: {got} (기대 {n}) {'OK' if got == n else '⚠ 불일치 — 파일 버전 확인'}")
```

---

## CELL 09 — Part 1-3a. 기존 4-category 분류 로드 + 값 검증

```python
# =============================================================================
# CELL 09 | Part 1-3a. 기존 BLASTp 4-category 분류 로드 (재계산 스킵)
#   - 기대값 BL21->MG1655: identical 2171 / high 1575 / low 125 / specific 217 (합 4088)
#   - 파일이 없으면 cls_BL = None 이 되고, CELL 10-11 fallback 으로 넘어간다
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

cls_BL = None
p = CATEGORY_TSV["BL21"]
if p.exists():
    cls_BL = load_category_tsv(p)
    vc = cls_BL.category.value_counts()
    print("=== BL21 -> MG1655 (로드) ===")
    print(vc.to_string())
    print("총:", len(cls_BL))
    expected = {"identical": 2171, "high-similarity": 1575,
                "low-similarity": 125, "strain-specific": 217}
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
#   - CELL 09 에서 cls_BL 이 로드됐으면 이 셀은 건너뛴다
# =============================================================================
if cls_BL is not None:
    print("CELL 09 에서 이미 로드됨 — 이 셀은 건너뜁니다.")
else:
    fa_bl = PROTEOMES["BL21DE3"]["faa"]
    fa_mg = PROTEOMES["MG1655"]["faa"]
    fmt = "query,target,fident,alnlen,evalue,bits,qlen,tlen,qcov,tcov"
    script = f"""
    cd "{BASE}"
    mmseqs easy-search "{fa_bl}" "{fa_mg}" search/BL21_vs_MG.m8 tmp/bm \\
      --format-output "{fmt}" -s 7.5 -e 1e-5 --max-seqs 5 --threads {THREADS}
    mmseqs easy-search "{fa_mg}" "{fa_bl}" search/MG_vs_BL21.m8 tmp/mb \\
      --format-output "{fmt}" -s 7.5 -e 1e-5 --max-seqs 5 --threads {THREADS}
    wc -l search/BL21_vs_MG.m8 search/MG_vs_BL21.m8
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

if cls_BL is not None:
    print("이미 분류가 로드되어 있어 건너뜁니다.")
else:
    bl_ids = list(PROTEOMES["BL21DE3"]["seqs"])
    mg_ids = list(PROTEOMES["MG1655"]["seqs"])
    cls_BL = classify(DIR["search"]/"BL21_vs_MG.m8", bl_ids)
    cls_MG = classify(DIR["search"]/"MG_vs_BL21.m8", mg_ids)
    print("=== BL21 기준 ===");   print(cls_BL.category.value_counts().to_string())
    print("\n=== MG1655 기준 (BL21 결손 = 억제자 후보) ===")
    print(cls_MG.category.value_counts().to_string())
    cls_MG.to_csv(DIR["results"]/"classify_MG1655_vs_BL21.csv", index=False)
    print("\n기대값(BL21): 2171 / 1575 / 125 / 217, (MG1655 역방향): 2199 / 1569 / 152 / 380")
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
assert cls_BL is not None, "분류 결과가 없습니다. CELL 09 또는 CELL 10-11 을 먼저 실행하세요."
cls_BL.to_csv(DIR["results"]/"classify_BL21_vs_MG1655.csv", index=False)

known = set(hdr2seq)
cls_BL["in_proteome"] = cls_BL.protein.isin(known)
if (~cls_BL.in_proteome).any():
    n = int((~cls_BL.in_proteome).sum())
    print(f"⚠ 분류에는 있으나 프로테옴 FASTA 에 없는 ID {n}개 — ID 체계(GenBank 버전 .1 등) 확인")
    print("  예시:", cls_BL.loc[~cls_BL.in_proteome, "protein"].head(5).tolist())

TRACK_A = cls_BL.query("category in ['identical','high-similarity'] and in_proteome").protein.tolist()
TRACK_B = cls_BL.query("category in ['low-similarity','strain-specific'] and in_proteome").protein.tolist()

pd.Series(TRACK_A).to_csv(DIR["results"]/"track_A_ids.txt", index=False, header=False)
pd.Series(TRACK_B).to_csv(DIR["results"]/"track_B_ids.txt", index=False, header=False)
print(f"Track A (RF2-PPI): {len(TRACK_A):5d}   (기대 ~3,746)")
print(f"Track B (Boltz-2): {len(TRACK_B):5d}   (기대 ~342)")
print("저장:", DIR["results"]/"track_A_ids.txt", "/", DIR["results"]/"track_B_ids.txt")
```

---

## CELL 13 — Part 2. Bacteria reference proteome 다운로드 (백그라운드)

```python
# =============================================================================
# CELL 13 | Part 2-1. UniProt Reference Proteomes (Bacteria) 다운로드
#   - paired MSA 는 "같은 genome 안의 bait/prey orthologue 짝짓기"라 이 DB 가 필요하다
#   - 약 12GB 다운로드 + 압축해제 ~40GB. ⚠ /mnt/af2results 여유 ~1.1T 확인하고 시작할 것
#   - 백그라운드 30~60분
# =============================================================================
script = f"""
cd "{BASE}/db"
REL=https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes
FN=$(curl -s "$REL/" | grep -o 'Reference_Proteomes_[0-9_]*\\.tar\\.gz' | sort -u | tail -1)
echo "release file = $FN"
[ -n "$FN" ] || {{ echo "릴리스 파일명을 못 찾음"; exit 1; }}
if command -v aria2c >/dev/null; then
  aria2c -x8 -s8 -c "$REL/$FN"
else
  curl -C - -O "$REL/$FN"
fi
tar -xzf "$FN" Bacteria/
du -sh Bacteria
echo DONE_download
"""
sh_bg("part2_download", script)
print("\n디스크 여유 먼저 확인:")
sh(f"df -h {BASE}", check=False)
```

---

## CELL 14 — Part 2b. 병합 + mmseqs createdb (백그라운드)

```python
# =============================================================================
# CELL 14 | Part 2-2. taxid 를 헤더에 박아 하나로 병합 -> mmseqs DB 생성
#   - 헤더를 ">taxid|accession" 으로 만들어 두면 pairing 시 organism key 가 공짜
#   - CELL 13 의 DONE_download 확인 후 실행. 1~2시간.
# =============================================================================
script = f"""
cd "{BASE}/db"
OUT=bacteria_ref.fasta
: > "$OUT"
find Bacteria -name '*.fasta.gz' ! -name '*_DNA*' ! -name '*additional*' | while read f; do
  tax=$(basename "$f" .fasta.gz | cut -d_ -f2)
  zcat "$f" | awk -v t="$tax" '/^>/{{split($1,a,"|"); acc=(length(a)>=2?a[2]:substr($1,2)); print ">"t"|"acc; next}}{{print}}'
done >> "$OUT"
echo "sequences: $(grep -c '^>' "$OUT")"
echo "proteomes: $(grep '^>' "$OUT" | cut -d'|' -f1 | tr -d '>' | sort -u | wc -l)"
mmseqs createdb "$OUT" bactDB
mmseqs createindex bactDB "{BASE}/tmp/idx" --threads {THREADS} --search-type 1 || true
ls -la
echo DONE_bactdb
"""
sh_bg("part2_bactdb", script, env=CONDA_ENV_RF2)
```

---

## CELL 15 — 백그라운드 작업 모니터 (수시 실행)

```python
# =============================================================================
# CELL 15 | 백그라운드 작업 상태 확인 — 긴 단계 돌릴 때마다 이 셀로 확인
# =============================================================================
for job in ["easy_search", "part2_download", "part2_bactdb",
            "part3_gate", "part4_prey", "part6_rf2ppi", "part7_boltz"]:
    if (DIR["logs"] / f"{job}.log").exists():
        bg_tail(job, n=8)
        print("-" * 70)

print("\n### 디스크 ###"); sh(f"df -h {BASE}", check=False)
print("### GPU ###");   sh("nvidia-smi --query-gpu=index,utilization.gpu,memory.used,memory.total "
                           "--format=csv,noheader", check=False)
```

---

## CELL 16 — Part 3. Feasibility gate 검색

```python
# =============================================================================
# CELL 16 | Part 3-1. Feasibility gate — bait orthologue 보유 genome 수 세기
#   - paired MSA 깊이 상한 = min(bait orthologue genome 수, prey orthologue genome 수)
#   - CooS 는 CO-oxidizing anaerobe 에만 분포 -> 이 값이 낮을 수 있다. 절대 건너뛰지 말 것.
# =============================================================================
fmt = "query,target,fident,evalue,bits,qstart,qend,qlen,tstart,tend,tlen,qaln,taln"
script = f"""
cd "{BASE}"
mmseqs createdb seq/baits.fasta db/baitDB
mmseqs search db/baitDB db/bactDB search/bait_res tmp/bs \\
  -s 7.5 --num-iterations 3 -e 1e-3 --max-seqs 20000 --threads {THREADS}
mmseqs convertalis db/baitDB db/bactDB search/bait_res search/bait_hits.m8 \\
  --format-output "{fmt}" --threads {THREADS}
wc -l search/bait_hits.m8
echo DONE_gate
"""
sh_bg("part3_gate", script, env=CONDA_ENV_RF2)
```

---

## CELL 17 — Part 3b. 깊이 판정 (★체크포인트 2)

```python
# =============================================================================
# CELL 17 | Part 3-2. Gate 판정 — Track A 를 진행할지 여기서 결정한다
#   >=500 : Track A 진행 / 200-500 : 진행하되 replicate 5회, 경계 점수 불신
#   <200  : Track A 포기, Track B(CELL 25~)에 자원 집중
# =============================================================================
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
MAX_PREY_LEN = 1000
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
    DIR["results"]/"trackA_skipped_long.csv", index=False)
```

---

## CELL 19 — Part 4b. prey homolog 검색 (최장 단계)

```python
# =============================================================================
# CELL 19 | Part 4-2. prey homolog 검색 — 이 파이프라인에서 가장 오래 걸린다 (4~10h)
#   - 반드시 백그라운드. 이 시간 동안 다른 GPU 에서 Track B(CELL 25~)를 먼저 돌려도 된다.
# =============================================================================
script = f"""
cd "{BASE}"
mmseqs createdb seq/prey_trackA.fasta db/preyDB
mmseqs search db/preyDB db/bactDB search/prey_res tmp/ps \\
  -s 7.5 --num-iterations 3 -e 1e-3 --max-seqs 20000 --threads {THREADS}
mmseqs convertalis db/preyDB db/bactDB search/prey_res search/prey_hits.m8 \\
  --format-output "{fmt}" --threads {THREADS}
wc -l search/prey_hits.m8
echo DONE_prey
"""
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
from collections import defaultdict

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

def best_hit_per_taxid(m8_path, chunksize=2_000_000):
    out = defaultdict(dict)
    for chunk in pd.read_csv(m8_path, sep="\t", names=AC, chunksize=chunksize):
        chunk["taxid"] = chunk.target.str.split("|").str[0]
        for r in chunk.itertuples(index=False):
            cur = out[r.query].get(r.taxid)
            if cur is None or r.bits > cur[0]:
                out[r.query][r.taxid] = (r.bits,
                                         to_query_frame(r.qaln, r.taln, r.qstart, r.qlen))
    return out

print("bait hits 정리 중...")
BAIT_MSA = best_hit_per_taxid(DIR["search"]/"bait_hits.m8")
print("prey hits 정리 중... (수 분 소요)")
PREY_MSA = best_hit_per_taxid(DIR["search"]/"prey_hits.m8")
print("bait depth:", {k: len(v) for k, v in BAIT_MSA.items()})
print("prey 개수:", len(PREY_MSA))
```

---

## CELL 21 — Part 5b. paired MSA 생성 + hhfilter

```python
# =============================================================================
# CELL 21 | Part 5-2. paired MSA 생성 (같은 organism 끼리 이어 붙이고 hhfilter 90%)
#   - MIN_PAIRED 미만이면 예측 자체를 스킵한다 (얕은 MSA 의 점수는 신뢰 불가)
#   - USE_SEGMENT_BAIT 이면 full-length + segment bait 를 모두 돌린다 (설계문서 대응 b/c)
# =============================================================================
MIN_PAIRED = 50
paired_dir = DIR["paired"]; paired_dir.mkdir(exist_ok=True)

BAITS_TO_RUN = [BAIT_KEY] + ([k for k in SEGMENTS] if USE_SEGMENT_BAIT else [])
BAITS_TO_RUN = [b for b in BAITS_TO_RUN if b in BAIT_MSA]
print("돌릴 bait:", len(BAITS_TO_RUN), "개")

input_lines, stats = [], []
for bkey in BAITS_TO_RUN:
    bseq, brows = BAIT_ALL[bkey], BAIT_MSA[bkey]
    L1 = len(bseq)
    for pid, prows in PREY_MSA.items():
        shared = set(brows) & set(prows)
        if len(shared) < MIN_PAIRED:
            stats.append((bkey, pid, len(shared), "skip")); continue
        raw = paired_dir / f"{bkey}__{pid.replace('|','_')}.raw.a3m"
        with open(raw, "w") as fh:
            fh.write(f">query\n{bseq}{sel[pid]}\n")
            for tx in sorted(shared):
                fh.write(f">{tx}\n{brows[tx][1]}{prows[tx][1]}\n")
        flt = str(raw).replace(".raw.a3m", ".a3m")
        subprocess.run(["hhfilter", "-i", str(raw), "-o", flt, "-id", "90", "-M", "first"],
                       check=True, capture_output=True)
        depth = sum(1 for l in open(flt) if l.startswith(">"))
        input_lines.append(f"{flt} {L1}")
        stats.append((bkey, pid, depth, "ok"))
        raw.unlink()

(DIR["rf2ppi"]/"input_file").write_text("\n".join(input_lines) + "\n")
sdf = pd.DataFrame(stats, columns=["bait","prey","paired_depth","status"])
sdf.to_csv(DIR["results"]/"paired_msa_stats.csv", index=False)
print(f"\n생성 {len(input_lines)}개 / 스킵 {(sdf.status=='skip').sum()}개")
print(sdf.query("status=='ok'").paired_depth.describe().to_string())
```

---

## CELL 22 — Part 6. RF2-PPI ×N replicate 실행

```python
# =============================================================================
# CELL 22 | Part 6-1. RF2-PPI 추론 — replicate 필수
#   - 비결정적 모델이다. README 기준 ~5% 쌍에서 SD>0.1, 중간 점수대 변동이 크다.
#   - N_REPLICATE 는 CELL 17 gate 결과에 따라 3 또는 5 (깊이 200-500 이면 5회)
# =============================================================================
reps = " ".join(str(i) for i in range(1, N_REPLICATE + 1))
script = f"""
cd "{BASE}/rf2ppi"
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
reps_df = []
for rep in range(1, N_REPLICATE + 1):
    log = DIR["rf2ppi"]/f"input_rep{rep}.log"
    if not log.exists():
        print(f"[없음] {log}"); continue
    d = pd.read_csv(log, sep=r"\s+", names=["msa","prob","sec"])
    stem = d.msa.apply(lambda p: Path(p).stem)
    d["bait"] = stem.apply(lambda s: s.split("__")[0])
    d["prey"] = stem.apply(lambda s: s.split("__")[-1])
    reps_df.append(d.set_index(["bait","prey"])["prob"].rename(f"rep{rep}"))

R = pd.concat(reps_df, axis=1)
R["mean"], R["sd"] = R.mean(axis=1), R.std(axis=1)
R = R.sort_values("mean", ascending=False)

# full-length / segment 를 prey 단위로 합집합(최대값)까지 같이 본다
R_best = R.reset_index().groupby("prey").agg(
    best_bait=("bait", lambda s: s.iloc[0]), mean=("mean","max"), sd=("sd","first")
).sort_values("mean", ascending=False)

meta = cls_BL.set_index("protein")
R_best = R_best.join(meta[["category","pident"]], how="left")
R_best = R_best.join(sdf.query("status=='ok'").groupby("prey").paired_depth.max(), how="left")
R_best["desc"] = [hdr2desc.get(i, "")[:70] for i in R_best.index]

R.to_csv(DIR["results"]/"trackA_RF2PPI_all_pairs.csv")
R_best.to_csv(DIR["results"]/"trackA_RF2PPI_ranked.csv")
print("=== 상위 40 후보 (bait 합집합) ===")
print(R_best.head(40).to_string())
print("\nmean>=0.74 (strict):", int((R_best["mean"] >= 0.74).sum()))
print("mean>=0.05 (loose) :", int((R_best["mean"] >= 0.05).sum()))
print("sd>0.1 (재현성 낮음, 재실행 권장):", int((R_best["sd"] > 0.1).sum()))
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

ctrl = [k for k in BAIT_SEQ if k.upper().startswith(("COOC","COOT","COOJ"))]
if ctrl:
    print("양성대조군 bait:", ctrl)
    print("  ※ CooC-CooS 쌍은 CELL 21 에서 bait-bait paired MSA 를 따로 만들어 돌려야 한다.")
    print("     (BAITS_TO_RUN 에 대조군 bait 를 넣고 prey 자리에 ChCODH2 를 두는 방식)")
else:
    print("⚠ 양성대조군 bait 가 없다 — 이 상태의 결과는 검증되지 않은 값이다.")

import matplotlib.pyplot as plt
d = R_best.dropna(subset=["paired_depth"])
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].scatter(d.paired_depth, d["mean"], s=6, alpha=.4)
ax[0].set_xlabel("paired MSA depth"); ax[0].set_ylabel("RF2-PPI prob (mean)")
ax[0].set_title("depth vs score (상관 있으면 artifact 의심)")
ax[1].hist(R_best["mean"], bins=60); ax[1].set_yscale("log")
ax[1].set_xlabel("RF2-PPI prob"); ax[1].set_title("score distribution")
plt.tight_layout(); plt.show()
print("Spearman r =", round(d[["paired_depth","mean"]].corr(method="spearman").iloc[0,1], 3))
```

---

## CELL 25 — Part 7. Boltz-2 설치

```python
# =============================================================================
# CELL 25 | Part 7-1. Boltz-2 설치 스크립트 (Track B — 여기가 주력일 가능성이 높다)
#   - 공진화 불필요. strain-specific 217 + low-sim 125 = 342개가 대상.
#   - Ni 을 ligand(CCD: NI)로 명시할 수 있는 게 이 문제에서 결정적 장점.
# =============================================================================
boltz_sh = f"""
source "$(conda info --base)/etc/profile.d/conda.sh"
conda create -y -n {CONDA_ENV_BOLTZ} python=3.11
conda activate {CONDA_ENV_BOLTZ}
pip install boltz -U
pip install pyyaml pandas
python -c "import boltz; print('boltz ok')"
"""
p = DIR["scripts"] / "install_boltz.sh"
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

n, long_skip = 0, []
for pid in TRACK_B:
    s = hdr2seq[pid]
    if len(s) > MAX_B_LEN:
        long_skip.append((pid, len(s))); continue
    write_boltz_yaml(f"{BAIT_KEY}__{pid.replace('|','_')}",
                     BAIT_SEQ[BAIT_KEY], s, with_ni=True)
    n += 1
print(f"Boltz-2 입력 {n}개 생성 -> {boltz_in}")
print(f"길이>{MAX_B_LEN} 제외: {len(long_skip)}개  {long_skip[:5]}")
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
BOLTZ_GPU = GPU_ID
script = f"""
cd "{BASE}/boltz"
CUDA_VISIBLE_DEVICES={BOLTZ_GPU} boltz predict inputs \\
  --out_dir out \\
  --use_msa_server \\
  --recycling_steps 3 \\
  --diffusion_samples 1 \\
  --output_format mmcif \\
  --num_workers 4
echo DONE_boltz
"""
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
    B["category"] = B.prey.map(cls_BL.set_index("protein")["category"])
    B["desc"] = B.prey.map(lambda i: hdr2desc.get(i, "")[:70])
B.to_csv(DIR["results"]/"trackB_boltz2_ranked.csv", index=False)
print(f"파싱 {len(B)}건")
print(B.head(40).to_string(index=False))
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
RESIDUES_METAL = "A112,A114"     # <- CELL 06 에서 확인한 체인으로 수정
RESIDUES_ATP   = ""              # <- Walker A 등 지정 후 run02 실행 (비우면 스킵)

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

## CELL 30 — Part 8b. ID crosswalk (구조 tid → GenBank protein ID)

```python
# =============================================================================
# CELL 30 | Part 8-2. ★ID crosswalk — Folddisco tid 를 GenBank ID 에 붙인다
#   문제: 구조 파일명이 AF-<UniProtID>-F1-model_v6.cif / cf_<UniParcID>.pdb 인데
#         우리 분류·프로테옴은 GenBank ID (QZI…, AAC…) 라 그냥은 조인이 안 된다.
#   해법: 구조 파일에서 CA 기준 서열을 뽑아 프로테옴 서열과 "완전일치" 매칭 (오프라인·정확)
#   - 5,000여 개 파싱에 수 분. 결과는 CSV 로 캐시하므로 한 번만 돌리면 된다.
# =============================================================================
XWALK_CSV = DIR["results"]/"id_crosswalk_struct_to_genbank.csv"

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

if XWALK_CSV.exists():
    XW = pd.read_csv(XWALK_CSV)
    print("캐시 로드:", XWALK_CSV, len(XW), "행")
else:
    seq2gb = {}
    for strain in ["BL21DE3", "Y19", "MG1655"]:
        for pid, s in PROTEOMES.get(strain, {}).get("seqs", {}).items():
            seq2gb.setdefault(s.rstrip("*"), (pid, strain))
    rows = []
    for tag, d in [("BL21", ASSET["struct_bl21"]), ("Y19", ASSET["struct_y19"])]:
        if not d.exists():
            print("[없음]", d); continue
        files = sorted(list(d.glob("*.cif")) + list(d.glob("*.pdb")))
        print(f"{tag}: {len(files)} 구조 파싱 중...")
        for i, f in enumerate(files):
            if i and i % 1000 == 0: print(f"  {i}/{len(files)}")
            tid = f.stem
            m = re.match(r"AF-([A-Z0-9]+)-F\d+", tid)
            acc_type, acc = ("UniProt", m.group(1)) if m else (
                ("UniParc", tid[3:]) if tid.startswith("cf_") else ("unknown", tid))
            s = struct_seq(f)
            gb, strain = seq2gb.get(s, (None, None))
            rows.append({"tid": tid, "db": tag, "acc_type": acc_type, "acc": acc,
                         "struct_len": len(s), "protein": gb, "matched_strain": strain})
    XW = pd.DataFrame(rows)
    XW.to_csv(XWALK_CSV, index=False)
    print("저장:", XWALK_CSV)

if len(XW):
    hit = XW.protein.notna().mean()
    print(f"\n매칭률 {hit:.1%}  ({XW.protein.notna().sum()}/{len(XW)})")
    print(XW.groupby(["db","acc_type"]).agg(n=("tid","size"),
                                            matched=("protein", lambda s: s.notna().sum())).to_string())
    print("\n미매칭이 많으면: 구조가 isoform/부분서열이거나 프로테옴 버전이 다른 것.")
    print("  -> UniProt idmapping API (UniProtKB_AC-ID -> EMBL-GenBank-DDBJ_CDS) 로 보완 가능")
```

---

## CELL 31 — Part 8c. Folddisco 결과 → motif CSV 변환

```python
# =============================================================================
# CELL 31 | Part 8-3. Folddisco TSV -> Part 8 이 기대하는 (protein, score) CSV 로 변환
#   - score 는 idf 를 기본으로 쓰고, min_rmsd 는 낮을수록 좋으므로 역수로 보정 결합
# =============================================================================
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
        frames.append(d[["protein","score","strain","tid","tid_stem"] +
                        [c for c in ["idf","min_rmsd","plddt","matching_residues"] if c in d.columns]])
    if not frames:
        print("변환할 결과 없음"); return pd.DataFrame(columns=["protein","score"])
    A = pd.concat(frames)
    out = (A.dropna(subset=["protein"]).groupby("protein", as_index=False)
             .score.max().sort_values("score", ascending=False))
    out.to_csv(out_csv, index=False)
    A.to_csv(str(out_csv).replace(".csv", "_detail.csv"), index=False)
    print(f"{out_csv.name}: {len(out)}개 단백질 (매칭 실패 {A.protein.isna().sum()}행 제외)")
    return out

FOLDDISCO_NI  = DIR["results"]/"folddisco_metal_motif.csv"
FOLDDISCO_ATP = DIR["results"]/"folddisco_atp_motif.csv"
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
# =============================================================================
def safe_series(df, key, col):
    try:
        return df.set_index(key)[col]
    except Exception:
        return pd.Series(dtype=float)

M = pd.DataFrame(index=sorted(set(TRACK_A) | set(TRACK_B)))
M["category"]      = cls_BL.set_index("protein")["category"]
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
M["n_evidence"] = ((M.rf2ppi > 0.3).fillna(False).astype(int)
                   + (M.boltz_iptm > 0.6).fillna(False).astype(int)
                   + M.ni_motif.notna().astype(int) + M.atp_motif.notna().astype(int))
M = M.sort_values(["priority", "n_evidence"], ascending=False)
M.to_csv(DIR["results"]/"INTEGRATED_candidate_ranking.csv")
print(M.head(50).to_string())
print("\n독립적 2개 이상 증거를 가진 후보:", int((M.n_evidence >= 2).sum()))
```

---

## CELL 33 — 최종 리포트 + 남은 TODO

```python
# =============================================================================
# CELL 33 | 최종 요약 + 남은 TODO
# =============================================================================
print("### 산출물 ###")
for f in sorted(DIR["results"].glob("*.csv")):
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
  [ ] BL21 proteome ID: 이 노트북은 UP000503272 로 통일했다. 협업팀 확인 남음
  [ ] Folddisco ATP-binding motif residue 지정 (CELL 29 RESIDUES_ATP)
  [ ] 3kji.pdb 체인 ID 확정 (CELL 06 출력으로 확인)
  [ ] crosswalk 미매칭분 UniProt idmapping 으로 보완 (CELL 30)
  [ ] apt-mark hold nvidia-* (드라이버 재발 방지)
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
[ ] CELL 32-33  통합 랭킹 + 리포트
```

**병렬화**: CELL 13-21(CPU)이 도는 동안 CELL 25-28(Boltz-2 Track B)을 다른 GPU에서 먼저 시작할 수 있다.
Track B는 MSA 서버를 쓰므로 Part 2가 필요 없다.
