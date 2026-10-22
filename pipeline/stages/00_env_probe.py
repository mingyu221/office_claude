#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""단계 0 — 재현성 고정: 환경 탐지 · 219 클러스터 재현 · 잃어버린 search 명령 복원.

인계문서(20261021_handoff_mmseqs_env.md) 5절 확인명령 + 7절 219 재현 +
8절 빈칸 채우기 + 2절 search 명령 복원을 하나의 스크립트로 자동화한다.
epel-af2 에서 **가장 먼저** 돌릴 스크립트다.

무엇을 읽는가
  (탐지)  PATH / CONDA_PREFIX / conda env 목록 / jupyter kernelspec / nproc / df,
          mmseqs 바이너리(--mmseqs 또는 common.find_mmseqs), readlink -f, ldd, file
  (재현)  --ipr-fasta   기본 $MNT/chcodh2/protein-matching-IPR010047.fasta
          --master-tsv  기본 $SPLIT/master_1827.tsv           (멤버십 대조용)
  (복원)  --scan-dir    기본 $DOCK 아래 모든 *.ipynb (셀 소스 + 셀 출력) 및
                        *.sh / *.log / *.txt / *.py / *.out / *.err
          --history     기본 ~/.bash_history ~/.zsh_history ~/.sh_history
                        ~/.histfile ~/.python_history ~/.local/share/fish/fish_history
                        ~/.ipython/profile_default/history.sqlite (표준 sqlite3 로 읽음)
          --align-tsv   기본 $OLDR/chcodh2_alignment.tsv — 헤더만 읽어 후보 채점에 쓴다

무엇을 쓰는가  (쓰기는 outdir 안에서만. 원본은 건드리지 않는다)
  <outdir>/env_probe.json   기계용 — 탐지·재현·복원 전체 결과
  <outdir>/env_probe.md     사람용 — 인계문서 8절 빈칸을 **채운 블록** 포함
  <outdir>/params.json      common.write_params (행 수·조인율·결측률·상한 포함)
  <outdir>/logs/*.log       실행한 외부 명령의 stdout/stderr
  <repro-dir>/id70_*        --reproduce 일 때만. mmseqs easy-cluster 산출물

무엇을 보장하는가
  1. 탐지 실패는 예외가 아니라 "없음" / "미확인(사유)" 로 기록된다. mmseqs·conda·
     jupyter·데이터가 하나도 없는 장비에서도 **끝까지 돌고** 두 산출물을 남긴다.
     그것이 이 스크립트의 자기검증이다.
  2. 숫자는 전부 세어서 stdout 과 params.json 양쪽에 적는다. 상한(파일 수·크기·
     후보 수)에 걸려 잘린 것은 "몇 개 처리 / 몇 개 버림" 을 명시한다.
  3. 219 재현에서 **"대표 ID 동일"** 과 **"분할 1:1 동일"** 을 구분해 보고한다.
     대표 ID 가 달라도 분할이 1:1 이면 실질적으로 같은 분할이다(인계문서 7절).
  4. 서열을 지우지 않는다. 이 단계는 서열 표를 쓰지 않으며 master_1827.tsv 는
     읽기 전용으로만 쓴다.
  5. 불일치는 discrepancies 로 모아 stdout·env_probe.json·env_probe.md·params.json
     네 곳에 적는다. --strict 면 산출물을 모두 쓴 뒤 종료코드 2 로 끝낸다.
     (조용히 넘어가지 않기 위한 장치이며, 산출물을 잃지 않기 위해 마지막에 끝낸다.)
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import traceback
import unicodedata
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# common 임포트 — 레포 구조 pipeline/stages/00_env_probe.py 기준 ../lib
# ---------------------------------------------------------------------------
_HERE = Path(__file__).resolve()
_LIB = _HERE.parents[1] / "lib"
if str(_LIB) not in sys.path:
    sys.path.insert(0, str(_LIB))
import common  # noqa: E402

import pandas as pd  # noqa: E402  (common 이 이미 요구한다)

NOT_FOUND = "없음"
UNKNOWN = "미확인"

# 인계문서 7절 — 이 숫자가 재현되어야 환경이 당시와 같다
EXPECTED_N_CLUSTERS = 219
EXPECTED_N_SEQUENCES = 1827

# 인계문서 1절 — selection_params.txt 에 기록된 클러스터링 옵션. 바꾸면 219 가 깨진다
CLUSTER_OPTS = ["--min-seq-id", "0.70", "-c", "0.8", "--cov-mode", "1", "-s", "7.5"]

# 인계문서 9절 — 복원 단서가 있을 가능성이 높은 노트북
KEY_NOTEBOOKS = ["ipr010047_mmseqs.ipynb", "homolog_split.ipynb", "cluster_sweep.ipynb"]

# 복원 후보 판정 — 'mmseqs' 와 아래 서브명령/옵션이 같은 줄에 나오면 후보
RE_MMSEQS = re.compile(r"mmseqs", re.IGNORECASE)
RE_SUBCMD = re.compile(r"\b(easy-search|search|convertalis|format-output)\b", re.IGNORECASE)
# 'research' 같은 단어는 \b 때문에 걸리지 않는다. 'easy-search' 는 걸린다.
RE_FORMAT_OUT = re.compile(
    r"--format-output[=\s]+(?:[\"']([^\"']+)[\"']|([^\s\"';|&)\]]+))", re.IGNORECASE
)
RE_PATHISH = re.compile(r"[\w./$~{}-]+\.(?:fasta|fa|faa|pdb|fas|seq)\b", re.IGNORECASE)
RE_LOOKS_CMD = re.compile(
    r"mmseqs[\"'\s,\]\[]+(?:easy-search|search|convertalis|easy-cluster|createdb)", re.IGNORECASE
)
RE_WS = re.compile(r"\s+")
# mmseqs 전용 토큰 — 이것이 보이면 'mmseqs' 문자열이 없어도 mmseqs 호출로 본다
RE_SPECIFIC = re.compile(r"\b(easy-search|convertalis|format-output)\b", re.IGNORECASE)
RE_ALIAS_ASSIGN = re.compile(r"(?:export\s+)?([A-Za-z_]\w{0,30})\s*=")
ALIAS_STOPWORDS = {"if", "for", "while", "print", "return", "echo", "python", "pip"}

# 상한 (전부 보고한다 — 7번 규칙)
CAP = dict(
    scan_file_mb=64,          # 파일 하나 최대 크기
    scan_notebooks=2000,      # 노트북 개수
    scan_text_files=5000,     # 텍스트 파일 개수
    scan_depth=12,            # 디렉터리 깊이
    history_lines=400_000,    # 히스토리 한 파일당 줄 수
    ipython_rows=200_000,     # IPython sqlite 행 수
    candidates_json=200,      # env_probe.json 에 남길 후보 수
    candidates_md=10,         # env_probe.md 에 남길 후보 수
    tmp_dirs=60,              # tmp 잔여물 개수
    tmp_du_entries=200_000,   # tmp 용량 계산 시 훑을 엔트리 수
    conda_list_envs=12,       # `conda list -n <env>` 를 돌릴 환경 수
    cmd_output_chars=8000,    # 외부 명령 출력 저장 길이
    cand_text_chars=2000,     # 후보 텍스트 저장 길이
)

SKIP_DIR_NAMES = {".git", "__pycache__", "node_modules", "site-packages", ".cache", ".mypy_cache"}
TEXT_EXTS = {".sh", ".log", ".txt", ".py", ".out", ".err"}


# ===========================================================================
# 작은 도구들 — 어떤 경우에도 예외를 올리지 않는다
# ===========================================================================
def now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def probe_cmd(cmd, *, timeout=120, dry_run=False, log_dir=None, label=None) -> dict:
    """외부 명령 실행. 실패를 예외로 올리지 않고 status 로 기록한다.

    status:  ok | 비정상종료 | 없음 | 시간초과 | 오류 | dry-run
    """
    cmd = [str(c) for c in cmd]
    rec = dict(cmd=cmd, cmd_str=" ".join(cmd), status=UNKNOWN, stdout="", stderr="", returncode=None)
    if dry_run:
        rec["status"] = "dry-run"
        return rec
    try:
        r = subprocess.run(
            cmd, capture_output=True, text=True, errors="replace", timeout=timeout
        )
    except FileNotFoundError:
        rec["status"] = NOT_FOUND
        rec["stderr"] = f"실행파일을 찾을 수 없음: {cmd[0]}"
    except subprocess.TimeoutExpired:
        rec["status"] = "시간초과"
        rec["stderr"] = f"{timeout}s 초과"
    except Exception as e:  # noqa: BLE001  탐지 단계에서는 무엇도 치명적이지 않다
        rec["status"] = "오류"
        rec["stderr"] = f"{type(e).__name__}: {e}"
    else:
        rec["returncode"] = r.returncode
        rec["status"] = "ok" if r.returncode == 0 else f"비정상종료(rc={r.returncode})"
        out, err = r.stdout or "", r.stderr or ""
        rec["stdout"] = out[: CAP["cmd_output_chars"]]
        rec["stderr"] = err[: CAP["cmd_output_chars"]]
        rec["stdout_truncated"] = len(out) > CAP["cmd_output_chars"]
        rec["stderr_truncated"] = len(err) > CAP["cmd_output_chars"]
        if log_dir:
            try:
                log_dir = Path(log_dir)
                log_dir.mkdir(parents=True, exist_ok=True)
                name = re.sub(r"[^\w.-]+", "_", label or Path(cmd[0]).name)[:80]
                (log_dir / f"{name}.log").write_text(
                    f"$ {rec['cmd_str']}\nrc={r.returncode}\n\n"
                    f"--- stdout ---\n{out}\n--- stderr ---\n{err}\n"
                )
            except OSError:
                pass
    return rec


def first_line(text: str, default: str = "") -> str:
    for ln in (text or "").splitlines():
        if ln.strip():
            return ln.strip()
    return default


def walk_limited(root: Path, maxdepth: int):
    """깊이 제한 + 잡음 디렉터리 제외 walk. 접근 실패는 조용히 건너뛴다."""
    root = Path(root)
    base = len(root.parts)
    for dirpath, dirnames, filenames in os.walk(root, topdown=True, onerror=lambda e: None):
        depth = len(Path(dirpath).parts) - base
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIR_NAMES]
        if depth >= maxdepth:
            dirnames[:] = []
        yield Path(dirpath), dirnames, filenames


def read_text_capped(path: Path, max_bytes: int):
    """(텍스트, 사유). 너무 크거나 못 읽으면 (None, 사유)."""
    try:
        size = path.stat().st_size
    except OSError as e:
        return None, f"stat 실패: {e}"
    if size > max_bytes:
        return None, f"크기 초과 ({size/1e6:.1f} MB > {max_bytes/1e6:.0f} MB)"
    try:
        return path.read_text(errors="replace"), None
    except OSError as e:
        return None, f"읽기 실패: {e}"


# ===========================================================================
# (a) 환경 탐지
# ===========================================================================
def detect_python_env(dry_run: bool) -> dict:
    """인계문서 5.4 의 파이썬 블록. 이 스크립트를 돌린 인터프리터가 기준이다."""
    path_env = os.environ.get("PATH", "")
    which_mmseqs = shutil.which("mmseqs")
    out = dict(
        python_executable=sys.executable or UNKNOWN,
        python_version=sys.version.split()[0],
        python_full=sys.version.replace("\n", " "),
        mmseqs_on_path=which_mmseqs or NOT_FOUND,
        PATH=path_env or NOT_FOUND,
        PATH_entries=len([p for p in path_env.split(os.pathsep) if p]),
        CONDA_PREFIX=os.environ.get("CONDA_PREFIX") or NOT_FOUND,
        CONDA_DEFAULT_ENV=os.environ.get("CONDA_DEFAULT_ENV") or NOT_FOUND,
        VIRTUAL_ENV=os.environ.get("VIRTUAL_ENV") or NOT_FOUND,
        LD_LIBRARY_PATH=os.environ.get("LD_LIBRARY_PATH") or NOT_FOUND,
        MPLBACKEND=os.environ.get("MPLBACKEND") or NOT_FOUND,
        user=os.environ.get("USER") or os.environ.get("LOGNAME") or UNKNOWN,
        host=os.uname().nodename,
        pandas=pd.__version__,
        numpy=__import__("numpy").__version__,
        inside_jupyter=bool(os.environ.get("JPY_PARENT_PID")),
    )
    out["uname"] = probe_cmd(["uname", "-a"], dry_run=dry_run, timeout=20)["stdout"].strip() or UNKNOWN
    return out


def infer_install_method(binary: Path, ldd_rec: dict, file_rec: dict) -> dict:
    """설치 방식 추정 — 경로 표식 + ldd + 배포 tarball 레이아웃 + 빌드 트리."""
    evidence, guesses = [], []
    low = str(binary).lower()
    if any(k in low for k in ("miniconda", "anaconda", "/conda/", "conda-bld", "mambaforge", "miniforge")):
        guesses.append("conda")
        evidence.append(f"경로에 conda 표식: {binary}")
    if "/envs/" in low:
        guesses.append("conda")
        evidence.append("경로에 '/envs/' — conda 환경 안")
    if ".pixi/envs" in low:
        guesses.append("기타(pixi)")
        evidence.append("경로에 '.pixi/envs' — pixi 환경")

    ldd_text = ((ldd_rec.get("stdout") or "") + (ldd_rec.get("stderr") or "")).lower()
    if "not a dynamic executable" in ldd_text or "statically linked" in ldd_text:
        guesses.append("정적바이너리")
        evidence.append("ldd: not a dynamic executable / statically linked")
    elif ldd_rec.get("status") == "ok" and ldd_text.strip():
        evidence.append(f"ldd: 동적 링크 ({len([l for l in ldd_text.splitlines() if '=>' in l])} 라이브러리)")
    elif ldd_rec.get("status") == NOT_FOUND:
        evidence.append("ldd 명령 없음 — 링크 방식 미확인")

    file_text = (file_rec.get("stdout") or "").lower()
    if "static" in file_text:
        guesses.append("정적바이너리")
        evidence.append(f"file: {first_line(file_rec.get('stdout',''))}")
    elif file_text.strip():
        evidence.append(f"file: {first_line(file_rec.get('stdout',''))}")

    # 공식 배포 tarball 레이아웃: <root>/bin/mmseqs 와 형제로 util/ matrices/ examples/
    try:
        root = binary.parent.parent
        sib = {p.name for p in root.iterdir()} if root.is_dir() else set()
        if binary.parent.name == "bin" and {"util", "matrices"} & sib:
            guesses.append("정적바이너리(공식 배포 tarball)")
            evidence.append(f"배포 tarball 레이아웃: {root}/{{bin,util,matrices}}")
    except OSError:
        pass

    # 소스 빌드 흔적
    try:
        for up in [binary.parent, *binary.parents[:4]]:
            if (up / "CMakeCache.txt").is_file():
                guesses.append("소스빌드")
                evidence.append(f"CMakeCache.txt 발견: {up/'CMakeCache.txt'}")
                break
    except (OSError, IndexError):
        pass

    if not guesses:
        guesses.append("기타/판정불가")
        evidence.append("경로·ldd·레이아웃 어디서도 표식을 찾지 못함")
    # 중복 제거(순서 유지). 더 구체적인 라벨이 있으면 포괄 라벨은 버린다
    uniq = list(dict.fromkeys(guesses))
    uniq = [g for g in uniq if not any(g != o and g in o for o in uniq)]
    return dict(guess=" 또는 ".join(uniq), candidates=uniq, evidence=evidence)


def detect_mmseqs(args, log_dir: Path) -> dict:
    """인계문서 5.1 / 5.2 / 6절 — 바이너리 실체와 버전."""
    out = dict(
        found=False,
        path=UNKNOWN,
        reason="",
        version=UNKNOWN,
        realpath=UNKNOWN,
        is_symlink=None,
        bin_dir_entries=[],
        ldd=None,
        file=None,
        install=dict(guess=UNKNOWN, candidates=[], evidence=[]),
        search_hits=[],
        equivalent_commands=[
            "which mmseqs", "type -a mmseqs", "mmseqs version",
            'readlink -f "$(which mmseqs)"',
        ],
    )

    # 1) 사용자 지정 → 2) common.find_mmseqs (CHCODH2_MMSEQS → PATH → 탐색경로)
    binary = None
    if args.mmseqs:
        p = Path(args.mmseqs).expanduser()
        if p.is_file() and os.access(p, os.X_OK):
            binary = p
        else:
            out["reason"] = f"--mmseqs 로 받은 경로가 실행 가능한 파일이 아니다: {p}"
    if binary is None:
        try:
            binary = Path(common.find_mmseqs())
        except FileNotFoundError as e:
            out["reason"] = (out["reason"] + " / " if out["reason"] else "") + str(e)
        except Exception as e:  # noqa: BLE001
            out["reason"] = f"탐색 중 오류: {type(e).__name__}: {e}"

    # 5.2 직접 탐색 — 못 찾았을 때 어디를 봤는지 남긴다
    hits, searched = [], []
    for root in [args.mnt / "soft", Path.home(), Path("/opt"), Path("/usr/local")]:
        searched.append(str(root))
        if not root.is_dir():
            continue
        try:
            for dirpath, _dn, filenames in walk_limited(root, 5):
                for fn in filenames:
                    if fn.startswith("mmseqs"):
                        hits.append(str(dirpath / fn))
                        if len(hits) >= 40:
                            break
                if len(hits) >= 40:
                    break
        except OSError:
            continue
        if len(hits) >= 40:
            break
    out["search_hits"] = hits
    out["searched_dirs"] = searched
    out["search_hits_capped"] = len(hits) >= 40

    if binary is None:
        why = out["reason"] or "PATH·탐색경로에 없음"
        out["path"] = f"{UNKNOWN}(사유: {why})"
        out["version"] = f"{UNKNOWN}(사유: 바이너리를 못 찾아 `mmseqs version` 을 돌릴 수 없음)"
        out["realpath"] = f"{UNKNOWN}(사유: 바이너리 없음)"
        out["install"] = dict(guess=f"{UNKNOWN}(사유: 바이너리 없음)", candidates=[],
                              evidence=[f"탐색 실패: {why}"])
        return out

    out.update(found=True, path=str(binary))
    try:
        out["is_symlink"] = Path(binary).is_symlink()
        out["realpath"] = str(Path(binary).resolve())
    except OSError as e:
        out["realpath"] = f"{UNKNOWN}(사유: {e})"

    # ls -la $(dirname) | grep mmseqs 의 파이썬 버전
    try:
        for e in sorted(Path(binary).parent.iterdir()):
            if "mmseqs" in e.name.lower():
                st = e.lstat()
                out["bin_dir_entries"].append(dict(
                    name=e.name, bytes=st.st_size,
                    mode=oct(st.st_mode & 0o7777),
                    mtime=datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
                    symlink_to=(os.readlink(e) if e.is_symlink() else None),
                ))
    except OSError as e:
        out["bin_dir_entries"] = [dict(error=str(e))]

    out["ldd"] = probe_cmd(["ldd", str(binary)], dry_run=args.dry_run, log_dir=log_dir,
                           label="ldd_mmseqs", timeout=60)
    out["file"] = probe_cmd(["file", "-b", str(binary)], dry_run=args.dry_run, log_dir=log_dir,
                            label="file_mmseqs", timeout=60)
    out["readlink_f"] = probe_cmd(["readlink", "-f", str(binary)], dry_run=args.dry_run, timeout=20)

    if args.dry_run:
        out["version"] = f"{UNKNOWN}(사유: --dry-run 으로 실행 안 함)"
    else:
        ver = probe_cmd([str(binary), "version"], dry_run=False, log_dir=log_dir,
                        label="mmseqs_version", timeout=120)
        out["version_cmd"] = ver
        if ver["status"] == "ok":
            out["version"] = first_line(ver["stdout"] or ver["stderr"], UNKNOWN)
        else:
            out["version"] = f"{UNKNOWN}(사유: `mmseqs version` {ver['status']})"
        # common 경유 값도 남겨 두 경로가 어긋나는지 본다
        try:
            out["version_via_common"] = common.mmseqs_version(str(binary))
        except Exception as e:  # noqa: BLE001
            out["version_via_common"] = f"{UNKNOWN}({e})"

    out["install"] = infer_install_method(Path(binary), out["ldd"] or {}, out["file"] or {})
    return out


def detect_conda(args, log_dir: Path) -> dict:
    """인계문서 5.3 — conda 환경 목록과 각 환경의 mmseqs 유무."""
    out = dict(conda_exe=NOT_FOUND, env_list=None, envs=[], n_envs=0,
               n_envs_conda_listed=0, conda_list_capped=False, note="")
    exe = shutil.which("conda") or shutil.which("mamba") or shutil.which("micromamba")
    if not exe:
        # CONDA_EXE / CONDA_PREFIX 로 한 번 더
        for cand in [os.environ.get("CONDA_EXE"),
                     (os.environ.get("CONDA_PREFIX", "") + "/bin/conda") if os.environ.get("CONDA_PREFIX") else None]:
            if cand and Path(cand).is_file():
                exe = cand
                break
    if not exe:
        out["note"] = f"{NOT_FOUND} — conda/mamba/micromamba 실행파일을 찾지 못함 (PATH·CONDA_EXE 모두)"
        return out
    out["conda_exe"] = exe

    env_rec = probe_cmd([exe, "env", "list", "--json"], dry_run=args.dry_run,
                        log_dir=log_dir, label="conda_env_list", timeout=180)
    out["env_list"] = env_rec
    prefixes = []
    if env_rec["status"] == "ok":
        try:
            prefixes = json.loads(env_rec["stdout"]).get("envs", [])
        except (json.JSONDecodeError, AttributeError):
            out["note"] = "conda env list --json 출력을 json 으로 못 읽음 — 평문으로 재시도"
    if not prefixes and not args.dry_run:
        plain = probe_cmd([exe, "env", "list"], dry_run=False, timeout=180)
        out["env_list_plain"] = plain
        for ln in (plain["stdout"] or "").splitlines():
            if ln.startswith("#") or not ln.strip():
                continue
            parts = ln.split()
            prefixes.append(parts[-1])

    for pref in prefixes:
        p = Path(pref)
        rec = dict(name=(p.name or str(p)), prefix=str(p),
                   has_bin_mmseqs=(p / "bin/mmseqs").is_file(),
                   active=(str(p) == os.environ.get("CONDA_PREFIX")))
        out["envs"].append(rec)
    out["n_envs"] = len(out["envs"])

    # `conda list -n <env> | grep mmseqs` — 느리므로 상한을 둔다 (상한은 반드시 보고)
    budget = min(args.max_conda_envs, len(out["envs"]))
    out["conda_list_capped"] = len(out["envs"]) > budget
    for rec in out["envs"][:budget]:
        if args.dry_run:
            rec["conda_list_mmseqs"] = "dry-run"
            continue
        r = probe_cmd([exe, "list", "-p", rec["prefix"]], dry_run=False, timeout=240)
        if r["status"] == "ok":
            lines = [l for l in r["stdout"].splitlines() if "mmseqs" in l.lower()]
            rec["conda_list_mmseqs"] = lines or NOT_FOUND
        else:
            rec["conda_list_mmseqs"] = f"{UNKNOWN}({r['status']})"
        out["n_envs_conda_listed"] += 1
    for rec in out["envs"][budget:]:
        rec["conda_list_mmseqs"] = f"{UNKNOWN}(사유: --max-conda-envs {args.max_conda_envs} 상한)"
    return out


def detect_jupyter(args, log_dir: Path) -> dict:
    """인계문서 5.4 — jupyter kernelspec list. 6절 함정(커널 PATH)의 근거."""
    out = dict(jupyter_exe=NOT_FOUND, kernelspec=None, kernels=[], note="")
    exe = shutil.which("jupyter")
    if not exe:
        out["note"] = f"{NOT_FOUND} — jupyter 실행파일이 PATH 에 없음"
        return out
    out["jupyter_exe"] = exe
    rec = probe_cmd([exe, "kernelspec", "list", "--json"], dry_run=args.dry_run,
                    log_dir=log_dir, label="jupyter_kernelspec", timeout=180)
    out["kernelspec"] = rec
    if rec["status"] == "ok":
        try:
            specs = json.loads(rec["stdout"]).get("kernelspecs", {})
            for name, meta in specs.items():
                spec = (meta or {}).get("spec", {}) or {}
                argv = spec.get("argv", []) or []
                out["kernels"].append(dict(
                    name=name,
                    display_name=spec.get("display_name", UNKNOWN),
                    resource_dir=(meta or {}).get("resource_dir", UNKNOWN),
                    python=(argv[0] if argv else UNKNOWN),
                ))
        except (json.JSONDecodeError, AttributeError):
            out["note"] = "kernelspec --json 출력을 json 으로 못 읽음"
    if not out["kernels"] and not args.dry_run:
        out["kernelspec_plain"] = probe_cmd([exe, "kernelspec", "list"], dry_run=False, timeout=180)
    return out


def dir_size_capped(root: Path, max_entries: int) -> dict:
    """tmp 잔여물 용량 — 엔트리 상한에 걸리면 '하한값' 임을 표시한다."""
    total, n = 0, 0
    truncated = False
    try:
        for dirpath, _dn, filenames in os.walk(root, onerror=lambda e: None):
            for fn in filenames:
                try:
                    total += (Path(dirpath) / fn).lstat().st_size
                except OSError:
                    pass
                n += 1
                if n >= max_entries:
                    truncated = True
                    break
            if truncated:
                break
    except OSError:
        pass
    return dict(bytes=total, entries=n, is_lower_bound=truncated)


def detect_resources(args, log_dir: Path) -> dict:
    """인계문서 5.5 — nproc / df / tmp 잔여물."""
    out = dict(cpu_count=os.cpu_count(), nproc=UNKNOWN, mnt=str(args.mnt),
               mnt_exists=args.mnt.is_dir(), disk=None, df=None,
               tmp_leftovers=[], n_tmp_found=0, tmp_capped=False)
    r = probe_cmd(["nproc"], dry_run=args.dry_run, timeout=30)
    out["nproc"] = first_line(r["stdout"], f"{UNKNOWN}({r['status']})")
    out["nproc_cmd"] = r

    if out["mnt_exists"]:
        try:
            du = shutil.disk_usage(args.mnt)
            out["disk"] = dict(total_bytes=du.total, used_bytes=du.used, free_bytes=du.free,
                               free_gib=round(du.free / 2**30, 1),
                               total_gib=round(du.total / 2**30, 1))
        except OSError as e:
            out["disk"] = dict(error=str(e))
        out["df"] = probe_cmd(["df", "-h", str(args.mnt)], dry_run=args.dry_run,
                              log_dir=log_dir, label="df_mnt", timeout=60)
    else:
        out["disk"] = dict(error=f"{args.mnt} 가 없음 — 이 장비가 epel-af2 가 아니다")

    # tmp 잔여물: MNT 아래 maxdepth 4, --tmp-since 이후 수정된 'tmp*' 디렉터리
    since_ts = None
    try:
        since_ts = datetime.strptime(args.tmp_since, "%Y-%m-%d").replace(tzinfo=timezone.utc).timestamp()
    except ValueError:
        out["tmp_since_error"] = f"--tmp-since 형식 오류({args.tmp_since}) — 날짜 필터 없이 훑는다"
    found = []
    if out["mnt_exists"]:
        for dirpath, dirnames, _fn in walk_limited(args.mnt, 4):
            for d in list(dirnames):
                if not d.lower().startswith("tmp"):
                    continue
                p = dirpath / d
                try:
                    mt = p.stat().st_mtime
                except OSError:
                    continue
                if since_ts and mt < since_ts:
                    continue
                found.append(p)
    out["n_tmp_found"] = len(found)
    out["tmp_capped"] = len(found) > CAP["tmp_dirs"]
    for p in found[: CAP["tmp_dirs"]]:
        size = dir_size_capped(p, CAP["tmp_du_entries"])
        out["tmp_leftovers"].append(dict(
            path=str(p),
            mtime=datetime.fromtimestamp(p.stat().st_mtime, timezone.utc).isoformat(),
            size_gib=round(size["bytes"] / 2**30, 2),
            size_is_lower_bound=size["is_lower_bound"],
            entries_scanned=size["entries"],
        ))
    # 재현 실행에 쓸 tmp 위치의 여유도 본다 (easy-cluster tmp 는 수 GB 가 될 수 있다)
    probe = args.repro_dir
    while not probe.is_dir() and probe != probe.parent:
        probe = probe.parent
    try:
        out["repro_dir_free_gib"] = round(shutil.disk_usage(probe).free / 2**30, 1)
        out["repro_dir_checked"] = str(probe)
    except OSError as e:
        out["repro_dir_free_gib"] = f"{UNKNOWN}({e})"
    return out


# ===========================================================================
# (b) 219 재현
# ===========================================================================
def _read_cluster_tsv(path: Path) -> pd.DataFrame:
    """mmseqs _cluster.tsv — 헤더 없는 2열. 헤더 추정에 맡기면 첫 줄을 잃는다."""
    df = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if df.shape[1] != 2:
        raise ValueError(
            f"{path} 는 2열이어야 한다 (대표\\t멤버). 실제 {df.shape[1]}열 — "
            "mmseqs 산출물이 아니거나 구분자가 탭이 아니다."
        )
    df.columns = ["cluster", "accession"]
    return df


def reproduce_clusters(args, log_dir: Path, disc: list) -> dict:
    """인계문서 7절 — easy-cluster 를 그대로 돌려 219/1827 을 확인하고 멤버십 대조."""
    out = dict(requested=True, status=UNKNOWN, reason="", command=None,
               n_clusters=None, n_sequences=None, n_rows=None,
               expected_clusters=args.expect_clusters, expected_sequences=args.expect_sequences,
               clusters_match=None, sequences_match=None, compare=None,
               outputs={}, threads=(args.threads or "미지정(전 코어)"))

    mm = args.mmseqs_resolved
    if args.dry_run:
        out.update(status="건너뜀", reason="--dry-run")
    elif not mm:
        out.update(status="건너뜀", reason="mmseqs 바이너리를 찾지 못함")
    elif not args.ipr_fasta.is_file():
        out.update(status="건너뜀", reason=f"입력 FASTA 없음: {args.ipr_fasta}")
    if out["status"] == "건너뜀":
        disc.append(f"219 재현을 건너뛰었다 (사유: {out['reason']}) — 재현성은 미확인 상태로 남는다")
        return out

    # 입력 FASTA 의 서열 수를 먼저 센다 (클러스터 tsv 와 대조할 기준)
    n_fasta, dup_acc = 0, 0
    seen = set()
    try:
        for header, _seq in common.read_fasta(args.ipr_fasta):
            n_fasta += 1
            acc = common.parse_accession(header)
            if acc in seen:
                dup_acc += 1
            seen.add(acc)
    except OSError as e:
        raise RuntimeError(f"입력 FASTA 를 읽을 수 없다: {args.ipr_fasta} ({e})") from e
    out["n_fasta_records"] = n_fasta
    out["n_fasta_unique_accessions"] = len(seen)
    out["n_fasta_duplicate_accessions"] = dup_acc
    if dup_acc:
        disc.append(f"입력 FASTA 에 중복 accession {dup_acc} 건 — 서열 수 {n_fasta} ≠ 고유 {len(seen)}")

    repro = args.repro_dir
    prefix = repro / "id70"
    tmpdir = repro / "tmp"
    repro.mkdir(parents=True, exist_ok=True)
    cmd = [mm, "easy-cluster", str(args.ipr_fasta), str(prefix), str(tmpdir), *CLUSTER_OPTS]
    if args.threads:
        cmd += ["--threads", str(args.threads)]
        disc.append(f"--threads {args.threads} 를 추가로 주었다 — 인계문서에 기록된 명령과 다르다")
    out["command"] = " ".join(cmd)

    rec = common.run_cmd(cmd, log_path=log_dir / "mmseqs_easy_cluster.log", check=False,
                         timeout=args.repro_timeout)
    out["run"] = {k: rec[k] for k in ("cmd_str", "returncode", "seconds", "stderr_tail")}
    if rec["returncode"] != 0:
        out.update(status="실패", reason=f"easy-cluster rc={rec['returncode']}")
        disc.append(f"easy-cluster 실패 (rc={rec['returncode']}) — stderr 는 logs/ 참조")
        return out

    tsv = Path(f"{prefix}_cluster.tsv")
    if not tsv.is_file():
        out.update(status="실패", reason=f"산출물이 없다: {tsv}")
        disc.append(f"easy-cluster 가 rc=0 인데 {tsv} 가 없다 — mmseqs 동작이 예상과 다르다")
        return out

    new = _read_cluster_tsv(tsv)
    out["n_rows"] = int(len(new))
    out["n_clusters"] = int(new["cluster"].nunique())
    out["n_sequences"] = int(new["accession"].nunique())
    out["clusters_match"] = (out["n_clusters"] == args.expect_clusters)
    out["sequences_match"] = (out["n_sequences"] == args.expect_sequences)
    out["outputs"] = {s: str(Path(f"{prefix}{s}")) for s in
                      ("_cluster.tsv", "_rep_seq.fasta", "_all_seqs.fasta")}
    if not out["clusters_match"]:
        disc.append(f"클러스터 수 {out['n_clusters']} ≠ 기대 {args.expect_clusters} "
                    "— mmseqs 버전 또는 입력이 달라졌다")
    if not out["sequences_match"]:
        disc.append(f"서열 수 {out['n_sequences']} ≠ 기대 {args.expect_sequences}")
    if out["n_sequences"] != len(seen):
        disc.append(f"클러스터 tsv 의 고유 멤버 {out['n_sequences']} ≠ FASTA 고유 accession {len(seen)}")
    out["status"] = "ok"

    out["compare"] = compare_with_master(args, new, disc)
    if not args.keep_tmp:
        shutil.rmtree(tmpdir, ignore_errors=True)
        out["tmp_removed"] = str(tmpdir)
    return out


def compare_with_master(args, new: pd.DataFrame, disc: list) -> dict:
    """기존 master_1827.tsv 와 멤버십 대조. 대표 ID 동일과 분할 1:1 동일을 구분한다."""
    cmp = dict(master=str(args.master_tsv), status=UNKNOWN)
    if not args.master_tsv.is_file():
        cmp.update(status="건너뜀", reason=f"master_1827.tsv 없음: {args.master_tsv}")
        return cmp

    old_raw = common.read_tsv(args.master_tsv)
    cols = common.resolve_columns(old_raw, ["accession", "cluster"], required=("accession", "cluster"))
    old = old_raw[[cols["accession"], cols["cluster"]]].copy()
    old.columns = ["accession", "cluster_old"]
    old["accession"] = common.norm_acc(old["accession"])
    old["cluster_old"] = common.norm_acc(old["cluster_old"])

    nw = new.copy()
    nw["accession"] = common.norm_acc(nw["accession"])
    nw["cluster_new"] = common.norm_acc(nw["cluster"])
    nw = nw[["accession", "cluster_new"]]

    cmp["master_columns_used"] = cols
    cmp["n_master_rows"] = int(len(old))
    cmp["n_new_rows"] = int(len(nw))

    # 중복은 조인 행 수를 부풀려 조용히 숫자를 망친다 — 세고, 버린 수를 보고한다
    d_old = int(old["accession"].duplicated().sum())
    d_new = int(nw["accession"].duplicated().sum())
    cmp["n_dup_accession_master"] = d_old
    cmp["n_dup_accession_new"] = d_new
    if d_old or d_new:
        disc.append(f"accession 중복: master {d_old} 건 / 신규 {d_new} 건 — 조인 전 첫 행만 남겼다")
        old = old.drop_duplicates("accession", keep="first")
        nw = nw.drop_duplicates("accession", keep="first")
    cmp["n_rows_dropped_for_join"] = d_old + d_new

    j = old.merge(nw, on="accession", how="outer")
    total = int(len(j))
    matched = j.dropna(subset=["cluster_old", "cluster_new"])
    cmp.update(
        n_join_rows=total,
        n_matched=int(len(matched)),
        join_rate=round(len(matched) / total, 6) if total else None,
        n_only_master=int(j["cluster_new"].isna().sum()),
        n_only_new=int(j["cluster_old"].isna().sum()),
        missing_rate_new=round(float(j["cluster_new"].isna().mean()), 6) if total else None,
        missing_rate_master=round(float(j["cluster_old"].isna().mean()), 6) if total else None,
        n_clusters_master=int(j["cluster_old"].nunique()),
        n_clusters_new=int(j["cluster_new"].nunique()),
        n_same_rep_id=int((j["cluster_old"] == j["cluster_new"]).sum()),
    )
    cmp["same_rep_id_rate"] = round(cmp["n_same_rep_id"] / total, 6) if total else None
    if cmp["n_only_master"] or cmp["n_only_new"]:
        disc.append(f"조인 불일치: master 전용 {cmp['n_only_master']} / 신규 전용 {cmp['n_only_new']}")

    if len(matched) == 0:
        cmp["status"] = "대조 불가(교집합 0)"
        disc.append("master 와 신규 클러스터링의 accession 교집합이 0 — 정규화 규칙이 어긋났을 수 있다")
        return cmp

    cm = common.crosstab_mapping(matched["cluster_old"], matched["cluster_new"])
    cmp["crosstab"] = cm
    cmp["ari"] = round(common.adjusted_rand_index(matched["cluster_old"], matched["cluster_new"]), 6)
    cmp["nmi"] = round(common.normalized_mutual_info(matched["cluster_old"], matched["cluster_new"]), 6)
    # 분할이 같다 = 쪼개짐도 합쳐짐도 없고 클러스터 개수가 같다 (대표 ID 와 무관)
    cmp["partition_identical"] = bool(
        cm["n_split"] == 0 and cm["n_merged"] == 0 and cm["n_old"] == cm["n_new"]
        and len(matched) == total
    )
    cmp["rep_ids_identical"] = bool(cmp["n_same_rep_id"] == total)
    cmp["status"] = "ok"
    if not cmp["partition_identical"]:
        disc.append(f"분할이 달라졌다: 쪼개진 클러스터 {cm['n_split']} / 합쳐진 클러스터 "
                    f"{cm['n_merged']} / old {cm['n_old']} vs new {cm['n_new']} (ARI={cmp['ari']})")
    return cmp


# ===========================================================================
# (c) search 명령 복원
# ===========================================================================
def align_columns(args) -> dict:
    """chcodh2_alignment.tsv 의 컬럼 구성 — 복원 후보 채점의 정답지."""
    out = dict(path=str(args.align_tsv), status=UNKNOWN, columns=None, logical=None, n_cols=None)
    if not args.align_tsv.is_file():
        out.update(status="없음", reason=f"파일 없음: {args.align_tsv}")
        return out
    try:
        df = common.read_tsv(args.align_tsv, nrows=5)
    except Exception as e:  # noqa: BLE001  채점용 보조 입력이므로 치명적이지 않다
        out.update(status="읽기 실패", reason=f"{type(e).__name__}: {e}")
        return out
    cols = [str(c) for c in df.columns]
    out.update(status="ok", columns=cols, n_cols=len(cols))
    out["headerless"] = all(re.fullmatch(r"col\d+", c) for c in cols)
    if out["headerless"]:
        out["note"] = ("헤더가 없는 파일로 읽혔다 — 컬럼명을 알 수 없어 "
                       "--format-output 대조는 '열 개수'만 쓴다")
    else:
        try:
            out["logical"] = common.resolve_columns(
                df, ["accession", "fident", "qstart", "qaln", "tstart", "taln", "qcov", "tcov"]
            )
        except Exception as e:  # noqa: BLE001
            out["logical_error"] = str(e)
    return out


def _mk_candidate(source, kind, locator, line_no, text, context, tier, looks_cmd, aliases):
    return dict(source=str(source), source_kind=kind, locator=locator, line_no=line_no,
                text=text[: CAP["cand_text_chars"]],
                text_truncated=len(text) > CAP["cand_text_chars"],
                context=context[: CAP["cand_text_chars"]],
                match_tier=tier, looks_like_command=looks_cmd,
                aliases_in_unit=sorted(aliases),
                format_output=None, score=0, confidence=UNKNOWN, reasons=[])


def collect_aliases(text: str):
    """`MM = "/절대/경로/mmseqs"` 처럼 mmseqs 를 변수에 담은 흔적을 찾는다.

    인계문서 6절이 **권장한 패턴**이라 실제 호출 줄에는 'mmseqs' 문자열이 없다
    (`!$MM convertalis ...`, `subprocess.run([MM, "search", ...])`).
    'mmseqs' 직접 매칭만 하면 복원해야 할 바로 그 줄을 놓친다.
    """
    aliases, assigns = set(), []
    for i, ln in enumerate(text.splitlines()):
        if "mmseqs" not in ln.lower():
            continue
        m = RE_ALIAS_ASSIGN.match(ln.lstrip("!%").lstrip())
        if not m:
            continue
        name = m.group(1)
        if len(name) < 2 or name.lower() == "mmseqs" or name in ALIAS_STOPWORDS:
            continue
        aliases.add(name)
        assigns.append(dict(line_no=i + 1, text=ln.strip()[: CAP["cand_text_chars"]], alias=name))
    return aliases, assigns


def _scan_unit(lines, source, kind, locator, cands, aliases, unit_has_mmseqs, *, ctx=3):
    """한 단위(파일 또는 셀)를 줄 단위로 훑는다. 명령이 이어질 수 있어 문맥도 담는다.

    세 가지 경로로 후보를 잡는다.
      직접        'mmseqs' + 서브명령이 같은 줄
      별칭        서브명령 + mmseqs 를 담은 변수명 (같은 파일/노트북에서 수집)
      전용옵션    서브명령이 mmseqs 전용(easy-search/convertalis/format-output)이고
                  그 파일/노트북이 mmseqs 를 언급한다
    """
    n0 = len(cands)
    alias_use = (re.compile(r"\b(?:" + "|".join(re.escape(a) for a in sorted(aliases)) + r")\b")
                 if aliases else None)
    names = ["mmseqs"] + sorted(aliases)
    cmd_re = re.compile(
        r"(?:\$\{?)?\b(?:" + "|".join(re.escape(n) for n in names) + r")\b\}?"
        r"[\"'\s,\]\[)]+(?:easy-search|search|convertalis|easy-cluster|createdb)",
        re.IGNORECASE,
    )
    for i, ln in enumerate(lines):
        if not RE_SUBCMD.search(ln):
            continue
        if RE_MMSEQS.search(ln):
            tier = "직접"
        elif alias_use and alias_use.search(ln):
            tier = "별칭"
        elif unit_has_mmseqs and RE_SPECIFIC.search(ln):
            tier = "전용옵션"
        else:
            continue
        lo, hi = max(0, i - ctx), min(len(lines), i + ctx + 1)
        cands.append(_mk_candidate(source, kind, locator, i + 1, ln.strip(),
                                   "\n".join(lines[lo:hi]), tier,
                                   bool(cmd_re.search(ln)), aliases))
    return len(cands) - n0


def scan_notebook(path: Path, cands: list, stats: dict, alias_log: list) -> None:
    """ipynb 는 json. nbformat 이 없으므로 표준 json 으로 읽는다. 셀 소스 + 셀 출력 모두."""
    text, why = read_text_capped(path, CAP["scan_file_mb"] * 10**6)
    if text is None:
        stats["notebooks_skipped"].append(dict(path=str(path), reason=why))
        return
    try:
        nb = json.loads(text)
    except (json.JSONDecodeError, UnicodeDecodeError) as e:
        stats["notebooks_json_error"].append(dict(path=str(path), reason=str(e)))
        return
    stats["n_notebooks_parsed"] += 1

    cells = nb.get("cells")
    if cells is None:  # nbformat 3 이하: worksheets[0].cells
        ws = nb.get("worksheets") or []
        cells = (ws[0].get("cells") if ws else []) or []

    # 별칭은 다른 셀(보통 맨 앞)에서 정의된다 — 노트북 전체에서 먼저 모은다
    sources = []
    for cell in cells:
        if not isinstance(cell, dict):
            continue
        src = cell.get("source") or cell.get("input") or ""
        sources.append("".join(src) if isinstance(src, list) else str(src))
    nb_text = "\n".join(sources)
    aliases, assigns = collect_aliases(nb_text)
    unit_has_mmseqs = "mmseqs" in nb_text.lower()
    for a in assigns:
        alias_log.append(dict(source=str(path), kind="notebook", **a))

    for idx, cell in enumerate(cells):
        if not isinstance(cell, dict):
            continue
        src = cell.get("source") or cell.get("input") or ""
        if isinstance(src, list):
            src = "".join(src)
        exec_cnt = cell.get("execution_count", cell.get("prompt_number"))
        loc = f"cell {idx}" + (f" (exec {exec_cnt})" if exec_cnt not in (None, "") else "")
        if src:
            _scan_unit(src.splitlines(), path, "notebook_cell", loc, cands,
                       aliases, unit_has_mmseqs)
        for oi, outp in enumerate(cell.get("outputs") or []):
            if not isinstance(outp, dict):
                continue
            chunks = []
            for key in ("text", "traceback"):
                v = outp.get(key)
                if isinstance(v, list):
                    chunks.append("".join(str(x) for x in v))
                elif isinstance(v, str):
                    chunks.append(v)
            data = outp.get("data") or {}
            for key in ("text/plain", "text/html"):
                v = data.get(key)
                if isinstance(v, list):
                    chunks.append("".join(str(x) for x in v))
                elif isinstance(v, str):
                    chunks.append(v)
            blob = "\n".join(chunks)
            if blob:
                _scan_unit(blob.splitlines(), path, "notebook_output",
                           f"{loc} output {oi}", cands, aliases, unit_has_mmseqs)


def scan_text(path: Path, cands: list, stats: dict, kind: str, alias_log: list) -> None:
    text, why = read_text_capped(path, CAP["scan_file_mb"] * 10**6)
    if text is None:
        stats["text_skipped"].append(dict(path=str(path), reason=why))
        return
    lines = text.splitlines()
    if len(lines) > CAP["history_lines"]:
        stats["text_truncated"].append(dict(path=str(path), n_lines=len(lines),
                                            kept=CAP["history_lines"]))
        lines = lines[-CAP["history_lines"]:]
    stats["n_text_scanned"] += 1
    stats["n_lines_scanned"] += len(lines)
    aliases, assigns = collect_aliases("\n".join(lines))
    for a in assigns:
        alias_log.append(dict(source=str(path), kind=kind, **a))
    _scan_unit(lines, path, kind, "line", cands, aliases,
               any("mmseqs" in l.lower() for l in lines))


def scan_ipython_sqlite(path: Path, cands: list, stats: dict, alias_log: list) -> None:
    """IPython 히스토리는 sqlite. 노트북 셀이 지워졌어도 여기 남아 있을 수 있다."""
    try:
        uri = f"file:{path}?mode=ro"
        con = sqlite3.connect(uri, uri=True, timeout=5)
        try:
            rows = con.execute(
                "SELECT session, line, source FROM history LIMIT ?", (CAP["ipython_rows"],)
            ).fetchall()
        finally:
            con.close()
    except sqlite3.Error as e:
        stats["text_skipped"].append(dict(path=str(path), reason=f"sqlite: {e}"))
        return
    stats["n_text_scanned"] += 1
    stats["n_ipython_rows"] = len(rows)
    if len(rows) >= CAP["ipython_rows"]:
        stats["text_truncated"].append(dict(path=str(path), n_lines=">=" + str(len(rows)),
                                            kept=CAP["ipython_rows"]))
    all_text = "\n".join((r[2] or "") for r in rows)
    aliases, assigns = collect_aliases(all_text)
    for a in assigns:
        alias_log.append(dict(source=str(path), kind="ipython_history", **a))
    unit_has_mmseqs = "mmseqs" in all_text.lower()
    for sess, line, source in rows:
        src = source or ""
        stats["n_lines_scanned"] += src.count("\n") + 1
        _scan_unit(src.splitlines(), path, "ipython_history",
                   f"session {sess} line {line}", cands, aliases, unit_has_mmseqs)


def score_candidate(c: dict, acols: dict) -> None:
    """후보 신뢰도. chcodh2_alignment.tsv 의 컬럼 구성과 맞는지가 가장 센 근거다."""
    blob = (c["text"] + "\n" + c["context"])
    low = blob.lower()
    score, reasons = 0, []

    if c.get("looks_like_command") or RE_LOOKS_CMD.search(blob):
        score += 10
        reasons.append("실행 이름 바로 뒤에 서브명령이 온다 — 실제 호출로 보인다")
    else:
        reasons.append("명령 꼴이 아니다(주석·산문일 수 있다)")
    if c.get("match_tier") == "별칭":
        score += 8
        reasons.append(f"'mmseqs' 문자열 없이 별칭 {', '.join(c.get('aliases_in_unit') or [])} "
                       "로 호출 (인계문서 6절이 권장한 패턴)")
    elif c.get("match_tier") == "전용옵션":
        score += 4
        reasons.append("mmseqs 전용 토큰(easy-search/convertalis/format-output)으로 잡았다")

    m_text = RE_FORMAT_OUT.search(c["text"])
    m = m_text or RE_FORMAT_OUT.search(blob)
    c["format_output_from"] = ("본문" if m_text else ("문맥" if m else None))
    if m_text:
        score += 5
        reasons.append("이 줄 자체가 --format-output 을 들고 있다 (문맥이 아니라 본문)")
    if m:
        raw = m.group(1) or m.group(2) or ""
        toks = [t.strip() for t in raw.split(",") if t.strip()]
        c["format_output"] = toks
        score += 15
        reasons.append(f"--format-output 이 있다 ({len(toks)} 개 필드, 출처 {c['format_output_from']})")
        if acols.get("status") == "ok" and acols.get("columns"):
            cols = [str(x) for x in acols["columns"]]
            if acols.get("headerless"):
                if len(toks) == acols["n_cols"]:
                    score += 20
                    reasons.append(f"필드 수가 정렬 파일 열 수와 같다 ({len(toks)})")
                else:
                    reasons.append(f"필드 수 {len(toks)} ≠ 정렬 파일 열 수 {acols['n_cols']}")
            else:
                st, sc = {t.lower() for t in toks}, {x.lower() for x in cols}
                if st == sc:
                    score += 40
                    reasons.append("--format-output 필드가 정렬 파일 컬럼과 정확히 일치")
                elif sc and sc <= st:
                    score += 30
                    reasons.append("--format-output 이 정렬 파일 컬럼을 모두 포함")
                elif sc and len(st & sc) / len(sc) >= 0.7:
                    score += 20
                    reasons.append(f"컬럼 겹침 {len(st & sc)}/{len(sc)}")
                else:
                    reasons.append(f"컬럼 겹침 낮음 {len(st & sc)}/{len(sc) or 1}")
        else:
            reasons.append("정렬 파일을 못 읽어 컬럼 대조는 못 했다")
    else:
        reasons.append("--format-output 이 없다 — qaln/taln 이 나오지 않는다(2절)")

    if any(k in low for k in ("q9f8a8", "chcodh2", "chcodh", "codh2")):
        score += 10
        reasons.append("질의(ChCODH2/Q9F8A8) 표식 있음")
    if "chcodh2_alignment" in low:
        score += 15
        reasons.append("산출 파일명 chcodh2_alignment 이 보인다")
    if "ipr010047" in low:
        score += 5
        reasons.append("IPR010047 입력 표식 있음")
    if "easy-search" in low:
        score += 5
        reasons.append("easy-search (createdb+search+convertalis 일체형)")
    if "convertalis" in low:
        score += 5
        reasons.append("convertalis — 포맷을 직접 지정한 흐름")
    if c["source_kind"].startswith("notebook"):
        score += 5
        reasons.append(f"출처가 노트북 ({c['locator']})")

    paths = RE_PATHISH.findall(blob)
    if paths:
        c["seq_paths"] = sorted(set(RE_PATHISH.findall(blob)))[:8]
        reasons.append(f"서열 파일 경로 후보 {len(set(c['seq_paths']))} 개")

    c["score"] = score
    c["confidence"] = "높음" if score >= 60 else ("중간" if score >= 35 else "낮음")
    c["reasons"] = reasons


def recover_search(args, acols: dict, disc: list) -> dict:
    """인계문서 2절 — 기록되지 않은 search 명령을 노트북·히스토리에서 복원한다."""
    stats = dict(n_notebooks_found=0, n_notebooks_parsed=0, n_text_scanned=0,
                 n_lines_scanned=0, notebooks_skipped=[], notebooks_json_error=[],
                 text_skipped=[], text_truncated=[])
    out = dict(requested=True, status=UNKNOWN, scan_dirs=[], history_files=[],
               key_notebooks={}, candidates=[], n_candidates=0, stats=stats,
               align_columns=acols)

    # --- 스캔 대상 수집
    notebooks, texts = [], []
    for d in args.scan_dir:
        d = Path(d)
        rec = dict(path=str(d), exists=d.is_dir())
        out["scan_dirs"].append(rec)
        if not d.is_dir():
            rec["reason"] = "디렉터리 없음"
            continue
        n_nb = n_tx = 0
        for dirpath, _dn, filenames in walk_limited(d, CAP["scan_depth"]):
            for fn in filenames:
                p = dirpath / fn
                suf = p.suffix.lower()
                if suf == ".ipynb":
                    notebooks.append(p)
                    n_nb += 1
                elif suf in TEXT_EXTS or (args.include_md and suf == ".md"):
                    texts.append(p)
                    n_tx += 1
        rec.update(n_notebooks=n_nb, n_text_files=n_tx)

    stats["n_notebooks_found"] = len(notebooks)
    stats["n_text_files_found"] = len(texts)
    stats["notebooks_capped"] = len(notebooks) > CAP["scan_notebooks"]
    stats["text_files_capped"] = len(texts) > CAP["scan_text_files"]
    stats["n_notebooks_dropped_by_cap"] = max(0, len(notebooks) - CAP["scan_notebooks"])
    stats["n_text_dropped_by_cap"] = max(0, len(texts) - CAP["scan_text_files"])

    # 핵심 노트북 존재 확인 (인계문서 9절)
    for name in KEY_NOTEBOOKS:
        hit = [str(p) for p in notebooks if p.name == name]
        out["key_notebooks"][name] = hit or NOT_FOUND
        if not hit:
            disc.append(f"핵심 노트북 {name} 을 스캔 범위에서 찾지 못했다")

    cands: list = []
    alias_log: list = []
    for p in notebooks[: CAP["scan_notebooks"]]:
        scan_notebook(p, cands, stats, alias_log)
    for p in texts[: CAP["scan_text_files"]]:
        scan_text(p, cands, stats, "text_file", alias_log)

    # --- 히스토리
    for h in args.history:
        h = Path(h).expanduser()
        rec = dict(path=str(h), exists=h.is_file())
        out["history_files"].append(rec)
        if not h.is_file():
            rec["reason"] = NOT_FOUND
            continue
        before = len(cands)
        if h.suffix == ".sqlite":
            scan_ipython_sqlite(h, cands, stats, alias_log)
        else:
            scan_text(h, cands, stats, "shell_history", alias_log)
        rec["n_candidates"] = len(cands) - before

    # --- 채점·정렬·중복 제거
    for c in cands:
        score_candidate(c, acols)
    seen, uniq = set(), []
    for c in sorted(cands, key=lambda x: (-x["score"], x["source"], x["line_no"] or 0)):
        key = (re.sub(r"\s+", " ", c["text"]).strip().lower(), c["source_kind"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(c)
    # mmseqs 를 담은 변수 대입 줄 — 당시 쓴 절대경로가 여기 남아 있을 수 있다
    seen_alias = set()
    out["alias_assignments"] = []
    for a in alias_log:
        key = (a["source"], a["text"])
        if key in seen_alias:
            continue
        seen_alias.add(key)
        out["alias_assignments"].append(a)
    out["n_alias_assignments"] = len(out["alias_assignments"])

    out["n_candidates"] = len(uniq)
    out["n_candidates_raw"] = len(cands)
    out["n_by_tier"] = {t: sum(1 for c in uniq if c["match_tier"] == t)
                        for t in ("직접", "별칭", "전용옵션")}
    out["n_duplicates_merged"] = len(cands) - len(uniq)
    out["candidates_capped"] = len(uniq) > CAP["candidates_json"]
    out["n_candidates_dropped_by_cap"] = max(0, len(uniq) - CAP["candidates_json"])
    out["candidates"] = uniq[: CAP["candidates_json"]]
    out["n_high"] = sum(1 for c in uniq if c["confidence"] == "높음")
    out["n_mid"] = sum(1 for c in uniq if c["confidence"] == "중간")
    out["status"] = "ok" if uniq else "후보 없음"
    if not uniq:
        disc.append("search 명령 후보를 하나도 찾지 못했다 — 2절의 명령은 여전히 미복원이다")
    elif not out["n_high"]:
        disc.append(f"search 명령 후보 {len(uniq)} 개 중 신뢰도 '높음' 이 없다 — 사람이 확인해야 한다")

    # 질의 서열 출처 추정 — 최상위 후보(와 그 문맥)에 나온 경로만. 후보를 섞으면 오히려 틀린다
    top = out["candidates"][0] if out["candidates"] else None
    out["query_source_guess"] = (sorted(set(top.get("seq_paths") or []))[:5] or None) if top else None
    out["n_candidates_with_paths"] = sum(1 for c in out["candidates"] if c.get("seq_paths"))
    return out


# ===========================================================================
# (d) 산출물
# ===========================================================================
def _pad_wide(s: str, width: int) -> str:
    """한글은 두 칸을 차지한다. 8절 블록을 눈으로 읽을 수 있게 폭을 맞춘다."""
    w = sum(2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1 for ch in s)
    return s + " " * max(1, width - w)


def _fill(value, reason: str) -> str:
    """빈칸은 빈칸으로 두지 않는다. 못 채우면 사유를 적는다."""
    if value in (None, "", [], {}) or str(value).startswith(UNKNOWN):
        return f"{UNKNOWN}(사유: {reason})" if not str(value or "").startswith(UNKNOWN) else str(value)
    return str(value)


def build_blank_block(res: dict) -> str:
    """인계문서 8절 블록을 채운 상태로 만든다. selection_params.txt 에 붙일 수 있게."""
    mm, cd, jp, rs = res["mmseqs"], res["conda"], res["jupyter"], res["resources"]
    py, rep, sr = res["python"], res["reproduce"], res["search_recovery"]

    path = mm["path"] if mm["found"] else _fill(None, mm.get("reason") or "PATH·탐색경로에 없음")
    version = _fill(mm["version"], "mmseqs 없음")
    install = _fill(mm["install"]["guess"], "바이너리를 못 찾아 판정 불가")
    if mm["install"].get("evidence"):
        install += f"  [근거: {'; '.join(mm['install']['evidence'][:2])}]"

    if cd["conda_exe"] == NOT_FOUND:
        conda_line = _fill(None, "conda 실행파일 없음")
    else:
        active = [e["name"] for e in cd["envs"] if e.get("active")]
        with_mm = [e["name"] for e in cd["envs"] if e.get("has_bin_mmseqs")]
        conda_line = (f"활성={active[0] if active else (py['CONDA_DEFAULT_ENV'])} · "
                      f"환경 {cd['n_envs']}개 · bin/mmseqs 보유={with_mm or '없음'}")

    if jp["jupyter_exe"] == NOT_FOUND:
        kernel_line = _fill(None, "jupyter 실행파일 없음")
    elif jp["kernels"]:
        kernel_line = " · ".join(f"{k['name']}({k['python']})" for k in jp["kernels"][:4])
        if len(jp["kernels"]) > 4:
            kernel_line += f" … 외 {len(jp['kernels'])-4}개"
    else:
        kernel_line = _fill(None, "kernelspec 목록을 읽지 못함")

    python_line = f"{py['python_version']} ({py['python_executable']})"
    threads_line = (f"nproc={rs['nproc']} · 기록된 명령에 --threads 없음 → 전 코어 사용"
                    if not res["args"]["threads"] else
                    f"nproc={rs['nproc']} · 이번 실행 --threads {res['args']['threads']}")

    if rs["tmp_leftovers"]:
        t = rs["tmp_leftovers"][0]
        tmp_line = (f"잔여물 {rs['n_tmp_found']}개 · 최상위 {t['path']} "
                    f"({t['size_gib']} GiB{'+' if t['size_is_lower_bound'] else ''}, {t['mtime'][:10]})")
    elif rs["mnt_exists"]:
        tmp_line = f"{rs['mnt']} 아래 tmp* 잔여물 없음 (재현 실행은 {res['args']['repro_dir']} 사용)"
    else:
        tmp_line = _fill(None, f"{rs['mnt']} 가 없어 확인 불가")

    best = sr["candidates"][0] if sr.get("candidates") else None
    if best:
        cmd_line1 = re.sub(r"\s+", " ", best["text"]).strip()[:300]
        cmd_line2 = (f"출처: {best['source']} / {best['locator']} "
                     f"(줄 {best['line_no']}) · 신뢰도 {best['confidence']}(점수 {best['score']})")
    else:
        cmd_line1 = _fill(None, sr.get("status") or "후보 없음")
        cmd_line2 = f"스캔: 노트북 {sr['stats']['n_notebooks_parsed']}개 / 텍스트·히스토리 {sr['stats']['n_text_scanned']}개"

    qsrc = sr.get("query_source_guess")
    query_line = (", ".join(qsrc) + "  [최상위 후보 기준 — 다른 후보의 경로는 6절]"
                  if qsrc else _fill(None, "복원 후보에서 질의 파일 경로를 못 찾음"))

    if rep["status"] == "ok":
        c = rep.get("compare") or {}
        verdict = (f"클러스터 {rep['n_clusters']} / 서열 {rep['n_sequences']} "
                   f"(기대 {rep['expected_clusters']}/{rep['expected_sequences']}: "
                   f"{'일치' if rep['clusters_match'] and rep['sequences_match'] else '불일치'})")
        if c.get("status") == "ok":
            verdict += (f" · 대표 ID 동일 {c['n_same_rep_id']}/{c['n_join_rows']}"
                        f" · 분할 1:1 {'동일' if c['partition_identical'] else '다름'}"
                        f" (ARI={c['ari']})")
        elif c.get("status"):
            verdict += f" · 멤버십 대조: {c['status']}"
        repro_line = verdict
    else:
        repro_line = _fill(None, f"{rep['status']}: {rep['reason']}")

    rows = [
        ("mmseqs 경로", path),
        ("mmseqs 버전", version),
        ("설치 방식", install),
        ("conda 환경", conda_line),
        ("Jupyter 커널", kernel_line),
        ("Python", python_line),
        ("스레드", threads_line),
        ("tmp 위치", tmp_line),
        ("", ""),
        ("search 명령 (2절)", cmd_line1),
        ("", cmd_line2),
        ("질의 서열 출처", query_line),
        ("", ""),
        ("219 재현 결과", repro_line),
    ]
    lines = []
    for k, v in rows:
        if not k and not v:
            lines.append("")
        else:
            lines.append(_pad_wide(k, 19) + v)
    return "\n".join(lines)


def render_md(res: dict, disc: list) -> str:
    mm, cd, jp, rs = res["mmseqs"], res["conda"], res["jupyter"], res["resources"]
    py, rep, sr = res["python"], res["reproduce"], res["search_recovery"]
    L = []
    A = L.append
    A("# 단계 0 — mmseqs 실행 환경 탐지 (env_probe)")
    A("")
    A("| | |")
    A("|---|---|")
    A(f"| 생성 | {res['generated_utc']} |")
    A(f"| 호스트 · 사용자 | `{py['host']}` · `{py['user']}` |")
    A(f"| 모드 | {'dry-run (외부 명령 실행 안 함)' if res['args']['dry_run'] else '실행'}"
      f" · reproduce={res['args']['reproduce']} · recover-search={res['args']['recover_search']} |")
    A(f"| 불일치/미확인 | **{len(disc)} 건** (아래 마지막 절) |")
    A(f"| 근거 문서 | 20261021_handoff_mmseqs_env.md 5·7·8절 및 2절 |")
    A("")
    A("## 1. 인계문서 8절 빈칸 — 채운 블록")
    A("")
    A("아래를 그대로 `$SPLIT/selection_params.txt` 에 붙이면 된다.")
    A("")
    A("```")
    A(build_blank_block(res))
    A("```")
    A("")

    A("## 2. mmseqs 바이너리")
    A("")
    A("| 항목 | 값 |")
    A("|---|---|")
    A(f"| 찾음 | {'예' if mm['found'] else '아니오'} |")
    A(f"| 경로 | `{mm['path']}` |")
    A(f"| readlink -f | `{mm['realpath']}` |")
    A(f"| 심볼릭 링크 | {mm['is_symlink']} |")
    A(f"| 버전 | {mm['version']} |")
    A(f"| 설치 방식 추정 | {mm['install']['guess']} |")
    for ev in mm["install"]["evidence"]:
        A(f"| └ 근거 | {ev} |")
    A(f"| PATH 상의 mmseqs | `{py['mmseqs_on_path']}` |")
    A(f"| 직접 탐색 결과 | {len(mm['search_hits'])} 건"
      f"{' (상한 40 에 걸림)' if mm.get('search_hits_capped') else ''} |")
    A("")
    if mm["search_hits"]:
        A("탐색된 mmseqs* 파일:")
        A("")
        for h in mm["search_hits"][:10]:
            A(f"- `{h}`")
        if len(mm["search_hits"]) > 10:
            A(f"- … 외 {len(mm['search_hits'])-10} 건")
        A("")
    else:
        A(f"탐색 대상 디렉터리: {', '.join('`'+d+'`' for d in mm.get('searched_dirs', []))} — 모두 적중 없음")
        A("")

    A("## 3. conda / Jupyter / Python")
    A("")
    A("| 항목 | 값 |")
    A("|---|---|")
    A(f"| conda 실행파일 | `{cd['conda_exe']}` |")
    A(f"| conda 환경 수 | {cd['n_envs']} (conda list 조회 {cd['n_envs_conda_listed']}개"
      f"{', 상한으로 일부 미조회' if cd.get('conda_list_capped') else ''}) |")
    A(f"| CONDA_PREFIX | `{py['CONDA_PREFIX']}` |")
    A(f"| CONDA_DEFAULT_ENV | `{py['CONDA_DEFAULT_ENV']}` |")
    A(f"| jupyter 실행파일 | `{jp['jupyter_exe']}` |")
    A(f"| 커널 수 | {len(jp['kernels'])} |")
    A(f"| python | `{py['python_executable']}` ({py['python_version']}) |")
    A(f"| pandas / numpy | {py['pandas']} / {py['numpy']} |")
    A(f"| PATH 항목 수 | {py['PATH_entries']} |")
    A(f"| LD_LIBRARY_PATH | `{py['LD_LIBRARY_PATH']}` |")
    A(f"| MPLBACKEND | `{py['MPLBACKEND']}` (ColabFold 누수 이력 — 인계문서 6절) |")
    A("")
    if cd["envs"]:
        A("| conda 환경 | prefix | bin/mmseqs | conda list |")
        A("|---|---|---|---|")
        for e in cd["envs"][:20]:
            cl = e.get("conda_list_mmseqs")
            cl = ("; ".join(cl) if isinstance(cl, list) else str(cl))[:60]
            A(f"| {e['name']}{' **(활성)**' if e.get('active') else ''} | `{e['prefix']}` | "
              f"{'있음' if e['has_bin_mmseqs'] else '없음'} | {cl} |")
        if len(cd["envs"]) > 20:
            A(f"| … 외 {len(cd['envs'])-20}개 | | | |")
        A("")
    if cd.get("note"):
        A(f"> conda: {cd['note']}")
        A("")
    if jp.get("note"):
        A(f"> jupyter: {jp['note']}")
        A("")
    if jp["kernels"]:
        A("| 커널 | 표시이름 | 파이썬 |")
        A("|---|---|---|")
        for k in jp["kernels"]:
            A(f"| {k['name']} | {k['display_name']} | `{k['python']}` |")
        A("")
    A("> 인계문서 6절 — `subprocess` 는 **커널의 PATH** 를 쓴다. 노트북에서는 위 "
      "`mmseqs 경로` 를 절대경로 상수로 박아 쓸 것.")
    A("")

    A("## 4. 자원 · tmp 잔여물")
    A("")
    A("| 항목 | 값 |")
    A("|---|---|")
    A(f"| nproc / os.cpu_count | {rs['nproc']} / {rs['cpu_count']} |")
    A(f"| {rs['mnt']} | {'존재' if rs['mnt_exists'] else '**없음**'} |")
    d = rs.get("disk") or {}
    A(f"| 디스크 | {('여유 '+str(d.get('free_gib'))+' GiB / 전체 '+str(d.get('total_gib'))+' GiB') if 'free_gib' in d else d.get('error', UNKNOWN)} |")
    A(f"| 재현 tmp 여유 | {rs.get('repro_dir_free_gib')} GiB (`{res['args']['repro_dir']}`) |")
    tmp_cap_note = f" (상한 {CAP['tmp_dirs']} 건까지만 상세)" if rs["tmp_capped"] else ""
    A(f"| tmp* 잔여물 | {rs['n_tmp_found']} 건{tmp_cap_note} |")
    A("")
    if rs["tmp_leftovers"]:
        A("| 경로 | 크기(GiB) | 수정시각 |")
        A("|---|---|---|")
        for t in rs["tmp_leftovers"]:
            A(f"| `{t['path']}` | {t['size_gib']}{'+' if t['size_is_lower_bound'] else ''} | {t['mtime'][:19]} |")
        A("")
        A("> `easy-cluster` 의 tmp 는 수 GB 가 될 수 있고 남아 있으면 다음 실행이 꼬인다(인계문서 4절).")
        A("")

    A("## 5. 219 재현 (인계문서 7절)")
    A("")
    if rep["status"] == "ok":
        A("| 항목 | 값 | 기대 |")
        A("|---|---|---|")
        A(f"| 클러스터 수 | **{rep['n_clusters']}** | {rep['expected_clusters']} |")
        A(f"| 서열 수 | **{rep['n_sequences']}** | {rep['expected_sequences']} |")
        A(f"| _cluster.tsv 행 수 | {rep['n_rows']} | = 서열 수 |")
        A(f"| 입력 FASTA 레코드 | {rep.get('n_fasta_records')} (고유 accession {rep.get('n_fasta_unique_accessions')}) | |")
        A(f"| 스레드 | {rep['threads']} | |")
        A("")
        A(f"실행한 명령: `{rep['command']}`")
        A("")
        c = rep.get("compare") or {}
        if c.get("status") == "ok":
            A("### master_1827.tsv 멤버십 대조")
            A("")
            A("| 항목 | 값 |")
            A("|---|---|")
            A(f"| 조인 행 수 | {c['n_join_rows']} |")
            A(f"| 양쪽에 있는 accession | {c['n_matched']} (조인 성공률 {c['join_rate']:.4f}) |")
            A(f"| master 전용 / 신규 전용 | {c['n_only_master']} / {c['n_only_new']} |")
            A(f"| 결측률 (신규 없음 / master 없음) | {c['missing_rate_new']:.4f} / {c['missing_rate_master']:.4f} |")
            A(f"| 클러스터 수 master / 신규 | {c['n_clusters_master']} / {c['n_clusters_new']} |")
            A(f"| **대표 ID 동일** | {c['n_same_rep_id']} / {c['n_join_rows']} ({c['same_rep_id_rate']:.4f}) |")
            A(f"| **분할 1:1 동일** | {'예 — 실질적으로 같은 분할' if c['partition_identical'] else '아니오'} |")
            A(f"| 쪼개짐 / 합쳐짐 | {c['crosstab']['n_split']} / {c['crosstab']['n_merged']} |")
            A(f"| ARI / NMI | {c['ari']} / {c['nmi']} |")
            A(f"| 중복 제거로 버린 행 | {c['n_rows_dropped_for_join']} |")
            A("")
            A("> **대표 ID 가 달라도 분할이 1:1 이면 동일하다.** 위 두 줄을 혼동하지 말 것 "
              "— mmseqs 가 어느 멤버를 대표로 고르는지는 버전에 따라 달라진다(인계문서 7절).")
            A("")
        elif c:
            A(f"멤버십 대조: {c.get('status')} — {c.get('reason', '')}")
            A("")
    else:
        A(f"**{rep['status']}** — {rep['reason']}")
        A("")
        A(f"재현하려면: `--reproduce` 와 함께 mmseqs 경로(`--mmseqs`), 입력 FASTA"
          f"(`--ipr-fasta {res['args']['ipr_fasta']}`) 가 필요하다.")
        A("")

    A("## 6. search 명령 복원 (인계문서 2절 — 가장 중요)")
    A("")
    st = sr["stats"]
    A("| 항목 | 값 |")
    A("|---|---|")
    A(f"| 상태 | {sr['status']} |")
    A(f"| 노트북 | 발견 {st['n_notebooks_found']} · 파싱 성공 {st['n_notebooks_parsed']}"
      f" · json 오류 {len(st['notebooks_json_error'])} · 건너뜀 {len(st['notebooks_skipped'])}"
      f" · 상한으로 버림 {st['n_notebooks_dropped_by_cap']} |")
    A(f"| 텍스트·히스토리 파일 | 스캔 {st['n_text_scanned']} · 줄 {st['n_lines_scanned']:,}"
      f" · 건너뜀 {len(st['text_skipped'])} · 꼬리만 읽음 {len(st['text_truncated'])}"
      f" · 상한으로 버림 {st['n_text_dropped_by_cap']} |")
    A(f"| 후보 | 원시 {sr['n_candidates_raw']} → 중복 제거 {sr['n_candidates']}"
      f" (높음 {sr['n_high']} · 중간 {sr['n_mid']}) · 상한으로 버림 {sr['n_candidates_dropped_by_cap']} |")
    tiers = sr.get("n_by_tier") or {}
    A(f"| 매칭 경로 | 직접 {tiers.get('직접', 0)} · 별칭 {tiers.get('별칭', 0)}"
      f" · 전용옵션 {tiers.get('전용옵션', 0)} |")
    A(f"| mmseqs 별칭 대입 줄 | {sr.get('n_alias_assignments', 0)} 건 |")
    ac = sr["align_columns"]
    A(f"| 채점 기준 파일 | `{ac['path']}` — {ac['status']}"
      f"{' (' + str(ac['n_cols']) + '열: ' + ', '.join(map(str, (ac['columns'] or [])[:12])) + ')' if ac.get('columns') else ''} |")
    A("")
    A("핵심 노트북(인계문서 9절):")
    A("")
    for k, v in sr["key_notebooks"].items():
        A(f"- `{k}` — {v if isinstance(v, str) else v[0]}")
    A("")
    unread = (list(st.get("notebooks_json_error") or []) + list(st.get("notebooks_skipped") or [])
              + list(st.get("text_skipped") or []) + list(st.get("text_truncated") or []))
    if unread:
        A("### 끝까지 읽지 못한 파일 — **여기에 그 명령이 있을 수 있다**")
        A("")
        for r in (st.get("notebooks_json_error") or [])[:20]:
            A(f"- json 파싱 실패: `{r['path']}` — {r['reason'][:120]}")
        for r in (st.get("notebooks_skipped") or [])[:20]:
            A(f"- 노트북 건너뜀: `{r['path']}` — {r['reason']}")
        for r in (st.get("text_skipped") or [])[:20]:
            A(f"- 텍스트 건너뜀: `{r['path']}` — {r['reason']}")
        for r in (st.get("text_truncated") or [])[:20]:
            A(f"- 꼬리 {r['kept']:,} 줄만 읽음 (전체 {r['n_lines']}): `{r['path']}`")
        A("")
        A("> 상한은 `CAP` 상수에 있다 (파일 "
          f"{CAP['scan_file_mb']} MB · 줄 {CAP['history_lines']:,} · 노트북 {CAP['scan_notebooks']}"
          f" · 텍스트 {CAP['scan_text_files']}). 위 파일은 직접 `grep -n mmseqs` 할 것.")
        A("")
    if sr.get("alias_assignments"):
        A("mmseqs 가 등장하는 대입 줄 — 별칭 매칭의 근거이자 **당시 쓴 절대경로의 단서**:")
        A("")
        for a in sr["alias_assignments"][:10]:
            A(f"- `{a['source']}` 줄 {a['line_no']} (`{a['alias']}`): `{a['text']}`")
        if sr["n_alias_assignments"] > 10:
            A(f"- … 외 {sr['n_alias_assignments'] - 10} 건 (env_probe.json)")
        A("")
    if sr["candidates"]:
        A(f"### 후보 상위 {min(CAP['candidates_md'], len(sr['candidates']))}"
          f" (전체 {sr['n_candidates']} 건은 env_probe.json)")
        A("")
        for i, c in enumerate(sr["candidates"][: CAP["candidates_md"]], 1):
            A(f"**{i}. 신뢰도 {c['confidence']} (점수 {c['score']}, 매칭 {c.get('match_tier')})** — "
              f"`{c['source']}` / {c['locator']} / 줄 {c['line_no']}")
            A("")
            A("```")
            A(c["text"])
            A("```")
            if c.get("format_output"):
                A(f"- `--format-output` ({c.get('format_output_from')}에서): "
                  f"{', '.join(c['format_output'])}")
            for r in c["reasons"]:
                A(f"- {r}")
            if c.get("seq_paths"):
                A(f"- 서열 경로 후보: {', '.join('`'+p+'`' for p in c['seq_paths'])}")
            A("")
            A("<details><summary>문맥</summary>")
            A("")
            A("```")
            A(c["context"])
            A("```")
            A("")
            A("</details>")
            A("")
    else:
        A("후보 없음. 아직 복원되지 않았다. 다음을 확인할 것:")
        A("")
        A("- `--scan-dir` 가 실제 노트북이 있는 디렉터리를 가리키는지 (현재: "
          + ", ".join(f"`{d['path']}`({'존재' if d['exists'] else '없음'})" for d in sr["scan_dirs"]) + ")")
        A("- `--history` 경로 (현재: "
          + ", ".join(f"`{h['path']}`({'존재' if h['exists'] else '없음'})" for h in sr["history_files"]) + ")")
        A("- Jupyter 서버 로그(`~/.local/share/jupyter/runtime`), 사용자 `unsit` 의 홈 디렉터리")
        A("")

    A("## 7. 불일치 · 미확인 목록")
    A("")
    if disc:
        for i, d_ in enumerate(disc, 1):
            A(f"{i}. {d_}")
    else:
        A("없음 — 탐지·재현·복원 모두 가정대로였다.")
    A("")
    A("## 8. 다음 할 일 (인계문서 10절)")
    A("")
    A("1. 위 1절 블록을 `$SPLIT/selection_params.txt` 에 추가한다.")
    A("2. `미확인(...)` 로 남은 칸은 사유를 보고 해당 명령을 직접 돌린다.")
    A("3. 6절 후보 중 신뢰도 '높음' 을 사람이 확인해 `search 명령` 칸을 확정한다.")
    A("4. 219 가 재현되지 않았다면 그 숫자를 새로 기록하고 분할 1:1 여부로 판단한다.")
    A("")
    return "\n".join(L)


# ===========================================================================
# CLI
# ===========================================================================
def default_history_paths() -> list:
    h = Path.home()
    return [h / ".bash_history", h / ".zsh_history", h / ".sh_history", h / ".histfile",
            h / ".python_history", h / ".local/share/fish/fish_history",
            h / ".ipython/profile_default/history.sqlite"]


def _empty_recovery(acols, status, reason, **extra) -> dict:
    """복원을 돌리지 못했을 때의 빈 결과. 렌더러가 기대하는 키를 모두 채운다."""
    return dict(requested=(status != "요청 안 함"), status=status, reason=reason,
                candidates=[], n_candidates=0, n_candidates_raw=0,
                n_candidates_dropped_by_cap=0, n_duplicates_merged=0, n_high=0, n_mid=0,
                n_by_tier={}, alias_assignments=[], n_alias_assignments=0,
                key_notebooks={}, scan_dirs=[], history_files=[], align_columns=acols,
                query_source_guess=None,
                stats=dict(n_notebooks_found=0, n_notebooks_parsed=0, n_text_files_found=0,
                           n_text_scanned=0, n_lines_scanned=0, notebooks_skipped=[],
                           notebooks_json_error=[], text_skipped=[], text_truncated=[],
                           n_notebooks_dropped_by_cap=0, n_text_dropped_by_cap=0),
                **extra)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="단계 0 — mmseqs 환경 탐지 · 219 재현 · search 명령 복원",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--outdir", type=Path, default=common.OUT_ROOT / "00_env_probe",
                   help="env_probe.json / env_probe.md / params.json 을 쓸 곳")
    p.add_argument("--mnt", type=Path, default=common.MNT, help="작업 루트 (디스크·tmp 탐지 기준)")
    p.add_argument("--dock", type=Path, default=common.DOCK, help="노트북 스캔 기본 루트")
    p.add_argument("--ipr-fasta", type=Path, default=common.IPR_FASTA, help="입력 FASTA (1,827 서열)")
    p.add_argument("--master-tsv", type=Path, default=common.MASTER_TSV, help="기존 기준표 (멤버십 대조)")
    p.add_argument("--align-tsv", type=Path, default=common.ALIGN_TSV,
                   help="chcodh2_alignment.tsv — 복원 후보 채점의 기준")
    p.add_argument("--mmseqs", default=os.environ.get("CHCODH2_MMSEQS"),
                   help="mmseqs 절대경로. 생략하면 common.find_mmseqs 로 탐색")
    p.add_argument("--scan-dir", type=Path, action="append", default=None,
                   help="복원 스캔 디렉터리 (반복 지정 가능). 기본 --dock")
    p.add_argument("--history", type=Path, action="append", default=None,
                   help="히스토리 파일 (반복 지정 가능). 기본 ~/.bash_history 등")
    p.add_argument("--include-md", action="store_true",
                   help="*.md 도 스캔 (인계문서 자체가 후보로 잡혀 잡음이 늘어난다)")
    p.add_argument("--reproduce", action="store_true", help="인계문서 7절 easy-cluster 재현 실행")
    p.add_argument("--recover-search", dest="recover_search", action="store_true",
                   help="인계문서 2절 search 명령 복원 스캔")
    p.add_argument("--no-recover-search", dest="recover_search", action="store_false",
                   help="복원 스캔 생략")
    p.set_defaults(recover_search=True)
    p.add_argument("--repro-dir", type=Path, default=None,
                   help="재현 산출물·tmp 위치. 기본 <outdir>/reproduce")
    p.add_argument("--repro-timeout", type=int, default=24 * 3600, help="easy-cluster 제한시간(초)")
    p.add_argument("--keep-tmp", action="store_true", help="재현 후 mmseqs tmp 를 남긴다")
    p.add_argument("--threads", type=int, default=None,
                   help="재현에 --threads 를 추가한다 (기록된 명령과 달라지므로 보고된다)")
    p.add_argument("--expect-clusters", type=int, default=EXPECTED_N_CLUSTERS)
    p.add_argument("--expect-sequences", type=int, default=EXPECTED_N_SEQUENCES)
    p.add_argument("--max-conda-envs", type=int, default=CAP["conda_list_envs"],
                   help="`conda list` 를 돌릴 환경 수 상한")
    p.add_argument("--tmp-since", default="2026-01-01", help="tmp 잔여물 수정시각 하한 (YYYY-MM-DD)")
    p.add_argument("--dry-run", action="store_true",
                   help="외부 명령을 하나도 실행하지 않는다. 데이터 없이 산출물 모양만 확인할 때")
    p.add_argument("--strict", action="store_true",
                   help="불일치가 있으면 산출물을 모두 쓴 뒤 종료코드 2 로 끝낸다")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    args.scan_dir = args.scan_dir or [args.dock]
    args.history = args.history or default_history_paths()
    args.repro_dir = args.repro_dir or (args.outdir / "reproduce")
    outdir = args.outdir
    log_dir = outdir / "logs"
    outdir.mkdir(parents=True, exist_ok=True)

    disc: list = []
    print(f"[00_env_probe] 시작 {now_utc()}  host={os.uname().nodename}  outdir={outdir}")
    if args.dry_run:
        print("[모드] --dry-run — 외부 명령을 실행하지 않는다")

    # --- (a) 환경 탐지
    py = detect_python_env(args.dry_run)
    mm = detect_mmseqs(args, log_dir)
    args.mmseqs_resolved = mm["path"] if mm["found"] else None
    cd = detect_conda(args, log_dir)
    jp = detect_jupyter(args, log_dir)
    rs = detect_resources(args, log_dir)

    print(f"[탐지] mmseqs   {'찾음: ' + mm['path'] if mm['found'] else NOT_FOUND + ' — ' + (mm['reason'] or '')}")
    print(f"[탐지] 버전     {mm['version']}")
    print(f"[탐지] 설치방식 {mm['install']['guess']}")
    print(f"[탐지] conda    {cd['conda_exe']} (환경 {cd['n_envs']}개)")
    print(f"[탐지] jupyter  {jp['jupyter_exe']} (커널 {len(jp['kernels'])}개)")
    print(f"[탐지] python   {py['python_executable']} ({py['python_version']})")
    print(f"[탐지] 자원     nproc={rs['nproc']} · {rs['mnt']} {'존재' if rs['mnt_exists'] else NOT_FOUND}"
          f" · tmp 잔여물 {rs['n_tmp_found']}건")
    if not mm["found"]:
        disc.append(f"mmseqs 바이너리를 찾지 못했다 ({mm['reason'] or 'PATH·탐색경로 모두 적중 없음'})")
    if str(mm["version"]).startswith(UNKNOWN):
        disc.append("mmseqs 버전을 확정하지 못했다 — 219 재현성이 버전에 묶여 있다(인계문서 4절)")
    if cd["conda_exe"] == NOT_FOUND:
        disc.append("conda 를 찾지 못했다 — 커널↔환경 대응이 미확인으로 남는다")
    if jp["jupyter_exe"] == NOT_FOUND:
        disc.append("jupyter 를 찾지 못했다 — kernelspec 미확인(인계문서 6절 함정 점검 불가)")
    if not rs["mnt_exists"]:
        disc.append(f"{rs['mnt']} 가 없다 — 이 장비는 epel-af2 가 아니다. 데이터 기반 항목은 전부 미확인")

    # --- (b) 219 재현
    # 가정이 깨지면 예외를 올려야 하지만(규칙 6), 탐지 결과를 잃으면 안 된다(규칙 5).
    # 그래서 여기서 붙잡아 두고 산출물을 모두 쓴 뒤 마지막에 그대로 다시 올린다.
    stage_exc = None
    if args.reproduce:
        try:
            rep = reproduce_clusters(args, log_dir, disc)
        except Exception as e:  # noqa: BLE001
            stage_exc = e
            rep = dict(requested=True, status="예외", reason=f"{type(e).__name__}: {e}",
                       traceback=traceback.format_exc(),
                       expected_clusters=args.expect_clusters,
                       expected_sequences=args.expect_sequences,
                       n_clusters=None, n_sequences=None, clusters_match=None,
                       sequences_match=None, compare=None, command=None,
                       threads=(args.threads or "미지정(전 코어)"))
            disc.append(f"219 재현 중 예외: {type(e).__name__}: {e}")
            print(f"[재현] 예외 — {type(e).__name__}: {e}")
        if rep["status"] == "ok":
            print(f"[재현] 클러스터 {rep['n_clusters']} / 서열 {rep['n_sequences']} "
                  f"(기대 {args.expect_clusters} / {args.expect_sequences}) "
                  f"→ {'일치' if rep['clusters_match'] and rep['sequences_match'] else '불일치'}")
            c = rep.get("compare") or {}
            if c.get("status") == "ok":
                print(f"[대조] 조인 {c['n_matched']}/{c['n_join_rows']} ({c['join_rate']:.4f}) · "
                      f"대표 ID 동일 {c['n_same_rep_id']} ({c['same_rep_id_rate']:.4f}) · "
                      f"분할 1:1 동일 {'예' if c['partition_identical'] else '아니오'} "
                      f"(ARI={c['ari']}, 쪼개짐 {c['crosstab']['n_split']}, 합쳐짐 {c['crosstab']['n_merged']})")
            else:
                print(f"[대조] {c.get('status')} — {c.get('reason', '')}")
        else:
            print(f"[재현] {rep['status']} — {rep['reason']}")
    else:
        rep = dict(requested=False, status="요청 안 함", reason="--reproduce 미지정",
                   expected_clusters=args.expect_clusters, expected_sequences=args.expect_sequences,
                   n_clusters=None, n_sequences=None, clusters_match=None, sequences_match=None,
                   compare=None, threads=(args.threads or "미지정(전 코어)"), command=None)
        print("[재현] 건너뜀 (--reproduce 미지정)")

    # --- (c) search 명령 복원
    acols = align_columns(args)
    if not args.recover_search:
        sr = _empty_recovery(acols, "요청 안 함", "--no-recover-search")
        print("[복원] 건너뜀 (--no-recover-search)")
    elif stage_exc is not None:
        sr = _empty_recovery(acols, "건너뜀", "앞 단계에서 예외가 발생해 중단했다")
        print("[복원] 건너뜀 — 앞 단계 예외로 중단 (산출물은 기록한다)")
    else:
        try:
            sr = recover_search(args, acols, disc)
        except Exception as e:  # noqa: BLE001
            stage_exc = e
            sr = _empty_recovery(acols, "예외", f"{type(e).__name__}: {e}",
                                 traceback=traceback.format_exc())
            disc.append(f"search 복원 중 예외: {type(e).__name__}: {e}")
            print(f"[복원] 예외 — {type(e).__name__}: {e}")
        if sr["status"] != "예외":
            st = sr["stats"]
            print(f"[복원] 노트북 {st['n_notebooks_found']}개 발견 / {st['n_notebooks_parsed']}개 파싱, "
                  f"텍스트·히스토리 {st['n_text_scanned']}개 / {st['n_lines_scanned']:,}줄 스캔")
            print(f"[복원] 후보 {sr['n_candidates']}건 (원시 {sr['n_candidates_raw']}, "
                  f"높음 {sr['n_high']}, 중간 {sr['n_mid']}, "
                  f"상한으로 버림 {sr['n_candidates_dropped_by_cap']})")
            tiers = sr.get("n_by_tier") or {}
            print(f"[복원] 매칭 경로 직접 {tiers.get('직접', 0)} · 별칭 {tiers.get('별칭', 0)} · "
                  f"전용옵션 {tiers.get('전용옵션', 0)} · "
                  f"별칭 대입 줄 {sr.get('n_alias_assignments', 0)}건")
            for c in sr["candidates"][:3]:
                flat = RE_WS.sub(" ", c["text"])[:150]
                print(f"       [{c['confidence']}/{c['score']}/{c.get('match_tier')}] "
                      f"{Path(c['source']).name} {c['locator']}: {flat}")

    # --- (d) 산출물
    res = dict(
        stage="00_env_probe",
        generated_utc=now_utc(),
        handoff="20261021_handoff_mmseqs_env.md (5·7·8절 + 2절)",
        args={k: (str(v) if isinstance(v, Path) else
                  [str(x) for x in v] if isinstance(v, list) else v)
              for k, v in vars(args).items()},
        caps=CAP,
        python=py, mmseqs=mm, conda=cd, jupyter=jp, resources=rs,
        reproduce=rep, search_recovery=sr,
        discrepancies=disc,
    )

    counts = dict(
        n_discrepancies=len(disc),
        mmseqs_found=bool(mm["found"]),
        mmseqs_version=mm["version"],
        n_conda_envs=cd["n_envs"],
        n_conda_envs_listed=cd["n_envs_conda_listed"],
        conda_list_capped=bool(cd.get("conda_list_capped")),
        n_jupyter_kernels=len(jp["kernels"]),
        n_tmp_leftovers_found=rs["n_tmp_found"],
        n_tmp_leftovers_detailed=len(rs["tmp_leftovers"]),
        reproduce_status=rep["status"],
        reproduce_n_clusters=rep.get("n_clusters"),
        reproduce_n_sequences=rep.get("n_sequences"),
        reproduce_clusters_match=rep.get("clusters_match"),
        search_status=sr.get("status"),
        n_notebooks_found=sr["stats"]["n_notebooks_found"],
        n_notebooks_parsed=sr["stats"]["n_notebooks_parsed"],
        n_notebooks_dropped_by_cap=sr["stats"].get("n_notebooks_dropped_by_cap", 0),
        n_text_files_scanned=sr["stats"]["n_text_scanned"],
        n_text_dropped_by_cap=sr["stats"].get("n_text_dropped_by_cap", 0),
        n_lines_scanned=sr["stats"]["n_lines_scanned"],
        n_candidates_raw=sr.get("n_candidates_raw", 0),
        n_candidates=sr.get("n_candidates", 0),
        n_candidates_high=sr.get("n_high", 0),
        n_candidates_dropped_by_cap=sr.get("n_candidates_dropped_by_cap", 0),
    )
    cmpres = (rep.get("compare") or {})
    if cmpres.get("status") == "ok":
        counts.update(
            join_total=cmpres["n_join_rows"], join_matched=cmpres["n_matched"],
            join_rate=cmpres["join_rate"],
            missing_rate_new=cmpres["missing_rate_new"],
            missing_rate_master=cmpres["missing_rate_master"],
            rep_ids_identical=cmpres["rep_ids_identical"],
            partition_identical=cmpres["partition_identical"],
            ari=cmpres["ari"], nmi=cmpres["nmi"],
            n_rows_dropped_for_join=cmpres["n_rows_dropped_for_join"],
        )
    res["counts"] = counts

    json_path = outdir / "env_probe.json"
    json_path.write_text(json.dumps(res, indent=2, ensure_ascii=False, default=str))
    md_path = outdir / "env_probe.md"
    md_path.write_text(render_md(res, disc))

    commands = [mm.get("version_cmd", {}).get("cmd_str"), (cd.get("env_list") or {}).get("cmd_str"),
                (jp.get("kernelspec") or {}).get("cmd_str"), (rep.get("run") or {}).get("cmd_str")]
    common.write_params(
        outdir, "00_env_probe",
        inputs=[args.ipr_fasta, args.master_tsv, args.align_tsv],
        outputs=[json_path, md_path],
        commands=[c for c in commands if c],
        counts=counts,
        discrepancies=disc,
        caps=CAP,
        mmseqs_path=mm["path"], mmseqs_version=mm["version"],
        install_method=mm["install"]["guess"],
        scan_dirs=[str(d) for d in args.scan_dir],
        history_files=[str(h) for h in args.history],
        dry_run=args.dry_run, strict=args.strict,
        reproduce_requested=bool(args.reproduce),
        recover_search_requested=bool(args.recover_search),
    )

    print(f"[쓰기] {json_path}")
    print(f"[쓰기] {md_path}")
    print(f"[쓰기] {outdir / 'params.json'}")
    print(f"[결과] 불일치·미확인 {len(disc)} 건")
    for i, d_ in enumerate(disc, 1):
        print(f"       {i}. {d_}")

    # 가정이 깨진 예외는 산출물을 남긴 뒤 그대로 올린다 (규칙 5 + 규칙 6)
    if stage_exc is not None:
        print("[예외] 가정이 깨졌다. 탐지 결과는 위 경로에 모두 기록했고, 아래 예외로 끝낸다.")
        raise stage_exc

    if disc and args.strict:
        print("[strict] 불일치가 있어 종료코드 2 로 끝낸다 (산출물은 모두 기록되었다)")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
