"""노트북에서 단계 스크립트를 불러 쓰기 위한 얇은 어댑터.

왜 필요한가
  단계 스크립트 이름이 숫자로 시작해(02_cluster_sweep.py) `import 02_cluster_sweep` 이 안 된다.
  여기서 importlib 로 감싸 노트북 첫 셀을 한 줄로 만든다. 로직은 전부 스크립트 쪽에 있고
  이 파일은 경로 해결과 Jupyter 함정 처리만 한다.

노트북 첫 셀
    import sys
    from pathlib import Path
    REPO = Path.home() / "office_claude"          # 여기만 고치면 된다
    sys.path.insert(0, str(REPO / "pipeline" / "lib"))
    import nb, common
    sweep = nb.load_stage("02_cluster_sweep")
"""
from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

LIB = Path(__file__).resolve().parent
PIPELINE = LIB.parent
REPO = PIPELINE.parent
STAGES = PIPELINE / "stages"


def load_stage(name: str):
    """단계 스크립트를 모듈로 불러온다. name 예: '02_cluster_sweep' (.py 생략 가능)."""
    name = name[:-3] if name.endswith(".py") else name
    path = STAGES / f"{name}.py"
    if not path.is_file():
        avail = sorted(p.stem for p in STAGES.glob("*.py"))
        raise FileNotFoundError(f"{path} 가 없다. 사용 가능: {avail}")
    mod_name = f"chcodh2_stage_{name}"
    if mod_name in sys.modules:
        return sys.modules[mod_name]
    if str(LIB) not in sys.path:
        sys.path.insert(0, str(LIB))  # 스크립트가 common 을 import 한다
    spec = importlib.util.spec_from_file_location(mod_name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


def find_mmseqs(explicit: str | None = None, *, verbose: bool = True) -> str | None:
    """mmseqs 절대경로. 인계문서 6절 — Jupyter 커널의 PATH 는 로그인 셸과 다르다.

    mmseqs 는 파이썬 패키지가 아니라 외부 실행파일이라 전용 커널이 필요 없다.
    어느 커널에서든 절대경로로 부르면 된다. 못 찾아도 예외를 던지지 않는다.
    """
    seen, cands = set(), []

    def add(path):
        if path and str(path) not in seen:
            seen.add(str(path))
            cands.append(Path(path))

    add(explicit)
    add(os.environ.get("CHCODH2_MMSEQS"))
    add(shutil.which("mmseqs"))

    # conda 환경 — 현재 커널의 환경, 그리고 형제 환경들
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        add(Path(prefix) / "bin/mmseqs")
        for env in sorted(Path(prefix).parent.glob("*/bin/mmseqs")):
            add(env)
    for root in (Path.home() / "anaconda3", Path.home() / "miniconda3",
                 Path.home() / "miniforge3", Path.home() / "mambaforge",
                 Path("/opt/conda")):
        add(root / "bin/mmseqs")
        if (root / "envs").is_dir():
            for env in sorted((root / "envs").glob("*/bin/mmseqs")):
                add(env)

    # 소프트웨어 트리와 표준 위치 (홈 전체 rglob 은 느려서 쓰지 않는다)
    for root in (Path("/mnt/af2results/mingyu/soft"), Path("/usr/local"), Path("/opt")):
        if root.is_dir():
            try:
                for d in ("mmseqs", "*/mmseqs", "*/bin/mmseqs", "*/*/bin/mmseqs"):
                    for hit in sorted(root.glob(d))[:5]:
                        add(hit)
            except (PermissionError, OSError):
                pass

    for c in cands:
        try:
            if c.is_file() and os.access(c, os.X_OK):
                path = str(c.resolve())
                if verbose:
                    print(f"mmseqs : {path}")
                    print(f"버전   : {_version(path)}")
                return path
        except OSError:
            continue

    if verbose:
        print("mmseqs 를 못 찾았다. 커널을 바꿀 필요는 없다 — 절대경로만 있으면 된다.")
        print('  터미널에서  readlink -f "$(which mmseqs)"  를 돌리고')
        print("  그 경로를 1절의 MMSEQS_PATH 에 적은 뒤 1절부터 다시 실행할 것.")
        print(f"  (이 커널의 CONDA_PREFIX = {os.environ.get('CONDA_PREFIX')})")
    return None


def _version(binary: str) -> str:
    try:
        r = subprocess.run([binary, "version"], capture_output=True, text=True, timeout=60)
        out = (r.stdout or r.stderr).strip()
        return out.splitlines()[0] if out else "unknown"
    except Exception as e:  # noqa: BLE001
        return f"unknown ({e})"


def env_report() -> None:
    """노트북 커널이 어디서 무엇으로 도는지. 첫 셀에서 한 번 찍어두면 나중에 헷갈리지 않는다."""
    import numpy as np
    import pandas as pd
    print(f"python       {sys.version.split()[0]}  ({sys.executable})")
    print(f"pandas       {pd.__version__}")
    print(f"numpy        {np.__version__}")
    print(f"CONDA_PREFIX {os.environ.get('CONDA_PREFIX')}")
    print(f"cwd          {os.getcwd()}")
    print(f"REPO         {REPO}")
    print(f"nproc        {os.cpu_count()}")
