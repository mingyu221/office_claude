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

    못 찾아도 예외를 던지지 않는다. 노트북에서는 None 을 받아 셀 안에서 직접 경로를
    지정하는 편이 낫다. 찾은 경우 버전까지 같이 출력한다.
    """
    cands = []
    if explicit:
        cands.append(Path(explicit))
    env = os.environ.get("CHCODH2_MMSEQS")
    if env:
        cands.append(Path(env))
    which = shutil.which("mmseqs")
    if which:
        cands.append(Path(which))
    for root in (Path("/mnt/af2results/mingyu/soft"), Path.home(), Path("/opt"), Path("/usr/local")):
        if not root.is_dir():
            continue
        try:
            cands.extend(sorted(root.rglob("mmseqs"))[:3])
        except (PermissionError, OSError):
            pass
    for c in cands:
        try:
            if c.is_file() and os.access(c, os.X_OK):
                p = str(c.resolve())
                if verbose:
                    print(f"mmseqs : {p}")
                    print(f"버전   : {_version(p)}")
                return p
        except OSError:
            continue
    if verbose:
        print("mmseqs 를 못 찾았다.")
        print("  로그인 셸에서  readlink -f \"$(which mmseqs)\"  를 돌려 절대경로를 확인하고")
        print("  이 셀에  MMSEQS = '/절대/경로/mmseqs'  로 직접 지정할 것.")
        print(f"  (참고) 이 커널의 PATH = {os.environ.get('PATH', '')[:200]}")
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
