"""ChCODH2 상동체 파이프라인 공용 모듈.

설계 원칙
  1. 의존성은 pandas + numpy + 표준 라이브러리뿐. (epel-af2 conda base 에서 설치 없이 실행)
  2. 서열을 지우지 않는다. QC 는 플래그 컬럼으로만 표현한다.
  3. 모든 위치는 ChCODH2 서열 번호(Q9F8A8, 636 aa, Met1 = 1)로 적는다.
  4. 모든 단계는 산출물 옆에 params.json(명령·버전·입력 해시·행 수)을 남긴다.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------
# 경로 — 환경변수로 덮어쓸 수 있다. 기본값은 인계문서(2026-10-21)의 epel-af2 경로.
# --------------------------------------------------------------------------
MNT = Path(os.environ.get("CHCODH2_MNT", "/mnt/af2results/mingyu"))
DOCK = Path(os.environ.get("CHCODH2_DOCK", MNT / "workspace/docking_simulation"))
SPLIT = Path(os.environ.get("CHCODH2_SPLIT", DOCK / "result/homolog_split"))
OLDR = Path(os.environ.get("CHCODH2_OLDR", DOCK / "result/homolog_selection"))

IPR_FASTA = Path(os.environ.get("CHCODH2_IPR_FASTA", MNT / "chcodh2/protein-matching-IPR010047.fasta"))
MASTER_TSV = Path(os.environ.get("CHCODH2_MASTER", SPLIT / "master_1827.tsv"))
ALIGN_TSV = Path(os.environ.get("CHCODH2_ALIGN", OLDR / "chcodh2_alignment.tsv"))

OUT_ROOT = Path(os.environ.get("CHCODH2_OUT", SPLIT / "cluster_sweep"))
REPO_LIB = Path(__file__).resolve().parent

# mmseqs 바이너리. 인계문서 6절(Jupyter PATH 함정) 때문에 탐색 후 검증한다.
MMSEQS_SEARCH_DIRS = [MNT / "soft", Path.home(), Path("/opt"), Path("/usr/local")]

# --------------------------------------------------------------------------
# ChCODH2 참조
# --------------------------------------------------------------------------
CHCODH2_ACC = "Q9F8A8"
CHCODH2_LEN = 636


def chcodh2_seq() -> str:
    """참조 서열(636 aa). lib/chcodh2.fasta 가 단일 출처다."""
    _, seq = next(iter(read_fasta(REPO_LIB / "chcodh2.fasta")))
    if len(seq) != CHCODH2_LEN:
        raise ValueError(f"참조 서열 길이가 {len(seq)} — {CHCODH2_LEN} 이어야 한다")
    return seq


# 논문에 보고된 위치. 규칙 검증용 정답지이며 설계 판정의 입력이 아니다.
#   ref = 1 Kim 2022 Nat Catal / 2 Kim 2024 Nat Commun / 3 Kim 2025 (FeS fortification)
KNOWN_POSITIONS = {
    41:  dict(aa="F", role="d_cluster_aromatic", note="F41A 활성 소실; 전자전달 관문", ref=2),
    43:  dict(aa="E", role="mediator_electrostatic", note="viologen 정전기 결합; 터널#1 입구", ref=2),
    57:  dict(aa="R", role="mediator_affinity", note="R57G/N59L 에서 EV 친화도 10배 개선", ref=2),
    59:  dict(aa="N", role="mediator_affinity", note="R57G/N59L 에서 EV 친화도 10배 개선", ref=2),
    206: dict(aa="Q", role="tunnel3_entrance", note="비선택 터널 #3 입구", ref=3),
    450: dict(aa="K", role="tunnel1_entrance", note="비선택 터널 #1 입구", ref=3),
    559: dict(aa="A", role="tunnel_bottleneck", note="A559W/H 로 O2 내성 61~148배; triad", ref=1),
    565: dict(aa="V", role="tunnel_bottleneck_triad", note="A559/V565/I580 triad", ref=1),
    580: dict(aa="I", role="tunnel_bottleneck_triad", note="A559/V565/I580 triad", ref=1),
    583: dict(aa="L", role="tunnel1_entrance", note="비선택 터널 #1 입구", ref=3),
    586: dict(aa="I", role="tunnel2_entrance", note="비선택 터널 #2 입구", ref=3),
    593: dict(aa="T", role="tunnel2_entrance", note="비선택 터널 #2 입구", ref=3),
    597: dict(aa="T", role="tunnel2_entrance", note="비선택 터널 #2 입구", ref=3),
    599: dict(aa="S", role="tunnel3_entrance", note="비선택 터널 #3 입구", ref=3),
    603: dict(aa="I", role="tunnel3_entrance", note="비선택 터널 #3 입구", ref=3),
    610: dict(aa="V", role="tunnel2_bottleneck", note="V610H 로 O2 내성 추가 개선", ref=3),
}


def verify_known_positions(seq: str | None = None) -> None:
    """참조 서열이 논문 번호와 맞는지 확인한다. 어긋나면 즉시 중단."""
    seq = seq or chcodh2_seq()
    bad = [(p, m["aa"], seq[p - 1]) for p, m in KNOWN_POSITIONS.items() if seq[p - 1] != m["aa"]]
    if bad:
        raise ValueError(f"번호 체계 불일치: {bad}")


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------
def read_fasta(path):
    """(header, sequence) 를 순서대로 돌려준다. 메모리에 전부 올리지 않는다."""
    header, chunks = None, []
    with open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:], []
            elif line:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def parse_accession(header: str) -> str:
    """FASTA 헤더에서 accession 을 뽑는다.

    InterPro 다운로드는 'A0A0A1|...' 또는 'tr|A0A0A1|NAME_ORG' 꼴이 섞인다.
    master_1827.tsv 와 조인하려면 양쪽에서 같은 규칙을 써야 하므로 여기 한 곳에만 둔다.
    """
    tok = header.split()[0]
    parts = tok.split("|")
    if len(parts) >= 3 and parts[0] in ("sp", "tr", "gi", "ref"):
        return parts[1].strip()
    return parts[0].strip()


def norm_acc(s):
    """Series/스칼라의 accession 을 parse_accession 과 같은 규칙으로 정규화."""
    if isinstance(s, pd.Series):
        return s.astype(str).map(lambda x: parse_accession(x))
    return parse_accession(str(s))


# --------------------------------------------------------------------------
# 컬럼 자동 탐지 — 인계문서가 'accession 또는 target' 처럼 불확실하게 적은 부분 대응
# --------------------------------------------------------------------------
def detect_column(df: pd.DataFrame, candidates, *, required=True, what=""):
    """후보 이름 중 실제로 존재하는 첫 컬럼명을 돌려준다(대소문자 무시)."""
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise KeyError(f"{what or candidates[0]} 컬럼을 못 찾음. 후보={candidates} 실제={list(df.columns)}")
    return None


COLSPEC = {
    "accession": ["accession", "target", "tid", "subject", "sseqid", "acc"],
    "query":     ["query", "qid", "qseqid"],
    "fident":    ["fident", "pident", "identity", "seq_id"],
    "qstart":    ["qstart", "q_start"],
    "qaln":      ["qaln", "qseq", "q_aln"],
    "tstart":    ["tstart", "t_start", "sstart"],
    "taln":      ["taln", "tseq", "t_aln"],
    "qcov":      ["qcov", "qcovs", "query_cov"],
    "tcov":      ["tcov", "tcovs", "target_cov"],
    "cluster":   ["cluster", "cluster_rep", "rep", "representative"],
    "evalue":    ["evalue", "e_value", "eval"],
    "bits":      ["bits", "bitscore", "score"],
}


def resolve_columns(df: pd.DataFrame, keys, *, required=()) -> dict:
    """COLSPEC 의 논리 이름 -> 실제 컬럼명. 없으면 None (required 면 예외)."""
    out = {}
    for k in keys:
        out[k] = detect_column(df, COLSPEC[k], required=(k in required), what=k)
    return out


def read_tsv(path, **kw) -> pd.DataFrame:
    """헤더가 있든 없든 읽는다. mmseqs convertalis 산출물은 헤더가 없을 수 있다."""
    df = pd.read_csv(path, sep="\t", dtype=str, **kw)
    # 첫 행이 데이터처럼 보이면 헤더 없는 파일이다 (컬럼명에 '|' 또는 숫자만)
    cols = list(df.columns)
    if any(("|" in c) or re.fullmatch(r"[\d.eE+-]+", c) for c in cols):
        df = pd.read_csv(path, sep="\t", dtype=str, header=None, **kw)
        df.columns = [f"col{i}" for i in range(df.shape[1])]
    return df


# --------------------------------------------------------------------------
# 재현성 기록
# --------------------------------------------------------------------------
def sha256_file(path, limit_mb=None) -> str:
    h = hashlib.sha256()
    n = 0
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 20), b""):
            h.update(blk)
            n += len(blk)
            if limit_mb and n > limit_mb * (1 << 20):
                break
    return h.hexdigest()


def find_mmseqs() -> str:
    """mmseqs 절대경로. 인계문서 6절 — PATH 를 믿지 않고 탐색 후 검증한다."""
    env = os.environ.get("CHCODH2_MMSEQS")
    if env and Path(env).is_file() and os.access(env, os.X_OK):
        return env
    which = shutil.which("mmseqs")
    if which:
        return str(Path(which).resolve())
    for root in MMSEQS_SEARCH_DIRS:
        if not root.is_dir():
            continue
        try:
            for cand in sorted(root.rglob("mmseqs")):
                if cand.is_file() and os.access(cand, os.X_OK):
                    return str(cand.resolve())
        except (PermissionError, OSError):
            continue
    raise FileNotFoundError(
        "mmseqs 를 못 찾음. CHCODH2_MMSEQS 환경변수에 절대경로를 지정할 것 (인계문서 6절)."
    )


def mmseqs_version(binary=None) -> str:
    binary = binary or find_mmseqs()
    try:
        r = subprocess.run([binary, "version"], capture_output=True, text=True, timeout=60)
        return (r.stdout or r.stderr).strip().splitlines()[0] if (r.stdout or r.stderr) else "unknown"
    except Exception as e:  # noqa: BLE001
        return f"unknown ({e})"


def run_cmd(cmd, *, log_path=None, check=True, timeout=None) -> dict:
    """명령을 돌리고 재현에 필요한 모든 것을 기록해 돌려준다."""
    cmd = [str(c) for c in cmd]
    t0 = datetime.now(timezone.utc)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    rec = dict(
        cmd=cmd,
        cmd_str=" ".join(cmd),
        returncode=r.returncode,
        started_utc=t0.isoformat(),
        seconds=round((datetime.now(timezone.utc) - t0).total_seconds(), 2),
        stdout_tail=(r.stdout or "")[-4000:],
        stderr_tail=(r.stderr or "")[-4000:],
    )
    if log_path:
        Path(log_path).parent.mkdir(parents=True, exist_ok=True)
        Path(log_path).write_text(
            f"$ {rec['cmd_str']}\nrc={r.returncode}\n\n--- stdout ---\n{r.stdout}\n--- stderr ---\n{r.stderr}\n"
        )
    if check and r.returncode != 0:
        raise RuntimeError(f"명령 실패 (rc={r.returncode}): {rec['cmd_str']}\n{rec['stderr_tail']}")
    return rec


def write_params(outdir, stage: str, *, inputs=None, outputs=None, commands=None, **extra) -> Path:
    """산출물 옆에 재현 정보를 남긴다. 인계문서가 존재하는 이유가 이것이다."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    def _file_rec(p):
        p = Path(p)
        if not p.exists():
            return dict(path=str(p), exists=False)
        rec = dict(path=str(p), exists=True, bytes=p.stat().st_size)
        if p.is_file() and p.stat().st_size < 2 * (1 << 30):
            rec["sha256"] = sha256_file(p)
            if p.suffix in (".tsv", ".csv", ".txt", ".fasta", ".fa"):
                with open(p, "rt", errors="replace") as fh:
                    rec["lines"] = sum(1 for _ in fh)
        return rec

    payload = dict(
        stage=stage,
        written_utc=datetime.now(timezone.utc).isoformat(),
        host=os.uname().nodename,
        user=os.environ.get("USER") or os.environ.get("LOGNAME"),
        cwd=os.getcwd(),
        python=sys.version.split()[0],
        python_executable=sys.executable,
        pandas=pd.__version__,
        numpy=np.__version__,
        conda_prefix=os.environ.get("CONDA_PREFIX"),
        inputs=[_file_rec(p) for p in (inputs or [])],
        outputs=[_file_rec(p) for p in (outputs or [])],
        commands=commands or [],
        **extra,
    )
    dest = outdir / "params.json"
    dest.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str))
    return dest


# --------------------------------------------------------------------------
# 분할 비교 지표 — sklearn 없이 numpy 로 구현 (의존성 최소화)
# --------------------------------------------------------------------------
def contingency(labels_a, labels_b):
    a = pd.Categorical(pd.Series(labels_a).astype(str))
    b = pd.Categorical(pd.Series(labels_b).astype(str))
    n = np.zeros((len(a.categories), len(b.categories)), dtype=np.int64)
    np.add.at(n, (a.codes, b.codes), 1)
    return n


def adjusted_rand_index(labels_a, labels_b) -> float:
    """ARI. 두 분할이 같으면 1, 무작위면 0 근처.

    퇴화 분할 주의 — 양쪽 모두 한 덩어리이거나 양쪽 모두 싱글턴이면 두 분할은 **완전히 같은데도**
    공식이 0/0 이 된다. sklearn 과 같은 규약으로 1.0 을 돌려준다. 이것을 NaN 으로 두면
    스윕에서 "분할이 전혀 바뀌지 않은 구간"이 "경계(급격히 재편됨)"로 뒤집혀 보고된다.
    """
    n = contingency(labels_a, labels_b)
    total = n.sum()
    if total < 2:
        return float("nan")
    n_a, n_b = n.shape
    if (n_a == n_b == 1) or (n_a == n_b == total):
        return 1.0
    comb2 = lambda x: x * (x - 1) / 2.0  # noqa: E731
    sum_ij = comb2(n).sum()
    sum_i = comb2(n.sum(axis=1)).sum()
    sum_j = comb2(n.sum(axis=0)).sum()
    exp = sum_i * sum_j / comb2(total)
    mx = (sum_i + sum_j) / 2.0
    return float("nan") if mx == exp else float((sum_ij - exp) / (mx - exp))


def normalized_mutual_info(labels_a, labels_b) -> float:
    """NMI (산술평균 정규화). 0~1."""
    n = contingency(labels_a, labels_b).astype(float)
    total = n.sum()
    if total == 0:
        return float("nan")
    p = n / total
    pi, pj = p.sum(axis=1), p.sum(axis=0)
    nz = p > 0
    mi = (p[nz] * np.log(p[nz] / np.outer(pi, pj)[nz])).sum()
    ent = lambda v: -(v[v > 0] * np.log(v[v > 0])).sum()  # noqa: E731
    hi, hj = ent(pi), ent(pj)
    denom = (hi + hj) / 2.0
    return 1.0 if denom == 0 else float(mi / denom)


def partition_summary(labels) -> dict:
    """한 분할의 모양. 스윕에서 임계값별로 비교하는 값들."""
    s = pd.Series(labels).astype(str)
    sizes = s.value_counts()
    total = int(len(s))
    return dict(
        n_sequences=total,
        n_clusters=int(sizes.size),
        n_singletons=int((sizes == 1).sum()),
        singleton_frac=float((sizes == 1).sum() / sizes.size) if sizes.size else float("nan"),
        largest_cluster=int(sizes.iloc[0]) if sizes.size else 0,
        largest_frac=float(sizes.iloc[0] / total) if total else float("nan"),
        median_size=float(sizes.median()) if sizes.size else float("nan"),
        mean_size=float(sizes.mean()) if sizes.size else float("nan"),
    )


def crosstab_mapping(labels_old, labels_new) -> dict:
    """두 분할이 1:1 로 대응하는지. 대표 ID 가 바뀌어도 분할이 같은지 보는 용도."""
    ct = contingency(labels_old, labels_new)
    rows_nonzero = (ct > 0).sum(axis=1)
    cols_nonzero = (ct > 0).sum(axis=0)
    return dict(
        n_old=int(ct.shape[0]),
        n_new=int(ct.shape[1]),
        old_maps_to_one=int((rows_nonzero == 1).sum()),
        new_maps_to_one=int((cols_nonzero == 1).sum()),
        one_to_one=int(((rows_nonzero == 1)[:, None] & (cols_nonzero == 1)[None, :] & (ct > 0)).sum()),
        n_split=int((rows_nonzero > 1).sum()),
        n_merged=int((cols_nonzero > 1).sum()),
    )
