"""ChCODH2 상동체 클러스터 스윕 — 노트북 전용 단일 모듈.

cluster_sweep.ipynb 옆에 두기만 하면 된다. git 도, 다른 파일도 필요 없다.
의존성은 pandas + numpy + 표준 라이브러리뿐.

이 파일은 pipeline/lib/common.py 와 pipeline/stages/02_cluster_sweep.py 에서
노트북이 쓰는 함수만 뽑아 합친 것이다. 원본은 git 레포에 있고 거기서 검증됐다.
직접 고치지 말 것 — 다음 갱신 때 덮어쓴다.
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
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "2026-10-25b"   # 노트북 1절이 이 값을 대조한다. 어긋나면 둘 중 하나가 낡은 것이다.

CHCODH2_ACC = "Q9F8A8"
CHCODH2_LEN = 636


def env_report() -> None:
    """커널이 어디서 무엇으로 도는지. 1절에서 한 번 찍어두면 나중에 헷갈리지 않는다."""
    print(f"python       {sys.version.split()[0]}  ({sys.executable})")
    print(f"pandas       {pd.__version__}")
    print(f"numpy        {np.__version__}")
    print(f"CONDA_PREFIX {os.environ.get('CONDA_PREFIX')}")
    print(f"모듈         {__file__}  (v{VERSION})")
    print(f"nproc        {os.cpu_count()}")


def find_mmseqs(explicit: str | None = None, *, verbose: bool = True) -> str | None:
    """mmseqs 절대경로. mmseqs 는 파이썬 패키지가 아니라 외부 실행파일이라
    전용 커널이 필요 없다. 못 찾아도 예외를 던지지 않는다."""
    seen, cands = set(), []

    def add(path):
        if path and str(path) not in seen:
            seen.add(str(path)); cands.append(Path(path))

    add(explicit)
    add(os.environ.get("CHCODH2_MMSEQS"))
    add(shutil.which("mmseqs"))
    prefix = os.environ.get("CONDA_PREFIX")
    if prefix:
        add(Path(prefix) / "bin/mmseqs")
        for env in sorted(Path(prefix).parent.glob("*/bin/mmseqs")): add(env)
    for root in (Path.home() / "mmseqs", Path.home() / "anaconda3",
                 Path.home() / "miniconda3", Path.home() / "miniforge3", Path("/opt/conda")):
        add(root / "bin/mmseqs")
        if (root / "envs").is_dir():
            for env in sorted((root / "envs").glob("*/bin/mmseqs")): add(env)
    for root in (Path("/mnt/af2results/mingyu/soft"), Path("/usr/local"), Path("/opt")):
        if root.is_dir():
            try:
                for pat in ("mmseqs", "*/mmseqs", "*/bin/mmseqs", "*/*/bin/mmseqs"):
                    for hit in sorted(root.glob(pat))[:5]: add(hit)
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
        print('  터미널에서  readlink -f "$(which mmseqs)"  를 돌려')
        print("  그 경로를 1절의 MMSEQS_PATH 에 적고 1절부터 다시 실행할 것.")
    return None


def _version(binary: str) -> str:
    try:
        r = subprocess.run([binary, "version"], capture_output=True, text=True, timeout=60)
        out = (r.stdout or r.stderr).strip()
        return out.splitlines()[0] if out else "unknown"
    except Exception as e:  # noqa: BLE001
        return f"unknown ({e})"


def dir_bytes(path) -> int:
    """디렉터리 총 바이트. tmp 용량 측정용."""
    total = 0
    for p in Path(path).rglob("*"):
        try:
            if p.is_file(): total += p.stat().st_size
        except OSError:
            pass
    return total


DEFAULT_THRESHOLDS = [round(p / 100.0, 2) for p in range(10, 91, 5)]

FIXED_COV = "0.8"

FIXED_COV_MODE = "1"

FIXED_SENSITIVITY = "7.5"

ADJACENT_COLUMNS = [
    "label_low", "label_high", "min_seq_id_low", "min_seq_id_high", "n_common",
    "ari", "nmi", "n_clusters_coarser", "n_clusters_finer", "delta_n_clusters",
    "finer_maps_to_one", "coarser_maps_to_one", "one_to_one", "n_finer_split",
    "n_coarser_merged", "is_stable", "ari_degenerate", "ref_jaccard",
    "ref_n_added", "ref_n_removed",
    "ref_n_source_clusters", "ref_merged_from", "ref_merged_truncated",
]

def cluster_members(ser: pd.Series, acc: str) -> list[str] | None:
    """acc 가 속한 클러스터의 멤버 목록. 없으면 None."""
    if acc not in ser.index:
        return None
    rep = ser.loc[acc]
    if isinstance(rep, pd.Series):          # 중복 index 방어 (정상 경로에서는 발생하지 않는다)
        rep = rep.iloc[0]
    return sorted(ser.index[ser.to_numpy() == rep].tolist())

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

def compare_partitions(a: pd.Series, b: pd.Series) -> dict:
    """공통 accession 위에서만 비교한다. 공통 크기와 피복률을 항상 같이 돌려준다."""
    common_idx = a.index.intersection(b.index)
    n_common = int(len(common_idx))
    out = dict(n_common=n_common,
               common_frac_a=float(n_common / len(a)) if len(a) else float("nan"),
               common_frac_b=float(n_common / len(b)) if len(b) else float("nan"))
    if n_common < 2:
        out.update(ari=float("nan"), nmi=float("nan"), ari_degenerate=False)
        return out
    la = a.reindex(common_idx).astype(str).tolist()
    lb = b.reindex(common_idx).astype(str).tolist()
    ari = adjusted_rand_index(la, lb)

    # ARI 의 퇴화 사례 보정. adjusted_rand_index 는 기대값과 최대값이 같아지면 NaN 을
    # 돌려주는데, 그런 경우가 정확히 둘이고 **둘 다 두 분할이 완전히 동일한 상황**이다.
    #   ① 양쪽이 전부 한 클러스터      (낮은 min-seq-id 에서 전부 한 묶음으로 뭉친다)
    #   ② 양쪽이 전부 싱글턴           (높은 min-seq-id 에서 아무것도 묶이지 않는다)
    # 동일한 분할의 ARI 는 1.0 이 맞다 (sklearn.adjusted_rand_score 도 이 두 경우를 1.0 으로
    # 특수 처리한다). 보정하지 않으면 is_stable 이 False 가 되어 **분할이 전혀 바뀌지 않은
    # 구간이 '경계(급격히 재편됨)' 로 뒤집혀 보고된다.** 기본 스윕 0.10~0.90 의 양쪽 끝이
    # 바로 그 구간이므로, 이 단계의 결론("안정 구간이 있다")이 정반대로 나온다.
    # 보정 사실은 숨기지 않고 ari_degenerate 컬럼으로 산출물에 남긴다.
    n_a, n_b = len(set(la)), len(set(lb))
    degenerate = (n_a == n_b == 1) or (n_a == n_b == n_common)
    if degenerate:
        ari = 1.0

    out.update(ari=ari, ari_degenerate=degenerate,
               nmi=normalized_mutual_info(la, lb),
               **crosstab_mapping(la, lb))
    return out

REFERENCE_THRESHOLD = 0.70            # master_1827.tsv 의 219 분할이 만들어진 값

def consecutive_runs(flags: list[bool]) -> list[tuple[int, int]]:
    """True 가 연속하는 구간의 (시작 index, 끝 index) 목록."""
    runs, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = i
        elif not f and start is not None:
            runs.append((start, i - 1))
            start = None
    if start is not None:
        runs.append((start, len(flags) - 1))
    return runs

def label_to_threshold(label: str) -> float:
    """'id70' -> 0.70, 'id12p5' -> 0.125. threshold_label 의 역함수."""
    m = re.fullmatch(r"id(\d+(?:p\d+)?)", label)
    if not m:
        raise ValueError(f"임계값 라벨을 해석할 수 없다: {label!r} (기대 형식 'id70' / 'id12p5')")
    return float(m.group(1).replace("p", ".")) / 100.0

def interpret_adjacent(df: pd.DataFrame, order: list[str], *, stable_ari: float) -> dict:
    """ARI 가 유지되는 '안정 구간' 과 급락하는 '경계' 를 뽑는다. 요약 md 의 해석 재료."""
    out = dict(stable_ari=stable_ari, n_pairs=int(len(df)))
    # 임계값 목록에서 바로 나오는 값이므로 인접 쌍이 없어도(임계값 1개) 정확히 적는다.
    # df.empty 분기에서 빠뜨리면 '0.70 하나만 돌린' 실행이 요약 1절에 「0.70 포함: 아니오」로 적힌다.
    thr = [label_to_threshold(lb) for lb in order]
    out["reference_threshold"] = REFERENCE_THRESHOLD
    out["reference_threshold_present"] = any(
        abs(t - REFERENCE_THRESHOLD) < 1e-9 for t in thr)

    if df.empty:
        out.update(stable_intervals=[], boundaries=[], ari_median=None, ari_min=None,
                   ari_max=None, reference_in_stable_interval=None,
                   reference_in_stable_interval_at_edge=None, ref_merge_steps=[],
                   n_ari_nan=0, n_ari_degenerate=0, ref_membership_intervals=[])
        return out

    # 중앙값·최솟값은 NaN 을 건너뛴 값이다. 몇 쌍이 NaN 이었는지 같이 적지 않으면
    # "n_pairs 쌍의 중앙값" 으로 오독된다.
    out["n_ari_nan"] = int(df["ari"].isna().sum())
    out["n_ari_degenerate"] = int(df["ari_degenerate"].sum()) if "ari_degenerate" in df.columns else 0
    out["ari_median"] = float(df["ari"].median())
    out["ari_min"] = float(df["ari"].min())
    out["ari_max"] = float(df["ari"].max())

    intervals = []
    for i, j in consecutive_runs(df["is_stable"].tolist()):
        intervals.append(dict(low=thr[i], high=thr[j + 1], n_steps=int(j - i + 1),
                              min_ari=float(df["ari"].iloc[i:j + 1].min()),
                              median_ari=float(df["ari"].iloc[i:j + 1].median())))
    out["stable_intervals"] = sorted(intervals, key=lambda d: -(d["high"] - d["low"]))

    bnd = df[~df["is_stable"]].sort_values("ari")
    out["boundaries"] = [dict(low=float(r.min_seq_id_low), high=float(r.min_seq_id_high),
                              ari=float(r.ari), nmi=float(r.nmi),
                              n_coarser_merged=int(r.n_coarser_merged))
                         for r in bnd.itertuples()]

    # 운용 임계값이 안정 구간 안에 있는가 — 이번 단계의 핵심 주장
    ref = REFERENCE_THRESHOLD
    hit = [iv for iv in out["stable_intervals"]
           if iv["low"] - 1e-9 <= ref <= iv["high"] + 1e-9]
    out["reference_threshold"] = ref
    out["reference_in_stable_interval"] = hit[0] if hit else None
    out["reference_threshold_present"] = any(abs(t - ref) < 1e-9 for t in thr)
    # 0.70 이 구간의 **끝점**이면 한쪽으로 한 스텝만 움직여도 분할이 바뀐다.
    # "구간 어디를 골라도 같다" 는 주장과 의미가 다르므로 구분해서 적는다.
    out["reference_in_stable_interval_at_edge"] = (
        bool(abs(hit[0]["low"] - ref) < 1e-9 or abs(hit[0]["high"] - ref) < 1e-9) if hit else None)

    # ChCODH2 클러스터 멤버가 전혀 바뀌지 않는 구간
    if "ref_jaccard" in df.columns:
        same = [(not np.isnan(v)) and abs(v - 1.0) < 1e-12 for v in df["ref_jaccard"].tolist()]
        out["ref_membership_intervals"] = [
            dict(low=thr[i], high=thr[j + 1], n_steps=int(j - i + 1))
            for i, j in consecutive_runs(same)]
        merge_steps = df[(df["ref_n_source_clusters"] > 1)]
        out["ref_merge_steps"] = [dict(low=float(r.min_seq_id_low), high=float(r.min_seq_id_high),
                                       n_source_clusters=int(r.ref_n_source_clusters),
                                       n_added=int(r.ref_n_added),
                                       merged_from=str(r.ref_merged_from))
                                  for r in merge_steps.itertuples()]
    else:
        out["ref_membership_intervals"] = []
        out["ref_merge_steps"] = []
    return out

def adjacent_table(partitions: dict[str, pd.Series], order: list[str], *,
                   stable_ari: float, ref_acc: str, max_list: int) -> tuple[pd.DataFrame, dict]:
    """인접 임계값 쌍의 안정성. 낮은 임계값 = coarser, 높은 임계값 = finer 로 고정 표기."""
    rows = []
    for lo_label, hi_label in zip(order, order[1:]):
        coarse = partitions[lo_label]        # 낮은 min-seq-id → 더 큰 묶음
        fine = partitions[hi_label]
        c = compare_partitions(fine, coarse)  # old=finer, new=coarser → n_merged 가 의미를 가진다
        row = dict(
            label_low=lo_label, label_high=hi_label,
            min_seq_id_low=label_to_threshold(lo_label),
            min_seq_id_high=label_to_threshold(hi_label),
            n_common=c["n_common"],
            ari=c["ari"], nmi=c["nmi"],
            n_clusters_coarser=int(coarse.nunique()), n_clusters_finer=int(fine.nunique()),
            delta_n_clusters=int(fine.nunique() - coarse.nunique()),
            finer_maps_to_one=c.get("old_maps_to_one", 0),
            coarser_maps_to_one=c.get("new_maps_to_one", 0),
            one_to_one=c.get("one_to_one", 0),
            n_finer_split=c.get("n_split", 0),
            n_coarser_merged=c.get("n_merged", 0),
            is_stable=bool(not np.isnan(c["ari"]) and c["ari"] >= stable_ari),
            ari_degenerate=bool(c.get("ari_degenerate", False)),
        )
        # ChCODH2 클러스터의 변화 — 높은 임계값에서 낮은 임계값으로 내려갈 때
        m_hi = cluster_members(fine, ref_acc)
        m_lo = cluster_members(coarse, ref_acc)
        if m_hi is None or m_lo is None:
            row.update(ref_jaccard=float("nan"), ref_n_added=0, ref_n_removed=0,
                       ref_n_source_clusters=0, ref_merged_from="", ref_merged_truncated=False)
        else:
            s_hi, s_lo = set(m_hi), set(m_lo)
            union = s_hi | s_lo
            src = sorted({str(fine.get(a, "")) for a in s_lo if a in fine.index})
            src = [s for s in src if s]
            row.update(
                ref_jaccard=float(len(s_hi & s_lo) / len(union)) if union else float("nan"),
                ref_n_added=len(s_lo - s_hi), ref_n_removed=len(s_hi - s_lo),
                ref_n_source_clusters=len(src),
                ref_merged_from=";".join(src[:max_list]),
                ref_merged_truncated=bool(len(src) > max_list),
            )
        rows.append(row)

    df = pd.DataFrame(rows, columns=ADJACENT_COLUMNS)
    missing = [c for c in ADJACENT_COLUMNS if c not in df.columns]
    if missing:
        raise AssertionError(f"sweep_adjacent.tsv 컬럼 누락: {missing}")
    interp = interpret_adjacent(df, order, stable_ari=stable_ari)
    return df, interp

def fmt_min_seq_id(t: float) -> str:
    """mmseqs 에 넘기는 문자열. 0.70 은 인계문서대로 '0.70' 으로 적는다."""
    s = f"{t:.2f}"
    if abs(float(s) - t) < 1e-12:
        return s
    return f"{t:.6f}".rstrip("0").rstrip(".")

def build_cluster_cmd(mmseqs: str, fasta: Path, prefix: Path, tmpdir: Path, t: float,
                      *, cov: str, cov_mode: str, sensitivity: str, threads: int) -> list[str]:
    """인계문서 1절의 명령 그대로. --min-seq-id 만 스윕 값으로 바뀐다."""
    return [str(mmseqs), "easy-cluster", str(fasta), str(prefix), str(tmpdir),
            "--min-seq-id", fmt_min_seq_id(t),
            "-c", str(cov), "--cov-mode", str(cov_mode), "-s", str(sensitivity),
            "--threads", str(threads)]

def dir_bytes(path: Path) -> int:
    total = 0
    try:
        for p in Path(path).rglob("*"):
            try:
                if p.is_file():
                    total += p.stat().st_size
            except OSError:
                continue
    except OSError:
        return -1
    return total

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

def safe_parse_accession(token) -> str:
    """parse_accession 과 같은 규칙. 빈 셀에서 IndexError 를 내지 않는다.

    (parse_accession 은 빈 문자열에서 터진다. 단계 1도 같은 래퍼를 들고 있으므로
     py 에 올릴 후보다.)
    """
    token = "" if token is None else str(token).strip()
    if not token:
        return ""
    return parse_accession(token)

def fasta_universe(fasta: Path) -> tuple[list[str], dict]:
    accs, n_records = [], 0
    for header, _seq in read_fasta(fasta):
        n_records += 1
        accs.append(safe_parse_accession(header))
    counts = pd.Series(accs).value_counts() if accs else pd.Series(dtype=int)
    uniq = sorted({a for a in accs if a})
    meta = dict(source="fasta", path=str(fasta), n_records=n_records,
                n_unique_accession=len(uniq),
                n_empty_accession=int(sum(1 for a in accs if not a)),
                n_duplicate_accession=int((counts > 1).sum()) if len(counts) else 0,
                duplicate_examples=sorted(counts[counts > 1].index.tolist())[:10] if len(counts) else [])
    return uniq, meta

def threshold_label(t: float) -> str:
    """0.70 -> 'id70', 0.125 -> 'id12p5'. 파일 이름과 열 이름의 단일 출처."""
    pct = t * 100.0
    if abs(pct - round(pct)) < 1e-9:
        return f"id{int(round(pct)):02d}"
    s = f"{pct:.4f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"id{s}"

def make_labels(thresholds) -> list[str]:
    """임계값 목록 -> 라벨 목록. 충돌과 왕복 손실을 여기서 막는다.

    라벨이 곧 파일 이름(id70_cluster.tsv)이고 열 이름(cluster@70)이므로, 서로 다른
    임계값이 같은 라벨을 받으면 **조용히 덮어쓰기**가 난다. 그래서 두 가지를 검사한다.
      ① 라벨 유일성       ② label_to_threshold(threshold_label(t)) == t (왕복 무손실)
    """
    labels = [threshold_label(t) for t in thresholds]
    if len(set(labels)) != len(labels):
        dup = sorted({lb for lb in labels if labels.count(lb) > 1})
        raise ValueError(
            f"임계값 라벨이 충돌한다: {dup}. 서로 다른 임계값이 같은 파일 이름(id<XX>_cluster.tsv)과 "
            f"같은 열 이름(cluster@XX)을 쓰게 되므로 조용히 덮어쓰기가 난다. 임계값을 조정할 것.")
    broken = [(t, lb) for t, lb in zip(thresholds, labels)
              if abs(label_to_threshold(lb) - t) > 1e-9]
    if broken:
        raise ValueError(
            f"임계값과 라벨의 왕복이 깨진다: {broken}. 라벨에서 임계값을 복원할 수 없으면 "
            f"--from-existing 이 잘못된 임계값을 읽는다.")
    return labels

def partition_from_cluster_df(df: pd.DataFrame, label: str, *,
                              allow_duplicate_accession: bool) -> tuple[pd.Series, dict]:
    """cluster.tsv -> accession 색인 Series. 중복 accession 을 조용히 넘기지 않는다."""
    dup_mask = df["accession"].duplicated(keep=False)
    n_dup_rows = int(dup_mask.sum())
    conflict = []
    if n_dup_rows:
        g = df[dup_mask].groupby("accession")["cluster"].nunique()
        conflict = sorted(g[g > 1].index.tolist())

    if conflict and not allow_duplicate_accession:
        raise ValueError(
            f"[{label}] 같은 accession 이 서로 다른 클러스터 둘 이상에 속한다 "
            f"({len(conflict)}개, 예: {conflict[:5]}). 와이드 테이블이 조용히 틀어지므로 중단한다. "
            f"입력 FASTA 에 중복 accession 이 있는지 먼저 확인하고, 그래도 진행하려면 "
            f"--allow-duplicate-accession (정렬 후 첫 배정만 사용).")

    ordered = df.sort_values(["accession", "cluster"], kind="mergesort")
    dedup = ordered.drop_duplicates("accession", keep="first")
    ser = pd.Series(dedup["cluster"].to_numpy(), index=dedup["accession"].to_numpy(), name=label)

    meta = dict(
        n_duplicate_accession_rows=n_dup_rows,
        n_conflicting_accession=len(conflict),
        conflicting_accession_examples=conflict[:10],
        conflicting_accession_truncated=bool(len(conflict) > 10),
        n_assignments=int(len(ser)),
    )
    return ser, meta

HEADER_TOKENS = {"cluster", "cluster_rep", "rep", "representative", "accession",
                 "member", "target", "query", "acc"}

def norm_acc_series(series: pd.Series) -> pd.Series:
    s = series.astype("object").fillna("").astype(str).str.strip()
    return s.map(safe_parse_accession)

def read_cluster_tsv(path: Path) -> tuple[pd.DataFrame, dict]:
    """easy-cluster 의 _cluster.tsv (헤더 없는 2열: 대표 \t 멤버) 를 읽는다.

    read_tsv 를 쓰지 않는 이유 — 그 함수는 '컬럼명에 | 또는 숫자' 로 헤더 유무를
    추정한다. accession 만 든 2열 파일은 그 조건에 걸리지 않아 **첫 서열이 헤더로 먹힌다.**
    여기서는 형식을 알고 있으므로 항상 header=None 으로 읽고, 사람이 붙인 헤더만 걸러낸다.
    """
    df = pd.read_csv(path, sep="\t", header=None, dtype=str,
                     keep_default_na=False, na_filter=False)
    if df.shape[1] != 2:
        raise ValueError(
            f"{path} 의 컬럼 수가 {df.shape[1]} — easy-cluster _cluster.tsv 는 2열이어야 한다. "
            f"파일이 다른 산출물이거나 구분자가 탭이 아닐 수 있다.")
    df.columns = ["cluster_raw", "accession_raw"]

    header_dropped = False
    if len(df) and all(str(v).strip().lower() in HEADER_TOKENS for v in df.iloc[0].tolist()):
        df = df.iloc[1:].reset_index(drop=True)
        header_dropped = True

    if df.empty:
        raise ValueError(f"{path} 에서 행을 하나도 읽지 못했다.")

    df["cluster"] = norm_acc_series(df["cluster_raw"])
    df["accession"] = norm_acc_series(df["accession_raw"])

    meta = dict(
        path=str(path),
        n_rows=int(len(df)),
        header_row_dropped=header_dropped,
        n_empty_accession=int((df["accession"] == "").sum()),
        n_empty_cluster=int((df["cluster"] == "").sum()),
        n_unique_accession=int(df["accession"].nunique()),
        n_unique_cluster=int(df["cluster"].nunique()),
    )
    if meta["n_empty_accession"] or meta["n_empty_cluster"]:
        raise ValueError(
            f"{path} 에 빈 accession {meta['n_empty_accession']}개 / 빈 cluster "
            f"{meta['n_empty_cluster']}개가 있다. 조인이 조용히 실패하므로 중단한다.")
    return df, meta

def read_markers(path: Path) -> tuple[list[dict], dict]:
    """마커 목록. 한 줄에 accession 하나, '#' 뒤 주석, 뒤따르는 토큰은 라벨."""
    markers, bad_lines, dup_accs = [], [], []
    n_lines = 0
    seen = set()
    for raw in Path(path).read_text(errors="replace").splitlines():
        n_lines += 1
        body = raw.split("#", 1)[0].strip()
        if not body:
            continue
        toks = body.split()
        acc = safe_parse_accession(toks[0])
        if not acc:
            bad_lines.append(raw[:120])
            continue
        label = " ".join(toks[1:]).strip()
        if acc in seen:
            dup_accs.append(acc)
            continue
        seen.add(acc)
        markers.append(dict(marker_accession=acc, marker_label=label or acc))
    meta = dict(
        path=str(path),
        n_lines=n_lines,
        n_markers=len(markers),
        n_unparsable_lines=len(bad_lines),
        unparsable_examples=bad_lines[:5],
        n_duplicates_dropped=len(dup_accs),
        duplicate_examples=sorted(set(dup_accs))[:5],
    )
    return markers, meta

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

def read_tsv(path, **kw) -> pd.DataFrame:
    """헤더가 있든 없든 읽는다. mmseqs convertalis 산출물은 헤더가 없을 수 있다."""
    df = pd.read_csv(path, sep="\t", dtype=str, **kw)
    # 첫 행이 데이터처럼 보이면 헤더 없는 파일이다 (컬럼명에 '|' 또는 숫자만)
    cols = list(df.columns)
    if any(("|" in c) or re.fullmatch(r"[\d.eE+-]+", c) for c in cols):
        df = pd.read_csv(path, sep="\t", dtype=str, header=None, **kw)
        df.columns = [f"col{i}" for i in range(df.shape[1])]
    return df

def detect_column(df: pd.DataFrame, candidates, *, required=True, what=""):
    """후보 이름 중 실제로 존재하는 첫 컬럼명을 돌려준다(대소문자 무시)."""
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise KeyError(f"{what or candidates[0]} 컬럼을 못 찾음. 후보={candidates} 실제={list(df.columns)}")
    return None

def resolve_columns(df: pd.DataFrame, keys, *, required=()) -> dict:
    """COLSPEC 의 논리 이름 -> 실제 컬럼명. 없으면 None (required 면 예외)."""
    out = {}
    for k in keys:
        out[k] = detect_column(df, COLSPEC[k], required=(k in required), what=k)
    return out

def read_master_partition(path: Path, *, allow_duplicate_accession: bool) -> tuple[pd.Series | None, dict]:
    """master_1827.tsv 의 cluster 컬럼을 참조 분할로 읽는다. 읽기 전용."""
    info = dict(path=str(path), exists=bool(Path(path).exists()), skipped=False)
    if not Path(path).exists():
        info.update(skipped=True, reason="파일 없음")
        return None, info

    m = read_tsv(path)
    info["columns"] = list(m.columns)
    info["n_rows"] = int(len(m))
    cols = resolve_columns(m, ["accession", "cluster"], required=("accession",))
    acc_col, cl_col = cols["accession"], cols["cluster"]
    info.update(accession_column=acc_col, cluster_column=cl_col)
    if cl_col is None:
        raise KeyError(
            f"master 에서 cluster 컬럼을 못 찾음. 후보={COLSPEC['cluster']} 실제={list(m.columns)}")

    acc = norm_acc_series(m[acc_col])
    cl = norm_acc_series(m[cl_col])
    keep = (acc != "") & (cl != "")
    info["n_rows_dropped_empty"] = int((~keep).sum())
    sub = pd.DataFrame({"accession": acc[keep], "cluster": cl[keep]})

    dup = sub["accession"].duplicated(keep=False)
    conflict = []
    if dup.any():
        g = sub[dup].groupby("accession")["cluster"].nunique()
        conflict = sorted(g[g > 1].index.tolist())
    info.update(n_duplicate_accession_rows=int(dup.sum()), n_conflicting_accession=len(conflict),
                conflicting_accession_examples=conflict[:10])
    if conflict and not allow_duplicate_accession:
        raise ValueError(
            f"master 에서 같은 accession 이 서로 다른 cluster 둘 이상에 속한다 "
            f"({len(conflict)}개, 예: {conflict[:5]}). 참조 분할이 모호하므로 중단한다. "
            f"--allow-duplicate-accession 으로 첫 배정만 쓸 수 있다.")

    sub = sub.sort_values(["accession", "cluster"], kind="mergesort").drop_duplicates("accession")
    ser = pd.Series(sub["cluster"].to_numpy(), index=sub["accession"].to_numpy(), name="master")
    info.update(n_assignments=int(len(ser)), n_clusters=int(ser.nunique()))
    return ser, info

def track_markers(partitions: dict[str, pd.Series], order: list[str], markers: list[dict],
                  *, max_list: int) -> tuple[pd.DataFrame, dict, dict]:
    rows, per_label = [], {}
    accs = [m["marker_accession"] for m in markers]
    labels = {m["marker_accession"]: m["marker_label"] for m in markers}

    for label in order:
        ser = partitions[label]
        present = [a for a in accs if a in ser.index]
        by_cluster: dict[str, list[str]] = {}
        for a in present:
            by_cluster.setdefault(str(ser.loc[a]), []).append(a)
        groups = {rep: i + 1 for i, rep in enumerate(sorted(by_cluster))}
        sizes = ser.value_counts()
        for a in accs:
            if a not in ser.index:
                rows.append(dict(label=label, min_seq_id=label_to_threshold(label),
                                 marker_accession=a, marker_label=labels[a], found=False,
                                 cluster="", cluster_size=0, co_marker_group=0,
                                 n_co_markers=0, co_markers="", co_markers_truncated=False))
                continue
            rep = str(ser.loc[a])
            co = sorted(x for x in by_cluster[rep] if x != a)
            rows.append(dict(label=label, min_seq_id=label_to_threshold(label),
                             marker_accession=a, marker_label=labels[a], found=True,
                             cluster=rep, cluster_size=int(sizes.get(rep, 0)),
                             co_marker_group=groups[rep], n_co_markers=len(co),
                             co_markers=";".join(co[:max_list]),
                             co_markers_truncated=bool(len(co) > max_list)))
        per_label[label] = dict(n_markers_found=len(present),
                                n_marker_clusters=len(by_cluster),
                                markers_all_together=("yes" if len(by_cluster) == 1 and len(present) > 1
                                                      else ("n/a" if len(present) < 2 else "no")))

    df = pd.DataFrame(rows)

    # 마커 쌍이 처음으로 한 묶음이 되는 가장 높은 임계값 — "어디서 갈라지나" 의 정답
    pair_rows = []
    for a, b in combinations(accs, 2):
        same_at = []
        for label in order:
            ser = partitions[label]
            if a in ser.index and b in ser.index and str(ser.loc[a]) == str(ser.loc[b]):
                same_at.append(label_to_threshold(label))
        pair_rows.append(dict(marker_a=a, marker_b=b,
                              label_a=labels[a], label_b=labels[b],
                              n_thresholds_together=len(same_at),
                              highest_together=max(same_at) if same_at else None,
                              lowest_together=min(same_at) if same_at else None))
    info = dict(n_markers=len(accs),
                markers=[dict(m) for m in markers],
                n_markers_in_universe=int(df.loc[df.found, "marker_accession"].nunique()) if len(df) else 0,
                missing_markers=sorted({a for a in accs
                                        if not any(a in partitions[lb].index for lb in order)}),
                pairs=pair_rows, n_pairs=len(pair_rows), max_list=max_list)
    return df, per_label, info

def track_reference(partitions: dict[str, pd.Series], order: list[str], ref_acc: str,
                    *, max_members: int) -> tuple[pd.DataFrame, dict, dict]:
    rows, per_label = [], {}
    found_any = False
    for label in order:
        ser = partitions[label]
        members = cluster_members(ser, ref_acc)
        if members is None:
            rows.append(dict(label=label, min_seq_id=label_to_threshold(label), ref_found=False,
                             ref_cluster="", ref_cluster_size=0, ref_is_mm_rep=False,
                             n_members_listed=0, members_truncated=False, members=""))
            per_label[label] = dict(ref_found=False, ref_cluster="", ref_cluster_size=0,
                                    ref_is_mm_rep=False)
            continue
        found_any = True
        rep = str(ser.loc[ref_acc])
        shown = members[:max_members]
        rows.append(dict(label=label, min_seq_id=label_to_threshold(label), ref_found=True,
                         ref_cluster=rep, ref_cluster_size=len(members),
                         ref_is_mm_rep=bool(rep == ref_acc),
                         n_members_listed=len(shown),
                         members_truncated=bool(len(members) > max_members),
                         members=";".join(shown)))
        per_label[label] = dict(ref_found=True, ref_cluster=rep, ref_cluster_size=len(members),
                                ref_is_mm_rep=bool(rep == ref_acc))
    df = pd.DataFrame(rows)
    info = dict(reference_accession=ref_acc, found_in_any_threshold=found_any,
                n_thresholds_found=int(df["ref_found"].sum()) if len(df) else 0,
                max_members_listed=max_members,
                size_min=int(df.loc[df.ref_found, "ref_cluster_size"].min()) if found_any else None,
                size_max=int(df.loc[df.ref_found, "ref_cluster_size"].max()) if found_any else None)
    return df, per_label, info

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
