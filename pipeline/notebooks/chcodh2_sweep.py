"""ChCODH2 상동체 클러스터 스윕 — 노트북 전용 단일 모듈.

cluster_sweep.ipynb 옆에 두기만 하면 된다. git 도, 다른 파일도 필요 없다.
의존성은 pandas + numpy + 표준 라이브러리뿐.

이 파일은 pipeline/lib/common.py 와 pipeline/stages/02_cluster_sweep.py 에서
노트북이 쓰는 함수만 뽑아 합친 것이다. 원본은 git 레포에 있고 거기서 검증됐다.
직접 고치지 말 것 — 다음 갱신 때 덮어쓴다.

2026-10-25c 에서 cluster_set_compare.ipynb 용 세트 비교 함수를 더했다 (파일 뒤쪽).
정렬 파일 읽기는 pipeline/stages/03_equiv_map.py 에서 검증된 로직을 옮긴 것이다.
앞쪽 함수(01 번 노트북이 쓰는 것)는 전혀 건드리지 않았다.
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

VERSION = "2026-10-25c"   # 노트북 1절이 이 값을 대조한다. 어긋나면 둘 중 하나가 낡은 것이다.

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


# ==========================================================================
# 보존도 분석 대상 세트 고르기 (2026-10-25c 추가)
#
# cluster_set_compare.ipynb 가 쓰는 부분. 두 방식을 **같은 잣대로** 재는 것이 목적이다.
#   A. 클러스터 기반 — ChCODH2 가 속한 mmseqs 클러스터의 멤버
#   B. identity 기반 — ChCODH2 와의 fident 가 기준 이상인 서열 전부
# 어느 쪽이 낫다고 여기서 정하지 않는다. 숫자만 같은 모양으로 돌려준다.
#
# 정렬 파일 읽기(헤더 판정·--columns·HSP 중복 처리)는 pipeline/stages/03_equiv_map.py
# 에서 검증된 로직을 그대로 옮겼다. 다만 이 노트북은 qaln/taln 을 쓰지 않으므로
# 필수 컬럼이 target·fident 두 개뿐이고, 헤더도 --columns 도 없을 때는 복원된
# --format-output 순서를 **추정으로 쓰되 그 사실을 반드시 보고한다.**
# ==========================================================================

# 복원된 easy-search 명령의 --format-output 순서. 헤더도 columns 도 없을 때의 추정값.
ALIGN_FORMAT_OUTPUT = ["query", "target", "fident", "alnlen", "qcov", "tcov",
                       "evalue", "bits", "qstart", "qend", "tstart", "tend",
                       "qaln", "taln"]

# 헤더 행 판정 어휘 (03_equiv_map.HEADER_FIELD_TOKENS 와 같은 목록).
ALIGN_HEADER_TOKENS = {
    "query", "target", "evalue", "gapopen", "pident", "fident", "nident", "qstart",
    "qend", "qlen", "tstart", "tend", "tlen", "alnlen", "raw", "bits", "cigar",
    "qseq", "tseq", "qheader", "theader", "qaln", "taln", "mismatch", "qcov", "tcov",
    "qset", "qsetid", "tset", "tsetid", "taxid", "taxname", "taxlineage", "accession",
    "cluster", "bitscore", "sseqid", "qseqid", "subject", "score", "length", "identity",
}

# master_1827.tsv 의 부가 컬럼 후보. 실제 컬럼명이 미확인이라 후보 탐지로 찾는다.
MASTER_FEATURE_SPEC = {
    "cluster":   ["cluster", "cluster_rep", "rep", "representative"],
    "is_mm_rep": ["is_mm_rep", "is_rep", "mm_rep", "is_representative", "rep_flag"],
    "cf":        ["cf", "cluster_fident", "cluster_median_fident", "median_fident",
                  "fident_cluster"],
    "cband":     ["cband", "c_band", "cluster_band", "band", "cf_band"],
    "dtype":     ["dtype", "d_type", "d_cluster", "dcluster", "d_cluster_type",
                  "dcluster_type", "d_site", "dcluster_sf", "sf_type"],
}

# 논문에서 확인된 잔기 번호. master 의 r<NNN> 컬럼과 대조한다 (r612 vs V610 문제).
PAPER_POSITIONS = (41, 43, 57, 59, 583, 610)

_POSITION_COL_RE = re.compile(r"^r(\d+)$", re.IGNORECASE)


def _to_float(series: pd.Series) -> pd.Series:
    """문자열 컬럼 -> float64. 빈 칸과 숫자가 아닌 값은 NaN. dtype 을 float64 로 고정해
    pandas 버전에 따라 nullable 로 올라오는 것을 막는다(마스크에 pd.NA 가 섞이면 터진다)."""
    s = pd.to_numeric(series.astype("object").astype(str).str.strip(), errors="coerce")
    return s.astype("float64")


def looks_like_align_header(tokens) -> tuple[bool, str]:
    """정렬 TSV 첫 행이 헤더인지 판정. 03_equiv_map.looks_like_header 와 같은 규칙.

    모든 토큰이 맨 식별자이고 그중 둘 이상이 mmseqs 필드 어휘에 있으면 헤더다.
    데이터 행은 fident/evalue 같은 숫자를 반드시 포함하므로 첫 조건에서 걸러진다.
    """
    toks = [str(t).strip() for t in tokens]
    if not toks:
        return False, "빈 행"
    ident = [bool(re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", t)) for t in toks]
    if not all(ident):
        bad = [t for t, k in zip(toks, ident) if not k][:5]
        return False, f"식별자가 아닌 토큰 {bad} 이 있어 데이터 행으로 본다"
    known = [t for t in toks if t.lower() in ALIGN_HEADER_TOKENS]
    if len(known) >= 2:
        return True, f"mmseqs 필드 어휘 일치 {len(known)}/{len(toks)}"
    return False, f"필드 어휘 일치가 {len(known)}개뿐이라 데이터 행으로 본다"


def _reject_duplicate_names(names, where: str) -> None:
    """컬럼명이 중복되면 df[name] 이 Series 가 아니라 DataFrame 이 되어 한참 뒤
    엉뚱한 곳에서 터진다. 어느 이름이 몇 번 나왔는지 여기서 바로 알려준다."""
    vc = pd.Series([str(n).strip().lower() for n in names]).value_counts()
    dup = {k: int(v) for k, v in vc[vc > 1].items()}
    if dup:
        raise ValueError(
            f"{where} 에 같은 컬럼명이 두 번 이상 나온다: {dup}. 어느 열이 fident 인지 "
            f"확정할 수 없으므로 중단한다. columns 인자로 열 순서를 직접 줄 것. "
            f"실제 이름 = {list(names)[:20]}")


def read_alignment(path, columns=None, *, ref_acc: str = CHCODH2_ACC,
                   assume_format_output: bool = True) -> tuple[pd.DataFrame, dict]:
    """chcodh2_alignment.tsv (ChCODH2 질의 easy-search 결과) 를 읽는다.

    반환 DataFrame 은 **target accession 색인, target 당 1행** (HSP 여러 개면 bits 최대
    1건만 남긴다). 컬럼은 fident(0~1 분수로 정규화) + 있으면 bits/evalue/qcov/tcov,
    그리고 몇 개의 HSP 에서 골랐는지 n_hsp.

    열 이름을 정하는 세 갈래 — 어느 쪽을 썼는지 meta['column_source'] 에 남는다.
      ① columns 인자가 있으면 그것 (--format-output 순서를 손으로 주는 경우)
      ② 파일에 헤더가 있으면 헤더
      ③ 둘 다 없고 열 수가 ALIGN_FORMAT_OUTPUT 과 같으면 **추정**.
         이때는 fident 열이 전부 숫자인지까지 확인하고, 하나라도 아니면 예외.
    """
    path = Path(path)
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str,
                      keep_default_na=False, na_filter=False)
    if raw.empty:
        raise ValueError(f"{path} 에서 행을 하나도 읽지 못했다. 파일이 비었거나 구분자가 탭이 아니다.")

    first = [str(v).strip() for v in raw.iloc[0].tolist()]
    hdr_guess, hdr_reason = looks_like_align_header(first)
    meta: dict = dict(path=str(path), n_raw_rows=int(len(raw)), n_fields=int(raw.shape[1]),
                      first_row=first[:20], header_guess=hdr_guess, header_reason=hdr_reason)

    if columns:
        names = [str(c).strip() for c in columns]
        if len(names) != raw.shape[1]:
            raise ValueError(
                f"columns 가 {len(names)}개인데 {path} 는 {raw.shape[1]}열이다. "
                f"--format-output 순서를 그대로 적었는지 확인할 것. 첫 행 = {first[:20]}")
        _reject_duplicate_names(names, "columns 인자")
        dropped = all(a.lower() == b.lower() for a, b in zip(first, names))
        body = raw.iloc[1:] if dropped else raw
        meta.update(column_source="columns 인자", header_row_dropped=bool(dropped))
    elif hdr_guess:
        names = first
        _reject_duplicate_names(names, f"{path} 의 헤더 행")
        body = raw.iloc[1:]
        meta.update(column_source="파일 헤더", header_row_dropped=True)
    elif assume_format_output and raw.shape[1] == len(ALIGN_FORMAT_OUTPUT):
        names = list(ALIGN_FORMAT_OUTPUT)
        body = raw
        meta.update(column_source="추정(복원된 --format-output 순서)", header_row_dropped=False)
    else:
        raise ValueError(
            f"{path} 에 헤더가 없다고 판정했다 ({hdr_reason}). {raw.shape[1]}열이며 첫 행 = "
            f"{first[:20]}.\n어느 열이 fident 인지 알 수 없으므로 columns 인자로 "
            f"--format-output 순서를 직접 줄 것. 예:\n  columns=cs.ALIGN_FORMAT_OUTPUT")

    df = body.reset_index(drop=True)
    df.columns = list(names)
    if df.empty:
        raise ValueError(f"{path} 에 헤더만 있고 데이터 행이 없다.")
    meta["columns"] = list(df.columns)

    cmap = resolve_columns(df, ["accession", "query", "fident", "qcov", "tcov",
                                "evalue", "bits"], required=("accession", "fident"))
    meta["column_map"] = {k: v for k, v in cmap.items()}

    # 질의가 둘 이상이면 fident 가 '무엇 대비' 인지 달라진다 — 조용히 틀리기 전에 멈춘다.
    if cmap.get("query"):
        qvals = sorted({str(v).strip() for v in df[cmap["query"]].tolist() if str(v).strip()})
        meta.update(n_query_values=len(qvals), query_values=qvals[:5])
        if len(qvals) > 1:
            raise ValueError(
                f"{path} 에 질의가 {len(qvals)}종류 있다 ({qvals[:5]}). 이 노트북은 ChCODH2 "
                f"단일 질의 결과를 전제하므로 fident 의 뜻이 달라진다. 중단한다.")
    else:
        meta.update(n_query_values=None, query_values=[])

    fid = _to_float(df[cmap["fident"]])
    n_bad = int(fid.isna().sum())
    strict = str(meta["column_source"]).startswith("추정")
    if n_bad and (strict or n_bad > max(1, int(0.01 * len(df)))):
        ex = df.loc[fid.isna(), cmap["fident"]].astype(str).head(5).tolist()
        raise ValueError(
            f"fident 컬럼('{cmap['fident']}') 에 숫자가 아닌 값이 {n_bad}/{len(df)}개 있다 "
            f"(예: {ex}). 열 이름을 {meta['column_source']} 로 정했는데 열이 어긋났을 수 있다. "
            f"columns 인자로 --format-output 순서를 직접 줄 것.")

    fmax_in = float(fid.max()) if int(fid.notna().sum()) else float("nan")
    is_percent = bool(fmax_in > 1.0)
    if is_percent:
        if fmax_in > 100.0 + 1e-6:
            raise ValueError(
                f"fident 최댓값이 {fmax_in} 이다. 0~1 분수도 0~100 퍼센트도 아니므로 "
                f"열이 어긋난 것으로 본다. columns 인자로 열 순서를 직접 줄 것.")
        fid = fid / 100.0
    if int(fid.notna().sum()) and float(fid.min()) < -1e-9:
        raise ValueError(f"fident 에 음수가 있다 (최솟값 {float(fid.min())}). 열이 어긋났다.")
    meta.update(fident_scale="percent(100으로 나눴다)" if is_percent else "fraction(0~1)",
                fident_is_percent=is_percent, fident_input_max=fmax_in,
                n_fident_missing=n_bad)

    acc = norm_acc_series(df[cmap["accession"]])
    n_empty = int((acc == "").sum())
    if n_empty:
        raise ValueError(
            f"target accession 이 빈 행이 {n_empty}개 있다. 조인이 조용히 실패하므로 중단한다.")

    work = pd.DataFrame({"accession": acc.to_numpy(), "fident": fid.to_numpy()})
    for key in ("bits", "evalue", "qcov", "tcov"):
        if cmap.get(key):
            work[key] = _to_float(df[cmap[key]]).to_numpy()
    work["_row"] = np.arange(len(work), dtype=np.int64)
    work["n_hsp"] = work.groupby("accession")["_row"].transform("size").to_numpy()

    keys, asc, used = [], [], []
    if "bits" in work.columns:
        keys.append("bits"); asc.append(False); used.append("bits 최대")
    if "evalue" in work.columns:
        keys.append("evalue"); asc.append(True); used.append("evalue 최소")
    keys += ["fident", "_row"]; asc += [False, True]
    used += ["fident 최대", "파일 등장 순서"]
    best = (work.sort_values(keys, ascending=asc, kind="mergesort")
                .drop_duplicates("accession", keep="first")
                .sort_values("accession", kind="mergesort")
                .set_index("accession")
                .drop(columns=["_row"]))

    n_rows, n_uniq = int(len(work)), int(work["accession"].nunique())
    ok = best["fident"].dropna()
    meta.update(
        n_rows=n_rows, n_unique_target=n_uniq, n_merged=n_rows - n_uniq,
        n_targets_multi_hsp=int((best["n_hsp"] > 1).sum()),
        hsp_tiebreak=" → ".join(used),
        fident_min=float(ok.min()) if len(ok) else None,
        fident_median=float(ok.median()) if len(ok) else None,
        fident_max=float(ok.max()) if len(ok) else None,
        reference_accession=ref_acc,
        reference_in_alignment=bool(ref_acc in best.index),
        reference_fident=(float(best.loc[ref_acc, "fident"])
                          if ref_acc in best.index and pd.notna(best.loc[ref_acc, "fident"])
                          else None),
    )
    meta["reference_fident_is_one"] = (None if meta["reference_fident"] is None
                                       else bool(abs(meta["reference_fident"] - 1.0) < 1e-6))
    return best, meta


def alignment_coverage(align_df: pd.DataFrame, universe, *, max_examples: int = 10) -> dict:
    """정렬 결과가 universe(FASTA 1,827) 를 몇 개 덮는지. 절단 여부를 여기서 센다."""
    uni = pd.Index(sorted(set(universe)))
    hit = pd.Index(align_df.index)
    inter = uni.intersection(hit)
    missing = uni.difference(hit)
    return dict(
        n_universe=int(len(uni)), n_aligned=int(len(hit)),
        n_in_universe=int(len(inter)), n_missing=int(len(missing)),
        frac_covered=float(len(inter) / len(uni)) if len(uni) else float("nan"),
        missing_examples=missing.tolist()[:max_examples],
        missing_truncated=bool(len(missing) > max_examples),
        n_aligned_outside_universe=int(len(hit.difference(uni))),
        outside_examples=hit.difference(uni).tolist()[:max_examples],
    )


def fident_summary(align_df: pd.DataFrame, *, universe=None) -> dict:
    """ChCODH2 대비 identity 분포. 기존 관측(중앙값 0.379, 95% 가 0.55 이하) 과 대조한다."""
    s = align_df["fident"]
    if universe is not None:
        s = s.reindex(pd.Index(sorted(set(universe)))).dropna()
    s = s.dropna().astype("float64")
    if not len(s):
        raise ValueError("fident 값이 하나도 없다 — 정렬 파일을 다시 확인할 것.")
    q = s.quantile([0.05, 0.25, 0.50, 0.75, 0.90, 0.95]).to_dict()
    return dict(
        n=int(len(s)), min=float(s.min()), max=float(s.max()), mean=float(s.mean()),
        q05=float(q[0.05]), q25=float(q[0.25]), median=float(q[0.50]),
        q75=float(q[0.75]), q90=float(q[0.90]), q95=float(q[0.95]),
        n_eq_one=int((s >= 1.0 - 1e-9).sum()),
        frac_le_0p55=float((s <= 0.55).mean()),
    )


def fident_histogram(align_df: pd.DataFrame, *, universe=None, step: float = 0.05) -> pd.DataFrame:
    """fident 히스토그램을 표로. matplotlib 이 없어도 분포를 눈으로 볼 수 있게."""
    s = align_df["fident"]
    if universe is not None:
        s = s.reindex(pd.Index(sorted(set(universe))))
    s = s.dropna().astype("float64")
    edges = np.round(np.arange(0.0, 1.0 + step / 2, step), 6)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        last = bool(abs(hi - edges[-1]) < 1e-12)
        sel = (s >= lo) & ((s <= hi) if last else (s < hi))
        n = int(sel.sum())
        rows.append(dict(bin_low=float(lo), bin_high=float(hi), n=n,
                         frac=float(n / len(s)) if len(s) else float("nan"),
                         n_at_or_above=int((s >= lo).sum())))
    df = pd.DataFrame(rows)
    nmax = int(df["n"].max()) if len(df) else 0
    df["bar"] = [("#" * int(round(40 * n / nmax))) if nmax else "" for n in df["n"]]
    return df


def identity_set(align_df: pd.DataFrame, cutoff: float, universe, *,
                 ref_acc: str = CHCODH2_ACC, include_reference: bool = True) -> set:
    """fident >= cutoff 인 target 집합 (universe 와의 교집합).

    include_reference=True 면 ChCODH2 자신을 **항상 넣는다** (설계 대상이 분석 세트에서
    빠지면 안 되므로). 자기 히트의 fident 는 1.0 이라 보통 어차피 통과한다.
    False 면 반대로 항상 뺀다. 어느 쪽이든 '포함됐는지' 는 set_profile 의
    has_reference 로 다시 확인된다.
    """
    cutoff = float(cutoff)
    if not (0.0 <= cutoff <= 1.0):
        raise ValueError(f"cutoff 는 0~1 분수다 (받은 값 {cutoff}). 퍼센트를 넣지 말 것.")
    if "fident" not in align_df.columns:
        raise KeyError(f"align_df 에 fident 컬럼이 없다. 실제 컬럼 = {list(align_df.columns)}")
    uni = set(universe)
    if not uni:
        raise ValueError("universe 가 비었다 — FASTA 를 먼저 읽을 것.")
    fid = align_df["fident"].astype("float64")
    members = {a for a in fid.index[fid >= cutoff - 1e-12].tolist() if a in uni}
    if include_reference:
        if ref_acc in uni:
            members.add(ref_acc)
    else:
        members.discard(ref_acc)
    return members


def identity_cutoff_table(align_df: pd.DataFrame, cutoffs, universe, *,
                          ref_acc: str = CHCODH2_ACC,
                          include_reference: bool = True) -> pd.DataFrame:
    """cutoff 마다 몇 개가 남는지. identity 기반 세트의 크기를 한눈에 본다."""
    uni = set(universe)
    rows = []
    for c in cutoffs:
        s = identity_set(align_df, c, uni, ref_acc=ref_acc,
                         include_reference=include_reference)
        rows.append(dict(cutoff=float(c), n=len(s),
                         frac_of_universe=float(len(s) / len(uni)) if uni else float("nan"),
                         has_reference=bool(ref_acc in s)))
    return pd.DataFrame(rows)


def cluster_set(cluster_series: pd.Series, ref_acc: str = CHCODH2_ACC) -> set:
    """ref_acc 가 속한 클러스터의 멤버 집합. 없으면 예외 (조용히 빈 세트를 돌려주지 않는다)."""
    members = cluster_members(cluster_series, ref_acc)
    if members is None:
        raise ValueError(
            f"{ref_acc} 가 분할 '{cluster_series.name}' 에 없다 ({len(cluster_series)}개 배정). "
            f"클러스터 파일이 잘렸거나 다른 입력으로 만든 것이다.")
    return set(members)


def load_existing_partition(clusters_dir, t: float, *,
                            allow_duplicate_accession: bool = False) -> tuple[pd.Series, dict]:
    """이미 만들어진 id<XX>_cluster.tsv 를 읽어 분할 Series 로. **mmseqs 를 돌리지 않는다.**"""
    label = threshold_label(t)
    path = Path(clusters_dir) / f"{label}_cluster.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} 가 없다. 이 노트북은 mmseqs 를 돌리지 않는다 — 01 번 묶음이 만든 "
            f"result/clusters/ 폴더를 가리키는지 1절의 CLUSTERS 를 확인할 것.")
    df, meta = read_cluster_tsv(path)
    ser, pmeta = partition_from_cluster_df(
        df, label, allow_duplicate_accession=allow_duplicate_accession)
    meta.update(pmeta, label=label, min_seq_id=float(t), n_clusters=int(ser.nunique()))
    return ser, meta


def link_inputs(targets, src: dict) -> list[str]:
    """input/ 에 없는 입력을 원본에서 심볼릭 링크로 건다 (원본은 읽기만 한다).

    무엇을 했는지/못 했는지 사람이 읽는 줄로 돌려준다. 원본이 없어도 예외를 던지지
    않는다 — 어느 파일이 없는지 2절에서 한 번에 보는 쪽이 낫기 때문이다.
    """
    out = []
    for dst in targets:
        dst = Path(dst)
        if dst.exists():
            out.append(f"있음      {dst.name}")
            continue
        origin = src.get(dst.name)
        if origin is None:
            out.append(f"⚠ 없음    {dst.name} — SRC 에 원본 위치가 없다. 1절에 적을 것")
            continue
        origin = Path(origin)
        if origin.exists():
            dst.parent.mkdir(parents=True, exist_ok=True)
            dst.symlink_to(origin)
            out.append(f"링크 생성  {dst.name}  →  {origin}")
        else:
            out.append(f"⚠ 원본을 못 찾음: {origin}  — 실제 위치를 1절 SRC 에 적을 것")
    return out


def read_master_table(path) -> tuple[pd.DataFrame | None, dict]:
    """master_1827.tsv 를 accession 색인 표로 읽는다 (읽기 전용).

    cf·cband·r41…·D-cluster 유형 컬럼은 실제 이름이 미확인이므로 후보 탐지로 찾고,
    **없으면 건너뛰되 meta['missing'] 에 반드시 남긴다.**
    """
    path = Path(path)
    info: dict = dict(path=str(path), exists=bool(path.exists()), skipped=False)
    if not path.exists():
        info.update(skipped=True, reason="파일 없음")
        return None, info

    m = read_tsv(path)
    info.update(n_rows=int(len(m)), columns=list(m.columns))
    acc_col = detect_column(m, COLSPEC["accession"], required=True, what="master accession")
    info["accession_column"] = acc_col

    found, missing = {}, []
    for key, cands in MASTER_FEATURE_SPEC.items():
        col = detect_column(m, cands, required=False, what=key)
        if col:
            found[key] = col
        else:
            missing.append(key)
    positions = {}
    for col in m.columns:
        hit = _POSITION_COL_RE.match(str(col).strip())
        if hit:
            positions[int(hit.group(1))] = str(col)
    info.update(found=found, missing=missing, position_columns=positions)

    acc = norm_acc_series(m[acc_col])
    keep = acc != ""
    info["n_rows_dropped_empty"] = int((~keep).sum())
    out = m.loc[keep].copy()
    out.index = pd.Index(acc[keep].to_numpy(), name="accession")
    dup = out.index.duplicated(keep=False)
    info["n_duplicate_accession_rows"] = int(dup.sum())
    if dup.any():
        raise ValueError(
            f"master 에 같은 accession 이 여러 행 있다 ({int(dup.sum())}행, 예: "
            f"{sorted(set(out.index[dup].tolist()))[:5]}). reindex 가 행을 불리므로 중단한다.")
    info["n_assignments"] = int(len(out))
    # set_profile 이 어느 컬럼을 썼는지 되짚을 수 있게 표에 붙여 둔다 (attrs 는 참고용이고,
    # set_profile 은 이것이 없어도 MASTER_FEATURE_SPEC 으로 다시 찾는다).
    out.attrs["found"] = dict(found)
    out.attrs["position_columns"] = dict(positions)
    return out, info


def check_position_numbering(position_columns: dict, *,
                             paper_positions=PAPER_POSITIONS, tolerance: int = 3) -> dict:
    """master 의 r<NNN> 컬럼 번호를 논문 번호와 대조한다.

    r612 가 논문 V610 과 2 어긋난다고 보고됐다. 질의가 pdb|1SU6|A 였으므로 PDB 번호가
    섞였을 가능성이 있다. **보존도 계산 전에 사람이 결정해야 하는 사안**이므로 여기서는
    판정만 하고 고치지 않는다.
    """
    found = sorted(int(p) for p in position_columns)
    paper = sorted(int(p) for p in paper_positions)
    exact = [p for p in paper if p in found]
    near, unmatched = [], []
    for p in paper:
        if p in found:
            continue
        cands = [(abs(f - p), f) for f in found if 0 < abs(f - p) <= tolerance
                 and f not in exact]
        if cands:
            d, f = sorted(cands)[0]
            near.append(dict(paper=p, master=f, delta=int(f - p),
                             column=position_columns[f]))
        else:
            unmatched.append(p)
    explained = {d["master"] for d in near} | set(exact)
    extra = [p for p in found if p not in explained]
    msgs = [f"master 위치 컬럼 {len(found)}개: {[position_columns[p] for p in found]}",
            f"논문 번호와 정확히 일치: {exact}"]
    for d in near:
        msgs.append(f"⚠ r{d['master']} 와 논문 {d['paper']} 가 {d['delta']:+d} 어긋난다 "
                    f"(컬럼 '{d['column']}')")
    if unmatched:
        msgs.append(f"⚠ 논문 번호인데 master 에 대응 컬럼이 없다: {unmatched}")
    if extra:
        msgs.append(f"⚠ 논문에 대응이 없는 master 컬럼: {[position_columns[p] for p in extra]}")
    needs = bool(near or unmatched or extra)
    if needs:
        msgs.append("→ UniProt 번호와 PDB(1SU6) 번호가 섞였을 수 있다. "
                    "**보존도 단계 전에 어느 번호 체계로 통일할지 정해야 한다.** "
                    "세트 선택(이 노트북)에는 영향이 없다.")
    return dict(found=found, paper=paper, exact_match=exact, near_miss=near,
                unmatched_paper=unmatched, extra_master=extra,
                needs_decision=needs, messages=msgs)


def classify_dcluster(value) -> str:
    """D-cluster 유형 값을 SF4 / SF2 / other / '' 로 정규화. 표기가 섞여 있어도 센다."""
    v = ("" if value is None else str(value)).strip().lower()
    if not v or v in ("nan", "na", "none", "-", "?"):
        return ""
    if "sf4" in v or "4fe" in v or "fe4" in v or v in ("4", "4fe4s", "sf4?"):
        return "SF4"
    if "sf2" in v or "2fe" in v or "fe2" in v or v in ("2", "2fe2s"):
        return "SF2"
    return "other"


def master_feature_columns(master: pd.DataFrame, keys=("dtype", "cband")) -> dict:
    """master 표에서 부가 컬럼의 실제 이름을 찾는다. read_master_table 이 이미 찾아
    attrs 에 남겨 두면 그것을 쓰고, 없으면 후보 탐지로 다시 찾는다. 못 찾으면 빠진다."""
    found = dict(master.attrs.get("found") or {})
    out = {}
    for key in keys:
        col = found.get(key)
        if not col or col not in master.columns:
            col = detect_column(master, MASTER_FEATURE_SPEC[key], required=False, what=key)
        if col:
            out[key] = col
    return out


def set_profile(members, *, universe, align_df, fine_partition, master=None,
                master_cols: dict | None = None,
                ref_acc: str = CHCODH2_ACC, name: str = "", method: str = "",
                cutoff=None, max_band_terms: int = 6) -> dict:
    """**한 세트를 재는 공통 잣대.** 두 방식을 같은 숫자로 비교하는 것이 이 함수의 목적이다.

      n              서열 수
      n_effective    fine_partition(0.90 분할) 기준 서로 다른 클러스터 수.
                     거의 같은 서열을 여러 개 세지 않은, **보존도의 실제 분모**다.
      redundancy     1 - n_effective / n
      has_reference  ChCODH2 자신이 들어 있는가 (대표만 쓰면 빠지는 일이 있다)
      fident_*       세트 내 ChCODH2 대비 identity 분포
      n_dtype_sf4/sf2  master 에 D-cluster 유형 컬럼이 있으면, 없으면 None
      cband_top      master 에 cband 컬럼이 있으면 분포
    """
    uni = set(universe)
    mem = sorted(set(members))
    if not mem:
        raise ValueError(f"[{name or '세트'}] 가 비었다. cutoff 가 너무 높거나 조인이 틀렸다.")
    outside = [m for m in mem if m not in uni]
    if outside:
        raise ValueError(
            f"[{name or '세트'}] 에 universe 밖 accession 이 {len(outside)}개 있다 "
            f"(예: {outside[:5]}). 조인 키가 어긋났으므로 중단한다.")
    if fine_partition is None:
        raise ValueError("fine_partition(0.90 분할) 이 없으면 n_effective 를 잴 수 없다.")
    missing_fine = [m for m in mem if m not in fine_partition.index]
    if missing_fine:
        raise ValueError(
            f"[{name or '세트'}] 의 {len(missing_fine)}개가 fine_partition "
            f"('{fine_partition.name}') 에 없다 (예: {missing_fine[:5]}). n_effective 의 "
            f"분모가 모호해지므로 중단한다 — 0.90 클러스터 파일이 잘렸는지 확인할 것.")

    groups = fine_partition.reindex(mem).astype(str)
    n = len(mem)
    n_eff = int(groups.nunique())
    rec: dict = dict(
        name=name or "", method=method or "", cutoff=(None if cutoff is None else float(cutoff)),
        n=n, n_effective=n_eff,
        redundancy=float(1.0 - n_eff / n) if n else float("nan"),
        has_reference=bool(ref_acc in set(mem)),
    )

    fid = align_df["fident"].reindex(mem).astype("float64") if align_df is not None else None
    ok = fid.dropna() if fid is not None else None
    rec.update(
        n_fident_missing=int(fid.isna().sum()) if fid is not None else None,
        fident_min=float(ok.min()) if ok is not None and len(ok) else None,
        fident_median=float(ok.median()) if ok is not None and len(ok) else None,
        fident_max=float(ok.max()) if ok is not None and len(ok) else None,
    )

    rec.update(n_master_missing=None, dtype_column=None, n_dtype_sf4=None,
               n_dtype_sf2=None, n_dtype_other=None, n_dtype_missing=None,
               cband_column=None, cband_top=None, n_cband_missing=None, cband_counts=None)
    if master is not None:
        in_master = [m for m in mem if m in master.index]
        rec["n_master_missing"] = int(n - len(in_master))
        sub = master.reindex(in_master)
        cols = dict(master_cols) if master_cols else master_feature_columns(master)
        dcol = cols.get("dtype")
        if dcol and dcol in sub.columns:
            kinds = [classify_dcluster(v) for v in sub[dcol].tolist()]
            rec.update(dtype_column=dcol,
                       n_dtype_sf4=int(sum(1 for k in kinds if k == "SF4")),
                       n_dtype_sf2=int(sum(1 for k in kinds if k == "SF2")),
                       n_dtype_other=int(sum(1 for k in kinds if k == "other")),
                       n_dtype_missing=int(sum(1 for k in kinds if not k)))
        bcol = cols.get("cband")
        if bcol and bcol in sub.columns:
            vals = sub[bcol].astype("object").fillna("").astype(str).str.strip()
            counts = vals[vals != ""].value_counts()
            rec.update(cband_column=bcol,
                       n_cband_missing=int((vals == "").sum()) + int(n - len(in_master)),
                       cband_counts={str(k): int(v) for k, v in counts.items()},
                       cband_top=";".join(f"{k}:{int(v)}"
                                          for k, v in counts.head(max_band_terms).items()))
    rec["members"] = mem
    return rec


def profile_table(profiles, *, sort_by: str = "n_effective") -> pd.DataFrame:
    """set_profile 결과 여러 개 -> 한 표. **n 이 아니라 n_effective 로 줄 세운다.**"""
    rows = [{k: v for k, v in p.items() if k not in ("members", "cband_counts")}
            for p in profiles]
    df = pd.DataFrame(rows)
    if sort_by and sort_by in df.columns:
        df = df.sort_values([sort_by, "n"], ascending=[False, True], kind="mergesort")
    return df.reset_index(drop=True)


def set_overlap(a, b, *, max_examples: int = 8) -> dict:
    """두 세트의 교집합/자카드/각 쪽에만 있는 수."""
    A, B = set(a), set(b)
    inter, union = A & B, A | B
    return dict(
        n_a=len(A), n_b=len(B), n_intersection=len(inter), n_union=len(union),
        jaccard=float(len(inter) / len(union)) if union else float("nan"),
        n_only_a=len(A - B), n_only_b=len(B - A),
        only_a_examples=sorted(A - B)[:max_examples],
        only_b_examples=sorted(B - A)[:max_examples],
    )


def overlap_matrix(sets_a: dict, sets_b: dict, *, value: str = "jaccard") -> pd.DataFrame:
    """행 = sets_a, 열 = sets_b 의 겹침 행렬. value 는 set_overlap 의 키 하나."""
    allowed = ("jaccard", "n_intersection", "n_only_a", "n_only_b", "n_union")
    if value not in allowed:
        raise ValueError(f"value 는 {allowed} 중 하나여야 한다 (받은 값 '{value}').")
    data = {kb: [set_overlap(va, vb)[value] for va in sets_a.values()]
            for kb, vb in sets_b.items()}
    return pd.DataFrame(data, index=pd.Index(list(sets_a), name="세트"))


def cutoff_stability(build_fn, cutoffs, *, ref_acc: str = CHCODH2_ACC,
                     stable_jaccard: float = 0.90, label: str = "") -> pd.DataFrame:
    """기준을 조금 움직였을 때 세트가 얼마나 바뀌는지 — 인접 cutoff 간 자카드.

    직전 분석에서 ChCODH2 의 클러스터 소속이 0.40~0.55 에서 불안정했다(대표가 바뀌며
    그룹을 건너뜀). identity 기반에서도 같은 구간이 흔들리는지 이 표로 본다.
    """
    vals = [float(c) for c in cutoffs]
    if len(vals) < 2:
        raise ValueError(f"cutoff 가 2개 이상이어야 인접 비교가 된다 (받은 값 {vals}).")
    if len(set(vals)) != len(vals):
        raise ValueError(f"cutoff 에 중복이 있다: {vals}")
    order = sorted(vals)
    built = {c: set(build_fn(c)) for c in order}
    rows = []
    for lo, hi in zip(order[:-1], order[1:]):
        ov = set_overlap(built[lo], built[hi], max_examples=5)
        rows.append(dict(
            label=label, cutoff_low=lo, cutoff_high=hi,
            n_low=ov["n_a"], n_high=ov["n_b"], n_common=ov["n_intersection"],
            jaccard=ov["jaccard"], n_only_low=ov["n_only_a"], n_only_high=ov["n_only_b"],
            delta_n=int(ov["n_b"] - ov["n_a"]),
            ref_in_low=bool(ref_acc in built[lo]), ref_in_high=bool(ref_acc in built[hi]),
            is_stable=bool(ov["jaccard"] >= stable_jaccard),
        ))
    return pd.DataFrame(rows)


def pick_by_effective(prof_df: pd.DataFrame, target_effective: int) -> pd.DataFrame:
    """방식마다 '실질 다양성 target 을 채우는 가장 좁은 후보' 를 하나씩 고른다.

    n_effective 가 큰 쪽을 그냥 고르면 **언제나 가장 느슨한 기준이 이긴다** (세트를
    넓히면 다양성은 자동으로 늘기 때문). 그래서 분모를 먼저 고정하고, 같은 분모를
    채우면서 ChCODH2 와 가장 가까운(fident 하한이 높은) 후보를 방식별로 뽑는다.
    목표를 못 채우는 방식은 그 방식의 최대 n_effective 후보로 대체하고 reached=False 로 표시한다.
    """
    rows = []
    for method, grp in prof_df.groupby("method", sort=True):
        ok = grp[grp["n_effective"] >= int(target_effective)]
        if len(ok):
            g = ok.sort_values(["n_effective", "n"], ascending=[True, True], kind="mergesort")
            rec = g.iloc[0].to_dict(); rec["reached"] = True
        else:
            g = grp.sort_values(["n_effective", "n"], ascending=[False, True], kind="mergesort")
            rec = g.iloc[0].to_dict(); rec["reached"] = False
        rec["target_effective"] = int(target_effective)
        rows.append(rec)
    return pd.DataFrame(rows)


def compare_summary(prof_df: pd.DataFrame, *, target_effective: int | None = None,
                    overlap: dict | None = None,
                    stability: pd.DataFrame | None = None) -> list[str]:
    """사람이 읽는 요약. **결론을 내지 않는다** — 후보와 잃는 것을 숫자로 적어 둔다.

    '최선' 을 n_effective 최대로 정의하지 않는다. 그렇게 하면 가장 느슨한 기준이 항상
    이기므로 비교가 되지 않는다. 대신 방식별 **범위**를 적고, target_effective 가 주어지면
    같은 분모를 채우는 가장 좁은 후보끼리 나란히 놓는다.
    """
    if not len(prof_df):
        raise ValueError("프로필 표가 비었다.")
    need = {"method", "name", "n", "n_effective", "redundancy", "has_reference"}
    missing = need - set(prof_df.columns)
    if missing:
        raise AssertionError(f"프로필 표 컬럼 누락: {sorted(missing)}")
    lines: list[str] = []

    for method, grp in prof_df.groupby("method", sort=True):
        fmin = pd.to_numeric(grp.get("fident_min"), errors="coerce")
        rng = ("" if fmin is None or not fmin.notna().any()
               else f", fident 하한 {float(fmin.min()):.3f}~{float(fmin.max()):.3f}")
        no_ref = grp.loc[~grp["has_reference"].astype(bool), "name"].tolist()
        warn = f"  ⚠ ChCODH2 가 빠진 후보 {no_ref}" if no_ref else ""
        lines.append(f"{method} 후보 {len(grp)}개: n {int(grp['n'].min())}~{int(grp['n'].max())}, "
                     f"n_effective {int(grp['n_effective'].min())}~{int(grp['n_effective'].max())}"
                     f"{rng}{warn}")

    if target_effective is not None:
        picks = pick_by_effective(prof_df, target_effective)
        parts = []
        for _, r in picks.iterrows():
            fr = ("" if pd.isna(r.get("fident_min")) or r.get("fident_min") is None
                  else f"/fident 하한 {float(r['fident_min']):.3f}")
            tag = "" if r["reached"] else "(⚠ 목표 미달, 이 방식의 최대치)"
            parts.append(f"{r['method']} → {r['name']} n={int(r['n'])}"
                         f"/eff={int(r['n_effective'])}"
                         f"/redun={float(r['redundancy']):.2f}{fr}{tag}")
        lines.append(f"n_effective>={int(target_effective)} 를 채우는 가장 좁은 후보: "
                     + "  |  ".join(parts)
                     + " — 서열 수가 많은 쪽은 그만큼 중복을, 하한이 낮은 쪽은 "
                       "그만큼 먼 서열을 안고 간다.")

    if overlap:
        lines.append(f"두 후보의 겹침: 공통 {overlap['n_intersection']}개, 자카드 "
                     f"{overlap['jaccard']:.2f} (한쪽에만 {overlap['n_only_a']} / "
                     f"{overlap['n_only_b']}) — 겹치지 않는 쪽이 보존도를 가르는 서열이다.")
    if stability is not None and len(stability):
        worst = stability.sort_values("jaccard", kind="mergesort").iloc[0]
        lines.append(f"기준 민감도: {worst['label']} 가 {float(worst['cutoff_low']):.2f}→"
                     f"{float(worst['cutoff_high']):.2f} 에서 자카드 {float(worst['jaccard']):.2f} "
                     f"로 가장 크게 바뀐다 (n {int(worst['n_low'])}→{int(worst['n_high'])}). "
                     f"그 구간을 기준으로 삼으면 0.05 차이로 세트가 달라진다.")
    lines.append("결론은 비워 둔다 — 넓히면 실질 다양성은 늘지만 다른 CODH 계열이 섞이고, "
                 "좁히면 순수하지만 보존도의 분모가 작아진다. 위 n_effective 와 "
                 "fident 하한을 보고 고를 것.")
    return lines
