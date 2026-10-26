"""ChCODH2 상동체 클러스터 스윕 — 노트북 전용 단일 모듈.

cluster_sweep.ipynb 옆에 두기만 하면 된다. git 도, 다른 파일도 필요 없다.
의존성은 pandas + numpy + 표준 라이브러리뿐.

이 파일은 pipeline/lib/common.py 와 pipeline/stages/02_cluster_sweep.py 에서
노트북이 쓰는 함수만 뽑아 합친 것이다. 원본은 git 레포에 있고 거기서 검증됐다.
직접 고치지 말 것 — 다음 갱신 때 덮어쓴다.

2026-10-25c 에서 cluster_set_compare.ipynb 용 세트 비교 함수를 더했다 (파일 뒤쪽).
2026-10-26a 에서 equiv_map.ipynb 용 좌표계 함수를 더했다 (파일 맨 뒤쪽).
정렬 파일 읽기와 qaln/taln 파싱은 pipeline/stages/03_equiv_map.py 에서 검증된 로직을
옮긴 것이다. 앞쪽 함수(01·02 번 노트북이 쓰는 것)는 전혀 건드리지 않았다.
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

VERSION = "2026-10-26a"   # 노트북 1절이 이 값을 대조한다. 어긋나면 둘 중 하나가 낡은 것이다.

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

def write_params(outdir, stage: str, *, inputs=None, outputs=None, commands=None,
                 filename: str = "params.json", **extra) -> Path:
    """산출물 옆에 재현 정보를 남긴다. 인계문서가 존재하는 이유가 이것이다.

    filename — 한 폴더를 여러 노트북이 공유하면 params.json 이 서로 덮어쓴다.
    노트북마다 params_<태그>.json 으로 나눠 쓸 것.
    """
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
    dest = outdir / filename
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


def assert_alignment_usable(align_df: pd.DataFrame, meta: dict, coverage: dict, *,
                            ref_acc: str = CHCODH2_ACC,
                            min_reference_fident: float = 0.99,
                            max_outside_frac: float = 0.01,
                            max_near_identical_frac: float = 0.01) -> list[str]:
    """정렬 파일이 **정말 ChCODH2 대비 fident 인지** 못 박는다. 어긋나면 경고가 아니라 예외.

    헤더도 --columns 도 없는 파일에서 열이 한 칸 밀리면, fident 자리에 qcov(0.95)나
    qstart(1) 가 들어와도 0~1 범위라 척도 판정을 그냥 통과한다. 그러면 identity 기반
    세트가 전부 universe 전체가 되는데 표만 보면 멀쩡해 보인다 — 가장 나쁜 실패다.
    그래서 easy-search 가 반드시 만족하는 것만 골라 검사한다.
      ① 자기 히트(ChCODH2 vs ChCODH2) 의 fident 는 1.0 이다.
      ② 1,827 서열의 fident 가 전부 같은 값일 수는 없다 (01번: 중앙 0.379, 연속 분포).
      ③ 이 정렬은 그 FASTA 를 DB 로 돌린 것이므로 universe 밖 target 이 나올 수 없다.
         나오면 절단이 아니라 **accession 정규화 불일치**다 (조인이 조용히 비는 경로).
      ④ ChCODH2 와 사실상 동일한(>=0.999) 서열이 여럿일 수는 없다. 자기 히트 1건이 정상.
    ①②④ 가 깨지면 열이 밀린 것이고, ③ 이 깨지면 조인 키가 어긋난 것이다. 어느 쪽이든
    아래 숫자 전부가 틀리므로 여기서 멈춘다. 임계값은 전부 키워드로 열어 두었다.
    """
    lines: list[str] = []
    hint = ("columns=cs.ALIGN_FORMAT_OUTPUT 를 1절 ALIGN_COLUMNS 에 넣어 열 순서를 "
            "직접 주거나, 실제 --format-output 순서를 확인할 것.")

    if not meta.get("reference_in_alignment"):
        raise ValueError(
            f"{ref_acc} 자기 히트가 정렬 결과에 없다. 질의가 ChCODH2 가 아니거나 "
            f"target 열이 어긋났다. {hint}")
    ref_fid = meta.get("reference_fident")
    if ref_fid is None or float(ref_fid) < min_reference_fident:
        raise ValueError(
            f"{ref_acc} 자기 히트의 fident 가 {ref_fid} 다. easy-search 자기 히트는 1.0 "
            f"이어야 하므로 fident 열이 아닌 것을 읽고 있다 "
            f"(열 이름 출처: {meta.get('column_source')}). {hint}")
    lines.append(f"자기 히트 fident {float(ref_fid):.4f} — 1.0 이어야 한다 (통과)")

    fmin, fmax = meta.get("fident_min"), meta.get("fident_max")
    if fmin is None or fmax is None:
        raise ValueError("fident 값이 하나도 없다 — 정렬 파일을 다시 확인할 것.")
    if float(fmax) - float(fmin) < 1e-9:
        raise ValueError(
            f"fident 가 전부 {float(fmin)} 로 같다. 1,827 서열의 ChCODH2 대비 identity 가 "
            f"상수일 수는 없으므로 qcov/tcov/qstart 같은 다른 열을 읽고 있다 "
            f"(열 이름 출처: {meta.get('column_source')}). {hint}")
    lines.append(f"fident 범위 {float(fmin):.4f}~{float(fmax):.4f} — 상수가 아니다 (통과)")

    n_aln = int(coverage.get("n_aligned") or 0)
    n_out = int(coverage.get("n_aligned_outside_universe") or 0)
    out_cap = max(1, int(max_outside_frac * n_aln))
    lines.append(f"universe 밖 target {n_out}/{n_aln}  (허용 {out_cap})"
                 + (f"  예: {coverage.get('outside_examples', [])[:5]}" if n_out else ""))
    if n_out > out_cap:
        raise ValueError(
            f"정렬 결과의 target 중 {n_out}/{n_aln} 개가 FASTA universe 에 없다 "
            f"(예: {coverage.get('outside_examples', [])[:5]}). 이 정렬은 그 FASTA 를 DB 로 "
            f"돌린 것이므로 절단이 아니라 accession 표기가 다른 것이다 (버전 접미사 .1, "
            f"sp|…| 접두사 등). 조인이 조용히 비면 identity 기반 세트가 ChCODH2 하나로 "
            f"쪼그라든다. parse_accession 규칙을 양쪽에 맞춘 뒤 다시 돌릴 것.")

    near = align_df["fident"].astype("float64")
    n_near = int((near >= 0.999).sum())
    near_cap = max(1, int(max_near_identical_frac * n_aln))
    lines.append(f"fident >= 0.999 인 target {n_near}개 — 자기 히트 1건이 정상 (허용 {near_cap})")
    if n_near > near_cap:
        raise ValueError(
            f"ChCODH2 와 사실상 동일한(fident>=0.999) target 이 {n_near}/{n_aln} 개다. "
            f"자기 히트 1건 말고는 나올 수 없으므로 fident 가 아닌 열(tcov/qcov 등) 을 읽고 "
            f"있을 가능성이 크다 (열 이름 출처: {meta.get('column_source')}). 중복 서열이 "
            f"실제로 그만큼 있다면 max_near_identical_frac 를 올려 다시 돌릴 것. {hint}")
    return lines


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
    claimed = set(exact)          # master 컬럼 하나를 두 논문 번호가 나눠 갖지 않게 한다
    for p in paper:
        if p in found:
            continue
        cands = [(abs(f - p), f) for f in found if 0 < abs(f - p) <= tolerance
                 and f not in claimed]
        if cands:
            d, f = sorted(cands)[0]
            claimed.add(f)
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
    넓히면 다양성은 자동으로 늘기 때문). 그래서 분모를 먼저 고정하고, target 을 넘는
    것 중 n_effective 가 가장 작은 후보(동률이면 n 이 작은 쪽)를 방식별로 뽑는다.
    즉 '목표 분모를 겨우 채우는 가장 좁은 후보' 다 — fident 하한으로 고르지는 않는다.
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
        fmin = (pd.to_numeric(grp["fident_min"], errors="coerce")
                if "fident_min" in grp.columns else None)
        rng = ("" if fmin is None or not bool(fmin.notna().any())
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


# ==========================================================================
# ChCODH2 좌표계 — 상동체 잔기를 1~636 번호로 (2026-10-26a 추가)
#
# equiv_map.ipynb 가 쓰는 부분. 파싱·검증 로직은 pipeline/stages/03_equiv_map.py
# (자기검증 79항목 통과) 에서 그대로 옮겼다. **새로 발명한 것은 없다.**
#   - 검증 1 (V1) qaln 갭 제거 == 참조서열[qstart..] 부분서열. 전 행 정확히 일치.
#   - 검증 2 (V2) len(qaln) == len(taln). **HSP 선택 전 전 행**에 적용.
#   - 검증 3 (V3) 자기 히트는 q_pos == t_pos 가 636개 전부 성립.
#   - 검증 4 (V4) base 자동 판별. 0-based 채택은 **파일에 색인 0 이 실제로 있을 때만.**
#   - 검증 5 (V5) qlen 컬럼이 있으면 636 과 대조.
#   - 검증 6 (V6) taln 소문자(소프트 마스킹) 검출 — 보존도만 조용히 과소계산된다.
#   - 검증 7 (V7) universe 축소 검출 — 행렬 행 수가 1,827 아래면 보고.
# 어느 것이든 깨지면 경고가 아니라 예외다. 뒤의 보존도·Fix/Bias 판정이 전부 조용히
# 틀리는 쪽이 훨씬 비싸기 때문이다.
#
# 번호 체계는 **참조 FASTA(pipeline/lib/chcodh2.fasta, Met1=1, 636 aa) 하나뿐**이다.
# master_1827.tsv 의 r* 컬럼은 대조 대상이지 좌표의 출처가 아니다
# (compare_master_positions 가 잔기 값으로 판정한다).
# ==========================================================================

GAP_CHAR = "-"
OUTSIDE_CHAR = "-"          # 정렬 블록 밖. 기본값은 갭과 같은 문자다(행렬 파일을 단순하게
                            # 유지하려고). 구분이 필요한 계산은 커버리지 마스크를 쓴다.
EXPECTED_N_SEQUENCES = 1827         # IPR010047 다운로드 서열 수 (universe)
STD_AA = "ACDEFGHIKLMNPQRSTVWY"

EQUIV_COLUMNS = ["accession", "q_pos", "q_aa", "t_pos", "t_aa"]
EQUIV_REQUIRED_KEYS = ("accession", "qstart", "qaln", "tstart", "taln")

# COLSPEC 에 없는 컬럼의 후보 이름. COLSPEC 은 01·02 가 쓰므로 건드리지 않는다.
EQUIV_EXTRA_COLSPEC = {
    "qend":   ["qend", "q_end"],
    "tend":   ["tend", "t_end", "send"],
    "qlen":   ["qlen", "q_len", "query_len"],
    "tlen":   ["tlen", "t_len", "target_len"],
    "alnlen": ["alnlen", "aln_len", "length", "alignment_length"],
}

# 논문에서 확인된 16개 위치 (pipeline/lib/common.KNOWN_POSITIONS 와 같은 내용).
# aa 는 참조서열에서 **반드시 재확인된다** (verify_known_positions).
KNOWN_POSITIONS = {
    41:  dict(aa="F", role="d_cluster_aromatic", note="F41A 활성 소실; 전자전달 관문"),
    43:  dict(aa="E", role="mediator_electrostatic", note="viologen 정전기 결합; 터널#1 입구"),
    57:  dict(aa="R", role="mediator_affinity", note="R57G/N59L 에서 EV 친화도 10배 개선"),
    59:  dict(aa="N", role="mediator_affinity", note="R57G/N59L 에서 EV 친화도 10배 개선"),
    206: dict(aa="Q", role="tunnel3_entrance", note="비선택 터널 #3 입구"),
    450: dict(aa="K", role="tunnel1_entrance", note="비선택 터널 #1 입구"),
    559: dict(aa="A", role="tunnel_bottleneck", note="A559W/H 로 O2 내성 61~148배; triad"),
    565: dict(aa="V", role="tunnel_bottleneck_triad", note="A559/V565/I580 triad"),
    580: dict(aa="I", role="tunnel_bottleneck_triad", note="A559/V565/I580 triad"),
    583: dict(aa="L", role="tunnel1_entrance", note="비선택 터널 #1 입구"),
    586: dict(aa="I", role="tunnel2_entrance", note="비선택 터널 #2 입구"),
    593: dict(aa="T", role="tunnel2_entrance", note="비선택 터널 #2 입구"),
    597: dict(aa="T", role="tunnel2_entrance", note="비선택 터널 #2 입구"),
    599: dict(aa="S", role="tunnel3_entrance", note="비선택 터널 #3 입구"),
    603: dict(aa="I", role="tunnel3_entrance", note="비선택 터널 #3 입구"),
    610: dict(aa="V", role="tunnel2_bottleneck", note="V610H 로 O2 내성 추가 개선"),
}

COVERAGE_COLUMNS = [
    "q_pos", "q_aa", "is_known_position", "known_role",
    "n_universe", "n_in_alignment", "n_not_in_alignment",
    "n_residue", "n_gap_within_aln", "n_outside_aln",
    "n_identical", "residue_frac_in_alignment", "residue_frac_universe",
    "identity_frac_of_residues",
]


def read_reference_seq(path, *, ref_acc: str = CHCODH2_ACC,
                       expect_len: int = CHCODH2_LEN) -> tuple[str, dict]:
    """참조서열을 읽는다. **이 서열 하나가 좌표계의 유일한 출처다.**

    레코드가 둘 이상이면 어느 것이 질의였는지 알 수 없으므로 예외. 길이가 636 이
    아니면 번호가 전부 밀리므로 예외 (expect_len=None 으로 끌 수 있다).
    """
    recs = list(read_fasta(path))
    if not recs:
        raise ValueError(f"{path} 에서 FASTA 레코드를 읽지 못했다.")
    if len(recs) > 1:
        raise ValueError(
            f"{path} 에 레코드가 {len(recs)}개 있다 (헤더 예: "
            f"{[h.split()[0] for h, _ in recs][:5]}). 어느 것이 search 질의였는지 "
            f"확정할 수 없으므로 중단한다 — 참조서열 1개 FASTA 를 줄 것.")
    header, seq = recs[0]
    seq = "".join(str(seq).split()).upper()
    acc = safe_parse_accession(header)
    bad = sorted({c for c in seq if c not in STD_AA})
    meta = dict(path=str(path), header=header, accession=acc, length=len(seq),
                expected_length=expect_len, nonstandard_chars=bad,
                accession_matches_reference=bool(acc == ref_acc))
    if expect_len is not None and len(seq) != int(expect_len):
        raise ValueError(
            f"참조서열 길이가 {len(seq)} 다 ({expect_len} 이어야 한다). 좌표 전체가 "
            f"밀리므로 중단한다. 파일 {path}, 헤더 '{header}'. 다른 아이소폼이거나 "
            f"태그가 붙은 서열일 수 있다.")
    if bad:
        raise ValueError(
            f"참조서열에 표준 20 아미노산이 아닌 문자가 있다: {bad}. 번호를 믿을 수 "
            f"없으므로 중단한다 (파일 {path}).")
    return seq, meta


def verify_known_positions(seq: str, positions: dict | None = None) -> dict:
    """논문 위치의 잔기가 참조서열과 맞는지 확인한다. 어긋나면 즉시 예외.

    **이것이 '번호 체계가 맞다' 의 유일한 직접 증거다.** 여기가 통과해야 아래의
    모든 위치 번호가 논문과 같은 뜻이 된다.
    """
    positions = positions or KNOWN_POSITIONS
    L = len(seq)
    out_of_range = [p for p in positions if not (1 <= int(p) <= L)]
    if out_of_range:
        raise ValueError(f"참조서열 {L} aa 범위를 벗어난 위치: {sorted(out_of_range)}")
    bad = [(int(p), m["aa"], seq[int(p) - 1]) for p, m in sorted(positions.items())
           if seq[int(p) - 1] != m["aa"]]
    if bad:
        raise ValueError(
            f"번호 체계 불일치 (위치, 논문 잔기, 참조서열 잔기) = {bad}. 참조 FASTA 가 "
            f"논문의 번호와 다른 서열이므로 중단한다.")
    return dict(n_positions=len(positions), positions=sorted(int(p) for p in positions),
                residues={int(p): seq[int(p) - 1] for p in sorted(positions)},
                passed=True)


def _resolve_equiv_columns(df: pd.DataFrame) -> dict:
    """정렬 표의 논리 이름 -> 실제 컬럼명. qaln/taln 이 없으면 좌표를 만들 수 없다."""
    spec = dict(COLSPEC)
    spec.update(EQUIV_EXTRA_COLSPEC)
    out = {}
    for key in ("accession", "query", "fident", "qstart", "qaln", "tstart", "taln",
                "qcov", "tcov", "evalue", "bits", "qend", "tend", "qlen", "tlen", "alnlen"):
        out[key] = detect_column(df, spec[key], required=False, what=key)
    missing = [k for k in EQUIV_REQUIRED_KEYS if not out.get(k)]
    if missing:
        raise KeyError(
            f"정렬 표에 {missing} 컬럼이 없다 (실제 컬럼 {list(df.columns)[:20]}). "
            f"좌표계는 qstart/qaln/tstart/taln 으로만 만들 수 있다. --format-output 에 "
            f"이 네 개가 들어갔는지 확인하고, 헤더가 없으면 columns 인자로 순서를 줄 것: "
            f"columns=cs.ALIGN_FORMAT_OUTPUT")
    return out


def read_alignment_rows(path, columns=None, *, ref_acc: str = CHCODH2_ACC,
                        assume_format_output: bool = True) -> tuple[pd.DataFrame, dict]:
    """정렬 파일을 **행 그대로** 읽는다 (HSP 병합 없음, qaln/taln 보존).

    read_alignment 은 target 당 1행으로 줄이고 qaln/taln 을 버리므로 좌표계에는 쓸 수
    없다. 그래서 같은 헤더 판정 규칙(looks_like_align_header)을 쓰되 모든 행을 남긴다.
    반환 표의 컬럼 이름은 논리 이름으로 정규화된다 (accession/qstart/qaln/tstart/taln/…).

    열 이름을 정하는 세 갈래는 read_alignment 과 같다 — ① columns 인자 ② 파일 헤더
    ③ 복원된 --format-output 순서(추정). 어느 쪽을 썼는지 meta['column_source'] 에 남는다.
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
            f"{first[:20]}.\n어느 열이 qaln/taln 인지 알 수 없으므로 columns 인자로 "
            f"--format-output 순서를 직접 줄 것. 예:\n  columns=cs.ALIGN_FORMAT_OUTPUT")

    df = body.reset_index(drop=True)
    df.columns = list(names)
    if df.empty:
        raise ValueError(f"{path} 에 헤더만 있고 데이터 행이 없다.")
    cmap = _resolve_equiv_columns(df)
    meta.update(columns=list(df.columns), column_map={k: v for k, v in cmap.items() if v})

    if cmap.get("query"):
        qvals = sorted({str(v).strip() for v in df[cmap["query"]].tolist() if str(v).strip()})
        meta.update(n_query_values=len(qvals), query_values=qvals[:5])
        if len(qvals) > 1:
            raise ValueError(
                f"{path} 에 질의가 {len(qvals)}종류 있다 ({qvals[:5]}). 좌표계는 ChCODH2 "
                f"단일 질의를 전제하므로 q_pos 의 뜻이 달라진다. 중단한다.")
    out = pd.DataFrame(index=pd.RangeIndex(len(df)))
    out["accession"] = norm_acc_series(df[cmap["accession"]]).to_numpy()
    for key in ("qaln", "taln"):
        out[key] = df[cmap[key]].astype("object").astype(str).str.strip().to_numpy()
    for key in ("qstart", "tstart", "qend", "tend", "qlen", "tlen", "alnlen"):
        if cmap.get(key):
            out[key] = _to_float(df[cmap[key]]).to_numpy()
    for key in ("fident", "bits", "evalue", "qcov", "tcov"):
        if cmap.get(key):
            out[key] = _to_float(df[cmap[key]]).to_numpy()
    if cmap.get("query"):
        out["query"] = df[cmap["query"]].astype("object").astype(str).str.strip().to_numpy()
    out["row_in_file"] = np.arange(len(out), dtype=np.int64)

    n_empty = int((out["accession"] == "").sum())
    if n_empty:
        raise ValueError(
            f"target accession 이 빈 행이 {n_empty}개 있다. 조인이 조용히 실패하므로 중단한다.")
    for key in ("qstart", "tstart"):
        n_bad = int(out[key].isna().sum())
        if n_bad:
            ex = df.loc[out[key].isna().to_numpy(), cmap[key]].astype(str).head(5).tolist()
            raise ValueError(
                f"{key} 를 정수로 읽을 수 없는 행이 {n_bad}/{len(out)}개 있다 (예 {ex}). "
                f"열 이름을 {meta['column_source']} 로 정했는데 열이 어긋났을 수 있다.")
    if "fident" in out.columns:
        fid = out["fident"]
        fmax = float(fid.max()) if int(fid.notna().sum()) else float("nan")
        if fmax > 1.0:
            if fmax > 100.0 + 1e-6:
                raise ValueError(f"fident 최댓값이 {fmax} 다 — 분수도 퍼센트도 아니다. 열이 어긋났다.")
            out["fident"] = fid / 100.0
        meta.update(fident_scale="percent(100으로 나눴다)" if fmax > 1.0 else "fraction(0~1)",
                    fident_input_max=fmax, n_fident_missing=int(fid.isna().sum()))
    meta.update(n_rows=int(len(out)), n_unique_target=int(out["accession"].nunique()),
                reference_accession=ref_acc,
                reference_in_alignment=bool((out["accession"] == ref_acc).any()))
    return out, meta


def check_aln_lengths(rows: pd.DataFrame, *, max_examples: int = 10) -> dict:
    """V2 — len(qaln) == len(taln). **HSP 선택 전 전 행**에 적용한다.

    선택 규칙이 bits/evalue 를 쓰므로, 깨진 행이 선택될 수도 있다. 그래서 선택 뒤가
    아니라 전 행을 먼저 본다. 쌍이 깨진 정렬에서는 잔기 대응을 만들 수 없다.
    """
    qa = rows["qaln"].astype(str)
    ta = rows["taln"].astype(str)
    empty = (qa.str.len() == 0) | (ta.str.len() == 0)
    if empty.any():
        ex = rows.loc[empty.to_numpy(), "accession"].head(max_examples).tolist()
        raise ValueError(
            f"[검증 2 실패] qaln 또는 taln 이 빈 행이 {int(empty.sum())}개 있다. 예시 {ex}. "
            f"--format-output 에 qaln/taln 이 실제로 포함되었는지 확인할 것.")
    bad = (qa.str.len() != ta.str.len()).to_numpy()
    if bad.any():
        ex = [f"{r.accession}: len(qaln)={len(str(r.qaln))} len(taln)={len(str(r.taln))}"
              for r in rows.loc[bad].head(max_examples).itertuples()]
        raise ValueError(
            f"[검증 2 실패] qaln 과 taln 의 길이가 다른 행이 {int(bad.sum())}/{len(rows)}개 "
            f"있다. 정렬 문자열 쌍이 깨졌으므로 잔기 대응을 만들 수 없다.\n  "
            + "\n  ".join(ex))
    return dict(check="V2", passed=True, n_checked=int(len(rows)),
                aln_len_min=int(qa.str.len().min()), aln_len_max=int(qa.str.len().max()),
                aln_len_median=float(qa.str.len().median()),
                detail=f"len(qaln)==len(taln) 전 행 {len(rows)}행 통과 "
                       f"(정렬 열 수 {int(qa.str.len().min())}~{int(qa.str.len().max())})")


def check_query_length(rows: pd.DataFrame, ref_seq: str, *, max_examples: int = 10) -> dict:
    """V5 — qlen 컬럼이 있으면 참조서열 길이와 대조한다.

    'search 질의가 정말 이 참조서열이었나' 의 가장 싼 직접 증거다. 1잔기 차이는
    0-based 오판으로 위장해 t_pos 를 조용히 +1 밀어버리므로 반드시 먼저 걸러야 한다.
    alnlen 은 정렬 열 수와 대조하되 판본마다 정의가 달라 불일치로만 보고한다.
    """
    L = len(ref_seq)
    info = dict(check="V5", reference_len=L, qlen_present=bool("qlen" in rows.columns),
                alnlen_present=bool("alnlen" in rows.columns), notes=[])
    if "qlen" not in rows.columns:
        info.update(passed=None, detail="qlen 컬럼이 없다 — 질의 서열 길이를 직접 확인할 수 "
                                        "없고 V1·V3 의 간접 증거에만 의존한다 "
                                        "(--format-output 에 qlen 을 넣으면 가장 싸게 확인된다).")
        info["notes"].append(info["detail"])
        return info
    ql = rows["qlen"]
    vals = sorted({int(v) for v in ql.dropna().tolist()})
    bad = (ql.notna() & (ql != L)).to_numpy()
    info.update(qlen_values=vals[:20], n_qlen_unparseable=int(ql.isna().sum()),
                n_qlen_mismatch=int(bad.sum()), n_rows=int(len(rows)))
    if int(ql.isna().sum()):
        info["notes"].append(f"qlen 을 숫자로 읽지 못한 행 {int(ql.isna().sum())}개 — 그 행은 건너뛴다.")
    if bad.any():
        ex = rows.loc[bad, "accession"].head(max_examples).tolist()
        raise ValueError(
            f"[검증 5 실패] 정렬 파일의 qlen 이 참조서열 길이 {L} 와 다른 행이 "
            f"{int(bad.sum())}/{len(rows)}개 있다 (qlen 값 {vals[:10]}). search 질의가 "
            f"참조 FASTA 의 {L} aa 가 아니었다는 직접 증거이며, 그러면 q_pos/t_pos 전체가 "
            f"다른 좌표계가 된다. 예시 accession {ex}\n"
            f"  특히 1잔기 차이는 0-based 오판으로 위장해 t_pos 를 조용히 밀어버린다.\n"
            f"  조치 — 질의 FASTA 를 찾아 참조서열과 대조하고, 다르면 올바른 질의로 "
            f"search 를 다시 돌릴 것.")
    if "alnlen" in rows.columns:
        al = rows["alnlen"]
        cols = rows["qaln"].astype(str).str.len()
        nbad = int((al.notna() & (al != cols)).sum())
        info["n_alnlen_mismatch"] = nbad
        if nbad:
            info["notes"].append(
                f"alnlen 이 len(qaln) 과 다른 행 {nbad}/{len(rows)}개 — 열 순서가 어긋났거나 "
                f"alnlen 정의가 다른 판본이다 (매핑은 qaln/taln 만 쓰므로 좌표에는 영향 없음).")
    info.update(passed=True,
                detail=f"qlen 전 행 {int(ql.notna().sum())}행이 참조서열 {L} aa 와 일치")
    return info


def check_soft_masking(rows: pd.DataFrame, *, upper_aln: bool = False) -> dict:
    """V6 — taln 소문자(소프트 마스킹) 검출.

    소문자는 V1~V3 를 **전부 통과하면서** 행렬에 그대로 들어간다. 그러면 'a' != 'A' 라서
    n_identical 과 하위 단계의 보존도가 조용히 과소계산된다. 그래서 경고가 아니라 예외다.
    """
    info = dict(check="V6", upper_aln=bool(upper_aln))
    counts: dict[str, int] = {}
    for key in ("qaln", "taln"):
        ser = rows[key].astype(str)
        n = int(ser.str.count(r"[a-z]").sum())
        info[f"n_lower_{key}"] = n
        if n:
            chars = sorted({c for s in ser for c in s if c.islower()})
            counts[key] = n
            info[f"lower_chars_{key}"] = chars[:12]
    if counts and not upper_aln:
        raise ValueError(
            f"[검증 6 실패] 정렬 문자열에 소문자 잔기가 있다 ({counts}). 소프트 마스킹된 "
            f"산출물이면 V1~V3 는 전부 통과하지만, 행렬에 소문자가 그대로 들어가면 "
            f"'a' != 'A' 가 되어 **n_identical 과 보존도가 조용히 과소계산된다.**\n"
            f"  조치 — parse_equiv(..., upper_aln=True) 로 다시 돌릴 것 "
            f"(qaln/taln 을 대문자로 올려 비교·매핑한다). 그 사실은 감사 json 에 남는다.")
    info.update(passed=True, n_lower_total=int(sum(counts.values())),
                detail=("소문자 없음 (대문자 정렬)" if not counts else
                        f"소문자 {sum(counts.values())}개를 upper_aln=True 로 대문자화했다 "
                        f"— 보존도 과소계산을 막았다"))
    return info


def select_best_hsp(rows: pd.DataFrame, *, max_examples: int = 5) -> tuple[pd.DataFrame, dict]:
    """target 당 1행으로 줄인다. 규칙은 03_equiv_map 과 같다 —
    bits 최대 → evalue 최소 → qcov 최대 → 파일 등장 순서."""
    work = rows.copy()
    work["n_hsp"] = work.groupby("accession")["row_in_file"].transform("size").to_numpy()
    keys, asc, used = [], [], []
    for key, ascending, label in (("bits", False, "bits 최대"),
                                  ("evalue", True, "evalue 최소"),
                                  ("qcov", False, "qcov 최대")):
        if key in work.columns:
            keys.append(key); asc.append(ascending); used.append(label)
    keys.append("row_in_file"); asc.append(True); used.append("파일 등장 순서")
    sel = (work.sort_values(keys, ascending=asc, kind="mergesort")
               .drop_duplicates("accession", keep="first")
               .sort_values("row_in_file", kind="mergesort")
               .reset_index(drop=True))
    multi = sel.loc[sel["n_hsp"] > 1, "accession"].head(max_examples).tolist()
    meta = dict(n_rows=int(len(rows)), n_selected=int(len(sel)),
                n_merged=int(len(rows) - len(sel)),
                n_targets_multi_hsp=int((sel["n_hsp"] > 1).sum()),
                multi_hsp_examples=multi, tiebreak=" → ".join(used))
    return sel, meta


def detect_qstart_base(sel: pd.DataFrame, ref_seq: str, *, assume_base=None,
                       max_examples: int = 10, upper_aln: bool = False) -> dict:
    """V1 의 전반부 + V4 — qstart 가 1-based 인지 0-based 인지 **데이터로** 판정한다.

    두 해석 모두로 전 행을 맞춰보고 이긴 쪽을 채택하되, **0-based 를 채택하려면 파일에
    색인 0 이 실제로 존재해야 한다.** 없으면 보류하고 예외다 — 질의가 참조와 1잔기
    어긋난 판이 0-based 로 오진되면 q_pos 는 맞게 나오고 V1 도 통과하면서 t_pos 만
    전 행 +1 밀려 **조용히 틀린다.** 가장 위험한 경로이므로 여기서 멈춘다.
    """
    gap = GAP_CHAR
    L = len(ref_seq)
    ref_up = ref_seq.upper()
    stat = {1: 0, 0: 0}
    case_only = {1: 0, 0: 0}
    out_of_range = {1: 0, 0: 0}
    examples: dict[int, list] = {1: [], 0: []}
    qstart_vals: list[int] = []
    tstart_vals: list[int] = []

    for r in sel.itertuples():
        qaln = str(r.qaln).upper() if upper_aln else str(r.qaln)
        qres = qaln.replace(gap, "")
        qstart = int(r.qstart)
        qstart_vals.append(qstart)
        tstart_vals.append(int(r.tstart))
        for base in (1, 0):
            start0 = qstart - 1 if base == 1 else qstart
            if start0 < 0 or start0 + len(qres) > L:
                out_of_range[base] += 1
                if len(examples[base]) < max_examples:
                    examples[base].append(dict(accession=r.accession, qstart=qstart,
                                               n_residues=len(qres), reason="범위 초과"))
                continue
            seg = ref_seq[start0:start0 + len(qres)]
            if qres == seg:
                stat[base] += 1
            elif qres.upper() == ref_up[start0:start0 + len(qres)]:
                case_only[base] += 1
                if len(examples[base]) < max_examples:
                    examples[base].append(dict(accession=r.accession, qstart=qstart,
                                               reason="대소문자만 다름"))
            elif len(examples[base]) < max_examples:
                found = ref_seq.find(qres) if len(qres) >= 8 else -1
                examples[base].append(dict(
                    accession=r.accession, qstart=qstart, n_residues=len(qres),
                    expected_prefix=seg[:30], observed_prefix=qres[:30],
                    found_at_1based=(found + 1) if found >= 0 else None,
                    reason="부분서열 불일치"))

    n = int(len(sel))
    n_qzero = int(sum(1 for v in qstart_vals if v == 0))
    n_tzero = int(sum(1 for v in tstart_vals if v == 0))
    zero = dict(qstart_min=(min(qstart_vals) if qstart_vals else None),
                tstart_min=(min(tstart_vals) if tstart_vals else None),
                n_qstart_zero=n_qzero, n_tstart_zero=n_tzero,
                zero_index_present=bool(n_qzero or n_tzero))

    if assume_base in (0, 1, "0", "1"):
        base = int(assume_base)
        basis = f"assume_base={base} 으로 강제 지정 (데이터 판정을 쓰지 않았다)"
    elif stat[1] == 0 and stat[0] == 0:
        raise ValueError(
            f"[검증 1 실패] 1-based / 0-based 어느 해석으로도 qaln 이 참조서열과 맞는 행이 "
            f"하나도 없다 ({n}행 검사). 질의 서열이 참조서열({L} aa)이 아니거나, "
            f"qaln/taln 컬럼이 뒤바뀌었거나, 참조 FASTA 가 다른 아이소폼이다.\n"
            f"  1-based 범위초과 {out_of_range[1]} · 대소문자만다름 {case_only[1]}\n"
            f"  0-based 범위초과 {out_of_range[0]} · 대소문자만다름 {case_only[0]}\n"
            f"  1-based 예시: {json.dumps(examples[1][:3], ensure_ascii=False)}\n"
            f"  0-based 예시: {json.dumps(examples[0][:3], ensure_ascii=False)}\n"
            f"  (대소문자만 다르면 upper_aln=True 로 대문자 변환 후 진행할 수 있다)")
    elif stat[1] > stat[0]:
        base = 1
        basis = (f"1-based 로 {stat[1]}/{n} 행 일치, 0-based 로 {stat[0]}/{n} 행 일치 → "
                 f"1-based 채택 (mmseqs 기본 가정과 일치)")
    elif stat[0] > stat[1]:
        if not zero["zero_index_present"]:
            raise ValueError(
                f"[검증 4 보류] 0-based 로 {stat[0]}/{n} 행이 맞고 1-based 로 {stat[1]}/{n} 행이 "
                f"맞지만, 파일 어디에도 색인 0 이 없다 (qstart 최소 {zero['qstart_min']} · "
                f"tstart 최소 {zero['tstart_min']} · qstart==0 {n_qzero}행 · "
                f"tstart==0 {n_tzero}행).\n"
                f"  진짜 0-based 산출물이라면 {n}개 국소 정렬 중 어딘가에 0 이 나와야 한다. "
                f"0 이 전혀 없다는 것은 base 가 아니라 **참조서열과 search 질의서열이 1잔기 "
                f"어긋났다**는 뜻일 가능성이 훨씬 높다 (개시 Met 유무·태그·다른 아이소폼).\n"
                f"  그대로 0-based 로 진행하면 q_pos 는 맞게 나오지만 t_pos 가 전부 +1 밀려 "
                f"**조용히 틀린다** (V1 은 통과해버리고, 자기 히트가 없으면 V3 도 못 한다).\n"
                f"  조치 — ① qlen 이 {L} 인지 본다(V5) ② search 질의 FASTA 를 참조서열과 "
                f"대조한다 ③ 진짜 0-based 임을 확인했으면 assume_base=0 으로 강제한다.")
        base = 0
        basis = (f"0-based 로 {stat[0]}/{n} 행 일치, 1-based 로 {stat[1]}/{n} 행 일치 → "
                 f"0-based 채택. mmseqs 의 1-based 가정이 깨졌다 (색인 0 증거: "
                 f"qstart==0 {n_qzero}행 · tstart==0 {n_tzero}행)")
    else:
        base = 1
        basis = (f"두 해석의 일치 행 수가 같다({stat[1]}/{n}). 1-based 를 채택하지만 "
                 f"근거가 약하다")
    return dict(check="V4", base=base, basis=basis, n_checked=n, **zero,
                n_match_1based=stat[1], n_match_0based=stat[0],
                n_case_only_1based=case_only[1], n_case_only_0based=case_only[0],
                n_out_of_range_1based=out_of_range[1], n_out_of_range_0based=out_of_range[0],
                ambiguous=bool(stat[1] == stat[0] and assume_base is None),
                assume_base=assume_base, upper_aln=bool(upper_aln),
                mismatch_examples_1based=examples[1][:max_examples],
                mismatch_examples_0based=examples[0][:max_examples],
                passed=True,
                detail=f"base={base} 채택 — {basis}")


def _build_equiv_arrays(sel: pd.DataFrame, ref_seq: str, base: int, *,
                        upper_aln: bool, max_examples: int) -> dict:
    """정렬 문자열을 ChCODH2 좌표로 투영한다. **행마다 V1 을 다시 확인한다.**

    qaln 에 잔기가 있는 **모든 열**이 결과에 들어간다 (상동체 쪽이 갭인 열도 t_aa='-',
    t_pos=0 으로 남는다). 그래야 '갭'(상동체에 그 잔기가 없다) 과 '정렬 블록 밖'
    (애초에 정렬되지 않았다) 이 구분되고, 보존도의 분모가 정확해진다.
    qaln 이 '-' 인 열(상동체에만 있는 삽입 잔기)은 ChCODH2 번호를 갖지 않으므로 버리고
    개수만 센다.
    """
    gap = GAP_CHAR
    L = len(ref_seq)
    acc_list: list[str] = []
    counts: list[int] = []
    q_pos_parts: list[np.ndarray] = []
    t_pos_parts: list[np.ndarray] = []
    q_aa_parts: list[np.ndarray] = []
    t_aa_parts: list[np.ndarray] = []
    per_acc: list[dict] = []
    mismatches: list[str] = []
    n_scanned = 0
    ins_cols = 0
    del_cols = 0
    end_mismatch = dict(qend=0, tend=0, tlen_exceeded=0)
    nonstd: dict[str, int] = {}
    seen: dict[str, int] = {}

    for r in sel.itertuples():
        n_scanned += 1
        acc = str(r.accession)
        if acc in seen:
            raise ValueError(
                f"{acc} 가 선택 후에도 두 번 이상 있다 (행 {seen[acc]} 와 {int(r.row_in_file)}). "
                f"HSP 선택이 실패했으므로 중단한다 — 같은 위치가 두 번 덮인다.")
        seen[acc] = int(r.row_in_file)
        qaln, taln = str(r.qaln), str(r.taln)
        if upper_aln:
            qaln, taln = qaln.upper(), taln.upper()
        qstart, tstart = int(r.qstart), int(r.tstart)
        qarr = np.array(list(qaln), dtype="<U1")
        tarr = np.array(list(taln), dtype="<U1")
        qnon = qarr != gap
        tnon = tarr != gap
        nq, nt = int(qnon.sum()), int(tnon.sum())

        start0 = qstart - 1 if base == 1 else qstart
        tstart0 = tstart - 1 if base == 1 else tstart
        if start0 < 0 or start0 + nq > L:
            raise ValueError(
                f"[검증 1 실패] {acc}: qstart={qstart} (base={base} → 0-based {start0}) 에서 "
                f"잔기 {nq}개는 참조서열 {L} aa 를 벗어난다. base 판정 또는 입력이 잘못됐다.")
        if tstart0 < 0:
            raise ValueError(f"{acc}: tstart={tstart} 가 base={base} 해석에서 음수 색인이 된다.")

        qres = "".join(qarr[qnon].tolist())
        seg = ref_seq[start0:start0 + nq]
        if qres != seg:
            diff = next((i for i in range(nq) if qres[i] != seg[i]), 0)
            found = ref_seq.find(qres) if nq >= 8 else -1
            mismatches.append(
                f"{acc} (파일 {int(r.row_in_file)}행): qstart={qstart}, 잔기 {nq}개. "
                f"첫 불일치 색인 {diff} — 참조 '{seg[diff]}'(위치 {start0 + diff + 1}) vs "
                f"qaln '{qres[diff]}'. 기대 '{seg[:30]}' vs 관측 '{qres[:30]}'"
                + (f". 이 부분서열은 참조의 {found + 1} 번에서 시작한다 "
                   f"(qstart 와 {found + 1 - (start0 + 1):+d} 어긋남)" if found >= 0 else ""))
            if len(mismatches) >= max_examples:
                break
            continue

        qcols = np.flatnonzero(qnon)
        q_pos = (start0 + np.cumsum(qnon) - 1)[qcols] + 1
        t_run = (tstart0 + np.cumsum(tnon) - 1)[qcols] + 1
        t_pos = np.where(tnon[qcols], t_run, 0)
        acc_list.append(acc)
        counts.append(int(qcols.size))
        q_pos_parts.append(q_pos.astype(np.int32))
        t_pos_parts.append(t_pos.astype(np.int32))
        q_aa_parts.append(qarr[qcols])
        t_aa_parts.append(tarr[qcols])

        n_del = int((qnon & ~tnon).sum())
        n_ins = int((~qnon & tnon).sum())
        del_cols += n_del
        ins_cols += n_ins
        for ch, cnt in zip(*np.unique(tarr[tnon], return_counts=True)):
            if str(ch) not in STD_AA:
                nonstd[str(ch)] = nonstd.get(str(ch), 0) + int(cnt)

        q_from, q_to = start0 + 1, start0 + nq
        t_from, t_to = tstart0 + 1, tstart0 + nt
        for key, col in (("qend", "qend"), ("tend", "tend")):
            val = getattr(r, col, None)
            if val is not None and not pd.isna(val):
                obs = int(val) if base == 1 else int(val) + 1
                if obs != (q_to if key == "qend" else t_to):
                    end_mismatch[key] += 1
        tlen = getattr(r, "tlen", None)
        if tlen is not None and not pd.isna(tlen) and t_to > int(tlen):
            end_mismatch["tlen_exceeded"] += 1
        per_acc.append(dict(accession=acc, q_from=q_from, q_to=q_to, t_from=t_from,
                            t_to=t_to, aln_cols=int(qarr.size), n_q_residues=nq,
                            n_t_residues=nt, n_mapped_positions=int(qcols.size),
                            n_deletion_cols=n_del, n_insertion_cols=n_ins))

    if mismatches:
        stopped = len(mismatches) >= max_examples and n_scanned < len(sel)
        raise ValueError(
            f"[검증 1 실패] qaln 에서 갭을 뺀 문자열이 참조서열의 qstart 부분서열과 다른 행이 "
            f"{len(mismatches)}개 있다 (base={base} 채택). 1-based/0-based 오프셋 오류이거나 "
            f"질의 서열이 참조서열이 아니다. 뒤의 모든 숫자가 조용히 틀리므로 중단한다. "
            f"검사한 행 {n_scanned}/{len(sel)}"
            + (f" — 불일치 {max_examples}건에서 탐색을 멈췄으므로 나머지 "
               f"{len(sel) - n_scanned}행은 검사하지 않았다. 실제 불일치는 더 많을 수 있다"
               if stopped else " (전 행 검사 완료)")
            + ":\n  " + "\n  ".join(mismatches))

    lower = {c: n for c, n in nonstd.items() if c.isalpha() and c.islower()}
    if lower:
        raise ValueError(
            f"[검증 6 실패] taln 에 소문자 잔기가 {sum(lower.values())}개 있다 "
            f"({dict(sorted(lower.items())[:8])}). 행렬에 소문자가 들어가면 보존도가 "
            f"조용히 과소계산된다. upper_aln=True 로 다시 돌릴 것.")
    notes: list[str] = []
    other = {c: n for c, n in nonstd.items() if not (c.isalpha() and c.islower())}
    if other:
        notes.append(
            f"상동체 쪽 비표준 문자 {sum(other.values())}개 ({json.dumps(other, ensure_ascii=False)}). "
            f"행렬에 그대로 들어가며 보존도의 분모(n_residue)에는 세지되 n_identical 에는 "
            f"들어가지 않는다. 04 번에서 X/*/U 를 어떻게 다룰지 정할 것.")
    for key, cnt in end_mismatch.items():
        if cnt:
            notes.append(f"{key} 가 매핑 결과와 어긋난 행 {cnt}개 (좌표는 qaln/taln 으로만 "
                         f"만들므로 영향 없음. 열 순서나 판본 차이를 의심할 것).")

    empty = np.array([], dtype=np.int32)
    equiv = pd.DataFrame({
        "accession": np.repeat(np.array(acc_list, dtype=object), counts) if acc_list
                     else np.array([], dtype=object),
        "q_pos": np.concatenate(q_pos_parts) if q_pos_parts else empty,
        "q_aa": np.concatenate(q_aa_parts) if q_aa_parts else np.array([], dtype="<U1"),
        "t_pos": np.concatenate(t_pos_parts) if t_pos_parts else empty,
        "t_aa": np.concatenate(t_aa_parts) if t_aa_parts else np.array([], dtype="<U1"),
    })[EQUIV_COLUMNS]
    return dict(equiv=equiv, per_acc=pd.DataFrame(per_acc), notes=notes,
                n_accessions=len(acc_list), insertion_cols=ins_cols, deletion_cols=del_cols,
                end_mismatch=end_mismatch, nonstandard_target_chars=nonstd)


def check_self_hit(equiv: pd.DataFrame, ref_seq: str, *, ref_acc: str = CHCODH2_ACC) -> dict:
    """V3 — 자기 자신 히트는 q_pos == t_pos 가 636개 전부 성립해야 한다.

    좌표계를 **독립적으로** 확인하는 유일한 그물이다. 자기 히트가 없으면 생략되고,
    그 사실 자체가 보고되어야 한다 (V1 하나만 남으므로 base 오판에 취약해진다).
    """
    info = dict(check="V3", reference_acc=ref_acc, present=False, expected=len(ref_seq))
    sub = equiv[equiv["accession"].to_numpy() == ref_acc]
    if not len(sub):
        info.update(passed=None, skipped=True,
                    detail=f"자기 히트({ref_acc})가 정렬에 없다 — 생략. V1 만 남으므로 "
                           f"base 오판을 교차 확인할 수단이 하나 줄었다.")
        return info
    qpos = sub["q_pos"].to_numpy()
    tpos = sub["t_pos"].to_numpy()
    qaa = sub["q_aa"].astype(str).to_numpy()
    taa = sub["t_aa"].astype(str).to_numpy()
    n_same = int((qpos == tpos).sum())
    n_aa = int((qaa == taa).sum())
    ok = (len(sub) == len(ref_seq) and n_same == len(ref_seq) and n_aa == len(ref_seq)
          and "".join(taa.tolist()) == ref_seq)
    info.update(present=True, skipped=False, n_pairs=int(len(sub)), n_q_eq_t=n_same,
                n_same_residue=n_aa,
                residues_equal_reference=bool("".join(taa.tolist()) == ref_seq), passed=bool(ok))
    if not ok:
        raise AssertionError(
            f"[검증 3 실패] 자기 자신 히트({ref_acc})가 완전 일치가 아니다: 대응 "
            f"{len(sub)}/{len(ref_seq)} · q_pos==t_pos {n_same}/{len(ref_seq)} · "
            f"q_aa==t_aa {n_aa}/{len(ref_seq)}. 좌표계가 어긋났거나 자기 히트가 국소 "
            f"정렬로 잘렸다. 자기 히트가 실제로 부분 정렬이라면 그 사실을 보고서에 적고 "
            f"이 검증을 사람이 판단할 것.")
    info["detail"] = (f"자기 히트 {ref_acc}: q_pos==t_pos {n_same}/{len(ref_seq)}, "
                      f"잔기도 {n_aa}/{len(ref_seq)} 일치 — 좌표계가 참조서열과 같다")
    return info


def parse_equiv(align_df: pd.DataFrame, ref_seq: str, *, assume_base=None,
                ref_acc: str = CHCODH2_ACC, upper_aln: bool = False,
                max_examples: int = 10) -> tuple[pd.DataFrame, dict]:
    """qstart/qaln 과 tstart/taln 으로 잔기 대응을 만든다. **이 파이프라인의 척추다.**

    align_df — read_alignment_rows 의 결과 (HSP 병합 전, 행 그대로).
    반환 (equiv, meta):
      equiv — accession, q_pos(1~636), q_aa(참조 잔기), t_pos(상동체 번호, 갭이면 0),
              t_aa(상동체 잔기, 갭이면 '-'). qaln 에 잔기가 있는 모든 열이 들어간다.
      meta  — checks 에 V1~V6 의 통과 여부와 근거가 순서대로 담긴다. base 를 무엇으로
              채택했고 색인 0 증거가 무엇인지 반드시 들어 있다.
    검증 순서가 중요하다 — V2(전 행) → V5 → V6 → HSP 선택 → V4/V1 → 매핑(V1 재확인) → V3.
    """
    if "qaln" not in align_df.columns or "taln" not in align_df.columns:
        raise KeyError(
            f"align_df 에 qaln/taln 이 없다 (컬럼 {list(align_df.columns)}). "
            f"read_alignment 가 아니라 **read_alignment_rows** 의 결과를 줄 것 "
            f"(read_alignment 은 target 당 1행으로 줄이고 정렬 문자열을 버린다).")
    checks: list[dict] = []
    checks.append(check_aln_lengths(align_df, max_examples=max_examples))
    checks.append(check_query_length(align_df, ref_seq, max_examples=max_examples))
    checks.append(check_soft_masking(align_df, upper_aln=upper_aln))
    sel, hsp = select_best_hsp(align_df)
    base_info = detect_qstart_base(sel, ref_seq, assume_base=assume_base,
                                   max_examples=max_examples, upper_aln=upper_aln)
    built = _build_equiv_arrays(sel, ref_seq, base_info["base"],
                                upper_aln=upper_aln, max_examples=max_examples)
    equiv = built["equiv"]
    n_bad_q = int((equiv["q_aa"].astype(str).to_numpy()
                   != np.array(list(ref_seq), dtype="<U1")[equiv["q_pos"].to_numpy() - 1]).sum())
    if n_bad_q:
        raise AssertionError(
            f"[검증 1 실패] 투영된 q_aa 가 참조서열과 다른 대응이 {n_bad_q}개 있다. "
            f"행별 검사를 통과했는데 전체 대조에서 어긋났다면 투영 색인 계산이 틀린 것이다.")
    checks.insert(0, dict(check="V1", passed=True, base=base_info["base"],
                          n_rows=int(len(sel)), n_pairs=int(len(equiv)),
                          detail=f"qaln 갭 제거 == 참조서열[qstart..] 부분서열: 선택된 "
                                 f"{len(sel)}행 전부 통과, 투영된 대응 {len(equiv):,}개의 "
                                 f"q_aa 도 전부 참조서열과 일치 (base={base_info['base']})"))
    checks.append(base_info)
    checks.append(check_self_hit(equiv, ref_seq, ref_acc=ref_acc))

    n_gap = int((equiv["t_aa"].astype(str).to_numpy() == GAP_CHAR).sum())
    meta = dict(
        reference_len=len(ref_seq), reference_acc=ref_acc, base=base_info["base"],
        base_basis=base_info["basis"], assume_base=assume_base, upper_aln=bool(upper_aln),
        n_input_rows=int(len(align_df)), n_selected_rows=int(len(sel)),
        n_accessions=int(built["n_accessions"]), n_pairs=int(len(equiv)),
        n_gap_pairs=n_gap, n_residue_pairs=int(len(equiv)) - n_gap,
        insertion_cols=int(built["insertion_cols"]), deletion_cols=int(built["deletion_cols"]),
        end_mismatch=built["end_mismatch"],
        nonstandard_target_chars=built["nonstandard_target_chars"],
        hsp=hsp, checks=checks, notes=list(built["notes"]),
        per_accession=built["per_acc"],
    )
    return equiv, meta


def equiv_check_lines(meta: dict) -> list[str]:
    """검증 결과를 한 줄씩. 노트북 3절이 **각 검증의 통과 여부를 명시적으로** 찍는다."""
    order = {"V1": "qaln↔참조 부분서열 일치", "V2": "len(qaln)==len(taln)",
             "V3": "자기 히트 q_pos==t_pos", "V4": "qstart base 판정",
             "V5": "qlen == 참조 길이", "V6": "taln 소문자(소프트 마스킹)"}
    lines = []
    for c in meta.get("checks", []):
        name = order.get(c.get("check"), c.get("check", "?"))
        mark = {True: "통과", False: "실패", None: "생략"}.get(c.get("passed"), "?")
        lines.append(f"  [{c.get('check')}] {mark:<4} {name} — {c.get('detail', '')}")
    for note in meta.get("notes", []):
        lines.append(f"  [참고] {note}")
    return lines


def build_residue_matrix(equiv: pd.DataFrame, universe, ref_seq: str, *,
                         outside_char: str = OUTSIDE_CHAR,
                         expect_n_rows: int = EXPECTED_N_SEQUENCES) -> pd.DataFrame:
    """행=accession(universe 전부), 열=1..636, 값=그 자리에 정렬된 상동체 잔기.

    **universe 는 줄지 않는다** — 정렬에 없는 서열도 행으로 남고 전부 outside_char 다
    (규칙 3: 1,827행은 끝까지 1,827행). 갭('-')과 정렬 블록 밖은 기본값에서 같은 문자가
    되므로, 그 둘을 구분해야 하는 계산은 attrs 의 커버리지 마스크를 쓴다
    (position_coverage 가 그렇게 한다).
    """
    rows = sorted(set(str(a) for a in universe if str(a)))
    if not rows:
        raise ValueError("universe 가 비었다 — FASTA 를 먼저 읽을 것.")
    L = len(ref_seq)
    index = pd.Index(rows, name="accession")
    acc = equiv["accession"].astype(str).to_numpy()
    ridx = index.get_indexer(acc)
    if len(ridx) and int((ridx < 0).sum()):
        missing = sorted(set(acc[ridx < 0].tolist()))
        raise ValueError(
            f"정렬에는 있고 universe 에는 없는 accession 이 {len(missing)}개 있다 "
            f"(예 {missing[:5]}). 행렬에 넣을 자리가 없으므로 서열이 조용히 사라진다. "
            f"accession 정규화 규칙(parse_accession)을 양쪽에 맞춘 뒤 다시 돌릴 것.")
    qpos = equiv["q_pos"].to_numpy()
    if len(qpos) and (int(qpos.min()) < 1 or int(qpos.max()) > L):
        raise ValueError(f"q_pos 가 1~{L} 범위를 벗어난다 (min {int(qpos.min())}, "
                         f"max {int(qpos.max())}). 좌표계가 참조서열과 다르다.")
    mat = np.full((len(rows), L), str(outside_char), dtype="<U1")
    covered = np.zeros((len(rows), L), dtype=bool)
    if len(ridx):
        mat[ridx, qpos - 1] = equiv["t_aa"].to_numpy(dtype="<U1")
        covered[ridx, qpos - 1] = True
    in_aln = covered.any(axis=1)

    out = pd.DataFrame(mat, index=index, columns=pd.RangeIndex(1, L + 1, name="q_pos"))
    notes: list[str] = []
    shrunk = bool(len(rows) < int(expect_n_rows))
    if shrunk:
        notes.append(
            f"[검증 7] universe 축소 — 행렬 행이 {len(rows)}개로 기대 {expect_n_rows} 보다 "
            f"적다. 어딘가에서 서열이 걸러졌다는 뜻이므로 보존도의 분모가 달라진다. "
            f"FASTA 와 master 의 행 수를 다시 확인할 것.")
    n_rows_no_aln = int((~in_aln).sum())
    if n_rows_no_aln:
        notes.append(f"정렬에 전혀 들어가지 않은 서열 {n_rows_no_aln}개 — 행은 남아 있고 "
                     f"값은 전부 '{outside_char}' 다 (n_not_in_alignment 로 센다).")
    out.attrs["meta"] = dict(
        check="V7", n_rows=len(rows), n_position_columns=L, expect_n_rows=int(expect_n_rows),
        universe_shrunk=shrunk, passed=(not shrunk),
        n_in_alignment=int(in_aln.sum()), n_not_in_alignment=n_rows_no_aln,
        n_covered_cells=int(covered.sum()), gap_char=GAP_CHAR, outside_char=str(outside_char),
        n_gap_cells=int((covered & (mat == GAP_CHAR)).sum()),
        memory_bytes=int(mat.nbytes + covered.nbytes), notes=notes,
        detail=(f"행렬 {len(rows)} x {L}, 정렬에 들어간 서열 {int(in_aln.sum())}개, "
                f"덮인 칸 {int(covered.sum()):,}개" + (" — 축소 감지" if shrunk else "")))
    out.attrs["covered"] = covered
    out.attrs["in_alignment"] = in_aln
    out.attrs["reference_seq"] = ref_seq
    return out


def _matrix_masks(matrix: pd.DataFrame, ref_seq: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """행렬에서 (값, 커버리지 마스크, 정렬 포함 여부) 를 꺼낸다.

    attrs 가 사라진 행렬(파일에서 다시 읽은 것 등)은 갭과 '정렬 블록 밖' 을 구분할 수
    없다. 그 상태로 커버리지를 계산하면 분모가 조용히 틀리므로 예외다.
    """
    vals = matrix.to_numpy(dtype="<U1")
    if vals.shape[1] != len(ref_seq):
        raise ValueError(f"행렬 열이 {vals.shape[1]}개인데 참조서열은 {len(ref_seq)} aa 다.")
    covered = matrix.attrs.get("covered")
    in_aln = matrix.attrs.get("in_alignment")
    if covered is None or in_aln is None:
        meta = matrix.attrs.get("meta") or {}
        outside = str(meta.get("outside_char", OUTSIDE_CHAR))
        if outside == GAP_CHAR:
            raise ValueError(
                "행렬에 커버리지 마스크(attrs['covered'])가 없다. 갭('-')과 정렬 블록 밖이 "
                "같은 문자이므로 n_gap_within_aln 과 n_outside_aln 을 구분할 수 없고, "
                "그대로 계산하면 보존도의 분모가 조용히 틀린다. build_residue_matrix 가 "
                "돌려준 행렬을 그대로 넘기거나(권장), 행렬을 다시 만들 것 "
                "(tsv 에서 읽은 행렬은 마스크를 잃는다).")
        covered = vals != outside
        in_aln = covered.any(axis=1)
    return vals, np.asarray(covered), np.asarray(in_aln)


def position_coverage(matrix: pd.DataFrame, ref_seq: str, *, members=None) -> pd.DataFrame:
    """636행. **보존도의 분모다.** 위치마다 몇 개가 잔기이고 갭이고 정렬 밖인지 센다.

    members 를 주면 그 부분집합만 센다 (티어별 커버리지가 이것을 쓴다).
    n_residue + n_gap_within_aln + n_outside_aln + n_not_in_alignment == 행 수 가
    모든 위치에서 성립해야 하며, 깨지면 마스크 계산이 틀린 것이므로 예외다.
    """
    vals, covered, in_aln = _matrix_masks(matrix, ref_seq)
    if members is not None:
        keep = matrix.index.isin(list(members))
        n_req = len({str(a) for a in members})
        if int(keep.sum()) != n_req:
            missing = sorted({str(a) for a in members} - set(matrix.index[keep].tolist()))
            raise ValueError(
                f"members 중 {len(missing)}개가 행렬에 없다 (예 {missing[:5]}). 분모가 "
                f"조용히 작아지므로 중단한다 — 같은 universe 에서 고른 세트인지 확인할 것.")
        vals, covered, in_aln = vals[keep], covered[keep], in_aln[keep]
    n_rows = int(vals.shape[0])
    n_in = int(in_aln.sum())
    n_not_in = n_rows - n_in
    residue = covered & (vals != GAP_CHAR)
    n_residue = residue.sum(axis=0).astype(np.int64)
    n_gap = (covered & (vals == GAP_CHAR)).sum(axis=0).astype(np.int64)
    n_outside = (in_aln[:, None] & ~covered).sum(axis=0).astype(np.int64)
    ref_arr = np.array(list(ref_seq), dtype="<U1")
    n_ident = (residue & (vals == ref_arr[None, :])).sum(axis=0).astype(np.int64)

    total = n_residue + n_gap + n_outside + n_not_in
    if not np.all(total == n_rows):
        bad = int(np.flatnonzero(total != n_rows)[0]) + 1
        raise AssertionError(
            f"position_coverage 분해 항등식 위반: 위치 {bad} 에서 {int(total[bad - 1])} != "
            f"행 수 {n_rows}. 커버리지 마스크 계산이 틀렸다.")
    pos = np.arange(1, len(ref_seq) + 1)
    nan = np.full(len(ref_seq), np.nan)
    df = pd.DataFrame({
        "q_pos": pos, "q_aa": ref_arr,
        "is_known_position": [int(p) in KNOWN_POSITIONS for p in pos],
        "known_role": [KNOWN_POSITIONS.get(int(p), {}).get("role", "") for p in pos],
        "n_universe": n_rows, "n_in_alignment": n_in, "n_not_in_alignment": n_not_in,
        "n_residue": n_residue, "n_gap_within_aln": n_gap, "n_outside_aln": n_outside,
        "n_identical": n_ident,
        "residue_frac_in_alignment": np.divide(n_residue, n_in, out=nan.copy(), where=(n_in > 0)),
        "residue_frac_universe": (n_residue / n_rows if n_rows else nan.copy()),
        "identity_frac_of_residues": np.divide(n_ident, n_residue, out=nan.copy(),
                                               where=(n_residue > 0)),
    })[COVERAGE_COLUMNS]
    if len(df) != len(ref_seq):
        raise AssertionError(f"position_coverage 행 수가 {len(df)} — {len(ref_seq)} 이어야 한다")
    return df


def coverage_summary(cov: pd.DataFrame, *, worst: int = 10,
                     low_frac: float = 0.5) -> dict:
    """커버리지 표를 사람이 읽는 숫자로. **보존도를 믿을 수 없는 위치를 미리 안다.**"""
    n_res = cov["n_residue"].to_numpy()
    frac = cov["residue_frac_universe"].to_numpy()
    low = cov.loc[frac < low_frac]
    known = cov[cov["is_known_position"]]
    order = cov.sort_values(["n_residue", "q_pos"], kind="mergesort").head(worst)
    return dict(
        n_positions=int(len(cov)), n_rows=int(cov["n_universe"].iloc[0]),
        n_in_alignment=int(cov["n_in_alignment"].iloc[0]),
        min_n_residue=int(n_res.min()), max_n_residue=int(n_res.max()),
        median_n_residue=float(np.median(n_res)),
        n_positions_zero_residue=int((n_res == 0).sum()),
        n_positions_below_low_frac=int(len(low)), low_frac=float(low_frac),
        low_runs_1based=[(a + 1, b + 1) for a, b in
                         consecutive_runs([bool(f < low_frac) for f in frac])][:20],
        worst_positions=[dict(q_pos=int(r.q_pos), q_aa=str(r.q_aa), n_residue=int(r.n_residue),
                              residue_frac=float(r.residue_frac_universe),
                              is_known=bool(r.is_known_position))
                         for r in order.itertuples()],
        known_min_n_residue=(int(known["n_residue"].min()) if len(known) else None),
        known_median_n_residue=(float(known["n_residue"].median()) if len(known) else None),
        known_worst=[dict(q_pos=int(r.q_pos), q_aa=str(r.q_aa), n_residue=int(r.n_residue),
                          role=str(r.known_role))
                     for r in known.sort_values("n_residue", kind="mergesort").head(5).itertuples()],
    )


def _equiv_at_positions(equiv: pd.DataFrame, positions, *, universe=None) -> pd.DataFrame:
    """equiv 에서 주어진 위치들만 뽑아 행=accession, 열=위치 로 돌린다.

    값: 상동체 잔기 / '-'(정렬 안에서의 갭) / ''(정렬 블록 밖 또는 정렬에 없음).
    세 가지를 구분해서 돌려주는 것이 요점이다 — 보존도의 분모가 거기서 갈린다.
    """
    pos = [int(p) for p in sorted(set(int(p) for p in positions))]
    if universe is None:
        rows = sorted(set(equiv["accession"].astype(str).tolist()))
    else:
        rows = sorted(set(str(a) for a in universe if str(a)))
    index = pd.Index(rows, name="accession")
    out = np.full((len(rows), len(pos)), "", dtype="<U1")
    qp = equiv["q_pos"].to_numpy()
    keep = np.isin(qp, np.array(pos, dtype=qp.dtype))
    if int(keep.sum()):
        sub_acc = equiv["accession"].astype(str).to_numpy()[keep]
        ridx = index.get_indexer(sub_acc)
        if int((ridx < 0).sum()):
            missing = sorted(set(sub_acc[ridx < 0].tolist()))
            raise ValueError(f"universe 에 없는 accession {len(missing)}개 (예 {missing[:5]}) — "
                             f"행렬과 같은 universe 를 줄 것.")
        lut = np.full(int(pos[-1]) + 1, -1, dtype=np.int64)   # 위치 -> 열 색인. 636개
        lut[np.array(pos, dtype=np.int64)] = np.arange(len(pos), dtype=np.int64)
        cidx = lut[qp[keep].astype(np.int64)]
        if int((cidx < 0).sum()):          # np.isin 과 lut 이 어긋나면 열이 밀린다
            raise AssertionError("위치->열 색인이 어긋났다 (_equiv_at_positions).")
        out[ridx, cidx] = equiv["t_aa"].to_numpy(dtype="<U1")[keep]
    return pd.DataFrame(out, index=index, columns=pos)


def known_position_equivalents(equiv: pd.DataFrame, positions=None, *, universe=None,
                               ref_seq: str | None = None) -> pd.DataFrame:
    """행=accession, 열=논문 위치. 각 상동체의 **그 자리 잔기**를 그대로 담는다.

    열 이름은 'pos41_F' 처럼 참조 잔기를 달고 나온다 (번호만 보면 무슨 자리인지
    잊어버리므로). ref_seq 를 주면 논문 잔기와 참조서열을 대조해 어긋나면 예외다.
    값의 뜻 — 잔기: 그 상동체의 잔기 / '-': 정렬 안에서 갭 / '': 정렬 블록 밖.
    """
    positions = sorted(int(p) for p in (positions or sorted(KNOWN_POSITIONS)))
    if ref_seq is not None:
        bad = [(p, KNOWN_POSITIONS[p]["aa"], ref_seq[p - 1]) for p in positions
               if p in KNOWN_POSITIONS and ref_seq[p - 1] != KNOWN_POSITIONS[p]["aa"]]
        if bad:
            raise AssertionError(f"참조 잔기와 KNOWN_POSITIONS 가 어긋난 위치: {bad}")
    wide = _equiv_at_positions(equiv, positions, universe=universe)
    aa = {p: (ref_seq[p - 1] if ref_seq is not None
              else KNOWN_POSITIONS.get(p, {}).get("aa", "?")) for p in positions}
    out = wide.rename(columns={p: f"pos{p}_{aa[p]}" for p in positions})
    vals = wide.to_numpy(dtype="<U1")
    out["n_residue_covered"] = ((vals != "") & (vals != GAP_CHAR)).sum(axis=1).astype(np.int32)
    out["n_gap"] = (vals == GAP_CHAR).sum(axis=1).astype(np.int32)
    out["n_outside"] = (vals == "").sum(axis=1).astype(np.int32)
    return out


def known_position_frequency(known: pd.DataFrame, ref_seq: str, *, top: int = 5) -> pd.DataFrame:
    """논문 위치마다 **어떤 잔기가 몇 번 나오는지.** F41 이 정말 보존되는지,
    A559 자리가 다른 CODH 에선 무엇인지 여기서 눈으로 확인된다."""
    rows = []
    for col in known.columns:
        m = re.fullmatch(r"pos(\d+)_(.)", str(col))
        if not m:
            continue
        p = int(m.group(1))
        ref_aa = ref_seq[p - 1]
        if m.group(2) != ref_aa:
            raise AssertionError(f"열 '{col}' 의 잔기가 참조서열 {p} 번({ref_aa}) 과 다르다.")
        ser = known[col].astype(str)
        n_out = int((ser == "").sum())
        n_gap = int((ser == GAP_CHAR).sum())
        res = ser[(ser != "") & (ser != GAP_CHAR)]
        vc = res.value_counts()
        n_res = int(len(res))
        top_s = "; ".join(f"{aa} {int(n)} ({n / n_res:.2f})" for aa, n in vc.head(top).items()) \
            if n_res else "(없음)"
        rows.append(dict(
            q_pos=p, q_aa=ref_aa, known_role=KNOWN_POSITIONS.get(p, {}).get("role", ""),
            n_residue=n_res, n_gap=n_gap, n_outside=n_out,
            n_identical=int(vc.get(ref_aa, 0)),
            identity_frac=(float(vc.get(ref_aa, 0)) / n_res if n_res else float("nan")),
            n_distinct_residue=int(len(vc)), top_residues=top_s))
    return pd.DataFrame(rows).sort_values("q_pos", kind="mergesort").reset_index(drop=True)


def _normalize_residue_series(ser: pd.Series) -> pd.Series:
    """master 의 잔기 컬럼을 한 글자 대문자로. 잔기가 아닌 값은 빈 문자열."""
    s = ser.astype("object").fillna("").astype(str).str.strip().str.upper()
    return s.where(s.str.fullmatch(r"[A-Z]"), "")


def compare_master_positions(master: pd.DataFrame, equiv: pd.DataFrame,
                             paper_positions=PAPER_POSITIONS, *, ref_seq: str,
                             window: int = 4, ref_acc: str = CHCODH2_ACC,
                             min_rows: int = 20, min_frac: float = 0.5,
                             margin: float = 0.10) -> dict:
    """master 의 r* 컬럼이 **어느 위치 번호인지 잔기 값으로 판정한다.**

    컬럼 이름(r612)이 아니라, 그 컬럼의 값이 이 좌표계의 어느 위치 잔기와 가장 많이
    일치하는지로 정한다. r612 가 논문 V610 과 +2 어긋난다고 보고됐고 질의가 pdb|1SU6|A
    였으므로 PDB 번호가 섞였을 수 있다 — 그 가설을 **숫자로** 검정하는 것이 이 함수다.
    좌표 자체는 고치지 않는다(참조서열 기준을 쓴다). 판정만 한다.

    후보는 **참조서열 전체 1..636 을 훑는다.** 이름 ±window 만 보면 어긋남이 그 창보다
    크면(예: PDB 번호가 10 어긋남) 창 안의 엉뚱한 위치가 1등이 되어 **확신에 찬 오답**이
    나온다. 실제로 참조서열에는 같은 잔기가 널려 있어(600번과 610번이 모두 V) 창 안의
    이웃이 0.38 쯤 맞는 일이 생긴다. window 는 표에 보여 줄 근접 후보를 고르는 데만 쓴다.

    판정은 **결정적일 때만** 한다 — 1등 일치율이 min_frac 이상이고 2등(인접 위치 제외)과
    margin 이상 벌어져야 한다. 아니면 '판정 불가' 로 남긴다. 애매한 것을 숫자처럼
    적어 두는 쪽이 조용히 틀리는 길이다.
    """
    L = len(ref_seq)
    position_columns = dict(master.attrs.get("position_columns") or {})
    if not position_columns:
        for col in master.columns:
            hit = _POSITION_COL_RE.match(str(col).strip())
            if hit:
                position_columns[int(hit.group(1))] = str(col)
    name_based = check_position_numbering(position_columns,
                                          paper_positions=paper_positions)
    out = dict(position_columns={int(k): str(v) for k, v in position_columns.items()},
               paper_positions=[int(p) for p in paper_positions], window=int(window),
               min_rows=int(min_rows), min_frac=float(min_frac), margin=float(margin),
               scanned_positions=int(L),
               name_based=name_based, columns={}, messages=[],
               n_decided=0, n_shifted=0, n_undecidable=0)
    if not position_columns:
        out["messages"].append("master 에 r<NNN> 꼴 위치 컬럼이 없다 — 대조할 것이 없다.")
        return out

    wide = _equiv_at_positions(equiv, range(1, L + 1))     # 전 위치를 훑는다 (n x 636)
    common = master.index.intersection(wide.index)
    out["n_master_rows"] = int(len(master))
    out["n_common_accession"] = int(len(common))
    if not len(common):
        out["messages"].append(
            "master 와 정렬의 accession 이 하나도 겹치지 않는다 — 조인 키가 다르다. "
            "대조를 할 수 없다.")
        return out

    ours_all = wide.loc[common].to_numpy(dtype="<U1")            # (n_common, L)
    covered_all = (ours_all != "") & (ours_all != GAP_CHAR)
    ref_arr = np.array(list(ref_seq), dtype="<U1")
    ref_row = int(common.get_indexer([ref_acc])[0]) if ref_acc in common else -1

    def _cand(p, named, n_cmp, n_eq, mref):
        """후보 한 칸. p 는 1-based 참조 위치."""
        return dict(position=int(p), ref_aa=str(ref_arr[p - 1]), delta=int(p - int(named)),
                    n_compared=int(n_cmp), n_equal=int(n_eq),
                    frac=(float(n_eq / n_cmp) if n_cmp else float("nan")),
                    reference_row_equal=(None if mref is None else
                                         bool(mref != "" and mref == str(ref_arr[p - 1]))))

    pos_all = np.arange(1, L + 1)
    for n in sorted(position_columns):
        col = position_columns[n]
        mval = _normalize_residue_series(master.loc[common, col])
        n_letter = int((mval != "").sum())
        rec = dict(column=col, named_position=int(n), n_rows=int(len(common)),
                   n_residue_values=n_letter,
                   n_distinct_values=int(mval[mval != ""].nunique()),
                   value_counts={str(k): int(v) for k, v in
                                 mval[mval != ""].value_counts().head(6).items()},
                   candidates=[], global_top=[], verdict="", named_position_frac=None,
                   best_position=None, decided=False)
        if n_letter < max(1, min_rows):
            rec["verdict"] = (f"판정 불가 — 잔기 값이 {n_letter}개뿐이다(최소 {min_rows}). "
                              f"이 컬럼은 잔기가 아니라 플래그/숫자일 수 있다 "
                              f"(값 분포 {rec['value_counts']}).")
            out["n_undecidable"] += 1
            out["columns"][col] = rec
            out["messages"].append(f"{col}: {rec['verdict']}")
            continue
        mv = mval.to_numpy(dtype="<U1")
        mref = (str(mv[ref_row]) if ref_row >= 0 else None)
        ok = covered_all & (mv != "")[:, None]
        n_cmp = ok.sum(axis=0).astype(np.int64)
        n_eq = (ok & (ours_all == mv[:, None])).sum(axis=0).astype(np.int64)
        frac = np.where(n_cmp > 0, n_eq / np.maximum(n_cmp, 1), np.nan)
        rec["candidates"] = [_cand(p, n, n_cmp[p - 1], n_eq[p - 1], mref)
                             for p in range(max(1, int(n) - window),
                                            min(L, int(n) + window) + 1)]
        rec["reference_row_matches"] = ([int(p) for p in pos_all.tolist()
                                         if mref and mref == str(ref_arr[p - 1])]
                                        if mref else [])
        rec["n_reference_row_matches"] = len(rec["reference_row_matches"])
        rec["named_in_range"] = bool(1 <= int(n) <= L)
        if not rec["named_in_range"]:      # r700 처럼 636 을 넘는 이름. PDB 번호의 흔적일 수 있다
            out["messages"].append(
                f"{col}: 이름의 번호 {n} 이 참조서열 1~{L} 밖이다 — 이 번호 체계는 "
                f"ChCODH2 서열 번호가 아니다. 값으로 어느 위치인지 아래에서 판정한다.")
        usable = (n_cmp >= min_rows)
        pool = usable if int(usable.sum()) else (n_cmp > 0)
        if not int(pool.sum()):
            rec["verdict"] = "판정 불가 — 이 좌표계에서 비교 가능한 행이 없다 (전부 갭/정렬 밖)."
            out["n_undecidable"] += 1
            out["columns"][col] = rec
            out["messages"].append(f"{col}: {rec['verdict']}")
            continue
        order = sorted(pos_all[pool].tolist(),
                       key=lambda p: (-frac[p - 1], abs(p - int(n)), p))
        rec["global_top"] = [_cand(p, n, n_cmp[p - 1], n_eq[p - 1], mref) for p in order[:8]]
        bp = int(order[0])
        rec["best_position"] = bp
        rec["best_frac"] = float(frac[bp - 1])
        rec["named_position_frac"] = (float(frac[int(n) - 1])
                                      if 1 <= int(n) <= L and int(n_cmp[int(n) - 1]) else None)
        # 2등은 1등에서 window 보다 멀리 떨어진 위치만 센다 (바로 옆이 같은 잔기라
        # 덩달아 높게 나오는 것을 2등으로 세면 멀쩡한 판정까지 '불가' 가 된다)
        rival = [p for p in order[1:] if abs(p - bp) > window]
        rec["runner_up_position"] = (int(rival[0]) if rival else None)
        rec["runner_up_frac"] = (float(frac[rival[0] - 1]) if rival else None)
        tag = (f"일치율 최고 {bp}({ref_arr[bp - 1]}) {rec['best_frac']:.3f} "
               f"[{int(n_eq[bp - 1])}/{int(n_cmp[bp - 1])}] — 참조서열 {L} 위치 전부에서 골랐다")
        named_txt = (f"이름 위치 {n} 은 참조서열 1~{L} 밖" if not rec["named_in_range"]
                     else "이름 위치 비교 불가" if rec["named_position_frac"] is None
                     else f"이름 위치 {n} 은 {rec['named_position_frac']:.3f}")
        rival_txt = ("먼 2등 없음" if rec["runner_up_position"] is None else
                     f"먼 2등 {rec['runner_up_position']} 은 {rec['runner_up_frac']:.3f}")
        ref_txt = (f"참조행 잔기 '{mref}' 와 같은 위치 {rec['reference_row_matches'][:8]}"
                   f"{'…' if rec['n_reference_row_matches'] > 8 else ''}"
                   if mref else "참조행 값이 없어 단서로 쓸 수 없다")
        weak = not np.isfinite(rec["best_frac"]) or rec["best_frac"] < float(min_frac)
        close = (rec["runner_up_frac"] is not None
                 and rec["best_frac"] - rec["runner_up_frac"] < float(margin))
        if weak or close:
            why = (f"1등 일치율이 하한 {float(min_frac):.2f} 에 못 미친다" if weak
                   else f"먼 2등과 {float(margin):.2f} 도 벌어지지 않았다")
            rec["verdict"] = (f"판정 불가 — {why}. {tag}, {named_txt}, {rival_txt}. "
                              f"이 컬럼이 어느 위치인지 **값으로 정할 수 없다** — 다른 서열 "
                              f"집합이나 다른 번호 체계에서 온 값이거나 잔기 컬럼이 아닐 수 "
                              f"있다. 04 번은 이 컬럼을 쓰지 말고 이 좌표계의 "
                              f"known_position_equivalents 를 쓸 것. {ref_txt}")
            out["n_undecidable"] += 1
            out["columns"][col] = rec
            out["messages"].append(f"{col}: {rec['verdict']}")
            continue
        rec["decided"] = True
        if bp == int(n):
            rec["verdict"] = (f"이름 그대로 {n} 이다 — {tag}, {rival_txt}. {ref_txt}")
            out["n_decided"] += 1
        else:
            rec["verdict"] = (f"**이름은 r{n} 인데 값은 {bp} 이다** ({bp - int(n):+d} 어긋남). "
                              f"{tag}, {named_txt}, {rival_txt}. {ref_txt}")
            out["n_shifted"] += 1
        out["columns"][col] = rec
        out["messages"].append(f"{col}: {rec['verdict']}")

    shifted = {c: r for c, r in out["columns"].items() if r.get("decided")
               and int(r["best_position"]) != int(r["named_position"])}
    out["shifted_columns"] = {c: int(r["best_position"]) for c, r in shifted.items()}
    out["needs_decision"] = bool(shifted or out["n_undecidable"])
    if shifted:
        out["messages"].append(
            "→ 위 컬럼들은 **이름의 번호가 아닌 다른 위치**의 잔기다. 이 노트북의 좌표는 "
            "참조서열(Met1=1) 기준이므로 영향받지 않지만, master 의 r* 를 그대로 쓰는 "
            "코드는 그 어긋남만큼 틀린 자리를 본다. 04 번에서 master 의 r* 를 쓰지 말고 "
            "이 좌표계의 known_position_equivalents 를 쓸 것.")
    elif out["n_undecidable"]:
        out["messages"].append(
            f"→ 값으로 위치가 **정해진** 컬럼은 이름 그대로였지만, 판정 불가가 "
            f"{out['n_undecidable']}개 남았다. 판정 불가는 '이름이 맞다' 가 아니다 — "
            f"그 컬럼은 어느 자리인지 모르는 것이므로 04 번에서 쓰지 말 것.")
    else:
        out["messages"].append(
            "→ master 의 r* 컬럼이 전부 이름 그대로의 위치였다. 그래도 04 번은 이 좌표계를 "
            "쓴다 (출처가 하나여야 조용히 틀리지 않는다).")
    return out


def master_compare_lines(mcmp: dict, *, focus=("r612",), top: int = 8) -> list[str]:
    """7절이 찍을 줄들 — 이름 기준 대조 → 잔기 값 기준 판정 → 문제 컬럼의 후보 표.

    표를 찍는 컬럼은 ① focus 에 적힌 것 ② 이름과 다른 위치로 판정된 것 ③ 판정 불가인 것.
    '판정 불가' 의 후보 표를 반드시 보여 준다 — 1등이 왜 결정적이지 않은지는 숫자를
    봐야 사람이 판단할 수 있고, 그 판단 없이 r* 를 쓰면 조용히 틀린 자리를 본다.
    """
    cols = ["position", "ref_aa", "delta", "n_compared", "n_equal", "frac",
            "reference_row_equal"]
    want = {str(f).lower() for f in (focus or ())}
    lines = ["", "■ 이름 기준 대조 (02 번과 같은 함수 — 컬럼 이름만 본다)"]
    lines += ["  " + m for m in (mcmp.get("name_based") or {}).get("messages", [])]
    lines += ["", "■ 잔기 값 기준 판정  ★이것이 답이다★"]
    lines += ["  " + m for m in mcmp.get("messages", [])]
    for col, rec in (mcmp.get("columns") or {}).items():
        named, best = rec.get("named_position"), rec.get("best_position")
        show = (str(col).lower() in want or not rec.get("decided")
                or (best is not None and int(best) != int(named)))
        if not show or not rec.get("global_top"):
            continue
        lines += ["", f"■ {col} — 참조서열 전 위치 상위 {min(top, len(rec['global_top']))} "
                      f"(이름 위치 {named}, 값 분포 {rec.get('value_counts')})",
                  pd.DataFrame(rec["global_top"][:top])[cols].round(3).to_string(index=False)]
        if rec.get("candidates"):
            lines += [f"  이름 ±{mcmp.get('window')} 근접 후보",
                      pd.DataFrame(rec["candidates"])[cols].round(3).to_string(index=False)]
        else:
            lines.append(f"  이름 위치 {named} 가 참조서열 1~{mcmp.get('scanned_positions')} "
                         f"밖이라 근접 후보가 없다 — 이름이 다른 번호 체계(PDB 등)다")
    return lines


def tier_coverage(matrix: pd.DataFrame, tiers: dict, ref_seq: str) -> pd.DataFrame:
    """티어마다 위치별 정렬 서열 수. **T1(31개)에서 커버리지가 낮은 위치는 보존도를
    낼 수 없다.** 04 번의 전제가 되는 표다. tiers = {라벨: accession 집합}."""
    frames = []
    for label, members in tiers.items():
        mem = sorted({str(a) for a in members})
        if not mem:
            raise ValueError(f"티어 '{label}' 이 비었다 — 세트를 먼저 고를 것.")
        cov = position_coverage(matrix, ref_seq, members=mem)
        cov.insert(0, "tier", str(label))
        cov.insert(1, "n_members", len(mem))
        frames.append(cov)
    out = pd.concat(frames, ignore_index=True)
    out["residue_frac_tier"] = out["n_residue"] / out["n_members"]
    return out


def tier_coverage_summary(tcov: pd.DataFrame, *, low_frac: float = 0.8,
                          worst: int = 5) -> pd.DataFrame:
    """티어별 한 줄 요약. low_frac 미만인 위치가 '보존도를 믿을 수 없는 자리' 다."""
    rows = []
    for label, g in tcov.groupby("tier", sort=False):
        low = g.loc[g["residue_frac_tier"] < low_frac]
        known = g[g["is_known_position"]]
        kl = known.loc[known["residue_frac_tier"] < low_frac]
        order = known.sort_values(["n_residue", "q_pos"], kind="mergesort").head(worst)
        rows.append(dict(
            tier=str(label), n_members=int(g["n_members"].iloc[0]),
            n_in_alignment=int(g["n_in_alignment"].iloc[0]),
            min_n_residue=int(g["n_residue"].min()),
            median_n_residue=float(g["n_residue"].median()),
            n_pos_zero=int((g["n_residue"] == 0).sum()),
            n_pos_below_frac=int(len(low)),
            known_min_n_residue=(int(known["n_residue"].min()) if len(known) else None),
            n_known_below_frac=int(len(kl)),
            known_below_frac=";".join(f"{int(r.q_pos)}{r.q_aa}" for r in kl.itertuples()),
            known_worst=";".join(f"{int(r.q_pos)}{r.q_aa}:{int(r.n_residue)}"
                                 for r in order.itertuples()),
        ))
    return pd.DataFrame(rows)


def write_residue_matrix(path, matrix: pd.DataFrame) -> dict:
    """행렬을 tsv 로. DataFrame.to_csv 를 쓰지 않고 직접 쓴다 (1,827 x 637 객체열 회피)."""
    path = Path(path)
    vals = matrix.to_numpy(dtype="<U1")
    with open(path, "wt", encoding="utf-8") as fh:
        fh.write("accession\t" + "\t".join(str(c) for c in matrix.columns) + "\n")
        for acc, row in zip(matrix.index.tolist(), vals):
            fh.write(str(acc) + "\t" + "\t".join(row.tolist()) + "\n")
    return dict(path=str(path), n_rows=int(vals.shape[0]),
                n_position_columns=int(vals.shape[1]), bytes=int(path.stat().st_size))


def write_equiv_long(path, equiv: pd.DataFrame, *, chunk: int = 200_000) -> dict:
    """long 테이블(100만행급)을 흘려 쓴다. 04 번이 필요할 때만 켠다 —
    보존도는 행렬만으로 계산된다."""
    path = Path(path)
    n = int(len(equiv))
    with open(path, "wt", encoding="utf-8") as fh:
        fh.write("\t".join(EQUIV_COLUMNS) + "\n")
        for start in range(0, n, int(chunk)):
            equiv.iloc[start:start + int(chunk)].to_csv(
                fh, sep="\t", header=False, index=False, lineterminator="\n")
    return dict(path=str(path), n_rows=n, columns=list(EQUIV_COLUMNS),
                bytes=int(path.stat().st_size))
