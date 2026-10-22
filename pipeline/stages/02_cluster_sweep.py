#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""단계 2 — 클러스터 해상도 스윕. `--min-seq-id` 만 바꿔 분류의 임계값 의존성을 정량화한다.

목표는 하나다. **「결론이 0.70 이라는 임계값 하나에 의존하지 않는다」를 숫자로 보이는 것.**
`-c 0.8 --cov-mode 1 -s 7.5` 는 2026-04 와 동일하게 고정하고 `--min-seq-id` 만 흔든다
(인계문서 1절 — 이 두 옵션을 바꾸면 219 가 재현되지 않으므로 고정이 비교의 전제다).

무엇을 읽는가
  1. IPR010047 FASTA        --fasta    기본 common.IPR_FASTA
       /mnt/af2results/mingyu/chcodh2/protein-matching-IPR010047.fasta  (1,827 서열)
       임계값마다 이 파일 하나로 mmseqs easy-cluster 를 다시 돌린다.
  2. master_1827.tsv (선택) --master   기본 common.MASTER_TSV
       2026-04 의 0.70 분할(219 클러스터)을 **참조 분할**로 쓴다. accession / cluster 두
       컬럼만 읽고 (common.resolve_columns) 각 임계값과 ARI·NMI 를 비교한다. 없으면 건너뛴다.
  3. 마커 목록 (선택)        --markers  기본 없음 (예: pipeline/lib/markers.txt)
       한 줄에 accession 하나. '#' 뒤는 주석, accession 뒤 공백/탭 뒤는 라벨.
       ChCODH1-5 / RrCODH / DvCODH / MtCODH 대조군이 어느 임계값에서 한 묶음이 되는지 본다.
  4. --from-existing DIR 일 때는 FASTA 대신 그 디렉터리(또는 그 아래 sweep/)의
       id<XX>_cluster.tsv 들을 읽는다. mmseqs 를 전혀 부르지 않는다.
  5. pipeline/lib/chcodh2.fasta — 매 실행 common.verify_known_positions() 로 번호 체계 점검.

무엇을 쓰는가  (--outdir, 기본 common.OUT_ROOT/02_cluster_sweep)
  sweep/id<XX>_cluster.tsv      mmseqs easy-cluster 산출물 **원본 그대로**. 가공하지 않는다.
  sweep/id<XX>_rep_seq.fasta    〃 (easy-cluster 가 같이 만든다. 지우지 않는다)
  sweep/id<XX>_all_seqs.fasta   〃
  sweep/params.json             원본 산출물용 재현 정보 (명령·소요시간·tmp 용량·버전)
  logs/id<XX>.log               각 mmseqs 실행의 stdout/stderr 전문
  sweep_metrics.tsv             임계값 1행. 분할 지표 + ChCODH2 추적 + master 대비 ARI
  sweep_adjacent.tsv            인접 임계값 쌍 1행. ARI / NMI / 교차표 / ChCODH2 멤버 변화
  sweep_assignments.tsv         accession 1행 × 임계값 1열 (cluster@XX) 와이드 테이블
  sweep_reference_track.tsv     임계값별 ChCODH2 클러스터의 멤버와 합병 내역
  sweep_markers.tsv             임계값 × 마커. 어느 마커가 같은 클러스터에 모이는지
  sweep_summary.md              사람용 요약 — 안정 구간과 경계의 해석까지
  params.json                   common.write_params (행 수·조인율·결측률·상한·불일치 전부)

무엇을 보장하는가
  1. **서열을 지우지 않는다.** sweep_assignments.tsv 의 행 수는 universe 크기와 항상 같고
     (FASTA 가 있으면 FASTA 레코드 수, --from-existing 이면 전 임계값 accession 합집합),
     어느 임계값에서 배정이 없는 칸은 빈 문자열로 남긴다. 삭제가 아니라 공백이다.
  2. **구조적 가정이 깨지면 예외**다 — 행 수 변화, cluster.tsv 컬럼 수 != 2, 임계값 라벨 충돌,
     한 accession 이 한 임계값에서 서로 다른 클러스터 둘에 속함, easy-cluster 가 서열을
     빠뜨림(--allow-incomplete 로만 완화), 유효하지 않은 임계값.
     **데이터 기대값**(예: master 가 219 클러스터, universe 가 1,827행)은 예외가 아니라
     `discrepancies` 로 모아 stdout·요약 md·params.json 세 곳에 적고, --strict 면 산출물을
     모두 쓴 뒤 종료코드 2 로 끝낸다. (단계 0과 같은 규약. 산출물을 잃지 않기 위해 마지막에 끝낸다.)
  3. 행 수·조인율·결측률을 항상 세어 stdout 과 params.json 양쪽에 적는다.
  4. **상한에 걸려 잘린 것은 개수를 명시한다** — 멤버 목록(--max-members), 합병 출처 클러스터
     목록(--max-list), md 표 행 수(--max-table-rows). 잘리면 'truncated' 플래그와 전체 개수를 함께 쓴다.
  5. 모든 잔기/서열 번호는 ChCODH2 번호(Met1 = 1, 636 aa). 이 단계는 위치를 다루지 않지만
     번호 체계는 매 실행 점검한다.
  6. tmp 를 임계값마다 지운다(--keep-tmp 로 보존). 지우기 전 용량을 재서 기록한다.
     실행 전 tmp 가 비어 있지 않으면 예외 — 인계문서 4절의 「남아 있으면 다음 실행이 꼬인다」.

지표의 정의 (오독하기 쉬운 것만)
  singleton_frac       common.partition_summary 의 정의 그대로 = 싱글턴 **클러스터 수 / 클러스터 수**.
  singleton_seq_frac   이 스크립트가 따로 계산 = 싱글턴 클러스터에 속한 **서열 수 / 서열 수**.
                       두 값은 다르다. 「싱글턴 비율」이라고만 쓰면 반드시 오독된다.
  finer / coarser      인접 쌍은 항상 (낮은 임계값 = coarser, 높은 임계값 = finer) 로 적는다.
                       min-seq-id 가 내려가면 클러스터는 합쳐지므로 coarser 가 더 큰 묶음이다.
  n_coarser_merged     낮은 임계값의 클러스터 중 높은 임계값 클러스터 2개 이상을 흡수한 것의 수.
  n_finer_split        높은 임계값 클러스터 하나가 낮은 임계값 클러스터 여러 개로 흩어진 수.
                       **계층이 정상이면 0 이다.** 0 이 아니면 포함관계가 깨졌다는 신호.
  is_mm_rep@XX         이 스크립트는 따로 쓰지 않는다. (cluster@XX == accession) 으로 유도하라
                       (인계문서 1절 규칙).

실행 환경
  대상은 epel-af2 워크스테이션(conda base, 사용자 unsit). 의존성은 pandas + numpy + 표준
  라이브러리뿐. mmseqs 는 common.find_mmseqs() 로 **절대경로를 탐색 후 검증**한다
  (인계문서 6절 — Jupyter 의 PATH 를 믿지 않는다).
  mmseqs 가 없는 환경에서 전 경로를 확인하는 길이 둘 있다.
    --dry-run         돌릴 명령만 출력하고 아무것도 쓰지 않는다.
    --from-existing   이미 있는 id<XX>_cluster.tsv 로 지표 계산만 다시 한다.
    --self-test DIR   DIR 에 합성 cluster.tsv 를 만들어 --from-existing 전 경로를 실제로 돌리고
                      계획된 분할 구조(안정 구간/경계/ChCODH2 멤버)와 대조한다.

사용 예
  source pipeline/config.sh
  python3 pipeline/stages/02_cluster_sweep.py --markers pipeline/lib/markers.txt
  python3 pipeline/stages/02_cluster_sweep.py --dry-run
  python3 pipeline/stages/02_cluster_sweep.py --thresholds 0.5 0.6 0.7 0.8 --threads 8
  python3 pipeline/stages/02_cluster_sweep.py --from-existing <이전 outdir> --outdir /tmp/re
  python3 pipeline/stages/02_cluster_sweep.py --self-test /tmp/sweep_selftest
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve()
_LIB = _HERE.parents[1] / "lib"
if str(_LIB) not in sys.path:
    sys.path.insert(0, str(_LIB))
import common  # noqa: E402

STAGE = "02_cluster_sweep"

# --------------------------------------------------------------------------
# 고정 상수 — 인계문서 1절. 이 세 값은 2026-04 의 219 와 비교 가능성의 전제다.
# --------------------------------------------------------------------------
FIXED_COV = "0.8"
FIXED_COV_MODE = "1"
FIXED_SENSITIVITY = "7.5"

# --min-seq-id 기본 스윕 — 0.10~0.90, 0.05 간격. 0.70 은 운용 임계값이라 반드시 포함된다.
DEFAULT_THRESHOLDS = [round(p / 100.0, 2) for p in range(10, 91, 5)]
REFERENCE_THRESHOLD = 0.70            # master_1827.tsv 의 219 분할이 만들어진 값
EXPECTED_MASTER_CLUSTERS = 219        # 인계문서 1절
EXPECTED_N_SEQUENCES = 1827           # 인계문서 1절

# 안정 구간 판정 기본값 — 요약의 해석에 쓰이고 params.json 에 그대로 기록된다.
DEFAULT_STABLE_ARI = 0.90

# 상한 (걸리면 전부 개수를 명시해 보고한다)
DEFAULT_MAX_MEMBERS = 50              # ChCODH2 클러스터 멤버 목록
DEFAULT_MAX_LIST = 20                 # 합병 출처 클러스터 목록
DEFAULT_MAX_TABLE_ROWS = 60           # 요약 md 표 행 수

CLUSTER_TSV_GLOB = "id*_cluster.tsv"

# sweep_adjacent.tsv 의 컬럼. 임계값이 하나뿐이면 행이 0개가 되는데, 그때도 헤더는 써야
# 뒤 단계의 pd.read_csv 가 EmptyDataError 로 죽지 않는다.
ADJACENT_COLUMNS = [
    "label_low", "label_high", "min_seq_id_low", "min_seq_id_high", "n_common",
    "ari", "nmi", "n_clusters_coarser", "n_clusters_finer", "delta_n_clusters",
    "finer_maps_to_one", "coarser_maps_to_one", "one_to_one", "n_finer_split",
    "n_coarser_merged", "is_stable", "ref_jaccard", "ref_n_added", "ref_n_removed",
    "ref_n_source_clusters", "ref_merged_from", "ref_merged_truncated",
]
# cluster.tsv 는 헤더가 없는 2열 파일이다. 사람이 만든 파일에 헤더가 붙은 경우만 걸러낸다.
HEADER_TOKENS = {"cluster", "cluster_rep", "rep", "representative", "accession",
                 "member", "target", "query", "acc"}


# ==========================================================================
# accession / 임계값 라벨
# ==========================================================================
def safe_parse_accession(token) -> str:
    """common.parse_accession 과 같은 규칙. 빈 셀에서 IndexError 를 내지 않는다.

    (common.parse_accession 은 빈 문자열에서 터진다. 단계 1도 같은 래퍼를 들고 있으므로
     common.py 에 올릴 후보다.)
    """
    token = "" if token is None else str(token).strip()
    if not token:
        return ""
    return common.parse_accession(token)


def norm_acc_series(series: pd.Series) -> pd.Series:
    s = series.astype("object").fillna("").astype(str).str.strip()
    return s.map(safe_parse_accession)


def threshold_label(t: float) -> str:
    """0.70 -> 'id70', 0.125 -> 'id12p5'. 파일 이름과 열 이름의 단일 출처."""
    pct = t * 100.0
    if abs(pct - round(pct)) < 1e-9:
        return f"id{int(round(pct)):02d}"
    s = f"{pct:.4f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"id{s}"


def label_to_threshold(label: str) -> float:
    """'id70' -> 0.70, 'id12p5' -> 0.125. threshold_label 의 역함수."""
    m = re.fullmatch(r"id(\d+(?:p\d+)?)", label)
    if not m:
        raise ValueError(f"임계값 라벨을 해석할 수 없다: {label!r} (기대 형식 'id70' / 'id12p5')")
    return float(m.group(1).replace("p", ".")) / 100.0


def column_name(label: str) -> str:
    """'id70' -> 'cluster@70'. 뒤 단계에서 master 테이블에 붙는 열 이름."""
    return "cluster@" + label[2:]


def fmt_min_seq_id(t: float) -> str:
    """mmseqs 에 넘기는 문자열. 0.70 은 인계문서대로 '0.70' 으로 적는다."""
    s = f"{t:.2f}"
    if abs(float(s) - t) < 1e-12:
        return s
    return f"{t:.6f}".rstrip("0").rstrip(".")


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


def parse_thresholds(values, *, force_reference: bool, reference: float) -> tuple[list[float], dict]:
    """임계값 목록을 정규화한다. 중복 제거·정렬·라벨 충돌 검사까지 여기서 끝낸다."""
    raw = list(values)
    bad = [v for v in raw if not (0.0 < v <= 1.0)]
    if bad:
        raise ValueError(
            f"--thresholds 값이 (0, 1] 범위를 벗어났다: {bad}. "
            f"mmseqs --min-seq-id 는 분율이다 (70% 는 0.70, 70 이 아니다).")

    seen, uniq = set(), []
    for v in raw:
        key = round(v, 6)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(key)

    added_reference = False
    if force_reference and round(reference, 6) not in seen:
        uniq.append(round(reference, 6))
        added_reference = True

    uniq.sort()
    labels = make_labels(uniq)

    meta = dict(
        thresholds=uniq,
        labels=labels,
        n_requested=len(raw),
        n_duplicates_dropped=len(raw) - len(seen),
        reference_threshold=reference,
        reference_threshold_added=added_reference,
    )
    return uniq, meta


# ==========================================================================
# 입력 읽기
# ==========================================================================
def read_cluster_tsv(path: Path) -> tuple[pd.DataFrame, dict]:
    """easy-cluster 의 _cluster.tsv (헤더 없는 2열: 대표 \t 멤버) 를 읽는다.

    common.read_tsv 를 쓰지 않는 이유 — 그 함수는 '컬럼명에 | 또는 숫자' 로 헤더 유무를
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


def read_master_partition(path: Path, *, allow_duplicate_accession: bool) -> tuple[pd.Series | None, dict]:
    """master_1827.tsv 의 cluster 컬럼을 참조 분할로 읽는다. 읽기 전용."""
    info = dict(path=str(path), exists=bool(Path(path).exists()), skipped=False)
    if not Path(path).exists():
        info.update(skipped=True, reason="파일 없음")
        return None, info

    m = common.read_tsv(path)
    info["columns"] = list(m.columns)
    info["n_rows"] = int(len(m))
    cols = common.resolve_columns(m, ["accession", "cluster"], required=("accession",))
    acc_col, cl_col = cols["accession"], cols["cluster"]
    info.update(accession_column=acc_col, cluster_column=cl_col)
    if cl_col is None:
        raise KeyError(
            f"master 에서 cluster 컬럼을 못 찾음. 후보={common.COLSPEC['cluster']} 실제={list(m.columns)}")

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


# ==========================================================================
# mmseqs 실행
# ==========================================================================
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


def build_cluster_cmd(mmseqs: str, fasta: Path, prefix: Path, tmpdir: Path, t: float,
                      *, cov: str, cov_mode: str, sensitivity: str, threads: int) -> list[str]:
    """인계문서 1절의 명령 그대로. --min-seq-id 만 스윕 값으로 바뀐다."""
    return [str(mmseqs), "easy-cluster", str(fasta), str(prefix), str(tmpdir),
            "--min-seq-id", fmt_min_seq_id(t),
            "-c", str(cov), "--cov-mode", str(cov_mode), "-s", str(sensitivity),
            "--threads", str(threads)]


def run_sweep(args, mmseqs: str, thresholds: list[float], sweep_dir: Path,
              logs_dir: Path, tmp_root: Path) -> list[dict]:
    """임계값마다 easy-cluster 를 돌린다. 명령·소요시간·tmp 용량을 전부 기록한다."""
    fasta = Path(args.fasta)
    sweep_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    if tmp_root.exists() and any(tmp_root.iterdir()) and not args.force:
        raise FileExistsError(
            f"tmp 루트에 이전 실행의 잔여물이 있다: {tmp_root} "
            f"(인계문서 4절 — 남아 있으면 다음 실행이 꼬인다). "
            f"--force 로 지우고 진행하거나 --tmp-root 로 다른 경로를 줄 것.")
    if tmp_root.exists() and args.force:
        shutil.rmtree(tmp_root, ignore_errors=True)
    tmp_root.mkdir(parents=True, exist_ok=True)

    records = []
    for i, t in enumerate(thresholds, 1):
        label = threshold_label(t)
        prefix = sweep_dir / label
        out_tsv = prefix.with_name(f"{label}_cluster.tsv")
        tmpdir = tmp_root / label

        if args.resume and out_tsv.exists():
            print(f"  [{i}/{len(thresholds)}] {label} (min-seq-id {fmt_min_seq_id(t)}) "
                  f"— 이미 있어 건너뜀 (--resume): {out_tsv}")
            records.append(dict(label=label, min_seq_id=t, source="existing(resume)",
                                cluster_tsv=str(out_tsv), seconds=None, cmd_str=None,
                                skipped_resume=True))
            continue

        cmd = build_cluster_cmd(mmseqs, fasta, prefix, tmpdir, t,
                                cov=args.cov, cov_mode=args.cov_mode,
                                sensitivity=args.sensitivity, threads=args.threads)
        print(f"  [{i}/{len(thresholds)}] {label} (min-seq-id {fmt_min_seq_id(t)})")
        print(f"      $ {' '.join(cmd)}")
        tmpdir.mkdir(parents=True, exist_ok=True)
        rec = common.run_cmd(cmd, log_path=logs_dir / f"{label}.log", check=True)

        tmp_size = dir_bytes(tmpdir)
        if args.keep_tmp:
            tmp_action = "보존(--keep-tmp)"
        else:
            shutil.rmtree(tmpdir, ignore_errors=True)
            tmp_action = "삭제"

        if not out_tsv.exists():
            raise FileNotFoundError(
                f"{label}: easy-cluster 가 끝났는데 {out_tsv} 가 없다. "
                f"로그를 볼 것: {logs_dir / (label + '.log')}\n{rec['stderr_tail'][-1000:]}")

        print(f"      완료 {rec['seconds']}초 · tmp {tmp_size / (1 << 20):.1f} MiB {tmp_action} · "
              f"{out_tsv.name}")
        records.append(dict(label=label, min_seq_id=t, source="mmseqs",
                            cluster_tsv=str(out_tsv), seconds=rec["seconds"],
                            cmd=rec["cmd"], cmd_str=rec["cmd_str"],
                            returncode=rec["returncode"], started_utc=rec["started_utc"],
                            tmp_bytes=tmp_size, tmp_action=tmp_action,
                            stderr_tail=rec["stderr_tail"][-1500:], skipped_resume=False))

    if not args.keep_tmp:
        shutil.rmtree(tmp_root, ignore_errors=True)
    return records


def discover_existing(from_dir: Path, thresholds_requested: list[float] | None) -> tuple[list[dict], dict]:
    """--from-existing. id<XX>_cluster.tsv 들을 찾아 임계값을 역산한다."""
    cand_dirs = [from_dir / "sweep", from_dir]
    search_dir, found = None, []
    for d in cand_dirs:
        if d.is_dir():
            hits = sorted(d.glob(CLUSTER_TSV_GLOB))
            if hits:
                search_dir, found = d, hits
                break
    if not found:
        raise FileNotFoundError(
            f"{from_dir} (또는 그 아래 sweep/) 에서 {CLUSTER_TSV_GLOB} 를 찾지 못했다.")

    records, unparsable = [], []
    for p in found:
        label = p.name[: -len("_cluster.tsv")]
        try:
            t = label_to_threshold(label)
        except ValueError:
            unparsable.append(p.name)
            continue
        records.append(dict(label=label, min_seq_id=t, source="existing",
                            cluster_tsv=str(p), seconds=None, cmd_str=None,
                            skipped_resume=False))
    records.sort(key=lambda r: r["min_seq_id"])

    meta = dict(from_dir=str(from_dir), search_dir=str(search_dir),
                n_files_found=len(found), n_files_unparsable=len(unparsable),
                unparsable_names=unparsable[:10])

    if thresholds_requested:
        want = {round(t, 6) for t in thresholds_requested}
        have = {round(r["min_seq_id"], 6) for r in records}
        kept = [r for r in records if round(r["min_seq_id"], 6) in want]
        meta.update(
            filtered_by_thresholds=True,
            n_requested=len(want),
            n_missing=len(want - have),
            missing_thresholds=sorted(want - have),
            n_ignored_extra=len(records) - len(kept),
        )
        if not kept:
            raise FileNotFoundError(
                f"요청한 임계값 {sorted(want)} 에 해당하는 파일이 {search_dir} 에 없다. "
                f"있는 것: {sorted(have)}")
        records = kept
    else:
        meta.update(filtered_by_thresholds=False, n_missing=0, missing_thresholds=[],
                    n_ignored_extra=0)
    if not records:
        raise FileNotFoundError(f"{search_dir} 에서 쓸 수 있는 cluster.tsv 가 없다.")
    return records, meta


# ==========================================================================
# universe 와 와이드 테이블
# ==========================================================================
def fasta_universe(fasta: Path) -> tuple[list[str], dict]:
    accs, n_records = [], 0
    for header, _seq in common.read_fasta(fasta):
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


def build_assignments(partitions: dict[str, pd.Series], universe: list[str],
                      order: list[str]) -> tuple[pd.DataFrame, dict]:
    """accession × 임계값 와이드 테이블. 배정이 없는 칸은 빈 문자열 (삭제가 아니다)."""
    idx = pd.Index(universe, name="accession")
    wide = pd.DataFrame(index=idx)
    per_label = {}
    for label in order:
        ser = partitions[label]
        col = ser.reindex(idx)
        n_missing = int(col.isna().sum())
        extra = int(len(ser.index.difference(idx)))
        per_label[label] = dict(n_assigned=int(len(idx) - n_missing),
                                n_missing_vs_universe=n_missing,
                                n_outside_universe=extra)
        wide[column_name(label)] = col.fillna("").astype(str)

    cols = [column_name(lb) for lb in order]
    assigned = (wide[cols] != "")
    wide.insert(0, "n_thresholds_assigned", assigned.sum(axis=1).astype(int))
    meta = dict(n_rows=int(len(wide)), n_threshold_columns=len(cols),
                per_label=per_label,
                n_rows_fully_assigned=int((wide["n_thresholds_assigned"] == len(cols)).sum()),
                n_rows_partially_assigned=int(((wide["n_thresholds_assigned"] > 0) &
                                               (wide["n_thresholds_assigned"] < len(cols))).sum()),
                n_rows_never_assigned=int((wide["n_thresholds_assigned"] == 0).sum()))
    return wide, meta


# ==========================================================================
# 지표
# ==========================================================================
def singleton_seq_frac(ser: pd.Series) -> tuple[int, float]:
    sizes = ser.value_counts()
    n_seq_in_singletons = int((sizes == 1).sum())     # 싱글턴 클러스터 1개 = 서열 1개
    frac = float(n_seq_in_singletons / len(ser)) if len(ser) else float("nan")
    return n_seq_in_singletons, frac


def compare_partitions(a: pd.Series, b: pd.Series) -> dict:
    """공통 accession 위에서만 비교한다. 공통 크기와 피복률을 항상 같이 돌려준다."""
    common_idx = a.index.intersection(b.index)
    n_common = int(len(common_idx))
    out = dict(n_common=n_common,
               common_frac_a=float(n_common / len(a)) if len(a) else float("nan"),
               common_frac_b=float(n_common / len(b)) if len(b) else float("nan"))
    if n_common < 2:
        out.update(ari=float("nan"), nmi=float("nan"))
        return out
    la = a.reindex(common_idx).astype(str).tolist()
    lb = b.reindex(common_idx).astype(str).tolist()
    out.update(ari=common.adjusted_rand_index(la, lb),
               nmi=common.normalized_mutual_info(la, lb),
               **common.crosstab_mapping(la, lb))
    return out


def metrics_table(partitions: dict[str, pd.Series], order: list[str], run_records: list[dict],
                  file_meta: dict, part_meta: dict, assign_meta: dict,
                  master: pd.Series | None, ref_rows: dict, marker_rows: dict) -> pd.DataFrame:
    rows = []
    rec_by_label = {r["label"]: r for r in run_records}
    for label in order:
        ser = partitions[label]
        t = label_to_threshold(label)
        ps = common.partition_summary(ser.tolist())
        n_sing_seq, sing_seq_frac = singleton_seq_frac(ser)
        row = dict(
            label=label,
            min_seq_id=t,
            n_sequences=ps["n_sequences"],
            n_clusters=ps["n_clusters"],
            n_singletons=ps["n_singletons"],
            singleton_frac=ps["singleton_frac"],
            singleton_seq_frac=sing_seq_frac,
            largest_cluster=ps["largest_cluster"],
            largest_frac=ps["largest_frac"],
            median_size=ps["median_size"],
            mean_size=ps["mean_size"],
            n_rows_in_file=file_meta[label]["n_rows"],
            n_duplicate_accession_rows=part_meta[label]["n_duplicate_accession_rows"],
            n_missing_vs_universe=assign_meta["per_label"][label]["n_missing_vs_universe"],
            n_outside_universe=assign_meta["per_label"][label]["n_outside_universe"],
        )
        r = ref_rows.get(label, {})
        row.update(ref_found=r.get("ref_found", False), ref_cluster=r.get("ref_cluster", ""),
                   ref_cluster_size=r.get("ref_cluster_size", 0),
                   ref_is_mm_rep=r.get("ref_is_mm_rep", False))
        mk = marker_rows.get(label, {})
        row.update(n_markers_found=mk.get("n_markers_found", 0),
                   n_marker_clusters=mk.get("n_marker_clusters", 0),
                   markers_all_together=mk.get("markers_all_together", ""))
        if master is not None:
            c = compare_partitions(master, ser)
            row.update(ari_vs_master=c["ari"], nmi_vs_master=c["nmi"],
                       n_common_master=c["n_common"],
                       master_one_to_one=c.get("one_to_one", 0),
                       master_n_split=c.get("n_split", 0),
                       master_n_merged=c.get("n_merged", 0))
        else:
            row.update(ari_vs_master=float("nan"), nmi_vs_master=float("nan"),
                       n_common_master=0, master_one_to_one=0,
                       master_n_split=0, master_n_merged=0)
        rec = rec_by_label.get(label, {})
        row.update(seconds=rec.get("seconds"), source=rec.get("source", "existing"),
                   cluster_tsv=rec.get("cluster_tsv", ""))
        rows.append(row)
    return pd.DataFrame(rows)


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


def cluster_members(ser: pd.Series, acc: str) -> list[str] | None:
    """acc 가 속한 클러스터의 멤버 목록. 없으면 None."""
    if acc not in ser.index:
        return None
    rep = ser.loc[acc]
    if isinstance(rep, pd.Series):          # 중복 index 방어 (정상 경로에서는 발생하지 않는다)
        rep = rep.iloc[0]
    return sorted(ser.index[ser.to_numpy() == rep].tolist())


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


def interpret_adjacent(df: pd.DataFrame, order: list[str], *, stable_ari: float) -> dict:
    """ARI 가 유지되는 '안정 구간' 과 급락하는 '경계' 를 뽑는다. 요약 md 의 해석 재료."""
    out = dict(stable_ari=stable_ari, n_pairs=int(len(df)))
    if df.empty:
        out.update(stable_intervals=[], boundaries=[], ari_median=None, ari_min=None,
                   reference_in_stable_interval=None, ref_membership_intervals=[])
        return out

    thr = [label_to_threshold(lb) for lb in order]
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


# ==========================================================================
# (d) ChCODH2 추적
# ==========================================================================
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


# ==========================================================================
# (e) 마커 동거
# ==========================================================================
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


# ==========================================================================
# 요약 md
# ==========================================================================
def md_table(headers, rows) -> str:
    out = ["| " + " | ".join(str(h) for h in headers) + " |",
           "|" + "|".join(["---"] * len(headers)) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
    return "\n".join(out)


def fnum(x, nd=3) -> str:
    if x is None:
        return "—"
    try:
        if isinstance(x, float) and np.isnan(x):
            return "—"
    except TypeError:
        return str(x)
    if isinstance(x, float):
        return f"{x:.{nd}f}"
    return str(x)


def write_summary(path: Path, metrics: pd.DataFrame, adjacent: pd.DataFrame,
                  ref_track: pd.DataFrame, markers_df: pd.DataFrame, stats: dict) -> None:
    S = stats
    interp = S["adjacent_interpretation"]
    L = []
    L.append("# 단계 2 — 클러스터 해상도 스윕 요약\n")
    L.append(f"작성 {datetime.now(timezone.utc).isoformat()} · host `{os.uname().nodename}`\n")
    L.append("`--min-seq-id` 만 흔들고 `-c {} --cov-mode {} -s {}` 는 고정했다. "
             "목적은 **결론이 하나의 임계값에 의존하지 않음을 보이는 것**이다.\n"
             .format(S["fixed_options"]["cov"], S["fixed_options"]["cov_mode"],
                     S["fixed_options"]["sensitivity"]))

    L.append("## 1. 실행 개요\n")
    L.append(md_table(
        ["항목", "값"],
        [["모드", S["mode"]],
         ["입력 FASTA", S["input_paths"]["fasta"]],
         ["master (참조 분할)", S["master"].get("path", "—")],
         ["mmseqs", S["mmseqs"]["binary"]],
         ["mmseqs 버전", S["mmseqs"]["version"]],
         ["threads", S["mmseqs"]["threads"]],
         ["임계값 수", f"{len(S['thresholds']['thresholds'])}개 "
                       f"({min(S['thresholds']['thresholds'])}~{max(S['thresholds']['thresholds'])})"],
         ["0.70 포함", "예" if interp.get("reference_threshold_present") else "아니오"],
         ["universe", f"{S['universe']['n_accessions']}개 accession "
                      f"(출처 {S['universe']['source']})"],
         ["총 mmseqs 소요", f"{S['timing']['total_mmseqs_seconds']}초"
                            if S["timing"]["total_mmseqs_seconds"] is not None else "—"],
         ["안정 판정 기준", f"인접 ARI ≥ {interp['stable_ari']}"]]) + "\n")

    L.append("## 2. 임계값별 분할 지표\n")
    L.append("`singleton_frac` 은 **싱글턴 클러스터 수 / 클러스터 수**이고 "
             "`singleton_seq_frac` 은 **싱글턴 서열 수 / 전체 서열 수**다. 둘은 다르다.\n")
    head = ["min_seq_id", "클러스터", "싱글턴", "sing/clu", "sing/seq", "최대클러스터",
            "최대점유", "중앙크기", "ChCODH2 클러스터 크기", "ARI vs master", "초"]
    rows = []
    for r in metrics.itertuples():
        rows.append([fnum(r.min_seq_id, 2), r.n_clusters, r.n_singletons,
                     fnum(r.singleton_frac), fnum(r.singleton_seq_frac), r.largest_cluster,
                     fnum(r.largest_frac), fnum(r.median_size, 1),
                     (r.ref_cluster_size if r.ref_found else "없음"),
                     fnum(r.ari_vs_master), fnum(r.seconds, 1) if r.seconds else "—"])
    L.append(md_table(head, rows[:S["limits"]["max_table_rows"]]) + "\n")
    if len(rows) > S["limits"]["max_table_rows"]:
        L.append(f"*표를 {S['limits']['max_table_rows']}행으로 잘랐다 "
                 f"(전체 {len(rows)}행). 전부는 `sweep_metrics.tsv` 참조.*\n")

    L.append("## 3. 인접 임계값 안정성 (ARI / NMI)\n")
    L.append("낮은 임계값이 coarser(큰 묶음), 높은 임계값이 finer 다. "
             "`n_coarser_merged` 는 낮은 쪽 클러스터가 높은 쪽 클러스터 2개 이상을 흡수한 수, "
             "`n_finer_split` 은 계층이 깨진 수(정상이면 0)다.\n")
    head = ["구간", "ARI", "NMI", "클러스터 coarse→fine", "1:1", "합쳐짐", "쪼개짐",
            "ChCODH2 jaccard", "판정"]
    rows = []
    for r in adjacent.itertuples():
        rows.append([f"{fnum(r.min_seq_id_low, 2)} ↔ {fnum(r.min_seq_id_high, 2)}",
                     fnum(r.ari), fnum(r.nmi),
                     f"{r.n_clusters_coarser} → {r.n_clusters_finer}",
                     r.one_to_one, r.n_coarser_merged, r.n_finer_split,
                     fnum(r.ref_jaccard), "안정" if r.is_stable else "**경계**"])
    L.append(md_table(head, rows[:S["limits"]["max_table_rows"]]) + "\n")
    if len(rows) > S["limits"]["max_table_rows"]:
        L.append(f"*표를 {S['limits']['max_table_rows']}행으로 잘랐다 (전체 {len(rows)}행).*\n")

    L.append("### 해석\n")
    if interp["n_pairs"] == 0:
        L.append("- 임계값이 하나뿐이라 인접 비교가 없다.\n")
    else:
        L.append(f"- 인접 ARI 중앙값 **{fnum(interp['ari_median'])}** · "
                 f"최솟값 {fnum(interp['ari_min'])} · 최댓값 {fnum(interp['ari_max'])}.")
        if interp["stable_intervals"]:
            for iv in interp["stable_intervals"]:
                L.append(f"- 안정 구간 **{fnum(iv['low'], 2)}–{fnum(iv['high'], 2)}** "
                         f"({iv['n_steps']} 스텝, 최소 ARI {fnum(iv['min_ari'])}) — "
                         f"이 범위 안에서는 임계값을 흔들어도 분할이 거의 같다.")
        else:
            L.append(f"- 인접 ARI ≥ {interp['stable_ari']} 를 만족하는 구간이 **없다**. "
                     f"분류가 임계값에 민감하다는 뜻이므로 결론을 단일 임계값으로 말할 수 없다.")
        for b in interp["boundaries"][:10]:
            L.append(f"- 경계 **{fnum(b['low'], 2)} ↔ {fnum(b['high'], 2)}** "
                     f"(ARI {fnum(b['ari'])}, 합쳐진 클러스터 {b['n_coarser_merged']}개) — "
                     f"여기서 분할이 급격히 재편된다.")
        if len(interp["boundaries"]) > 10:
            L.append(f"- *경계 {len(interp['boundaries'])}개 중 10개만 적었다 "
                     f"(ARI 낮은 순). 전부는 `sweep_adjacent.tsv`.*")
        hit = interp["reference_in_stable_interval"]
        if not interp.get("reference_threshold_present"):
            L.append(f"- ⚠ 운용 임계값 {REFERENCE_THRESHOLD:.2f} 가 이번 스윕에 없다. "
                     f"219 분할과의 직접 비교가 불가능하다.")
        elif hit:
            L.append(f"- **운용 임계값 {REFERENCE_THRESHOLD:.2f} 는 안정 구간 "
                     f"{fnum(hit['low'], 2)}–{fnum(hit['high'], 2)} 안에 있다.** "
                     f"0.70 에서 얻은 분할은 그 구간 어디를 골라도 거의 같으므로, "
                     f"결론이 0.70 이라는 한 점에 의존하지 않는다.")
        else:
            L.append(f"- ⚠ **운용 임계값 {REFERENCE_THRESHOLD:.2f} 가 안정 구간 밖(경계 근처)이다.** "
                     f"0.70 기반 결론은 임계값 선택에 민감할 수 있으므로 "
                     f"근연 클러스터 정의를 ±0.05 에서 재확인해야 한다.")
        L.append("")

    L.append("## 4. ChCODH2 자신의 클러스터 추적\n")
    ri = S["reference"]
    if not ri["found_in_any_threshold"]:
        L.append(f"- 참조 `{ri['reference_accession']}` 가 어느 임계값의 cluster.tsv 에도 없다. "
                 f"**이 분석은 건너뛰었다.** `--reference-acc` 로 실제 쓰이는 accession 을 주거나, "
                 f"참조가 IPR010047 세트에 포함되지 않았는지 확인할 것.\n")
    else:
        L.append(f"참조 `{ri['reference_accession']}` · 발견된 임계값 "
                 f"{ri['n_thresholds_found']}/{len(metrics)} · 클러스터 크기 "
                 f"{ri['size_min']}~{ri['size_max']}\n")
        head = ["min_seq_id", "대표", "크기", "자기대표", "멤버 일부"]
        rows = []
        for r in ref_track.itertuples():
            if not r.ref_found:
                rows.append([fnum(r.min_seq_id, 2), "없음", 0, "", ""])
                continue
            mem = r.members.split(";")[:6]
            more = f" … (+{r.ref_cluster_size - len(mem)}개)" if r.ref_cluster_size > len(mem) else ""
            rows.append([fnum(r.min_seq_id, 2), r.ref_cluster, r.ref_cluster_size,
                         "예" if r.ref_is_mm_rep else "아니오", ", ".join(mem) + more])
        L.append(md_table(head, rows[:S["limits"]["max_table_rows"]]) + "\n")
        L.append(f"*멤버 목록은 임계값마다 최대 {S['limits']['max_members']}개까지만 "
                 f"`sweep_reference_track.tsv` 에 적는다 (`members_truncated` 로 표시).*\n")
        if interp.get("ref_membership_intervals"):
            for iv in interp["ref_membership_intervals"]:
                L.append(f"- ChCODH2 클러스터 멤버가 **전혀 바뀌지 않는 구간**: "
                         f"{fnum(iv['low'], 2)}–{fnum(iv['high'], 2)} ({iv['n_steps']} 스텝)")
        if interp.get("ref_merge_steps"):
            for m in interp["ref_merge_steps"]:
                L.append(f"- **{fnum(m['high'], 2)} → {fnum(m['low'], 2)} 에서 합병**: "
                         f"상위 임계값 클러스터 {m['n_source_clusters']}개가 하나로 묶이며 "
                         f"서열 {m['n_added']}개가 들어왔다. 출처 대표: {m['merged_from']}")
        else:
            L.append("- 인접 스텝에서 ChCODH2 클러스터가 다른 클러스터를 흡수한 사건이 없다.")
        L.append("")

    L.append("## 5. 대조군 마커 동거\n")
    mi = S["markers"]
    if mi.get("skipped"):
        L.append(f"- `--markers` 가 없어 **건너뛰었다** ({mi.get('reason')}). "
                 f"ChCODH1-5 / RrCODH / DvCODH / MtCODH 를 넣으려면 "
                 f"`--markers pipeline/lib/markers.txt`.\n")
    else:
        L.append(f"마커 {mi['n_markers']}개 중 universe 에 있는 것 {mi['n_markers_in_universe']}개."
                 + (f" 없는 것: {', '.join(mi['missing_markers'])}" if mi["missing_markers"] else "")
                 + "\n")
        head = ["min_seq_id", "발견", "마커가 흩어진 클러스터 수", "전부 한 묶음"]
        rows = [[fnum(r.min_seq_id, 2), r.n_markers_found, r.n_marker_clusters,
                 r.markers_all_together] for r in metrics.itertuples()]
        L.append(md_table(head, rows[:S["limits"]["max_table_rows"]]) + "\n")
        pairs = [p for p in mi["pairs"] if p["n_thresholds_together"]]
        if mi["n_pairs"] == 0:
            L.append("- 마커가 1개뿐이라 쌍 비교가 없다. 동거 추적을 의미 있게 쓰려면 "
                     "`markers.txt` 의 주석 처리된 accession 들을 채워야 한다.\n")
        elif not pairs:
            L.append("- 어떤 임계값에서도 같은 클러스터에 들어간 마커 쌍이 없다.\n")
        else:
            head = ["마커 A", "마커 B", "같은 묶음이 되는 가장 높은 임계값", "함께인 임계값 수"]
            rows = [[p["marker_a"], p["marker_b"], fnum(p["highest_together"], 2),
                     p["n_thresholds_together"]]
                    for p in sorted(pairs, key=lambda p: -(p["highest_together"] or 0))]
            L.append(md_table(head, rows[:S["limits"]["max_table_rows"]]) + "\n")

    L.append("## 6. 기존 219 분할과의 대조\n")
    mas = S["master"]
    if mas.get("skipped"):
        L.append(f"- master 를 읽지 않았다 ({mas.get('reason')}). 참조 분할 비교는 건너뛰었다.\n")
    else:
        L.append(md_table(["항목", "값"],
                          [["파일", mas["path"]],
                           ["accession 컬럼", mas.get("accession_column")],
                           ["cluster 컬럼", mas.get("cluster_column")],
                           ["행 수", mas.get("n_rows")],
                           ["클러스터 수", mas.get("n_clusters")],
                           ["기대 클러스터 수", EXPECTED_MASTER_CLUSTERS],
                           ["universe 와 공통 accession",
                            f"{mas.get('n_common_universe')} "
                            f"({fnum(mas.get('join_rate_universe'))})"]]) + "\n")
        best = metrics.loc[metrics["ari_vs_master"].idxmax()] if metrics["ari_vs_master"].notna().any() else None
        if best is not None:
            L.append(f"- master 분할과 ARI 가 가장 높은 임계값은 "
                     f"**{fnum(float(best['min_seq_id']), 2)}** (ARI {fnum(float(best['ari_vs_master']))}).")
            ref_row = metrics[np.isclose(metrics["min_seq_id"], REFERENCE_THRESHOLD)]
            if len(ref_row):
                v = float(ref_row.iloc[0]["ari_vs_master"])
                L.append(f"- 0.70 재실행 vs master ARI = **{fnum(v)}** "
                         f"(1.000 이면 2026-04 분할이 그대로 재현된 것이다. "
                         f"대표 ID 가 달라도 ARI 는 1.000 이 될 수 있다 — 인계문서 7절).")
        L.append("")

    L.append("## 7. 숫자 점검\n")
    u, a = S["universe"], S["assignments"]
    L.append(md_table(
        ["항목", "값"],
        [["universe accession 수", u["n_accessions"]],
         ["기대 서열 수(--expect-n)", S["expect_n"] or "검사 안 함"],
         ["sweep_assignments.tsv 행 수", a["n_rows"]],
         ["임계값 열 수", a["n_threshold_columns"]],
         ["전 임계값에 배정된 행", a["n_rows_fully_assigned"]],
         ["일부만 배정된 행", a["n_rows_partially_assigned"]],
         ["한 번도 배정되지 않은 행", a["n_rows_never_assigned"]],
         ["중복 accession 배정 행 합", sum(v["n_duplicate_accession_rows"]
                                          for v in S["partition_meta"].values())],
         ["universe 밖 accession 합", sum(v["n_outside_universe"]
                                         for v in a["per_label"].values())]]) + "\n")
    L.append("상한에 걸려 잘린 것\n")
    L.append(md_table(["항목", "상한", "잘림"],
                      [["ChCODH2 멤버 목록", S["limits"]["max_members"],
                        "예" if (len(ref_track) and bool(ref_track["members_truncated"].any())) else "아니오"],
                       ["합병 출처 클러스터 목록", S["limits"]["max_list"],
                        "예" if (len(adjacent) and bool(adjacent.get("ref_merged_truncated",
                                                                     pd.Series(dtype=bool)).any())) else "아니오"],
                       ["마커 동거 목록", S["limits"]["max_list"],
                        "예" if (len(markers_df) and bool(markers_df["co_markers_truncated"].any())) else "아니오"],
                       ["md 표 행 수", S["limits"]["max_table_rows"], "표마다 위에 명시"]]) + "\n")

    L.append("## 8. 불일치\n")
    if S["discrepancies"]:
        for d in S["discrepancies"]:
            L.append(f"- {d}")
        L.append("")
    else:
        L.append("- 없음\n")

    L.append("## 9. 한계\n")
    L.append("- ARI/NMI 는 공통 accession 위에서만 계산한다. `n_common` 이 universe 보다 "
             "작으면 그 비교는 부분집합에 대한 것이다 (표에 수치로 적혀 있다).")
    L.append("- 안정/경계 판정은 `--stable-ari` 한 값에 의존한다. 그 값이 params.json 에 "
             "기록되어 있으니 바꿔서 재실행해 결론이 유지되는지 확인할 것.")
    L.append("- `-c`/`--cov-mode`/`-s` 를 바꾸면 219 와의 비교 가능성이 깨진다. "
             "이번 실행값은 1절에 적혀 있다.")
    L.append("- mmseqs 가 어느 멤버를 대표로 고르는지는 버전에 따라 다르다. 대표 ID 가 바뀌어도 "
             "분할이 같으면 ARI 는 1.000 이다 (인계문서 7절).")
    Path(path).write_text("\n".join(L) + "\n")


# ==========================================================================
# 본체
# ==========================================================================
def run(args) -> dict:
    t_start = datetime.now(timezone.utc)
    fasta = Path(args.fasta)
    master_path = Path(args.master)
    outdir = Path(args.outdir)
    sweep_dir = outdir / "sweep"
    logs_dir = outdir / "logs"
    tmp_root = Path(args.tmp_root) if args.tmp_root else (outdir / "tmp")

    discrepancies: list[str] = []

    def note(msg: str) -> None:
        discrepancies.append(msg)
        print(f"[불일치] {msg}")

    common.verify_known_positions()     # 번호 체계 점검 (매 실행, 단계 1과 동일)

    thresholds, thr_meta = parse_thresholds(
        args.thresholds, force_reference=not args.no_force_reference_threshold,
        reference=REFERENCE_THRESHOLD)
    order = thr_meta["labels"]
    if thr_meta["reference_threshold_added"]:
        print(f"[주의] 운용 임계값 {REFERENCE_THRESHOLD:.2f} 가 --thresholds 에 없어 자동으로 추가했다 "
              f"(--no-force-reference-threshold 로 끌 수 있다).")
    if thr_meta["n_duplicates_dropped"]:
        print(f"[주의] 중복 임계값 {thr_meta['n_duplicates_dropped']}개를 제거했다.")

    # ---------------- dry-run ----------------
    if args.dry_run:
        print("=== --dry-run : mmseqs 를 실행하지 않고 돌릴 명령만 출력한다 (아무것도 쓰지 않는다) ===")
        try:
            mmseqs = common.find_mmseqs()
            mm_note = "탐색 성공"
        except FileNotFoundError as e:
            mmseqs = "mmseqs"
            mm_note = f"못 찾음 → 'mmseqs' 로 표기 ({e})"
        print(f"mmseqs      {mmseqs}  [{mm_note}]")
        for label, p in [("FASTA", fasta), ("master", master_path),
                         ("markers", Path(args.markers) if args.markers else None),
                         ("참조 FASTA", common.REPO_LIB / "chcodh2.fasta")]:
            if p is None:
                print(f"  {label:10s} 지정 안 됨 (해당 분석은 건너뛴다)")
            else:
                print(f"  {label:10s} {'있음' if p.exists() else '없음'}  {p}")
        print(f"고정 옵션   -c {args.cov} --cov-mode {args.cov_mode} -s {args.sensitivity} "
              f"--threads {args.threads}")
        print(f"임계값      {len(thresholds)}개 {thresholds}")
        print(f"tmp 루트    {tmp_root}  (임계값마다 생성 후 "
              f"{'보존' if args.keep_tmp else '삭제'})")
        print("--- 실행 예정 명령 ---")
        for t in thresholds:
            label = threshold_label(t)
            cmd = build_cluster_cmd(mmseqs, fasta, sweep_dir / label, tmp_root / label, t,
                                    cov=args.cov, cov_mode=args.cov_mode,
                                    sensitivity=args.sensitivity, threads=args.threads)
            print("  $ " + " ".join(cmd))
        print("--- 산출 예정 파일 ---")
        for t in thresholds:
            print(f"  {sweep_dir / (threshold_label(t) + '_cluster.tsv')}")
        for name in ["sweep_metrics.tsv", "sweep_adjacent.tsv", "sweep_assignments.tsv",
                     "sweep_reference_track.tsv", "sweep_markers.tsv", "sweep_summary.md",
                     "params.json"]:
            print(f"  {outdir / name}")
        print(f"  {sweep_dir / 'params.json'}")
        print(f"  {logs_dir / 'id<XX>.log'}")
        print(f"와이드 열 이름  {[column_name(lb) for lb in order]}")
        print(f"안정 판정       인접 ARI ≥ {args.stable_ari}")
        print(f"참조 서열 {common.CHCODH2_LEN} aa · 논문 위치 "
              f"{len(common.KNOWN_POSITIONS)}개 번호 체계 점검 통과")
        return dict(dry_run=True, thresholds=thresholds,
                    commands=[" ".join(build_cluster_cmd(
                        mmseqs, fasta, sweep_dir / threshold_label(t), tmp_root / threshold_label(t),
                        t, cov=args.cov, cov_mode=args.cov_mode,
                        sensitivity=args.sensitivity, threads=args.threads)) for t in thresholds])

    # ---------------- 클러스터링 또는 기존 산출물 ----------------
    mm_info = dict(binary=None, version=None, threads=args.threads)
    if args.from_existing:
        mode = "from-existing"
        print(f"=== 단계 2 클러스터 스윕 ({mode}) — mmseqs 를 부르지 않는다 ===")
        explicit = args.thresholds if args.thresholds_explicit else None
        run_records, existing_meta = discover_existing(Path(args.from_existing), explicit)
        order = [r["label"] for r in run_records]
        thresholds = [r["min_seq_id"] for r in run_records]
        thr_meta.update(thresholds=thresholds, labels=order, from_existing=existing_meta)
        print(f"cluster.tsv {existing_meta['n_files_found']}개 발견 "
              f"({existing_meta['search_dir']}) → 사용 {len(run_records)}개")
        if existing_meta["n_files_unparsable"]:
            note(f"라벨을 해석할 수 없는 파일 {existing_meta['n_files_unparsable']}개를 무시했다: "
                 f"{existing_meta['unparsable_names']}")
        if existing_meta.get("n_missing"):
            note(f"--thresholds 로 요청했으나 파일이 없는 임계값 {existing_meta['n_missing']}개: "
                 f"{existing_meta['missing_thresholds']}")
        if existing_meta.get("n_ignored_extra"):
            print(f"[주의] --thresholds 밖의 파일 {existing_meta['n_ignored_extra']}개는 무시했다.")
        mm_info.update(binary="(사용 안 함 — from-existing)", version="(사용 안 함)")
    else:
        mode = "mmseqs"
        if not fasta.exists():
            raise FileNotFoundError(
                f"입력 FASTA 가 없다: {fasta}. --fasta 로 지정하거나 --dry-run / "
                f"--from-existing 을 쓸 것.")
        mmseqs = args.mmseqs or common.find_mmseqs()
        version = common.mmseqs_version(mmseqs)
        mm_info.update(binary=str(mmseqs), version=version)
        print(f"=== 단계 2 클러스터 스윕 ({mode}) ===")
        print(f"mmseqs      {mmseqs}  [{version}]")
        print(f"입력 FASTA  {fasta}")
        print(f"고정 옵션   -c {args.cov} --cov-mode {args.cov_mode} -s {args.sensitivity} "
              f"--threads {args.threads}")
        if (args.cov, args.cov_mode, args.sensitivity) != (FIXED_COV, FIXED_COV_MODE, FIXED_SENSITIVITY):
            note(f"고정 옵션이 인계문서 1절 값(-c {FIXED_COV} --cov-mode {FIXED_COV_MODE} "
                 f"-s {FIXED_SENSITIVITY})과 다르다. 219 분할과의 비교 가능성이 깨진다.")
        print(f"임계값      {len(thresholds)}개 {thresholds}")
        run_records = run_sweep(args, mmseqs, thresholds, sweep_dir, logs_dir, tmp_root)

    # ---------------- cluster.tsv 읽기 ----------------
    partitions: dict[str, pd.Series] = {}
    file_meta, part_meta = {}, {}
    for rec in run_records:
        label = rec["label"]
        df, fmeta = read_cluster_tsv(Path(rec["cluster_tsv"]))
        ser, pmeta = partition_from_cluster_df(
            df, label, allow_duplicate_accession=args.allow_duplicate_accession)
        partitions[label], file_meta[label], part_meta[label] = ser, fmeta, pmeta
        if fmeta["header_row_dropped"]:
            note(f"{label}: cluster.tsv 첫 행이 헤더처럼 보여 제외했다 "
                 f"(easy-cluster 산출물은 헤더가 없다).")
        if pmeta["n_conflicting_accession"]:
            note(f"{label}: accession {pmeta['n_conflicting_accession']}개가 여러 클러스터에 "
                 f"걸쳐 있어 정렬 후 첫 배정만 썼다 (--allow-duplicate-accession).")

    # ---------------- universe ----------------
    if fasta.exists():
        universe, uni_meta = fasta_universe(fasta)
        if uni_meta["n_empty_accession"]:
            note(f"FASTA 헤더에서 accession 을 뽑지 못한 레코드 {uni_meta['n_empty_accession']}개 "
                 f"— 그 행들은 universe 에서 빠진다.")
        if uni_meta["n_duplicate_accession"]:
            note(f"FASTA 에 중복 accession {uni_meta['n_duplicate_accession']}개 "
                 f"(예: {uni_meta['duplicate_examples'][:5]}) — universe 는 유일 accession 기준이다.")
    else:
        all_acc = set()
        for ser in partitions.values():
            all_acc |= set(ser.index.tolist())
        universe = sorted(all_acc)
        uni_meta = dict(source="union(cluster.tsv)", path=None, n_records=len(universe),
                        n_unique_accession=len(universe), n_empty_accession=0,
                        n_duplicate_accession=0, duplicate_examples=[])
        print(f"[주의] FASTA 가 없어 universe 를 cluster.tsv {len(partitions)}개의 "
              f"accession 합집합({len(universe)}개)으로 잡았다.")
    uni_meta["n_accessions"] = len(universe)
    print(f"universe    {len(universe)}개 accession (출처 {uni_meta['source']})")

    if args.expect_n and len(universe) != args.expect_n:
        note(f"universe 가 {len(universe)}개 — 기대값 {args.expect_n} 과 다르다 "
             f"(--expect-n 0 으로 검사를 끌 수 있다).")

    # easy-cluster 는 모든 서열을 배정해야 한다. 빠진 것은 기본적으로 예외.
    uni_idx = pd.Index(universe)
    for label in order:
        miss = int(len(uni_idx.difference(partitions[label].index)))
        if miss:
            msg = (f"{label}: universe {len(universe)}개 중 {miss}개가 cluster.tsv 에 없다. "
                   f"easy-cluster 는 모든 입력 서열을 배정하므로 입력이 다른 파일일 가능성이 높다.")
            if args.allow_incomplete:
                note(msg + " (--allow-incomplete 로 계속하며 해당 칸은 빈 문자열로 둔다)")
            else:
                raise ValueError(msg + " 의도한 것이면 --allow-incomplete 를 줄 것.")

    # ---------------- master 참조 분할 ----------------
    master, master_info = read_master_partition(
        master_path, allow_duplicate_accession=args.allow_duplicate_accession)
    if master is None:
        print(f"[주의] master 를 읽지 않았다 ({master_info['reason']}) — 참조 분할 비교를 건너뛴다.")
    else:
        n_common = int(len(uni_idx.intersection(master.index)))
        master_info.update(n_common_universe=n_common,
                           join_rate_universe=float(n_common / len(universe)) if universe else float("nan"),
                           n_master_outside_universe=int(len(master.index.difference(uni_idx))))
        print(f"master      {master_info['n_rows']}행 · 클러스터 {master_info['n_clusters']}개 · "
              f"universe 와 공통 {n_common}개 "
              f"({master_info['join_rate_universe'] * 100:.2f}%)")
        if master_info["n_clusters"] != args.expect_master_clusters and args.expect_master_clusters:
            note(f"master 클러스터 수가 {master_info['n_clusters']} — 기대값 "
                 f"{args.expect_master_clusters}(인계문서 1절) 와 다르다.")
        if n_common == 0:
            note("master 와 universe 의 공통 accession 이 0 — accession 정규화 규칙이 다를 수 있다. "
                 "ARI vs master 는 전부 NaN 이 된다.")
        if master_info["join_rate_universe"] < args.min_join_rate:
            raise ValueError(
                f"master 조인율 {master_info['join_rate_universe']:.4f} < --min-join-rate "
                f"{args.min_join_rate}. 참조 비교가 부분집합에만 성립하므로 중단한다.")

    # ---------------- 마커 ----------------
    if args.markers:
        mpath = Path(args.markers)
        if not mpath.exists():
            raise FileNotFoundError(f"--markers 파일이 없다: {mpath}")
        markers, marker_file_meta = read_markers(mpath)
        if not markers:
            note(f"마커 파일 {mpath} 에서 accession 을 하나도 읽지 못했다 "
                 f"(전부 주석일 수 있다) — 동거 추적을 건너뛴다.")
        if marker_file_meta["n_unparsable_lines"]:
            note(f"마커 파일에서 해석 못한 줄 {marker_file_meta['n_unparsable_lines']}개 "
                 f"(예: {marker_file_meta['unparsable_examples'][:2]})")
    else:
        markers, marker_file_meta = [], dict(skipped=True, reason="--markers 지정 안 됨")
        print("[주의] --markers 가 없어 대조군 동거 추적을 건너뛴다.")

    # ---------------- (d) ChCODH2 ----------------
    ref_acc = safe_parse_accession(args.reference_acc)
    ref_track, ref_rows, ref_info = track_reference(partitions, order, ref_acc,
                                                    max_members=args.max_members)
    if not ref_info["found_in_any_threshold"]:
        msg = (f"참조 {ref_acc} 가 어느 임계값의 cluster.tsv 에도 없다 — ChCODH2 클러스터 추적을 "
               f"건너뛴다. --reference-acc 로 실제 accession 을 지정할 것.")
        if args.require_reference:
            raise ValueError(msg + " (--require-reference 지정됨)")
        note(msg)
    else:
        print(f"ChCODH2     {ref_acc} — {ref_info['n_thresholds_found']}/{len(order)} 임계값에서 발견, "
              f"클러스터 크기 {ref_info['size_min']}~{ref_info['size_max']}")

    # ---------------- (e) 마커 ----------------
    if markers:
        markers_df, marker_rows, marker_info = track_markers(partitions, order, markers,
                                                             max_list=args.max_list)
        marker_info.update(marker_file_meta)
        marker_info["skipped"] = False
        if marker_info["missing_markers"]:
            note(f"universe 에 없는 마커 {len(marker_info['missing_markers'])}개: "
                 f"{marker_info['missing_markers'][:10]}")
        print(f"마커        {marker_info['n_markers']}개 중 "
              f"{marker_info['n_markers_in_universe']}개 발견 · 쌍 {marker_info['n_pairs']}개")
    else:
        markers_df = pd.DataFrame(columns=["label", "min_seq_id", "marker_accession",
                                           "marker_label", "found", "cluster", "cluster_size",
                                           "co_marker_group", "n_co_markers", "co_markers",
                                           "co_markers_truncated"])
        marker_rows = {}
        marker_info = dict(marker_file_meta, skipped=True,
                           reason=marker_file_meta.get("reason", "마커 없음"),
                           n_markers=0, n_markers_in_universe=0, n_pairs=0,
                           missing_markers=[], pairs=[])

    # ---------------- 표 만들기 ----------------
    wide, assign_meta = build_assignments(partitions, universe, order)
    if len(wide) != len(universe):
        raise AssertionError(f"와이드 테이블 행 수가 {len(wide)} — universe {len(universe)} 와 다르다.")
    if assign_meta["n_rows_never_assigned"]:
        note(f"어느 임계값에서도 배정되지 않은 accession {assign_meta['n_rows_never_assigned']}개")

    metrics = metrics_table(partitions, order, run_records, file_meta, part_meta,
                            assign_meta, master, ref_rows, marker_rows)
    adjacent, interp = adjacent_table(partitions, order, stable_ari=args.stable_ari,
                                      ref_acc=ref_acc, max_list=args.max_list)

    # master 쪽 cluster 열을 와이드 테이블에 함께 둔다 (뒤 단계의 조인 편의)
    if master is not None:
        wide.insert(1, "master_cluster", master.reindex(pd.Index(universe)).fillna("").astype(str))
    wide.insert(1, "is_reference", [a == ref_acc for a in universe])
    marker_set = {m["marker_accession"] for m in markers}
    wide.insert(2, "is_marker", [a in marker_set for a in universe])

    print(f"지표        임계값 {len(metrics)}행 · 인접 쌍 {len(adjacent)}행 · "
          f"와이드 {len(wide)}행 × {assign_meta['n_threshold_columns']}열")
    if len(adjacent):
        print(f"안정성      인접 ARI 중앙 {interp['ari_median']:.4f} · 최소 {interp['ari_min']:.4f} · "
              f"안정 구간 {len(interp['stable_intervals'])}개 · 경계 {len(interp['boundaries'])}개")
        hit = interp["reference_in_stable_interval"]
        if hit:
            print(f"            0.70 은 안정 구간 {hit['low']:.2f}–{hit['high']:.2f} 안에 있다 "
                  f"(최소 ARI {hit['min_ari']:.4f})")
        elif interp.get("reference_threshold_present"):
            print("            0.70 이 안정 구간 밖이다 — 요약 md 3절 해석 참조")

    # ---------------- 쓰기 ----------------
    outdir.mkdir(parents=True, exist_ok=True)
    paths = dict(
        metrics=outdir / "sweep_metrics.tsv",
        adjacent=outdir / "sweep_adjacent.tsv",
        assignments=outdir / "sweep_assignments.tsv",
        reference_track=outdir / "sweep_reference_track.tsv",
        markers=outdir / "sweep_markers.tsv",
        summary=outdir / "sweep_summary.md",
    )
    metrics.to_csv(paths["metrics"], sep="\t", index=False)
    adjacent.to_csv(paths["adjacent"], sep="\t", index=False)
    wide.reset_index().to_csv(paths["assignments"], sep="\t", index=False)
    ref_track.to_csv(paths["reference_track"], sep="\t", index=False)
    markers_df.to_csv(paths["markers"], sep="\t", index=False)

    total_seconds = [r.get("seconds") for r in run_records if r.get("seconds") is not None]
    stats = dict(
        mode=mode,
        input_paths=dict(fasta=str(fasta), master=str(master_path),
                         markers=str(args.markers) if args.markers else None,
                         from_existing=str(args.from_existing) if args.from_existing else None,
                         reference_fasta=str(common.REPO_LIB / "chcodh2.fasta")),
        output_paths={k: str(v) for k, v in paths.items()},
        mmseqs=mm_info,
        fixed_options=dict(cov=args.cov, cov_mode=args.cov_mode, sensitivity=args.sensitivity,
                           handoff_cov=FIXED_COV, handoff_cov_mode=FIXED_COV_MODE,
                           handoff_sensitivity=FIXED_SENSITIVITY),
        thresholds=thr_meta,
        expect_n=args.expect_n,
        expect_master_clusters=args.expect_master_clusters,
        limits=dict(max_members=args.max_members, max_list=args.max_list,
                    max_table_rows=args.max_table_rows),
        universe=uni_meta,
        file_meta=file_meta,
        partition_meta=part_meta,
        assignments=assign_meta,
        master=master_info,
        reference=ref_info,
        markers=marker_info,
        adjacent_interpretation=interp,
        metrics_records=metrics.replace({np.nan: None}).to_dict(orient="records"),
        adjacent_records=adjacent.replace({np.nan: None}).to_dict(orient="records"),
        runs=run_records,
        timing=dict(total_mmseqs_seconds=round(sum(total_seconds), 2) if total_seconds else None,
                    n_mmseqs_runs=len(total_seconds),
                    wall_seconds=round((datetime.now(timezone.utc) - t_start).total_seconds(), 2)),
        tmp=dict(root=str(tmp_root), kept=bool(args.keep_tmp),
                 bytes_per_run={r["label"]: r.get("tmp_bytes") for r in run_records}),
        discrepancies=discrepancies,
        column_naming=dict(wide_columns=[column_name(lb) for lb in order],
                           is_mm_rep_rule="is_mm_rep@XX = (cluster@XX == accession)"),
    )
    write_summary(paths["summary"], metrics, adjacent, ref_track, markers_df, stats)

    inputs = [p for p in [fasta, master_path,
                          Path(args.markers) if args.markers else None,
                          common.REPO_LIB / "chcodh2.fasta"] if p and Path(p).exists()]
    params = common.write_params(
        outdir, STAGE, inputs=inputs, outputs=list(paths.values()),
        commands=[" ".join(sys.argv)] + [r["cmd_str"] for r in run_records if r.get("cmd_str")],
        **stats)

    if mode == "mmseqs" and sweep_dir.exists():
        raw_outputs = sorted(sweep_dir.glob("id*_cluster.tsv"))
        common.write_params(
            sweep_dir, STAGE + ".raw", inputs=[fasta], outputs=raw_outputs,
            commands=[r["cmd_str"] for r in run_records if r.get("cmd_str")],
            mmseqs=mm_info, fixed_options=stats["fixed_options"], thresholds=thr_meta,
            runs=run_records, timing=stats["timing"], tmp=stats["tmp"],
            note="easy-cluster 원본 산출물. 가공하지 않는다.")

    print("출력  " + "\n      ".join(str(v) for v in paths.values()))
    print(f"      {params}")
    if discrepancies:
        print(f"[요약] 불일치 {len(discrepancies)}건 — params.json 과 요약 md 8절 참조")
    return dict(stats, metrics=metrics, adjacent=adjacent, assignments=wide,
                reference_track=ref_track, markers_table=markers_df,
                partitions=partitions, paths=paths)


# ==========================================================================
# 합성 데이터 자기검증 — mmseqs 없이 --from-existing 전 경로를 실제로 돌린다
# ==========================================================================
# 계획된 분할: 임계값이 내려갈수록 묶음이 커지는 완전 중첩 계층.
#   0.90 / 0.85 / 0.80 → 24 클러스터 (크기 5)      ← 평탄 구간
#   0.75 / 0.70 / 0.65 → 12 클러스터 (크기 10)     ← 평탄 구간 (0.70 이 여기)
#   0.60 / 0.55        →  4 클러스터 (크기 30)
#   0.50               →  1 클러스터 (크기 120)
SYNTH_PLAN = {
    0.90: 24, 0.85: 24, 0.80: 24,
    0.75: 12, 0.70: 12, 0.65: 12,
    0.60: 4, 0.55: 4,
    0.50: 1,
}
SYNTH_N = 120
SYNTH_REF_INDEX = 3          # Q9F8A8 로 바꿔 끼울 accession 위치 (0-based)


def synth_accessions() -> list[str]:
    accs = [f"ACC{i:04d}" for i in range(1, SYNTH_N + 1)]
    accs[SYNTH_REF_INDEX] = common.CHCODH2_ACC
    return accs


def synth_partition(accs: list[str], n_clusters: int) -> list[tuple[str, str]]:
    """중첩 계층을 보장하려면 블록을 '연속 구간' 으로 잘라야 한다."""
    size = SYNTH_N // n_clusters
    rows = []
    for b in range(n_clusters):
        block = accs[b * size:(b + 1) * size]
        rep = block[0]                        # 대표는 블록의 첫 accession
        for a in block:
            rows.append((rep, a))
    return rows


def make_synthetic(dirpath: Path) -> dict:
    """합성 cluster.tsv / master.tsv / markers.txt / FASTA 를 만든다."""
    sweep = dirpath / "sweep"
    sweep.mkdir(parents=True, exist_ok=True)
    accs = synth_accessions()

    for t, k in sorted(SYNTH_PLAN.items()):
        rows = synth_partition(accs, k)
        # mmseqs 산출물처럼 멤버 순서를 섞어도 결과가 같아야 한다
        rows = rows[1:] + rows[:1]
        path = sweep / f"{threshold_label(t)}_cluster.tsv"
        path.write_text("".join(f"{rep}\t{acc}\n" for rep, acc in rows))

    # master = 0.70 과 같은 분할이지만 대표 ID 를 블록의 '마지막' accession 으로 바꿨다.
    # 대표 ID 가 달라도 ARI 는 1.000 이어야 한다 (인계문서 7절).
    size = SYNTH_N // SYNTH_PLAN[0.70]
    mrows = []
    for b in range(SYNTH_PLAN[0.70]):
        block = accs[b * size:(b + 1) * size]
        rep = block[-1]
        for a in block:
            mrows.append((a, rep, int(a == rep)))
    master = dirpath / "master_synth.tsv"
    master.write_text("accession\tcluster\tis_mm_rep\n"
                      + "".join(f"{a}\t{c}\t{r}\n" for a, c, r in mrows))

    # 마커 4개 — 0.70 블록 안에서 2개, 다른 블록에서 2개
    markers = dirpath / "markers_synth.txt"
    markers.write_text(
        "# 합성 마커\n"
        f"{common.CHCODH2_ACC}\tChCODH2 synthetic\n"
        f"{accs[9]}\tsame_block_at_0.70\n"
        f"{accs[10]}\tnext_block_at_0.70\n"
        f"{accs[60]}\tfar_block\n"
        "ZZZZZZ\tuniverse 에 없는 마커\n"
        "#\t주석만 있는 줄\n")

    # FASTA — universe 를 FASTA 에서 잡는 경로도 함께 검증한다
    fasta = dirpath / "synth.fasta"
    fasta.write_text("".join(f">tr|{a}|SYN_{i}\nMKV{'A' * 20}\n" for i, a in enumerate(accs)))

    return dict(sweep=sweep, master=master, markers=markers, fasta=fasta, accessions=accs)


def self_test(dirpath: Path, base_args) -> None:
    dirpath = Path(dirpath)
    if dirpath.exists():
        shutil.rmtree(dirpath)
    dirpath.mkdir(parents=True)
    syn = make_synthetic(dirpath)
    accs = syn["accessions"]
    print(f"=== 자기검증 : 합성 데이터 {dirpath} ===")
    print(f"  cluster.tsv {len(SYNTH_PLAN)}개 · accession {SYNTH_N}개 · "
          f"계획된 클러스터 수 {dict(sorted(SYNTH_PLAN.items()))}")

    outdir = dirpath / "out"
    argv = ["--from-existing", str(dirpath),
            "--fasta", str(syn["fasta"]),
            "--master", str(syn["master"]),
            "--markers", str(syn["markers"]),
            "--outdir", str(outdir),
            "--expect-n", str(SYNTH_N),
            "--expect-master-clusters", str(SYNTH_PLAN[0.70]),
            "--no-force-reference-threshold"]
    print("  $ 02_cluster_sweep.py " + " ".join(argv))
    res = run(build_parser().parse_args(argv))

    checks: list[tuple[str, bool, str]] = []

    def chk(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    m = res["metrics"].set_index("label")
    adj = res["adjacent"]
    wide = res["assignments"]
    interp = res["adjacent_interpretation"]

    # (b) 임계값별 지표
    chk("임계값 행 수", len(m) == len(SYNTH_PLAN), f"{len(m)} vs {len(SYNTH_PLAN)}")
    for t, k in SYNTH_PLAN.items():
        lb = threshold_label(t)
        chk(f"{lb} 클러스터 수 == {k}", int(m.loc[lb, "n_clusters"]) == k,
            f"실제 {m.loc[lb, 'n_clusters']}")
    chk("id50 최대점유 == 1.0", abs(float(m.loc["id50", "largest_frac"]) - 1.0) < 1e-12)
    chk("id90 중앙크기 == 5", float(m.loc["id90", "median_size"]) == 5.0)
    chk("싱글턴 0개 (계획상 없음)", int(m["n_singletons"].sum()) == 0)

    # (c) 인접 안정성 — 평탄 구간은 ARI 1.0, 경계는 1.0 미만
    a = adj.set_index(["label_low", "label_high"])
    chk("id80↔id85 ARI == 1.0", abs(float(a.loc[("id80", "id85"), "ari"]) - 1.0) < 1e-12)
    chk("id85↔id90 ARI == 1.0", abs(float(a.loc[("id85", "id90"), "ari"]) - 1.0) < 1e-12)
    chk("id65↔id70 ARI == 1.0", abs(float(a.loc[("id65", "id70"), "ari"]) - 1.0) < 1e-12)
    chk("id75↔id80 ARI < 1.0 (경계)", float(a.loc[("id75", "id80"), "ari"]) < 1.0,
        f"ARI={a.loc[('id75', 'id80'), 'ari']:.4f}")
    chk("id50↔id55 ARI 가 가장 낮다",
        float(a["ari"].idxmin() == 0 or True) and
        float(a.loc[("id50", "id55"), "ari"]) == float(a["ari"].min()),
        f"min={a['ari'].min():.4f}")
    chk("계층 중첩이라 n_finer_split == 0", int(adj["n_finer_split"].sum()) == 0,
        f"합 {int(adj['n_finer_split'].sum())}")
    chk("경계마다 n_coarser_merged > 0",
        bool((adj.loc[~adj["is_stable"], "n_coarser_merged"] > 0).all()))
    chk("평탄 구간은 one_to_one == 클러스터 수",
        int(a.loc[("id85", "id90"), "one_to_one"]) == SYNTH_PLAN[0.90])

    # 안정 구간 / 경계 해석
    ivs = {(round(iv["low"], 2), round(iv["high"], 2)) for iv in interp["stable_intervals"]}
    chk("안정 구간 0.80–0.90 검출", (0.80, 0.90) in ivs, str(sorted(ivs)))
    chk("안정 구간 0.65–0.75 검출", (0.65, 0.75) in ivs, str(sorted(ivs)))
    chk("0.70 이 안정 구간 안", interp["reference_in_stable_interval"] is not None,
        str(interp["reference_in_stable_interval"]))
    chk("경계가 3개 (0.50/0.55/0.60/0.75 전이)", len(interp["boundaries"]) == 3,
        str([(b["low"], b["high"]) for b in interp["boundaries"]]))

    # (d) ChCODH2 추적
    rt = res["reference_track"].set_index("label")
    chk("ChCODH2 전 임계값에서 발견", bool(rt["ref_found"].all()))
    chk("id90 ChCODH2 클러스터 크기 5", int(rt.loc["id90", "ref_cluster_size"]) == 5)
    chk("id70 ChCODH2 클러스터 크기 10", int(rt.loc["id70", "ref_cluster_size"]) == 10)
    chk("id50 ChCODH2 클러스터 크기 120", int(rt.loc["id50", "ref_cluster_size"]) == SYNTH_N)
    chk("ChCODH2 대표는 ACC0001 (블록 첫 accession)",
        str(rt.loc["id90", "ref_cluster"]) == "ACC0001", str(rt.loc["id90", "ref_cluster"]))
    chk("ChCODH2 멤버 불변 구간 0.80–0.90 검출",
        (0.80, 0.90) in {(round(iv["low"], 2), round(iv["high"], 2))
                         for iv in interp["ref_membership_intervals"]},
        str(interp["ref_membership_intervals"]))
    chk("ChCODH2 합병 사건 3건", len(interp["ref_merge_steps"]) == 3,
        str([(s["high"], s["low"], s["n_source_clusters"]) for s in interp["ref_merge_steps"]]))
    # 0.80 → 0.75 에서 크기 5 클러스터 둘이 하나(크기 10)로 합쳐진다
    step = a.loc[("id75", "id80")]
    chk("0.80→0.75 합병: 상위 클러스터 2개 흡수", int(step["ref_n_source_clusters"]) == 2,
        str(step["ref_n_source_clusters"]))
    chk("0.80→0.75 합병: 서열 5개 추가", int(step["ref_n_added"]) == 5, str(step["ref_n_added"]))
    chk("0.80→0.75 합병 출처가 2개 기록",
        len(str(step["ref_merged_from"]).split(";")) == 2, str(step["ref_merged_from"]))
    chk("ChCODH2 멤버 추가 수가 단조", bool((adj["ref_n_removed"] == 0).all()),
        f"removed 합 {int(adj['ref_n_removed'].sum())}")

    # (e) 마커
    mk = res["markers_table"]
    mi = res["markers"]
    chk("마커 5개 파싱 (주석 줄 제외, universe 에 없는 1개 포함)",
        mi["n_markers"] == 5, f"{mi['n_markers']}")
    chk("universe 에 없는 마커 ZZZZZZ 보고", mi["missing_markers"] == ["ZZZZZZ"],
        str(mi["missing_markers"]))
    r70 = mk[(mk["label"] == "id70") & (mk["marker_accession"] == common.CHCODH2_ACC)].iloc[0]
    chk("id70 에서 ChCODH2 와 ACC0010 동거", accs[9] in str(r70["co_markers"]).split(";"),
        str(r70["co_markers"]))
    chk("id70 에서 ACC0011 은 다른 그룹",
        accs[10] not in str(r70["co_markers"]).split(";"))
    chk("id50 에서 마커 전부 한 묶음",
        str(m.loc["id50", "markers_all_together"]) == "yes",
        str(m.loc["id50", "markers_all_together"]))
    # 0.90 에서는 크기 5 블록이므로 마커 4개가 서로 다른 4개 클러스터에 흩어진다
    chk("id90 마커 클러스터 수 4", int(m.loc["id90", "n_marker_clusters"]) == 4,
        str(m.loc["id90", "n_marker_clusters"]))
    chk("id70 마커 클러스터 수 3", int(m.loc["id70", "n_marker_clusters"]) == 3,
        str(m.loc["id70", "n_marker_clusters"]))
    pair = next(p for p in mi["pairs"]
                if {p["marker_a"], p["marker_b"]} == {common.CHCODH2_ACC, accs[10]})
    chk("ChCODH2–ACC0011 쌍은 0.60 에서 처음 한 묶음",
        abs((pair["highest_together"] or 0) - 0.60) < 1e-9, str(pair))

    # (f) master 대비
    chk("id70 ARI vs master == 1.0 (대표 ID 가 달라도)",
        abs(float(m.loc["id70", "ari_vs_master"]) - 1.0) < 1e-12,
        str(m.loc["id70", "ari_vs_master"]))
    chk("id90 ARI vs master < 1.0", float(m.loc["id90", "ari_vs_master"]) < 1.0)
    chk("master 조인율 1.0", abs(res["master"]["join_rate_universe"] - 1.0) < 1e-12)

    # (g) 와이드 테이블
    chk("와이드 행 수 == universe", len(wide) == SYNTH_N, f"{len(wide)}")
    chk("와이드 열 이름 cluster@70 포함", "cluster@70" in wide.columns, str(list(wide.columns)))
    chk("와이드에 빈 칸 없음", int((wide[[c for c in wide.columns
                                        if c.startswith('cluster@')]] == "").sum().sum()) == 0)
    chk("n_thresholds_assigned 전부 9", bool((wide["n_thresholds_assigned"] == len(SYNTH_PLAN)).all()))
    chk("is_reference 1행", int(wide["is_reference"].sum()) == 1)
    chk("master_cluster 열 존재", "master_cluster" in wide.columns)

    # 산출물 존재
    for k, p in res["paths"].items():
        chk(f"산출물 {k}", Path(p).exists(), str(p))
    chk("params.json", (outdir / "params.json").exists())
    summary_text = (outdir / "sweep_summary.md").read_text()
    chk("요약에 안정 구간 해석 포함", "안정 구간" in summary_text)
    chk("요약에 ChCODH2 추적 포함", "ChCODH2 자신의 클러스터 추적" in summary_text)
    # 합성 입력에는 universe 에 없는 마커 ZZZZZZ 하나만 심어 두었다. 그 1건만 보고되어야 한다.
    chk("불일치는 ZZZZZZ 1건뿐",
        len(res["discrepancies"]) == 1 and "ZZZZZZ" in res["discrepancies"][0],
        str(res["discrepancies"]))
    chk("요약 8절에 불일치가 적혔다", "ZZZZZZ" in summary_text)

    # --- 결측 경로: 한 파일에서 서열을 빼면 기본은 예외, --allow-incomplete 면 공백 ---
    short_dir = dirpath / "incomplete"
    (short_dir / "sweep").mkdir(parents=True)
    for p in sorted((dirpath / "sweep").glob("id*_cluster.tsv")):
        lines = p.read_text().splitlines(keepends=True)
        tgt = short_dir / "sweep" / p.name
        tgt.write_text("".join(lines[:-1]) if p.name == "id70_cluster.tsv" else "".join(lines))
    base = ["--from-existing", str(short_dir), "--fasta", str(syn["fasta"]),
            "--master", str(syn["master"]), "--expect-n", str(SYNTH_N),
            "--no-force-reference-threshold"]
    raised = False
    try:
        run(build_parser().parse_args(base + ["--outdir", str(short_dir / "out_strict")]))
    except ValueError as e:
        raised = "allow-incomplete" in str(e)
    chk("서열 누락 시 기본은 예외", raised)
    res2 = run(build_parser().parse_args(
        base + ["--allow-incomplete", "--outdir", str(short_dir / "out_allow")]))
    chk("--allow-incomplete 면 행 수 유지", len(res2["assignments"]) == SYNTH_N)
    chk("--allow-incomplete 면 빈 칸 1개",
        int((res2["assignments"]["cluster@70"] == "").sum()) == 1)
    chk("--allow-incomplete 면 불일치 보고", len(res2["discrepancies"]) >= 1,
        str(res2["discrepancies"]))

    # --- 라벨 충돌 / 범위 / 중복 accession 가드 ---
    # 라벨 충돌 가드 — 0.7 과 0.7000000001 은 둘 다 'id70' 이 된다 (parse_thresholds 의
    # 6자리 중복 제거를 통과하지 못하는 값이므로 make_labels 를 직접 때린다)
    try:
        make_labels([0.7, 0.7000000001])
        chk("라벨 충돌 예외", False)
    except ValueError as e:
        chk("라벨 충돌 예외", "충돌" in str(e), str(e)[:60])
    chk("라벨 왕복 무손실", make_labels([0.1, 0.7, 0.125, 0.9]) ==
        ["id10", "id70", "id12p5", "id90"], str(make_labels([0.1, 0.7, 0.125, 0.9])))
    try:
        parse_thresholds([70.0], force_reference=False, reference=0.7)
        chk("범위 밖 임계값 예외", False)
    except ValueError as e:
        chk("범위 밖 임계값 예외", "분율" in str(e))
    dup_dir = dirpath / "dupacc"
    (dup_dir / "sweep").mkdir(parents=True)
    (dup_dir / "sweep" / "id70_cluster.tsv").write_text(
        "R1\tACC0001\nR1\tACC0002\nR2\tACC0001\nR2\tACC0003\n")
    try:
        run(build_parser().parse_args(["--from-existing", str(dup_dir),
                                       "--fasta", "/nonexistent.fasta",
                                       "--master", "/nonexistent.tsv",
                                       "--outdir", str(dup_dir / "out"),
                                       "--expect-n", "0",
                                       "--no-force-reference-threshold"]))
        chk("중복 accession 예외", False)
    except ValueError as e:
        chk("중복 accession 예외", "여러 클러스터" in str(e) or "클러스터 둘" in str(e), str(e)[:80])

    # --- dry-run 경로 ---
    dr = run(build_parser().parse_args(
        ["--dry-run", "--fasta", str(syn["fasta"]), "--outdir", str(dirpath / "dry")]))
    chk("dry-run 명령 생성", len(dr["commands"]) == len(DEFAULT_THRESHOLDS),
        f"{len(dr['commands'])}")
    chk("dry-run 은 아무것도 쓰지 않는다", not (dirpath / "dry").exists())
    chk("dry-run 명령에 고정 옵션 포함",
        all(("--cov-mode" in c and "-s 7.5" in c and "--min-seq-id" in c) for c in dr["commands"]))

    n_fail = sum(1 for _, ok, _ in checks if not ok)
    print(f"\n=== 자기검증 결과 : {len(checks) - n_fail}/{len(checks)} 통과 ===")
    for name, ok, detail in checks:
        if not ok:
            print(f"  [실패] {name}  {detail}")
    if n_fail:
        raise AssertionError(f"자기검증 {n_fail}건 실패")
    print("모든 검증 통과. 합성 산출물은 " + str(outdir))


# ==========================================================================
# CLI
# ==========================================================================
def default_threads() -> int:
    """다른 작업(RF2-PPI 등)과 충돌하지 않게 2코어를 남긴다 — 인계문서 4절."""
    return max(1, (os.cpu_count() or 1) - 2)


class _ThresholdsAction(argparse.Action):
    """--thresholds 를 사용자가 직접 줬는지 기록한다 (--from-existing 필터링에 쓴다)."""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, [float(v) for v in values])
        setattr(namespace, "thresholds_explicit", True)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="02_cluster_sweep.py",
        description="단계 2 — mmseqs easy-cluster 의 --min-seq-id 를 스윕해 클러스터 분류의 "
                    "임계값 의존성을 정량화한다 (-c / --cov-mode / -s 는 고정).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="mmseqs 가 없는 환경에서는 --dry-run / --from-existing / --self-test 를 쓸 것.")
    p.add_argument("--fasta", default=str(common.IPR_FASTA), help="입력 IPR010047 FASTA")
    p.add_argument("--master", default=str(common.MASTER_TSV),
                   help="참조 분할 (master_1827.tsv). 없으면 비교를 건너뛴다")
    p.add_argument("--markers", default=None,
                   help="마커 accession 목록 파일 (예: pipeline/lib/markers.txt). "
                        "없으면 대조군 동거 추적을 건너뛴다")
    p.add_argument("--outdir", default=str(common.OUT_ROOT / STAGE), help="산출물 디렉터리")
    p.add_argument("--thresholds", nargs="+", type=float, action=_ThresholdsAction,
                   default=list(DEFAULT_THRESHOLDS),
                   help=f"--min-seq-id 스윕 값 (기본 {DEFAULT_THRESHOLDS})")
    p.add_argument("--no-force-reference-threshold", action="store_true",
                   help=f"운용 임계값 {REFERENCE_THRESHOLD:.2f} 자동 추가를 끈다")
    p.add_argument("--mmseqs", default=None,
                   help="mmseqs 절대경로. 기본은 common.find_mmseqs() (인계문서 6절)")
    p.add_argument("--threads", type=int, default=default_threads(),
                   help=f"mmseqs --threads (기본 nproc-2 = {default_threads()}). "
                        f"인계문서 4절 — 다른 작업과 충돌을 피한다")
    p.add_argument("--cov", default=FIXED_COV, help=f"-c (고정값 {FIXED_COV}. 바꾸면 219 와 비교 불가)")
    p.add_argument("--cov-mode", default=FIXED_COV_MODE, help=f"--cov-mode (고정값 {FIXED_COV_MODE})")
    p.add_argument("--sensitivity", default=FIXED_SENSITIVITY,
                   help=f"-s (고정값 {FIXED_SENSITIVITY})")
    p.add_argument("--tmp-root", default=None,
                   help="easy-cluster tmp 의 부모 (기본 <outdir>/tmp). 수 GB 가 될 수 있다")
    p.add_argument("--keep-tmp", action="store_true", help="실행 후 tmp 를 지우지 않는다")
    p.add_argument("--force", action="store_true", help="tmp 루트에 잔여물이 있으면 지우고 진행")
    p.add_argument("--resume", action="store_true",
                   help="이미 있는 id<XX>_cluster.tsv 는 다시 돌리지 않는다")
    p.add_argument("--from-existing", default=None, metavar="DIR",
                   help="DIR (또는 DIR/sweep) 의 id<XX>_cluster.tsv 로 지표만 다시 계산한다. "
                        "mmseqs 를 부르지 않는다")
    p.add_argument("--reference-acc", default=common.CHCODH2_ACC,
                   help=f"추적할 참조 accession (기본 {common.CHCODH2_ACC})")
    p.add_argument("--require-reference", action="store_true",
                   help="참조 accession 이 없으면 예외 (기본은 보고 후 건너뜀)")
    p.add_argument("--stable-ari", type=float, default=DEFAULT_STABLE_ARI,
                   help=f"인접 ARI 가 이 값 이상이면 '안정' (기본 {DEFAULT_STABLE_ARI})")
    p.add_argument("--expect-n", type=int, default=EXPECTED_N_SEQUENCES,
                   help=f"기대 서열 수. 다르면 불일치로 보고 (기본 {EXPECTED_N_SEQUENCES}, 0 이면 검사 안 함)")
    p.add_argument("--expect-master-clusters", type=int, default=EXPECTED_MASTER_CLUSTERS,
                   help=f"master 의 기대 클러스터 수 (기본 {EXPECTED_MASTER_CLUSTERS}, 0 이면 검사 안 함)")
    p.add_argument("--min-join-rate", type=float, default=0.0,
                   help="master 조인율이 이 값 미만이면 예외 (기본 0.0 = 보고만)")
    p.add_argument("--allow-incomplete", action="store_true",
                   help="cluster.tsv 가 universe 의 일부만 담아도 계속한다 (기본은 예외)")
    p.add_argument("--allow-duplicate-accession", action="store_true",
                   help="한 accession 이 여러 클러스터에 걸쳐도 첫 배정만 쓰고 계속한다")
    p.add_argument("--max-members", type=int, default=DEFAULT_MAX_MEMBERS,
                   help=f"ChCODH2 멤버 목록 상한 (기본 {DEFAULT_MAX_MEMBERS}, 잘리면 명시)")
    p.add_argument("--max-list", type=int, default=DEFAULT_MAX_LIST,
                   help=f"합병 출처·동거 마커 목록 상한 (기본 {DEFAULT_MAX_LIST})")
    p.add_argument("--max-table-rows", type=int, default=DEFAULT_MAX_TABLE_ROWS,
                   help=f"요약 md 표의 행 수 상한 (기본 {DEFAULT_MAX_TABLE_ROWS})")
    p.add_argument("--strict", action="store_true",
                   help="불일치가 있으면 산출물을 모두 쓴 뒤 종료코드 2 로 끝낸다")
    p.add_argument("--dry-run", action="store_true",
                   help="mmseqs 를 실행하지 않고 돌릴 명령만 출력한다 (아무것도 쓰지 않는다)")
    p.add_argument("--self-test", metavar="DIR",
                   help="DIR 에 합성 cluster.tsv 를 만들어 --from-existing 전 경로를 실제로 돌려 검증")
    p.set_defaults(thresholds_explicit=False)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test:
        self_test(Path(args.self_test), args)
        return 0
    res = run(args)
    if args.strict and res.get("discrepancies"):
        print(f"[--strict] 불일치 {len(res['discrepancies'])}건 — 종료코드 2")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
