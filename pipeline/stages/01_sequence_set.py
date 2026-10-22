#!/usr/bin/env python3
"""단계 1 — 서열 세트 정의. IPR010047 FASTA 에 QC 플래그를 붙인 단일 기준표를 만든다.

무엇을 읽는가
  1. IPR010047 FASTA           기본값 common.IPR_FASTA
       /mnt/af2results/mingyu/chcodh2/protein-matching-IPR010047.fasta  (1,827 서열)
  2. master_1827.tsv (선택)    기본값 common.MASTER_TSV
       2026-04 mmseqs easy-cluster(--min-seq-id 0.70 -c 0.8 --cov-mode 1 -s 7.5) 결과인
       cluster / is_mm_rep 컬럼을 accession 으로 조인해 가져온다. 없으면 그 두 컬럼은 빈값.
  3. pipeline/lib/chcodh2.fasta  참조 ChCODH2(Q9F8A8, 636 aa) — common.chcodh2_seq()

무엇을 쓰는가  (--outdir, 기본값 common.OUT_ROOT/01_sequence_set)
  1. sequence_set.tsv          한 행 = FASTA 레코드 하나. 입력 레코드 수와 항상 같다.
  2. sequence_set_summary.md   길이 분포 · 플래그 집계 · 중복 · 조인율 · 참조 확인 표
  3. params.json               common.write_params() — 입력 해시, 행 수, 모든 집계값, 임계값

무엇을 보장하는가
  1. 서열을 절대 삭제하지 않는다. QC 는 불리언/범주 플래그 컬럼으로만 표현한다.
     출력 행 수 == FASTA 레코드 수. --expect-n(기본 1827)과 다르면 예외를 던진다.
  2. master 조인은 반드시 행 수를 보존한다(many_to_one 검증). master 쪽 accession 이
     중복이면 예외 — 조용히 행이 불어나는 것을 막는다.
  3. 조인율·결측률·중복 수를 stdout 과 params.json 과 요약 md 세 곳에 모두 적는다.
     조인율이 100% 가 아니면 매칭 실패 accession 예시를 출력한다(기본 10개, 잘림을 명시).
  4. 참조 Q9F8A8 이 FASTA 에 있으면 common.chcodh2_seq() 와 전문 대조한다. 불일치면
     예외 — 잔기 번호 체계가 어긋났다는 뜻이므로 하위 단계 전체가 무효이기 때문이다.
     (FASTA 에 없는 것은 예외가 아니다. --require-reference 로 강제할 수 있다.)
  5. 위치 번호는 전부 ChCODH2 서열 번호(Met1 = 1, 636 aa). 이 단계는 위치를 직접 다루지
     않지만 common.verify_known_positions() 로 번호 체계를 매 실행 점검한다.
  6. 아래 네 가지는 '조용히 틀린 숫자'로 끝나므로 전부 예외다.
     a. FASTA 첫 줄이 '>' 가 아니거나 UTF-8 BOM 으로 시작  — common.read_fasta 가 첫 '>' 앞의
        줄을 조용히 버려 레코드가 통째로 사라진다(1827 -> 1826). assert_fasta_header_first().
     b. --len-core-min > --len-core-max  — core 밴드가 사라져 qc_pass 가 전부 False 가 된다.
     c. --min-join-rate > 0 인데 master 파일이 없음  — 조인을 안 하면 게이트가 조용히 통과된다.
     d. master 의 is_mm_rep 컬럼이 참 토큰 0개인데 (cluster == accession) 인 행은 있음
        — 토큰 어휘 불일치로 is_mm_rep 이 전부 False 가 된다. 실제 값 분포를 예외에 적는다.
  7. TSV 안전성. 탭/개행/수직탭이 든 **모든 문자열 셀**을 공백으로 치환하고 치환한 셀 수와
     행 수를 보고한다. raw_header 만 막으면 source_organism(원본 헤더에서 뽑으므로 탭이
     남을 수 있다) 때문에 cut/awk 에서 컬럼이 한 칸 밀린다.

컬럼 정의 (출력 순서대로)
  accession          common.parse_accession(header). master 조인 키.
  raw_header         '>' 를 뗀 원본 헤더. 탭은 공백으로 치환(TSV 보호, 치환 건수 보고).
  seq_len            아미노산 수 (대소문자 무관, 공백 제거 후)
  n_X                'X' 개수 (대문자 변환 후)
  n_ambiguous        X/B/Z/J/U/O 개수 합
  n_cys              'C' 개수 — 금속 클러스터 리간드 후보 수의 1차 지표
  n_lowercase        소문자 문자 수. 0 이 아니면 소문자 마스킹(저복잡도 영역)이 섞인 입력이다.
                     ※ n_X / n_ambiguous / n_cys / nonstandard_chars / seq_sha256 은 모두
                       **대문자로 변환한 서열**에 대해 계산한다. 따라서 소문자 마스킹만
                       다른 서열은 완전 중복으로 묶이고, 소문자 자체는 qc_nonstandard 를
                       켜지 않는다. 마스킹 여부는 이 컬럼으로만 판단하라.
  nonstandard_chars  20종 표준 아미노산 밖의 문자를 정렬해 이어붙인 문자열 (없으면 빈값)
  seq_sha256         대문자 변환 서열의 sha256 (완전 중복 판정의 근거)
  qc_len_band        'short' (<620) / 'core' (620~700) / 'long' (>700)   ← Kim 2022 Nat Catal
                     경계값은 --len-core-min/--len-core-max 로 바꿀 수 있다.
  qc_has_ambiguous   X/B/Z/J/U/O 를 하나라도 포함
  qc_nonstandard     20종 표준 아미노산 밖 문자를 하나라도 포함
                     ※ 정의상 qc_nonstandard ⊇ qc_has_ambiguous 다. 어떤 문자였는지는
                       nonstandard_chars 를 보라.
  is_putative_fusion seq_len > 700. CdhA/ACS 융합 의심. **길이 기반 추정일 뿐이며 도메인
                     조성을 확인한 것이 아니다.** 단계 2 이후에서 도메인으로 재판정할 것.
  dup_group          완전히 동일한 서열(대문자 변환 후 일치)끼리 묶은 ID 'dup_0001'.
                     중복이 없는 서열은 빈값.
  dup_group_size     그 그룹의 행 수 (중복 없으면 1)
  is_dup_rep         dup_group 대표 여부. **중복이 없는 서열도 True** 이므로
                     df[df.is_dup_rep] 가 곧 유일 서열 집합이다.
                     대표 선정 규칙: is_reference 우선 → accession 사전순 → 등장 순서.
  is_dup_accession   같은 accession 이 FASTA 에 2회 이상 등장. **빈 accession 은 제외**
                     (빈값끼리 묶어 세면 중복 수가 부풀려진다 — 빈 행은 n_empty_accession 으로
                      따로 보고한다).
  accession_count    그 accession 의 FASTA 등장 횟수. 빈 accession 행은 1 로 둔다.
  source_organism    헤더의 'OS=...' 또는 말미 '[...]' 에서만 추출. 못 뽑으면 빈값.
                     (UniProt 엔트리 니모닉 'CODH2_CARHZ' 의 _CARHZ 는 쓰지 않는다 —
                      전체 학명과 섞이면 조용히 틀린 집계가 되기 때문)
  source_taxid       헤더의 'OX=' (없으면 빈값)
  source_db          'sp|'→UniProtKB/Swiss-Prot, 'tr|'→UniProtKB/TrEMBL, 'ref|'→RefSeq,
                     'gi|'→NCBI, 그 외 'status=' 값, 없으면 빈값
  is_reference       accession 이 Q9F8A8 인 행
  cluster            master 의 클러스터 대표 accession을 common.norm_acc 로 정규화한 값.
                     조인 실패 행은 빈값.
  cluster_raw        master 에 적혀 있던 원본 문자열
  is_mm_rep          master 에 컬럼이 있으면 그 값, 없으면 (cluster == accession) 으로
                     유도(인계문서 1절 규칙). 어느 쪽이었는지는 params.json 참조.
  joined_master      master 조인 성공 여부
  qc_pass            종합 판정 =
                       (not qc_has_ambiguous) and (not qc_nonstandard)
                       and (qc_len_band == 'core') and is_dup_rep
                     즉 「표준 20종 문자만 · 길이 코어 구간 · 유일 서열 대표」.
                     is_putative_fusion 은 core 밴드가 이미 배제하므로 식에 넣지 않았다.
                     ※ qc_pass 는 **보존도 계산용 기본 집합**을 뜻할 뿐이고, 행 삭제의
                       근거가 아니다. 기준을 바꾸려면 임계값 옵션만 바꿔 재실행하면 된다.
  seq                서열 원문 (--no-sequence 로 생략 가능)

실행 환경
  대상은 epel-af2 워크스테이션(conda base). 의존성은 pandas + numpy + 표준 라이브러리뿐.
  실제 데이터가 없는 곳에서도 점검할 수 있도록 --dry-run 과 --self-test 를 제공한다.

사용 예
  python pipeline/stages/01_sequence_set.py                      # 기본 경로로 실행
  python pipeline/stages/01_sequence_set.py --dry-run            # 경로·존재만 확인
  python pipeline/stages/01_sequence_set.py --self-test /tmp/st   # 합성 데이터 생성+검증
  python pipeline/stages/01_sequence_set.py --fasta X.fasta --master M.tsv --outdir OUT \
         --expect-n 0 --min-join-rate 0.99
"""
from __future__ import annotations

import argparse
import hashlib
import random
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "pipeline" / "lib"))
import common  # noqa: E402

STAGE = "01_sequence_set"

# --------------------------------------------------------------------------
# 상수 — 아미노산 문자 집합과 기본 길이 밴드
# --------------------------------------------------------------------------
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_AA = set("XBZJUO")          # X 모호, B=D/N, Z=E/Q, J=I/L, U=Sec, O=Pyl
DEFAULT_CORE_MIN = 620                # 이상이면 core (Kim 2022 Nat Catal 의 단량체 길이대)
DEFAULT_CORE_MAX = 700                # 초과면 long + 융합 의심

# master 의 is_mm_rep 컬럼 이름 후보. common.COLSPEC 에는 없어서 여기 둔다.
IS_MM_REP_CANDIDATES = ["is_mm_rep", "mm_rep", "is_rep", "is_representative", "rep_flag", "is_mmseqs_rep"]
TRUE_TOKENS = {"1", "1.0", "true", "t", "yes", "y", "rep"}

DUP_GROUPS_IN_PARAMS = 200     # params.json 에 적는 중복 그룹 상세의 상한 (잘리면 명시)
DUP_GROUPS_IN_SUMMARY = 20     # 요약 md 의 중복 그룹 표 상한 (잘리면 개수를 함께 적는다)

QC_PASS_FORMULA = ("(not qc_has_ambiguous) and (not qc_nonstandard) "
                   "and (qc_len_band == 'core') and is_dup_rep")

OUT_COLUMNS = [
    "accession", "raw_header", "seq_len", "n_X", "n_ambiguous", "n_cys", "n_lowercase",
    "nonstandard_chars", "seq_sha256",
    "qc_len_band", "qc_has_ambiguous", "qc_nonstandard", "is_putative_fusion",
    "dup_group", "dup_group_size", "is_dup_rep", "is_dup_accession", "accession_count",
    "source_organism", "source_taxid", "source_db", "is_reference",
    "cluster", "cluster_raw", "is_mm_rep", "joined_master", "qc_pass", "seq",
]

# 헤더 메타데이터 추출 정규식
_OS_RE = re.compile(r"\bOS=(.+?)(?=\s+[A-Z][A-Za-z]{1,3}=|$)")
_OX_RE = re.compile(r"\bOX=(\d+)")
_BRACKET_RE = re.compile(r"\[([^\[\]]+)\]\s*$")
_STATUS_RE = re.compile(r"\bstatus=([A-Za-z]+)")
_DB_BY_PREFIX = {"sp": "UniProtKB/Swiss-Prot", "tr": "UniProtKB/TrEMBL",
                 "ref": "RefSeq", "gi": "NCBI", "up": "UniParc"}

# TSV 를 깨뜨리는 문자. 탭이 든 셀은 pandas 가 따옴표로 감싸 주지만
# awk / cut / sort 로 읽으면 컬럼이 한 칸 밀린다 (인계문서 7절이 cut -f 를 쓴다).
_TSV_UNSAFE = re.compile(r"[\t\r\n\v\f]")


# --------------------------------------------------------------------------
# 입력 인자 검증 — 잘못된 임계값은 '조용히 틀린 숫자'로 끝나므로 예외로 막는다
# --------------------------------------------------------------------------
def validate_args(args) -> None:
    if args.len_core_min > args.len_core_max:
        raise ValueError(
            f"--len-core-min({args.len_core_min}) 이 --len-core-max({args.len_core_max}) 보다 크다. "
            f"이대로 돌리면 core 밴드가 한 행도 없이 사라져 qc_pass 가 전부 False 가 되고 "
            f"is_putative_fusion 이 {args.len_core_max} aa 초과에서 켜진다 — 두 값을 바꿔 쓴 것이 "
            f"아닌지 확인할 것")
    if args.len_core_min < 1:
        raise ValueError(f"--len-core-min 은 1 이상이어야 한다 (받은 값 {args.len_core_min})")
    if args.unmatched_examples < 0:
        raise ValueError(
            f"--unmatched-examples 는 0 이상이어야 한다 (받은 값 {args.unmatched_examples}). "
            f"음수는 파이썬 슬라이스가 뒤에서 잘라내므로 예시 목록이 조용히 엉뚱해진다")
    if not (0.0 <= args.min_join_rate <= 1.0):
        raise ValueError(f"--min-join-rate 는 0.0~1.0 이어야 한다 (받은 값 {args.min_join_rate})")
    if args.expect_n < 0:
        raise ValueError(f"--expect-n 은 0(검사 끄기) 이상이어야 한다 (받은 값 {args.expect_n})")


# --------------------------------------------------------------------------
# FASTA 머리 검사 — common.read_fasta 는 첫 '>' 앞의 줄을 조용히 버린다.
# BOM 한 바이트 때문에 첫 레코드가 통째로 사라지는 경로를 예외로 바꾼다.
# --------------------------------------------------------------------------
def assert_fasta_header_first(path: Path) -> None:
    with open(path, "rb") as fh:
        head = fh.read(1 << 16)
    if not head.strip():
        raise ValueError(f"FASTA 가 비어 있거나 공백뿐이다: {path}")
    if head.startswith(b"\xef\xbb\xbf"):
        raise ValueError(
            f"FASTA 가 UTF-8 BOM 으로 시작한다: {path}\n"
            "  common.read_fasta 는 BOM 이 붙은 첫 줄을 '>' 헤더로 보지 않으므로 "
            "**첫 레코드가 조용히 사라진다** (1827 -> 1826).\n"
            "  제거: sed -i '1s/^\\xEF\\xBB\\xBF//' <fasta>")
    first = next(line for line in head.split(b"\n") if line.strip())
    if not first.startswith(b">"):
        raise ValueError(
            f"FASTA 첫 줄이 '>' 로 시작하지 않는다: {path}\n"
            f"  첫 줄={first[:80]!r}\n"
            "  common.read_fasta 는 첫 '>' 앞의 모든 줄을 조용히 버리므로 레코드가 "
            "통째로 사라진다. 머리쓰레기를 지우고 다시 돌릴 것")


# --------------------------------------------------------------------------
# accession 정규화 — common.parse_accession 은 빈 문자열에서 IndexError 를 낸다.
# 결측 셀이 섞인 입력에서 그것이 의미불명 예외로 터지는 것을 막는다(정규화 규칙은 동일).
# --------------------------------------------------------------------------
def safe_parse_accession(token: str) -> str:
    token = "" if token is None else str(token).strip()
    return common.parse_accession(token) if token else ""


def safe_norm_acc(series: pd.Series) -> pd.Series:
    s = series.astype("object").fillna("").astype(str).str.strip()
    return s.map(safe_parse_accession)


# --------------------------------------------------------------------------
# 헤더 파싱
# --------------------------------------------------------------------------
def parse_header_meta(header: str) -> dict:
    """헤더에서 생물종 / taxid / 출처 DB 를 뽑는다. 못 뽑으면 빈 문자열."""
    organism = ""
    m = _OS_RE.search(header)
    if m:
        organism = m.group(1).strip()
    else:
        m = _BRACKET_RE.search(header)
        if m:
            organism = m.group(1).strip()

    taxid = ""
    m = _OX_RE.search(header)
    if m:
        taxid = m.group(1)

    db = ""
    tok = header.split()[0] if header.split() else ""
    parts = tok.split("|")
    if len(parts) >= 3 and parts[0].lower() in _DB_BY_PREFIX:
        db = _DB_BY_PREFIX[parts[0].lower()]
    else:
        m = _STATUS_RE.search(header)
        if m:
            db = m.group(1).strip().lower()
    return dict(source_organism=organism, source_taxid=taxid, source_db=db)


# --------------------------------------------------------------------------
# FASTA -> 레코드 표
# --------------------------------------------------------------------------
def build_records(fasta_path: Path, *, core_min: int, core_max: int) -> tuple[pd.DataFrame, dict]:
    """FASTA 를 한 번 훑어 행 하나씩 만든다. 삭제·병합 없음."""
    assert_fasta_header_first(fasta_path)
    rows = []
    n_tab_sanitized = 0
    n_cells_sanitized = 0
    n_rows_cells_sanitized = 0
    for header, raw_seq in common.read_fasta(fasta_path):
        seq = re.sub(r"\s+", "", raw_seq)          # 공백/개행은 서열이 아니다
        up = seq.upper()
        chars = set(up)
        nonstd = sorted(chars - STANDARD_AA)
        length = len(seq)

        safe_header = _TSV_UNSAFE.sub(" ", header)
        if safe_header != header:
            n_tab_sanitized += 1

        if length < core_min:
            band = "short"
        elif length <= core_max:
            band = "core"
        else:
            band = "long"

        rec = dict(
            accession=safe_parse_accession(header),
            raw_header=safe_header,
            seq_len=length,
            n_X=up.count("X"),
            n_ambiguous=sum(up.count(c) for c in AMBIGUOUS_AA),
            n_cys=up.count("C"),
            n_lowercase=sum(1 for ch in seq if ch.islower()),
            nonstandard_chars="".join(nonstd),
            seq_sha256=hashlib.sha256(up.encode()).hexdigest(),
            qc_len_band=band,
            qc_has_ambiguous=bool(chars & AMBIGUOUS_AA),
            qc_nonstandard=bool(nonstd),
            is_putative_fusion=bool(length > core_max),
            seq=seq,
        )
        rec.update(parse_header_meta(header))
        # parse_header_meta 는 '원본' 헤더에서 뽑으므로 source_organism 에 탭이 남을 수 있다.
        # raw_header 만 막아 두면 TSV 가 cut/awk 에서 한 칸 밀린다 — 모든 문자열 셀을 같은
        # 규칙으로 정리하고 몇 셀을 고쳤는지 센다 (절대 규칙 7).
        n_bad = 0
        for k, v in rec.items():
            if isinstance(v, str) and _TSV_UNSAFE.search(v):
                rec[k] = _TSV_UNSAFE.sub(" ", v).strip()
                n_bad += 1
        if n_bad:
            n_cells_sanitized += n_bad
            n_rows_cells_sanitized += 1
        rows.append(rec)

    if not rows:
        raise ValueError(f"FASTA 에서 레코드를 하나도 읽지 못했다: {fasta_path}")

    df = pd.DataFrame(rows).reset_index(drop=True)

    # accession 중복 — 조인 자체는 보존되지만 하위 단계가 조용히 틀릴 수 있어 플래그로 남긴다
    acc_counts = df["accession"].value_counts()
    df["accession_count"] = df["accession"].map(acc_counts).astype(int)
    # 빈 accession 은 '같은 accession' 이 아니라 '없는 accession' 이다. 빈값끼리 묶어 세면
    # 중복 accession 수가 조용히 부풀려진다 — 빈 행은 n_empty_accession 으로만 보고한다.
    empty_acc = df["accession"] == ""
    df.loc[empty_acc, "accession_count"] = 1
    df["is_dup_accession"] = (df["accession_count"] > 1) & ~empty_acc

    n_empty_acc = int(empty_acc.sum())
    meta = dict(
        n_records=int(len(df)),
        n_empty_accession=n_empty_acc,
        # 빈값을 뺀 실제 accession 종 수. nunique() 는 "" 도 한 종으로 센다.
        n_unique_accession=int(df.loc[~empty_acc, "accession"].nunique()),
        n_unique_accession_incl_empty=int(df["accession"].nunique()),
        n_unique_sequence=int(df["seq_sha256"].nunique()),
        n_header_tab_sanitized=int(n_tab_sanitized),
        n_tsv_cells_sanitized=int(n_cells_sanitized),
        n_rows_tsv_cells_sanitized=int(n_rows_cells_sanitized),
        n_empty_sequence=int((df["seq_len"] == 0).sum()),
        n_rows_with_lowercase=int((df["n_lowercase"] > 0).sum()),
    )
    return df, meta


# --------------------------------------------------------------------------
# 참조 ChCODH2 대조
# --------------------------------------------------------------------------
def check_reference(df: pd.DataFrame, *, require: bool = False) -> dict:
    """Q9F8A8 행을 표시하고 서열을 전문 대조한다. 불일치면 예외."""
    ref_seq = common.chcodh2_seq()
    common.verify_known_positions(ref_seq)          # 번호 체계 점검 (매 실행)

    target = common.CHCODH2_ACC.upper()
    hits = [i for i, a in enumerate(df["accession"].astype("object")) if str(a).upper() == target]
    df["is_reference"] = False

    info = dict(
        reference_accession=common.CHCODH2_ACC,
        reference_len=len(ref_seq),
        known_positions_verified=True,
        found_in_fasta=bool(hits),
        n_reference_rows=len(hits),
        reference_row_index=[int(i) for i in hits],
        sequence_identical=None,
        case_difference_only=None,
    )
    if not hits:
        msg = (f"참조 {common.CHCODH2_ACC} 가 FASTA 에 없다. "
               f"InterPro 세트에 참조가 포함되지 않았을 수 있다 — 하위 단계는 "
               f"lib/chcodh2.fasta 를 질의 서열로 써야 한다.")
        if require:
            raise ValueError(msg + " (--require-reference 지정됨)")
        print(f"[주의] {msg}")
        return info

    identical, case_only = True, False
    for i in hits:
        df.loc[i, "is_reference"] = True
        got = str(df.loc[i, "seq"])
        if got.upper() != ref_seq:
            diff_at = next((k + 1 for k in range(min(len(got), len(ref_seq)))
                            if got.upper()[k] != ref_seq[k]), min(len(got), len(ref_seq)) + 1)
            raise ValueError(
                "참조 서열 불일치 — 잔기 번호 체계가 어긋났으므로 파이프라인 전체가 무효다.\n"
                f"  행 {i} accession={df.loc[i, 'accession']} header={df.loc[i, 'raw_header'][:120]}\n"
                f"  FASTA 길이={len(got)}  lib/chcodh2.fasta 길이={len(ref_seq)}\n"
                f"  첫 불일치 위치(ChCODH2 번호)={diff_at} "
                f"FASTA={got.upper()[diff_at-1:diff_at+9]!r} "
                f"REF={ref_seq[diff_at-1:diff_at+9]!r}"
            )
        if got != ref_seq:
            case_only = True
    info.update(sequence_identical=identical, case_difference_only=case_only)
    common.verify_known_positions(str(df.loc[hits[0], "seq"]).upper())
    print(f"[확인] 참조 {common.CHCODH2_ACC} FASTA 행 {hits} — 서열 {len(ref_seq)} aa 전문 일치, "
          f"논문 위치 {len(common.KNOWN_POSITIONS)}개 대조 통과")
    return info


# --------------------------------------------------------------------------
# 완전 중복 서열 그룹
# --------------------------------------------------------------------------
def assign_dup_groups(df: pd.DataFrame) -> dict:
    """대문자 변환 서열이 완전히 같은 행끼리 dup_group 을 부여한다.

    대표는 is_reference 우선 → accession 사전순 → 등장 순서. (재현 가능하게 결정적)
    중복이 없는 행도 is_dup_rep = True 로 둔다 (유일 서열 집합을 바로 뽑기 위함).
    """
    by_hash: dict[str, list[int]] = {}
    for pos, h in enumerate(df["seq_sha256"].astype("object")):
        by_hash.setdefault(str(h), []).append(pos)

    dup_group = [""] * len(df)
    dup_size = [1] * len(df)
    is_rep = [True] * len(df)

    is_ref = list(df["is_reference"]) if "is_reference" in df.columns else [False] * len(df)
    accs = [str(a) for a in df["accession"]]

    gid = 0
    group_detail = []
    for h, idxs in by_hash.items():
        if len(idxs) == 1:
            continue
        gid += 1
        name = f"dup_{gid:04d}"
        rep = sorted(idxs, key=lambda i: (0 if is_ref[i] else 1, accs[i], i))[0]
        for i in idxs:
            dup_group[i] = name
            dup_size[i] = len(idxs)
            is_rep[i] = (i == rep)
        group_detail.append(dict(dup_group=name, size=len(idxs), rep_accession=accs[rep],
                                 members=[accs[i] for i in idxs]))

    df["dup_group"] = dup_group
    df["dup_group_size"] = dup_size
    df["is_dup_rep"] = is_rep

    n_rows_in_groups = int(sum(1 for g in dup_group if g))
    n_pairs = int(sum(d["size"] * (d["size"] - 1) // 2 for d in group_detail))
    stats = dict(
        n_dup_groups=gid,
        n_rows_in_dup_groups=n_rows_in_groups,
        n_redundant_rows=n_rows_in_groups - gid,              # 대표를 뺀 잉여 행 수
        n_identical_pairs=n_pairs,                            # 완전중복 '쌍' 의 수
        largest_dup_group_size=max((d["size"] for d in group_detail), default=0),
        n_unique_sequence=int(df["seq_sha256"].nunique()),
        dup_groups=sorted(group_detail, key=lambda d: (-d["size"], d["dup_group"])),
    )
    return stats


# --------------------------------------------------------------------------
# master_1827.tsv 조인
# --------------------------------------------------------------------------
def join_master(df: pd.DataFrame, master_path: Path, *, n_examples: int = 10,
                min_join_rate: float = 0.0, allow_master_dup: bool = False) -> dict:
    """accession 으로 cluster / is_mm_rep 을 가져온다. 행 수를 반드시 보존한다."""
    n_before = len(df)
    info = dict(master_path=str(master_path), master_exists=bool(Path(master_path).exists()))

    if not Path(master_path).exists():
        # --min-join-rate 를 요구했는데 조인을 아예 하지 않으면 그 게이트는 조용히 통과된다.
        # 사용자가 명시적으로 건 검사이므로 예외로 막는다.
        if min_join_rate > 0:
            raise FileNotFoundError(
                f"--min-join-rate {min_join_rate} 를 요구했는데 master 파일이 없다: {master_path}\n"
                "  조인을 하지 않으면 조인율을 확인할 수 없고 cluster/is_mm_rep 이 전부 빈값이 된다.\n"
                "  --master 로 올바른 경로를 주거나, --min-join-rate 0 으로 검사를 명시적으로 끌 것")
        df["cluster"] = ""
        df["cluster_raw"] = ""
        df["is_mm_rep"] = False
        df["joined_master"] = False
        print(f"[주의] master 파일이 없어 cluster/is_mm_rep 을 비워 둔다: {master_path}")
        print(f"       조인율 해당없음 · 결측률 100.00% ({n_before}/{n_before} 행이 cluster 없음)")
        info.update(join_rate=None, missing_rate=1.0, n_matched=0, n_unmatched=n_before,
                    n_clusters=0, n_is_mm_rep=0, skipped=True)
        return info

    m = common.read_tsv(master_path)
    info["master_columns"] = list(m.columns)
    info["master_rows"] = int(len(m))

    cols = common.resolve_columns(m, ["accession", "cluster"], required=("accession",))
    acc_col, cl_col = cols["accession"], cols["cluster"]
    rep_col = common.detect_column(m, IS_MM_REP_CANDIDATES, required=False, what="is_mm_rep")
    info.update(master_accession_column=acc_col, master_cluster_column=cl_col,
                master_is_mm_rep_column=rep_col)

    if cl_col is None:
        raise KeyError(f"master 에서 cluster 컬럼을 못 찾음. 후보={common.COLSPEC['cluster']} "
                       f"실제={list(m.columns)}")

    mm = pd.DataFrame({
        "accession": safe_norm_acc(m[acc_col]),
        "cluster_raw": m[cl_col].astype("object").fillna("").astype(str).str.strip(),
    })
    mm["cluster"] = safe_norm_acc(mm["cluster_raw"])
    n_empty_acc = int((mm["accession"] == "").sum())
    n_empty_cluster = int((mm["cluster"] == "").sum())
    info["n_master_empty_accession"] = n_empty_acc
    info["n_master_empty_cluster"] = n_empty_cluster
    if n_empty_acc:
        raise ValueError(
            f"master 의 '{acc_col}' 컬럼이 빈 행이 {n_empty_acc}개다. 조인 키가 없으므로 "
            f"조인율이 조용히 떨어진다 — master 를 먼저 고칠 것")
    if n_empty_cluster:
        print(f"[주의] master 의 '{cl_col}' 이 빈 행 {n_empty_cluster}개 — 해당 행은 "
              f"joined_master=True 이지만 cluster 가 빈값이 된다")
    info["n_cluster_label_normalized"] = int((mm["cluster"] != mm["cluster_raw"]).sum())

    if rep_col is not None:
        raw = m[rep_col].astype("object").fillna("").astype(str).str.strip().str.lower()
        mm["is_mm_rep"] = raw.isin(TRUE_TOKENS)
        info["is_mm_rep_source"] = f"master 컬럼 '{rep_col}'"
        info["is_mm_rep_raw_value_counts"] = {str(k): int(v) for k, v in
                                              raw.value_counts().head(10).items()}
        # 토큰 어휘가 어긋나면 is_mm_rep 이 전부 False 가 되고도 아무 일 없이 끝난다.
        # mmseqs 분할에는 반드시 대표가 있으므로(cluster == accession 인 행) 그 모순을 예외로 막는다.
        n_parsed_true = int(mm["is_mm_rep"].sum())
        n_derived_true = int((mm["cluster"] == mm["accession"]).sum())
        info["n_is_mm_rep_parsed_true"] = n_parsed_true
        info["n_is_mm_rep_derived_true"] = n_derived_true
        if n_parsed_true == 0 and n_derived_true > 0:
            raise ValueError(
                f"master 의 '{rep_col}' 컬럼에서 참으로 해석된 행이 0개인데 "
                f"(cluster == accession) 인 행은 {n_derived_true}개다. 참 토큰 어휘가 어긋났다는 "
                f"뜻이며, 이대로면 is_mm_rep 이 조용히 전부 False 가 되어 하위 단계의 대표 선정이 "
                f"무효가 된다.\n"
                f"  실제 값 분포(상위 10) = {info['is_mm_rep_raw_value_counts']}\n"
                f"  인식하는 참 토큰 = {sorted(TRUE_TOKENS)}\n"
                f"  master 를 고치거나 스크립트의 TRUE_TOKENS 를 넓힐 것")
    else:
        mm["is_mm_rep"] = mm["cluster"] == mm["accession"]
        info["is_mm_rep_source"] = "유도: cluster == accession (인계문서 1절)"
    # 주어진 값과 유도 규칙이 어긋나는 행 수 — 조용한 불일치를 드러낸다
    info["n_is_mm_rep_vs_derived_mismatch"] = int(
        (mm["is_mm_rep"] != (mm["cluster"] == mm["accession"])).sum())

    dup_acc = mm["accession"][mm["accession"].duplicated(keep=False)]
    if len(dup_acc) and not allow_master_dup:
        raise ValueError(
            f"master accession 이 중복이다 ({mm['accession'].duplicated().sum()}개 중복 행). "
            f"조인에서 행이 불어나 1,827 이 깨진다. 예시={sorted(set(dup_acc))[:10]} "
            f"(--allow-master-dup-accession 으로 첫 행만 쓰도록 강제할 수 있다)")
    if len(dup_acc):
        info["master_dup_accession_dropped"] = int(mm["accession"].duplicated().sum())
        mm = mm.drop_duplicates(subset="accession", keep="first")
        print(f"[주의] master 중복 accession {info['master_dup_accession_dropped']}행을 "
              f"첫 행만 남기고 버렸다 (--allow-master-dup-accession)")

    out = df.merge(mm[["accession", "cluster", "cluster_raw", "is_mm_rep"]],
                   on="accession", how="left", validate="many_to_one")
    if len(out) != n_before:
        raise AssertionError(f"조인이 행 수를 바꿨다: {n_before} -> {len(out)}")

    # 컬럼별로 dtype 을 지켜 대입한다. 블록 대입은 object 로 뭉개져 ~(bool) 이 정수가 된다.
    joined = out["cluster"].notna().to_numpy(dtype=bool)
    df["cluster"] = out["cluster"].astype("object").fillna("").astype(str).to_numpy()
    df["cluster_raw"] = out["cluster_raw"].astype("object").fillna("").astype(str).to_numpy()
    df["is_mm_rep"] = out["is_mm_rep"].astype("object").fillna(False).to_numpy(dtype=bool)
    df["joined_master"] = joined

    n_matched = int(df["joined_master"].sum())
    rate = n_matched / n_before if n_before else float("nan")
    unmatched = [str(a) for a in df.loc[~df["joined_master"], "accession"]]
    fasta_accs = set(str(a) for a in df["accession"])
    orphans = [a for a in mm["accession"].astype("object") if str(a) not in fasta_accs]

    info.update(
        n_matched=n_matched,
        n_unmatched=n_before - n_matched,
        join_rate=round(float(rate), 6),
        missing_rate=round(1.0 - float(rate), 6),
        n_master_orphans=len(orphans),
        unmatched_examples=unmatched[:n_examples],
        unmatched_examples_truncated=bool(len(unmatched) > n_examples),
        master_orphan_examples=[str(a) for a in orphans[:n_examples]],
        master_orphan_examples_truncated=bool(len(orphans) > n_examples),
        # 빈 cluster 라벨("")은 클러스터가 아니다 — 세면 개수가 1 부풀려진다
        n_clusters=int(df.loc[df["joined_master"] & (df["cluster"] != ""), "cluster"].nunique()),
        n_is_mm_rep=int(df["is_mm_rep"].sum()),
        skipped=False,
    )

    print(f"[조인] master {len(mm)}행 ↔ FASTA {n_before}행 : 성공 {n_matched} "
          f"({rate*100:.2f}%) · 실패 {n_before - n_matched} · master 쪽 미사용 {len(orphans)}")
    print(f"       클러스터 {info['n_clusters']}개 · is_mm_rep {info['n_is_mm_rep']}개 "
          f"({info['is_mm_rep_source']})")
    if rate < 1.0:
        shown = unmatched[:n_examples]
        print(f"[주의] 조인율이 100% 가 아니다. 미매칭 accession {len(unmatched)}개 중 "
              f"{len(shown)}개만 표시:")
        for a in shown:
            print(f"         {a if a else '<빈 accession — 헤더에서 못 뽑음>'}")
        if orphans:
            o = orphans[:n_examples]
            print(f"       master 쪽에만 있는 accession {len(orphans)}개 중 {len(o)}개: "
                  f"{[str(x) for x in o]}")
        print("       → 양쪽 accession 정규화 규칙(common.parse_accession)이 다른지 확인할 것")
    if n_matched == 0:
        raise ValueError(
            "master 와 한 건도 매칭되지 않았다. accession 정규화 규칙이 서로 다르다는 뜻이다. "
            f"FASTA 예시={sorted(fasta_accs)[:5]} master 예시={sorted(set(mm['accession']))[:5]}")
    if rate < min_join_rate:
        raise ValueError(f"조인율 {rate:.4f} < --min-join-rate {min_join_rate}")
    return info


# --------------------------------------------------------------------------
# 요약
# --------------------------------------------------------------------------
def md_table(headers, rows) -> str:
    out = ["| " + " | ".join(str(h) for h in headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if v is None else str(v) for v in r) + " |")
    return "\n".join(out)


def length_histogram(lengths, bin_width=50, max_bins=40) -> tuple[list, int]:
    """텍스트 히스토그램용 (구간, 개수) 목록. 구간이 너무 많으면 폭을 2배씩 늘린다."""
    lo_raw, hi_raw = int(min(lengths)), int(max(lengths))
    while (hi_raw // bin_width - lo_raw // bin_width + 1) > max_bins:
        bin_width *= 2
    lo = (lo_raw // bin_width) * bin_width
    edges = list(range(lo, hi_raw + bin_width + 1, bin_width))
    counts, _ = np.histogram(np.asarray(lengths, dtype=float), bins=edges)
    return [(edges[i], edges[i + 1] - 1, int(counts[i])) for i in range(len(counts))], bin_width


def flag_rows(df: pd.DataFrame) -> list:
    n = len(df)
    flags = ["qc_has_ambiguous", "qc_nonstandard", "is_putative_fusion",
             "is_dup_accession", "is_dup_rep", "is_reference", "joined_master",
             "is_mm_rep", "qc_pass"]
    rows = []
    for f in flags:
        if f in df.columns:
            k = int(df[f].sum())
            rows.append([f, k, f"{k / n * 100:.2f}%", n - k])
    return rows


def write_summary(path: Path, df: pd.DataFrame, stats: dict) -> None:
    n = len(df)
    L = df["seq_len"].to_numpy()
    q = lambda p: int(np.percentile(L, p))  # noqa: E731
    fasta_rec = dict(stats["fasta"])
    band_counts = df["qc_len_band"].value_counts()
    hist, bw = length_histogram(L)
    mx = max((c for *_, c in hist), default=1) or 1
    dup = stats["duplicates"]
    jn = stats["join"]
    ref = stats["reference"]

    parts = []
    parts.append(f"# 단계 1 — 서열 세트 정의 요약\n")
    parts.append(f"생성 시각(UTC) `{stats['written_utc']}` · 스테이지 `{STAGE}`\n")

    parts.append("## 1. 개요\n")
    parts.append(md_table(["항목", "값"], [
        ["입력 FASTA", f"`{stats['input_paths']['fasta']}`"],
        ["FASTA sha256", (fasta_rec.get("sha256") or "")[:16]],
        ["행 수 (= FASTA 레코드 수)", n],
        ["기대 행 수 (--expect-n)", stats["expect_n"] or "검사 안 함"],
        ["유일 accession (빈값 제외)", fasta_rec["n_unique_accession"]],
        ["**accession 을 못 뽑은 행**", fasta_rec["n_empty_accession"]],
        ["유일 서열 (대문자 기준)", fasta_rec["n_unique_sequence"]],
        ["빈 서열 행", fasta_rec["n_empty_sequence"]],
        ["소문자 포함 행", fasta_rec["n_rows_with_lowercase"]],
        ["헤더 탭/개행 치환 행", fasta_rec["n_header_tab_sanitized"]],
        ["문자열 셀 탭/개행 치환 (셀 / 행)",
         f"{fasta_rec['n_tsv_cells_sanitized']} / {fasta_rec['n_rows_tsv_cells_sanitized']}"],
        ["길이 밴드 경계", f"short &lt; {stats['core_min']} ≤ core ≤ {stats['core_max']} &lt; long"],
    ]) + "\n")
    if fasta_rec["n_empty_accession"]:
        parts.append(f"⚠ accession 을 못 뽑은 행이 {fasta_rec['n_empty_accession']}개다. "
                     "이 행들은 master 조인에서 **반드시 실패**하므로 조인율을 그만큼 끌어내린다. "
                     "헤더 형식을 확인할 것.\n")
    parts.append("※ `sequence_set.tsv` 의 빈 셀(`cluster` · `dup_group` · `source_organism` 등)은 "
                 "`common.read_tsv`(= `dtype=str`, `keep_default_na` 기본값)로 다시 읽으면 "
                 "**빈 문자열이 아니라 NaN** 이 된다. 하위 단계에서 `!= \"\"` 로 거르면 조용히 "
                 "틀린다 — `.fillna(\"\")` 를 먼저 하거나 `.notna()` 를 쓸 것.\n")

    parts.append("## 2. 길이 분포\n")
    parts.append(md_table(["통계", "값 (aa)"], [
        ["최소", int(L.min())], ["1 퍼센타일", q(1)], ["1사분위", q(25)],
        ["중앙값", q(50)], ["평균", f"{L.mean():.1f}"], ["3사분위", q(75)],
        ["99 퍼센타일", q(99)], ["최대", int(L.max())],
    ]) + "\n")
    parts.append(md_table(["qc_len_band", "정의", "행 수", "비율"], [
        ["short", f"&lt; {stats['core_min']}", int(band_counts.get("short", 0)),
         f"{band_counts.get('short', 0) / n * 100:.2f}%"],
        ["core", f"{stats['core_min']} ~ {stats['core_max']}", int(band_counts.get("core", 0)),
         f"{band_counts.get('core', 0) / n * 100:.2f}%"],
        ["long", f"&gt; {stats['core_max']}", int(band_counts.get("long", 0)),
         f"{band_counts.get('long', 0) / n * 100:.2f}%"],
    ]) + "\n")
    parts.append(f"히스토그램 (구간 폭 {bw} aa)\n\n```")
    for lo, hi, c in hist:
        parts.append(f"{lo:>5}-{hi:<5} {c:>5} {'#' * int(round(40 * c / mx))}")
    parts.append("```\n")

    parts.append("## 3. 플래그 집계\n")
    parts.append(md_table(["플래그", "True", "비율", "False"], flag_rows(df)) + "\n")
    parts.append(f"`qc_pass` 정의 = `{QC_PASS_FORMULA}`\n")
    parts.append("`is_putative_fusion` 은 **길이(&gt; %d aa) 기반 추정일 뿐**이며 "
                 "도메인 조성을 확인한 것이 아니다. CdhA/ACS 융합 여부는 단계 2 이후 "
                 "도메인 수준에서 재판정해야 한다.\n" % stats["core_max"])
    parts.append("`qc_nonstandard` 는 정의상 `qc_has_ambiguous` 를 포함한다 "
                 "(X/B/Z/J/U/O 도 표준 20종이 아니다). 실제 문자는 `nonstandard_chars` 참조.\n")
    nsc = df.loc[df["qc_nonstandard"], "nonstandard_chars"].value_counts()
    if len(nsc):
        parts.append("비표준 문자 조성\n")
        parts.append(md_table(["nonstandard_chars", "행 수"],
                              [[k, int(v)] for k, v in nsc.items()]) + "\n")

    parts.append("## 4. 완전 중복 서열\n")
    parts.append(md_table(["항목", "값"], [
        ["중복 그룹 수", dup["n_dup_groups"]],
        ["중복에 속한 행 수", dup["n_rows_in_dup_groups"]],
        ["잉여 행 수 (대표 제외)", dup["n_redundant_rows"]],
        ["**완전중복 쌍 수**", dup["n_identical_pairs"]],
        ["가장 큰 중복 그룹 크기", dup["largest_dup_group_size"]],
        ["유일 서열 수", dup["n_unique_sequence"]],
        ["is_dup_rep = True 행 수", int(df["is_dup_rep"].sum())],
        ["중복 accession 행 수", int(df["is_dup_accession"].sum())],
    ]) + "\n")
    parts.append("동일 서열이 여러 accession 으로 들어 있으면 **보존도가 부풀려진다.** "
                 "보존도 계산에는 `is_dup_rep` 가 True 인 행만 쓰라.\n")
    shown = dup["dup_groups"][:DUP_GROUPS_IN_SUMMARY]
    if shown:
        parts.append(f"중복 그룹 (전체 {dup['n_dup_groups']}개 중 {len(shown)}개 표시)\n")
        parts.append(md_table(["dup_group", "크기", "대표", "멤버"],
                              [[d["dup_group"], d["size"], d["rep_accession"],
                                ", ".join(d["members"][:8]) +
                                (f" … (+{len(d['members']) - 8})" if len(d["members"]) > 8 else "")]
                               for d in shown]) + "\n")

    parts.append("## 5. master_1827.tsv 조인\n")
    if jn.get("skipped"):
        parts.append(f"master 파일이 없어 건너뜀: `{jn['master_path']}` — "
                     "`cluster` / `is_mm_rep` 은 빈값/False 다.\n")
        parts.append(md_table(["항목", "값"], [
            ["조인 성공률", "해당없음 (조인 미수행)"],
            ["결측률", f"100.00% ({n}/{n} 행이 cluster 없음)"],
            ["클러스터 수", 0],
        ]) + "\n")
    else:
        parts.append(md_table(["항목", "값"], [
            ["master 경로", f"`{jn['master_path']}`"],
            ["master 행 수", jn["master_rows"]],
            ["master 컬럼", ", ".join(f"`{c}`" for c in jn["master_columns"])],
            ["사용한 accession 컬럼", f"`{jn['master_accession_column']}`"],
            ["사용한 cluster 컬럼", f"`{jn['master_cluster_column']}`"],
            ["is_mm_rep 출처", jn["is_mm_rep_source"]],
            ["**조인 성공률**", f"{jn['join_rate'] * 100:.2f}% ({jn['n_matched']}/{n})"],
            ["조인 실패 행", jn["n_unmatched"]],
            ["결측률", f"{jn['missing_rate'] * 100:.2f}%"],
            ["master 쪽에만 있는 accession", jn["n_master_orphans"]],
            ["cluster 라벨 정규화로 바뀐 행", jn["n_cluster_label_normalized"]],
            ["is_mm_rep 과 (cluster==accession) 불일치", jn["n_is_mm_rep_vs_derived_mismatch"]],
            ["클러스터 수", jn["n_clusters"]],
            ["is_mm_rep True", jn["n_is_mm_rep"]],
        ]) + "\n")
        if jn["n_unmatched"]:
            ex = jn["unmatched_examples"]
            parts.append(f"미매칭 accession (전체 {jn['n_unmatched']}개 중 {len(ex)}개 표시)\n\n```")
            parts.extend(ex)
            parts.append("```\n")
            parts.append("조인율이 100% 가 아니면 **양쪽 accession 정규화 규칙이 다른 것**을 "
                         "먼저 의심하라 (`common.parse_accession` 한 곳에서만 정규화한다).\n")

    parts.append("## 6. 참조 ChCODH2 확인\n")
    parts.append(md_table(["항목", "값"], [
        ["참조 accession", ref["reference_accession"]],
        ["참조 길이", f"{ref['reference_len']} aa"],
        ["FASTA 에 존재", "예" if ref["found_in_fasta"] else "**아니오**"],
        ["해당 행 수", ref["n_reference_rows"]],
        ["행 인덱스", ref["reference_row_index"]],
        ["서열 전문 일치", "예" if ref["sequence_identical"] else ("-" if ref["sequence_identical"] is None else "아니오")],
        ["대소문자만 다름", ref["case_difference_only"]],
        ["논문 위치 16개 대조", "통과" if ref["known_positions_verified"] else "실패"],
    ]) + "\n")
    parts.append("모든 잔기 위치는 ChCODH2 서열 번호(Met1 = 1, 636 aa)다. PDB 번호를 쓰지 않는다.\n")

    parts.append("## 7. 산출물\n")
    parts.append(md_table(["파일", "설명"], [
        [f"`{Path(stats['output_paths']['tsv']).name}`", f"{n}행 × {len(stats['out_columns'])}열 기준표"],
        [f"`{Path(stats['output_paths']['summary']).name}`", "이 문서"],
        ["`params.json`", "입력 해시 · 행 수 · 모든 집계값 · 임계값"],
    ]) + "\n")

    path.write_text("\n".join(parts), encoding="utf-8")


# --------------------------------------------------------------------------
# 본체
# --------------------------------------------------------------------------
def run(args) -> dict:
    validate_args(args)
    fasta = Path(args.fasta)
    master = Path(args.master)
    outdir = Path(args.outdir)

    if args.dry_run:
        print("=== --dry-run : 경로와 존재 여부만 확인하고 아무것도 쓰지 않는다 ===")
        for label, p in [("FASTA", fasta), ("master", master),
                         ("참조 FASTA", common.REPO_LIB / "chcodh2.fasta")]:
            print(f"  {label:10s} {'있음' if p.exists() else '없음'}  {p}")
        print(f"  출력 예정 {outdir / 'sequence_set.tsv'}")
        print(f"            {outdir / 'sequence_set_summary.md'}")
        print(f"            {outdir / 'params.json'}")
        print(f"  길이 밴드 short < {args.len_core_min} ≤ core ≤ {args.len_core_max} < long")
        print(f"  qc_pass = {QC_PASS_FORMULA}")
        print(f"  기대 행 수 {args.expect_n or '검사 안 함'} · 최소 조인율 {args.min_join_rate}")
        common.verify_known_positions()
        print(f"  참조 서열 {common.CHCODH2_LEN} aa 및 논문 위치 "
              f"{len(common.KNOWN_POSITIONS)}개 번호 체계 점검 통과")
        return dict(dry_run=True)

    if not fasta.exists():
        raise FileNotFoundError(f"입력 FASTA 가 없다: {fasta} (--fasta 로 지정하거나 --dry-run)")

    print(f"=== 단계 1 서열 세트 정의 ===")
    print(f"입력 FASTA   {fasta}")
    df, fmeta = build_records(fasta, core_min=args.len_core_min, core_max=args.len_core_max)
    print(f"레코드 {fmeta['n_records']}행 · 유일 accession {fmeta['n_unique_accession']} · "
          f"유일 서열 {fmeta['n_unique_sequence']}")
    if args.expect_n and fmeta["n_records"] != args.expect_n:
        raise ValueError(
            f"FASTA 레코드 수가 {fmeta['n_records']} — 기대값 {args.expect_n} 과 다르다. "
            f"입력이 바뀌었으면 --expect-n 을 명시적으로 고치고, 0 을 주면 검사를 끈다.")
    n_dup_acc_rows = int(df["is_dup_accession"].sum())
    if n_dup_acc_rows:
        print(f"[주의] 중복 accession 행 {n_dup_acc_rows}개: 행 {fmeta['n_records']} vs 유일 "
              f"{fmeta['n_unique_accession']}(빈값 제외) — is_dup_accession 플래그 참조")
    if fmeta["n_tsv_cells_sanitized"]:
        print(f"[주의] 탭/개행이 든 문자열 셀 {fmeta['n_tsv_cells_sanitized']}개"
              f"({fmeta['n_rows_tsv_cells_sanitized']}행)를 공백으로 치환했다 — "
              f"TSV 컬럼이 밀리는 것을 막기 위함")
    if fmeta["n_empty_sequence"]:
        print(f"[주의] 빈 서열 {fmeta['n_empty_sequence']}행 (삭제하지 않고 short 밴드로 둔다)")
    if fmeta["n_empty_accession"]:
        print(f"[주의] accession 을 뽑지 못한 행 {fmeta['n_empty_accession']}개 — "
              f"master 조인에서 반드시 실패한다 (헤더 형식을 확인할 것)")

    ref_info = check_reference(df, require=args.require_reference)
    dup_stats = assign_dup_groups(df)
    print(f"[중복] 그룹 {dup_stats['n_dup_groups']}개 · 완전중복 쌍 "
          f"{dup_stats['n_identical_pairs']}쌍 · 잉여 행 {dup_stats['n_redundant_rows']}개 · "
          f"유일 서열 {dup_stats['n_unique_sequence']}개")

    join_info = join_master(df, master, n_examples=args.unmatched_examples,
                            min_join_rate=args.min_join_rate,
                            allow_master_dup=args.allow_master_dup_accession)

    df["qc_pass"] = (~df["qc_has_ambiguous"] & ~df["qc_nonstandard"]
                     & (df["qc_len_band"] == "core") & df["is_dup_rep"])

    out_cols = [c for c in OUT_COLUMNS if c in df.columns]
    if args.no_sequence:
        out_cols = [c for c in out_cols if c != "seq"]
    missing = [c for c in OUT_COLUMNS if c not in df.columns]
    if missing:
        raise AssertionError(f"출력 컬럼 누락: {missing}")
    if len(df) != fmeta["n_records"]:
        raise AssertionError(f"행 수가 변했다: {fmeta['n_records']} -> {len(df)}")

    outdir.mkdir(parents=True, exist_ok=True)
    tsv_path = outdir / "sequence_set.tsv"
    md_path = outdir / "sequence_set_summary.md"
    df[out_cols].to_csv(tsv_path, sep="\t", index=False)

    n_pass = int(df["qc_pass"].sum())
    print(f"[판정] qc_pass {n_pass}/{len(df)} ({n_pass / len(df) * 100:.2f}%)  "
          f"= {QC_PASS_FORMULA}")

    stats = dict(
        written_utc=datetime.now(timezone.utc).isoformat(),
        input_paths=dict(fasta=str(fasta), master=str(master),
                         reference_fasta=str(common.REPO_LIB / "chcodh2.fasta")),
        output_paths=dict(tsv=str(tsv_path), summary=str(md_path)),
        out_columns=out_cols,
        expect_n=args.expect_n,
        core_min=args.len_core_min,
        core_max=args.len_core_max,
        qc_pass_formula=QC_PASS_FORMULA,
        n_rows=int(len(df)),
        fasta=dict(fmeta, sha256=common.sha256_file(fasta)),
        band_counts={str(k): int(v) for k, v in df["qc_len_band"].value_counts().items()},
        flag_counts={f: int(df[f].sum()) for f in
                     ["qc_has_ambiguous", "qc_nonstandard", "is_putative_fusion",
                      "is_dup_accession", "is_dup_rep", "is_reference", "joined_master",
                      "is_mm_rep", "qc_pass"]},
        length_stats=dict(min=int(df["seq_len"].min()), max=int(df["seq_len"].max()),
                          mean=round(float(df["seq_len"].mean()), 2),
                          median=float(df["seq_len"].median())),
        duplicates=dict(
            {k: v for k, v in dup_stats.items() if k != "dup_groups"},
            dup_groups=dup_stats["dup_groups"][:DUP_GROUPS_IN_PARAMS],
            dup_groups_listed=min(len(dup_stats["dup_groups"]), DUP_GROUPS_IN_PARAMS),
            dup_groups_truncated=bool(len(dup_stats["dup_groups"]) > DUP_GROUPS_IN_PARAMS),
        ),
        join=join_info,
        reference=ref_info,
    )
    write_summary(md_path, df, stats)

    params = common.write_params(
        outdir, STAGE,
        inputs=[fasta] + ([master] if master.exists() else []) + [common.REPO_LIB / "chcodh2.fasta"],
        outputs=[tsv_path, md_path],
        commands=[" ".join(sys.argv)],
        **{k: v for k, v in stats.items() if k not in ("written_utc",)},
    )
    print(f"출력  {tsv_path}\n      {md_path}\n      {params}")
    return dict(stats, dataframe=df)


# --------------------------------------------------------------------------
# 합성 데이터 자기검증 — 실제 데이터 없이 각 플래그가 의도대로 켜지는지 확인
# --------------------------------------------------------------------------
def make_synthetic(dirpath: Path) -> dict:
    """합성 FASTA + 합성 master tsv 를 만든다. 결정적(seed 고정)."""
    dirpath.mkdir(parents=True, exist_ok=True)
    rnd = random.Random(20261022)
    alpha = sorted(STANDARD_AA)
    gen = lambda n: "M" + "".join(rnd.choice(alpha) for _ in range(n - 1))  # noqa: E731

    ref = common.chcodh2_seq()
    core_a, core_b = gen(636), gen(690)
    short_s, long_s = gen(400), gen(760)
    amb = gen(650)
    amb = amb[:100] + "XXX" + amb[103:]                 # X 3개
    nonstd = gen(645)
    nonstd = nonstd[:50] + "U" + nonstd[51:] + "*"      # U(모호집합) + *(그 외)
    lower = gen(640)
    lower = lower[:30] + lower[30:60].lower() + lower[60:]

    recs = [
        ("tr|SYN00001|SYN1_ORGA normal core OS=Organismus alpha OX=111 GN=cooS PE=4 SV=1", core_a),
        ("tr|SYN00002|SYN2_ORGB normal core OS=Organismus beta OX=222 GN=cooS PE=4 SV=1", core_b),
        ("tr|SYN00003|SYN3_ORGC short fragment OS=Organismus gamma OX=333", short_s),
        ("tr|SYN00004|SYN4_ORGD long putative fusion OS=Organismus delta OX=444", long_s),
        ("tr|SYN00005|SYN5_ORGE ambiguous X OS=Organismus epsilon OX=555", amb),
        ("tr|SYN00006|SYN6_ORGF nonstandard U and stop OS=Organismus zeta OX=666", nonstd),
        ("tr|SYN00007|SYN7_ORGG lowercase masked OS=Organismus eta OX=777", lower),
        # 완전 중복 쌍 (core_a 와 동일)
        ("tr|SYN00008|SYN8_ORGH exact duplicate of SYN00001 OS=Organismus theta OX=888", core_a),
        # 완전 중복 3중 그룹 (core_b 와 동일)
        ("tr|SYN00009|SYN9_ORGI exact duplicate of SYN00002 OS=Organismus iota OX=999", core_b),
        ("sp|SYN00010|SY10_ORGJ exact duplicate of SYN00002 [Organismus kappa]", core_b),
        # 중복 accession (서열은 다름)
        ("tr|SYN00003|SYN3_ORGC2 same accession different sequence OS=Organismus gamma OX=333",
         gen(660)),
        # 참조
        ("tr|Q9F8A8|CODH2_CARHZ carbon monoxide dehydrogenase 2 "
         "OS=Carboxydothermus hydrogenoformans OX=246194 GN=cooS2 PE=1 SV=1", ref),
    ]
    fasta = dirpath / "synthetic_ipr.fasta"
    with open(fasta, "wt") as fh:
        for h, s in recs:
            fh.write(f">{h}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # 합성 master — 일부러 2개 accession(SYN00006, SYN00007)을 빼서 조인율 100% 미만을 만든다.
    #              cluster 라벨은 'tr|X|Y' 꼴로 둬서 정규화 경로도 태운다.
    rows = [
        ("tr|SYN00001|SYN1_ORGA", "tr|SYN00001|SYN1_ORGA", "1"),
        ("tr|SYN00002|SYN2_ORGB", "tr|SYN00002|SYN2_ORGB", "1"),
        ("tr|SYN00003|SYN3_ORGC", "tr|SYN00001|SYN1_ORGA", "0"),
        ("tr|SYN00004|SYN4_ORGD", "tr|SYN00004|SYN4_ORGD", "1"),
        ("tr|SYN00005|SYN5_ORGE", "tr|SYN00002|SYN2_ORGB", "0"),
        ("tr|SYN00008|SYN8_ORGH", "tr|SYN00001|SYN1_ORGA", "0"),
        ("tr|SYN00009|SYN9_ORGI", "tr|SYN00002|SYN2_ORGB", "0"),
        ("sp|SYN00010|SY10_ORGJ", "tr|SYN00002|SYN2_ORGB", "0"),
        ("tr|Q9F8A8|CODH2_CARHZ", "tr|Q9F8A8|CODH2_CARHZ", "1"),
        ("tr|SYN99999|ORPHAN", "tr|SYN99999|ORPHAN", "1"),   # FASTA 에 없는 고아 행
    ]
    master = dirpath / "synthetic_master.tsv"
    with open(master, "wt") as fh:
        fh.write("accession\tcluster\tis_mm_rep\n")
        for a, c, r in rows:
            fh.write(f"{a}\t{c}\t{r}\n")
    return dict(fasta=fasta, master=master, n_records=len(recs), n_master_rows=len(rows))


def self_test(dirpath: Path, base_args) -> None:
    syn = make_synthetic(dirpath)
    print(f"=== --self-test : 합성 데이터 {syn['n_records']} 서열 / master "
          f"{syn['n_master_rows']}행 ===")
    print(f"  {syn['fasta']}\n  {syn['master']}\n")

    args = argparse.Namespace(**vars(base_args))
    args.self_test = None
    args.fasta = str(syn["fasta"])
    args.master = str(syn["master"])
    args.outdir = str(dirpath / "out")
    args.expect_n = syn["n_records"]
    args.min_join_rate = 0.0
    args.dry_run = False
    res = run(args)
    def one(acc):
        sub = res["dataframe"][res["dataframe"]["accession"] == acc]
        return sub.iloc[0]

    checks = []
    checks.append(("행 수 보존", len(res["dataframe"]) == syn["n_records"]))
    checks.append(("short 밴드", one("SYN00003")["qc_len_band"] == "short"))
    checks.append(("core 밴드", one("SYN00001")["qc_len_band"] == "core"))
    checks.append(("long 밴드 + 융합 의심",
                   one("SYN00004")["qc_len_band"] == "long" and bool(one("SYN00004")["is_putative_fusion"])))
    checks.append(("X 3개 + 모호 플래그",
                   one("SYN00005")["n_X"] == 3 and bool(one("SYN00005")["qc_has_ambiguous"])))
    checks.append(("비표준 문자 플래그(U,*)",
                   bool(one("SYN00006")["qc_nonstandard"])
                   and set(one("SYN00006")["nonstandard_chars"]) == {"U", "*"}))
    checks.append(("소문자 30자 집계", one("SYN00007")["n_lowercase"] == 30))
    checks.append(("중복 그룹 2개", res["duplicates"]["n_dup_groups"] == 2))
    checks.append(("완전중복 쌍 4쌍 (2쌍짓기 1 + 3중 3)",
                   res["duplicates"]["n_identical_pairs"] == 4))
    checks.append(("잉여 행 3개", res["duplicates"]["n_redundant_rows"] == 3))
    checks.append(("중복 대표는 1그룹당 1행",
                   int(res["dataframe"]["is_dup_rep"].sum()) == res["duplicates"]["n_unique_sequence"]))
    checks.append(("중복 accession 플래그", int(res["dataframe"]["is_dup_accession"].sum()) == 2))
    checks.append(("참조 발견 + 전문 일치",
                   res["reference"]["found_in_fasta"] and res["reference"]["sequence_identical"]))
    checks.append(("참조 행 qc_pass", bool(one("Q9F8A8")["qc_pass"])))
    checks.append(("조인 매칭 10행 / 미매칭 2행",
                   res["join"]["n_matched"] == 10 and res["join"]["n_unmatched"] == 2))
    checks.append(("미매칭 예시 출력됨",
                   res["join"]["unmatched_examples"] == ["SYN00006", "SYN00007"]))
    checks.append(("중복 accession 1건이 master 1행에 many_to_one 으로 붙음",
                   int((res["dataframe"]["accession"] == "SYN00003").sum()) == 2
                   and int(res["dataframe"].loc[res["dataframe"]["accession"] == "SYN00003",
                                                "joined_master"].sum()) == 2))
    checks.append(("master 고아 1개", res["join"]["n_master_orphans"] == 1))
    checks.append(("cluster 라벨 정규화 수행", res["join"]["n_cluster_label_normalized"] == 10))
    checks.append(("cluster 정규화 결과가 accession 과 비교 가능",
                   one("SYN00003")["cluster"] == "SYN00001"))
    checks.append(("is_mm_rep master 컬럼 사용",
                   res["join"]["master_is_mm_rep_column"] == "is_mm_rep"))
    checks.append(("클러스터 4개 (SYN00001/SYN00002/SYN00004/Q9F8A8)",
                   res["join"]["n_clusters"] == 4))
    expect_pass = set()
    for _, r in res["dataframe"].iterrows():
        if (not r["qc_has_ambiguous"] and not r["qc_nonstandard"]
                and r["qc_len_band"] == "core" and r["is_dup_rep"]):
            expect_pass.add(r["raw_header"])
    got_pass = set(res["dataframe"].loc[res["dataframe"]["qc_pass"], "raw_header"])
    checks.append(("qc_pass 식 일치", expect_pass == got_pass))
    checks.append(("TSV 행 수 = 1 + n",
                   sum(1 for _ in open(Path(args.outdir) / "sequence_set.tsv")) == syn["n_records"] + 1))
    checks.append(("요약 md 존재", (Path(args.outdir) / "sequence_set_summary.md").exists()))
    checks.append(("params.json 존재", (Path(args.outdir) / "params.json").exists()))

    checks.append(("문자열 셀 탭 치환 집계 존재",
                   "n_tsv_cells_sanitized" in res["fasta"]))

    # ---------------- 음성 경로 (예외가 떠야 하는 것들) ----------------
    # 전부 '조용히 틀린 숫자'로 끝나던 경로다. 되돌아가면 여기서 잡힌다.
    neg_dir = dirpath / "neg"
    neg_dir.mkdir(parents=True, exist_ok=True)

    def must_raise(name, mutate, exc=Exception):
        a = argparse.Namespace(**vars(args))
        a.outdir = str(neg_dir / name)
        mutate(a)
        try:
            run(a)
        except exc as e:
            print(f"  (음성경로) {name} -> {type(e).__name__}: {str(e).splitlines()[0][:110]}")
            return True
        except Exception as e:  # noqa: BLE001
            print(f"  (음성경로) {name} -> 예상과 다른 예외 {type(e).__name__}: {e}")
            return False
        return False

    def _inverted(a):
        a.len_core_min, a.len_core_max = 700, 620
    checks.append(("밴드 역전은 예외", must_raise("inverted_bands", _inverted, ValueError)))

    def _neg_examples(a):
        a.unmatched_examples = -1
    checks.append(("음수 --unmatched-examples 는 예외",
                   must_raise("neg_examples", _neg_examples, ValueError)))

    def _bad_rate(a):
        a.min_join_rate = 2.0
    checks.append(("범위 밖 --min-join-rate 는 예외",
                   must_raise("bad_rate", _bad_rate, ValueError)))

    def _gate_without_master(a):
        a.master = str(neg_dir / "does_not_exist.tsv")
        a.min_join_rate = 1.0
    checks.append(("master 없는데 --min-join-rate 요구하면 예외",
                   must_raise("gate_no_master", _gate_without_master, FileNotFoundError)))

    bom = neg_dir / "bom.fasta"
    bom.write_bytes(b"\xef\xbb\xbf" + Path(syn["fasta"]).read_bytes())

    def _bom(a):
        a.fasta = str(bom)
        a.expect_n = 0
    checks.append(("BOM FASTA 는 예외 (첫 레코드 조용히 소실 방지)",
                   must_raise("bom", _bom, ValueError)))

    junk = neg_dir / "junk.fasta"
    junk.write_bytes(b"# comment line\n" + Path(syn["fasta"]).read_bytes())

    def _junk(a):
        a.fasta = str(junk)
        a.expect_n = 0
    checks.append(("'>' 아닌 머리줄은 예외", must_raise("junk_head", _junk, ValueError)))

    # master 의 is_mm_rep 토큰이 참으로 해석되지 않는 경우
    bad_master = neg_dir / "master_bad_rep.tsv"
    src = Path(syn["master"]).read_text().splitlines()
    with open(bad_master, "wt") as fh:
        fh.write(src[0] + "\n")
        for line in src[1:]:
            a_, c_, _r = line.split("\t")
            fh.write(f"{a_}\t{c_}\t{c_}\n")      # 참 토큰 대신 대표 accession 문자열

    def _bad_rep(a):
        a.master = str(bad_master)
    checks.append(("is_mm_rep 참 토큰 0개는 예외", must_raise("bad_rep_token", _bad_rep, ValueError)))

    # 탭이 든 헤더 -> source_organism 까지 치환되는지 (TSV 컬럼 밀림 방지)
    tabf = neg_dir / "tab_header.fasta"
    with open(tabf, "wt") as fh:
        fh.write(">tr|TAB001|T_ORG OS=Organismus\talpha OX=1\n" + "M" * 640 + "\n")
    a2 = argparse.Namespace(**vars(args))
    a2.fasta, a2.expect_n = str(tabf), 0
    a2.master = str(neg_dir / "does_not_exist.tsv")
    a2.min_join_rate = 0.0
    a2.outdir = str(neg_dir / "tab_out")
    r2 = run(a2)
    tab_line = (Path(a2.outdir) / "sequence_set.tsv").read_text().splitlines()[1]
    hdr_line = (Path(a2.outdir) / "sequence_set.tsv").read_text().splitlines()[0]
    checks.append(("탭 헤더: 출력 TSV 필드 수가 헤더와 같다",
                   tab_line.count("\t") == hdr_line.count("\t")))
    checks.append(("탭 헤더: source_organism 에 탭이 남지 않는다",
                   "\t" not in str(r2["dataframe"].iloc[0]["source_organism"])))
    checks.append(("탭 헤더: 치환 셀 수 보고됨", r2["fasta"]["n_tsv_cells_sanitized"] >= 1))

    # 빈 accession 은 중복 accession 으로 세지 않는다
    emptyf = neg_dir / "empty_acc.fasta"
    with open(emptyf, "wt") as fh:
        for i in range(3):
            fh.write(">\n" + "M" + "A" * (630 + i) + "\n")
    a3 = argparse.Namespace(**vars(args))
    a3.fasta, a3.expect_n = str(emptyf), 0
    a3.master = str(neg_dir / "does_not_exist.tsv")
    a3.min_join_rate = 0.0
    a3.outdir = str(neg_dir / "empty_acc_out")
    r3 = run(a3)
    checks.append(("빈 accession 3행은 is_dup_accession 0",
                   int(r3["dataframe"]["is_dup_accession"].sum()) == 0
                   and r3["fasta"]["n_empty_accession"] == 3))
    checks.append(("빈 accession 은 유일 accession 수에서 제외",
                   r3["fasta"]["n_unique_accession"] == 0))

    print("\n=== 자기검증 결과 ===")
    bad = []
    for name, ok in checks:
        print(f"  {'OK  ' if ok else 'FAIL'} {name}")
        if not ok:
            bad.append(name)
    if bad:
        raise AssertionError(f"자기검증 실패: {bad}")
    print(f"\n합성 데이터 자기검증 {len(checks)}항목 전부 통과. 산출물: {args.outdir}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="01_sequence_set.py",
        description="단계 1 — IPR010047 FASTA 에 QC 플래그를 붙여 sequence_set.tsv 를 만든다 "
                    "(행 삭제 없음).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="qc_pass = " + QC_PASS_FORMULA)
    p.add_argument("--fasta", default=str(common.IPR_FASTA), help="입력 IPR010047 FASTA")
    p.add_argument("--master", default=str(common.MASTER_TSV),
                   help="master_1827.tsv (없으면 cluster/is_mm_rep 을 비워 둔다)")
    p.add_argument("--outdir", default=str(common.OUT_ROOT / STAGE), help="산출물 디렉터리")
    p.add_argument("--expect-n", type=int, default=1827,
                   help="기대 레코드 수. 다르면 예외. 0 이면 검사 안 함 (기본 1827)")
    p.add_argument("--len-core-min", type=int, default=DEFAULT_CORE_MIN,
                   help=f"core 밴드 하한 (기본 {DEFAULT_CORE_MIN})")
    p.add_argument("--len-core-max", type=int, default=DEFAULT_CORE_MAX,
                   help=f"core 밴드 상한 (기본 {DEFAULT_CORE_MAX})")
    p.add_argument("--min-join-rate", type=float, default=0.0,
                   help="master 조인율이 이 값 미만이면 예외 (기본 0.0 = 보고만. "
                        "단 0건 매칭은 항상 예외)")
    p.add_argument("--unmatched-examples", type=int, default=10,
                   help="미매칭 accession 을 몇 개까지 보여줄지 (기본 10)")
    p.add_argument("--require-reference", action="store_true",
                   help="참조 Q9F8A8 이 FASTA 에 없으면 예외")
    p.add_argument("--allow-master-dup-accession", action="store_true",
                   help="master accession 중복 시 첫 행만 쓰고 계속한다 (기본은 예외)")
    p.add_argument("--no-sequence", action="store_true", help="출력 TSV 에서 seq 컬럼을 뺀다")
    p.add_argument("--dry-run", action="store_true", help="경로·존재만 확인하고 쓰지 않는다")
    p.add_argument("--self-test", metavar="DIR",
                   help="DIR 에 합성 데이터를 만들어 전 경로를 실행하고 플래그를 검증한다")
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    if args.self_test:
        self_test(Path(args.self_test), args)
        return 0
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
