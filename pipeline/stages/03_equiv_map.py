#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""단계 3 — 좌표계 구축. 상동체의 모든 잔기를 ChCODH2 번호(1~636)로 옮긴다.

이 단계는 파이프라인의 척추다. 여기서 1잔기만 밀리면 보존도·Fix/Bias 대상·도킹 박스
중심까지 **전부 조용히** 틀린다. 그래서 산출물보다 검증이 먼저다.

무엇을 읽는가
  1. chcodh2_alignment.tsv   --alignment  기본 common.ALIGN_TSV
       /mnt/af2results/mingyu/.../homolog_selection/chcodh2_alignment.tsv
       ChCODH2 를 질의로 1,827 서열에 돌린 mmseqs search 산출물.
       **그 search 명령은 기록이 남아있지 않다**(인계문서 2절). 그러므로 이 파일은
       「신뢰하는 입력」이 아니라 「감사 대상」이다. (A) 절이 그 감사다.
       qaln/taln 이 들어 있으므로 `convertalis --format-output` (또는 easy-search) 산출물이며,
       헤더가 없을 수 있다. 그 경우 --columns 로 --format-output 순서를 직접 준다.
  2. master_1827.tsv (선택)  --master   기본 common.MASTER_TSV
       accession 컬럼만 읽어 **universe(전체 서열 집합)** 로 쓴다. 정렬 파일에 빠진
       accession 을 세려면 분모가 필요하기 때문이다.
  3. IPR010047 FASTA (선택)  --fasta    기본 common.IPR_FASTA
       master 가 없을 때의 universe 대체 출처. 둘 다 없으면 정렬 파일의 accession 만으로
       돌아가고, 「빠진 서열을 셀 수 없다」는 사실을 불일치로 보고한다.
  4. pipeline/lib/chcodh2.fasta — common.chcodh2_seq(). 검증 1의 정답지이고
     common.verify_known_positions() 로 매 실행 번호 체계를 점검한다.

무엇을 쓰는가  (--outdir, 기본 common.OUT_ROOT/03_equiv_map)
  alignment_audit.md             (A) 입력 감사 + (B) 검증 증거 + 해석.
                                 **파싱 전에 한 번 쓰고, 파싱 성공 후 덮어쓴다.**
                                 파싱이 예외로 죽어도 감사 결과는 디스크에 남는다.
  equiv_map.tsv                  long. accession · q_pos(1~636) · q_aa · t_pos · t_aa
                                 한 행 = 잔기 대응 한 쌍 (양쪽 모두 잔기인 열만).
  residue_matrix.tsv             wide. 행=accession, 열=1..636, 값=정렬된 상동체 잔기.
                                 보존도 계산의 직접 입력.
  known_position_equivalents.tsv 행=accession, 열=KNOWN_POSITIONS 16개 ('pos41_F' 꼴).
                                 클러스터 특성화에 바로 쓰는 표.
  position_coverage.tsv          636행. 위치별 잔기 수/갭 수/미커버 수 — 보존도의 분모.
  params.json                    common.write_params. 행 수·조인율·결측률·상한·판정 근거 전부.
                                 **실패해도 status="failed" 로 남긴다.**

무엇을 보장하는가 — 깨지면 경고가 아니라 예외다
  검증 1 (오프셋)  qaln 에서 갭을 뺀 문자열 == 참조서열[qstart..] 부분서열. 전 행 정확히 일치.
                   mmseqs 의 qstart 는 1-based 라고 가정하지만, 그 가정을 **믿지 않는다.**
                   1-based / 0-based 양쪽으로 전 행을 맞춰보고 이긴 쪽을 채택하며,
                   채택 근거(각 base 의 일치 행 수)를 stdout·md·params.json 에 적는다.
                   채택한 base 로도 안 맞는 행이 하나라도 있으면 예외.
                   이것이 1-based/0-based 오프셋 오류를 잡는 그물이다.
  검증 2 (길이)    len(qaln) == len(taln). 다르면 예외. (base 판정보다 먼저 한다)
  검증 3 (자기히트) 자기 자신(Q9F8A8) 히트가 있으면 equiv_map 에서 q_pos == t_pos 가
                   636개 전부 성립하고 q_aa == t_aa 여야 한다. 아니면 예외
                   (--allow-partial-self-hit 로만 불일치 보고로 완화).
  구조 가정        단일 질의(query 컬럼이 두 값 이상이면 예외), 필수 컬럼 존재,
                   qstart 범위(start0 >= 0, start0 + 잔기수 <= 636), 행 수 보존,
                   position_coverage 의 분해 항등식.
  데이터 기대값    정렬 행 수 < --expected-n, 자기 히트 부재, universe 부재, qend/tend 불일치
                   같은 것은 예외가 아니라 `discrepancies` 로 모아 stdout·md·params.json
                   세 곳에 적고, --strict 면 산출물을 모두 쓴 뒤 종료코드 2 로 끝낸다.
                   (단계 0/2와 같은 규약. 산출물을 잃지 않기 위해 마지막에 끝낸다.)
  서열 보존        행을 지우지 않는다. residue_matrix / known_position_equivalents 의 행은
                   기본적으로 **universe 전체**(1,827)이고, 정렬에 없는 accession 은 전부
                   갭 행으로 남는다. 삭제가 아니라 갭이다. --matrix-aligned-only 로만 줄인다.
                   universe 에 없고 정렬에만 있는 accession 도 버리지 않고 뒤에 붙이며,
                   그 사실을 불일치로 보고한다.

좌표 규약 — 오독하면 전부 틀린다
  q_pos   **ChCODH2 서열 번호. Met1 = 1, 1~636.** PDB 번호가 아니다.
  t_pos   상동체 자신의 1-based 번호. 입력 base 판정 결과를 1-based 로 환산해 적는다
          (입력이 0-based 였다면 +1 해서 쓴다. 어느 쪽이었는지는 params.json 참조).
  갭      residue_matrix 의 '-' 는 두 가지다 —
            (가) 정렬 블록 안의 갭(상동체에 그 잔기가 없음, taln 의 '-')
            (나) 정렬 블록 밖(그 위치가 애초에 정렬되지 않음) 또는 정렬 파일에 없는 accession
          둘은 position_coverage.tsv 의 n_gap_within_aln / n_outside_aln /
          n_not_in_alignment 로 구분된다. 행 단위로는 known_position_equivalents.tsv 의
          in_alignment 컬럼으로 구분하라. (--outside-char 로 (나)만 다른 문자로 쓸 수 있지만
          기본값은 '-' 다. 바꾸면 행렬의 알파벳이 달라지므로 하위 단계와 합의할 것.)
  삽입    상동체에만 있는 잔기(qaln 이 '-' 인 열)는 ChCODH2 번호를 갖지 않으므로 equiv_map 과
          residue_matrix 양쪽에서 **구조적으로 표현되지 않는다.** 참조 좌표계 투영의 본질적
          한계다. 몇 개가 그렇게 사라졌는지는 params.json 의 insertion_* 와 감사 md 7절에
          숫자로 적는다.

실행 환경
  대상은 epel-af2 워크스테이션(conda base). 의존성은 pandas + numpy + 표준 라이브러리뿐.
  실제 데이터가 없는 곳에서도 점검할 수 있도록 --dry-run 과 --self-test 를 제공한다.

common 에 없어서 이 스크립트가 들고 있는 것 (common.py 는 수정하지 않았다)
  EXTRA_COLSPEC        qend/tend/qlen/tlen/alnlen 후보 이름 (common.COLSPEC 에 없다)
  HEADER_FIELD_TOKENS  헤더 행 판정용 mmseqs 필드 이름 어휘
  read_alignment()     common.read_tsv 를 쓰지 않는다. 이유는 그 함수 주석 참조(아래 구현).
  safe_parse_accession 빈 셀에서 IndexError 를 내지 않는 래퍼 (단계 1/2와 동일)
  ensure_utf8_io()     LC_ALL=C 같은 ASCII 로케일에서 **전 과정을 돈 뒤** params.json 쓰기만
                       죽는 것을 막는다 (common.write_params 는 encoding 을 지정하지 않는다).
                       미리 크게 실패시키고 PYTHONUTF8=1 처방을 알려준다. 자세한 설명은 구현 주석.

사용 예
  python pipeline/stages/03_equiv_map.py                       # 기본 경로로 실행
  python pipeline/stages/03_equiv_map.py --dry-run             # 경로·규약만 확인
  python pipeline/stages/03_equiv_map.py --self-test /tmp/st3   # 합성 데이터 생성+검증
  python pipeline/stages/03_equiv_map.py --alignment A.tsv --outdir OUT \
         --columns query,target,fident,alnlen,qstart,qend,tstart,tend,qaln,taln,evalue,bits,qcov,tcov
"""
from __future__ import annotations

import argparse
import json
import locale
import os
import random
import re
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "pipeline" / "lib"))
import common  # noqa: E402

STAGE = "03_equiv_map"

# --------------------------------------------------------------------------
# 상수
# --------------------------------------------------------------------------
EXPECTED_N_SEQUENCES = 1827          # 인계문서 1절 (IPR010047 다운로드 서열 수)
MMSEQS_DEFAULT_MAX_SEQS = 300        # search 기본값. 행 수가 이 근처면 절단 의심
MMSEQS_DEFAULT_EVALUE = 0.001        # search 기본 e-value

DEFAULT_GAP_CHAR = "-"
DEFAULT_MAX_EXAMPLES = 10            # 누락/불일치 예시 상한 (잘리면 개수를 명시)
DEFAULT_MAX_TABLE_ROWS = 40          # 감사 md 표 행 수 상한
DEFAULT_MAX_DETAIL = 200             # params.json 상세 목록 상한

# 필수/선택 컬럼 — 존재 확인은 common.resolve_columns 로 한다.
REQUIRED_KEYS = ("accession", "qstart", "qaln", "tstart", "taln")
OPTIONAL_KEYS = ("query", "fident", "qcov", "tcov", "evalue", "bits")

# common.COLSPEC 에 없는 컬럼의 후보 이름. common.py 를 고치지 않기 위해 여기 둔다.
EXTRA_COLSPEC = {
    "qend":   ["qend", "q_end"],
    "tend":   ["tend", "t_end", "send"],
    "qlen":   ["qlen", "q_len", "query_len"],
    "tlen":   ["tlen", "t_len", "target_len"],
    "alnlen": ["alnlen", "aln_len", "length", "alignment_length"],
}

# 헤더 행 판정 어휘. mmseqs --format-output 이 쓰는 이름 + COLSPEC/EXTRA_COLSPEC 의 후보.
HEADER_FIELD_TOKENS = {
    "query", "target", "evalue", "gapopen", "pident", "fident", "nident", "qstart",
    "qend", "qlen", "tstart", "tend", "tlen", "alnlen", "raw", "bits", "cigar",
    "qseq", "tseq", "qheader", "theader", "qaln", "taln", "mismatch", "qcov", "tcov",
    "qset", "qsetid", "tset", "tsetid", "taxid", "taxname", "taxlineage", "accession",
    "cluster", "bitscore", "sseqid", "qseqid", "subject", "score", "length", "identity",
}

# 등가 매핑 순서(중복 HSP 중 하나를 고르는 규칙). 과제 지정: bits 최대 -> evalue 최소 -> qcov 최대.
HSP_TIEBREAK = [("bits", False), ("evalue", True), ("qcov", False)]
HSP_TIEBREAK_DESC = "bits 최대 → evalue 최소 → qcov 최대 → 파일 등장 순서"

EQUIV_COLUMNS = ["accession", "q_pos", "q_aa", "t_pos", "t_aa"]
COVERAGE_COLUMNS = [
    "q_pos", "q_aa", "is_known_position", "known_role",
    "n_universe", "n_in_alignment", "n_not_in_alignment",
    "n_residue", "n_gap_within_aln", "n_outside_aln",
    "n_identical", "residue_frac_in_alignment", "residue_frac_universe",
    "identity_frac_of_residues",
]


# --------------------------------------------------------------------------
# 인코딩 점검 — 산출물을 잃지 않기 위해 가장 먼저 한다.
#
# epel-af2 에 LC_ALL=C / LANG 미설정으로 들어오면 파이썬의 기본 텍스트 인코딩이
# ASCII 가 된다. 그러면 (1) 이 스크립트의 한국어 print 가 첫 줄부터
# UnicodeEncodeError 로 죽고, (2) common.write_params 가 Path.write_text 로
# ensure_ascii=False JSON 을 쓰므로 **전 과정을 다 돌린 뒤 마지막에** params.json
# 쓰기가 죽는다. (2)가 더 나쁘다 — 산출 TSV 는 남았는데 재현 정보만 없는 상태가 된다.
# common.py 는 수정하지 않기로 했으므로, 여기서 미리 크게 실패시키고 처방을 알려준다.
# --------------------------------------------------------------------------
def ensure_utf8_io() -> dict:
    """stdout/stderr 를 UTF-8 로 올리고, 파일 쓰기 인코딩이 ASCII 면 예외."""
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")
        except (AttributeError, ValueError, OSError):
            pass
    enc = (locale.getpreferredencoding(False) or "").lower().replace("-", "").replace("_", "")
    info = dict(preferred_encoding=locale.getpreferredencoding(False),
                utf8_mode=bool(sys.flags.utf8_mode),
                lang=os.environ.get("LANG"), lc_all=os.environ.get("LC_ALL"))
    if enc != "utf8" and not sys.flags.utf8_mode:
        raise RuntimeError(
            f"파일 쓰기 기본 인코딩이 '{locale.getpreferredencoding(False)}' 다 (UTF-8 이 아니다). "
            f"이대로 돌리면 전 과정을 마친 뒤 params.json 쓰기에서 UnicodeEncodeError 로 죽어 "
            f"재현 정보만 사라진다. 다음 중 하나로 다시 실행할 것:\n"
            f"  PYTHONUTF8=1 python {Path(sys.argv[0]).name} ...\n"
            f"  python -X utf8 {Path(sys.argv[0]).name} ...\n"
            f"  LANG=en_US.UTF-8 python {Path(sys.argv[0]).name} ...\n"
            f"(현재 LANG={info['lang']} LC_ALL={info['lc_all']})")
    return info


# ==========================================================================
# 작은 유틸
# ==========================================================================
def safe_parse_accession(token) -> str:
    """common.parse_accession 과 같은 규칙. 빈 셀에서 IndexError 를 내지 않는다.

    (common.parse_accession 은 빈 문자열에서 header.split()[0] 으로 터진다. 단계 1/2도
     같은 래퍼를 들고 있으므로 규칙은 한 곳(common)에만 있고 방어만 여기서 한다.)
    """
    if token is None:
        return ""
    s = str(token).strip()
    if not s or s.lower() == "nan":
        return ""
    return common.parse_accession(s)


def norm_acc_series(series: pd.Series) -> pd.Series:
    return series.astype(str).map(safe_parse_accession)


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
    if isinstance(x, (float, np.floating)):
        return f"{float(x):.{nd}f}"
    return str(x)


def pct(num, den, nd=2) -> str:
    if not den:
        return "—"
    return f"{num / den * 100:.{nd}f}%"


def cap_list(values, limit):
    """목록을 상한으로 자르고 (잘린목록, 전체수, 잘림여부) 를 돌려준다. 규칙 7."""
    values = list(values)
    return values[:limit], len(values), len(values) > limit


def numeric_summary(series: pd.Series | None, name: str) -> dict:
    """수치 컬럼의 분포. 파싱 실패 개수를 반드시 같이 센다."""
    if series is None:
        return dict(column=name, present=False)
    raw = series.astype(str).str.strip()
    vals = pd.to_numeric(raw, errors="coerce")
    blank = int((raw == "").sum() + (raw.str.lower() == "nan").sum())
    bad = int(vals.isna().sum()) - blank
    ok = vals.dropna().to_numpy(dtype=float)
    rec = dict(column=name, present=True, n=int(len(raw)),
               n_blank=blank, n_unparseable=max(bad, 0), n_valid=int(ok.size))
    if ok.size:
        q = np.percentile(ok, [0, 25, 50, 75, 100])
        rec.update(min=float(q[0]), p25=float(q[1]), median=float(q[2]),
                   p75=float(q[3]), max=float(q[4]), mean=float(ok.mean()),
                   n_gt_one=int((ok > 1.0).sum()), n_eq_one=int((ok == 1.0).sum()))
    return rec


# ==========================================================================
# 정렬 파일 읽기 — 헤더가 있을 수도, 없을 수도 있다
# ==========================================================================
def looks_like_header(tokens) -> tuple[bool, str]:
    """첫 행이 헤더인지 판정한다.

    common.read_tsv 를 쓰지 않는 이유 — 그 함수는 '컬럼명에 | 또는 숫자' 로 헤더 유무를
    추정한다. convertalis 산출물은 보통 숫자를 포함하므로 대개 맞지만, 맞춘 뒤에는 컬럼이
    col0..colN 이 되어 어느 컬럼이 qaln 인지 알 수 없다. 여기서는 (가) 헤더가 있으면
    그대로 쓰고 (나) 없으면 --columns 를 요구하는 두 갈래를 **명시적으로** 가른다.

    판정 규칙 — 모든 토큰이 맨 식별자([A-Za-z_][A-Za-z0-9_]*) 이고, 그중 둘 이상이
    mmseqs 필드 어휘에 있으면 헤더다. 데이터 행은 fident/qstart 같은 숫자를 반드시
    포함하므로 첫 조건에서 걸러진다.
    """
    toks = [str(t).strip() for t in tokens]
    if not toks:
        return False, "빈 행"
    ident = [bool(re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", t)) for t in toks]
    if not all(ident):
        bad = [t for t, k in zip(toks, ident) if not k][:5]
        return False, f"식별자가 아닌 토큰 {bad} 이 있어 데이터 행으로 본다"
    known = [t for t in toks if t.lower() in HEADER_FIELD_TOKENS]
    if len(known) >= 2:
        return True, f"mmseqs 필드 어휘 일치 {len(known)}/{len(toks)}"
    return False, f"필드 어휘 일치가 {len(known)}개뿐이라 데이터 행으로 본다"


def read_alignment(path: Path, columns: list[str] | None) -> tuple[pd.DataFrame, dict]:
    """정렬 TSV 를 읽는다. 헤더 유무를 가려내고 --columns 를 우선한다."""
    raw = pd.read_csv(path, sep="\t", header=None, dtype=str,
                      keep_default_na=False, na_filter=False)
    if raw.empty:
        raise ValueError(f"{path} 에서 행을 하나도 읽지 못했다. 파일이 비었거나 구분자가 탭이 아니다.")

    first = [str(v).strip() for v in raw.iloc[0].tolist()]
    hdr_guess, hdr_reason = looks_like_header(first)
    meta = dict(path=str(path), n_raw_rows=int(len(raw)), n_fields=int(raw.shape[1]),
                first_row=first[:30], header_guess=hdr_guess, header_reason=hdr_reason)

    def _reject_duplicate_names(names, where: str) -> None:
        """컬럼명이 중복되면 df[name] 이 Series 가 아니라 DataFrame 이 되어 한참 뒤
        `'DataFrame' object has no attribute 'str'` 같은 엉뚱한 곳에서 터진다.
        어느 이름이 몇 번 나왔는지 여기서 바로 알려준다."""
        vc = pd.Series([str(n).strip().lower() for n in names]).value_counts()
        dup = vc[vc > 1]
        if len(dup):
            raise ValueError(
                f"{where} 에 같은 컬럼명이 두 번 이상 나온다: "
                f"{ {k: int(v) for k, v in dup.items()} }. "
                f"어느 열이 qaln/taln 인지 확정할 수 없으므로 중단한다. "
                f"mmseqs --format-output 에 같은 필드를 중복으로 적었는지 확인하고, "
                f"--columns 로 열 순서를 직접 줄 것. 실제 이름 = {list(names)[:30]}")

    if columns:
        if len(columns) != raw.shape[1]:
            raise ValueError(
                f"--columns 가 {len(columns)}개인데 {path} 는 {raw.shape[1]}열이다. "
                f"--format-output 순서를 그대로 적었는지 확인할 것. "
                f"첫 행 = {first[:30]}")
        _reject_duplicate_names(columns, "--columns")
        df = raw.copy()
        df.columns = list(columns)
        dropped = False
        if all(a.lower() == b.lower() for a, b in zip(first, columns)):
            df = df.iloc[1:].reset_index(drop=True)      # 헤더가 붙어 있었다
            dropped = True
        meta.update(column_source="--columns", header_row_dropped=dropped,
                    columns=list(columns))
    else:
        if not hdr_guess:
            raise ValueError(
                f"{path} 에 헤더가 없다고 판정했다 ({hdr_reason}). {raw.shape[1]}열이며 "
                f"첫 행 = {first[:30]}.\n"
                f"어느 컬럼이 qaln/taln 인지 알 수 없으므로 --columns 로 mmseqs "
                f"--format-output 순서를 직접 줄 것. 예:\n"
                f"  --columns query,target,fident,alnlen,qstart,qend,tstart,tend,"
                f"qaln,taln,evalue,bits,qcov,tcov")
        _reject_duplicate_names(first, f"{path} 의 헤더 행")
        df = raw.iloc[1:].reset_index(drop=True)
        df.columns = first
        meta.update(column_source="파일 헤더", header_row_dropped=True, columns=first)

    if df.empty:
        raise ValueError(f"{path} 에 헤더만 있고 데이터 행이 없다.")
    df = df.reset_index(drop=True)
    df["_row"] = np.arange(len(df), dtype=np.int64)      # 파일 등장 순서 (최종 tiebreak)
    meta["n_rows"] = int(len(df))
    return df, meta


def resolve(df: pd.DataFrame) -> dict:
    """논리 이름 -> 실제 컬럼명. 필수 컬럼 확인은 common.resolve_columns 가 한다."""
    colmap = common.resolve_columns(df, list(REQUIRED_KEYS) + list(OPTIONAL_KEYS),
                                    required=REQUIRED_KEYS)
    for key, cands in EXTRA_COLSPEC.items():
        colmap[key] = common.detect_column(df, cands, required=False, what=key)
    return colmap


# ==========================================================================
# universe (분모) — 정렬 파일에서 몇 개가 빠졌는지 세려면 전체 집합이 필요하다
# ==========================================================================
def load_universe(master: Path, fasta: Path, *, max_examples: int) -> tuple[list[str] | None, dict]:
    """master_1827.tsv 또는 IPR FASTA 에서 accession 순서를 읽는다."""
    meta: dict = dict(source=None, master=str(master), fasta=str(fasta),
                      master_exists=master.exists(), fasta_exists=fasta.exists())
    if master.exists():
        m = common.read_tsv(master)
        col = common.detect_column(m, common.COLSPEC["accession"], required=False,
                                   what="accession")
        if col is None:
            meta["master_error"] = (f"master 에 accession 컬럼이 없다 (실제 컬럼 {list(m.columns)}). "
                                    f"FASTA 로 넘어간다.")
        else:
            accs = norm_acc_series(m[col]).tolist()
            dups = pd.Series(accs).duplicated()
            uniq = [a for a, d in zip(accs, dups) if not d and a]
            meta.update(source="master", source_path=str(master), column=col,
                        n_rows=int(len(m)), n_unique=len(uniq),
                        n_duplicate_rows=int(dups.sum()),
                        n_empty_accession=int(sum(1 for a in accs if not a)))
            return uniq, meta
    if fasta.exists():
        accs = [safe_parse_accession(h) for h, _ in common.read_fasta(fasta)]
        dups = pd.Series(accs).duplicated()
        uniq = [a for a, d in zip(accs, dups) if not d and a]
        meta.update(source="fasta", source_path=str(fasta), n_rows=len(accs),
                    n_unique=len(uniq), n_duplicate_rows=int(dups.sum()),
                    n_empty_accession=int(sum(1 for a in accs if not a)))
        return uniq, meta
    meta["error"] = ("master 도 FASTA 도 없다. universe 를 모르므로 '정렬 파일에서 빠진 "
                     "accession' 을 셀 수 없다 — 핵심 발견을 놓치는 상태다.")
    return None, meta


# ==========================================================================
# (A) 감사 — 파싱 전에 입력을 의심한다
# ==========================================================================
def audit_alignment(df: pd.DataFrame, colmap: dict, universe: list[str] | None,
                    *, args, read_meta: dict, universe_meta: dict) -> dict:
    """정렬 파일 자체의 감사. 여기서 나오는 숫자가 (B) 이후 해석의 전제다."""
    discrepancies: list[str] = []
    acc = norm_acc_series(df[colmap["accession"]])
    df["accession"] = acc

    n_rows = int(len(df))
    uniq = acc.drop_duplicates().tolist()
    n_uniq = len(uniq)
    n_empty_acc = int((acc == "").sum())
    if n_empty_acc:
        raise ValueError(
            f"정렬 파일에 accession 을 뽑지 못한 행이 {n_empty_acc}개 있다. "
            f"조인이 조용히 실패하므로 중단한다. accession 컬럼 = {colmap['accession']}")

    # --- 질의 단일성. 여러 질의가 섞인 파일이면 ChCODH2 좌표계 가정이 깨진다.
    query_info: dict = dict(column=colmap.get("query"))
    if colmap.get("query"):
        q = norm_acc_series(df[colmap["query"]])
        qvals = sorted(q.unique().tolist())
        query_info.update(n_distinct=len(qvals), values=qvals[:20],
                          values_truncated=len(qvals) > 20)
        if args.query_acc:
            keep = q == safe_parse_accession(args.query_acc)
            query_info.update(filtered_to=args.query_acc, n_kept=int(keep.sum()),
                              n_dropped=int((~keep).sum()))
            if not keep.any():
                raise ValueError(
                    f"--query-acc {args.query_acc} 에 해당하는 행이 없다. "
                    f"파일의 질의 값 = {qvals[:20]}")
            df = df[keep].reset_index(drop=True)
            acc = df["accession"]
            n_rows = int(len(df))
            uniq = acc.drop_duplicates().tolist()
            n_uniq = len(uniq)
        elif len(qvals) > 1:
            raise ValueError(
                f"query 컬럼({colmap['query']})에 질의가 {len(qvals)}종 섞여 있다: {qvals[:10]}. "
                f"ChCODH2 단일 질의 좌표계 가정이 깨지므로 중단한다. "
                f"--query-acc {args.reference_acc} 로 골라낼 것.")
        if len(qvals) == 1 and qvals[0] != safe_parse_accession(args.reference_acc):
            discrepancies.append(
                f"query 값이 '{qvals[0]}' 로 참조 accession '{args.reference_acc}' 와 다르다. "
                f"같은 단백질의 다른 표기일 수 있으나 확인이 필요하다 "
                f"(검증 1이 통과하면 서열은 ChCODH2 와 동일하다는 뜻이다).")
    else:
        query_info.update(n_distinct=None,
                          note="query 컬럼이 없다. 단일 질의라고 가정한다 (검증 1이 그 가정을 검사한다).")

    # --- 행 수. 과제 (A) 의 핵심 발견 지점.
    truncation: dict = dict(expected_n=args.expected_n, n_rows=n_rows, n_unique=n_uniq)
    missing_examples: list[str] = []
    if universe is not None:
        uset = set(universe)
        aset = set(uniq)
        missing = [a for a in universe if a not in aset]
        extra = [a for a in uniq if a not in uset]
        ex, n_missing, trunc = cap_list(missing, args.max_examples)
        missing_examples = ex
        truncation.update(
            n_universe=len(universe), universe_source=universe_meta.get("source"),
            n_missing=n_missing, missing_examples=ex, missing_examples_truncated=trunc,
            n_only_in_alignment=len(extra),
            only_in_alignment_examples=cap_list(extra, args.max_examples)[0],
            join_rate=(n_uniq - len(extra)) / len(universe) if universe else float("nan"),
        )
        if n_missing:
            discrepancies.append(
                f"**정렬 파일에 {n_missing}/{len(universe)} accession 이 없다** "
                f"({pct(n_missing, len(universe))} 결측). mmseqs search 기본값 "
                f"--max-seqs {MMSEQS_DEFAULT_MAX_SEQS} / e-value {MMSEQS_DEFAULT_EVALUE} "
                f"에서 상동체가 잘렸을 가능성이 높다 (인계문서 2절 — search 명령 미기록). "
                f"예시 {ex}")
            discrepancies.append(
                f"빠진 {n_missing}개는 fident 가 없으므로 **cf(클러스터 median fident)와 "
                f"cband 도 같이 빠진다.** 그 클러스터의 cf 는 남은 멤버만으로 계산되어 "
                f"조용히 편향되고, 멤버가 전부 빠진 클러스터는 cf 가 아예 없다. "
                f"단계 4 이후에서 cf 의 분모를 반드시 함께 보고할 것.")
        if extra:
            discrepancies.append(
                f"universe({universe_meta.get('source')})에 없고 정렬 파일에만 있는 accession "
                f"{len(extra)}개. 버리지 않고 행렬 뒤에 붙인다. 예시 "
                f"{cap_list(extra, args.max_examples)[0]}")
    else:
        discrepancies.append(universe_meta.get("error", "universe 를 읽지 못했다."))

    # --- universe 자체의 행 수. 규칙 3 — 1,827행은 끝까지 1,827행이어야 한다.
    #     master 에 정규화 후 중복이나 빈 accession 이 있으면 universe 가 조용히 줄어들고,
    #     그대로 residue_matrix 의 행 수와 position_coverage 의 n_universe(=보존도의 분모)가 된다.
    #     그 분모가 틀리면 뒤의 모든 보존도 비율이 조용히 틀리므로 숫자로 못박아 둔다.
    if universe is not None:
        n_uni_rows = universe_meta.get("n_rows")
        n_uni_dup = int(universe_meta.get("n_duplicate_rows") or 0)
        n_uni_empty = int(universe_meta.get("n_empty_accession") or 0)
        truncation.update(universe_n_source_rows=n_uni_rows,
                          universe_n_duplicate_rows=n_uni_dup,
                          universe_n_empty_accession=n_uni_empty)
        if n_uni_dup or n_uni_empty:
            discrepancies.append(
                f"universe 출처({universe_meta.get('source')}, {n_uni_rows}행)에서 "
                f"정규화 후 중복 {n_uni_dup}행 · 빈 accession {n_uni_empty}행을 제외해 "
                f"universe 가 {len(universe)}개로 줄었다. **residue_matrix 의 행 수와 "
                f"position_coverage 의 n_universe(보존도의 분모)가 그만큼 작아진다.** "
                f"master 의 accession 컬럼을 확인할 것 (규칙 3).")
        if args.expected_n and len(universe) != args.expected_n:
            truncation["universe_shortfall"] = args.expected_n - len(universe)
            discrepancies.append(
                f"universe 가 {len(universe)}개로 기대 {args.expected_n}개와 다르다 "
                f"({args.expected_n - len(universe):+d}). 보존도의 분모가 되는 수이므로 "
                f"master({universe_meta.get('source_path')})가 1,827 서열 전체인지 확인할 것.")

    if args.expected_n and n_uniq < args.expected_n:
        truncation["shortfall"] = args.expected_n - n_uniq
        discrepancies.append(
            f"유일 accession {n_uniq} < 기대 {args.expected_n} — {args.expected_n - n_uniq}개 부족.")
    if args.expected_n and n_uniq > args.expected_n:
        discrepancies.append(
            f"유일 accession {n_uniq} > 기대 {args.expected_n}. 정렬 DB 가 IPR010047 1,827 "
            f"서열이 아니었을 수 있다.")
    if n_uniq in (MMSEQS_DEFAULT_MAX_SEQS, MMSEQS_DEFAULT_MAX_SEQS + 1) or \
       n_rows in (MMSEQS_DEFAULT_MAX_SEQS, MMSEQS_DEFAULT_MAX_SEQS + 1):
        truncation["max_seqs_suspicion"] = True
        discrepancies.append(
            f"행 수/유일 accession 수가 {MMSEQS_DEFAULT_MAX_SEQS} 와 거의 같다. "
            f"mmseqs search 의 --max-seqs 기본값에 정확히 걸린 모양이다. "
            f"search 를 --max-seqs 2000 이상으로 다시 돌릴 것을 강력히 권한다.")

    # --- HSP 중복
    vc = acc.value_counts()
    dup_accs = vc[vc > 1]
    hsp = dict(
        n_accession_with_multiple_hsp=int(len(dup_accs)),
        n_rows_in_duplicate_accession=int(dup_accs.sum()),
        n_rows_to_merge=int(dup_accs.sum() - len(dup_accs)),
        max_hsp_per_accession=int(vc.max()) if len(vc) else 0,
        tiebreak_rule=HSP_TIEBREAK_DESC,
        examples=[dict(accession=a, n_hsp=int(c)) for a, c in dup_accs.head(args.max_examples).items()],
        examples_truncated=bool(len(dup_accs) > args.max_examples),
    )
    if len(dup_accs):
        discrepancies.append(
            f"한 accession 에 HSP 가 여러 개인 경우 {len(dup_accs)}개 "
            f"(행 {int(dup_accs.sum())}개 → {len(dup_accs)}개로 병합, "
            f"{hsp['n_rows_to_merge']}행 버림). 규칙: {HSP_TIEBREAK_DESC}")

    # --- 수치 분포
    dist = {k: numeric_summary(df[colmap[k]] if colmap.get(k) else None, k)
            for k in ("fident", "qcov", "tcov", "evalue", "bits")}
    fid = dist["fident"]
    if fid.get("present") and fid.get("n_gt_one"):
        discrepancies.append(
            f"fident 에 1.0 을 넘는 값이 {fid['n_gt_one']}개 있다. mmseqs fident 는 0~1 "
            f"분수이므로, 이 컬럼은 퍼센트(pident)일 가능성이 높다. cf/cband 의 축이 "
            f"100배 달라지므로 단계 4 이전에 반드시 확정할 것.")
    for k in ("fident", "qcov", "tcov"):
        d = dist[k]
        if d.get("present") and (d.get("n_unparseable") or d.get("n_blank")):
            discrepancies.append(
                f"{k} 에 파싱 불가 {d['n_unparseable']}개 / 빈값 {d['n_blank']}개.")
        if not d.get("present"):
            discrepancies.append(f"{k} 컬럼이 없다 (품질 참고/선정 근거가 줄어든다).")
    if dist["evalue"].get("present") and dist["evalue"].get("n_valid"):
        mx = dist["evalue"]["max"]
        dist["evalue"]["at_default_evalue_cutoff"] = bool(mx <= MMSEQS_DEFAULT_EVALUE * 1.000001)
        if dist["evalue"]["at_default_evalue_cutoff"]:
            discrepancies.append(
                f"최대 e-value 가 {mx:g} 로 기본 절단값 {MMSEQS_DEFAULT_EVALUE} 이하다. "
                f"-e 를 올리지 않은 기본 search 였다는 증거이며, 먼 상동체가 버려졌을 수 있다.")

    # --- 자기 자신 히트
    ref_acc = safe_parse_accession(args.reference_acc)
    self_rows = df[df["accession"] == ref_acc]
    self_info = dict(reference_acc=ref_acc, n_rows=int(len(self_rows)))
    if len(self_rows):
        if colmap.get("fident"):
            fvals = pd.to_numeric(self_rows[colmap["fident"]], errors="coerce").tolist()
            self_info["fident_values"] = [None if pd.isna(v) else float(v) for v in fvals]
            self_info["fident_is_one"] = bool(any(
                (not pd.isna(v)) and (abs(v - 1.0) < 1e-9 or abs(v - 100.0) < 1e-9) for v in fvals))
            if not self_info["fident_is_one"]:
                discrepancies.append(
                    f"자기 자신 히트({ref_acc})의 fident 가 {self_info['fident_values']} 로 "
                    f"1.0(또는 100) 이 아니다. 질의 서열이 ChCODH2 가 아닐 수 있다 "
                    f"(인계문서 2절 — 질의 서열 출처 미확인).")
        for k in ("qcov", "tcov"):
            if colmap.get(k):
                self_info[f"{k}_values"] = [
                    None if pd.isna(v) else float(v) for v in
                    pd.to_numeric(self_rows[colmap[k]], errors="coerce").tolist()]
    else:
        discrepancies.append(
            f"자기 자신 히트({ref_acc})가 정렬 파일에 없다. 검증 3(q_pos == t_pos 636개)을 "
            f"할 수 없으므로 좌표계의 가장 강한 그물이 하나 빠진다. "
            f"search DB 에 질의 자신이 포함되지 않았다는 뜻이다.")

    return dict(
        read=read_meta, universe=universe_meta, columns=colmap, query=query_info,
        truncation=truncation, hsp=hsp, distributions=dist, self_hit=self_info,
        n_rows=n_rows, n_unique_accession=n_uniq, missing_examples=missing_examples,
        discrepancies=discrepancies, dataframe=df,
    )


# ==========================================================================
# HSP 선택 — 한 accession 당 한 정렬
# ==========================================================================
def select_hsps(df: pd.DataFrame, colmap: dict, *, max_detail: int) -> tuple[pd.DataFrame, dict]:
    """accession 당 HSP 하나를 고른다. 어떤 키로 골랐는지와 몇 건을 버렸는지 보고한다."""
    keys, asc, used = [], [], []
    work = df.copy()
    for key, ascending in HSP_TIEBREAK:
        if colmap.get(key):
            work[f"_k_{key}"] = pd.to_numeric(work[colmap[key]], errors="coerce")
            keys.append(f"_k_{key}")
            asc.append(ascending)
            used.append(f"{key}({'asc' if ascending else 'desc'})")
    if not keys:
        # bits/evalue/qcov 가 전부 없다 — 정렬된 잔기 수가 많은 쪽을 쓴다 (명시적 대체 규칙)
        work["_k_nres"] = work[colmap["qaln"]].astype(str).str.replace(
            DEFAULT_GAP_CHAR, "", regex=False).str.len()
        keys, asc, used = ["_k_nres"], [False], ["qaln_residues(desc)"]

    ordered = work.sort_values(by=keys + ["_row"], ascending=asc + [True], kind="stable")
    first = ordered.drop_duplicates(subset=["accession"], keep="first")
    chosen_rows = set(first["_row"].tolist())
    dropped = ordered[~ordered["_row"].isin(chosen_rows)]

    detail = []
    if len(dropped):
        for a, grp in ordered[ordered["accession"].isin(set(dropped["accession"]))].groupby(
                "accession", sort=True):
            rec = dict(accession=a, n_hsp=int(len(grp)),
                       chosen_row=int(grp.iloc[0]["_row"]),
                       dropped_rows=[int(v) for v in grp.iloc[1:]["_row"].tolist()])
            for key in ("bits", "evalue", "qcov"):
                if colmap.get(key):
                    rec[key] = [None if pd.isna(v) else float(v)
                                for v in pd.to_numeric(grp[colmap[key]], errors="coerce").tolist()]
            detail.append(rec)
            if len(detail) >= max_detail:
                break

    sel = first.sort_values("_row", kind="stable").reset_index(drop=True)
    meta = dict(
        n_rows_in=int(len(df)), n_rows_selected=int(len(sel)),
        n_rows_dropped=int(len(df) - len(sel)),
        tiebreak_keys_used=used, tiebreak_rule=HSP_TIEBREAK_DESC,
        detail=detail, detail_listed=len(detail),
        detail_truncated=bool(len(dropped) and len(detail) >= max_detail),
    )
    if len(sel) != sel["accession"].nunique():
        raise AssertionError("HSP 선택 후에도 accession 중복이 남았다 — 선택 로직 결함")
    return sel, meta


# ==========================================================================
# (B) 파싱 — 좌표 기준(base) 판정과 세 가지 검증
# ==========================================================================
def parse_int(value, *, what: str, row: int) -> int:
    s = str(value).strip()
    try:
        return int(float(s))
    except (TypeError, ValueError):
        raise ValueError(f"{what} 를 정수로 읽을 수 없다 (행 {row}, 값 '{s}')") from None


def check_lengths(sel: pd.DataFrame, colmap: dict, *, max_examples: int) -> dict:
    """검증 2 — len(qaln) == len(taln). base 판정보다 먼저 한다."""
    qa = sel[colmap["qaln"]].astype(str).str.strip()
    ta = sel[colmap["taln"]].astype(str).str.strip()
    bad = qa.str.len() != ta.str.len()
    empty = (qa.str.len() == 0) | (ta.str.len() == 0)
    if empty.any():
        ex = sel.loc[empty, "accession"].head(max_examples).tolist()
        raise ValueError(
            f"qaln 또는 taln 이 빈 행이 {int(empty.sum())}개 있다. 예시 {ex}. "
            f"--format-output 에 qaln/taln 이 실제로 포함되었는지 확인할 것.")
    if bad.any():
        rows = []
        for _, r in sel[bad].head(max_examples).iterrows():
            rows.append(f"{r['accession']}: len(qaln)={len(str(r[colmap['qaln']]).strip())} "
                        f"len(taln)={len(str(r[colmap['taln']]).strip())}")
        raise ValueError(
            f"[검증 2 실패] qaln 과 taln 의 길이가 다른 행이 {int(bad.sum())}개 있다. "
            f"정렬 문자열 쌍이 깨졌으므로 잔기 대응을 만들 수 없다.\n  " + "\n  ".join(rows))
    return dict(n_checked=int(len(sel)), aln_len_min=int(qa.str.len().min()),
                aln_len_max=int(qa.str.len().max()),
                aln_len_median=float(qa.str.len().median()))


def check_query_length(df: pd.DataFrame, colmap: dict, ref: str, *,
                       max_examples: int) -> tuple[dict, list[str]]:
    """qlen 대조 — 「search 질의가 정말 이 참조서열이었나」의 가장 싼 직접 증거.

    인계문서 2절이 '질의 서열 출처 미확인' 으로 남긴 빈칸을 정렬 파일 자신이 메워 준다.
    qlen != 636 이면 좌표계의 전제가 깨진 것이므로(규칙 6) 예외다. 특히 1잔기 차이는
    0-based 오판으로 위장해 t_pos 를 조용히 +1 밀어버리므로 반드시 먼저 걸러야 한다.
    alnlen 은 정렬 열 수와 대조하되, 버전에 따라 정의가 달라 불일치로만 보고한다.
    """
    disc: list[str] = []
    info: dict = dict(reference_len=len(ref), qlen_present=bool(colmap.get("qlen")),
                      alnlen_present=bool(colmap.get("alnlen")))
    if colmap.get("qlen"):
        ql = pd.to_numeric(df[colmap["qlen"]].astype(str).str.strip(), errors="coerce")
        vals = sorted({int(v) for v in ql.dropna().tolist()})
        bad = ql.notna() & (ql != len(ref))
        info.update(qlen_column=colmap["qlen"], qlen_values=vals[:20],
                    n_qlen_unparseable=int(ql.isna().sum()),
                    n_qlen_mismatch=int(bad.sum()), n_rows=int(len(df)))
        if int(ql.isna().sum()):
            disc.append(f"qlen 을 숫자로 읽지 못한 행 {int(ql.isna().sum())}개 — 대조를 건너뛴다.")
        if bad.any():
            ex = df.loc[bad, colmap["accession"]].astype(str).head(max_examples).tolist()
            raise ValueError(
                f"정렬 파일의 qlen 이 참조서열 길이 {len(ref)} 와 다른 행이 {int(bad.sum())}/"
                f"{len(df)}개 있다 (qlen 값 {vals[:10]}). "
                f"search 질의가 pipeline/lib/chcodh2.fasta 의 {common.CHCODH2_ACC} "
                f"{len(ref)} aa 가 아니었다는 직접 증거이며, 그러면 q_pos/t_pos 전체가 "
                f"다른 좌표계가 된다 (인계문서 2절 '질의 서열 출처 미확인').\n"
                f"  특히 {len(ref) - vals[0] if vals else '?'} 잔기 차이는 0-based 오판으로 "
                f"위장해 t_pos 를 조용히 밀어버린다. 예시 accession {ex}\n"
                f"  조치 — 질의 FASTA 를 찾아 lib/chcodh2.fasta 와 대조하고, 다르면 올바른 "
                f"질의로 search 를 다시 돌릴 것.")
    else:
        disc.append("qlen 컬럼이 없다. 'search 질의가 정말 ChCODH2 636 aa 였나' 를 직접 "
                    "확인할 수 없고 검증 1·3 의 간접 증거에만 의존한다 "
                    "(--format-output 에 qlen 을 넣어 다시 뽑으면 가장 싸게 확인된다).")
    if colmap.get("alnlen"):
        al = pd.to_numeric(df[colmap["alnlen"]].astype(str).str.strip(), errors="coerce")
        cols = df[colmap["qaln"]].astype(str).str.strip().str.len()
        bad = al.notna() & (al != cols)
        info.update(alnlen_column=colmap["alnlen"], n_alnlen_mismatch=int(bad.sum()))
        if bad.any():
            disc.append(
                f"alnlen 이 len(qaln) 과 다른 행 {int(bad.sum())}/{len(df)}개. "
                f"--columns 순서가 어긋났거나 alnlen 정의가 다른 판본이다 "
                f"(매핑은 qaln/taln 만 쓰므로 좌표에는 영향 없음).")
    return info, disc


def detect_qstart_base(sel: pd.DataFrame, colmap: dict, ref: str, *,
                       assume: str, max_examples: int, upper_aln: bool) -> dict:
    """검증 1 의 전반부 — qstart 가 1-based 인지 0-based 인지 **데이터로** 판정한다.

    mmseqs 의 qstart 는 1-based 라고 가정하지만 그 가정을 믿지 않는다. 두 해석 모두로
    전 행을 맞춰보고 이긴 쪽을 채택하며, 각 base 의 일치 행 수를 채택 근거로 남긴다.
    """
    gap = DEFAULT_GAP_CHAR
    L = len(ref)
    ref_up = ref.upper()
    stat = {1: 0, 0: 0}
    case_only = {1: 0, 0: 0}
    out_of_range = {1: 0, 0: 0}
    examples: dict[int, list] = {1: [], 0: []}
    # 0-based 의 직접 증거 = 어딘가에 색인 0 이 있다는 것. 숫자로 세어 채택 근거에 넣는다.
    qstart_vals: list[int] = []
    tstart_vals: list[int] = []

    for _, r in sel.iterrows():
        qaln = str(r[colmap["qaln"]]).strip()
        if upper_aln:
            qaln = qaln.upper()
        qres = qaln.replace(gap, "")
        qstart = parse_int(r[colmap["qstart"]], what="qstart", row=int(r["_row"]))
        qstart_vals.append(qstart)
        tstart_vals.append(parse_int(r[colmap["tstart"]], what="tstart", row=int(r["_row"])))
        for base in (1, 0):
            start0 = qstart - 1 if base == 1 else qstart
            if start0 < 0 or start0 + len(qres) > L:
                out_of_range[base] += 1
                if len(examples[base]) < max_examples:
                    examples[base].append(dict(accession=r["accession"], qstart=qstart,
                                               n_residues=len(qres), reason="범위 초과"))
                continue
            seg = ref[start0:start0 + len(qres)]
            if qres == seg:
                stat[base] += 1
            elif qres.upper() == ref_up[start0:start0 + len(qres)]:
                case_only[base] += 1
                if len(examples[base]) < max_examples:
                    examples[base].append(dict(accession=r["accession"], qstart=qstart,
                                               reason="대소문자만 다름"))
            elif len(examples[base]) < max_examples:
                found = ref.find(qres) if len(qres) >= 8 else -1
                examples[base].append(dict(
                    accession=r["accession"], qstart=qstart, n_residues=len(qres),
                    expected_prefix=seg[:30], observed_prefix=qres[:30],
                    found_at_1based=(found + 1) if found >= 0 else None,
                    reason="부분서열 불일치"))

    n = int(len(sel))
    n_qstart_zero = int(sum(1 for v in qstart_vals if v == 0))
    n_tstart_zero = int(sum(1 for v in tstart_vals if v == 0))
    zero_stat = dict(
        qstart_min=(min(qstart_vals) if qstart_vals else None),
        tstart_min=(min(tstart_vals) if tstart_vals else None),
        n_qstart_zero=n_qstart_zero, n_tstart_zero=n_tstart_zero,
        zero_index_present=bool(n_qstart_zero or n_tstart_zero),
    )

    if assume in ("1", "0"):
        base = int(assume)
        basis = f"--assume-qstart-base {assume} 로 강제 지정"
    elif stat[1] == 0 and stat[0] == 0:
        raise ValueError(
            f"[검증 1 실패] 1-based / 0-based 어느 해석으로도 qaln 이 참조서열과 맞는 행이 "
            f"하나도 없다 ({n}행 검사). 질의 서열이 ChCODH2(Q9F8A8, 636 aa)가 아니거나, "
            f"qaln/taln 컬럼이 뒤바뀌었거나, 참조 FASTA 가 다른 아이소폼이다.\n"
            f"  1-based 범위초과 {out_of_range[1]} · 대소문자만다름 {case_only[1]}\n"
            f"  0-based 범위초과 {out_of_range[0]} · 대소문자만다름 {case_only[0]}\n"
            f"  1-based 예시: {json.dumps(examples[1][:3], ensure_ascii=False)}\n"
            f"  0-based 예시: {json.dumps(examples[0][:3], ensure_ascii=False)}\n"
            f"  (대소문자만 다르면 --upper-aln 으로 대문자 변환 후 진행할 수 있다)")
    elif stat[1] > stat[0]:
        base = 1
        basis = (f"1-based 로 {stat[1]}/{n} 행 일치, 0-based 로 {stat[0]}/{n} 행 일치 → "
                 f"1-based 채택 (mmseqs 기본 가정과 일치)")
    elif stat[0] > stat[1]:
        # 0-based 채택은 t_pos 에도 그대로 적용된다(tstart0 = tstart). 그런데 파일 어디에도
        # 색인 0 이 없으면 「진짜 0-based」보다 「참조서열과 질의서열이 1잔기 어긋났다」가
        # 훨씬 흔한 원인이다 (개시 Met 유무, 태그, 다른 아이소폼).
        # 그 경우 q_pos 는 맞게 나오지만 **t_pos 가 전부 +1 밀린 채 조용히 틀린다.**
        # tend/tlen 보조 확인도 같은 base 를 쓰므로 자기모순 없이 통과해 잡아주지 못하고,
        # 자기 히트(검증 3)가 없으면 그물이 하나도 남지 않는다. 그래서 여기서 중단한다.
        if not zero_stat["zero_index_present"]:
            raise ValueError(
                f"[검증 1 보류] 0-based 로 {stat[0]}/{n} 행이 맞고 1-based 로 {stat[1]}/{n} 행이 "
                f"맞지만, 파일 어디에도 색인 0 이 없다 "
                f"(qstart 최소 {zero_stat['qstart_min']} · tstart 최소 {zero_stat['tstart_min']} · "
                f"qstart==0 {n_qstart_zero}행 · tstart==0 {n_tstart_zero}행).\n"
                f"  진짜 0-based 산출물이라면 {n}개 국소 정렬 중 어딘가에 0 이 나와야 한다. "
                f"0 이 전혀 없다는 것은 base 가 아니라 **참조서열과 search 질의서열이 1잔기 "
                f"어긋났다**는 뜻일 가능성이 훨씬 높다 (개시 Met 유무·태그·다른 아이소폼. "
                f"인계문서 2절 '질의 서열 출처 미확인').\n"
                f"  그대로 0-based 로 진행하면 q_pos 는 맞게 나오지만 t_pos 가 전부 +1 밀려 "
                f"**조용히 틀린다**(tend/tlen 보조 확인은 같은 base 를 쓰므로 잡아내지 못하고, "
                f"자기 히트가 없으면 검증 3도 못 한다).\n"
                f"  조치 — ① 정렬 파일의 qlen 이 {L} 인지 본다 ② search 질의 FASTA 를 "
                f"pipeline/lib/chcodh2.fasta 와 대조한다 ③ 진짜 0-based 임을 확인했으면 "
                f"--assume-qstart-base 0 으로 강제한다.")
        base = 0
        basis = (f"**0-based 로 {stat[0]}/{n} 행 일치, 1-based 로 {stat[1]}/{n} 행 일치 → "
                 f"0-based 채택.** mmseqs 의 1-based 가정이 깨졌다 — 산출물의 출처를 확인할 것 "
                 f"(색인 0 증거: qstart==0 {n_qstart_zero}행 · tstart==0 {n_tstart_zero}행)")
    else:
        base = 1
        basis = (f"두 해석의 일치 행 수가 같다({stat[1]}/{n}). 1-based 를 채택하지만 "
                 f"근거가 약하다")

    return dict(base=base, basis=basis, n_checked=n, **zero_stat,
                n_match_1based=stat[1], n_match_0based=stat[0],
                n_case_only_1based=case_only[1], n_case_only_0based=case_only[0],
                n_out_of_range_1based=out_of_range[1], n_out_of_range_0based=out_of_range[0],
                ambiguous=bool(stat[1] == stat[0] and assume == "auto"),
                mismatch_examples_1based=examples[1], mismatch_examples_0based=examples[0],
                assume=assume, upper_aln=bool(upper_aln))


def build_mapping(sel: pd.DataFrame, colmap: dict, ref: str, base: int, *,
                  rows: list[str], row_index: dict, args) -> dict:
    """정렬 문자열을 ChCODH2 좌표로 투영한다. 검증 1을 행마다 다시 확인한다."""
    gap = DEFAULT_GAP_CHAR
    outside = args.outside_char
    L = len(ref)
    n = len(rows)

    mat = np.full((n, L), outside, dtype="<U1")
    cov = np.zeros((n, L), dtype=bool)
    in_aln = np.zeros(n, dtype=bool)

    equiv: list[tuple] = []
    per_acc: list[dict] = []
    mismatches: list[str] = []
    n_scanned = 0
    ins_cols_total = 0
    ins_res_total = 0
    del_cols_total = 0
    end_mismatch = dict(qend=0, tend=0, tlen_exceeded=0)
    nonstd_tchars: dict[str, int] = {}

    for _, r in sel.iterrows():
        n_scanned += 1
        acc = r["accession"]
        row = int(r["_row"])
        qaln = str(r[colmap["qaln"]]).strip()
        taln = str(r[colmap["taln"]]).strip()
        if args.upper_aln:
            qaln, taln = qaln.upper(), taln.upper()
        qstart = parse_int(r[colmap["qstart"]], what="qstart", row=row)
        tstart = parse_int(r[colmap["tstart"]], what="tstart", row=row)

        qarr = np.array(list(qaln), dtype="<U1")
        tarr = np.array(list(taln), dtype="<U1")
        qnon = qarr != gap
        tnon = tarr != gap
        nq = int(qnon.sum())
        nt = int(tnon.sum())

        start0 = qstart - 1 if base == 1 else qstart
        tstart0 = tstart - 1 if base == 1 else tstart
        if start0 < 0 or start0 + nq > L:
            raise ValueError(
                f"[검증 1 실패] {acc}: qstart={qstart} (base={base} → 0-based {start0}) 에서 "
                f"잔기 {nq}개는 참조서열 {L} aa 를 벗어난다. base 판정 또는 입력이 잘못됐다.")
        if tstart0 < 0:
            raise ValueError(f"{acc}: tstart={tstart} 가 base={base} 해석에서 음수 색인이 된다.")

        qres = "".join(qarr[qnon])
        seg = ref[start0:start0 + nq]
        if qres != seg:
            diff = next((i for i in range(nq) if qres[i] != seg[i]), 0)
            found = ref.find(qres) if nq >= 8 else -1
            mismatches.append(
                f"{acc} (파일 {row}행): qstart={qstart}, 잔기 {nq}개. "
                f"첫 불일치 색인 {diff} — 참조 '{seg[diff]}'(위치 {start0 + diff + 1}) vs "
                f"qaln '{qres[diff]}'. 기대 '{seg[:30]}' vs 관측 '{qres[:30]}'" +
                (f". 이 부분서열은 참조의 {found + 1} 번 위치에서 시작한다 "
                 f"(qstart 와 {found + 1 - (start0 + 1)} 만큼 어긋남)" if found >= 0 else ""))
            if len(mismatches) >= args.max_examples:
                break
            continue

        ridx = row_index[acc]
        in_aln[ridx] = True
        qidx0 = start0 + np.cumsum(qnon) - 1
        tidx0 = tstart0 + np.cumsum(tnon) - 1
        qcols = np.flatnonzero(qnon)
        mat[ridx, qidx0[qcols]] = tarr[qcols]
        cov[ridx, qidx0[qcols]] = True

        both = np.flatnonzero(qnon & tnon)
        equiv.append((acc, (qidx0[both] + 1).astype(np.int32), qarr[both],
                      (tidx0[both] + 1).astype(np.int32), tarr[both]))

        n_del = int((qnon & ~tnon).sum())
        n_ins_cols = int((~qnon & tnon).sum())
        del_cols_total += n_del
        ins_cols_total += n_ins_cols
        ins_res_total += n_ins_cols
        for ch, cnt in zip(*np.unique(tarr[tnon], return_counts=True)):
            if str(ch) not in "ACDEFGHIKLMNPQRSTVWY":
                nonstd_tchars[str(ch)] = nonstd_tchars.get(str(ch), 0) + int(cnt)

        q_from, q_to = start0 + 1, start0 + nq
        t_from, t_to = tstart0 + 1, tstart0 + nt
        if colmap.get("qend"):
            qe = parse_int(r[colmap["qend"]], what="qend", row=row)
            if (qe if base == 1 else qe + 1) != q_to:
                end_mismatch["qend"] += 1
        if colmap.get("tend"):
            te = parse_int(r[colmap["tend"]], what="tend", row=row)
            if (te if base == 1 else te + 1) != t_to:
                end_mismatch["tend"] += 1
        if colmap.get("tlen"):
            tl = parse_int(r[colmap["tlen"]], what="tlen", row=row)
            if t_to > tl:
                end_mismatch["tlen_exceeded"] += 1

        per_acc.append(dict(accession=acc, q_from=q_from, q_to=q_to,
                            t_from=t_from, t_to=t_to, aln_cols=int(len(qarr)),
                            n_q_residues=nq, n_t_residues=nt,
                            n_equiv_pairs=int(both.size), n_deletion_cols=n_del,
                            n_insertion_cols=n_ins_cols))

    if mismatches:
        # 상한에 걸려 탐색을 멈췄을 수 있다 — 검사한 행 수와 남은 행 수를 반드시 적는다 (규칙 7).
        stopped_early = len(mismatches) >= args.max_examples and n_scanned < len(sel)
        raise ValueError(
            f"[검증 1 실패] qaln 에서 갭을 뺀 문자열이 참조서열의 qstart 부분서열과 다른 행이 "
            f"{len(mismatches)}개 있다 (base={base} 채택). 1-based/0-based 오프셋 오류이거나 "
            f"질의 서열이 ChCODH2 가 아니다. 뒤의 모든 숫자가 조용히 틀리므로 중단한다. "
            f"검사한 행 {n_scanned}/{len(sel)}"
            + (f" — 불일치 {args.max_examples}건(--max-examples)에서 탐색을 멈췄으므로 "
               f"나머지 {len(sel) - n_scanned}행은 검사하지 않았다. 실제 불일치는 더 많을 수 있다"
               if stopped_early else " (전 행 검사 완료)")
            + ":\n  " + "\n  ".join(mismatches))

    # --- 상동체 쪽 알파벳 점검. 행렬의 알파벳은 하위 단계(보존도)와의 계약이다.
    #     taln 만 소문자인 입력(소프트 마스킹)은 검증 1·2·3 을 전부 통과하면서
    #     행렬에 소문자를 그대로 심는다. 그러면 'a' != 'A' 라서 n_identical 과
    #     하위 단계의 보존도가 **조용히 과소계산**된다. 그래서 경고가 아니라 예외다(규칙 6).
    disc: list[str] = []
    lower_chars = {c: n for c, n in nonstd_tchars.items() if c.isalpha() and c.islower()}
    if lower_chars:
        tot = sum(lower_chars.values())
        raise ValueError(
            f"taln 에 소문자 잔기가 {tot}개 있다 (문자 {len(lower_chars)}종, "
            f"{dict(sorted(lower_chars.items())[:8])}...). qaln 은 대문자이므로 검증 1~3 은 "
            f"전부 통과하지만, 행렬에 소문자가 그대로 들어가면 'a' != 'A' 가 되어 "
            f"**n_identical 과 하위 단계의 보존도가 조용히 과소계산된다.**\n"
            f"  소프트 마스킹된 산출물이면 --upper-aln 을 주고 다시 돌릴 것 "
            f"(qaln/taln 을 대문자로 올려 비교·매핑한다).")
    other_nonstd = {c: n for c, n in nonstd_tchars.items() if not (c.isalpha() and c.islower())}
    if other_nonstd:
        disc.append(
            f"상동체 쪽 비표준 문자 {sum(other_nonstd.values())}개 "
            f"({json.dumps(other_nonstd, ensure_ascii=False)}). 행렬에 그대로 들어가며 "
            f"보존도의 분모(n_residue)에는 세지되 n_identical 에는 들어가지 않는다. "
            f"하위 단계에서 X/*/U 등을 어떻게 다룰지 합의할 것.")

    return dict(mat=mat, cov=cov, in_aln=in_aln, equiv=equiv, per_acc=per_acc,
                insertion_cols=ins_cols_total, insertion_residues=ins_res_total,
                deletion_cols=del_cols_total, end_mismatch=end_mismatch,
                nonstandard_target_chars=nonstd_tchars, discrepancies=disc)


def check_self_hit(mapping: dict, ref: str, *, row_index: dict, ref_acc: str,
                   allow_partial: bool) -> tuple[dict, list[str]]:
    """검증 3 — 자기 자신 히트는 q_pos == t_pos 가 636개 전부 성립해야 한다."""
    disc: list[str] = []
    info = dict(reference_acc=ref_acc, present=False)
    if ref_acc not in row_index:
        disc.append(f"[검증 3 생략] 자기 자신 히트({ref_acc})가 정렬에 없다.")
        return info, disc
    ridx = row_index[ref_acc]
    if not mapping["in_aln"][ridx]:
        disc.append(f"[검증 3 생략] {ref_acc} 행이 정렬 매핑에 들어가지 않았다.")
        return info, disc

    pairs = next((e for e in mapping["equiv"] if e[0] == ref_acc), None)
    if pairs is None:
        disc.append(f"[검증 3 생략] {ref_acc} 의 등가 쌍이 비었다.")
        return info, disc

    _, qpos, qaa, tpos, taa = pairs
    n_pairs = int(qpos.size)
    n_same = int((qpos == tpos).sum())
    n_aa_same = int((qaa == taa).sum())
    matrix_row = "".join(mapping["mat"][ridx])
    info.update(present=True, n_equiv_pairs=n_pairs, n_q_eq_t=n_same,
                n_same_residue=n_aa_same, expected=len(ref),
                matrix_row_equals_reference=bool(matrix_row == ref))
    ok = (n_pairs == len(ref) and n_same == len(ref) and n_aa_same == len(ref)
          and matrix_row == ref)
    info["passed"] = bool(ok)
    if not ok:
        msg = (f"[검증 3 실패] 자기 자신 히트({ref_acc})가 완전 일치가 아니다: "
               f"등가 쌍 {n_pairs}/{len(ref)} · q_pos==t_pos {n_same}/{len(ref)} · "
               f"q_aa==t_aa {n_aa_same}/{len(ref)} · 행렬행==참조서열 "
               f"{matrix_row == ref}. 좌표계가 어긋났거나 자기 히트가 국소 정렬로 잘렸다.")
        if allow_partial:
            disc.append(msg + " (--allow-partial-self-hit 로 계속 진행)")
        else:
            raise AssertionError(
                msg + "\n  자기 히트가 실제로 부분 정렬이라면 --allow-partial-self-hit 를 "
                      "주고, 그 사실을 보고서에 반드시 적을 것.")
    return info, disc


# ==========================================================================
# (C) 산출 테이블
# ==========================================================================
def write_equiv_map(path: Path, equiv: list[tuple]) -> dict:
    """long 테이블. 1,827 x ~630 = 100만행대라서 accession 단위로 흘려 쓴다."""
    n = 0
    with open(path, "wt", encoding="utf-8") as fh:
        fh.write("\t".join(EQUIV_COLUMNS) + "\n")
        for acc, qpos, qaa, tpos, taa in equiv:
            fh.writelines([f"{acc}\t{q}\t{a}\t{t}\t{b}\n"
                           for q, a, t, b in zip(qpos.tolist(), qaa.tolist(),
                                                 tpos.tolist(), taa.tolist())])
            n += int(qpos.size)
    return dict(path=str(path), n_rows=n, columns=list(EQUIV_COLUMNS))


def write_residue_matrix(path: Path, mat: np.ndarray, rows: list[str]) -> dict:
    """wide 테이블. DataFrame 을 만들지 않고 직접 쓴다 (1,827 x 637 객체열 회피)."""
    labels = [str(i) for i in range(1, mat.shape[1] + 1)]
    with open(path, "wt", encoding="utf-8") as fh:
        fh.write("accession\t" + "\t".join(labels) + "\n")
        for i, acc in enumerate(rows):
            fh.write(acc + "\t" + "\t".join(mat[i].tolist()) + "\n")
    return dict(path=str(path), n_rows=len(rows), n_position_columns=len(labels))


def write_known_positions(path: Path, mat: np.ndarray, cov: np.ndarray,
                          in_aln: np.ndarray, rows: list[str], ref: str) -> dict:
    """클러스터 특성화에 바로 쓰는 표. 열 이름에 참조 잔기를 넣는다 ('pos41_F')."""
    positions = sorted(common.KNOWN_POSITIONS)
    cols = [f"pos{p}_{common.KNOWN_POSITIONS[p]['aa']}" for p in positions]
    idx = np.array([p - 1 for p in positions], dtype=np.int64)
    sub = mat[:, idx]
    covered = (cov[:, idx] & (sub != DEFAULT_GAP_CHAR)).sum(axis=1)

    out = pd.DataFrame(sub, columns=cols)
    out.insert(0, "accession", rows)
    out["in_alignment"] = in_aln
    out["n_known_covered"] = covered.astype(np.int32)
    out.to_csv(path, sep="\t", index=False)

    ref_mismatch = [p for p in positions if ref[p - 1] != common.KNOWN_POSITIONS[p]["aa"]]
    if ref_mismatch:
        raise AssertionError(f"참조 잔기와 KNOWN_POSITIONS 가 어긋난 위치: {ref_mismatch}")
    return dict(path=str(path), n_rows=int(len(out)), positions=positions, columns=cols,
                extra_columns=["in_alignment", "n_known_covered"],
                n_rows_with_all_known=int((covered == len(positions)).sum()),
                n_rows_with_no_known=int((covered == 0).sum()))


def build_coverage(mat: np.ndarray, cov: np.ndarray, in_aln: np.ndarray,
                   ref: str) -> tuple[pd.DataFrame, dict]:
    """636행. 보존도의 분모. 분해 항등식이 깨지면 예외."""
    gap = DEFAULT_GAP_CHAR
    n_rows = mat.shape[0]
    n_in = int(in_aln.sum())
    n_not_in = n_rows - n_in

    residue_mask = cov & (mat != gap)
    n_residue = residue_mask.sum(axis=0).astype(np.int64)
    n_gap_within = (cov & (mat == gap)).sum(axis=0).astype(np.int64)
    n_outside = (in_aln[:, None] & ~cov).sum(axis=0).astype(np.int64)
    ref_arr = np.array(list(ref), dtype="<U1")
    n_ident = (residue_mask & (mat == ref_arr[None, :])).sum(axis=0).astype(np.int64)

    total = n_residue + n_gap_within + n_outside + n_not_in
    if not np.all(total == n_rows):
        bad = int(np.flatnonzero(total != n_rows)[0]) + 1
        raise AssertionError(
            f"position_coverage 분해 항등식 위반: 위치 {bad} 에서 "
            f"{int(total[bad - 1])} != 행 수 {n_rows}. 커버리지 마스크 계산이 틀렸다.")

    pos = np.arange(1, mat.shape[1] + 1)
    df = pd.DataFrame({
        "q_pos": pos,
        "q_aa": ref_arr,
        "is_known_position": [p in common.KNOWN_POSITIONS for p in pos],
        "known_role": [common.KNOWN_POSITIONS.get(p, {}).get("role", "") for p in pos],
        "n_universe": n_rows,
        "n_in_alignment": n_in,
        "n_not_in_alignment": n_not_in,
        "n_residue": n_residue,
        "n_gap_within_aln": n_gap_within,
        "n_outside_aln": n_outside,
        "n_identical": n_ident,
        "residue_frac_in_alignment": np.divide(n_residue, n_in,
                                               out=np.full(mat.shape[1], np.nan),
                                               where=(n_in > 0)),
        "residue_frac_universe": n_residue / n_rows if n_rows else np.nan,
        "identity_frac_of_residues": np.divide(n_ident, n_residue, out=np.full(mat.shape[1], np.nan),
                                               where=(n_residue > 0)),
    })[COVERAGE_COLUMNS]
    if len(df) != len(ref):
        raise AssertionError(f"position_coverage 행 수가 {len(df)} — {len(ref)} 이어야 한다")

    known = df[df["is_known_position"]]
    meta = dict(
        n_rows=int(len(df)), n_in_alignment=n_in, n_not_in_alignment=n_not_in,
        residue_total=int(n_residue.sum()),
        min_n_residue=int(n_residue.min()), max_n_residue=int(n_residue.max()),
        median_n_residue=float(np.median(n_residue)),
        n_positions_zero_residue=int((n_residue == 0).sum()),
        worst_positions=[dict(q_pos=int(p), q_aa=str(ref[p - 1]), n_residue=int(n_residue[p - 1]))
                         for p in np.argsort(n_residue)[:5] + 1],
        known_min_n_residue=int(known["n_residue"].min()) if len(known) else None,
        known_median_n_residue=float(known["n_residue"].median()) if len(known) else None,
    )
    return df, meta


# ==========================================================================
# 감사 md
# ==========================================================================
def write_audit(path: Path, audit: dict, *, args, parse: dict | None,
                outputs: dict | None, discrepancies: list[str]) -> None:
    P: list[str] = []
    stamp = datetime.now(timezone.utc).isoformat()
    P.append(f"# 단계 3 정렬 감사 — chcodh2_alignment.tsv\n")
    P.append(f"작성 {stamp} · stage `{STAGE}`\n")
    if parse is None:
        P.append("> **파싱 전에 먼저 쓴 판본이다.** 파싱이 성공하면 이 파일은 (B) 검증 결과를\n"
                 "> 포함한 완본으로 덮어쓰인다. 이 문구가 남아 있으면 파싱이 예외로 중단된 것이다.\n")
    P.append("정렬 파일을 만든 mmseqs search 명령은 **기록이 남아있지 않다**(인계문서 2절).\n"
             "그래서 이 문서는 산출물 설명이 아니라 입력 감사다.\n")

    rd, tr, hsp = audit["read"], audit["truncation"], audit["hsp"]
    cm = audit["columns"]

    P.append("## 1. 입력과 컬럼 해석\n")
    P.append(md_table(["항목", "값"], [
        ["정렬 파일", f"`{rd['path']}`"],
        ["원시 행 수", rd["n_raw_rows"]],
        ["데이터 행 수", audit["n_rows"]],
        ["열 수", rd["n_fields"]],
        ["컬럼 출처", rd.get("column_source")],
        ["헤더 판정", f"{rd['header_guess']} — {rd['header_reason']}"],
        ["헤더 행 제거", rd.get("header_row_dropped")],
        ["universe 출처", f"{audit['universe'].get('source')} "
                          f"(`{audit['universe'].get('source_path', '—')}`)"],
    ]) + "\n")
    P.append("컬럼 해석 (common.resolve_columns + 이 스크립트의 EXTRA_COLSPEC)\n")
    P.append(md_table(["논리 이름", "실제 컬럼", "필수"], [
        [k, cm.get(k) or "— 없음 —", "필수" if k in REQUIRED_KEYS else "선택"]
        for k in list(REQUIRED_KEYS) + list(OPTIONAL_KEYS) + list(EXTRA_COLSPEC)
    ]) + "\n")
    q = audit["query"]
    P.append(f"질의 단일성 — query 컬럼 `{q.get('column')}` · 서로 다른 값 "
             f"{q.get('n_distinct')} · {q.get('values', q.get('note', ''))}\n")

    P.append("## 2. 행 수와 누락 accession — 핵심 발견 지점\n")
    rows = [["정렬 파일 데이터 행", audit["n_rows"]],
            ["유일 accession", audit["n_unique_accession"]],
            ["기대 서열 수 (--expected-n)", args.expected_n or "검사 안 함"]]
    if "n_universe" in tr:
        rows += [["universe 서열 수", tr["n_universe"]],
                 ["정렬에 없는 accession", f"{tr['n_missing']} "
                                           f"({pct(tr['n_missing'], tr['n_universe'])})"],
                 ["정렬에만 있는 accession", tr["n_only_in_alignment"]],
                 ["조인율", fnum(tr.get("join_rate"), 4)]]
    if "shortfall" in tr:
        rows += [["기대치 대비 부족", tr["shortfall"]]]
    P.append(md_table(["항목", "값"], rows) + "\n")
    if tr.get("n_missing"):
        ex = tr["missing_examples"]
        P.append(f"누락 예시 ({len(ex)}개 표시"
                 f"{' · 상한에 걸려 잘림' if tr.get('missing_examples_truncated') else ''}, "
                 f"전체 {tr['n_missing']}개): `{', '.join(ex)}`\n")
        P.append(f"**해석.** mmseqs search 는 기본값이 `--max-seqs {MMSEQS_DEFAULT_MAX_SEQS}` 에 "
                 f"e-value {MMSEQS_DEFAULT_EVALUE} 다. 질의 하나당 보고되는 타깃이 "
                 f"{MMSEQS_DEFAULT_MAX_SEQS} 개로 잘리므로, 1,827 서열 전체가 들어오지 않는 것이 "
                 f"정상 동작이다. 명령이 기록되지 않았으므로 이 값이 쓰였는지는 알 수 없다.\n")
        P.append(f"**경고 — 같이 빠지는 것.** 빠진 {tr['n_missing']}개는 `fident` 가 없다. "
                 f"`fident → cf(클러스터 median) → cband` 사슬이 그대로 끊기므로 "
                 f"**cf 와 cband 도 같이 빠진다.** 멤버 일부만 빠진 클러스터는 cf 가 남은 멤버로만 "
                 f"계산되어 조용히 편향되고, 멤버가 전부 빠진 클러스터는 cf 가 아예 없다. "
                 f"단계 4 이후에서 cf 를 쓸 때는 반드시 분모(멤버 중 fident 가 있는 수)를 "
                 f"함께 보고할 것.\n")
    if tr.get("max_seqs_suspicion"):
        P.append(f"⚠ 행 수가 `--max-seqs` 기본값 {MMSEQS_DEFAULT_MAX_SEQS} 에 정확히 걸린 모양이다. "
                 f"search 를 `--max-seqs 2000` 이상으로 다시 돌릴 것.\n")

    P.append("## 3. 중복 HSP 병합\n")
    P.append(md_table(["항목", "값"], [
        ["HSP 가 여러 개인 accession", hsp["n_accession_with_multiple_hsp"]],
        ["그 accession 들의 총 행 수", hsp["n_rows_in_duplicate_accession"]],
        ["병합으로 버린 행", hsp["n_rows_to_merge"]],
        ["한 accession 최대 HSP 수", hsp["max_hsp_per_accession"]],
        ["선택 규칙", hsp["tiebreak_rule"]],
    ]) + "\n")
    if parse and parse.get("hsp_select"):
        hs = parse["hsp_select"]
        P.append(f"실제 사용한 정렬 키: `{', '.join(hs['tiebreak_keys_used'])}` · "
                 f"입력 {hs['n_rows_in']}행 → 선택 {hs['n_rows_selected']}행 "
                 f"(버림 {hs['n_rows_dropped']}행)\n")
        if hs["detail"]:
            P.append(md_table(["accession", "HSP 수", "채택 행", "버린 행", "bits"], [
                [d["accession"], d["n_hsp"], d["chosen_row"],
                 ",".join(str(x) for x in d["dropped_rows"]),
                 ",".join(fnum(v, 1) for v in d.get("bits", []))]
                for d in hs["detail"][:args.max_table_rows]
            ]) + "\n")
            if len(hs["detail"]) > args.max_table_rows:
                P.append(f"(표는 {args.max_table_rows}행만 표시 · 상세 전체 "
                         f"{hs['detail_listed']}건"
                         f"{' · params.json 상한에서 잘림' if hs['detail_truncated'] else ''})\n")
    if hsp["examples"]:
        P.append(f"예시: " + ", ".join(f"{e['accession']}×{e['n_hsp']}" for e in hsp["examples"])
                 + (f" (상한 표시)" if hsp["examples_truncated"] else "") + "\n")

    P.append("## 4. 수치 분포\n")
    rows = []
    for k, d in audit["distributions"].items():
        if not d.get("present"):
            rows.append([k, "— 컬럼 없음 —", "", "", "", "", "", ""])
            continue
        rows.append([k, d["n_valid"], d.get("n_blank"), d.get("n_unparseable"),
                     fnum(d.get("min"), 4), fnum(d.get("median"), 4),
                     fnum(d.get("max"), 4), fnum(d.get("mean"), 4)])
    P.append(md_table(["컬럼", "유효", "빈값", "파싱불가", "min", "median", "max", "mean"], rows) + "\n")
    fid = audit["distributions"]["fident"]
    if fid.get("present"):
        P.append(f"fident = 1.0 인 행 {fid.get('n_eq_one')}개 · 1.0 초과 {fid.get('n_gt_one')}개. "
                 f"초과가 있으면 이 컬럼은 분수가 아니라 퍼센트(pident)다.\n")
    ev = audit["distributions"]["evalue"]
    if ev.get("present") and ev.get("n_valid"):
        P.append(f"e-value 최대 {ev['max']:g} — 기본 절단값 {MMSEQS_DEFAULT_EVALUE} "
                 f"{'이하 (기본 search 로 보인다)' if ev.get('at_default_evalue_cutoff') else '초과 (-e 를 올려 돌렸다)'}\n")

    P.append("## 5. 자기 자신 히트\n")
    sh = audit["self_hit"]
    P.append(md_table(["항목", "값"], [
        ["참조 accession", sh["reference_acc"]],
        ["행 수", sh["n_rows"]],
        ["fident 값", sh.get("fident_values")],
        ["fident == 1.0", sh.get("fident_is_one")],
        ["qcov / tcov", f"{sh.get('qcov_values')} / {sh.get('tcov_values')}"],
    ]) + "\n")
    if sh["n_rows"] == 0:
        P.append("자기 히트가 없다. **검증 3(q_pos == t_pos 636개)을 할 수 없다** — "
                 "좌표계의 가장 강한 그물이 하나 빠진 상태다.\n")

    if parse is None:
        P.append("## 6. 파싱 — 미수행\n")
        P.append("이 판본은 파싱 전에 쓰였다.\n")
    else:
        base = parse["base"]
        P.append("## 6. 좌표 기준(base) 판정과 검증\n")
        P.append(f"**가정.** mmseqs 의 `qstart` 는 1-based 다. 이 스크립트는 그 가정을 믿지 않고 "
                 f"두 해석을 모두 시험한다.\n")
        P.append(md_table(["해석", "qaln↔참조 일치 행", "대소문자만 다름", "범위 초과"], [
            ["1-based (qstart-1)", base["n_match_1based"], base["n_case_only_1based"],
             base["n_out_of_range_1based"]],
            ["0-based (qstart)", base["n_match_0based"], base["n_case_only_0based"],
             base["n_out_of_range_0based"]],
        ]) + "\n")
        P.append(f"**채택: {base['base']}-based.** 근거 — {base['basis']}. "
                 f"검사 행 {base['n_checked']}개"
                 + (" · ⚠ 두 해석이 동률이라 근거가 약하다" if base.get("ambiguous") else "") + "\n")
        P.append(f"색인 0 의 직접 증거 — qstart 최소 {base.get('qstart_min')} · tstart 최소 "
                 f"{base.get('tstart_min')} · qstart==0 {base.get('n_qstart_zero')}행 · "
                 f"tstart==0 {base.get('n_tstart_zero')}행. "
                 f"0 이 하나도 없는데 0-based 가 이기면 base 문제가 아니라 **참조서열과 질의서열이 "
                 f"1잔기 어긋난 것**일 가능성이 높고, 그 경우 t_pos 가 전부 +1 밀린다 "
                 f"(자동 판정은 이 조합을 예외로 거절한다).\n")
        ql = parse.get("query_length", {})
        if ql.get("qlen_present"):
            P.append(f"질의 길이(qlen) 대조 — 컬럼 `{ql.get('qlen_column')}` · 값 "
                     f"{ql.get('qlen_values')} · 참조 {ql.get('reference_len')} aa · "
                     f"불일치 {ql.get('n_qlen_mismatch')}행 (불일치가 있으면 예외다). "
                     f"인계문서 2절의 '질의 서열 출처 미확인' 을 이 컬럼이 직접 메워 준다.\n")
        else:
            P.append("질의 길이(qlen) 컬럼이 없다 — search 질의가 정말 참조 636 aa 였는지 "
                     "직접 확인할 수 없고 검증 1·3 의 간접 증거에만 의존한다.\n")
        P.append(md_table(["검증", "내용", "결과"], [
            ["검증 2", "len(qaln) == len(taln)",
             f"{parse['lengths']['n_checked']}행 전부 통과 "
             f"(정렬 열 수 {parse['lengths']['aln_len_min']}~{parse['lengths']['aln_len_max']}, "
             f"중앙값 {fnum(parse['lengths']['aln_len_median'], 1)})"],
            ["검증 1", "qaln 갭 제거 == 참조[qstart..] 부분서열",
             f"{parse['n_mapped']}행 전부 정확히 일치 (base={base['base']})"],
            ["검증 3", f"자기 히트에서 q_pos == t_pos 가 {parse['reference_seq_len']}개",
             ("미수행 — 자기 히트 없음" if not parse["self_hit"].get("present")
              else f"통과 ({parse['self_hit']['n_q_eq_t']}/{parse['reference_seq_len']})"
              if parse["self_hit"].get("passed")
              else f"**실패했으나 --allow-partial-self-hit 로 진행** "
                   f"({parse['self_hit'].get('n_q_eq_t')}/{parse['reference_seq_len']})")],
        ]) + "\n")
        em = parse["end_mismatch"]
        P.append(f"보조 확인 — qend 불일치 {em['qend']}행 · tend 불일치 {em['tend']}행 · "
                 f"tlen 초과 {em['tlen_exceeded']}행 (해당 컬럼이 있을 때만 센다). "
                 f"0 이 아니면 base 판정이나 정렬 문자열을 다시 볼 것.\n")

        P.append("## 7. 매핑 요약\n")
        P.append(md_table(["항목", "값"], [
            ["행렬 행 수 (= universe)", parse["n_matrix_rows"]],
            ["정렬이 있는 accession", parse["n_in_alignment"]],
            ["정렬이 없어 전부 갭인 행", parse["n_not_in_alignment"]],
            ["등가 쌍 (equiv_map 행 수)", parse["n_equiv_rows"]],
            ["결손 열 (q 잔기, t 갭)", parse["deletion_cols"]],
            ["삽입 열 (q 갭, t 잔기)", parse["insertion_cols"]],
            ["위치당 잔기 수 min/median/max",
             f"{parse['coverage']['min_n_residue']} / "
             f"{fnum(parse['coverage']['median_n_residue'], 1)} / "
             f"{parse['coverage']['max_n_residue']}"],
            ["잔기가 하나도 없는 위치", parse["coverage"]["n_positions_zero_residue"]],
            ["16개 보고 위치의 잔기 수 min/median",
             f"{parse['coverage']['known_min_n_residue']} / "
             f"{fnum(parse['coverage']['known_median_n_residue'], 1)}"],
        ]) + "\n")
        P.append(f"**삽입은 표현되지 않는다.** 상동체에만 있는 잔기 {parse['insertion_residues']}개는 "
                 f"ChCODH2 번호를 가질 수 없어 equiv_map 과 residue_matrix 양쪽에서 사라진다. "
                 f"참조 좌표계 투영의 본질적 한계이며, 보존도 계산에는 영향이 없지만 "
                 f"'상동체에 추가 루프가 있다' 는 관찰은 이 산출물로 할 수 없다.\n")
        if parse["nonstandard_target_chars"]:
            P.append(f"상동체 쪽 비표준 문자: "
                     f"{json.dumps(parse['nonstandard_target_chars'], ensure_ascii=False)} "
                     f"(X 등 모호 문자는 보존도에서 분모에 들어가되 일치로 세지 않는다)\n")

    P.append("## 8. 산출물\n")
    if outputs:
        P.append(md_table(["파일", "행", "설명"], [
            [f"`{Path(outputs['equiv_map']['path']).name}`", outputs["equiv_map"]["n_rows"],
             "long. 한 행 = 잔기 대응 한 쌍 (q_pos 는 ChCODH2 번호 1~636)"],
            [f"`{Path(outputs['residue_matrix']['path']).name}`",
             outputs["residue_matrix"]["n_rows"],
             f"wide. 열 = 1..{outputs['residue_matrix']['n_position_columns']}. 보존도의 직접 입력"],
            [f"`{Path(outputs['known_positions']['path']).name}`",
             outputs["known_positions"]["n_rows"],
             f"보고 위치 {len(outputs['known_positions']['positions'])}개 "
             f"({', '.join(outputs['known_positions']['columns'][:4])}, …) + in_alignment"],
            [f"`{Path(outputs['coverage']['path']).name}`", outputs["coverage"]["n_rows"],
             "위치별 잔기/갭/미커버 수 — 보존도 신뢰도의 분모"],
            ["`alignment_audit.md`", "—", "이 문서"],
            ["`params.json`", "—", "입력 해시 · 행 수 · 판정 근거 · 상한 전부"],
        ]) + "\n")
    else:
        P.append("아직 쓰이지 않았다 (파싱 전 판본).\n")

    P.append("## 9. 불일치 목록\n")
    if discrepancies:
        P.append(f"{len(discrepancies)}건. 예외가 아니라 **보고**이며, `--strict` 면 종료코드 2 다.\n")
        for i, d in enumerate(discrepancies, 1):
            P.append(f"{i}. {d}")
        P.append("")
    else:
        P.append("없음.\n")

    path.write_text("\n".join(P), encoding="utf-8")


# ==========================================================================
# 본체
# ==========================================================================
def run(args) -> dict:
    alignment = Path(args.alignment)
    master = Path(args.master)
    fasta = Path(args.fasta)
    outdir = Path(args.outdir)

    # 행렬은 한 칸에 한 글자라는 계약이다. 두 글자 이상이면 numpy '<U1' 이 조용히 잘라내고,
    # 빈 문자열이면 TSV 에 빈 칸이 생겨 하위 단계에서 NaN 이 된다. 미리 거절한다.
    if len(args.outside_char) != 1 or args.outside_char in ("\t", "\n", "\r"):
        raise ValueError(
            f"--outside-char 는 탭/줄바꿈이 아닌 **한 글자**여야 한다 (받은 값 "
            f"{args.outside_char!r}, {len(args.outside_char)}글자). 두 글자 이상이면 행렬이 "
            f"조용히 첫 글자로 잘리고, 빈 문자열이면 TSV 의 빈 칸이 하위 단계에서 NaN 이 된다.")

    if args.dry_run:
        print("=== --dry-run : 경로와 규약만 확인하고 아무것도 쓰지 않는다 ===")
        for label, p in [("정렬 TSV", alignment), ("master", master), ("IPR FASTA", fasta),
                         ("참조 FASTA", common.REPO_LIB / "chcodh2.fasta")]:
            print(f"  {label:10s} {'있음' if p.exists() else '없음'}  {p}")
        for name in ("alignment_audit.md", "equiv_map.tsv", "residue_matrix.tsv",
                     "known_position_equivalents.tsv", "position_coverage.tsv", "params.json"):
            print(f"  출력 예정 {outdir / name}")
        common.verify_known_positions()
        print(f"  참조 {common.CHCODH2_ACC} {common.CHCODH2_LEN} aa · 보고 위치 "
              f"{len(common.KNOWN_POSITIONS)}개 번호 체계 점검 통과")
        print(f"  좌표 규약 q_pos = ChCODH2 번호 (Met1 = 1, 1~{common.CHCODH2_LEN}) · "
              f"t_pos = 상동체 1-based")
        print(f"  base 판정 {args.assume_qstart_base} · HSP 선택 {HSP_TIEBREAK_DESC}")
        print(f"  기대 서열 수 {args.expected_n or '검사 안 함'} · 갭 '{DEFAULT_GAP_CHAR}' · "
              f"미커버 '{args.outside_char}'")
        return dict(dry_run=True)

    if not alignment.exists():
        raise FileNotFoundError(
            f"정렬 TSV 가 없다: {alignment} (--alignment 로 지정하거나 --dry-run / --self-test)")

    print("=== 단계 3 좌표계 구축 (equiv_map) ===")
    common.verify_known_positions()
    ref = common.chcodh2_seq()
    print(f"참조       {common.CHCODH2_ACC} {len(ref)} aa · 보고 위치 "
          f"{len(common.KNOWN_POSITIONS)}개 번호 체계 점검 통과")
    print(f"정렬 TSV   {alignment}")

    aln, read_meta = read_alignment(alignment, args.columns)
    print(f"[읽기] {read_meta['n_raw_rows']}행 × {read_meta['n_fields']}열 → 데이터 "
          f"{read_meta['n_rows']}행 · 컬럼 출처 {read_meta['column_source']}")
    colmap = resolve(aln)
    print(f"[컬럼] " + " · ".join(f"{k}={colmap[k]}" for k in REQUIRED_KEYS))
    print(f"        선택 " + " · ".join(
        f"{k}={colmap.get(k) or '없음'}" for k in list(OPTIONAL_KEYS) + list(EXTRA_COLSPEC)))

    universe, uni_meta = load_universe(master, fasta, max_examples=args.max_examples)
    if universe is not None:
        print(f"[universe] {uni_meta['source']} {len(universe)}개 "
              f"({uni_meta.get('source_path')})")
    else:
        print(f"[universe] 없음 — {uni_meta.get('error')}")

    audit = audit_alignment(aln, colmap, universe, args=args, read_meta=read_meta,
                            universe_meta=uni_meta)
    aln = audit.pop("dataframe")
    discrepancies = list(audit["discrepancies"])
    print(f"[감사] 데이터 {audit['n_rows']}행 · 유일 accession {audit['n_unique_accession']} · "
          f"기대 {args.expected_n or '검사 안 함'}")
    tr = audit["truncation"]
    if "n_universe" in tr:
        print(f"       universe {tr['n_universe']} · 정렬에 없는 accession {tr['n_missing']} "
              f"({pct(tr['n_missing'], tr['n_universe'])}) · 조인율 {fnum(tr.get('join_rate'), 4)} · "
              f"정렬에만 있는 것 {tr['n_only_in_alignment']}")
        if tr["n_missing"]:
            print(f"       누락 예시 {tr['missing_examples']}"
                  f"{' (상한에서 잘림)' if tr.get('missing_examples_truncated') else ''}")
    print(f"       HSP 중복 accession {audit['hsp']['n_accession_with_multiple_hsp']}개 · "
          f"병합으로 버릴 행 {audit['hsp']['n_rows_to_merge']}개")

    outdir.mkdir(parents=True, exist_ok=True)
    md_path = outdir / "alignment_audit.md"
    write_audit(md_path, audit, args=args, parse=None, outputs=None,
                discrepancies=discrepancies)
    print(f"[감사] 파싱 전 판본 기록 {md_path}")

    try:
        # 검증 2 는 HSP 선택보다 **먼저**, 전 행에 대고 한다. 선택 뒤에 하면 버려진 HSP 의
        # 깨진 정렬 문자열이 검사되지 않고 통과하며, bits 가 없을 때의 대체 선택키가
        # 바로 그 qaln 길이이므로 깨진 행이 선택될 수도 있다.
        lengths = check_lengths(aln, colmap, max_examples=args.max_examples)
        print(f"[검증 2] len(qaln) == len(taln) 전 행 {lengths['n_checked']}행 통과 "
              f"(HSP 선택 전)")

        qlen_info, qlen_disc = check_query_length(aln, colmap, ref,
                                                  max_examples=args.max_examples)
        discrepancies += qlen_disc
        if qlen_info.get("qlen_present"):
            print(f"[qlen] 질의 길이 {qlen_info['qlen_values']} == 참조 {len(ref)} 확인 "
                  f"(불일치 {qlen_info['n_qlen_mismatch']}행)")
        else:
            print(f"[qlen] 컬럼 없음 — 질의 서열 동일성을 직접 확인할 수 없다")

        sel, hsp_meta = select_hsps(aln, colmap, max_detail=args.max_detail)
        print(f"[HSP] {hsp_meta['n_rows_in']}행 → {hsp_meta['n_rows_selected']}행 "
              f"(버림 {hsp_meta['n_rows_dropped']}행) · 키 "
              f"{', '.join(hsp_meta['tiebreak_keys_used'])}")

        base = detect_qstart_base(sel, colmap, ref, assume=args.assume_qstart_base,
                                  max_examples=args.max_examples, upper_aln=args.upper_aln)
        print(f"[검증 1] base 판정 — 1-based {base['n_match_1based']}/{base['n_checked']} 일치 · "
              f"0-based {base['n_match_0based']}/{base['n_checked']} 일치")
        print(f"         채택 {base['base']}-based · {base['basis']}")
        print(f"         색인 0 증거 — qstart 최소 {base['qstart_min']} · tstart 최소 "
              f"{base['tstart_min']} · qstart==0 {base['n_qstart_zero']}행 · "
              f"tstart==0 {base['n_tstart_zero']}행")
        if base["base"] == 0 and not base["zero_index_present"]:
            # --assume-qstart-base 0 으로 강제한 경로. 자동 판정은 이 조합을 거부한다.
            discrepancies.append(
                f"0-based 를 강제 채택했는데 파일에 색인 0 이 하나도 없다 "
                f"(qstart 최소 {base['qstart_min']} · tstart 최소 {base['tstart_min']}). "
                f"**t_pos 가 전부 +1 밀려 있을 수 있다** — 참조서열과 질의서열이 1잔기 "
                f"어긋난 경우가 0-based 로 위장하는 전형적인 모양이다. qlen 과 자기 히트로 "
                f"반드시 확인할 것.")
        if base.get("ambiguous"):
            discrepancies.append(
                f"base 판정이 동률이다 (1-based {base['n_match_1based']} = 0-based "
                f"{base['n_match_0based']}). 1-based 를 채택했으나 근거가 약하다.")

        aligned = sel["accession"].tolist()
        if args.matrix_aligned_only or universe is None:
            rows = aligned
            row_source = "정렬에 있는 accession" + (" (--matrix-aligned-only)"
                                                   if args.matrix_aligned_only else "")
        else:
            extra = [a for a in aligned if a not in set(universe)]
            rows = list(universe) + extra
            row_source = f"universe({uni_meta['source']}) + 정렬에만 있는 {len(extra)}개"
        if len(rows) != len(set(rows)):
            raise AssertionError("행렬 행 accession 에 중복이 있다 — universe 정규화 결함")
        row_index = {a: i for i, a in enumerate(rows)}
        print(f"[행렬] 행 {len(rows)}개 ({row_source})")

        mapping = build_mapping(sel, colmap, ref, base["base"], rows=rows,
                                row_index=row_index, args=args)
        discrepancies += mapping["discrepancies"]
        n_in = int(mapping["in_aln"].sum())
        print(f"[검증 1] qaln↔참조 부분서열 일치 {n_in}행 전부 통과 (base={base['base']})")

        self_info, self_disc = check_self_hit(
            mapping, ref, row_index=row_index, ref_acc=safe_parse_accession(args.reference_acc),
            allow_partial=args.allow_partial_self_hit)
        discrepancies += self_disc
        if self_info.get("present"):
            print(f"[검증 3] 자기 히트 등가 쌍 {self_info['n_equiv_pairs']}/{len(ref)} · "
                  f"q_pos==t_pos {self_info['n_q_eq_t']}/{len(ref)} · "
                  f"행렬행==참조서열 {self_info['matrix_row_equals_reference']} → "
                  f"{'통과' if self_info['passed'] else '실패(완화됨)'}")
        else:
            print(f"[검증 3] 자기 히트 없음 — 생략")

        em = mapping["end_mismatch"]
        for k, v in em.items():
            if v:
                discrepancies.append(
                    f"{k} 보조 확인 불일치 {v}행 — base 판정({base['base']}-based) 또는 정렬 "
                    f"문자열을 재확인할 것 (매핑 자체는 qaln/qstart 만 쓰므로 영향 없음).")

        cov_df, cov_meta = build_coverage(mapping["mat"], mapping["cov"],
                                          mapping["in_aln"], ref)

        outputs = dict(
            equiv_map=write_equiv_map(outdir / "equiv_map.tsv", mapping["equiv"]),
            residue_matrix=write_residue_matrix(outdir / "residue_matrix.tsv",
                                                mapping["mat"], rows),
            known_positions=write_known_positions(outdir / "known_position_equivalents.tsv",
                                                  mapping["mat"], mapping["cov"],
                                                  mapping["in_aln"], rows, ref),
            coverage=dict(path=str(outdir / "position_coverage.tsv"), **cov_meta),
        )
        cov_df.to_csv(outdir / "position_coverage.tsv", sep="\t", index=False)

        n_equiv = outputs["equiv_map"]["n_rows"]
        expect_equiv = sum(p["n_equiv_pairs"] for p in mapping["per_acc"])
        if n_equiv != expect_equiv:
            raise AssertionError(f"equiv_map 행 수 {n_equiv} != 매핑 합계 {expect_equiv}")
        if outputs["residue_matrix"]["n_rows"] != len(rows):
            raise AssertionError("residue_matrix 행 수가 변했다")
        if outputs["known_positions"]["n_rows"] != len(rows):
            raise AssertionError("known_position_equivalents 행 수가 변했다")

        print(f"[산출] equiv_map {n_equiv}행 · residue_matrix {len(rows)}×{len(ref)} · "
              f"known_positions {len(rows)}×{len(outputs['known_positions']['columns'])} · "
              f"position_coverage {cov_meta['n_rows']}행")
        print(f"       정렬 있는 행 {n_in} · 전부 갭인 행 {len(rows) - n_in} · "
              f"결손 열 {mapping['deletion_cols']} · 삽입 열 {mapping['insertion_cols']} "
              f"(삽입 잔기는 ChCODH2 번호가 없어 산출물에 표현되지 않는다)")
        print(f"       위치당 잔기 수 min {cov_meta['min_n_residue']} / median "
              f"{fnum(cov_meta['median_n_residue'], 1)} / max {cov_meta['max_n_residue']} · "
              f"잔기 0인 위치 {cov_meta['n_positions_zero_residue']}개")

        parse = dict(
            hsp_select=hsp_meta, lengths=lengths, query_length=qlen_info,
            base=base, self_hit=self_info, reference_seq_len=len(ref),
            n_matrix_rows=len(rows), matrix_row_source=row_source,
            n_in_alignment=n_in, n_not_in_alignment=len(rows) - n_in,
            n_mapped=n_in, n_equiv_rows=n_equiv,
            deletion_cols=mapping["deletion_cols"],
            insertion_cols=mapping["insertion_cols"],
            insertion_residues=mapping["insertion_residues"],
            end_mismatch=em, nonstandard_target_chars=mapping["nonstandard_target_chars"],
            coverage=cov_meta,
            per_accession=mapping["per_acc"][:args.max_detail],
            per_accession_listed=min(len(mapping["per_acc"]), args.max_detail),
            per_accession_truncated=bool(len(mapping["per_acc"]) > args.max_detail),
        )
        write_audit(md_path, audit, args=args, parse=parse, outputs=outputs,
                    discrepancies=discrepancies)

        stats = dict(
            status="ok",
            input_paths=dict(alignment=str(alignment), master=str(master), fasta=str(fasta),
                             reference_fasta=str(common.REPO_LIB / "chcodh2.fasta")),
            output_paths={k: v["path"] for k, v in outputs.items()} | dict(audit=str(md_path)),
            coordinate_convention=dict(
                q_pos="ChCODH2 1-based (Met1 = 1, 1..636)",
                t_pos="상동체 1-based (입력 base 를 1-based 로 환산)",
                gap_char=DEFAULT_GAP_CHAR, outside_char=args.outside_char,
                qstart_base_adopted=base["base"], qstart_base_option=args.assume_qstart_base),
            thresholds=dict(expected_n=args.expected_n, max_examples=args.max_examples,
                            max_table_rows=args.max_table_rows, max_detail=args.max_detail,
                            matrix_aligned_only=bool(args.matrix_aligned_only),
                            allow_partial_self_hit=bool(args.allow_partial_self_hit),
                            upper_aln=bool(args.upper_aln),
                            query_acc=args.query_acc),
            audit={k: v for k, v in audit.items() if k != "discrepancies"},
            parse=parse,
            output_tables=outputs,
            discrepancies=discrepancies,
            n_discrepancies=len(discrepancies),
        )
        params = common.write_params(
            outdir, STAGE,
            inputs=[alignment] + [p for p in (master, fasta) if p.exists()]
                   + [common.REPO_LIB / "chcodh2.fasta"],
            outputs=[outdir / n for n in ("equiv_map.tsv", "residue_matrix.tsv",
                                          "known_position_equivalents.tsv",
                                          "position_coverage.tsv", "alignment_audit.md")],
            commands=[" ".join(sys.argv)], **stats)
        print(f"출력  " + "\n      ".join(
            [outputs[k]["path"] for k in ("equiv_map", "residue_matrix", "known_positions",
                                          "coverage")] + [str(md_path), str(params)]))
        if discrepancies:
            print(f"[요약] 불일치 {len(discrepancies)}건 — 감사 md 9절과 params.json 참조")
            for d in discrepancies:
                print(f"  - {d}")
        return dict(stats, coverage_df=cov_df, mapping=mapping, rows=rows,
                    row_index=row_index, selected=sel, outputs=outputs,
                    discrepancies=discrepancies)

    except Exception as exc:  # noqa: BLE001
        # 규칙 5 — 실패해도 산출 디렉터리에 params.json 을 남긴다. 실패가 기록되지 않으면
        # 다음 사람이 「돌렸는데 왜 산출물이 없나」를 처음부터 다시 조사해야 한다.
        common.write_params(
            outdir, STAGE,
            inputs=[alignment], outputs=[md_path], commands=[" ".join(sys.argv)],
            status="failed", error_type=type(exc).__name__, error=str(exc),
            traceback=traceback.format_exc()[-6000:],
            audit={k: v for k, v in audit.items() if k != "discrepancies"},
            discrepancies=discrepancies, n_discrepancies=len(discrepancies))
        print(f"[실패] {type(exc).__name__} — params.json 에 status=failed 로 기록했다: "
              f"{outdir / 'params.json'}", file=sys.stderr)
        raise


# ==========================================================================
# 합성 데이터 자기검증
# ==========================================================================
ALN_FIELDS = ["query", "target", "fident", "alnlen", "qstart", "qend", "tstart", "tend",
              "qaln", "taln", "evalue", "bits", "qcov", "tcov"]


def _mutate(seq: str, rnd: random.Random, rate: float = 0.25) -> str:
    alpha = "ACDEFGHIKLMNPQRSTVWY"
    out = []
    for ch in seq:
        if rnd.random() < rate:
            c = rnd.choice(alpha)
            while c == ch:
                c = rnd.choice(alpha)
            out.append(c)
        else:
            out.append(ch)
    return "".join(out)


def _row(acc, qaln, taln, qstart, tstart, *, fident, bits, evalue, qcov, tcov,
         base: int = 1) -> list[str]:
    """ALN_FIELDS 순서의 한 행. base=0 이면 qstart/tstart/qend/tend 를 0-based 로 적는다."""
    nq = len(qaln.replace("-", ""))
    nt = len(taln.replace("-", ""))
    off = 0 if base == 1 else -1
    return [common.CHCODH2_ACC, acc, f"{fident:.3f}", str(len(qaln)),
            str(qstart + off), str(qstart + nq - 1 + off),
            str(tstart + off), str(tstart + nt - 1 + off),
            qaln, taln, f"{evalue:g}", f"{bits:.1f}", f"{qcov:.3f}", f"{tcov:.3f}"]


def naive_map(ref: str, qstart: int, qaln: str, tstart: int, taln: str) -> list[tuple]:
    """검증용 **독립 구현**. 벡터화 없이 한 글자씩 세어 등가 쌍을 만든다.

    build_mapping 의 numpy 누적합 구현과 결과가 같아야 한다. 같은 버그를 두 번 쓰지
    않도록 일부러 다른 방식으로 적었다.
    """
    out = []
    q = qstart
    t = tstart
    for cq, ct in zip(qaln, taln):
        if cq != "-" and ct != "-":
            out.append((q, cq, t, ct))
        if cq != "-":
            q += 1
        if ct != "-":
            t += 1
    return out


def make_synthetic(dirpath: Path) -> dict:
    """합성 정렬 TSV 들을 만든다. 결정적(seed 고정)."""
    dirpath.mkdir(parents=True, exist_ok=True)
    rnd = random.Random(20261022)
    ref = common.chcodh2_seq()
    L = len(ref)

    # --- (1) 정상 정렬: 전장, 갭 없음. 보고 위치에 특정 잔기를 심어 둔다.
    t1 = list(_mutate(ref, rnd, 0.20))
    t1[40] = "Y"      # pos41 F -> Y
    t1[558] = "W"     # pos559 A -> W
    t1[564] = ref[564]
    syn1 = _row("SYN00001", ref, "".join(t1), 1, 1,
                fident=0.80, bits=900.0, evalue=1e-180, qcov=1.0, tcov=1.0)

    # --- (2) 양쪽 갭: q 10..309 (300 잔기). 결손 3 + 삽입 2.
    seg = ref[9:309]
    mt = list(_mutate(seg, rnd, 0.25))
    mt[31] = "H"      # q=41 (열 31) -> H
    qaln2 = seg[:53] + "--" + seg[53:]
    taln2 = "".join(mt[:50]) + "---" + "KK" + "".join(mt[53:])
    syn2 = _row("SYN00002", qaln2, taln2, 10, 5,
                fident=0.55, bits=320.0, evalue=1e-80, qcov=0.47, tcov=0.52)

    # --- (3) 자기 자신 히트: 완전 일치
    selfrow = _row(common.CHCODH2_ACC, ref, ref, 1, 1,
                   fident=1.0, bits=1300.0, evalue=0.0, qcov=1.0, tcov=1.0)

    # --- (4) 같은 accession 의 HSP 2건. bits 가 큰 쪽(A)이 채택되어야 한다.
    hsp_a = _row("SYN00003", ref[99:199], _mutate(ref[99:199], rnd, 0.3), 100, 1,
                 fident=0.70, bits=250.0, evalue=1e-40, qcov=0.157, tcov=0.90)
    hsp_b = _row("SYN00003", ref[299:349], _mutate(ref[299:349], rnd, 0.3), 300, 400,
                 fident=0.60, bits=80.0, evalue=1e-5, qcov=0.079, tcov=0.10)

    # --- (5) 오른쪽 끝 국소 정렬: 터널 위치만 커버 (41 은 미커버)
    t5 = list(_mutate(ref[549:620], rnd, 0.3))
    t5[9] = "W"       # q=559 -> W
    t5[60] = "H"      # q=610 -> H
    syn4 = _row("SYN00004", ref[549:620], "".join(t5), 550, 10,
                fident=0.68, bits=140.0, evalue=1e-25, qcov=0.112, tcov=0.20)

    def dump(path: Path, rows, header=True):
        with open(path, "wt", encoding="utf-8") as fh:
            if header:
                fh.write("\t".join(ALN_FIELDS) + "\n")
            for r in rows:
                fh.write("\t".join(r) + "\n")

    good_rows = [syn1, syn2, selfrow, hsp_a, hsp_b, syn4]
    ok_path = dirpath / "aln_good.tsv"
    dump(ok_path, good_rows)
    noheader = dirpath / "aln_good_noheader.tsv"
    dump(noheader, good_rows, header=False)

    # --- 0-based 판본: 같은 정렬을 qstart/tstart 0-based 로 적는다.
    zero_rows = [
        _row("SYN00001", ref, "".join(t1), 1, 1, fident=0.80, bits=900.0, evalue=1e-180,
             qcov=1.0, tcov=1.0, base=0),
        _row("SYN00002", qaln2, taln2, 10, 5, fident=0.55, bits=320.0, evalue=1e-80,
             qcov=0.47, tcov=0.52, base=0),
        _row(common.CHCODH2_ACC, ref, ref, 1, 1, fident=1.0, bits=1300.0, evalue=0.0,
             qcov=1.0, tcov=1.0, base=0),
        _row("SYN00003", ref[99:199], _mutate(ref[99:199], rnd, 0.3), 100, 1,
             fident=0.70, bits=250.0, evalue=1e-40, qcov=0.157, tcov=0.90, base=0),
        _row("SYN00004", ref[549:620], "".join(t5), 550, 10, fident=0.68, bits=140.0,
             evalue=1e-25, qcov=0.112, tcov=0.20, base=0),
    ]
    zero_path = dirpath / "aln_zerobased.tsv"
    dump(zero_path, zero_rows)

    # --- 실패 케이스 ① qaln 이 참조와 어긋남 (qstart=200 인데 내용은 ref[0:80])
    bad_q = dirpath / "aln_bad_qaln.tsv"
    dump(bad_q, [syn1, selfrow,
                 _row("SYNBAD1", ref[0:80], _mutate(ref[0:80], rnd, 0.2), 200, 1,
                      fident=0.75, bits=200.0, evalue=1e-30, qcov=0.126, tcov=0.5)])

    # --- 실패 케이스 ② len(qaln) != len(taln)
    bad_len_row = _row("SYNBAD2", ref[0:100], _mutate(ref[0:100], rnd, 0.2), 1, 1,
                       fident=0.8, bits=200.0, evalue=1e-30, qcov=0.157, tcov=0.5)
    bad_len_row[9] = bad_len_row[9][:-7]          # taln 만 7자 짧게
    bad_len = dirpath / "aln_bad_len.tsv"
    dump(bad_len, [syn1, selfrow, bad_len_row])

    # --- 실패 케이스 ③ 자기 히트가 부분 정렬 (검증 3)
    bad_self = dirpath / "aln_bad_self.tsv"
    dump(bad_self, [syn1,
                    _row(common.CHCODH2_ACC, ref[0:600], ref[0:600], 1, 1,
                         fident=1.0, bits=1200.0, evalue=0.0, qcov=0.94, tcov=0.94)])

    # --- 실패 케이스 ④ 질의가 여러 종 섞임
    multi = dirpath / "aln_multi_query.tsv"
    other = list(syn1)
    other[0] = "P99999"
    dump(multi, [syn1, selfrow, other])

    # --- 실패 케이스 ⑤ 질의서열이 참조와 1잔기 어긋남 (개시 Met 가 빠진 질의로 search).
    #     0-based 로 위장해 q_pos 는 맞게 나오지만 t_pos 가 전부 +1 밀린다. 자기 히트도 없다.
    #     색인 0 이 하나도 없다는 점이 유일한 단서이고, 그것으로 거절해야 한다.
    shifted_q = ref[1:]                      # 635 aa
    shift_rows = []
    for i, (s, n) in enumerate([(10, 300), (200, 400)]):
        qa = shifted_q[s - 1:s - 1 + n]
        shift_rows.append(_row(f"SYNSHIFT{i}", qa, _mutate(qa, rnd, 0.25), s, 7 + i,
                               fident=0.7, bits=600.0, evalue=1e-99, qcov=0.9, tcov=0.9))
    shift = dirpath / "aln_shifted_query.tsv"
    dump(shift, shift_rows)

    # --- 실패 케이스 ⑥ qlen 이 참조 길이와 다름 (좌표는 1-based 로 정상)
    qlen_fields = ALN_FIELDS + ["qlen"]
    qlen_bad = dirpath / "aln_bad_qlen.tsv"
    with open(qlen_bad, "wt", encoding="utf-8") as fh:
        fh.write("\t".join(qlen_fields) + "\n")
        for r in (syn1, selfrow):
            fh.write("\t".join(list(r) + [str(L - 1)]) + "\n")

    # --- 실패 케이스 ⑦ taln 만 소문자 (소프트 마스킹). 검증 1~3 은 통과하지만
    #     행렬에 소문자가 들어가 n_identical 이 조용히 과소계산된다.
    lower_tgt = dirpath / "aln_lower_taln.tsv"
    low_row = _row("SYNLOW1", ref, "".join(t1).lower(), 1, 1,
                   fident=0.80, bits=890.0, evalue=1e-170, qcov=1.0, tcov=1.0)
    dump(lower_tgt, [selfrow, low_row])

    # --- 실패 케이스 ⑧ 헤더에 같은 컬럼명이 두 번
    dupcol = dirpath / "aln_dup_column.tsv"
    with open(dupcol, "wt", encoding="utf-8") as fh:
        f2 = list(ALN_FIELDS)
        f2[3] = "qaln"                       # alnlen -> qaln (중복)
        fh.write("\t".join(f2) + "\n")
        for r in (syn1, selfrow):
            fh.write("\t".join(r) + "\n")

    # --- master 에 정규화 후 중복 accession 이 있는 판본 (universe 가 조용히 줄어드는 경로)
    master_dup = dirpath / "synthetic_master_dup.tsv"
    with open(master_dup, "wt", encoding="utf-8") as fh:
        fh.write("accession\tcluster\tis_mm_rep\n")
        for a in ["tr|SYN00001|S1_ORGA", "sp|SYN00001|S1_ALT", "tr|SYN00002|S2_ORGB",
                  f"tr|{common.CHCODH2_ACC}|CODH2_CARHZ"]:
            fh.write(f"{a}\tc1\t0\n")

    # --- universe (master 대역). SYN00005/SYN00006 은 일부러 정렬에 없다.
    master = dirpath / "synthetic_master.tsv"
    uni = ["tr|SYN00001|S1_ORGA", "tr|SYN00002|S2_ORGB", "tr|SYN00003|S3_ORGC",
           "tr|SYN00004|S4_ORGD", "tr|SYN00005|S5_ORGE", "tr|SYN00006|S6_ORGF",
           f"tr|{common.CHCODH2_ACC}|CODH2_CARHZ"]
    with open(master, "wt", encoding="utf-8") as fh:
        fh.write("accession\tcluster\tis_mm_rep\n")
        for a in uni:
            fh.write(f"{a}\t{uni[0]}\t0\n")

    return dict(
        dir=dirpath, ref=ref, ref_len=L, master=master,
        good=ok_path, noheader=noheader, zerobased=zero_path,
        bad_qaln=bad_q, bad_len=bad_len, bad_self=bad_self, multi_query=multi,
        shifted=shift, bad_qlen=qlen_bad, lower_taln=lower_tgt, dup_column=dupcol,
        master_dup=master_dup,
        n_universe=len(uni), n_good_rows=len(good_rows),
        qaln2=qaln2, taln2=taln2, t1="".join(t1), t5="".join(t5),
        n_equiv_expected=L + 297 + L + 100 + 71,
    )


def _args_for(base_args, **over) -> argparse.Namespace:
    a = argparse.Namespace(**vars(base_args))
    a.self_test = None
    a.dry_run = False
    for k, v in over.items():
        setattr(a, k, v)
    return a


def expect_raises(label: str, fn) -> tuple[str, bool, str, str]:
    """예외가 나야 하는 케이스. (라벨, 발생여부, 한줄요약, 전문) 을 돌려준다.

    전문을 따로 돌려주는 이유 — 검증 1/2 의 예외 메시지는 여러 줄이고 핵심 증거(어느
    accession 이 어긋났는지)가 둘째 줄 이후에 있다. 한 줄만 보고 판정하면 그 증거를 놓친다.
    """
    try:
        fn()
    except Exception as exc:  # noqa: BLE001
        full = f"{type(exc).__name__}: {exc}"
        return label, True, f"{type(exc).__name__}: {str(exc).splitlines()[0][:200]}", full
    return label, False, "예외가 발생하지 않았다 (실패)", ""


def self_test(dirpath: Path, base_args) -> None:
    syn = make_synthetic(dirpath)
    ref = syn["ref"]
    print(f"=== --self-test : 합성 정렬 {syn['n_good_rows']}행 / universe "
          f"{syn['n_universe']}개 ===")
    for k in ("good", "noheader", "zerobased", "bad_qaln", "bad_len", "bad_self",
              "multi_query", "shifted", "bad_qlen", "lower_taln", "dup_column",
              "master", "master_dup"):
        print(f"  {k:12s} {syn[k]}")
    print()

    checks: list[tuple[str, bool, str]] = []

    def chk(name, ok, detail=""):
        checks.append((name, bool(ok), str(detail)))

    # ---------------- 정상 경로 ----------------
    out1 = dirpath / "out_good"
    args = _args_for(base_args, alignment=str(syn["good"]), master=str(syn["master"]),
                     fasta=str(dirpath / "does_not_exist.fasta"), outdir=str(out1),
                     expected_n=syn["n_universe"])
    res = run(args)

    eq = pd.read_csv(out1 / "equiv_map.tsv", sep="\t", dtype={"accession": str})
    mx = pd.read_csv(out1 / "residue_matrix.tsv", sep="\t", dtype=str,
                     keep_default_na=False, na_filter=False)
    kp = pd.read_csv(out1 / "known_position_equivalents.tsv", sep="\t", dtype=str,
                     keep_default_na=False, na_filter=False)
    cv = pd.read_csv(out1 / "position_coverage.tsv", sep="\t")
    params = json.loads((out1 / "params.json").read_text())

    chk("산출 5개 + params.json 생성",
        all((out1 / f).exists() for f in ("alignment_audit.md", "equiv_map.tsv",
                                          "residue_matrix.tsv",
                                          "known_position_equivalents.tsv",
                                          "position_coverage.tsv", "params.json")))
    chk("params status=ok", params["status"] == "ok", params["status"])

    # (A) 감사 — 누락 accession
    tr = res["audit"]["truncation"]
    chk("누락 accession 2개 (SYN00005/6) 탐지",
        tr["n_missing"] == 2 and set(tr["missing_examples"]) == {"SYN00005", "SYN00006"},
        str(tr.get("missing_examples")))
    chk("누락이 cf/cband 경고로 보고됨",
        any("cband" in d for d in res["discrepancies"]))
    chk("universe 출처 master", res["audit"]["universe"]["source"] == "master")
    chk("조인율 5/7", abs(tr["join_rate"] - 5 / 7) < 1e-9, fnum(tr["join_rate"], 4))

    # (A) 감사 — HSP 중복
    hs = res["parse"]["hsp_select"]
    chk("HSP 중복 1건 탐지 (SYN00003 ×2)",
        res["audit"]["hsp"]["n_accession_with_multiple_hsp"] == 1
        and res["audit"]["hsp"]["n_rows_to_merge"] == 1)
    chk("HSP 6행 → 5행 병합", hs["n_rows_in"] == 6 and hs["n_rows_selected"] == 5)
    chk("bits 최대 HSP 채택 (q 100..199)",
        sorted(eq.loc[eq.accession == "SYN00003", "q_pos"].agg(["min", "max"]).tolist())
        == [100, 199],
        str(eq.loc[eq.accession == "SYN00003", "q_pos"].agg(["min", "max"]).tolist()))
    chk("HSP 선택 키 기록", hs["tiebreak_keys_used"][0].startswith("bits"),
        str(hs["tiebreak_keys_used"]))

    # (A) 감사 — 자기 히트 / 분포
    chk("자기 히트 fident=1.0 확인", res["audit"]["self_hit"]["fident_is_one"] is True)
    chk("fident 분포 기록", res["audit"]["distributions"]["fident"]["n_valid"] == 6)
    chk("qcov/tcov 분포 기록",
        res["audit"]["distributions"]["qcov"]["present"]
        and res["audit"]["distributions"]["tcov"]["present"])

    # (B) 검증 1 — base 판정
    bs = res["parse"]["base"]
    chk("1-based 채택 + 전 행 일치",
        bs["base"] == 1 and bs["n_match_1based"] == 5 and bs["n_match_0based"] == 0,
        f"1based={bs['n_match_1based']} 0based={bs['n_match_0based']}")
    chk("검증 2 는 HSP 선택 전 전 행(6행)에 적용", res["parse"]["lengths"]["n_checked"] == 6,
        str(res["parse"]["lengths"]["n_checked"]))
    chk("qlen 컬럼 없음이 불일치로 보고됨", any("qlen 컬럼이 없다" in d for d in res["discrepancies"]))
    chk("색인 0 증거가 기록됨 (1-based 이므로 0 없음)",
        bs["n_qstart_zero"] == 0 and bs["n_tstart_zero"] == 0 and bs["qstart_min"] == 1,
        f"qmin={bs['qstart_min']} qzero={bs['n_qstart_zero']}")

    # (B) 검증 3 — 자기 히트 완전 일치
    sh = res["parse"]["self_hit"]
    chk(f"검증 3: q_pos==t_pos {len(ref)}개 전부",
        sh["passed"] and sh["n_q_eq_t"] == len(ref) and sh["n_equiv_pairs"] == len(ref)
        and sh["matrix_row_equals_reference"], json.dumps(sh, ensure_ascii=False))
    selfm = mx.loc[mx.accession == common.CHCODH2_ACC].iloc[0, 1:].tolist()
    chk("자기 히트 행렬행 == 참조서열", "".join(selfm) == ref)
    selfeq = eq[eq.accession == common.CHCODH2_ACC]
    chk("자기 히트 equiv_map 636행 · q_aa==t_aa",
        len(selfeq) == len(ref) and (selfeq.q_pos == selfeq.t_pos).all()
        and (selfeq.q_aa == selfeq.t_aa).all())

    # (B) 독립 구현과의 대조 — naive_map
    naive2 = naive_map(ref, 10, syn["qaln2"], 5, syn["taln2"])
    got2 = eq[eq.accession == "SYN00002"][["q_pos", "q_aa", "t_pos", "t_aa"]]
    chk("양쪽 갭 정렬이 독립 구현(naive_map)과 일치",
        len(naive2) == len(got2)
        and all(tuple(r) == n for r, n in zip(got2.itertuples(index=False), naive2)),
        f"naive {len(naive2)}행 vs 산출 {len(got2)}행")
    chk("양쪽 갭: 등가 쌍 297개 (q 잔기 300 − 결손 3)", len(got2) == 297, str(len(got2)))
    chk("양쪽 갭: 결손 3열 · 삽입 2열",
        res["parse"]["deletion_cols"] == 3 and res["parse"]["insertion_cols"] == 2,
        f"del={res['parse']['deletion_cols']} ins={res['parse']['insertion_cols']}")
    chk("양쪽 갭: q=41 -> t=36 'H'",
        got2[(got2.q_pos == 41)].iloc[0].tolist() == [41, "F", 36, "H"],
        str(got2[(got2.q_pos == 41)].iloc[0].tolist()))
    chk("양쪽 갭: 삽입 2잔기 뒤 q=63 -> t=57",
        got2[(got2.q_pos == 63)].iloc[0].tolist() == [63, ref[62], 57, got2[(got2.q_pos == 63)].iloc[0].t_aa],
        str(got2[(got2.q_pos == 63)].iloc[0].tolist()))
    chk("양쪽 갭: 결손 q=60,61,62 는 equiv_map 에 없다",
        got2[got2.q_pos.isin([60, 61, 62])].empty)
    chk("양쪽 갭: 결손 위치는 행렬에서 '-'",
        mx.loc[mx.accession == "SYN00002", ["60", "61", "62"]].iloc[0].tolist() == ["-"] * 3)

    # (C) equiv_map 총 행 수
    chk(f"equiv_map 총 {syn['n_equiv_expected']}행",
        len(eq) == syn["n_equiv_expected"], f"{len(eq)}행")
    chk("equiv_map q_pos 범위 1..636",
        int(eq.q_pos.min()) == 1 and int(eq.q_pos.max()) == len(ref))
    chk("equiv_map q_aa 가 전부 참조 잔기와 일치",
        all(ref[p - 1] == a for p, a in zip(eq.q_pos.tolist(), eq.q_aa.tolist())))

    # (C) residue_matrix
    chk(f"residue_matrix 행 = universe {syn['n_universe']}", len(mx) == syn["n_universe"],
        f"{len(mx)}행")
    chk("residue_matrix 열 = accession + 636", mx.shape[1] == len(ref) + 1, str(mx.shape))
    chk("정렬에 없는 SYN00005 는 전부 갭 (삭제되지 않음)",
        set(mx.loc[mx.accession == "SYN00005"].iloc[0, 1:].tolist()) == {"-"})
    chk("SYN00001 행렬행 == 심어둔 상동체 서열", "".join(
        mx.loc[mx.accession == "SYN00001"].iloc[0, 1:].tolist()) == syn["t1"])
    chk("SYN00004 는 550..620 만 잔기", (
        set(mx.loc[mx.accession == "SYN00004"].iloc[0, 1:550].tolist()) == {"-"}
        and "-" not in mx.loc[mx.accession == "SYN00004"].iloc[0, 550:621].tolist()
        and set(mx.loc[mx.accession == "SYN00004"].iloc[0, 621:].tolist()) == {"-"}))

    # (C) known_position_equivalents
    poscols = [c for c in kp.columns if c.startswith("pos")]
    chk(f"known 열 {len(common.KNOWN_POSITIONS)}개 · 이름에 참조 잔기 포함",
        len(poscols) == len(common.KNOWN_POSITIONS) and poscols[0] == "pos41_F"
        and "pos559_A" in poscols, str(poscols[:4]))
    chk("known 행 = universe", len(kp) == syn["n_universe"])
    k1 = kp.loc[kp.accession == "SYN00001"].iloc[0]
    chk("SYN00001 pos41_F = Y · pos559_A = W",
        k1["pos41_F"] == "Y" and k1["pos559_A"] == "W",
        f"{k1['pos41_F']}/{k1['pos559_A']}")
    k2 = kp.loc[kp.accession == "SYN00002"].iloc[0]
    chk("SYN00002 pos41_F = H · pos559_A = '-' (미커버)",
        k2["pos41_F"] == "H" and k2["pos559_A"] == "-")
    k4 = kp.loc[kp.accession == "SYN00004"].iloc[0]
    chk("SYN00004 pos559_A = W · pos610_V = H · pos41_F = '-'",
        k4["pos559_A"] == "W" and k4["pos610_V"] == "H" and k4["pos41_F"] == "-",
        f"{k4['pos559_A']}/{k4['pos610_V']}/{k4['pos41_F']}")
    kself = kp.loc[kp.accession == common.CHCODH2_ACC].iloc[0]
    chk("자기 히트의 known 열 = 참조 잔기",
        all(kself[f"pos{p}_{m['aa']}"] == m["aa"] for p, m in common.KNOWN_POSITIONS.items()))
    chk("in_alignment 플래그: 5 True / 2 False",
        (kp.in_alignment.astype(str).str.lower() == "true").sum() == 5)
    chk("정렬에 없는 행은 n_known_covered = 0",
        int(kp.loc[kp.accession == "SYN00006", "n_known_covered"].iloc[0]) == 0)

    # (C) position_coverage
    chk(f"position_coverage {len(ref)}행", len(cv) == len(ref), f"{len(cv)}행")
    chk("분해 항등식 (잔기+블록내갭+블록외+정렬없음 == universe)",
        bool(((cv.n_residue + cv.n_gap_within_aln + cv.n_outside_aln
               + cv.n_not_in_alignment) == cv.n_universe).all()))
    c41 = cv[cv.q_pos == 41].iloc[0]
    chk("pos41: 잔기 3 · 블록외 2 · 정렬없음 2 · 동일 1",
        (c41.n_residue, c41.n_gap_within_aln, c41.n_outside_aln,
         c41.n_not_in_alignment, c41.n_identical) == (3, 0, 2, 2, 1),
        str([c41.n_residue, c41.n_gap_within_aln, c41.n_outside_aln,
             c41.n_not_in_alignment, c41.n_identical]))
    c61 = cv[cv.q_pos == 61].iloc[0]
    chk("pos61(결손): 잔기 2 · 블록내갭 1 · 블록외 2 · 정렬없음 2",
        (c61.n_residue, c61.n_gap_within_aln, c61.n_outside_aln,
         c61.n_not_in_alignment) == (2, 1, 2, 2),
        str([c61.n_residue, c61.n_gap_within_aln, c61.n_outside_aln, c61.n_not_in_alignment]))
    c559 = cv[cv.q_pos == 559].iloc[0]
    chk("pos559: 잔기 3 (SYN00001/SYN00004/자기)",
        c559.n_residue == 3 and c559.is_known_position, str(c559.n_residue))
    chk("known_role 주석 존재",
        cv[cv.q_pos == 559].iloc[0].known_role == "tunnel_bottleneck")
    chk("q_aa 열이 참조서열과 완전히 일치", "".join(cv.q_aa.tolist()) == ref)

    # --- 헤더 없는 convertalis 산출물 + --columns
    out2 = dirpath / "out_noheader"
    res2 = run(_args_for(base_args, alignment=str(syn["noheader"]), master=str(syn["master"]),
                         fasta=str(dirpath / "nope.fasta"), outdir=str(out2),
                         expected_n=syn["n_universe"], columns=list(ALN_FIELDS)))
    eq2 = pd.read_csv(out2 / "equiv_map.tsv", sep="\t", dtype={"accession": str})
    chk("헤더 없는 파일 + --columns 가 같은 결과", len(eq2) == len(eq),
        f"{len(eq2)} vs {len(eq)}")
    chk("--columns 경로가 params 에 기록",
        res2["audit"]["read"]["column_source"] == "--columns")
    lbl, ok, msg, full = expect_raises(
        "헤더 없는 파일을 --columns 없이 주면 예외",
        lambda: run(_args_for(base_args, alignment=str(syn["noheader"]),
                              master=str(syn["master"]), outdir=str(dirpath / "out_nh_fail"),
                              expected_n=0, columns=None)))
    chk(lbl, ok and "--columns" in full and "헤더가 없다" in full, msg)

    # --- 0-based 자동 판별
    out3 = dirpath / "out_zero"
    res3 = run(_args_for(base_args, alignment=str(syn["zerobased"]), master=str(syn["master"]),
                         fasta=str(dirpath / "nope.fasta"), outdir=str(out3), expected_n=0))
    bs3 = res3["parse"]["base"]
    chk("0-based 입력을 자동 판별",
        bs3["base"] == 0 and bs3["n_match_0based"] == 5 and bs3["n_match_1based"] == 0,
        f"채택 {bs3['base']}-based · 1based={bs3['n_match_1based']} 0based={bs3['n_match_0based']}")
    chk("채택 근거가 출력에 남는다", "0-based 채택" in bs3["basis"], bs3["basis"])
    eq3 = pd.read_csv(out3 / "equiv_map.tsv", sep="\t", dtype={"accession": str})
    sub_a = eq.loc[eq.accession == "SYN00002", ["q_pos", "t_pos", "q_aa", "t_aa"]].reset_index(drop=True)
    sub_b = eq3.loc[eq3.accession == "SYN00002", ["q_pos", "t_pos", "q_aa", "t_aa"]].reset_index(drop=True)
    chk("0-based 판본의 매핑이 1-based 판본과 동일", sub_a.equals(sub_b))
    chk("0-based 판본도 검증 3 통과", res3["parse"]["self_hit"]["passed"])

    # --- --matrix-aligned-only
    out4 = dirpath / "out_alignedonly"
    res4 = run(_args_for(base_args, alignment=str(syn["good"]), master=str(syn["master"]),
                         fasta=str(dirpath / "nope.fasta"), outdir=str(out4),
                         expected_n=0, matrix_aligned_only=True))
    chk("--matrix-aligned-only 면 행 5개", res4["parse"]["n_matrix_rows"] == 5,
        str(res4["parse"]["n_matrix_rows"]))

    # ---------------- 예외가 나야 하는 케이스 ----------------
    fail_evidence: list[tuple[str, bool, str]] = []   # (라벨, 발생여부, 예외 전문)
    lbl, ok, msg, full = expect_raises(
        "qaln 이 참조와 어긋나면 예외 (검증 1)",
        lambda: run(_args_for(base_args, alignment=str(syn["bad_qaln"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_bad_qaln"), expected_n=0)))
    chk(lbl, ok and "검증 1 실패" in full and "SYNBAD1" in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    lbl, ok, msg, full = expect_raises(
        "len(qaln) != len(taln) 면 예외 (검증 2)",
        lambda: run(_args_for(base_args, alignment=str(syn["bad_len"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_bad_len"), expected_n=0)))
    chk(lbl, ok and "검증 2 실패" in full and "SYNBAD2" in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    lbl, ok, msg, full = expect_raises(
        "자기 히트가 부분 정렬이면 예외 (검증 3)",
        lambda: run(_args_for(base_args, alignment=str(syn["bad_self"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_bad_self"), expected_n=0)))
    chk(lbl, ok and "검증 3 실패" in full and "600/636" in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    lbl, ok, msg, full = expect_raises(
        "질의가 여러 종 섞이면 예외",
        lambda: run(_args_for(base_args, alignment=str(syn["multi_query"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_multi"), expected_n=0)))
    chk(lbl, ok and "질의가" in full and "P99999" in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    # 질의서열이 참조와 1잔기 어긋나 0-based 로 위장하는 경로 — t_pos 가 조용히 +1 밀린다.
    lbl, ok, msg, full = expect_raises(
        "질의서열 1잔기 어긋남(0-based 위장)이면 예외 — t_pos 조용한 +1 차단",
        lambda: run(_args_for(base_args, alignment=str(syn["shifted"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_shifted"), expected_n=0)))
    chk(lbl, ok and "색인 0 이 없다" in full and "t_pos" in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    # 같은 파일을 --assume-qstart-base 0 으로 강제하면 진행하되 불일치로 크게 적는다.
    res_sh = run(_args_for(base_args, alignment=str(syn["shifted"]), master=str(syn["master"]),
                           fasta=str(dirpath / "nope.fasta"),
                           outdir=str(dirpath / "out_shifted_forced"), expected_n=0,
                           assume_qstart_base="0"))
    chk("--assume-qstart-base 0 강제 시 t_pos +1 위험을 불일치로 보고",
        any("t_pos 가 전부 +1" in d for d in res_sh["discrepancies"]),
        str(res_sh["discrepancies"])[:200])

    lbl, ok, msg, full = expect_raises(
        "qlen 이 참조 길이와 다르면 예외",
        lambda: run(_args_for(base_args, alignment=str(syn["bad_qlen"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_bad_qlen"), expected_n=0,
                              columns=None)))
    chk(lbl, ok and "qlen" in full and str(len(ref)) in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    lbl, ok, msg, full = expect_raises(
        "taln 만 소문자면 예외 (n_identical 조용한 과소계산 차단)",
        lambda: run(_args_for(base_args, alignment=str(syn["lower_taln"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_lower"), expected_n=0)))
    chk(lbl, ok and "소문자" in full and "--upper-aln" in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    out6 = dirpath / "out_lower_upper"
    res6 = run(_args_for(base_args, alignment=str(syn["lower_taln"]), master=str(syn["master"]),
                         fasta=str(dirpath / "nope.fasta"), outdir=str(out6), expected_n=0,
                         upper_aln=True))
    cv6 = pd.read_csv(out6 / "position_coverage.tsv", sep="\t")
    mx6 = pd.read_csv(out6 / "residue_matrix.tsv", sep="\t", dtype=str,
                      keep_default_na=False, na_filter=False)
    chk("--upper-aln 이면 소문자 입력도 대문자로 매핑된다",
        "".join(mx6.loc[mx6.accession == "SYNLOW1"].iloc[0, 1:].tolist()) == syn["t1"].upper())
    chk("--upper-aln 후 n_identical 이 0 이 아니다 (과소계산 해소)",
        int(cv6.n_identical.sum()) > len(ref), str(int(cv6.n_identical.sum())))

    lbl, ok, msg, full = expect_raises(
        "헤더에 같은 컬럼명이 두 번이면 명확한 예외",
        lambda: run(_args_for(base_args, alignment=str(syn["dup_column"]),
                              master=str(syn["master"]), fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_dupcol"), expected_n=0)))
    chk(lbl, ok and "같은 컬럼명" in full and "qaln" in full, full[:400])
    fail_evidence.append((lbl, ok, full))

    lbl, ok, msg, full = expect_raises(
        "--outside-char 가 한 글자가 아니면 예외",
        lambda: run(_args_for(base_args, alignment=str(syn["good"]), master=str(syn["master"]),
                              fasta=str(dirpath / "nope.fasta"),
                              outdir=str(dirpath / "out_oc"), expected_n=0,
                              outside_char="..")))
    chk(lbl, ok and "한 글자" in full, full[:300])

    # master 의 중복 accession 으로 universe 가 줄어드는 경로 — 보존도의 분모가 작아진다.
    out7 = dirpath / "out_master_dup"
    res7 = run(_args_for(base_args, alignment=str(syn["good"]), master=str(syn["master_dup"]),
                         fasta=str(dirpath / "nope.fasta"), outdir=str(out7),
                         expected_n=4))
    chk("master 중복으로 universe 가 줄면 불일치로 보고",
        any("중복" in d and "n_universe" in d for d in res7["discrepancies"]),
        str(res7["discrepancies"])[:300])
    chk("universe 가 기대 수와 다르면 불일치로 보고",
        any("universe 가" in d and "기대" in d for d in res7["discrepancies"]),
        str(res7["discrepancies"])[:300])

    # 실패해도 감사 md 와 params.json(status=failed) 이 남는지
    fail_dir = dirpath / "out_bad_qaln"
    chk("실패해도 alignment_audit.md 가 남는다", (fail_dir / "alignment_audit.md").exists())
    if (fail_dir / "params.json").exists():
        fp = json.loads((fail_dir / "params.json").read_text())
        chk("실패 params.json status=failed", fp.get("status") == "failed", str(fp.get("status")))
        chk("실패 params.json 에 감사 결과 보존", "truncation" in fp.get("audit", {}))
    else:
        chk("실패 params.json 생성", False, "없음")

    # 완화 플래그
    out5 = dirpath / "out_self_relaxed"
    res5 = run(_args_for(base_args, alignment=str(syn["bad_self"]), master=str(syn["master"]),
                         fasta=str(dirpath / "nope.fasta"), outdir=str(out5), expected_n=0,
                         allow_partial_self_hit=True))
    chk("--allow-partial-self-hit 면 불일치 보고로 완화",
        any("검증 3 실패" in d for d in res5["discrepancies"]),
        str(res5["discrepancies"])[:200])

    # 감사 md 내용
    md = (out1 / "alignment_audit.md").read_text()
    chk("감사 md 에 누락·cband 경고 포함", "cband" in md and "max-seqs" in md)
    chk("감사 md 에 base 채택 근거 포함", "채택: 1-based" in md)
    chk("감사 md 에 검증 1/2/3 표 포함",
        "검증 1" in md and "검증 2" in md and "검증 3" in md)

    print("\n=== 예외가 나야 하는 케이스 — 실제 예외 증거 ===")
    for lbl, ok, msg in fail_evidence:
        body = "\n".join("        " + ln for ln in msg.splitlines())
        print(f"  {'OK  ' if ok else 'FAIL'} {lbl}\n{body}")

    print("\n=== 자기검증 결과 ===")
    bad = []
    for name, ok, detail in checks:
        print(f"  {'OK  ' if ok else 'FAIL'} {name}" + (f"   [{detail}]" if not ok else ""))
        if not ok:
            bad.append(f"{name} ({detail})")
    if bad:
        raise AssertionError(f"자기검증 실패 {len(bad)}건: {bad}")
    print(f"\n합성 데이터 자기검증 {len(checks)}항목 전부 통과. 산출물: {dirpath}")


# ==========================================================================
# CLI
# ==========================================================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="03_equiv_map.py",
        description="단계 3 — chcodh2_alignment.tsv 의 qaln/taln 을 파싱해 모든 상동체 잔기를 "
                    "ChCODH2 번호(1~636)로 옮긴다.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="q_pos = ChCODH2 번호 (Met1 = 1, 1~636). PDB 번호가 아니다.\n"
               "검증 1 qaln 갭제거 == 참조[qstart..] · 검증 2 len(qaln) == len(taln) · "
               "검증 3 자기히트 q_pos == t_pos 636개.\n"
               "데이터가 없는 곳에서는 --dry-run / --self-test 를 쓸 것.")
    p.add_argument("--alignment", default=str(common.ALIGN_TSV),
                   help="chcodh2_alignment.tsv (mmseqs search/convertalis 산출물)")
    p.add_argument("--master", default=str(common.MASTER_TSV),
                   help="universe(분모) 출처. accession 컬럼만 읽는다")
    p.add_argument("--fasta", default=str(common.IPR_FASTA),
                   help="master 가 없을 때의 universe 대체 출처")
    p.add_argument("--outdir", default=str(common.OUT_ROOT / STAGE), help="산출물 디렉터리")
    p.add_argument("--columns", default=None,
                   type=lambda s: [c.strip() for c in re.split(r"[,\s]+", s.strip()) if c.strip()],
                   help="헤더 없는 convertalis 산출물의 컬럼 순서를 직접 준다 "
                        "(mmseqs --format-output 순서. 쉼표 또는 공백 구분)")
    p.add_argument("--expected-n", "--expect-n", dest="expected_n", type=int,
                   default=EXPECTED_N_SEQUENCES,
                   help=f"기대 서열 수. 적으면 불일치로 보고 (기본 {EXPECTED_N_SEQUENCES}, "
                        f"0 이면 검사 안 함)")
    p.add_argument("--reference-acc", default=common.CHCODH2_ACC,
                   help=f"질의/자기 히트 accession (기본 {common.CHCODH2_ACC})")
    p.add_argument("--query-acc", default=None,
                   help="query 컬럼에 여러 질의가 섞인 파일에서 이 accession 의 행만 쓴다")
    p.add_argument("--assume-qstart-base", choices=["auto", "1", "0"], default="auto",
                   help="qstart 의 base. 기본 auto = 1-based/0-based 양쪽을 맞춰보고 채택")
    p.add_argument("--upper-aln", action="store_true",
                   help="qaln/taln 을 대문자로 변환해 비교·매핑한다 (소문자 마스킹 입력 대응)")
    p.add_argument("--outside-char", default=DEFAULT_GAP_CHAR,
                   help=f"정렬 블록 밖/정렬 없음을 나타낼 문자 (기본 '{DEFAULT_GAP_CHAR}'. "
                        f"바꾸면 행렬 알파벳이 달라지므로 하위 단계와 합의할 것)")
    p.add_argument("--matrix-aligned-only", action="store_true",
                   help="행렬 행을 정렬에 있는 accession 으로 제한 (기본은 universe 전체)")
    p.add_argument("--allow-partial-self-hit", action="store_true",
                   help="자기 히트가 완전 일치가 아니어도 불일치 보고로 계속한다 (기본은 예외)")
    p.add_argument("--max-examples", type=int, default=DEFAULT_MAX_EXAMPLES,
                   help=f"누락·불일치 예시 상한 (기본 {DEFAULT_MAX_EXAMPLES}, 잘리면 명시)")
    p.add_argument("--max-table-rows", type=int, default=DEFAULT_MAX_TABLE_ROWS,
                   help=f"감사 md 표 행 수 상한 (기본 {DEFAULT_MAX_TABLE_ROWS})")
    p.add_argument("--max-detail", type=int, default=DEFAULT_MAX_DETAIL,
                   help=f"params.json 상세 목록 상한 (기본 {DEFAULT_MAX_DETAIL})")
    p.add_argument("--strict", action="store_true",
                   help="불일치가 있으면 산출물을 모두 쓴 뒤 종료코드 2 로 끝낸다")
    p.add_argument("--dry-run", action="store_true", help="경로·규약만 확인하고 쓰지 않는다")
    p.add_argument("--self-test", metavar="DIR",
                   help="DIR 에 합성 정렬 TSV 를 만들어 전 경로를 실행하고 검증한다 "
                        "(예외가 나야 하는 케이스까지)")
    return p


def main(argv=None) -> int:
    enc = ensure_utf8_io()      # 한국어 메시지와 params.json 을 잃지 않기 위해 가장 먼저
    args = build_parser().parse_args(argv)
    if enc.get("utf8_mode") and enc.get("preferred_encoding", "").lower() != "utf-8":
        print(f"[인코딩] UTF-8 모드로 보정됨 (locale 기본값 {enc['preferred_encoding']})")
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
