"""Label sign handling.

Mirrors the helper surface of the project's ``sklib.py`` (handoff §14.3) so this
package drops into the existing notebooks. The point of routing everything
through helpers is that ``nsmallest``/``min``/``rank`` scattered across cells is
exactly where a silently reversed ranking comes from.

Raw ``ad_med_ds`` is negative and lower means stronger binding.
``score = -ad_med_ds`` is positive and higher means stronger binding.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import BETTER, LABEL_SIGN, NOISE_FLOOR, RAW_LABEL

DIR_WORD = "높을수록 강함" if BETTER == "high" else "낮을수록 강함"
DIR_SORT = False if BETTER == "high" else True   # ascending= for sort_values
LABEL_DESC = f"score = {LABEL_SIGN:+.0f} x {RAW_LABEL} ({DIR_WORD})"


def to_score(raw) -> np.ndarray | float:
    """Raw ``ad_med_ds`` -> training label."""
    if isinstance(raw, pd.Series):
        return LABEL_SIGN * raw
    return LABEL_SIGN * np.asarray(raw, dtype=float)


def to_raw(score) -> np.ndarray | float:
    """Training label -> raw ``ad_med_ds``, for reporting in kcal/mol."""
    if isinstance(score, pd.Series):
        return score / LABEL_SIGN
    return np.asarray(score, dtype=float) / LABEL_SIGN


def assert_score_sign(y) -> None:
    """Catch both a missing inversion and a double inversion.

    Raw values are all negative, so a correct label array is strictly positive.
    """
    arr = np.asarray(y, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        raise ValueError("label array is empty or all-NaN")
    if BETTER == "high" and arr.min() <= 0:
        raise AssertionError(
            f"score must be positive (got min {arr.min():.3f}) -- "
            "sign inversion missing or applied twice"
        )
    if BETTER == "low" and arr.max() >= 0:
        raise AssertionError(f"raw label must be negative (got max {arr.max():.3f})")


def is_better(value, reference) -> bool:
    return value > reference if BETTER == "high" else value < reference


def rank_better(s: pd.Series) -> pd.Series:
    """Rank with the strongest binder at 1."""
    return s.rank(ascending=(BETTER == "low"), method="min")


def sort_better(df: pd.DataFrame, col: str) -> pd.DataFrame:
    return df.sort_values(col, ascending=DIR_SORT)


def n_best(df: pd.DataFrame, col: str, n: int) -> pd.DataFrame:
    return sort_better(df, col).head(n)


def strongest(a) -> float:
    arr = np.asarray(a, dtype=float)
    return float(np.nanmax(arr) if BETTER == "high" else np.nanmin(arr))


def weakest(a) -> float:
    arr = np.asarray(a, dtype=float)
    return float(np.nanmin(arr) if BETTER == "high" else np.nanmax(arr))


def strong_pct(a, q: float) -> float:
    """Boundary value of the strongest ``q`` percent."""
    arr = np.asarray(a, dtype=float)
    pct = 100.0 - q if BETTER == "high" else q
    return float(np.nanpercentile(arr, pct))


def n_better_than(arr, value) -> int:
    a = np.asarray(arr, dtype=float)
    return int((a > value).sum() if BETTER == "high" else (a < value).sum())


def resolvable(delta: float) -> bool:
    """Is a score difference larger than the assay's own noise?

    Differences at or below :data:`NOISE_FLOOR` cannot be claimed -- the ruler
    has no finer graduation (methods §1).
    """
    return abs(float(delta)) > NOISE_FLOOR
