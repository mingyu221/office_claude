"""Evaluation, noise accounting and covariate control.

Three rules from the handoffs are enforced here rather than left to discipline:

1. **Report rho per group, never pooled across groups with different label
   centres.** Pooling inflates it for free off the between-group offset. This
   project has the precedent: a pooled rho of +0.730 reversed to -0.539 inside
   one subgroup.
2. **Nothing below the noise floor is claimable.** The docking assay resolves
   0.280 kcal/mol; a score difference under that reads as the same value.
3. **cband is never residualised away.** Similarity is a legitimate predictor.
   Technical covariates (RMSD, center_dev, clashes, pLDDT, source, mode) are
   fair game -- those are artefacts, and they are confounded with cband, which
   is exactly why they have to be in the regression.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C

# Technical covariates. cband is deliberately absent -- see module docstring.
TECH_COVARIATES = ("RMSD", "center_dev", "crit_clash", "pLDDT",
                   "struct_source", "mode")

# Bands with too few test members to carry their own confidence interval.
THIN_BANDS = ("≥60", ">=60", "<30")


def spearman(a, b) -> float:
    a, b = np.asarray(a, float).ravel(), np.asarray(b, float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3:
        return float("nan")
    return float(stats.spearmanr(a[ok], b[ok]).statistic)


def spearman_ci(a, b, n_boot: int = 2000, alpha: float = 0.05,
                seed: int = C.SEED) -> tuple[float, float, float]:
    """Point estimate and bootstrap percentile interval."""
    a, b = np.asarray(a, float).ravel(), np.asarray(b, float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    n = a.size
    if n < 4:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = spearman(a[idx], b[idx])
    boots = boots[np.isfinite(boots)]
    lo, hi = np.percentile(boots, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return spearman(a, b), float(lo), float(hi)


def report_by_band(df: pd.DataFrame, y: str, pred: str,
                   band: str = "cband") -> pd.DataFrame:
    """Per-band rho, with the thin tail bands pooled.

    ``>=60`` has 2 test members and ``<30`` has 3. A rank correlation on two
    points can only be +-1, so those are pooled into one 'tails' row instead of
    being given a meaningless coefficient of their own.
    """
    rows = []
    thin = df[df[band].astype(str).isin(THIN_BANDS)]
    for name, grp in df.groupby(band, observed=True):
        if str(name) in THIN_BANDS:
            continue
        r, lo, hi = spearman_ci(grp[y], grp[pred])
        rows.append({"band": str(name), "n": len(grp), "rho": r,
                     "lo": lo, "hi": hi})
    if len(thin):
        r, lo, hi = spearman_ci(thin[y], thin[pred])
        rows.append({"band": "tails (pooled)", "n": len(thin), "rho": r,
                     "lo": lo, "hi": hi})
    r, lo, hi = spearman_ci(df[y], df[pred])
    rows.append({"band": "ALL", "n": len(df), "rho": r, "lo": lo, "hi": hi})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- OLS
@dataclass
class OLSFit:
    names: list[str]
    beta: np.ndarray
    se: np.ndarray
    t: np.ndarray
    p: np.ndarray
    r2: float
    resid: np.ndarray
    dof: int
    mask: np.ndarray          # rows actually used (finite y and design)

    def table(self) -> pd.DataFrame:
        return pd.DataFrame({"term": self.names, "beta": self.beta,
                             "se": self.se, "t": self.t, "p": self.p})


def _design(df: pd.DataFrame, terms: list[str]) -> tuple[np.ndarray, list[str]]:
    cols = [np.ones(len(df))]
    names = ["intercept"]
    for term in terms:
        if term not in df.columns:
            continue
        s = df[term]
        # `is_numeric_dtype` rather than `== object`: pandas 3 gives string
        # columns a dedicated str dtype, so the object check silently treats
        # every categorical covariate as numeric and NaNs the whole design.
        if not pd.api.types.is_numeric_dtype(s):
            levels = sorted(str(x) for x in s.dropna().unique())
            for lv in levels[1:]:                 # first level is the reference
                cols.append((s.astype(str) == lv).to_numpy(float))
                names.append(f"{term}[{lv}]")
        else:
            v = pd.to_numeric(s, errors="coerce").to_numpy(float)
            fill = np.nanmean(v) if np.isfinite(v).any() else 0.0
            v = np.where(np.isfinite(v), v, fill)
            cols.append(v)
            names.append(term)
    return np.column_stack(cols), names


def ols(df: pd.DataFrame, y: str, terms: list[str]) -> OLSFit:
    X, names = _design(df, terms)
    yv = pd.to_numeric(df[y], errors="coerce").to_numpy(float)
    ok = np.isfinite(yv) & np.isfinite(X).all(axis=1)
    X, yv = X[ok], yv[ok]
    beta, *_ = np.linalg.lstsq(X, yv, rcond=None)
    resid = yv - X @ beta
    n, p = X.shape
    dof = max(1, n - p)
    sigma2 = float(resid @ resid) / dof
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.where(se > 0, beta / se, np.nan)
    pv = 2 * stats.t.sf(np.abs(t), dof)
    ss = float(((yv - yv.mean()) ** 2).sum()) if yv.size else 0.0
    r2 = 1.0 - float(resid @ resid) / ss if ss > 0 else float("nan")
    return OLSFit(names, beta, se, t, pv, r2, resid, dof, ok)


def covariate_regression(df: pd.DataFrame, y: str = C.LABEL,
                         include_cband: bool = True) -> OLSFit:
    """``score ~ cband + technical covariates``.

    RMSD and center_dev are required, not optional: assembly quality degrades
    with falling identity (RMSD 0.61 -> 2.51 across the bands), so without them
    the cband coefficient mixes real binding differences with assembly
    artefacts and cannot be read at all.
    """
    terms = list(TECH_COVARIATES)
    if include_cband and "cband" in df.columns:
        terms = ["cband", *terms]
    missing = [c for c in ("RMSD", "center_dev") if c not in df.columns]
    if missing:
        raise KeyError(
            f"{missing} absent -- the cband coefficient is uninterpretable "
            "without them (handoff §2.4)"
        )
    return ols(df, y, terms)


def residualize(df: pd.DataFrame, y: str = C.LABEL) -> np.ndarray:
    """Label with technical artefacts regressed out -- for sensitivity only.

    The primary analysis uses the raw label; this is the companion that shows
    whether a trend survives once assembly quality is accounted for. If the two
    disagree, that disagreement is itself the result worth reporting.
    """
    fit = ols(df, y, [c for c in TECH_COVARIATES if c in df.columns])
    out = np.full(len(df), np.nan)
    out[fit.mask] = fit.resid
    return out


# ------------------------------------------------------------- noise floor
@dataclass
class NoiseFloor:
    median_within: float
    mean_within: float
    p90_within: float
    n_pairs: int
    var_between: float
    snr_canonical: float
    snr_simple: float

    def verdict(self) -> str:
        m = self.median_within
        if m <= 0.1:
            return "라벨 그대로 사용 가능"
        if m <= 0.45:
            return "288 세트와 비슷한 수준. 진행 가능"
        if m <= 1.0:
            return "경계. 구조당 여러 번 조립해 평균이 필요"
        return "라벨이 잡음에 묻힘. 재설계 필요"


def noise_floor_from_pairs(scores: pd.DataFrame, cluster: str = "cluster",
                           value: str = C.LABEL,
                           tier: str = "tier") -> NoiseFloor:
    """Measure the assay's own resolution from representative/sibling pairs.

    Each expansion member is the second medoid of a training cluster, >=70%
    identical to the representative, predicted and assembled independently. The
    pair's score difference therefore carries run-to-run noise *plus* structure
    prediction *plus* assembly -- the two noise sources the 288-variant set never
    had, which is why its SNR cannot be reused here.

    This is an upper bound: real sequence differences are folded in too.
    """
    df = scores.dropna(subset=[value, cluster])
    within: list[float] = []
    between: list[float] = []
    for _, grp in df.groupby(cluster, observed=True):
        if tier in grp.columns:
            rep = grp[grp[tier] == "core"]
            exp = grp[grp[tier] == "expansion"]
            if len(rep) >= 1 and len(exp) >= 1:
                within.append(abs(float(rep[value].iloc[0]) - float(exp[value].iloc[0])))
                between.append(float(rep[value].iloc[0]))
        elif len(grp) >= 2:
            v = grp[value].to_numpy(float)
            within.append(abs(v[0] - v[1]))
            between.append(v[0])

    if not within:
        raise ValueError("no representative/expansion pairs found")
    w = np.asarray(within, float)
    var_between = float(np.var(np.asarray(between, float), ddof=1))
    mean_w2 = float((w ** 2).mean())
    med = float(np.median(w))
    return NoiseFloor(
        median_within=med,
        mean_within=float(w.mean()),
        p90_within=float(np.percentile(w, 90)),
        n_pairs=w.size,
        var_between=var_between,
        snr_canonical=var_between / mean_w2 if mean_w2 > 0 else float("inf"),
        snr_simple=(float(np.std(between, ddof=1)) / med) ** 2 if med > 0 else float("inf"),
    )


def claimable(rho_a: float, rho_b: float, n: int, alpha: float = 0.05) -> dict:
    """Is a rho difference bigger than its own sampling error?

    Uses the Fisher-z two-sample test on the same n. Two models evaluated on the
    same 55 structures are correlated, so this is conservative in the wrong
    direction -- treat it as a screen, not a verdict.
    """
    if n < 4:
        return {"delta": float("nan"), "p": float("nan"), "claimable": False}
    za, zb = np.arctanh(np.clip([rho_a, rho_b], -0.999999, 0.999999))
    se = np.sqrt(2.0 / (n - 3))
    z = (za - zb) / se
    p = float(2 * stats.norm.sf(abs(z)))
    return {"delta": float(rho_a - rho_b), "z": float(z), "p": p,
            "claimable": p < alpha}
