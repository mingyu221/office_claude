"""The Km link.

The modelling chain has two legs, and only the first one is actually fitted:

    sequence/structure --(field QSAR)--> predicted docking score
    predicted docking score --(validated proxy, n=19)--> Km

Leg 2 is *not* a trained model and must not be presented as one. There are 19
measured Km values. What they support is (a) confirming the docking score tracks
Km at all, and (b) a calibration that maps a predicted score onto a Km scale with
openly wide intervals. Fitting a Km regressor on 19 points, then reporting its
cross-validated performance, would be reporting the noise.

Two documented traps are handled explicitly.

**Group-mixing inflation.** A pooled rho across groups with different label
centres gains correlation for free from the between-group offset. The precedent
in this project is stark: pooled rho +0.730 on the Km set reversed to -0.539
inside the ``dchg <= -1`` subgroup. So ``validate_proxy`` reports per-group rho
whenever a group column is supplied, and labels the pooled value as reference
only.

**dchg as a confounder.** Charge change correlates with the docking score on its
own, so a raw correlation can be the charge term rediscovered.
``partial_spearman`` gives the correlation with dchg partialled out -- the
reported figure in the handoffs is the dchg-residual rho +0.600 (p=0.0066), and
that, not the raw +0.913, is the defensible number.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C
from .evaluate import spearman, spearman_ci


def partial_spearman(x, y, control) -> dict:
    """Spearman correlation of ``x`` and ``y`` with ``control`` partialled out.

    Rank-transform everything first, then correlate the residuals of each
    against the control. Rank-based so it stays robust and comparable to the
    project's primary metric.
    """
    x, y, c = (np.asarray(v, float).ravel() for v in (x, y, control))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(c)
    x, y, c = x[ok], y[ok], c[ok]
    n = x.size
    if n < 5:
        return {"rho": float("nan"), "p": float("nan"), "n": int(n)}

    rx, ry, rc = (stats.rankdata(v) for v in (x, y, c))
    A = np.column_stack([np.ones(n), rc])

    def resid(v):
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    ex, ey = resid(rx), resid(ry)
    if np.allclose(ex, 0) or np.allclose(ey, 0):
        return {"rho": float("nan"), "p": float("nan"), "n": int(n)}
    r = float(np.corrcoef(ex, ey)[0, 1])
    dof = n - 3
    if dof <= 0 or abs(r) >= 1:
        return {"rho": r, "p": float("nan"), "n": int(n)}
    t = r * np.sqrt(dof / max(1e-12, 1 - r * r))
    return {"rho": r, "p": float(2 * stats.t.sf(abs(t), dof)), "n": int(n)}


@dataclass
class ProxyReport:
    pooled_rho: float
    pooled_lo: float
    pooled_hi: float
    n: int
    per_group: pd.DataFrame | None
    partial: dict | None

    def headline(self) -> str:
        if self.partial and np.isfinite(self.partial.get("rho", np.nan)):
            return (f"dchg-partialled rho = {self.partial['rho']:+.3f} "
                    f"(p={self.partial['p']:.4f}, n={self.partial['n']}); "
                    f"pooled rho = {self.pooled_rho:+.3f} [reference only]")
        return (f"pooled rho = {self.pooled_rho:+.3f} "
                f"[{self.pooled_lo:+.3f}, {self.pooled_hi:+.3f}], n={self.n}")


def validate_proxy(df: pd.DataFrame, km: str = "km",
                   score: str = C.RAW_LABEL,
                   group: str | None = None,
                   control: str | None = "dchg") -> ProxyReport:
    """Does the docking score track measured Km?

    Note the sign carefully. Lower ``ad_med_ds`` means stronger binding, and
    lower Km means tighter apparent affinity, so a *positive* rho between the raw
    values is the expected direction. Pass ``score=C.LABEL`` instead and the
    expected sign flips.
    """
    d = df.dropna(subset=[km, score])
    r, lo, hi = spearman_ci(d[km], d[score])

    per_group = None
    if group and group in d.columns:
        rows = []
        for name, grp in d.groupby(group, observed=True):
            if len(grp) >= 4:
                g, glo, ghi = spearman_ci(grp[km], grp[score])
            else:
                g = spearman(grp[km], grp[score]); glo = ghi = float("nan")
            rows.append({"group": str(name), "n": len(grp), "rho": g,
                         "lo": glo, "hi": ghi})
        per_group = pd.DataFrame(rows)

    partial = None
    if control and control in d.columns:
        partial = partial_spearman(d[km], d[score], d[control])

    return ProxyReport(r, lo, hi, len(d), per_group, partial)


@dataclass
class KmCalibration:
    """Monotone map from predicted score to Km. Indicative, not a trained model."""

    slope: float
    intercept: float
    log10: bool
    n: int
    resid_sd: float
    loo_mae: float
    rho: float
    x_range: tuple[float, float]

    def predict(self, pred_score) -> np.ndarray:
        x = np.asarray(pred_score, float).ravel()
        y = self.intercept + self.slope * x
        return 10.0 ** y if self.log10 else y

    def interval(self, pred_score, k: float = 1.96) -> tuple[np.ndarray, np.ndarray]:
        """Crude prediction band from the residual spread. Wide on purpose."""
        x = np.asarray(pred_score, float).ravel()
        y = self.intercept + self.slope * x
        lo, hi = y - k * self.resid_sd, y + k * self.resid_sd
        return ((10.0 ** lo, 10.0 ** hi) if self.log10 else (lo, hi))

    def extrapolating(self, pred_score) -> np.ndarray:
        """Flag inputs outside the calibration's own score range."""
        x = np.asarray(pred_score, float).ravel()
        return (x < self.x_range[0]) | (x > self.x_range[1])

    def summary(self) -> str:
        unit = "log10(Km)" if self.log10 else "Km"
        return (f"{unit} = {self.intercept:+.4f} {self.slope:+.4f}*score  "
                f"(n={self.n}, rho={self.rho:+.3f}, LOO MAE={self.loo_mae:.4f} "
                f"in {unit} units)")


def calibrate_km(pred_score, km, log10: bool = True) -> KmCalibration:
    """Least-squares calibration with leave-one-out error.

    ``log10=True`` because Km spans orders of magnitude and its error is
    multiplicative. LOO MAE is reported rather than R^2: with n=19, R^2 invites
    the reader to treat this as a fitted predictor, and it is not one.
    """
    x = np.asarray(pred_score, float).ravel()
    y = np.asarray(km, float).ravel()
    ok = np.isfinite(x) & np.isfinite(y)
    if log10:
        ok &= y > 0
    x, y = x[ok], y[ok]
    n = x.size
    if n < 5:
        raise ValueError(f"need at least 5 paired points, got {n}")
    yt = np.log10(y) if log10 else y

    A = np.column_stack([np.ones(n), x])
    beta, *_ = np.linalg.lstsq(A, yt, rcond=None)
    resid = yt - A @ beta

    errs = []
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        b, *_ = np.linalg.lstsq(A[m], yt[m], rcond=None)
        errs.append(abs(yt[i] - (b[0] + b[1] * x[i])))

    return KmCalibration(
        slope=float(beta[1]), intercept=float(beta[0]), log10=log10, n=int(n),
        resid_sd=float(resid.std(ddof=2)) if n > 2 else float("nan"),
        loo_mae=float(np.mean(errs)), rho=spearman(x, yt),
        x_range=(float(x.min()), float(x.max())),
    )
