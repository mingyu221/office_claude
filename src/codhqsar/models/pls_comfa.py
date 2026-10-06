"""PLS field model -- the CoMFA/CoMSIA regressor.

PLS and not ridge/RF, for the reason CoMFA has always used it: there are tens of
thousands of correlated descriptors and ~157 rows, so an ordinary least-squares
family model has no unique solution and a tree ensemble cannot produce a
coefficient field to contour. PLS extracts a handful of latent directions and
its coefficients map straight back onto grid nodes, which is where the whole
interpretive payoff lives.

Component count is chosen by cluster-grouped cross-validation on the training
set only. Three guards travel with the model, because a field model with 59k
descriptors can fit anything:

``q2``              cross-validated R^2. The number a CoMFA reviewer looks for.
``rho_xval``        cross-validated Spearman -- this project's primary metric,
                    since the question is which homolog binds better, not the
                    absolute kcal/mol.
``y_randomization`` refit on shuffled labels. Tells you how much of the score is
                    the selection procedure rather than signal. This project has
                    already measured that selecting the best of 34 ESM layers on
                    noise alone buys rho 0.12-0.19, so the check is not optional.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.stats import spearmanr
from sklearn.cross_decomposition import PLSRegression

from .. import config as C
from ..features.fieldset import BlockScaler, ColumnFilter


def _spearman(a, b) -> float:
    a, b = np.asarray(a, float).ravel(), np.asarray(b, float).ravel()
    if a.size < 3 or np.allclose(a, a[0]) or np.allclose(b, b[0]):
        return float("nan")
    return float(spearmanr(a, b).statistic)


def _q2(y_true, y_pred) -> float:
    """Cross-validated R^2: 1 - PRESS/SS. Can go negative, and should be shown
    when it does."""
    y_true = np.asarray(y_true, float).ravel()
    y_pred = np.asarray(y_pred, float).ravel()
    press = float(((y_true - y_pred) ** 2).sum())
    ss = float(((y_true - y_true.mean()) ** 2).sum())
    return 1.0 - press / ss if ss > 0 else float("nan")


@dataclass
class CVResult:
    n_components: int
    q2: float
    rho_xval: float
    rho_folds: list[float] = field(default_factory=list)


class CoMFAModel:
    """Column filter -> block scaling -> PLS, with grouped component selection."""

    def __init__(
        self,
        blocks: np.ndarray,
        max_components: int = 8,
        min_sigma: float = 0.05,
        select_on: str = "rho",
        n_splits: int = 5,
        seed: int = C.SEED,
    ) -> None:
        if select_on not in ("rho", "q2"):
            raise ValueError("select_on must be 'rho' or 'q2'")
        self.blocks_in = np.asarray(blocks)
        self.max_components = max_components
        self.min_sigma = min_sigma
        self.select_on = select_on
        self.n_splits = n_splits
        self.seed = seed

        self.filter_: ColumnFilter | None = None
        self.scaler_: BlockScaler | None = None
        self.pls_: PLSRegression | None = None
        self.cv_table_: list[CVResult] = []
        self.chosen_: CVResult | None = None
        self.y_mean_: float = 0.0
        self.dmodx_threshold_: float | None = None

    # ------------------------------------------------------------ internals
    def _prep_fit(self, X: np.ndarray):
        self.filter_ = ColumnFilter(min_sigma=self.min_sigma).fit(X)
        Xf = self.filter_.transform(X)
        blocks = self.blocks_in[self.filter_.keep_]
        self.scaler_ = BlockScaler(blocks).fit(Xf)
        return self.scaler_.transform(Xf), blocks

    def _prep_apply(self, X: np.ndarray) -> np.ndarray:
        if self.filter_ is None or self.scaler_ is None:
            raise RuntimeError("model is not fitted")
        return self.scaler_.transform(self.filter_.transform(X))

    def _cap(self, n_rows: int, n_cols: int) -> int:
        return max(1, min(self.max_components, n_rows - 2, n_cols))

    # ------------------------------------------------------------- fitting
    def cross_validate(self, X: np.ndarray, y: np.ndarray,
                       folds: list[tuple[np.ndarray, np.ndarray]]) -> list[CVResult]:
        """Grouped CV over component counts.

        The filter and the scaler are refitted inside every fold. Fitting them
        once outside would let held-out rows set the descriptor set and the
        scaling -- a quiet leak that inflates every number downstream.
        """
        X = np.asarray(X, float)
        y = np.asarray(y, float).ravel()
        cap = self._cap(X.shape[0], X.shape[1])
        table: list[CVResult] = []

        for k in range(1, cap + 1):
            pred = np.full(y.shape, np.nan)
            rho_folds: list[float] = []
            for tr, te in folds:
                cf = ColumnFilter(min_sigma=self.min_sigma).fit(X[tr])
                Xtr, Xte = cf.transform(X[tr]), cf.transform(X[te])
                bs = BlockScaler(self.blocks_in[cf.keep_]).fit(Xtr)
                Xtr, Xte = bs.transform(Xtr), bs.transform(Xte)
                kk = max(1, min(k, Xtr.shape[0] - 2, Xtr.shape[1]))
                pls = PLSRegression(n_components=kk, scale=False)
                pls.fit(Xtr, y[tr] - y[tr].mean())
                pred[te] = pls.predict(Xte).ravel() + y[tr].mean()
                rho_folds.append(_spearman(y[te], pred[te]))
            table.append(CVResult(k, _q2(y, pred), _spearman(y, pred), rho_folds))

        self.cv_table_ = table
        return table

    def fit(self, X: np.ndarray, y: np.ndarray,
            folds: list[tuple[np.ndarray, np.ndarray]] | None = None,
            n_components: int | None = None) -> "CoMFAModel":
        X = np.asarray(X, float)
        y = np.asarray(y, float).ravel()

        if n_components is None:
            if folds is None:
                raise ValueError("pass folds (cluster-grouped) or n_components")
            table = self.cross_validate(X, y, folds)
            key = (lambda r: r.rho_xval) if self.select_on == "rho" else (lambda r: r.q2)
            scored = [r for r in table if np.isfinite(key(r))]
            if not scored:
                raise RuntimeError("every component count gave a non-finite score")
            self.chosen_ = max(scored, key=key)
            n_components = self.chosen_.n_components

        Xs, _ = self._prep_fit(X)
        k = max(1, min(n_components, Xs.shape[0] - 2, Xs.shape[1]))
        self.y_mean_ = float(y.mean())
        self.pls_ = PLSRegression(n_components=k, scale=False)
        self.pls_.fit(Xs, y - self.y_mean_)
        return self

    # ----------------------------------------------------------- inference
    def predict(self, X: np.ndarray) -> np.ndarray:
        if self.pls_ is None:
            raise RuntimeError("model is not fitted")
        return self.pls_.predict(self._prep_apply(X)).ravel() + self.y_mean_

    def latent(self, X: np.ndarray) -> np.ndarray:
        """PLS scores -- the space the applicability domain is measured in."""
        if self.pls_ is None:
            raise RuntimeError("model is not fitted")
        return self.pls_.transform(self._prep_apply(X))

    # ---------------------------------------------------------- diagnostics
    def coefficients(self) -> np.ndarray:
        """Regression coefficients expanded back to the full column set.

        Filtered-out columns get exactly 0, so the vector lines up with the
        original descriptor names and can be reshaped onto grid nodes.
        """
        if self.pls_ is None or self.filter_ is None:
            raise RuntimeError("model is not fitted")
        full = np.zeros(self.filter_.keep_.shape[0], dtype=float)
        coef = np.asarray(self.pls_.coef_, dtype=float).ravel()
        # undo block scaling so coefficients are per raw field unit
        blocks = self.blocks_in[self.filter_.keep_]
        scaled = coef.copy()
        for b, s in self.scaler_.block_scale_.items():
            scaled[blocks == b] /= s
        full[self.filter_.keep_] = scaled
        return full

    def x_residual_distance(self, X: np.ndarray) -> np.ndarray:
        """DModX -- how far a sample sits from the model's X-plane.

        Leverage and nearest-neighbour distance are measured in the latent
        space, so they only see the directions the model actually uses. A
        structure can be grossly unlike anything in the training set and still
        project near the centre of that space; the test fixture does exactly
        that with a double-density scaffold. DModX is the complementary measure:
        it is the part of X the latent components fail to explain, which is what
        catches "nothing like this was ever fitted".

        Returned as the per-sample RMS residual over the retained columns.
        """
        if self.pls_ is None:
            raise RuntimeError("model is not fitted")
        Xs = self._prep_apply(X)
        T = self.pls_.transform(Xs)
        P = np.asarray(self.pls_.x_loadings_, dtype=float)
        E = Xs - T @ P.T
        dof = max(1, Xs.shape[1] - P.shape[1])
        return np.sqrt((E ** 2).sum(axis=1) / dof)

    def fit_dmodx_threshold(self, X_train: np.ndarray,
                            percentile: float = 95.0) -> float:
        """Threshold for :meth:`x_residual_distance`, from training rows."""
        d = self.x_residual_distance(X_train)
        self.dmodx_threshold_ = float(np.percentile(d, percentile))
        return self.dmodx_threshold_

    def y_randomization(self, X: np.ndarray, y: np.ndarray,
                        folds: list[tuple[np.ndarray, np.ndarray]],
                        n_rep: int = 10, seed: int | None = None) -> dict[str, float]:
        """Refit on shuffled labels and report what selection alone buys.

        The reported model is only interesting to the extent it clears this.
        """
        rng = np.random.default_rng(self.seed if seed is None else seed)
        y = np.asarray(y, float).ravel()
        best: list[float] = []
        key = (lambda r: r.rho_xval) if self.select_on == "rho" else (lambda r: r.q2)
        for _ in range(n_rep):
            ys = rng.permutation(y)
            probe = CoMFAModel(self.blocks_in, self.max_components,
                               self.min_sigma, self.select_on, self.n_splits,
                               self.seed)
            table = probe.cross_validate(X, ys, folds)
            vals = [key(r) for r in table if np.isfinite(key(r))]
            best.append(max(vals) if vals else float("nan"))
        arr = np.asarray(best, float)
        return {"mean": float(np.nanmean(arr)), "max": float(np.nanmax(arr)),
                "values": arr}
