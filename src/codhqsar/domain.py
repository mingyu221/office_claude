"""Applicability domain.

A QSAR prediction is only defensible inside the chemical/sequence space the
model was fitted on. This project will eventually score ~1,500 remaining
IPR010047 sequences, metagenomic hits, and de novo designs -- the last of those
is far outside anything in the training set, so an AD check is not optional
there, it is the difference between a ranked candidate list and noise.

Two complementary measures:

``leverage``   classical QSAR h_i = x_i' (X'X)^-1 x_i with the 3p/n warning
               line. Cheap, but assumes the design matrix is well conditioned,
               so it is computed in the model's latent (PLS) space.
``knn_dist``   mean distance to the k nearest training points. Robust and it
               is what actually catches a de novo sequence that sits in a void.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class DomainModel:
    """Fitted applicability domain. Fit on training latents only."""

    train_latent: np.ndarray
    cov_inv: np.ndarray
    h_star: float
    knn_k: int
    knn_threshold: float

    def leverage(self, latent: np.ndarray) -> np.ndarray:
        x = np.atleast_2d(np.asarray(latent, dtype=float))
        return np.einsum("ij,jk,ik->i", x, self.cov_inv, x)

    def knn_distance(self, latent: np.ndarray) -> np.ndarray:
        x = np.atleast_2d(np.asarray(latent, dtype=float))
        d = np.linalg.norm(x[:, None, :] - self.train_latent[None, :, :], axis=2)
        k = min(self.knn_k, d.shape[1])
        return np.sort(d, axis=1)[:, :k].mean(axis=1)

    def inside(self, latent: np.ndarray) -> np.ndarray:
        """True where both criteria pass. Conservative: needs both."""
        return (self.leverage(latent) <= self.h_star) & (
            self.knn_distance(latent) <= self.knn_threshold
        )

    def report(self, latent: np.ndarray) -> dict[str, np.ndarray]:
        h = self.leverage(latent)
        d = self.knn_distance(latent)
        return {
            "leverage": h,
            "leverage_ratio": h / self.h_star,
            "knn_dist": d,
            "knn_ratio": d / self.knn_threshold,
            "inside": (h <= self.h_star) & (d <= self.knn_threshold),
        }


def fit_domain(train_latent: np.ndarray, knn_k: int = 5,
               knn_percentile: float = 95.0) -> DomainModel:
    x = np.atleast_2d(np.asarray(train_latent, dtype=float))
    n, p = x.shape
    # Ridge-regularised inverse: latent dims can still be near-collinear.
    cov = x.T @ x
    cov_inv = np.linalg.pinv(cov + 1e-9 * np.eye(p))
    h_star = 3.0 * p / n

    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)
    k = min(knn_k, max(1, n - 1))
    own = np.sort(d, axis=1)[:, :k].mean(axis=1)
    knn_threshold = float(np.percentile(own, knn_percentile))
    return DomainModel(x, cov_inv, float(h_star), k, knn_threshold)
