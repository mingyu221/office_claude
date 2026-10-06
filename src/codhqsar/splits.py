"""Cluster-wise splitting and cross-validation.

Rule 3 of the selection protocol: the split unit is the MMseqs2 cluster, never
the sequence. Members of one cluster are >=70% identical, so their docking
scores are nearly the same; letting siblings straddle train and test turns the
evaluation into a memorisation check.

SeekRank's stock random 5-fold violates this and inflates rho. Everything here
groups by cluster, including the inner folds used to pick a layer or a PLS
component count -- selection leakage through the inner loop is just as real as
through the test set (methods §9).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold

from . import config as C


def assert_no_cluster_leakage(train: pd.DataFrame, test: pd.DataFrame) -> None:
    shared = set(train["cluster"]) & set(test["cluster"])
    if shared:
        raise AssertionError(
            f"{len(shared)} clusters appear in both train and test: "
            f"{sorted(shared)[:5]}"
        )


def assert_nested(frames: list[set[str]]) -> None:
    """Experiment A/C/D1 curves must be nested so they read monotonically."""
    for i, (a, b) in enumerate(zip(frames, frames[1:])):
        if not a <= b:
            raise AssertionError(f"step {i} is not a subset of step {i + 1}")


def stratified_cluster_holdout(
    df: pd.DataFrame,
    frac: float = 0.25,
    strata: str = "cband",
    seed: int = C.SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Hold out ``frac`` of clusters per stratum.

    Only used to regenerate a split from scratch. The project's live split is
    already fixed in ``dock_targets.tsv`` -- prefer :func:`split_from_targets`,
    because a re-derived split is a different split even at the same seed.
    """
    rng = np.random.default_rng(seed)
    core = df[df.get("tier", "core") == "core"]
    test_clusters: list = []
    for _, grp in core.groupby(strata, observed=True):
        clusters = np.array(sorted(grp["cluster"].unique()))
        k = int(round(len(clusters) * frac))
        if k:
            test_clusters += list(rng.choice(clusters, size=k, replace=False))
    tc = set(test_clusters)
    test = df[df["cluster"].isin(tc)]
    train = df[~df["cluster"].isin(tc)]
    assert_no_cluster_leakage(train, test)
    return train, test


def split_from_targets(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Use the frozen ``split`` column. This is the canonical path.

    Expansion rows (tier == 'expansion') are second cluster members and belong
    to *training* clusters only, so they may join train (experiment D2) but can
    never enter test.
    """
    if "split" not in df.columns:
        raise KeyError("dock_targets.tsv has no `split` column")
    train = df[df["split"] == "train"].copy()
    test = df[df["split"] == "test"].copy()
    if "tier" in df.columns:
        stray = test[test["tier"] == "expansion"]
        if len(stray):
            raise AssertionError(
                f"{len(stray)} expansion rows are in test -- their cluster "
                "representatives are in train, so this is leakage"
            )
    assert_no_cluster_leakage(train, test)
    return train, test


def cluster_folds(df: pd.DataFrame, n_splits: int = 5):
    """GroupKFold over clusters. Use for every inner selection loop."""
    n_groups = df["cluster"].nunique()
    if n_groups < n_splits:
        raise ValueError(f"{n_groups} clusters < {n_splits} folds")
    gkf = GroupKFold(n_splits=n_splits)
    return list(gkf.split(df, groups=df["cluster"].to_numpy()))


def position_holdout(variants: pd.DataFrame, n_positions: int = 3,
                     seed: int = C.SEED) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Hold out whole mutated positions from the 288-variant set.

    The variants are a position x substitution grid. Splitting at random puts
    other substitutions of the same position into training, so the model
    memorises the position effect and the generalisation test collapses.
    """
    rng = np.random.default_rng(seed)
    positions = np.array(sorted(variants["pos"].astype(str).unique()))
    held = set(rng.choice(positions, size=min(n_positions, len(positions)),
                          replace=False))
    mask = variants["pos"].astype(str).isin(held)
    return variants[~mask].copy(), variants[mask].copy()
