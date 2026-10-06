"""Readers for the canonical project tables.

``master_1827.tsv`` is the single source of truth; everything else is derived.
Readers are deliberately strict: a missing column fails loudly rather than
silently producing NaN, because the one serious data incident in this project
(23 structures assembled onto the wrong reference dimer) was a silent NaN fill
of ``arch`` and ``f41_pos``.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import config as C
from .labels import assert_score_sign, to_score


def _read(path: Path, required: list[str], sep: str = "\t") -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Set CODHQSAR_BASE to the project root "
            f"(currently {C.BASE})."
        )
    df = pd.read_csv(path, sep=sep)
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise KeyError(f"{path.name} is missing columns {missing}; has {list(df.columns)}")
    return df


def load_targets() -> pd.DataFrame:
    """321 docking targets: accession, tier, cband, split, cluster, f41_pos."""
    df = _read(C.T_TARGETS, ["accession", "tier", "cband", "f41_pos"])
    df["f41_pos"] = df["f41_pos"].astype("int64")
    if df["accession"].duplicated().any():
        dup = df.loc[df["accession"].duplicated(), "accession"].tolist()
        raise ValueError(f"duplicate accessions in dock_targets.tsv: {dup}")
    return df


def load_manifest() -> pd.DataFrame:
    """Assembly covariates. Korean column names are kept as written on disk."""
    df = _read(C.T_MANIFEST, ["accession"])
    # Normalise the two Korean covariate names to ASCII for modelling code.
    rename = {"임계충돌": "crit_clash", "충돌": "clash", "쌍": "pairs"}
    return df.rename(columns={k: v for k, v in rename.items() if k in df.columns})


def load_scores() -> pd.DataFrame:
    """Docking results. Adds the sign-corrected ``score`` column and checks it.

    Also verifies every row really is a 50-run median: a row with
    ``ad_n_ds != 50`` means ``--heuristics 0 --autostop 0`` was dropped, which
    breaks both the 50-run convention and seed reproducibility.
    """
    df = _read(C.T_SCORES, ["accession", C.RAW_LABEL])
    if "ad_n_ds" in df.columns:
        bad = df.loc[df["ad_n_ds"] != 50, "accession"].tolist()
        if bad:
            raise ValueError(
                f"{len(bad)} rows have ad_n_ds != 50 (first: {bad[:5]}). "
                "Re-dock those with --heuristics 0 --autostop 0."
            )
    df[C.LABEL] = to_score(df[C.RAW_LABEL])
    assert_score_sign(df[C.LABEL])
    return df


def load_master() -> pd.DataFrame:
    return _read(C.T_MASTER, ["accession", "cluster"])


def load_experiment_sets() -> pd.DataFrame:
    return _read(C.T_EXPERIMENTS, ["exp", "step", "rep", "accession"])


def load_variants() -> pd.DataFrame:
    """288 ChCODH2 point variants. Separate noise regime -- never pooled naively."""
    df = _read(C.T_VARIANTS, ["pos", "wt", "mut", C.RAW_LABEL], sep=",")
    df[C.LABEL] = to_score(df[C.RAW_LABEL])
    return df


def build_dataset(require_scores: bool = True) -> pd.DataFrame:
    """Join targets + covariates + labels into the single modelling frame.

    ``require_scores=False`` lets the dataset be assembled before docking
    finishes, so that splits and features can be built while labels are still
    blind -- which is the order the project's pre-registration demands.
    """
    df = load_targets()
    man = load_manifest()
    cov = [c for c in ("RMSD", "center_dev", "crit_clash", "clash", "mode",
                       "struct_source", "pLDDT", "d_arom") if c in man.columns]
    df = df.merge(man[["accession", *cov]], on="accession", how="left")

    if require_scores:
        sc = load_scores()
        keep = ["accession", C.RAW_LABEL, C.LABEL]
        keep += [c for c in ("ad_iqr_ds", "ad_best_ds") if c in sc.columns]
        df = df.merge(sc[keep], on="accession", how="left")
        missing = df[C.LABEL].isna().sum()
        if missing:
            raise ValueError(f"{missing} targets have no docking score yet")
        assert_score_sign(df[C.LABEL])

    if "cluster" not in df.columns:
        df = df.merge(load_master()[["accession", "cluster"]], on="accession", how="left")
    if df["cluster"].isna().any():
        n = int(df["cluster"].isna().sum())
        raise ValueError(f"{n} targets have no cluster -- cluster-wise splitting impossible")
    return df
