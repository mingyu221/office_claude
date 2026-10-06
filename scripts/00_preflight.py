#!/usr/bin/env python3
"""Pre-flight checks. Run before any label is looked at.

Mirrors the checklist in the experiment-design brief §6.3 and §12: dataset
integrity, nested curves, no test leakage, and -- once docking scores exist --
the noise floor, which has to be measured before any correlation is reported.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from codhqsar import config as C, evaluate as E, splits as S, tables as T


def check_experiment_sets(targets: pd.DataFrame) -> list[str]:
    problems: list[str] = []
    try:
        exp = T.load_experiment_sets()
    except FileNotFoundError as err:
        return [f"experiment_sets.tsv unreadable: {err}"]

    test = set(targets.loc[targets.get("split") == "test", "accession"])
    bleed = set(exp["accession"]) & test
    if bleed:
        problems.append(f"{len(bleed)} test accessions appear in a training set: "
                        f"{sorted(bleed)[:5]}")

    for (name, rep), grp in exp[exp["exp"].isin(["A", "C", "D1"])].groupby(["exp", "rep"]):
        steps = [set(g["accession"]) for _, g in grp.groupby("step")]
        try:
            S.assert_nested(steps)
        except AssertionError as err:
            problems.append(f"{name} rep{rep}: {err}")

    n_unique = exp["accession"].nunique()
    if n_unique != 266:
        problems.append(f"unique training sequences = {n_unique}, expected 266")

    print("\nstep sizes per experiment:")
    print(exp.groupby(["exp", "rep", "step"]).size()
             .unstack("step").fillna(0).astype(int).to_string())
    return problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--with-labels", action="store_true",
                    help="docking has finished; also measure the noise floor")
    args = ap.parse_args()

    print(f"BASE = {C.BASE}")
    problems: list[str] = []

    targets = T.load_targets()
    print(f"\ntargets: {len(targets)} rows "
          f"({targets.get('tier', pd.Series(dtype=str)).value_counts().to_dict()})")

    if "split" in targets.columns:
        train, test = S.split_from_targets(targets)
        print(f"split: train {len(train)} / test {len(test)}; "
              f"clusters {train['cluster'].nunique()}/{test['cluster'].nunique()}")
    else:
        problems.append("dock_targets.tsv has no `split` column")

    problems += check_experiment_sets(targets)

    if args.with_labels:
        df = T.build_dataset(require_scores=True)
        print(f"\nlabels: n={len(df)}  "
              f"score {df[C.LABEL].min():.3f}..{df[C.LABEL].max():.3f}  "
              f"IQR {df[C.LABEL].quantile(.75) - df[C.LABEL].quantile(.25):.3f}")
        nf = E.noise_floor_from_pairs(df)
        print(f"\nnoise floor from {nf.n_pairs} representative/sibling pairs")
        print(f"  median |delta|  {nf.median_within:.3f} kcal/mol "
              f"(documented AutoDock value {C.NOISE_FLOOR})")
        print(f"  mean / p90      {nf.mean_within:.3f} / {nf.p90_within:.3f}")
        print(f"  SNR canonical   {nf.snr_canonical:.2f}   (simple {nf.snr_simple:.2f})")
        print(f"  verdict         {nf.verdict()}")
        if nf.median_within > 1.0:
            problems.append("noise floor above 1.0 -- labels are buried in noise")

        fit = E.covariate_regression(df)
        print("\ncovariate regression  score ~ cband + technical covariates")
        print(fit.table().round(4).to_string(index=False))
        print(f"  R2 = {fit.r2:.4f}")
        print("  (cband is uninterpretable unless RMSD and center_dev are in "
              "this model -- they are confounded with it)")

    print("\n" + "=" * 60)
    if problems:
        print(f"{len(problems)} PROBLEM(S):")
        for p in problems:
            print(f"  - {p}")
        return 1
    print("all pre-flight checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
