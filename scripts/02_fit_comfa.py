#!/usr/bin/env python3
"""Fit the field model and put it against the baselines that matter.

The comparison set is not arbitrary. Composition at rho_xval 0.4218 is the
measured ceiling of what this label has yielded so far -- ESM-1b's best layer
beat it by one sixth of an error bar. So the field model is reported against
composition and against shell descriptors, and the y-randomisation number is
printed next to it, because selecting the best of several component counts buys
correlation from noise alone.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from codhqsar import config as C, evaluate as E, splits as S, tables as T
from codhqsar.domain import fit_domain
from codhqsar.features.composition import composition, read_fasta
from codhqsar.features.fieldset import load_field_matrix
from codhqsar.features.site import shell_descriptors
from codhqsar.models.pls_comfa import CoMFAModel
from codhqsar.pdbio import read_pdb


def _fit_report(name, X, y, df, train_idx, test_idx, blocks, folds, args):
    model = CoMFAModel(blocks, max_components=args.max_components,
                       min_sigma=args.min_sigma, select_on=args.select_on)
    model.fit(X[train_idx], y[train_idx], folds=folds)
    pred_test = model.predict(X[test_idx])

    rho, lo, hi = E.spearman_ci(y[test_idx], pred_test)
    chosen = model.chosen_
    print(f"\n--- {name} ---")
    print(f"  columns kept        {model.filter_.n_kept} / {X.shape[1]}")
    print(f"  components          {chosen.n_components}")
    print(f"  rho_xval (train)    {chosen.rho_xval:+.4f}")
    print(f"  q2       (train)    {chosen.q2:+.4f}")
    print(f"  rho_test             {rho:+.4f}  [{lo:+.4f}, {hi:+.4f}]  "
          f"(n={len(test_idx)})")

    if args.y_random:
        yr = model.y_randomization(X[train_idx], y[train_idx], folds,
                                   n_rep=args.y_random)
        print(f"  y-randomisation     mean {yr['mean']:+.4f}  max {yr['max']:+.4f}"
              f"   <- the model must clear this")
        if chosen.rho_xval <= yr["max"]:
            print("  !! rho_xval does not exceed the noise maximum: "
                  "this result is selection bias, not signal")

    return model, pred_test, rho


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", type=Path,
                    default=C.CACHE / "fields_comsia_SELF.npz")
    ap.add_argument("--max-components", type=int, default=8)
    ap.add_argument("--min-sigma", type=float, default=0.05)
    ap.add_argument("--select-on", choices=["rho", "q2"], default="rho")
    ap.add_argument("--y-random", type=int, default=10,
                    help="y-randomisation repeats; 0 to skip")
    ap.add_argument("--residualized", action="store_true",
                    help="sensitivity run on the artefact-residualised label")
    ap.add_argument("--out", type=Path, default=C.CACHE / "fit_report")
    args = ap.parse_args()

    df = T.build_dataset(require_scores=True)
    fm = load_field_matrix(args.fields)
    df = df[df["accession"].isin(fm.accessions)].copy()
    fm = fm.subset_rows(df["accession"].tolist())

    y = df[C.LABEL].to_numpy(float)
    if args.residualized:
        y = E.residualize(df)
        print("NOTE: sensitivity run on residualised labels. cband is NOT "
              "residualised -- similarity is a legitimate predictor.")
        ok = np.isfinite(y)
        df, y = df[ok].copy(), y[ok]
        fm = fm.subset_rows(df["accession"].tolist())

    train_mask = (df["split"] == "train").to_numpy()
    test_mask = (df["split"] == "test").to_numpy()
    train_idx = np.flatnonzero(train_mask)
    test_idx = np.flatnonzero(test_mask)
    S.assert_no_cluster_leakage(df[train_mask], df[test_mask])
    print(f"train {len(train_idx)} / test {len(test_idx)}  "
          f"({df.loc[train_mask, 'cluster'].nunique()} vs "
          f"{df.loc[test_mask, 'cluster'].nunique()} clusters, no overlap)")

    folds = S.cluster_folds(df[train_mask].reset_index(drop=True),
                            n_splits=min(5, df.loc[train_mask, "cluster"].nunique()))

    results: dict[str, float] = {}
    model, pred_field, rho_field = _fit_report(
        "3D field (CoMSIA/PLS)", fm.X, y, df, train_idx, test_idx,
        fm.blocks, folds, args)
    results["field"] = rho_field

    # ---- baseline: shell descriptors (alignment-free, cheap) -------------
    try:
        shells = np.vstack([
            shell_descriptors(read_pdb(Path(C.ASSEMBLED) / f"{a}_dimer.pdb"),
                              fm.centers[a])
            for a in fm.accessions
        ])
        _, _, rho_shell = _fit_report(
            "shell descriptors", shells, y, df, train_idx, test_idx,
            np.zeros(shells.shape[1], int), folds, args)
        results["shell"] = rho_shell
    except (FileNotFoundError, ValueError) as err:
        print(f"\n[shell baseline skipped: {err}]")

    # ---- baseline: 24-dim composition (the measured ceiling) ------------
    train_fa = C.SPLIT / "homolog_train.fasta"
    test_fa = C.SPLIT / "homolog_test.fasta"
    if train_fa.exists() and test_fa.exists():
        seqs = {**read_fasta(train_fa), **read_fasta(test_fa)}
        have = [a for a in fm.accessions if a in seqs]
        if len(have) == len(fm.accessions):
            comp = np.vstack([composition(seqs[a]) for a in fm.accessions])
            _, _, rho_comp = _fit_report(
                "composition (24-dim)", comp, y, df, train_idx, test_idx,
                np.zeros(comp.shape[1], int), folds, args)
            results["composition"] = rho_comp
        else:
            print(f"\n[composition skipped: {len(fm.accessions) - len(have)} "
                  "accessions missing from the FASTA files]")

    # ---- is the field model's edge claimable? ---------------------------
    print("\n=== comparison on the fixed 55-structure test set ===")
    for k, v in results.items():
        print(f"  {k:<14} rho_test {v:+.4f}")
    if "composition" in results:
        verdict = E.claimable(results["field"], results["composition"], len(test_idx))
        print(f"\n  field - composition = {verdict['delta']:+.4f}  "
              f"p = {verdict['p']:.4f}  ->  "
              f"{'claimable' if verdict['claimable'] else 'NOT claimable at n=55'}")
        print("  (a difference this test cannot resolve is not a result; say so)")

    # ---- banded report + applicability domain --------------------------
    out = df.loc[test_mask, ["accession", "cband", "cluster", C.LABEL]].copy()
    out["pred"] = pred_field
    print("\n=== per-band test rho (thin tails pooled) ===")
    print(E.report_by_band(out, C.LABEL, "pred").round(4).to_string(index=False))

    dom = fit_domain(model.latent(fm.X[train_idx]))
    rep = dom.report(model.latent(fm.X[test_idx]))
    thr = model.fit_dmodx_threshold(fm.X[train_idx])
    dmodx = model.x_residual_distance(fm.X[test_idx])
    out["leverage_ratio"] = rep["leverage_ratio"]
    out["knn_ratio"] = rep["knn_ratio"]
    out["dmodx_ratio"] = dmodx / thr
    # Latent-space leverage/kNN only see the directions the model uses, so a
    # structure unlike anything fitted can still project near the centre.
    # DModX catches that, and a point has to pass both to count as in-domain.
    out["in_domain"] = rep["inside"] & (dmodx <= thr)
    print(f"\napplicability domain: {int(out['in_domain'].sum())}/{len(test_idx)} "
          f"test structures inside")
    print(f"  latent h* = {dom.h_star:.4f}; DModX threshold = {thr:.4f}")
    print(f"  latent-only would pass {int(rep['inside'].sum())}; "
          f"DModX-only {int((dmodx <= thr).sum())}")

    args.out.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out / "test_predictions.tsv", sep="\t", index=False)
    pd.DataFrame([{"n_components": r.n_components, "q2": r.q2,
                   "rho_xval": r.rho_xval} for r in model.cv_table_]
                 ).to_csv(args.out / "cv_table.tsv", sep="\t", index=False)
    np.save(args.out / "coefficients.npy", model.coefficients())
    np.save(args.out / "train_rows.npy", train_idx)
    print(f"\nwritten -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
