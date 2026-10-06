#!/usr/bin/env python3
"""Score candidate structures with the fitted field model, gated by domain.

This is the step that makes the model reusable rather than a one-off report. It
takes any set of assembled dimers -- the ~1,500 remaining IPR010047 sequences,
ColabFoldDB hits, or de novo designs -- computes the same fields, applies the
same filter/scaler/PLS, and refuses to present an out-of-domain prediction as
if it were comparable.

A de novo sequence is, by construction, far from the training set. Without the
domain columns the ranking of such a list is decoration.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from codhqsar import config as C, tables as T
from codhqsar.domain import fit_domain
from codhqsar.features.fieldset import build_field_matrix, load_field_matrix
from codhqsar.km_bridge import calibrate_km
from codhqsar.labels import n_best, rank_better, to_raw
from codhqsar.models.pls_comfa import CoMFAModel
from codhqsar.splits import cluster_folds


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True,
                    help="TSV with accession and f41_pos columns")
    ap.add_argument("--assembled", type=Path, default=C.ASSEMBLED)
    ap.add_argument("--fields", type=Path, default=C.CACHE / "fields_comsia_SELF.npz")
    ap.add_argument("--km-table", type=Path, default=None,
                    help="optional km_merge.csv for an indicative Km column")
    ap.add_argument("--top", type=int, default=50)
    ap.add_argument("--out", type=Path, default=C.CACHE / "candidate_scores.tsv")
    args = ap.parse_args()

    # Refit on the full labelled set: for scoring candidates there is no reason
    # to hold 55 structures back, and the test-set numbers are already reported.
    train = T.build_dataset(require_scores=True)
    fm = load_field_matrix(args.fields)
    train = train[train["accession"].isin(fm.accessions)].copy()
    fm = fm.subset_rows(train["accession"].tolist())
    y = train[C.LABEL].to_numpy(float)

    folds = cluster_folds(train.reset_index(drop=True), n_splits=5)
    model = CoMFAModel(fm.blocks).fit(fm.X, y, folds=folds)
    print(f"refit on all {len(train)} labelled structures, "
          f"{model.chosen_.n_components} components "
          f"(rho_xval {model.chosen_.rho_xval:+.4f})")

    dom = fit_domain(model.latent(fm.X))
    dmodx_thr = model.fit_dmodx_threshold(fm.X)

    cand = pd.read_csv(args.candidates, sep="\t")
    for col in ("accession", "f41_pos"):
        if col not in cand.columns:
            raise KeyError(f"{args.candidates} needs a `{col}` column")
    print(f"scoring {len(cand)} candidates ...")

    cfm = build_field_matrix(cand, grid=fm.grid, source="comsia",
                             fields=tuple(fm.block_names),
                             assembled=args.assembled,
                             verify_center_dev=False, progress=True)

    out = pd.DataFrame({"accession": cfm.accessions})
    out["pred_score"] = model.predict(cfm.X)
    out[f"pred_{C.RAW_LABEL}"] = to_raw(out["pred_score"])
    rep = dom.report(model.latent(cfm.X))
    dmodx = model.x_residual_distance(cfm.X)
    out["leverage_ratio"] = rep["leverage_ratio"].round(3)
    out["knn_ratio"] = rep["knn_ratio"].round(3)
    out["dmodx_ratio"] = (dmodx / dmodx_thr).round(3)
    # Both gates, because they fail on different things: the latent measures
    # catch a point far along the directions the model uses, DModX catches a
    # structure the model's X-plane does not describe at all. A de novo design
    # typically trips the second, not the first.
    out["in_domain"] = rep["inside"] & (dmodx <= dmodx_thr)
    out["rank"] = rank_better(out["pred_score"]).astype(int)

    if args.km_table and args.km_table.exists():
        km = pd.read_csv(args.km_table)
        if {"km", C.RAW_LABEL} <= set(km.columns):
            cal = calibrate_km(-km[C.RAW_LABEL], km["km"])
            print(f"\nKm calibration: {cal.summary()}")
            print("  indicative only -- n=19 supports a scale, not a predictor")
            out["km_indicative"] = cal.predict(out["pred_score"])
            lo, hi = cal.interval(out["pred_score"])
            out["km_lo"], out["km_hi"] = lo, hi
            out["km_extrapolated"] = cal.extrapolating(out["pred_score"])

    out = out.sort_values("rank")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)

    inside = out[out["in_domain"]]
    print(f"\n{len(inside)}/{len(out)} candidates inside the applicability domain")
    print(f"\ntop {args.top} in-domain candidates:")
    print(n_best(inside, "pred_score", args.top).to_string(index=False))
    if len(out) - len(inside):
        print(f"\n{len(out) - len(inside)} out-of-domain predictions are in the "
              "file but excluded from this ranking -- the model is extrapolating "
              "there and their order is not comparable")
    print(f"\nwritten -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
