"""End-to-end run on synthetic structures.

The real data lives on another machine, so this is what establishes that the
wiring works: fake dimers with a *planted geometric* signal (acidic residues
crowded toward the box centre bind a dication better), then the whole chain --
fields, cluster-wise split, PLS with grouped component selection, held-out
evaluation, y-randomisation, contours, domain gating.

Planting the signal in geometry rather than in a descriptor matters. The model is
only credited if it finds the signal through the field representation, and the
y-randomisation check confirms the score is not just what component selection
buys from noise.
"""

from __future__ import annotations

import numpy as np
import pytest

from codhqsar import evaluate as E
from codhqsar import splits as S
from codhqsar.contours import export_contours, split_by_field, to_volume
from codhqsar.domain import fit_domain
from codhqsar.features.fields_autogrid import read_map, write_gpf
from codhqsar.features.fieldset import build_field_matrix, load_field_matrix
from codhqsar.features.grid import make_grid
from codhqsar.models.pls_comfa import CoMFAModel
from synth import make_series, write_fake_dimer

GRID = dict(stride=8, radius=9.0)      # coarse grid keeps the test quick


@pytest.fixture(scope="module")
def series(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("series")
    df, assembled = make_series(tmp, n_clusters=30)
    grid = make_grid(**GRID)
    fm = build_field_matrix(df, grid=grid, source="comsia",
                            assembled=assembled, verify_center_dev=False)
    return df, fm, assembled, tmp


def test_field_matrix_shape_and_finiteness(series):
    df, fm, _, _ = series
    assert fm.X.shape == (len(df), fm.grid.n_nodes * len(fm.block_names))
    assert np.isfinite(fm.X).all()
    assert len(fm.columns) == fm.X.shape[1]


def test_electrostatic_field_tracks_the_planted_bias(series):
    """Sanity check on the descriptor itself, before any model is fitted."""
    df, fm, _, _ = series
    blocks = split_by_field(np.arange(fm.X.shape[1]), fm.grid, fm.block_names)
    elec_cols = blocks["electrostatic"].astype(int)
    near = np.linalg.norm(fm.grid.offsets, axis=1) <= 6.0
    mean_elec = fm.X[:, elec_cols[near]].mean(axis=1)
    rho = E.spearman(df["bias"].to_numpy(), mean_elec)
    # more acidic crowding -> more negative potential at the site
    assert rho < -0.4, f"electrostatic field does not track the bias (rho={rho})"


def test_model_recovers_signal_on_held_out_clusters(series):
    df, fm, _, _ = series
    y = df["score"].to_numpy(float)
    tr = np.flatnonzero((df["split"] == "train").to_numpy())
    te = np.flatnonzero((df["split"] == "test").to_numpy())
    S.assert_no_cluster_leakage(df.iloc[tr], df.iloc[te])

    folds = S.cluster_folds(df.iloc[tr].reset_index(drop=True), n_splits=5)
    model = CoMFAModel(fm.blocks, max_components=4, min_sigma=1e-6)
    model.fit(fm.X[tr], y[tr], folds=folds)

    rho_test = E.spearman(y[te], model.predict(fm.X[te]))
    assert model.chosen_.rho_xval > 0.3, f"train CV failed: {model.chosen_}"
    assert rho_test > 0.3, f"held-out rho too low: {rho_test}"

    noise = model.y_randomization(fm.X[tr], y[tr], folds, n_rep=3)
    assert model.chosen_.rho_xval > noise["max"], (
        f"rho_xval {model.chosen_.rho_xval} does not clear the noise ceiling "
        f"{noise['max']} -- that would be selection bias, not signal"
    )


def test_filter_and_scaler_are_fitted_on_training_rows_only(series):
    """A refit on a subset must not inherit the full set's column choice."""
    df, fm, _, _ = series
    y = df["score"].to_numpy(float)
    tr = np.flatnonzero((df["split"] == "train").to_numpy())
    folds = S.cluster_folds(df.iloc[tr].reset_index(drop=True), n_splits=5)

    full = CoMFAModel(fm.blocks, max_components=3, min_sigma=1e-6)
    full.fit(fm.X, y, folds=S.cluster_folds(df.reset_index(drop=True), 5))
    part = CoMFAModel(fm.blocks, max_components=3, min_sigma=1e-6)
    part.fit(fm.X[tr], y[tr], folds=folds)
    assert full.scaler_.mean_.shape != part.scaler_.mean_.shape or not np.allclose(
        full.scaler_.mean_, part.scaler_.mean_
    ), "scaler statistics are identical -- it was not refitted on the subset"


def test_coefficients_align_with_columns_and_contour(series, tmp_path):
    df, fm, _, _ = series
    y = df["score"].to_numpy(float)
    folds = S.cluster_folds(df.reset_index(drop=True), n_splits=5)
    model = CoMFAModel(fm.blocks, max_components=3, min_sigma=1e-6).fit(
        fm.X, y, folds=folds)

    coef = model.coefficients()
    assert coef.shape[0] == fm.X.shape[1] == len(fm.columns)
    assert np.isfinite(coef).all()
    # filtered-out columns must be exactly zero
    assert np.all(coef[~model.filter_.keep_] == 0.0)

    out = export_contours(tmp_path / "ctr", coef, fm.X, fm.grid,
                          fm.block_names, fm.centers[fm.accessions[0]])
    for name in fm.block_names:
        assert out[name].exists()
    assert out["pymol"].read_text().startswith("#")

    vol = to_volume(split_by_field(coef, fm.grid, fm.block_names)["steric"], fm.grid)
    assert vol.shape == fm.grid.shape()


def test_domain_flags_a_distant_structure(series, tmp_path):
    """A structure built on a different scaffold should fall outside."""
    import pandas as pd

    df, fm, assembled, _ = series
    y = df["score"].to_numpy(float)
    folds = S.cluster_folds(df.reset_index(drop=True), n_splits=5)
    model = CoMFAModel(fm.blocks, max_components=3, min_sigma=1e-6).fit(
        fm.X, y, folds=folds)
    dom = fit_domain(model.latent(fm.X))

    odd_dir = tmp_path / "odd"
    odd_dir.mkdir()
    write_fake_dimer(odd_dir / "ODD1_dimer.pdb", seed=999, n_res=120,
                     charge_bias=0.0, center=(48.884, -6.486, 31.065))
    odd = build_field_matrix(
        pd.DataFrame([{"accession": "ODD1", "f41_pos": 41}]),
        grid=fm.grid, source="comsia", fields=tuple(fm.block_names),
        assembled=odd_dir, verify_center_dev=False)
    rep = dom.report(model.latent(odd.X))
    thr = model.fit_dmodx_threshold(fm.X)
    dmodx = model.x_residual_distance(odd.X)[0]

    # Latent-space leverage/kNN can miss this: the odd structure may project
    # near the centre of the few directions the model uses. DModX is the measure
    # that sees it, which is why both are reported.
    assert dmodx > thr, (
        f"DModX {dmodx:.4f} did not exceed the training threshold {thr:.4f}; "
        "a double-density scaffold read as in-domain"
    )
    assert np.isfinite(rep["knn_ratio"][0])


def test_field_matrix_cache_roundtrip(series, tmp_path):
    df, fm, _, _ = series
    path = tmp_path / "fields.npz"
    np.savez_compressed(
        path, X=fm.X.astype(np.float32), accessions=np.asarray(fm.accessions),
        columns=np.asarray(fm.columns), blocks=fm.blocks,
        block_names=np.asarray(fm.block_names), grid_index=fm.grid.index,
        grid_offsets=fm.grid.offsets, npts=fm.grid.npts,
        spacing=fm.grid.spacing, stride=fm.grid.stride, radius=fm.grid.radius,
        centers=np.vstack([fm.centers[a] for a in fm.accessions]),
    )
    back = load_field_matrix(path)
    assert back.accessions == fm.accessions
    assert back.block_names == fm.block_names
    assert np.allclose(back.X, fm.X, atol=1e-4)
    assert np.array_equal(back.grid.index, fm.grid.index)


def test_autogrid_roundtrip(tmp_path):
    """Map parsing must preserve AutoGrid's x-fastest ordering."""
    n = 4
    vals = [float(x * 100 + y * 10 + z)
            for z in range(n + 1) for y in range(n + 1) for x in range(n + 1)]
    header = ["GRID_PARAMETER_FILE t.gpf", "GRID_DATA_FILE t.maps.fld",
              "MACROMOLECULE t.pdbqt", "SPACING 0.375",
              f"NELEMENTS {n} {n} {n}", "CENTER 1.000 2.000 3.000"]
    p = tmp_path / "t.C.map"
    p.write_text("\n".join(header + [f"{v:.3f}" for v in vals]) + "\n")
    arr, hdr = read_map(p)
    assert arr.shape == (n + 1,) * 3
    assert hdr["center"] == (1.0, 2.0, 3.0)
    for x in range(n + 1):
        for y in range(n + 1):
            for z in range(n + 1):
                assert arr[x, y, z] == pytest.approx(x * 100 + y * 10 + z)


def test_gpf_has_no_parameter_file_line(tmp_path):
    """A GPF generator emits `parameter_file None` and autogrid4 segfaults."""
    text = write_gpf(tmp_path / "t.gpf", "t", (1.0, 2.0, 3.0)).read_text()
    assert "parameter_file" not in text
    assert "npts 80 80 80" in text
    assert "gridcenter 1.000 2.000 3.000" in text
    # map lines must follow ligand_types order
    lig = [l.split()[1:] for l in text.splitlines() if l.startswith("ligand_types")][0]
    maps = [l.split()[1].split(".")[-2] for l in text.splitlines()
            if l.startswith("map ")]
    assert maps == lig
