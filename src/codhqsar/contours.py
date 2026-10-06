"""Coefficient contours -- the reason to run a field model at all.

A PLS coefficient lives on a grid node, so it can be put back in space and
looked at. That is what turns a ranking model into a design instruction: instead
of "this homolog scores 5.2" you get "adding bulk here raises the score, negative
charge there lowers it", which is directly actionable for the variant and de novo
tracks.

Plotted quantity is ``coefficient x field standard deviation`` (the usual
``StDev*Coeff`` map), not the bare coefficient. A node that barely varies across
the series can carry a huge coefficient and mean nothing; weighting by the
observed spread shows where the model is actually responding to real variation.

Sign convention follows the project label: ``score = -ad_med_ds`` and higher is
stronger binding, so a **positive** contour marks a place where more of that
field property is predicted to bind EV better.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from . import config as C
from .features.grid import GridSpec


def stdev_coeff(coefficients: np.ndarray, X_train: np.ndarray) -> np.ndarray:
    """``coef * sigma`` per column, computed on training rows only."""
    coef = np.asarray(coefficients, float).ravel()
    X = np.asarray(X_train, float)
    if coef.shape[0] != X.shape[1]:
        raise ValueError(f"{coef.shape[0]} coefficients for {X.shape[1]} columns")
    sigma = X.std(axis=0, ddof=1) if X.shape[0] > 1 else np.zeros(X.shape[1])
    return coef * sigma


def to_volume(values: np.ndarray, grid: GridSpec) -> np.ndarray:
    """Scatter one field's masked node values back into the dense strided box."""
    values = np.asarray(values, float).ravel()
    if values.shape[0] != grid.n_nodes:
        raise ValueError(f"{values.shape[0]} values for {grid.n_nodes} nodes")
    half = grid.npts // 2
    vol = np.zeros(grid.shape(), dtype=np.float32)
    idx = (grid.index + half) // grid.stride
    vol[idx[:, 0], idx[:, 1], idx[:, 2]] = values
    return vol


def split_by_field(values: np.ndarray, grid: GridSpec,
                   block_names: list[str]) -> dict[str, np.ndarray]:
    """Cut a concatenated descriptor vector back into per-field blocks."""
    values = np.asarray(values, float).ravel()
    n = grid.n_nodes
    if values.shape[0] != n * len(block_names):
        raise ValueError(
            f"{values.shape[0]} values != {n} nodes x {len(block_names)} fields"
        )
    return {name: values[i * n:(i + 1) * n] for i, name in enumerate(block_names)}


def write_dx(path: str | Path, volume: np.ndarray, origin: np.ndarray,
             spacing: float) -> Path:
    """Write an OpenDX scalar volume (readable by PyMOL `load` and ChimeraX).

    DX stores values with the last axis varying fastest, which is the C order
    numpy already uses, so the array is written as-is.
    """
    vol = np.asarray(volume, dtype=np.float32)
    nx, ny, nz = vol.shape
    ox, oy, oz = (float(v) for v in origin)
    flat = vol.ravel(order="C")

    lines = [
        f"object 1 class gridpositions counts {nx} {ny} {nz}",
        f"origin {ox:.6f} {oy:.6f} {oz:.6f}",
        f"delta {spacing:.6f} 0.000000 0.000000",
        f"delta 0.000000 {spacing:.6f} 0.000000",
        f"delta 0.000000 0.000000 {spacing:.6f}",
        f"object 2 class gridconnections counts {nx} {ny} {nz}",
        f"object 3 class array type double rank 0 items {flat.size} data follows",
    ]
    for i in range(0, flat.size, 3):
        lines.append(" ".join(f"{v:.6e}" for v in flat[i:i + 3]))
    lines.append('attribute "dep" string "positions"')
    lines.append("object \"regular positions regular connections\" class field")
    lines.append("component \"positions\" value 1")
    lines.append("component \"connections\" value 2")
    lines.append("component \"data\" value 3")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return path


def export_contours(
    out_dir: str | Path,
    coefficients: np.ndarray,
    X_train: np.ndarray,
    grid: GridSpec,
    block_names: list[str],
    center: np.ndarray,
    reference_pdb: str | Path | None = None,
    level_percentile: float = 95.0,
) -> dict[str, Path]:
    """Write one .dx per field plus a PyMOL script that loads them.

    ``center`` should be the structure ``reference_pdb`` belongs to -- normally
    the ChCODH2 self-control, whose residue numbering really is 41, which makes
    it the natural backdrop for reading the maps.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sc = stdev_coeff(coefficients, X_train)
    per_field = split_by_field(sc, grid, block_names)

    half = grid.npts // 2
    corner = np.asarray(center, float) - half * grid.spacing
    step = grid.spacing * grid.stride

    written: dict[str, Path] = {}
    levels: dict[str, float] = {}
    for name, vals in per_field.items():
        vol = to_volume(vals, grid)
        written[name] = write_dx(out_dir / f"{name}.dx", vol, corner, step)
        nz = np.abs(vals[vals != 0])
        levels[name] = float(np.percentile(nz, level_percentile)) if nz.size else 0.0

    pml = [
        "# CoMFA/CoMSIA coefficient contours (StDev*Coeff).",
        f"# score = {C.LABEL_SIGN:+.0f} * {C.RAW_LABEL}; higher = stronger binding.",
        "# Positive (blue) : more of this property here -> predicted stronger binding.",
        "# Negative (red)  : more of this property here -> predicted weaker binding.",
        "bg_color white",
    ]
    if reference_pdb:
        pml += [f"load {reference_pdb}, ref", "hide everything, ref",
                "show cartoon, ref", "color grey80, ref"]
    for name, path in written.items():
        lv = levels[name] or 1e-6
        pml += [
            f"load {path.name}, map_{name}",
            f"isosurface pos_{name}, map_{name}, {lv:.6g}",
            f"isosurface neg_{name}, map_{name}, {-lv:.6g}",
            f"color skyblue, pos_{name}",
            f"color firebrick, neg_{name}",
            f"set transparency, 0.35, pos_{name}",
            f"set transparency, 0.35, neg_{name}",
            f"disable pos_{name}", f"disable neg_{name}",
        ]
    first = next(iter(written), None)
    if first:
        pml += [f"enable pos_{first}", f"enable neg_{first}"]
    pml.append("# enable/disable one field at a time to read them separately")
    (out_dir / "contours.pml").write_text("\n".join(pml) + "\n")
    written["pymol"] = out_dir / "contours.pml"
    return written


def top_nodes(coefficients: np.ndarray, X_train: np.ndarray, grid: GridSpec,
              block_names: list[str], n: int = 20):
    """The strongest contour nodes as a table, with coordinates in Angstrom.

    Useful as text next to the pictures: a node offset is what tells you which
    residue positions to look at.
    """
    import pandas as pd

    sc = stdev_coeff(coefficients, X_train)
    per_field = split_by_field(sc, grid, block_names)
    rows = []
    for name, vals in per_field.items():
        order = np.argsort(-np.abs(vals))[:n]
        for j in order:
            if vals[j] == 0:
                continue
            dx, dy, dz = grid.offsets[j]
            rows.append({
                "field": name, "value": float(vals[j]),
                "dx": float(dx), "dy": float(dy), "dz": float(dz),
                "dist_from_center": float(np.linalg.norm(grid.offsets[j])),
                "direction": "stronger" if vals[j] > 0 else "weaker",
            })
    return pd.DataFrame(rows).sort_values("value", key=np.abs, ascending=False)
