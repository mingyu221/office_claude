"""Classic CoMFA fields, read from the project's own autogrid4 maps.

These are the physically honest counterpart to the CoMSIA fields: the very
grids AutoDock searched to produce the labels. ``*.C.map`` / ``*.A.map`` /
``*.N.map`` are Lennard-Jones interaction energies for an aliphatic carbon,
aromatic carbon and nitrogen probe; ``*.e.map`` is the Coulomb potential and
``*.d.map`` the desolvation term.

Two practical notes.

**The maps may have been deleted.** The docking protocol deletes them right
after docking (27 MB x 321 = 8 GB) and keeps only the ``.dlg``. Regenerating
them is cheap -- autogrid4 is grid computation, not search -- but it must use
the hand-written GPF from the docking brief. Never a GPF generator: those emit
a ``parameter_file None`` line and autogrid4 segfaults on it.

**Clamping is unavoidable here.** The 1/r^12 term diverges inside the protein,
so values are clipped to +-``ENERGY_CUTOFF``. Every buried node then shares one
constant, which is exactly the artefact CoMSIA was designed to avoid -- so
expect the two field families to behave differently in the buried shell, and
read that as a property of the descriptor, not a bug.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .. import config as C
from .grid import GridSpec

PROBE_SUFFIX = {
    "lj_C": "C.map", "lj_A": "A.map", "lj_N": "N.map",
    "elec": "e.map", "desolv": "d.map",
}
AUTOGRID_FIELDS = tuple(PROBE_SUFFIX)


def read_map(path: str | Path) -> tuple[np.ndarray, dict]:
    """Parse one autogrid map into an ``[x, y, z]`` array plus its header.

    AutoGrid writes values with x varying fastest, then y, then z, so the flat
    block reshapes as ``[z, y, x]`` and is transposed back.
    """
    path = Path(path)
    header: dict = {}
    values: list[float] = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            key = s.split()[0]
            if key == "SPACING":
                header["spacing"] = float(s.split()[1])
            elif key == "NELEMENTS":
                header["nelements"] = tuple(int(v) for v in s.split()[1:4])
            elif key == "CENTER":
                header["center"] = tuple(float(v) for v in s.split()[1:4])
            elif key in ("GRID_PARAMETER_FILE", "GRID_DATA_FILE", "MACROMOLECULE"):
                header[key.lower()] = s.split(maxsplit=1)[1] if " " in s else ""
            else:
                try:
                    values.append(float(s))
                except ValueError:
                    continue

    for required in ("spacing", "nelements", "center"):
        if required not in header:
            raise ValueError(f"{path.name} header has no {required.upper()} line")

    nx, ny, nz = header["nelements"]
    expected = (nx + 1) * (ny + 1) * (nz + 1)
    if len(values) != expected:
        raise ValueError(
            f"{path.name}: got {len(values)} values, expected {expected} "
            f"for NELEMENTS {nx} {ny} {nz}"
        )
    arr = np.asarray(values, dtype=np.float32).reshape(nz + 1, ny + 1, nx + 1)
    return arr.transpose(2, 1, 0), header


def extract_nodes(arr: np.ndarray, grid: GridSpec, npts: int = C.NPTS,
                  cutoff: float | None = C.ENERGY_CUTOFF) -> np.ndarray:
    """Pull the masked node subset out of a dense map and clamp it."""
    half = npts // 2
    idx = grid.index + half
    if idx.min() < 0 or (idx.max() >= np.asarray(arr.shape)).any():
        raise ValueError("grid nodes fall outside the map extent")
    vals = arr[idx[:, 0], idx[:, 1], idx[:, 2]].astype(np.float64)
    if cutoff is not None:
        vals = np.clip(vals, -cutoff, cutoff)
    return vals


def load_fields(tag: str, grid: GridSpec, maps_dir: str | Path = C.MAPS_DIR,
                fields: tuple[str, ...] = AUTOGRID_FIELDS,
                expect_center: np.ndarray | None = None,
                tol: float = 1e-3) -> dict[str, np.ndarray]:
    """Read every requested probe map for one receptor tag.

    ``expect_center`` is checked against the map header. A mismatch means the
    maps were generated for a different box than the one being modelled -- the
    single most likely way to produce a field matrix that silently describes
    the wrong volume.
    """
    maps_dir = Path(maps_dir)
    out: dict[str, np.ndarray] = {}
    for name in fields:
        if name not in PROBE_SUFFIX:
            raise ValueError(f"unknown autogrid field {name!r}")
        path = maps_dir / f"{tag}.{PROBE_SUFFIX[name]}"
        arr, header = read_map(path)
        if abs(header["spacing"] - C.SPACING) > 1e-6:
            raise ValueError(
                f"{path.name}: spacing {header['spacing']} != {C.SPACING}"
            )
        if header["nelements"] != (C.NPTS, C.NPTS, C.NPTS):
            raise ValueError(
                f"{path.name}: NELEMENTS {header['nelements']} != "
                f"{(C.NPTS,) * 3}"
            )
        if expect_center is not None:
            d = float(np.linalg.norm(
                np.asarray(header["center"]) - np.asarray(expect_center, float)))
            if d > tol:
                raise ValueError(
                    f"{path.name}: map centre is {d:.3f} A from the expected "
                    "box centre -- these maps are for a different box"
                )
        out[name] = extract_nodes(arr, grid)
    return out


def write_gpf(path: str | Path, tag: str, center: tuple[float, float, float],
              npts: int = C.NPTS, spacing: float = C.SPACING) -> Path:
    """Write the hand-rolled GPF from the docking brief.

    `map` lines must follow `ligand_types` order, and there must be no
    `parameter_file` line.
    """
    x, y, z = center
    lig = C.AUTOGRID_LIGAND_TYPES
    lines = [
        f"npts {npts} {npts} {npts}",
        f"gridfld {tag}.maps.fld",
        f"spacing {spacing}",
        f"receptor_types {C.AUTOGRID_RECEPTOR_TYPES}",
        f"ligand_types {' '.join(lig)}",
        f"receptor {tag}.pdbqt",
        f"gridcenter {x:.3f} {y:.3f} {z:.3f}",
        "smooth 0.5",
        *[f"map {tag}.{t}.map" for t in lig],
        f"elecmap {tag}.e.map",
        f"dsolvmap {tag}.d.map",
        "dielectric -0.1465",
    ]
    path = Path(path)
    path.write_text("\n".join(lines) + "\n")
    return path
