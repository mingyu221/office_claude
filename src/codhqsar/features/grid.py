"""Grid geometry shared by every field descriptor.

The grid is the docking box: npts 80, spacing 0.375 A, half-edge 15.000 A.
Reusing the docking box exactly is the point -- a field computed over a
different volume would not describe the volume the labels came from. (The box
edge is not free either: a half-edge of 18 A collapsed within-position
correlation from 0.93 to 0.78, which is why npts is 80 and not 108.)

Alignment, the thing CoMFA usually lives or dies on, is already solved here:

* **Orientation** is shared. All 321 assembled dimers were superposed onto
  1SU6, and 6OND was itself first moved into the 1SU6 frame, so the coordinate
  axes mean the same thing in every structure.
* **Origin** is each structure's own F41-equivalent CA, exactly as the docking
  box centre was defined. So descriptor column *j* means "this offset from this
  protein's own binding-site anchor" in every row.

That self-referenced origin is a real choice, not a detail, so both options are
implemented and the primary one is pre-registered:

``SELF``   origin = each structure's own F41-equivalent CA. **Primary.**
           Matches the docking box rule, so the field describes the same volume
           that was scored. ``center_dev`` (up to 9.69 A) then carries the axis
           shift as a covariate rather than smearing the grid.
``REF``    origin = 1SU6 PHE41 CA for every structure. Sensitivity analysis.
           More CoMFA-orthodox (one absolute frame) but a structure whose axis
           has moved 9 A gets its pocket pushed off-grid.

Pick before looking at labels. Choosing the one with the better test rho is the
same error as choosing an ESM layer on the test set, which this project already
measured: pure noise yields rho 0.283 that way.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .. import config as C

Origin = Literal["SELF", "REF"]


@dataclass(frozen=True)
class GridSpec:
    """A masked, strided subset of the docking grid, expressed as offsets."""

    offsets: np.ndarray        # (n_nodes, 3) offsets from the origin, Angstrom
    index: np.ndarray          # (n_nodes, 3) integer node indices, for export
    npts: int
    spacing: float
    stride: int
    radius: float

    @property
    def n_nodes(self) -> int:
        return self.offsets.shape[0]

    def coords(self, origin: np.ndarray) -> np.ndarray:
        """Absolute coordinates of every node for one structure."""
        return self.offsets + np.asarray(origin, dtype=float).reshape(1, 3)

    def shape(self) -> tuple[int, int, int]:
        """Dense shape of the strided box, for writing a volume back out."""
        n = self.npts // self.stride + 1
        return (n, n, n)


def make_grid(npts: int = C.NPTS, spacing: float = C.SPACING,
              stride: int = C.FIELD_STRIDE,
              radius: float | None = C.FIELD_RADIUS) -> GridSpec:
    """Build the node set.

    A full 81^3 box is 531,441 nodes; at 321 structures and six fields that is
    terabytes. Striding by 2 and keeping a sphere of ``radius`` around the
    centre brings it to ~10k nodes per field, which is both tractable and
    honest -- nodes 15 A away in the corner of the box are solvent.
    """
    if npts % stride:
        raise ValueError(f"stride {stride} must divide npts {npts}")
    half = npts // 2
    axis = np.arange(-half, half + 1, stride, dtype=np.int32)
    ix, iy, iz = np.meshgrid(axis, axis, axis, indexing="ij")
    idx = np.stack([ix.ravel(), iy.ravel(), iz.ravel()], axis=1)
    off = idx.astype(np.float64) * spacing
    if radius is not None:
        keep = np.linalg.norm(off, axis=1) <= radius
        idx, off = idx[keep], off[keep]
    return GridSpec(off, idx, npts, spacing, stride, radius or np.inf)


def box_center(structure, f41_pos: int, chain: str = "A",
               origin: Origin = "SELF") -> np.ndarray:
    """Resolve the grid origin for one structure.

    ``f41_pos`` comes from ``dock_targets.tsv``. Never search for residue 41 by
    number: homolog numbering differs, and ``resSeq == 41`` picks up an
    unrelated residue.
    """
    if origin == "REF":
        return np.asarray(C.REF_CENTER, dtype=float)
    ca = structure.ca(chain, int(f41_pos))
    if ca is None:
        raise ValueError(f"chain {chain} residue {f41_pos} has no CA atom")
    return np.asarray(ca, dtype=float)


def center_dev(center: np.ndarray) -> float:
    """Distance from the 1SU6 origin -- the published check value.

    Verifying this against ``dock_targets.tsv`` confirms the right residue was
    picked. The magnitude itself is not a pass/fail criterion: a large value is
    a property of that homolog, not an error (docking brief §5).
    """
    return float(np.linalg.norm(np.asarray(center, float) - np.asarray(C.REF_CENTER, float)))


def field_column_names(grid: GridSpec, fields: list[str]) -> list[str]:
    """Stable descriptor names so a fitted model can be reapplied later."""
    names: list[str] = []
    for f in fields:
        for i, j, k in grid.index:
            names.append(f"{f}:{i:+d},{j:+d},{k:+d}")
    return names
