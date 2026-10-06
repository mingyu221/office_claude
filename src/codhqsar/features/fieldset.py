"""Assembling the field matrix and pruning its columns.

A raw field matrix is ~59k columns for 157 rows. Most of those columns carry
nothing: nodes deep inside the protein are identical in every structure, and
nodes out in solvent are zero everywhere. Classic CoMFA prunes them with a
minimum-sigma filter before PLS, and so does this.

The filter is fitted on training rows only. That is not pedantry: deciding which
columns to keep by looking at all rows lets the test structures influence the
descriptor set, which is the same class of leakage that produced rho 0.283 from
pure noise when an ESM layer was chosen on the test set.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .. import config as C
from ..pdbio import read_pdb
from . import fields_autogrid as FA
from . import fields_comsia as FC
from .grid import GridSpec, Origin, box_center, center_dev, make_grid


@dataclass
class FieldMatrix:
    """Field descriptors for a set of structures."""

    X: np.ndarray                 # (n_rows, n_cols) float64
    accessions: list[str]
    columns: list[str]
    blocks: np.ndarray            # (n_cols,) int, which field each column is
    block_names: list[str]
    grid: GridSpec
    centers: dict[str, np.ndarray]

    def subset_rows(self, accessions: list[str]) -> "FieldMatrix":
        pos = {a: i for i, a in enumerate(self.accessions)}
        missing = [a for a in accessions if a not in pos]
        if missing:
            raise KeyError(f"{len(missing)} accessions not in matrix: {missing[:5]}")
        rows = [pos[a] for a in accessions]
        return FieldMatrix(self.X[rows], list(accessions), self.columns,
                           self.blocks, self.block_names, self.grid,
                           {a: self.centers[a] for a in accessions})


def structure_path(accession: str, assembled: Path = C.ASSEMBLED) -> Path:
    """Locate the assembled dimer.

    Only one naming rule applies to assembled structures --
    ``{accession}_dimer.pdb`` -- unlike the monomer directories, where AFDB
    files carry the accession in the middle of the filename.
    """
    direct = assembled / f"{accession}_dimer.pdb"
    if direct.exists():
        return direct
    hits = sorted(assembled.glob(f"{accession}*_dimer.pdb"))
    if hits:
        return hits[0]
    raise FileNotFoundError(f"no assembled dimer for {accession} under {assembled}")


def build_field_matrix(
    targets,                      # DataFrame with accession + f41_pos
    grid: GridSpec | None = None,
    source: str = "comsia",
    fields: tuple[str, ...] | None = None,
    origin: Origin = "SELF",
    assembled: Path = C.ASSEMBLED,
    maps_dir: Path = C.MAPS_DIR,
    verify_center_dev: bool = True,
    progress: bool = False,
) -> FieldMatrix:
    """Compute or read fields for every row of ``targets``.

    ``verify_center_dev`` re-derives ``center_dev`` and compares it to the
    published column. That is the one check the docking brief asks for: it
    confirms the right residue was picked. The magnitude is never a pass/fail
    criterion -- a large value is a property of that homolog.
    """
    grid = grid or make_grid()
    if source == "comsia":
        fields = fields or FC.FIELD_NAMES
    elif source == "autogrid":
        fields = fields or FA.AUTOGRID_FIELDS
    else:
        raise ValueError("source must be 'comsia' or 'autogrid'")

    rows: list[np.ndarray] = []
    accs: list[str] = []
    centers: dict[str, np.ndarray] = {}
    mismatches: list[tuple[str, float, float]] = []

    for n, rec in enumerate(targets.itertuples(index=False), start=1):
        acc = str(rec.accession)
        st = read_pdb(structure_path(acc, assembled))
        ctr = box_center(st, int(rec.f41_pos), origin=origin)
        centers[acc] = ctr

        if verify_center_dev and hasattr(rec, "center_dev"):
            published = getattr(rec, "center_dev")
            if published is not None and np.isfinite(published):
                got = center_dev(ctr)
                if abs(got - float(published)) > 0.05:
                    mismatches.append((acc, got, float(published)))

        if source == "comsia":
            fmap = FC.compute_fields(st, grid, ctr, fields=fields)
        else:
            fmap = FA.load_fields(acc, grid, maps_dir=maps_dir, fields=fields,
                                  expect_center=ctr)
        rows.append(np.concatenate([fmap[f] for f in fields]))
        accs.append(acc)
        if progress and n % 25 == 0:
            print(f"  {n}/{len(targets)}", flush=True)

    if mismatches:
        lines = "\n".join(f"  {a}: computed {g:.3f} vs published {p:.3f}"
                          for a, g, p in mismatches[:10])
        raise ValueError(
            f"{len(mismatches)} structures disagree with the published "
            f"center_dev -- the wrong residue was picked:\n{lines}"
        )

    X = np.vstack(rows)
    blocks = np.concatenate([np.full(grid.n_nodes, i) for i in range(len(fields))])
    columns = [f"{f}:{i:+d},{j:+d},{k:+d}"
               for f in fields for i, j, k in grid.index]
    return FieldMatrix(X, accs, columns, blocks, list(fields), grid, centers)


class ColumnFilter:
    """Minimum-sigma column filter, fitted on training rows only.

    Drops columns whose standard deviation across the training series is below
    ``min_sigma``. A column that does not vary cannot explain a label that does,
    and keeping tens of thousands of them only dilutes the PLS weights.
    """

    def __init__(self, min_sigma: float = 0.05,
                 max_const_frac: float = 0.98) -> None:
        self.min_sigma = min_sigma
        self.max_const_frac = max_const_frac
        self.keep_: np.ndarray | None = None

    def fit(self, X: np.ndarray) -> "ColumnFilter":
        X = np.asarray(X, dtype=float)
        sigma = X.std(axis=0, ddof=1) if X.shape[0] > 1 else np.zeros(X.shape[1])
        keep = sigma >= self.min_sigma
        # also drop columns pinned to a single value in nearly every row
        if X.shape[0] > 2:
            mode_frac = np.array([
                (np.isclose(col, np.median(col))).mean() for col in X.T
            ])
            keep &= mode_frac < self.max_const_frac
        if not keep.any():
            raise ValueError(
                f"min_sigma={self.min_sigma} removed every column; "
                "lower it or check the field computation"
            )
        self.keep_ = keep
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.keep_ is None:
            raise RuntimeError("ColumnFilter is not fitted")
        return np.asarray(X, dtype=float)[:, self.keep_]

    def fit_transform(self, X: np.ndarray) -> np.ndarray:
        return self.fit(X).transform(X)

    @property
    def n_kept(self) -> int:
        return int(self.keep_.sum()) if self.keep_ is not None else 0


class BlockScaler:
    """CoMFA block scaling.

    Steric and electrostatic fields carry different units and very different
    magnitudes. Without block scaling the larger-variance block dominates the
    PLS latent directions for a purely dimensional reason. Each block is divided
    by its own root-mean-variance so blocks enter on comparable footing, then
    columns are mean-centred -- the standard CoMFA_STD treatment.
    """

    def __init__(self, blocks: np.ndarray) -> None:
        self.blocks = np.asarray(blocks)
        self.mean_: np.ndarray | None = None
        self.block_scale_: dict[int, float] = {}

    def fit(self, X: np.ndarray) -> "BlockScaler":
        X = np.asarray(X, dtype=float)
        if self.blocks.shape[0] != X.shape[1]:
            raise ValueError(
                f"blocks has {self.blocks.shape[0]} entries for {X.shape[1]} columns"
            )
        self.mean_ = X.mean(axis=0)
        Xc = X - self.mean_
        for b in np.unique(self.blocks):
            m = self.blocks == b
            v = float(np.sqrt((Xc[:, m] ** 2).mean())) if m.any() else 1.0
            self.block_scale_[int(b)] = v if v > 1e-12 else 1.0
        return self

    def transform(self, X: np.ndarray) -> np.ndarray:
        if self.mean_ is None:
            raise RuntimeError("BlockScaler is not fitted")
        Xc = np.asarray(X, dtype=float) - self.mean_
        out = Xc.copy()
        for b, s in self.block_scale_.items():
            out[:, self.blocks == b] /= s
        return out

    def fit_transform(self, X: np.ndarray) -> np.ndarray:
        return self.fit(X).transform(X)


def load_field_matrix(path: str | Path) -> FieldMatrix:
    """Reload a cached field matrix, grid geometry included.

    The grid is stored alongside the matrix on purpose: column *j* only means
    anything in combination with the node layout that produced it, and a model
    reapplied against a differently strided grid would be silently wrong.
    """
    z = np.load(Path(path), allow_pickle=False)
    grid = GridSpec(
        offsets=z["grid_offsets"], index=z["grid_index"],
        npts=int(z["npts"]), spacing=float(z["spacing"]),
        stride=int(z["stride"]), radius=float(z["radius"]),
    )
    accs = [str(a) for a in z["accessions"]]
    centers = {a: z["centers"][i] for i, a in enumerate(accs)}
    return FieldMatrix(
        X=z["X"].astype(np.float64), accessions=accs,
        columns=[str(c) for c in z["columns"]], blocks=z["blocks"],
        block_names=[str(b) for b in z["block_names"]], grid=grid,
        centers=centers,
    )
