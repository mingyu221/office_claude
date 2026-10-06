#!/usr/bin/env python3
"""Compute the 3D field matrix for all docking targets and cache it.

~2 s per structure for the CoMSIA fields at the production grid, so roughly
10-15 minutes for 321 single-threaded. The autogrid route instead reads existing
maps; regenerate those with the hand-written GPF if they were deleted after
docking.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from codhqsar import config as C, tables as T
from codhqsar.features.fieldset import build_field_matrix
from codhqsar.features.grid import make_grid


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=["comsia", "autogrid"], default="comsia")
    ap.add_argument("--origin", choices=["SELF", "REF"], default="SELF",
                    help="SELF is the pre-registered primary; REF is sensitivity")
    ap.add_argument("--radius", type=float, default=C.FIELD_RADIUS)
    ap.add_argument("--stride", type=int, default=C.FIELD_STRIDE)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--no-verify-center", action="store_true")
    args = ap.parse_args()

    grid = make_grid(stride=args.stride, radius=args.radius)
    targets = T.build_dataset(require_scores=False)
    print(f"{len(targets)} targets, grid {grid.n_nodes} nodes "
          f"(stride {args.stride}, radius {args.radius} A)")

    t0 = time.time()
    fm = build_field_matrix(
        targets, grid=grid, source=args.source, origin=args.origin,
        verify_center_dev=not args.no_verify_center, progress=True,
    )
    print(f"X = {fm.X.shape} in {time.time() - t0:.1f}s")

    out = args.out or (C.CACHE / f"fields_{args.source}_{args.origin}.npz")
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out, X=fm.X.astype(np.float32),
        accessions=np.asarray(fm.accessions),
        columns=np.asarray(fm.columns), blocks=fm.blocks,
        block_names=np.asarray(fm.block_names),
        grid_index=fm.grid.index, grid_offsets=fm.grid.offsets,
        npts=fm.grid.npts, spacing=fm.grid.spacing,
        stride=fm.grid.stride, radius=fm.grid.radius,
        centers=np.vstack([fm.centers[a] for a in fm.accessions]),
    )
    print(f"cached -> {out}  ({out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
