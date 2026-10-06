#!/usr/bin/env python3
"""Export coefficient contours for PyMOL.

Read against the ChCODH2 self-control structure, whose residue numbering really
is 41, so the maps sit on a backbone you can identify residues on.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np

from codhqsar import config as C
from codhqsar.contours import export_contours, top_nodes
from codhqsar.features.fieldset import load_field_matrix


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", type=Path, default=C.CACHE / "fields_comsia_SELF.npz")
    ap.add_argument("--fit", type=Path, default=C.CACHE / "fit_report")
    ap.add_argument("--out", type=Path, default=C.CACHE / "contours")
    ap.add_argument("--reference", type=Path,
                    default=C.ASSEMBLED / "Q9F8A8_ChCODH2_assembled.pdb")
    ap.add_argument("--level-percentile", type=float, default=95.0)
    ap.add_argument("--top", type=int, default=20)
    args = ap.parse_args()

    fm = load_field_matrix(args.fields)
    coef = np.load(args.fit / "coefficients.npy")
    train_rows = np.load(args.fit / "train_rows.npy")
    X_train = fm.X[train_rows]

    # Contours are drawn in the frame of whichever structure backs them.
    ref_acc = next((a for a in fm.accessions if "ChCODH2" in a or a == "Q9F8A8"),
                   fm.accessions[0])
    center = fm.centers[ref_acc]
    print(f"contours centred on {ref_acc} at {np.round(center, 3)}")

    written = export_contours(
        args.out, coef, X_train, fm.grid, fm.block_names, center,
        reference_pdb=args.reference if args.reference.exists() else None,
        level_percentile=args.level_percentile,
    )
    for name, path in written.items():
        print(f"  {name:<14} {path}")

    table = top_nodes(coef, X_train, fm.grid, fm.block_names, n=args.top)
    table.to_csv(args.out / "top_nodes.tsv", sep="\t", index=False)
    print("\nstrongest contour nodes (offsets in A from the box centre):")
    print(table.head(args.top).round(3).to_string(index=False))
    print("\npositive = more of this property here is predicted to bind EV "
          "more strongly")
    print(f"open with:  pymol {written['pymol']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
