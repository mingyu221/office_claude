"""CoMSIA-style similarity-index fields computed straight from a structure.

Why CoMSIA rather than textbook CoMFA Lennard-Jones/Coulomb fields:

1. **No singularity, so no cutoff.** A 1/r^12 term blows up inside the protein
   and has to be clamped at +-30 kcal/mol. The clamp is an arbitrary constant
   that every descriptor inside the protein then shares, and it makes the field
   discontinuous. A Gaussian stays finite everywhere.
2. **No autogrid dependency.** Fields come from the assembled PDB directly, so
   any candidate -- a remaining IPR010047 sequence, a metagenomic hit, a de novo
   design -- can be scored without running receptor prep and grid computation.
   That is what makes the fitted model reusable rather than a one-off.
3. **More than two properties.** Hydrophobic, H-bond donor/acceptor and an
   aromatic field can be carried alongside steric and electrostatic. For this
   system that matters: EV is a divalent bipyridinium cation, so the aromatic
   stacking term is chemistry, not decoration.

The classic Lennard-Jones/Coulomb fields are available too, read from the
project's own autogrid maps -- see :mod:`codhqsar.features.fields_autogrid`.
Those are physically identical to what produced the labels, so the two together
make a fair comparison rather than a preference.

Field definition, for property *k* at node *j*:

    A_k(j) = sum_i  w_ik * exp(-alpha * r_ij^2)

with ``alpha = 0.3`` (Klebe's default) and distances in Angstrom.
"""

from __future__ import annotations

import warnings

import numpy as np

from .. import config as C
from ..pdbio import (Structure, aromatic_centroids, formal_charges,
                     hbond_weights, hydrophobic_weights, steric_weights)
from .grid import GridSpec

FIELD_NAMES = ("steric", "electrostatic", "hydrophobic", "hbd", "hba", "aromatic")

# exp(-0.3 * 8^2) = 7e-9: beyond this an atom contributes nothing measurable.
_CUTOFF = 8.0
_CHUNK = 1024


# Fields whose probe points are the atoms themselves. Only `aromatic` differs:
# it is evaluated at ring centroids.
_ATOM_FIELDS = ("steric", "electrostatic", "hydrophobic", "hbd", "hba")


def _gaussian_weighted(nodes: np.ndarray, centers: np.ndarray,
                       weights: np.ndarray, alpha: float,
                       cutoff: float = _CUTOFF) -> np.ndarray:
    """``sum_i w_ik exp(-alpha r_ij^2)`` for every property column *k* at once.

    ``weights`` is ``(n_centers, n_properties)``. Sharing one distance matrix
    across properties is what keeps this affordable: the six fields differ only
    in their weights, so recomputing distances per field would cost 6x for
    nothing.

    Dense chunking rather than a sparse distance matrix is deliberate -- a node
    landing exactly on an atom gives distance 0, and a sparse container treats
    that as an absent entry, dropping the single largest contribution.
    """
    n_props = weights.shape[1]
    out = np.zeros((nodes.shape[0], n_props), dtype=np.float64)
    if centers.shape[0] == 0:
        return out

    keep = np.flatnonzero(np.any(weights != 0.0, axis=1))
    if keep.size == 0:
        return out
    centers, weights = centers[keep], weights[keep]

    c2 = cutoff * cutoff
    for start in range(0, nodes.shape[0], _CHUNK):
        chunk = nodes[start:start + _CHUNK]
        d2 = ((chunk[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        contrib = np.where(d2 <= c2, np.exp(-alpha * d2), 0.0)
        out[start:start + chunk.shape[0]] = contrib @ weights
    return out


def compute_fields(
    st: Structure,
    grid: GridSpec,
    center: np.ndarray,
    fields: tuple[str, ...] = FIELD_NAMES,
    alpha: float = C.COMSIA_ALPHA,
    charges: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Similarity indices for one structure.

    Parameters
    ----------
    charges
        Per-atom charges aligned to ``st``. ``None`` uses the coarse formal
        charge model; pass a PQR-derived array for the accurate field.
    """
    unknown = set(fields) - set(FIELD_NAMES)
    if unknown:
        raise ValueError(f"unknown fields {sorted(unknown)}")

    nodes = grid.coords(center)

    # Only atoms that can reach the grid matter.
    reach = grid.radius + _CUTOFF
    keep = np.linalg.norm(st.xyz - np.asarray(center, float), axis=1) <= reach
    near = st.select(keep)
    if len(near) == 0:
        raise ValueError("no atoms within reach of the grid -- wrong box centre?")

    out: dict[str, np.ndarray] = {}

    atom_fields = [f for f in fields if f in _ATOM_FIELDS]
    if atom_fields:
        donor, acceptor = hbond_weights(near)
        if "electrostatic" in atom_fields:
            if charges is None:
                q = formal_charges(near)
            else:
                charges = np.asarray(charges, dtype=float)
                if charges.shape[0] != st.xyz.shape[0]:
                    raise ValueError(
                        f"charges has {charges.shape[0]} entries "
                        f"for {st.xyz.shape[0]} atoms"
                    )
                q = charges[keep]
        else:
            q = None

        table = {
            "steric": lambda: steric_weights(near),
            "electrostatic": lambda: q,
            "hydrophobic": lambda: hydrophobic_weights(near),
            "hbd": lambda: donor,
            "hba": lambda: acceptor,
        }
        W = np.stack([table[f]() for f in atom_fields], axis=1)
        vals = _gaussian_weighted(nodes, near.xyz, W, alpha)
        for i, f in enumerate(atom_fields):
            out[f] = vals[:, i]

    if "aromatic" in fields:
        pts, w = aromatic_centroids(near)
        out["aromatic"] = _gaussian_weighted(
            nodes, pts, w.reshape(-1, 1), alpha
        )[:, 0]

    # An identically zero field is almost always a structure problem, not a
    # real observation: the electrostatic and aromatic fields key on specific
    # side-chain atom names, so a backbone-only or poly-ALA model silently
    # produces a column block of zeros. Those columns then get dropped by the
    # variance filter and the whole run quietly proceeds one field short.
    for name, vals in out.items():
        if not np.any(vals):
            warnings.warn(
                f"{name} field is identically zero -- the structure probably "
                "lacks the side-chain atoms this field reads "
                "(charged/aromatic side chains, or an all-backbone model)",
                RuntimeWarning, stacklevel=2,
            )

    return out


def flatten(fieldmap: dict[str, np.ndarray],
            fields: tuple[str, ...] = FIELD_NAMES) -> np.ndarray:
    """Concatenate fields in a fixed order -- column order must never drift."""
    return np.concatenate([fieldmap[f] for f in fields if f in fieldmap])
