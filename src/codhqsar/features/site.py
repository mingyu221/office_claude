"""Binding-site residue descriptors.

Two flavours, both intended as the middle rung between composition (no spatial
information at all) and the full 3D field.

``site_descriptors``   the documented S3 feature: four descriptors -- charge,
                       volume, aromaticity, hydrophobicity -- at K mapped
                       positions. Needs the homolog's residue numbers for those
                       positions, which come from the existing MSA anchored on
                       ``f41_pos``.
``shell_descriptors``  alignment-free. Aggregates the same four properties over
                       concentric shells around the box centre. Nothing to map,
                       so it applies unchanged to de novo designs where there is
                       no alignment to the family to anchor on.

Both exist because the 288-variant problem is a *feature* problem, not a sample
problem: mean-pooled ESM moves almost not at all when 1 of 633 residues changes,
so the variants collapse onto one coordinate. Selecting which variants to keep by
their label would be informative censoring; changing the representation is the
legitimate fix.
"""

from __future__ import annotations

import numpy as np

from ..pdbio import KD_HYDRO, Structure

# Side-chain volumes, A^3 (Zamyatnin).
RES_VOLUME = {
    "GLY": 60.1, "ALA": 88.6, "SER": 89.0, "CYS": 108.5, "ASP": 111.1,
    "PRO": 112.7, "ASN": 114.1, "THR": 116.1, "GLU": 138.4, "VAL": 140.0,
    "GLN": 143.8, "HIS": 153.2, "MET": 162.9, "ILE": 166.7, "LEU": 166.7,
    "LYS": 168.6, "ARG": 173.4, "PHE": 189.9, "TYR": 193.6, "TRP": 227.8,
}
RES_CHARGE = {"ASP": -1.0, "GLU": -1.0, "LYS": 1.0, "ARG": 1.0, "HIS": 0.1}
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}

PROPERTIES = ("charge", "volume", "aromatic", "hydrophobic")


def residue_properties(res_name: str) -> np.ndarray:
    rn = str(res_name).upper()
    return np.array([
        RES_CHARGE.get(rn, 0.0),
        RES_VOLUME.get(rn, 120.0),
        1.0 if rn in AROMATIC_RES else 0.0,
        KD_HYDRO.get(rn, 0.0),
    ], dtype=float)


def site_descriptors(st: Structure, positions: list[int],
                     chain: str = "A") -> np.ndarray:
    """4 descriptors per mapped position, flattened in the given order.

    A position absent from the structure yields zeros, which is distinguishable
    from a real residue because volume is never 0 for one that is present.
    """
    out = []
    for pos in positions:
        rn = st.residue_name(chain, int(pos))
        out.append(np.zeros(len(PROPERTIES)) if rn is None
                   else residue_properties(rn))
    return np.concatenate(out) if out else np.zeros(0)


def site_feature_names(positions: list[int]) -> list[str]:
    return [f"pos{p}_{prop}" for p in positions for prop in PROPERTIES]


def shell_descriptors(st: Structure, center: np.ndarray,
                      edges: tuple[float, ...] = (4.0, 6.0, 8.0, 10.0, 12.0),
                      chains: tuple[str, ...] | None = None) -> np.ndarray:
    """Property sums per concentric shell around the box centre.

    One residue is counted once, at its side-chain anchor (CB, or CA for
    glycine), so a large residue does not get extra weight for having more
    atoms. Both chains count by default: the site is on the dimer interface.
    """
    anchor = (st.atom_name == "CB") | ((st.res_name == "GLY") & (st.atom_name == "CA"))
    if chains is not None:
        anchor &= np.isin(st.chain, np.asarray(chains, dtype="<U1"))
    idx = np.flatnonzero(anchor)
    if idx.size == 0:
        return np.zeros(len(edges) * (len(PROPERTIES) + 1))

    d = np.linalg.norm(st.xyz[idx] - np.asarray(center, float), axis=1)
    props = np.vstack([residue_properties(st.res_name[i]) for i in idx])

    blocks = []
    lo = 0.0
    for hi in edges:
        m = (d >= lo) & (d < hi)
        blocks.append(np.concatenate([props[m].sum(axis=0) if m.any()
                                      else np.zeros(len(PROPERTIES)),
                                      [float(m.sum())]]))
        lo = hi
    return np.concatenate(blocks)


def shell_feature_names(edges: tuple[float, ...] = (4.0, 6.0, 8.0, 10.0, 12.0)) -> list[str]:
    names, lo = [], 0.0
    for hi in edges:
        tag = f"{lo:g}-{hi:g}A"
        names += [f"shell_{tag}_{p}" for p in PROPERTIES] + [f"shell_{tag}_count"]
        lo = hi
    return names


def apply_point_mutation(st: Structure, chain: str, res_seq: int,
                         new_res: str) -> Structure:
    """Relabel one residue in place -- a cheap way to describe a point variant.

    This changes the residue *identity* only; side-chain atoms are untouched, so
    it is valid for identity-based descriptors (site/shell) and NOT for the 3D
    fields, which need real side-chain coordinates. Generating 288 variant field
    sets properly means threading each sequence onto the backbone first.
    """
    m = (st.chain == chain) & (st.res_seq == int(res_seq))
    if not m.any():
        raise ValueError(f"chain {chain} residue {res_seq} not present")
    res_name = st.res_name.copy()
    res_name[m] = str(new_res).upper()[:3]
    return Structure(st.xyz, st.element, st.atom_name, res_name,
                     st.res_seq, st.chain, st.bfactor)
