"""Minimal fixed-column PDB reader plus the per-atom property tables the field
descriptors need.

Deliberately dependency-free: these assembled dimers are plain ATOM records
written by the assembly code, and pulling in a full structure library would add
a parsing layer that can silently renumber residues. Residue numbering is
load-bearing here -- the box centre is found by residue number.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

# --------------------------------------------------------------- properties
# van der Waals radii (A), Bondi. Used as the CoMSIA steric probe weight.
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "H": 1.20,
       "FE": 2.00, "P": 1.80, "SE": 1.90}

# Kyte-Doolittle hydropathy, per residue.
KD_HYDRO = {
    "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5, "MET": 1.9,
    "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8, "TRP": -0.9,
    "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5, "GLN": -3.5,
    "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5,
}

# Coarse formal charges at pH 8, placed on the ionising atoms.
#
# Deliberately coarse, and the reason is documented: the dominant determinants
# of EV binding in this system are charge, volume and aromaticity (methods
# §13.2), so a formal-charge model captures the leading term. For the accurate
# field, pass a PQR produced by the project's own pdb2pqr/PROPKA pH 8 setup --
# see `charges_from_pqr`.
#
# HIS is left neutral (pKa ~6 << 8). CYS is left neutral too: in the holoenzyme
# these cysteines coordinate the Fe-S clusters and would be thiolate, but the
# clusters are absent from every structure in this series (limitation §8.1), so
# neutral is the self-consistent choice.
FORMAL_CHARGE = {
    ("ASP", "OD1"): -0.5, ("ASP", "OD2"): -0.5,
    ("GLU", "OE1"): -0.5, ("GLU", "OE2"): -0.5,
    ("LYS", "NZ"): 1.0,
    ("ARG", "NE"): 0.2, ("ARG", "NH1"): 0.4, ("ARG", "NH2"): 0.4,
}

AROMATIC_RING = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}

HB_DONOR_O = {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH")}
HB_ACCEPTOR_N = {("HIS", "ND1"), ("HIS", "NE2")}


@dataclass
class Structure:
    """Heavy atoms of one structure."""

    xyz: np.ndarray        # (n, 3) float64
    element: np.ndarray    # (n,) '<U2'
    atom_name: np.ndarray  # (n,) '<U4'
    res_name: np.ndarray   # (n,) '<U3'
    res_seq: np.ndarray    # (n,) int32
    chain: np.ndarray      # (n,) '<U1'
    bfactor: np.ndarray    # (n,) float32  -- pLDDT for predicted models

    def __len__(self) -> int:
        return self.xyz.shape[0]

    def select(self, mask: np.ndarray) -> "Structure":
        return Structure(self.xyz[mask], self.element[mask], self.atom_name[mask],
                         self.res_name[mask], self.res_seq[mask], self.chain[mask],
                         self.bfactor[mask])

    def ca(self, chain: str, res_seq: int) -> np.ndarray | None:
        m = ((self.chain == chain) & (self.res_seq == res_seq)
             & (self.atom_name == "CA"))
        idx = np.flatnonzero(m)
        return self.xyz[idx[0]] if idx.size else None

    def residue_name(self, chain: str, res_seq: int) -> str | None:
        m = (self.chain == chain) & (self.res_seq == res_seq)
        idx = np.flatnonzero(m)
        return str(self.res_name[idx[0]]) if idx.size else None


def read_pdb(path: str | Path, heavy_only: bool = True,
             chains: tuple[str, ...] | None = None) -> Structure:
    """Parse ATOM/HETATM records.

    Both chains are read by default and that matters: the EV site sits on the
    dimer interface, with the opposite protomer's cysteines 10-12 A away
    (structure brief §3). A monomer-only field would describe a pocket with one
    wall missing.
    """
    xyz, el, an, rn, rs, ch, bf = [], [], [], [], [], [], []
    alt_seen: set[tuple] = set()
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            name = line[12:16].strip()
            altloc = line[16].strip()
            resn = line[17:20].strip()
            chain = line[21].strip() or "A"
            if chains is not None and chain not in chains:
                continue
            try:
                seq = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            # element column 77-78; fall back to the atom name's first letter
            elem = line[76:78].strip().upper()
            if not elem:
                elem = "".join(c for c in name if c.isalpha())[:1].upper()
            if heavy_only and elem == "H":
                continue
            if altloc and altloc not in ("A", "1"):
                continue   # keep one altloc only
            key = (chain, seq, name)
            if key in alt_seen:
                continue
            alt_seen.add(key)
            try:
                b = float(line[60:66])
            except ValueError:
                b = 0.0
            xyz.append((x, y, z)); el.append(elem); an.append(name)
            rn.append(resn); rs.append(seq); ch.append(chain); bf.append(b)

    if not xyz:
        raise ValueError(f"no ATOM records parsed from {path}")
    return Structure(
        np.asarray(xyz, dtype=np.float64),
        np.asarray(el, dtype="<U2"), np.asarray(an, dtype="<U4"),
        np.asarray(rn, dtype="<U3"), np.asarray(rs, dtype=np.int32),
        np.asarray(ch, dtype="<U1"), np.asarray(bf, dtype=np.float32),
    )


def charges_from_pqr(path: str | Path) -> np.ndarray:
    """Read the per-atom charge column of a PQR file (pdb2pqr output)."""
    q = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                parts = line.split()
                if len(parts) >= 2:
                    q.append(float(parts[-2]))
    return np.asarray(q, dtype=float)


def formal_charges(st: Structure) -> np.ndarray:
    q = np.zeros(len(st), dtype=float)
    for i in range(len(st)):
        q[i] = FORMAL_CHARGE.get((str(st.res_name[i]), str(st.atom_name[i])), 0.0)
    return q


def steric_weights(st: Structure) -> np.ndarray:
    return np.asarray([VDW.get(str(e), 1.70) for e in st.element], dtype=float)


def hydrophobic_weights(st: Structure) -> np.ndarray:
    """Hydropathy placed on the side-chain beta carbon (CA for glycine)."""
    w = np.zeros(len(st), dtype=float)
    anchor = (st.atom_name == "CB") | ((st.res_name == "GLY") & (st.atom_name == "CA"))
    for i in np.flatnonzero(anchor):
        w[i] = KD_HYDRO.get(str(st.res_name[i]), 0.0)
    return w


def aromatic_centroids(st: Structure) -> tuple[np.ndarray, np.ndarray]:
    """Ring centroids and unit weights.

    Stacking is a ring-face interaction, so one point at the ring centroid is a
    better probe centre than six atom positions -- and it avoids weighting a
    six-membered ring six times more than a five-membered one.
    """
    pts: list[np.ndarray] = []
    for chain in np.unique(st.chain):
        cm = st.chain == chain
        for seq in np.unique(st.res_seq[cm]):
            m = cm & (st.res_seq == seq)
            resn = str(st.res_name[np.flatnonzero(m)[0]])
            ring = AROMATIC_RING.get(resn)
            if ring is None:
                continue
            sel = m & np.isin(st.atom_name, np.asarray(ring, dtype="<U4"))
            if sel.sum() >= 3:
                pts.append(st.xyz[sel].mean(axis=0))
    if not pts:
        return np.zeros((0, 3)), np.zeros(0)
    arr = np.asarray(pts, dtype=float)
    return arr, np.ones(len(arr), dtype=float)


def hbond_weights(st: Structure) -> tuple[np.ndarray, np.ndarray]:
    """(donor, acceptor) indicator weights, coarse by element and residue."""
    donor = np.zeros(len(st), dtype=float)
    acceptor = np.zeros(len(st), dtype=float)
    for i in range(len(st)):
        e, rn, an = str(st.element[i]), str(st.res_name[i]), str(st.atom_name[i])
        if e == "N":
            donor[i] = 1.0
            if (rn, an) in HB_ACCEPTOR_N:
                acceptor[i] = 1.0
        elif e == "O":
            acceptor[i] = 1.0
            if (rn, an) in HB_DONOR_O:
                donor[i] = 1.0
    return donor, acceptor
