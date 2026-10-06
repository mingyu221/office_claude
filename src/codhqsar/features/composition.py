"""24-dimensional composition baseline.

Not a straw man -- the measured baseline. Train-set cluster cross-validation put
composition at rho_xval 0.4218 +- 0.1471 against ESM-1b layer 2 at 0.4438 +-
0.1383, a gap of one sixth of the error bar, and ESM's top layer 33 came in
*below* composition. ESM layer 0 is a linear transform of composition and was
the best layer on the test set.

So the honest reading is that the predictable part of this label saturates at
composition level. Any new feature family -- the 3D fields included -- has to
beat *this*, not beat ESM, or it has shown nothing.

Composition is also the practical scorer: 24 dimensions, no GPU, no 7.3 GB
checkpoint, seconds for tens of thousands of candidates.
"""

from __future__ import annotations

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA)}

# Net charge at neutral pH, per residue.
CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
AROMATIC = set("FWYH")
HYDROPHOBIC = set("AVLIMFWC")

FEATURE_NAMES = (
    [f"frac_{a}" for a in AA]
    + ["length", "net_charge_per_residue", "aromatic_frac", "hydrophobic_frac"]
)


def composition(seq: str) -> np.ndarray:
    """24-vector for one sequence.

    Non-standard characters (``X``, ``*``) are counted in the length but not in
    any amino-acid bin. They are deliberately not deleted: SeekRank's stock
    pipeline strips them and silently shortens the sequence, which matters for
    metagenomic hits where hundreds of ``X`` occur.
    """
    s = str(seq).strip().upper()
    n = len(s)
    if n == 0:
        raise ValueError("empty sequence")
    counts = np.zeros(len(AA), dtype=float)
    charge = 0.0
    arom = hyd = 0
    for ch in s:
        idx = AA_INDEX.get(ch)
        if idx is not None:
            counts[idx] += 1.0
        charge += CHARGE.get(ch, 0.0)
        arom += ch in AROMATIC
        hyd += ch in HYDROPHOBIC
    return np.concatenate([counts / n, [float(n), charge / n, arom / n, hyd / n]])


def composition_matrix(seqs: dict[str, str]) -> tuple[np.ndarray, list[str]]:
    accs = list(seqs)
    return np.vstack([composition(seqs[a]) for a in accs]), accs


def read_fasta(path) -> dict[str, str]:
    """Read a FASTA keyed by the first whitespace-delimited header token.

    The numeric field in these headers is a *placeholder* holding ``fident``.
    It is not a label. Training on it teaches the model to predict sequence
    identity, so labels are always joined from ``dock_scores.tsv`` by accession
    and never read out of a header.
    """
    seqs: dict[str, str] = {}
    name, buf = None, []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0] if len(line) > 1 else "", []
            elif line:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs
