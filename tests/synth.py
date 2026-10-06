"""Synthetic structures for testing the pipeline without the real dataset.

Each fake homolog is a residue cloud around a known box centre with one anchor
residue at ``f41_pos``. The label is generated from a *geometric* property of the
cloud -- acidic side chains crowded toward the centre -- so a working field model
has something real to find through the representation, and a broken one cannot
fake it.

Side-chain atoms are emitted for the charged and aromatic residues. Backbone-only
structures would leave the electrostatic and aromatic fields identically zero,
since those read specific side-chain atom names.
"""

from __future__ import annotations

import numpy as np

RESIDUES = ["ALA", "ASP", "GLU", "LYS", "ARG", "PHE", "TYR", "LEU", "GLY", "SER"]

BACKBONE = (("N", (-0.6, 0.3, 0.0)), ("CA", (0.0, 0.0, 0.0)),
            ("C", (0.7, 0.4, 0.0)), ("O", (1.1, 1.2, 0.0)),
            ("CB", (0.2, -1.2, 0.5)))

# Side-chain atoms that the property tables actually key on.
SIDECHAIN = {
    "ASP": (("CG", (0.4, -2.4, 0.8)), ("OD1", (1.4, -2.9, 1.3)),
            ("OD2", (-0.4, -3.1, 0.3))),
    "GLU": (("CG", (0.4, -2.4, 0.8)), ("CD", (0.6, -3.6, 1.2)),
            ("OE1", (1.6, -4.1, 1.7)), ("OE2", (-0.2, -4.3, 0.7))),
    "LYS": (("CG", (0.4, -2.4, 0.8)), ("CD", (0.6, -3.6, 1.2)),
            ("CE", (0.8, -4.8, 1.6)), ("NZ", (1.0, -6.0, 2.0))),
    "ARG": (("CG", (0.4, -2.4, 0.8)), ("CD", (0.6, -3.6, 1.2)),
            ("NE", (0.8, -4.8, 1.6)), ("CZ", (1.0, -5.6, 2.2)),
            ("NH1", (1.9, -6.4, 2.0)), ("NH2", (0.3, -5.7, 3.2))),
    "PHE": (("CG", (0.4, -2.4, 0.8)), ("CD1", (1.5, -3.0, 1.3)),
            ("CD2", (-0.6, -3.1, 0.6)), ("CE1", (1.6, -4.3, 1.6)),
            ("CE2", (-0.5, -4.4, 0.9)), ("CZ", (0.6, -5.0, 1.4))),
    "TYR": (("CG", (0.4, -2.4, 0.8)), ("CD1", (1.5, -3.0, 1.3)),
            ("CD2", (-0.6, -3.1, 0.6)), ("CE1", (1.6, -4.3, 1.6)),
            ("CE2", (-0.5, -4.4, 0.9)), ("CZ", (0.6, -5.0, 1.4)),
            ("OH", (0.7, -6.3, 1.7))),
    "SER": (("OG", (0.4, -2.2, 0.9)),),
}


def write_fake_dimer(path, seed: int, f41_pos: int = 41, n_res: int = 60,
                     center=(48.884, -6.486, 31.065), charge_bias: float = 0.0):
    """Write a two-chain PDB. ``charge_bias`` pulls acidic residues inward.

    That bias is the planted signal: more negative charge near the box centre
    should mean stronger binding for a dication like ethyl viologen.
    """
    rng = np.random.default_rng(seed)
    ctr = np.asarray(center, float)
    lines, serial = [], 1

    for chain, shift in (("A", np.zeros(3)), ("B", np.array([14.0, 0.0, 0.0]))):
        for i in range(n_res):
            seq = f41_pos + i - n_res // 2
            if chain == "A" and seq == f41_pos:
                res, pos = "PHE", ctr.copy()
            else:
                res = RESIDUES[rng.integers(len(RESIDUES))]
                r = rng.uniform(2.0, 13.0)
                if res in ("ASP", "GLU") and charge_bias:
                    r = max(2.0, r - charge_bias * rng.uniform(0, 6.0))
                v = rng.normal(size=3)
                pos = ctr + shift + r * v / np.linalg.norm(v)

            atoms = [a for a in BACKBONE if not (res == "GLY" and a[0] == "CB")]
            atoms += list(SIDECHAIN.get(res, ()))
            for atom, off in atoms:
                p = pos + np.asarray(off)
                # Fixed-column PDB. The blank after the atom name is altLoc
                # (col 17); dropping it shifts every later field by one and the
                # reader then finds no CA.
                lines.append(
                    f"ATOM  {serial:5d} {atom:<4s} {res:>3s} {chain}{seq:4d}    "
                    f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  1.00 95.00"
                    f"          {atom[0]:>2s}"
                )
                serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return path


def make_series(tmp_path, n_clusters: int = 30, test_every: int = 5):
    """A fake series mirroring the real tier/split structure.

    Expansion members are second members of *training* clusters only. Test
    clusters contribute their representative and nothing else -- that is the
    actual rule, and generating it any other way would plant leakage in the
    fixture rather than in the code under test.
    """
    import pandas as pd

    assembled = tmp_path / "assembled"
    assembled.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in range(n_clusters):
        bias = (c % 5) / 4.0
        is_test = (c % test_every == test_every - 1)
        members = [("core", 0)] if is_test else [("core", 0), ("expansion", 1)]
        for tier, m in members:
            acc = f"FAKE{c:03d}{m}"
            write_fake_dimer(assembled / f"{acc}_dimer.pdb",
                             seed=1000 + 7 * c + m,
                             charge_bias=bias + 0.03 * m)
            rows.append({
                "accession": acc, "cluster": f"C{c:03d}", "f41_pos": 41,
                "tier": tier, "cband": ["40-50", "30-40", "50-60"][c % 3],
                "split": "test" if is_test else "train", "bias": bias,
            })

    df = pd.DataFrame(rows)
    rng = np.random.default_rng(0)
    df["score"] = 4.0 + 1.8 * df["bias"] + rng.normal(0, 0.12, len(df))
    df["ad_med_ds"] = -df["score"]
    df["RMSD"] = rng.uniform(0.5, 2.5, len(df))
    df["center_dev"] = rng.uniform(0, 1.2, len(df))
    df["crit_clash"] = rng.poisson(1.0, len(df))
    df["pLDDT"] = rng.uniform(93, 97, len(df))
    df["struct_source"] = rng.choice(["AFDB", "ColabFold"], len(df), p=[0.9, 0.1])
    df["mode"] = rng.choice(["global", "local"], len(df))
    return df, assembled
