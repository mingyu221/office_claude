"""Fixed constants and paths.

Every value here is pinned by a handoff document. Do not change one without
recording why -- several of them (label sign, grid geometry, noise floor) are
load-bearing for comparability with the 288-variant run and with the docking
labels themselves.
"""

from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------- paths
# The data lives on workstation `epel-af2`. Override with CODHQSAR_BASE.
BASE = Path(os.environ.get(
    "CODHQSAR_BASE",
    "/mnt/af2results/mingyu/workspace/docking_simulation",
))

INPUT = BASE / "input"
SPLIT = BASE / "result" / "homolog_split"
ASSEMBLED = INPUT / "assembled"
MAPS_DIR = Path(os.environ.get("CODHQSAR_MAPS", str(BASE / "scratch" / "maps")))
CACHE = Path(os.environ.get("CODHQSAR_CACHE", str(BASE / "result" / "qsar_cache")))

# Canonical tables (handoff 9.1)
T_MASTER = SPLIT / "master_1827.tsv"
T_TARGETS = SPLIT / "dock_targets.tsv"
T_MANIFEST = SPLIT / "assembly_manifest_all.tsv"
T_SCORES = SPLIT / "dock_scores.tsv"
T_EXPERIMENTS = SPLIT / "experiment_sets.tsv"
T_EXPANSION = SPLIT / "expansion_cap2.tsv"
T_SITE = SPLIT / "site_check.tsv"
T_VARIANTS = BASE / "result" / "Veri_Ca_diameter" / "engines_288_CA_r15.csv"

# ---------------------------------------------------------------- label
# handoff §14: score = -ad_med_ds, higher == stronger binding.
RAW_LABEL = "ad_med_ds"
LABEL = "score"
LABEL_SIGN = -1.0
BETTER = "high"
assert (LABEL_SIGN, BETTER) in ((-1.0, "high"), (1.0, "low"))

# Empirical noise floor of the docking assay, kcal/mol (methods §3).
# Sign-invariant: it is a median |delta|.
NOISE_FLOOR = 0.280

# Reference points measured on ChCODH2 (handoff §6)
WT_AD_MED_DS = -4.38
WT_IQR = 0.09
ASSEMBLY_OFFSET = 0.080          # assembled dimer vs crystal, < WT_IQR
PREP_OFFSET = 0.075              # pipeline systematic, < NOISE_FLOOR

# ---------------------------------------------------------------- grid
# Must match the docking box exactly (docking brief §2) or the field
# descriptors would not describe the volume the labels came from.
NPTS = 80                        # points per axis -> 81 nodes
SPACING = 0.375                  # Angstrom
HALF_EDGE = NPTS * SPACING / 2.0  # 15.000 A
AUTOGRID_RECEPTOR_TYPES = "A C HD N NA OA SA"
AUTOGRID_LIGAND_TYPES = ["A", "C", "N"]

# 1SU6 chain A PHE41 CA -- the origin `center_dev` is measured against.
REF_CENTER = (48.884, -6.486, 31.065)

# Field extraction: a full 81^3 grid is 531441 nodes per map. We only ever
# read a masked, strided subset -- see features/grid.py.
FIELD_RADIUS = 10.0              # A from box centre
FIELD_STRIDE = 2                 # every 2nd node -> 0.75 A effective
ENERGY_CUTOFF = 30.0             # kcal/mol, standard CoMFA clamp

# CoMSIA Gaussian attenuation (Klebe 1994 default)
COMSIA_ALPHA = 0.3

# ---------------------------------------------------------------- misc
SEED = 20260916                  # selection / split / experiment-set seed
CLUSTER_MIN_SEQ_ID = 0.70
ESM_MODEL = "esm1b_t33_650M_UR50S"
ESM_LAYER = 2                    # chosen by train-set cluster CV (methods §13.4)

# Km external validation (methods §10)
KM_N = 19
