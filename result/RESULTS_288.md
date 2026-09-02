# CODH/EV docking verification — 288 variants, 2 engines, 2 ligands (2026-09-02)

Supersedes the n=172 / n=19 numbers in the earlier handoff.

## What was run

288 ChCODH2 A559W-background structures (`input/1su6_varient/`, 15 positions x 19
substitutions + N59N + A559W background + Wildtype), docked with **AutoDock-GPU**
(`autodock_gpu_128wi`, AD4 scoring function) and **RxDock 0.1.0**, each with **two
ligands** that differ only in starting conformer:

| ligand | source | heavy | H | formal charge | rot | inter-ring dihedral |
|---|---|---|---|---|---|---|
| `ev_ds`  | `oxEV_prepared.sd` — the ligand C-Docker used | 16 | 18 | +2 | 3 | −30.0° |
| `ev_lig` | `ev_ref.sd` — crystal reference conformer | 16 | 0 | +2 | 3 | +2.02° (planar) |

Both give a PDBQT of 16 atoms, charge sum +1.998, types `A C N`, TORSDOF 3.

Every stage was validated against an existing artifact before the batch ran:
receptor prep reproduced byte-for-byte; AutoDock reproduced `D63A_A559W.dlg` with
all 50 sorted energies identical; RxDock reproduced `RX_D63A_A559W_out.sd` with
all 50 sorted scores identical.

Runtime: AutoDock 288 grids + 576 dockings ~15 min; RxDock 288 cavities + 576
dockings 6.9 min.

## Results

n = 285 against C-Docker; n = 21 against measured kinetics.

| engine | metric | ligand | rho C-Docker | group-centred | rho Km | rho kcat/Km | SNR | rho ligand-swap |
|---|---|---|---|---|---|---|---|---|
| AutoDock | best | ds  | +0.540 | +0.260 | +0.876 | −0.807 | 6.35 | 0.954 |
| AutoDock | best | lig | +0.521 | +0.242 | +0.882 | −0.822 | 5.74 | 0.954 |
| AutoDock | med  | ds  | +0.550 | +0.262 | +0.856 | −0.801 | 6.62 | 0.938 |
| AutoDock | med  | lig | +0.518 | +0.249 | **+0.903** | **−0.844** | 5.85 | 0.938 |
| RxDock | best | ds  | −0.131 | −0.073 | −0.049 | +0.148 | 0.44 | 0.453 |
| RxDock | best | lig | −0.161 | −0.090 | −0.251 | +0.295 | 0.44 | 0.453 |
| RxDock | med  | ds  | −0.125 | −0.086 | −0.101 | +0.125 | 0.35 | 0.505 |
| RxDock | med  | lig | −0.116 | −0.089 | −0.268 | +0.265 | 0.33 | 0.505 |

Reference points: `cdock` vs Km +0.719; `dchg` (charge change, computable from
sequence alone) vs Km +0.815.

### AutoDock beats C-Docker at predicting measured Km

+0.903 vs +0.719, and above the sequence-only charge descriptor (+0.815). The
n=19 finding in the earlier handoff reproduces at n=21.

### The information is beyond charge for AutoDock, not for C-Docker

Partial (rank-residual) correlation with Km after controlling `dchg`:

| metric | residual rho | p |
|---|---|---|
| `ad_med_lig` | **+0.668** | **0.00093** |
| `ad_best_lig` | +0.551 | 0.0097 |
| `ad_best_ds` | +0.527 | 0.014 |
| `ad_med_ds` | +0.434 | 0.049 |
| `cdock` | +0.173 | 0.45 |

Eight tests in the block; Bonferroni threshold 0.00625. `ad_med_lig` passes.

### The ligand-conformer question is closed

AutoDock ranking is essentially invariant to the starting conformer
(ligand-swap rho 0.94–0.95, mean score shift +0.05 kcal/mol). The earlier
handoff recorded "the inter-ring dihedral is relaxed during search" as an
explicit **assumption**; it is now measured. Which ligand is used does not change
any conclusion — `ds` matches C-Docker, `lig` scores marginally better on Km, and
at n=21 that gap is not discriminating.

### RxDock fails all three checks — but judgement stays deferred

SNR 0.33–0.44 (below 1: between-variant spread is smaller than pose noise), no
correlation with C-Docker after group centring, none with Km, and a ligand-swap
rho of only 0.45–0.51 — swapping an input detail that should not matter
reshuffles half its ranking. `rx_best` has Pearson +0.996 against `rx_best_lig`
while Spearman is +0.453: a few extreme variants carry the linear fit, the bulk
sits in a narrow band whose internal order is noise.

Per the standing policy this is **not** a rejection. RxDock runs alongside
AutoDock through the 350-homolog stage; these are the gate-1 and gate-2 records.
The ligand-swap rho is a useful ceiling: a scale that moves this much under an
irrelevant perturbation cannot carry a larger genuine correlation.

## Caveats

- **C-Docker agreement fell from +0.720 (n=172) to +0.550 (n=286).** The six new
  positions lowered it. By position (ad_best, ds): N59 +0.65, G62 +0.65, P60 +0.64,
  F61 +0.61, L46 +0.45 … I58 −0.18, K66 −0.18, R57 +0.03.
- **AutoDock score correlates +0.87 with `dchg`.** It is close to a charge
  function; the residual signal is real but small next to that.
- **kcat carries no signal** (all metrics n.s.). The signal is in Km only, which
  is what a binding score should predict.
- **The 21 kinetic measurements sit in 4 positions** (R57 8, N59 6, F61 5, G62 2),
  are not a random sample, and the dchg=0 stratum is n=7. The earlier "C-Docker
  sees a steric signal at P60" question cannot be settled here — no P60 variant
  has measured Km.
- **RxDock cavity volume varies 371–1723 A^3** across the 288 (4.6x). Always stored
  as `cav_vol` for QC.

## Corrections to the earlier handoff

1. The execution binary is **AutoDock-GPU** (`autodock_gpu_128wi`); "AutoDock4" is
   the scoring function. Same grids, same `.dlg`.
2. **`--heuristics 0 --autostop 0` are load-bearing.** Without them AutoDock-GPU
   adapts the run count and stops early, breaking the fixed 50-run sample and
   seed reproducibility.
3. Hydrogen counts were recorded backwards: **DS 18, reference 0** (not the
   reverse). H=0 is correct, not a defect — EV is diethyl viologen
   C14H18N2(2+), all 18 hydrogens are non-polar C–H and are merged under the
   united-atom convention.
4. The ligand preparation tool was missing from the record: **Meeko**.
   `ev_ref_meeko.pdbqt` == `ev_lig.pdbqt`. Regenerating the ligand with obabel or
   prepare_ligand4 yields charge sum 0.000 instead of +2 — the numbers still look
   normal, so the loss is easy to miss.
5. Dihedral −178.16° and +2.02° are the same planar state (ring symmetry picks the
   opposite ortho carbon).
6. RxDock electrostatics: no long-range Coulomb term, but ionic interactions do
   exist as a short-range geometric term. `ionic-atom-charges.json` is not a
   contradiction.
7. The 15 A figure means different things per engine — a 30 A cube half-width for
   AutoDock (fixed for all 288), a cavity-detection sphere for RxDock (search
   space then varies per structure).

## Open

- `ev_lig` is plausibly the crystal reference conformer (its centroid is at a
  site coordinate and the file is named `ev_ref.sd`). If so, redocking accuracy
  measured against that crystal is flattered, and the 8X9F 1.99 A / MDW 3.79 A gap
  may have a second cause besides overfitting. **Unverified in this session** —
  test by comparing against the crystal ligand directly, and control with a
  neutral conformer (RDKit ETKDG) docked into both crystals.
- Both engines share `ev_lig` as input, so they share that bias too.

## Outputs

`result/veri_engines/` — `manifest.csv`, `engines_288.csv`, `analysis_288.csv`,
`km_merge_288.csv`, `summary_288.csv`.
Raw data: `prep/ds/*.dlg`, `prep/lig/*.dlg`, `RXDS_*_out.sd`, `RXLG_*_out.sd`
(288 each). The pre-existing `prep/*.dlg` and `RX_*_out.sd` from the n=172 run
were left untouched.
