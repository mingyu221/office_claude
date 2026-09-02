# CODH/EV docking verification

Record of the docking pipeline run against the ChCODH2 A559W variant set
(288 structures) with AutoDock-GPU and RxDock, and the correlation analysis
against C-Docker scores and measured Km / kcat.

- `script/docking_pipeline.py` — the pipeline as it was run, with the fixed
  conditions and the traps that cost time, annotated inline.
- `result/RESULTS_288.md` — results, caveats, and corrections to the earlier
  handoff document.

The data itself lives on the workstation under
`/mnt/af2results/mingyu/workspace/docking_simulation/`; this repository holds
the procedure and the findings, not the structures.
