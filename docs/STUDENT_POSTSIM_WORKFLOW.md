# Student Guide: Post-Simulation Analysis (Production-Only)

This guide teaches how to run and interpret **post-simulation analyses** for alanine dipeptide using only completed production trajectories.

Default policy in this repository:
- no FES smoothing (`hist` density, `fes_smooth_sigma_bins=0`)
- crop-aware common-support comparison (`crop_policy=either`)

## 1) What You Need

Required per comparison case:
- one NPBC production dump (`*_prod.dump`)
- one PBC production dump (`*_prod.dump`)
- optional logs (`*.log`) for completion checks

Provided in this repository:
- case mapping: `configs/cases.json`
- trajectory links: `data/trajectories/`

## 2) Environment (No `mace_new` Needed)

```bash
conda env create -f environment/conda_postsim.yml
conda activate postsim-analysis
```

If conda is unavailable:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r environment/requirements.txt
```

## 3) Script Order (Use These Names)

Run from repository root.

1. `./scripts/00_check_inputs.sh`
2. `./scripts/05_run_case_pipeline.sh <case_name> [bootstrap_reps] [high_value_frame_stride]`

For all cases:

```bash
./scripts/run_all_cases.sh [bootstrap_reps] [high_value_frame_stride]
```

## 4) What Each Script Produces

### `00_check_inputs.sh`
Purpose:
- validates that configured production trajectories exist
- prints frame counts, timestep range, stride, and log completion lines

Output:
- terminal summary only (validation report)

### `01_split_production_dump.py`
Purpose:
- splits each production dump into `pseudo_eq` and `pseudo_prod` chunks only to satisfy legacy interfaces
- does **not** use real equilibration trajectories

Output (per case workspace):
- `workspaces/<case>/npbc_eq.dump`
- `workspaces/<case>/npbc_prod.dump`
- `workspaces/<case>/pbc_eq.dump`
- `workspaces/<case>/pbc_prod.dump`
- `workspaces/<case>/split_manifest.json`

### `02_run_main_compare.py`
Purpose (comprehensive FES + basin treatment):
- block-bootstrap uncertainty (time-correlation-aware)
- coarse/fine FES grids (50x50 and 100x100)
- histogram density estimation (no smoothing)
- FES-defined basin masks with `fes_smooth_sigma_bins=0` (no mask smoothing)
- crop-aware FES handling with common-support policy (`either`)

Output:
- `results/<case>/01_main_compare/metrics_summary.json`
- `results/<case>/01_main_compare/report.md`
- FES, marginals, transitions, basin CSV/plots

### `03_run_high_value_observables.py`
Purpose (high-value observables):
- hydration shell occupancy
- H-bond statistics
- radial density and coordination curves
- torsion autocorrelation curves

Output:
- `results/<case>/02_high_value_observables/high_value_metrics_summary.json`
- `results/<case>/02_high_value_observables/high_value_report.md`
- observable CSV and plots

### `04_assess_fes_significance.py`
Purpose:
- turns bootstrap metrics into an explicit significance decision

Criterion used:
- significant if `JSD 95% CI lower > 0` **and** `Overlap 95% CI upper < 1`

Output:
- `results/<case>/03_significance/fes_significance_summary.json`
- `results/<case>/03_significance/fes_significance_summary.md`

### `06_rare_basin_diagnostics.py`
Purpose:
- checks whether rare-basin NPBC visits may be fortuitous
- combines basin populations, first-entry timing, and transition-event counts

Output:
- `results/<case>/03_significance/rare_basin_diagnostics.json`
- `results/<case>/03_significance/rare_basin_diagnostics.md`

## 5) Which Comparisons Are Configured

Current production-only cases in `configs/cases.json`:
- `npbc_corr_vs_pbc_corr_prod_only`
- `npbc_oldbias_vs_pbc_nvt_anchor_prod_only`

These are the target simulation pairs to compare in this training repository.

## 6) Suggested Teaching Workflow

1. Validate inputs with `00_check_inputs.sh`.
2. Run one case with `bootstrap_reps=300` for a classroom demo.
3. Re-run the same case with `bootstrap_reps=1000` for final reporting.
4. Use `high_value_frame_stride=10` for quick teaching demos and `1` for final-quality observables.
5. Ask student to explain consistency between:
   - global FES significance,
   - basin-level significance counts,
   - high-value observables (hydration/H-bond/coordination/ACF).

## 7) Expected Runtime Notes

- main compare (`02`) is the heaviest step (bootstrap + FES computations)
- high-value observables (`03`) can also be long for full trajectories
- prefer long unattended runs for final (`bootstrap_reps=1000`)
