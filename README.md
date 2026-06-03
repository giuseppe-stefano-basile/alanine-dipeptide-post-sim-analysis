# Alanine Dipeptide Post-Simulation Analysis Repository

This repository is a student-facing, production-only analysis workflow for OFF23 alanine dipeptide trajectories, including the Leonardo NPT -> NPBC/PBC production workflow.

Goals:
- Use **completed production simulations only**.
- Keep workflow names simple and didactic.
- Run a **comprehensive** analysis stack:
  - main FES/dihedral/structural comparison with bootstrap + basin analysis,
  - high-value observables (hydration shells, H-bonds, radial densities, torsion ACF),
  - explicit FES significance summary,
  - rare-basin diagnostics.

Default methodological policy (teaching + rigor):
- **No FES smoothing** (`density_method=hist`, `fes_smooth_sigma_bins=0`).
- **Common-region FES comparison** via crop-aware masking with `crop_policy=either`
  (only shared sampled support is kept for crop-aware basin and marginal comparisons).

## Repository Layout

- `configs/`: comparison case definitions.
- `data/`: optional local scratch area only; input discovery does not depend on it.
- `environment/`: environment setup for students without `mace_new`.
- `scripts/`: clean pipeline scripts.
- `docs/`: student-oriented didactic guide.
- `reference_core/`: copied core OFF23 analysis scripts used internally.
- `workspaces/`: generated production-only split inputs per case.
- `results/`: final outputs per case.

## Quick Start

```bash
cd /path/to/post_sim_analysis_repo

# create env (no mace_new required)
conda env create -f environment/conda_postsim.yml
conda activate postsim-analysis

# choose where this repo should search for completed trajectories/logs
SEARCH_ROOT=/path/to/npt_bulk_equilibration_workflow

# verify input trajectories
./scripts/00_validate_production_inputs.sh \
  --case leonardo_npt_prod_only \
  --search-root "$SEARCH_ROOT"

# fast classroom demo (lighter bootstrap and subsampled high-value observables)
./scripts/05_run_single_comparison_case.sh \
  --case leonardo_npt_prod_only \
  --search-root "$SEARCH_ROOT" \
  --bootstrap-reps 300 \
  --high-value-frame-stride 10

# run selected configured cases; omit --case only when the search root contains all inputs
./scripts/run_all_comparison_cases.sh \
  --case leonardo_npt_prod_only \
  --search-root "$SEARCH_ROOT" \
  --bootstrap-reps 1000 \
  --high-value-frame-stride 1
```

Student teaching guide:

- `docs/MASTER_STUDENT_TRAJECTORY_ANALYSIS_WORKFLOW.md`
- `DATA_REQUIREMENTS.md`

## Inputs Required

Each case requires:
- `npbc_prod_dump`
- `pbc_prod_dump`
- Optional: corresponding production logs (used for completion checks)

The config stores only file names, relative paths, or glob patterns. The user must choose one or more directories to search with `--search-root`; no local machine path or repo-local trajectory directory is assumed.

Before running the scripts on a new machine, read `DATA_REQUIREMENTS.md` and set `SEARCH_ROOT` to the directory containing the transferred Leonardo production outputs.

All configured cases are mapped in:
- `configs/production_comparison_cases.json`

Included case names:
- `leonardo_npt_prod_only`: reads the expected Leonardo workflow outputs (`npbc_production/traj_alanine_nbpc_prod.dump` and `pbc_production/traj_alanine_pbc_prod.dump`, with alternate `logs/dump.*.lammpstrj` names).
- `npbc_corr_vs_pbc_corr_prod_only`
- `npbc_oldbias_vs_pbc_nvt_anchor_prod_only`

## Outputs Expected (per case)

Under `results/<case_name>/`:
- `01_main_compare/`:
  - `report.md`
  - `metrics_summary.json`
  - `basin_populations*.csv`
  - `*_fes_*.png`, transition plots, torsion ACF, structural plots.
- `02_high_value_observables/`:
  - `high_value_report.md`
  - `high_value_metrics_summary.json`
  - hydration/H-bond/radial density CSV+plots.
- `03_significance/`:
  - `fes_significance_summary.json`
  - `fes_significance_summary.md`
  - `rare_basin_diagnostics.json`
  - `rare_basin_diagnostics.md`

## Production-Only Guarantee

The workflow splits each production dump into two contiguous chunks (`pseudo_eq` + `pseudo_prod`) only to satisfy legacy script interfaces that expect eq+prod inputs.
No equilibration trajectory is used.
The union of chunks equals the full production trajectory exactly once.
