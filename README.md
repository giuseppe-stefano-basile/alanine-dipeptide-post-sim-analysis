# Alanine Dipeptide Post-Simulation Analysis Repository

This repository is a student-facing, production-only analysis workflow for OFF23 stage13 trajectories.

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
- `data/`: symlinks to source trajectories and logs.
- `environment/`: environment setup for students without `mace_new`.
- `scripts/`: clean pipeline scripts.
- `docs/`: student-oriented didactic guide.
- `reference_core/`: copied core OFF23 analysis scripts used internally.
- `workspaces/`: generated production-only split inputs per case.
- `results/`: final outputs per case.

## Quick Start

```bash
cd /home/utente/giuseppe/ML_Embedding/MACE/Alanine_dipeptide/post_sim_analysis_repo

# create env (no mace_new required)
conda env create -f environment/conda_postsim.yml
conda activate postsim-analysis

# verify input trajectories
./scripts/00_check_inputs.sh

# run all configured production-only cases
./scripts/run_all_cases.sh 1000 1

# fast classroom demo (lighter bootstrap and subsampled high-value observables)
./scripts/run_all_cases.sh 300 10
```

Student teaching guide:

- `docs/STUDENT_POSTSIM_WORKFLOW.md`

## Inputs Required

Each case requires:
- `npbc_prod_dump`
- `pbc_prod_dump`
- Optional: corresponding production logs (used for completion checks)

All configured defaults are already mapped in:
- `configs/cases.json`

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
