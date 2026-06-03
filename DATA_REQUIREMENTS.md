# Data Requirements

This repository contains one corrected OFF23 tutorial comparison through Git LFS. Final Leonardo production trajectories are not committed to normal Git history and must be transferred separately, then passed to the scripts with `--search-root`.

## Bundled Tutorial Reference

After cloning, fetch the Git LFS tutorial files:

```bash
git lfs install
git lfs pull
```

Bundled files:

```text
data/tutorial_reference/
  npbc_corr_prod.dump
  pbc_corr_prod.dump
  npbc_corr_prod.log
  pbc_corr_prod.log
  CHECKSUMS.sha256
```

Run the corrected tutorial comparison with:

```bash
SEARCH_ROOT=data/tutorial_reference

./scripts/00_validate_production_inputs.sh \
  --case npbc_corr_vs_pbc_corr_prod_only \
  --search-root "$SEARCH_ROOT"

./scripts/05_run_single_comparison_case.sh \
  --case npbc_corr_vs_pbc_corr_prod_only \
  --search-root "$SEARCH_ROOT" \
  --bootstrap-reps 300 \
  --high-value-frame-stride 10
```

## Official Student Case

The final student comparison is:

- case name: `leonardo_npt_prod_only`
- comparison: final Leonardo NPBC production vs final Leonardo PBC production

## Preferred Bundle Layout

Ask for, download, or extract the trajectory bundle with this structure:

```text
trajectory_bundle/
  npbc_production/
    traj_alanine_nbpc_prod.dump
    npbc_production.log
  pbc_production/
    traj_alanine_pbc_prod.dump
    pbc_production.log
```

The logs are useful for completion checks, but the two production dumps are the required files.

## Accepted Alternate Dump Names

The resolver also accepts these Leonardo-style alternatives:

```text
trajectory_bundle/
  npbc_production/logs/dump.npbc_prod.lammpstrj
  pbc_production/logs/dump.pbc_prod.lammpstrj
```

## How To Choose `SEARCH_ROOT`

Set `SEARCH_ROOT` to the directory that contains `npbc_production/` and `pbc_production/`:

```bash
SEARCH_ROOT=/path/to/trajectory_bundle
```

Then validate:

```bash
./scripts/00_validate_production_inputs.sh \
  --case leonardo_npt_prod_only \
  --search-root "$SEARCH_ROOT"
```

Run the fast teaching analysis:

```bash
./scripts/05_run_single_comparison_case.sh \
  --case leonardo_npt_prod_only \
  --search-root "$SEARCH_ROOT" \
  --bootstrap-reps 300 \
  --high-value-frame-stride 10
```

Run the final-quality analysis:

```bash
./scripts/05_run_single_comparison_case.sh \
  --case leonardo_npt_prod_only \
  --search-root "$SEARCH_ROOT" \
  --bootstrap-reps 1000 \
  --high-value-frame-stride 1
```

## Transfer Check

After transfer, it is good practice to record checksums:

```bash
sha256sum npbc_production/traj_alanine_nbpc_prod.dump
sha256sum pbc_production/traj_alanine_pbc_prod.dump
```

Compare these values with the checksums generated on Leonardo when the files were packaged.

## Important

Do not commit new raw `.dump`, `.lammpstrj`, `.dcd`, `.xtc`, `.trr`, `.nc`, or large `.log` files to normal Git history. The only committed dump files should be intentional Git LFS tutorial-reference files. Future Leonardo dumps should remain external unless they are deliberately added through Git LFS or a separate release bundle.
