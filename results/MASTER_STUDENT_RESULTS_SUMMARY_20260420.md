# FES Significance Mini-Report (2026-04-20, No-Smoothing Rerun)

## Method update applied in this rerun

- Density estimator: `hist` (no KDE smoothing)
- Basin mask smoothing: `fes_smooth_sigma_bins=0`
- Crop-aware common-support comparison: `crop_unsampled_fes=true`, `crop_policy=either`, `crop_min_prob=0`
- Bootstrap: block bootstrap (`400` reps in this rerun)
- High-value observables: frame stride `20` (fast rerun mode)

## Cases (production-only)

1. `npbc_corr_vs_pbc_corr_prod_only`
2. `npbc_oldbias_vs_pbc_nvt_anchor_prod_only`

## Global FES significance (95%)

Criterion: significant if `JSD 95% CI lower > 0` AND `Overlap 95% CI upper < 1`.

### 1) npbc_corr_vs_pbc_corr_prod_only

- Coarse (50x50): JSD CI `[0.04902, 0.22907]`, Overlap CI `[0.49758, 0.80829]` -> **Significant**
- Fine (100x100): JSD CI `[0.12262, 0.29182]`, Overlap CI `[0.42928, 0.67410]` -> **Significant**
- Basin-level significant deltas: `0/8`
- Common-support kept fraction:
  - coarse: `1 - 0.7540 = 0.2460`
  - fine: `1 - 0.8552 = 0.1448`

### 2) npbc_oldbias_vs_pbc_nvt_anchor_prod_only

- Coarse (50x50): JSD CI `[0.06313, 0.26001]`, Overlap CI `[0.48246, 0.79112]` -> **Significant**
- Fine (100x100): JSD CI `[0.12745, 0.31274]`, Overlap CI `[0.42551, 0.68534]` -> **Significant**
- Basin-level significant deltas: `0/8`
- Common-support kept fraction:
  - coarse: `1 - 0.7768 = 0.2232`
  - fine: `1 - 0.8805 = 0.1195`

## Rare-basin interpretation (fortuitous-transition question)

- Under strict no-smoothing + common-support, **global FES differences remain significant** for both comparisons.
- However, basin-level deltas are not significant (`0/8` in both cases), and several extremely low-population basins are classified as **possibly fortuitous** in the rare-basin diagnostics.
- Practical readout: yes, NPBC visits to very rare basins can plausibly include fortuitous excursions, while the overall FES difference between lanes is still robust.

See per-case files:
- `results/<case>/03_significance/fes_significance_summary.md`
- `results/<case>/03_significance/rare_basin_diagnostics.md`
