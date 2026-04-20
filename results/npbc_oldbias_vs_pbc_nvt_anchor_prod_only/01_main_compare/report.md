# Full Solute Dynamics Comparison (NPBC eq+prod vs PBC eq+prod)

- Frames NPBC/PBC: **10001 / 10001** (balanced N=10001)
- Grid coarse/fine: **50x50 / 100x100**
- Raw JSD(coarse)/Overlap(coarse): **0.12006 / 0.72053**
- Balanced JSD(coarse) mean [95% CI]: **0.14730 [0.06313, 0.26001]**
- Balanced Overlap(coarse) mean [95% CI]: **0.66103 [0.48246, 0.79112]**

## Analysis settings
- density_method: `hist`
- kde_bandwidth_deg: `18`
- basin_definition: `fes`
- fes_smooth_sigma_bins: `0`
- bootstrap_mode: `block` (block_len_used=1250, iact_est=1307.97)

## FES cropping
- crop_unsampled_fes: `True`
- crop_policy: `either`
- crop_min_prob: `0`
- masked_fraction_raw72: `0.7768`
- masked_fraction_raw144: `0.8805`

## Torsion autocorrelation
- phi IACT NPBC/PBC: 1307.97 / 58.24 frames
- psi IACT NPBC/PBC: 121.19 / 432.69 frames
- phi IACT NPBC/PBC: 653.98 / 29.12 ps
- psi IACT NPBC/PBC: 60.59 / 216.35 ps
- files: `torsion_autocorrelation.csv`, `torsion_autocorrelation_centered.csv`, `19_torsion_autocorrelation.png`, `20_torsion_autocorrelation_centered.png`

## Structural metrics
- RMSD(all) mean NPBC/PBC: 3.063 / 3.014 A
- RMSD(heavy) mean NPBC/PBC: 2.363 / 2.301 A
- RMSF |delta| mean/max: 0.079 / 0.179 A
- Rg(all) mean NPBC/PBC: 2.761 / 2.782 A
- Transition JSD/Overlap: 0.14924 / 0.70331
- Transition events NPBC/PBC: 502 / 541 (see transition_events.csv, state_first_entry.csv)
- Significant basin deltas (95% CI excludes 0): 0 / 8
