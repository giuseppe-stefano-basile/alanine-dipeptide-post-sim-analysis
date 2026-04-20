# Full Solute Dynamics Comparison (NPBC eq+prod vs PBC eq+prod)

- Frames NPBC/PBC: **10001 / 10001** (balanced N=10001)
- Grid coarse/fine: **50x50 / 100x100**
- Raw JSD(coarse)/Overlap(coarse): **0.09294 / 0.74243**
- Balanced JSD(coarse) mean [95% CI]: **0.11991 [0.04902, 0.22907]**
- Balanced Overlap(coarse) mean [95% CI]: **0.68143 [0.49758, 0.80829]**

## Analysis settings
- density_method: `hist`
- kde_bandwidth_deg: `18`
- basin_definition: `fes`
- fes_smooth_sigma_bins: `0`
- bootstrap_mode: `block` (block_len_used=1250, iact_est=1286.97)

## FES cropping
- crop_unsampled_fes: `True`
- crop_policy: `either`
- crop_min_prob: `0`
- masked_fraction_raw72: `0.7540`
- masked_fraction_raw144: `0.8552`

## Torsion autocorrelation
- phi IACT NPBC/PBC: 1286.97 / 87.95 frames
- psi IACT NPBC/PBC: 341.86 / 507.80 frames
- phi IACT NPBC/PBC: 643.48 / 43.97 ps
- psi IACT NPBC/PBC: 170.93 / 253.90 ps
- files: `torsion_autocorrelation.csv`, `torsion_autocorrelation_centered.csv`, `19_torsion_autocorrelation.png`, `20_torsion_autocorrelation_centered.png`

## Structural metrics
- RMSD(all) mean NPBC/PBC: 3.101 / 3.113 A
- RMSD(heavy) mean NPBC/PBC: 2.397 / 2.421 A
- RMSF |delta| mean/max: 0.022 / 0.062 A
- Rg(all) mean NPBC/PBC: 2.802 / 2.793 A
- Transition JSD/Overlap: 0.13623 / 0.67624
- Transition events NPBC/PBC: 857 / 932 (see transition_events.csv, state_first_entry.csv)
- Significant basin deltas (95% CI excludes 0): 0 / 8
