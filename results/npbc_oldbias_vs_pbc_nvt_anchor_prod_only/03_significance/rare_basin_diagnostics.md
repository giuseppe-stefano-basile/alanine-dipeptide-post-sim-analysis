# Rare Basin Diagnostics

- Rare threshold (balanced population): `0.0200`
- Heuristic: possibly fortuitous if basin delta is non-significant and either both lane populations are tiny (<0.005) or NPBC first-entry is late (>=0.7 ns) with very few entries (<=2).

- Rare basins found: **6**
- Possibly fortuitous: **3**
- Supported difference: **0**
- Uncertain: **3**
- Interpretation tip: with strict no-smoothing basin masks, transition counts can inflate from jagged boundaries; prioritize CI significance + occupancy magnitude.

| Basin | NPBC bal | PBC bal | Delta 95% CI | NPBC first (ns) | PBC first (ns) | NPBC to-events | PBC to-events | Label |
|---|---:|---:|---|---:|---:|---:|---:|---|
| alphaR | 0.00937 | 0.01117 | [-0.01164, 0.00947] | 0.065 | 0.073 | 30 | 64 | uncertain |
| alpha_prime | 0.00304 | 0.00544 | [-0.00644, 0.00140] | 0.196 | 0.604 | 27 | 60 | uncertain |
| PPII | 0.02096 | 0.01809 | [-0.00606, 0.01205] | 0.015 | 0.014 | 106 | 129 | uncertain |
| alphaL | 0.00013 | 0.00011 | [-0.00031, 0.00030] | 1.106 | 0.665 | 28 | 3 | possibly_fortuitous |
| alphaL_C7ax | 0.00000 | 0.00000 | [0.00000, 0.00000] | 0.103 | nan | 49 | 0 | possibly_fortuitous |
| C7ax | 0.00190 | 0.00106 | [-0.00167, 0.00430] | 0.066 | 2.454 | 65 | 7 | possibly_fortuitous |
