# Rare Basin Diagnostics

- Rare threshold (balanced population): `0.0200`
- Heuristic: possibly fortuitous if basin delta is non-significant and either both lane populations are tiny (<0.005) or NPBC first-entry is late (>=0.7 ns) with very few entries (<=2).

- Rare basins found: **5**
- Possibly fortuitous: **3**
- Supported difference: **0**
- Uncertain: **2**
- Interpretation tip: with strict no-smoothing basin masks, transition counts can inflate from jagged boundaries; prioritize CI significance + occupancy magnitude.

| Basin | NPBC bal | PBC bal | Delta 95% CI | NPBC first (ns) | PBC first (ns) | NPBC to-events | PBC to-events | Label |
|---|---:|---:|---|---:|---:|---:|---:|---|
| alphaR | 0.02474 | 0.01665 | [-0.01448, 0.03744] | 0.138 | 0.412 | 60 | 39 | uncertain |
| PPII | 0.01896 | 0.02147 | [-0.01255, 0.00791] | 0.024 | 0.035 | 105 | 158 | uncertain |
| alphaL | 0.00063 | 0.00092 | [-0.00141, 0.00061] | 0.176 | 0.417 | 23 | 19 | possibly_fortuitous |
| alphaL_C7ax | 0.00000 | 0.00000 | [0.00000, 0.00000] | 2.126 | nan | 23 | 0 | possibly_fortuitous |
| C7ax | 0.00000 | 0.00000 | [0.00000, 0.00000] | 1.134 | nan | 43 | 0 | possibly_fortuitous |
