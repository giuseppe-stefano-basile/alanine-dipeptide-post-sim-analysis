# Data Scratch Area

This directory contains the Git LFS-backed corrected tutorial reference data and can also be used as optional local scratch space.

For the bundled tutorial case, use:

```bash
SEARCH_ROOT=data/tutorial_reference
```

For Leonardo or other student runs, place or extract trajectory bundles anywhere convenient and pass that directory with `--search-root`.

Raw trajectory files are ignored by Git unless they are intentionally allowed and tracked through Git LFS.
