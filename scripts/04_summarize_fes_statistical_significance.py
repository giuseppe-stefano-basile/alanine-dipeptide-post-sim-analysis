#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Assess FES significance from main comparison outputs.")
    p.add_argument("--main-compare-dir", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def read_significant_basins(path: Path) -> list[dict]:
    if not path.exists():
        return []
    rows = []
    with path.open() as f:
        r = csv.DictReader(f)
        for row in r:
            rows.append(row)
    return rows


def main() -> None:
    args = parse_args()
    main_dir = Path(args.main_compare_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    metrics_file = main_dir / "metrics_summary.json"
    if not metrics_file.exists():
        raise FileNotFoundError(f"Missing {metrics_file}")

    metrics = json.loads(metrics_file.read_text())
    boot = metrics["dihedral_balanced_bootstrap"]

    jsd72_ci = boot["jsd_72_ci"]
    ov72_ci = boot["overlap_72_ci"]
    jsd144_ci = boot["jsd_144_ci"]
    ov144_ci = boot["overlap_144_ci"]

    sig_coarse = (float(jsd72_ci[0]) > 0.0) and (float(ov72_ci[1]) < 1.0)
    sig_fine = (float(jsd144_ci[0]) > 0.0) and (float(ov144_ci[1]) < 1.0)

    basin_sig = metrics.get("basin_significance", {})
    n_sig = int(basin_sig.get("n_significant_95", 0))
    n_total = int(basin_sig.get("n_total", 0))

    sig_rows = read_significant_basins(main_dir / "basin_populations_significant_only.csv")

    summary = {
        "criteria": {
            "global_fes_significant_if": "JSD 95% CI lower bound > 0 and Overlap 95% CI upper bound < 1"
        },
        "coarse_grid": {
            "jsd_mean": float(boot["jsd_72_mean"]),
            "jsd_95ci": [float(jsd72_ci[0]), float(jsd72_ci[1])],
            "overlap_mean": float(boot["overlap_72_mean"]),
            "overlap_95ci": [float(ov72_ci[0]), float(ov72_ci[1])],
            "global_fes_difference_significant_95": bool(sig_coarse),
        },
        "fine_grid": {
            "jsd_mean": float(boot["jsd_144_mean"]),
            "jsd_95ci": [float(jsd144_ci[0]), float(jsd144_ci[1])],
            "overlap_mean": float(boot["overlap_144_mean"]),
            "overlap_95ci": [float(ov144_ci[0]), float(ov144_ci[1])],
            "global_fes_difference_significant_95": bool(sig_fine),
        },
        "basin_population_significance": {
            "n_significant_95": n_sig,
            "n_total": n_total,
            "significant_rows_from_csv": sig_rows,
        },
    }

    (outdir / "fes_significance_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = []
    lines.append("# FES Significance Summary")
    lines.append("")
    lines.append("Criterion used: JSD 95% CI lower bound > 0 AND Overlap 95% CI upper bound < 1.")
    lines.append("")
    lines.append("## Coarse grid")
    lines.append(f"- JSD mean: {summary['coarse_grid']['jsd_mean']:.5f}")
    lines.append(f"- JSD 95% CI: [{summary['coarse_grid']['jsd_95ci'][0]:.5f}, {summary['coarse_grid']['jsd_95ci'][1]:.5f}]")
    lines.append(f"- Overlap mean: {summary['coarse_grid']['overlap_mean']:.5f}")
    lines.append(f"- Overlap 95% CI: [{summary['coarse_grid']['overlap_95ci'][0]:.5f}, {summary['coarse_grid']['overlap_95ci'][1]:.5f}]")
    lines.append(f"- Significant at 95%: {summary['coarse_grid']['global_fes_difference_significant_95']}")
    lines.append("")
    lines.append("## Fine grid")
    lines.append(f"- JSD mean: {summary['fine_grid']['jsd_mean']:.5f}")
    lines.append(f"- JSD 95% CI: [{summary['fine_grid']['jsd_95ci'][0]:.5f}, {summary['fine_grid']['jsd_95ci'][1]:.5f}]")
    lines.append(f"- Overlap mean: {summary['fine_grid']['overlap_mean']:.5f}")
    lines.append(f"- Overlap 95% CI: [{summary['fine_grid']['overlap_95ci'][0]:.5f}, {summary['fine_grid']['overlap_95ci'][1]:.5f}]")
    lines.append(f"- Significant at 95%: {summary['fine_grid']['global_fes_difference_significant_95']}")
    lines.append("")
    lines.append("## Basin-level significance")
    lines.append(f"- Significant basin deltas (95%): {n_sig}/{n_total}")
    lines.append(f"- Significant-basin rows found in CSV: {len(sig_rows)}")

    (outdir / "fes_significance_summary.md").write_text("\n".join(lines) + "\n")
    print(str(outdir))


if __name__ == "__main__":
    main()
