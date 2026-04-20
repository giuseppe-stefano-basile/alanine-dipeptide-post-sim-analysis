#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Diagnose whether NPBC rare-basin visits might be fortuitous by combining "
            "basin populations, transition events, and first-entry timing."
        )
    )
    p.add_argument("--main-compare-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--rare-threshold", type=float, default=0.02, help="Rare basin threshold on NPBC balanced mean population")
    return p.parse_args()


def read_csv_dicts(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open() as f:
        return list(csv.DictReader(f))


def to_float(v: str, default: float = float("nan")) -> float:
    try:
        return float(v)
    except Exception:
        return default


def to_int(v: str, default: int = 0) -> int:
    try:
        return int(float(v))
    except Exception:
        return default


def classify_fortuitous(
    npbc_first_time_ns: float,
    npbc_to_events: int,
    delta_sig: bool,
    npbc_bal_mean: float,
    pbc_bal_mean: float,
    delta_ci_lo: float,
    delta_ci_hi: float,
) -> str:
    # Teaching-oriented heuristic:
    # - supported_difference: statistically significant basin delta
    # - possibly_fortuitous: very small occupancy (both lanes) with non-significant delta
    #   OR late/very-few NPBC entries and non-significant delta
    # - uncertain: everything else
    if delta_sig:
        return "supported_difference"

    both_tiny = (npbc_bal_mean < 0.005) and (pbc_bal_mean < 0.005)
    ci_straddles_zero = (delta_ci_lo <= 0.0) and (delta_ci_hi >= 0.0)
    if both_tiny and ci_straddles_zero:
        return "possibly_fortuitous"

    if (npbc_to_events <= 2) and (npbc_first_time_ns >= 0.7) and ci_straddles_zero:
        return "possibly_fortuitous"

    return "uncertain"


def main() -> None:
    args = parse_args()
    main_dir = Path(args.main_compare_dir).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    basin_rows = read_csv_dicts(main_dir / "basin_populations.csv")
    first_entry = read_csv_dicts(main_dir / "state_first_entry.csv")
    transitions = read_csv_dicts(main_dir / "transition_events.csv")
    metrics_path = main_dir / "metrics_summary.json"
    if not metrics_path.exists():
        raise FileNotFoundError(f"Missing {metrics_path}")
    metrics = json.loads(metrics_path.read_text())

    npbc_first = {}
    pbc_first = {}
    for row in first_entry:
        st = str(row.get("state", ""))
        lane = str(row.get("lane", ""))
        tns = to_float(str(row.get("time_ns", "nan")))
        if lane == "npbc" and st not in npbc_first:
            npbc_first[st] = tns
        if lane == "pbc" and st not in pbc_first:
            pbc_first[st] = tns

    to_counts = defaultdict(lambda: {"npbc": 0, "pbc": 0})
    for row in transitions:
        lane = str(row.get("lane", ""))
        to_state = str(row.get("to_state", ""))
        if to_state == "other":
            continue
        if lane in ("npbc", "pbc"):
            to_counts[to_state][lane] += 1

    rare_rows = []
    for row in basin_rows:
        basin = str(row.get("basin", ""))
        npbc_bal = to_float(str(row.get("npbc_bal_mean", "nan")))
        pbc_bal = to_float(str(row.get("pbc_bal_mean", "nan")))
        d_lo = to_float(str(row.get("delta_bal_ci_lo", "nan")))
        d_hi = to_float(str(row.get("delta_bal_ci_hi", "nan")))
        delta_sig = to_int(str(row.get("delta_significant_95", "0"))) == 1
        if not (npbc_bal < float(args.rare_threshold) or pbc_bal < float(args.rare_threshold)):
            continue
        npbc_t = float(npbc_first.get(basin, float("nan")))
        pbc_t = float(pbc_first.get(basin, float("nan")))
        rare_rows.append(
            {
                "basin": basin,
                "npbc_bal_mean": npbc_bal,
                "pbc_bal_mean": pbc_bal,
                "delta_ci_lo": d_lo,
                "delta_ci_hi": d_hi,
                "delta_significant_95": delta_sig,
                "npbc_first_entry_time_ns": npbc_t,
                "pbc_first_entry_time_ns": pbc_t,
                "npbc_to_events": int(to_counts[basin]["npbc"]),
                "pbc_to_events": int(to_counts[basin]["pbc"]),
            }
        )

    for r in rare_rows:
        r["fortuitous_label"] = classify_fortuitous(
            npbc_first_time_ns=float(r["npbc_first_entry_time_ns"]) if r["npbc_first_entry_time_ns"] == r["npbc_first_entry_time_ns"] else 1.0,
            npbc_to_events=int(r["npbc_to_events"]),
            delta_sig=bool(r["delta_significant_95"]),
            npbc_bal_mean=float(r["npbc_bal_mean"]),
            pbc_bal_mean=float(r["pbc_bal_mean"]),
            delta_ci_lo=float(r["delta_ci_lo"]),
            delta_ci_hi=float(r["delta_ci_hi"]),
        )

    summary = {
        "criterion": {
            "rare_threshold_balanced_population": float(args.rare_threshold),
            "fortuitous_heuristic": (
                "possibly_fortuitous if basin delta is not significant and either (both NPBC/PBC balanced populations <0.005) "
                "or (NPBC first-entry >=0.7 ns and NPBC to-events <=2)"
            ),
        },
        "global_context": {
            "transition_events": metrics.get("transition_events", {}),
            "basin_significance": metrics.get("basin_significance", {}),
        },
        "rare_basins": rare_rows,
        "counts": {
            "n_rare_basins": int(len(rare_rows)),
            "n_possibly_fortuitous": int(sum(1 for r in rare_rows if r["fortuitous_label"] == "possibly_fortuitous")),
            "n_supported_difference": int(sum(1 for r in rare_rows if r["fortuitous_label"] == "supported_difference")),
            "n_uncertain": int(sum(1 for r in rare_rows if r["fortuitous_label"] == "uncertain")),
        },
    }

    (outdir / "rare_basin_diagnostics.json").write_text(json.dumps(summary, indent=2) + "\n")

    md = []
    md.append("# Rare Basin Diagnostics")
    md.append("")
    md.append(f"- Rare threshold (balanced population): `{args.rare_threshold:.4f}`")
    md.append(
        "- Heuristic: possibly fortuitous if basin delta is non-significant and either both lane populations are tiny (<0.005) "
        "or NPBC first-entry is late (>=0.7 ns) with very few entries (<=2)."
    )
    md.append("")
    md.append(f"- Rare basins found: **{summary['counts']['n_rare_basins']}**")
    md.append(f"- Possibly fortuitous: **{summary['counts']['n_possibly_fortuitous']}**")
    md.append(f"- Supported difference: **{summary['counts']['n_supported_difference']}**")
    md.append(f"- Uncertain: **{summary['counts']['n_uncertain']}**")
    md.append("- Interpretation tip: with strict no-smoothing basin masks, transition counts can inflate from jagged boundaries; prioritize CI significance + occupancy magnitude.")
    md.append("")
    if rare_rows:
        md.append("| Basin | NPBC bal | PBC bal | Delta 95% CI | NPBC first (ns) | PBC first (ns) | NPBC to-events | PBC to-events | Label |")
        md.append("|---|---:|---:|---|---:|---:|---:|---:|---|")
        for r in rare_rows:
            md.append(
                "| {b} | {n:.5f} | {p:.5f} | [{lo:.5f}, {hi:.5f}] | {tn:.3f} | {tp:.3f} | {en} | {ep} | {lab} |".format(
                    b=r["basin"],
                    n=float(r["npbc_bal_mean"]),
                    p=float(r["pbc_bal_mean"]),
                    lo=float(r["delta_ci_lo"]),
                    hi=float(r["delta_ci_hi"]),
                    tn=float(r["npbc_first_entry_time_ns"]) if r["npbc_first_entry_time_ns"] == r["npbc_first_entry_time_ns"] else float("nan"),
                    tp=float(r["pbc_first_entry_time_ns"]) if r["pbc_first_entry_time_ns"] == r["pbc_first_entry_time_ns"] else float("nan"),
                    en=int(r["npbc_to_events"]),
                    ep=int(r["pbc_to_events"]),
                    lab=r["fortuitous_label"],
                )
            )
    else:
        md.append("No rare basins under the selected threshold.")
    (outdir / "rare_basin_diagnostics.md").write_text("\n".join(md) + "\n")
    print(str(outdir))


if __name__ == "__main__":
    main()
