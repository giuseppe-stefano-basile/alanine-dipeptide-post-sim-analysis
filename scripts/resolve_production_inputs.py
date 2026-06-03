#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from postsim_paths import (
    InputResolutionError,
    load_cases,
    normalize_search_roots,
    resolve_case,
    select_cases,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Resolve configured production trajectory specs against user-selected search roots. "
            "No repository-local trajectory directory is assumed."
        )
    )
    parser.add_argument(
        "--cases-json",
        default=str(Path(__file__).resolve().parents[1] / "configs" / "production_comparison_cases.json"),
        help="Comparison case config file.",
    )
    parser.add_argument(
        "--case",
        action="append",
        default=[],
        help="Case name to resolve. Repeat to resolve multiple cases. Defaults to all cases.",
    )
    parser.add_argument(
        "--search-root",
        action="append",
        default=[],
        help=(
            "Directory chosen by the user to search for trajectories/logs. "
            "Repeat for multiple roots, or pass a colon-separated list."
        ),
    )
    parser.add_argument(
        "--format",
        choices=("json", "tsv", "names"),
        default="json",
        help="Output format for resolved cases.",
    )
    parser.add_argument(
        "--strict-logs",
        action="store_true",
        help="Require optional production logs to resolve.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        cases_file = Path(args.cases_json).expanduser().resolve()
        cases = select_cases(load_cases(cases_file), args.case)
        roots = normalize_search_roots(args.search_root)

        resolved = []
        warnings = []
        for case in cases:
            item, item_warnings = resolve_case(case, roots, strict_logs=args.strict_logs)
            resolved.append(item)
            warnings.extend(item_warnings)

        if args.format == "json":
            print(
                json.dumps(
                    {
                        "cases_json": str(cases_file),
                        "search_roots": [str(root) for root in roots],
                        "cases": resolved,
                        "warnings": warnings,
                    },
                    indent=2,
                )
            )
        elif args.format == "tsv":
            for item in resolved:
                print(
                    "\t".join(
                        [
                            item["name"],
                            item.get("description", ""),
                            item["npbc_prod_dump"],
                            item["pbc_prod_dump"],
                            item.get("npbc_prod_log", ""),
                            item.get("pbc_prod_log", ""),
                        ]
                    )
                )
        else:
            for item in resolved:
                print(item["name"])

        if warnings:
            for warning in warnings:
                print(f"WARNING: {warning}", file=sys.stderr)
        return 0

    except InputResolutionError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
