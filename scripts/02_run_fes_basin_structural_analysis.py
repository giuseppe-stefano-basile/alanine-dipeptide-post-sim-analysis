#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Run comprehensive main comparison (FES + basins + structural) on production-only split inputs. "
            "Defaults are no-smoothing + common-support FES comparison."
        )
    )
    p.add_argument("--workspace", required=True, help="Workspace with npbc_eq.dump/npbc_prod.dump/pbc_eq.dump/pbc_prod.dump")
    p.add_argument("--outdir", required=True)
    p.add_argument("--bootstrap-reps", type=int, default=1000)
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--density-method", choices=("hist", "kde"), default="hist")
    p.add_argument("--kde-bandwidth-deg", type=float, default=18.0)
    p.add_argument("--basin-definition", choices=("fixed", "fes"), default="fes")
    p.add_argument("--fes-smooth-sigma-bins", type=float, default=0.0)
    p.add_argument("--crop-policy", choices=("either", "both"), default="either")
    p.add_argument("--crop-min-prob", type=float, default=0.0)
    return p.parse_args()


def display_path(path: Path, base: Path) -> str:
    absolute = path.expanduser().absolute()
    try:
        return str(absolute.relative_to(base.resolve()))
    except ValueError:
        return str(absolute)


def main() -> None:
    args = parse_args()

    script_dir = Path(__file__).resolve().parent
    repo_dir = script_dir.parent
    core = (repo_dir / "reference_core" / "run_stage13corr_full_solute_compare.py").resolve()
    analysis_module = (repo_dir / "reference_core" / "analyze_stage13_dihedrals.py").resolve()
    workspace = Path(args.workspace).resolve()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    for need in ["npbc_eq.dump", "npbc_prod.dump", "pbc_eq.dump", "pbc_prod.dump"]:
        if not (workspace / need).exists():
            raise FileNotFoundError(f"Missing split dump: {workspace / need}")

    cmd = [
        "python3", str(core),
        "--stage-root", str(workspace),
        "--analysis-module", str(analysis_module),
        "--npbc-eq-dump", "npbc_eq.dump",
        "--npbc-prod-dump", "npbc_prod.dump",
        "--pbc-eq-dump", "pbc_eq.dump",
        "--pbc-prod-dump", "pbc_prod.dump",
        "--bootstrap-reps", str(args.bootstrap_reps),
        "--seed", str(args.seed),
        "--rolling-window", "500",
        "--fes-bins-coarse", "50",
        "--fes-bins-fine", "100",
        "--density-method", str(args.density_method),
        "--kde-bandwidth-deg", str(args.kde_bandwidth_deg),
        "--basin-definition", str(args.basin_definition),
        "--fes-smooth-sigma-bins", str(args.fes_smooth_sigma_bins),
        "--bootstrap-mode", "block",
        "--bootstrap-block-len", "0",
        "--crop-unsampled-fes",
        "--crop-policy", str(args.crop_policy),
        "--crop-min-prob", str(args.crop_min_prob),
        "--outdir", str(outdir),
    ]

    subprocess.run(cmd, check=True)

    manifest_paths = {
        str(core): display_path(core, repo_dir),
        str(analysis_module): display_path(analysis_module, repo_dir),
        str(workspace): display_path(workspace, repo_dir),
        str(outdir): display_path(outdir, repo_dir),
    }
    manifest_cmd = [manifest_paths.get(token, token) for token in cmd]
    manifest = {
        "script": display_path(core, repo_dir),
        "workspace": display_path(workspace, repo_dir),
        "outdir": display_path(outdir, repo_dir),
        "command": manifest_cmd,
        "notes": (
            "Comprehensive configuration: no-smoothing histogram density + FES-defined basins "
            "(sigma=0 by default) + block bootstrap + common-support FES crop."
        ),
    }
    (outdir / "main_compare_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(str(outdir))


if __name__ == "__main__":
    main()
