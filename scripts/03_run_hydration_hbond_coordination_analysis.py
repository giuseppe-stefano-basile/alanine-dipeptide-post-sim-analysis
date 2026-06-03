#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run high-value observables on production-only split inputs.")
    p.add_argument("--workspace", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--frame-stride", type=int, default=1)
    p.add_argument("--rmax", type=float, default=10.0)
    p.add_argument("--dr", type=float, default=0.25)
    p.add_argument("--coordination-cutoff", type=float, default=3.5)
    p.add_argument("--shell1-cutoff", type=float, default=3.5)
    p.add_argument("--shell2-cutoff", type=float, default=5.0)
    p.add_argument("--hbond-distance-cutoff", type=float, default=2.5)
    p.add_argument("--hbond-angle-cutoff-deg", type=float, default=150.0)
    p.add_argument("--donor-bond-cutoff", type=float, default=1.35)
    p.add_argument("--torsion-max-lag", type=int, default=600)
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
    core = (repo_dir / "reference_core" / "run_stage13corr_high_value_observables.py").resolve()
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
        "--frame-stride", str(args.frame_stride),
        "--rmax", str(args.rmax),
        "--dr", str(args.dr),
        "--coordination-cutoff", str(args.coordination_cutoff),
        "--shell1-cutoff", str(args.shell1_cutoff),
        "--shell2-cutoff", str(args.shell2_cutoff),
        "--hbond-distance-cutoff", str(args.hbond_distance_cutoff),
        "--hbond-angle-cutoff-deg", str(args.hbond_angle_cutoff_deg),
        "--donor-bond-cutoff", str(args.donor_bond_cutoff),
        "--torsion-max-lag", str(args.torsion_max_lag),
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
        "notes": "High-value observables: torsion ACF, hydration-shell occupancy, H-bond statistics, radial density/CN",
    }
    (outdir / "high_value_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(str(outdir))


if __name__ == "__main__":
    main()
