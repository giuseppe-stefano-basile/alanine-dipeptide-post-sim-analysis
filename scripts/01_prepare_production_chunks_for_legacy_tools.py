#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Split each production dump into pseudo_eq + pseudo_prod chunks (production-only workflow helper)."
    )
    p.add_argument("--npbc-prod", required=True)
    p.add_argument("--pbc-prod", required=True)
    p.add_argument("--workspace", required=True, help="Workspace directory where split dumps are written")
    p.add_argument("--eq-frames", type=int, default=1001, help="Number of initial production frames to write to pseudo_eq")
    return p.parse_args()


def read_one_frame(fin):
    while True:
        line = fin.readline()
        if not line:
            return None, None
        if line.startswith("ITEM: TIMESTEP"):
            break

    frame_lines = [line]
    step_line = fin.readline()
    if not step_line:
        return None, None
    frame_lines.append(step_line)
    step = int(step_line.strip())

    # NUMBER OF ATOMS
    line = fin.readline()
    if not line:
        return None, None
    frame_lines.append(line)

    nat_line = fin.readline()
    if not nat_line:
        return None, None
    frame_lines.append(nat_line)
    nat = int(nat_line.strip())

    # BOX BOUNDS header + 3 rows
    line = fin.readline()
    if not line:
        return None, None
    frame_lines.append(line)
    for _ in range(3):
        x = fin.readline()
        if not x:
            return None, None
        frame_lines.append(x)

    # ATOMS header
    line = fin.readline()
    if not line:
        return None, None
    frame_lines.append(line)

    for _ in range(nat):
        x = fin.readline()
        if not x:
            return None, None
        frame_lines.append(x)

    return step, frame_lines


def split_dump(src: Path, dst_eq: Path, dst_prod: Path, eq_frames: int) -> dict:
    if eq_frames < 1:
        raise ValueError("eq_frames must be >= 1")

    n_total = 0
    n_eq = 0
    n_prod = 0
    first_step = None
    last_step = None

    dst_eq.parent.mkdir(parents=True, exist_ok=True)
    dst_prod.parent.mkdir(parents=True, exist_ok=True)

    with src.open("r") as fin, dst_eq.open("w") as feq, dst_prod.open("w") as fprod:
        while True:
            step, frame_lines = read_one_frame(fin)
            if frame_lines is None:
                break
            n_total += 1
            if first_step is None:
                first_step = step
            last_step = step

            if n_total <= eq_frames:
                feq.writelines(frame_lines)
                n_eq += 1
            else:
                fprod.writelines(frame_lines)
                n_prod += 1

    if n_total == 0:
        raise RuntimeError(f"No frames parsed from {src}")
    if n_prod == 0:
        raise RuntimeError(
            f"eq_frames={eq_frames} consumes all frames for {src}. Reduce --eq-frames."
        )

    return {
        "source": str(src),
        "frames_total": n_total,
        "frames_pseudo_eq": n_eq,
        "frames_pseudo_prod": n_prod,
        "first_step": first_step,
        "last_step": last_step,
        "pseudo_eq_file": str(dst_eq),
        "pseudo_prod_file": str(dst_prod),
    }


def display_path(path: Path, base: Path) -> str:
    absolute = path.expanduser().absolute()
    try:
        return str(absolute.relative_to(base.resolve()))
    except ValueError:
        return str(absolute)


def main() -> None:
    args = parse_args()

    repo_dir = Path(__file__).resolve().parents[1]
    npbc_src = Path(args.npbc_prod).expanduser().absolute()
    pbc_src = Path(args.pbc_prod).expanduser().absolute()
    workspace = Path(args.workspace).expanduser().absolute()
    workspace.mkdir(parents=True, exist_ok=True)

    npbc_eq = workspace / "npbc_eq.dump"
    npbc_prod = workspace / "npbc_prod.dump"
    pbc_eq = workspace / "pbc_eq.dump"
    pbc_prod = workspace / "pbc_prod.dump"

    npbc_summary = split_dump(npbc_src, npbc_eq, npbc_prod, args.eq_frames)
    pbc_summary = split_dump(pbc_src, pbc_eq, pbc_prod, args.eq_frames)
    for item in (npbc_summary, pbc_summary):
        item["source"] = display_path(Path(item["source"]), repo_dir)
        item["pseudo_eq_file"] = display_path(Path(item["pseudo_eq_file"]), repo_dir)
        item["pseudo_prod_file"] = display_path(Path(item["pseudo_prod_file"]), repo_dir)

    summary = {
        "mode": "production_only_split",
        "eq_frames": int(args.eq_frames),
        "npbc": npbc_summary,
        "pbc": pbc_summary,
    }
    (workspace / "split_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
