#!/usr/bin/env python3
"""Comprehensive stage13 alanine NPBC/PBC dihedral + structural analysis.

This script produces:
- Raw and balanced (bootstrap-matched) phi/psi FES comparisons
- Basin and top-state populations with uncertainty
- Time-block stability diagnostics
- Solute structural metrics (RMSD, RMSF, radius of gyration, key distance)
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAVE_PLOT = True
except Exception:
    HAVE_PLOT = False


KB_KJ_MOL_K = 0.00831446261815324
TYPE_TO_ELEMENT = {1: "C", 2: "N", 3: "O", 4: "H"}

# (name, phi_center_deg, psi_center_deg)
BASIN_CENTERS = [
    ("alphaR", -65.0, -40.0),
    ("alpha_prime", -100.0, 0.0),
    ("C7eq", -75.0, 75.0),
    ("PPII", -60.0, 145.0),
    ("C5", -155.0, 160.0),
    ("alphaL", 60.0, 50.0),
    ("alphaL_C7ax", 80.0, -60.0),
    ("C7ax", 50.0, -150.0),
]


@dataclass
class LaneParse:
    source_name: str
    step: np.ndarray
    phi: np.ndarray
    psi: np.ndarray
    solute_coords: np.ndarray  # (nframes, nat_solute, 3)
    solute_atom_ids: np.ndarray
    solute_atom_types: np.ndarray
    box_lengths: np.ndarray  # (nframes, 3)
    boundary_flags: tuple[str, str, str]
    total_frames: int
    frames_with_required_atoms: int
    finite_frames: int
    columns_ok: bool


def parse_args() -> argparse.Namespace:
    root_default = Path(
        "Alanine_dipeptide/MACE-OFF_2023_small/stage13_alanine_prepared/stage13postopt_s1901000_20260316_001735"
    )
    ap = argparse.ArgumentParser(description="Comprehensive stage13 NPBC/PBC alanine analysis.")
    ap.add_argument("--stage13-root", default=str(root_default))
    ap.add_argument("--npbc-eq-dump", default="npbc/traj_alanine_nbpc_off23_stage13post_eq.dump")
    ap.add_argument("--npbc-prod-dump", default="npbc/traj_alanine_nbpc_off23_stage13post_prod.dump")
    ap.add_argument("--pbc-eq-dump", default="pbc/traj_alanine_pbc_off23_stage13post_eq.dump")
    ap.add_argument("--pbc-prod-dump", default="pbc/traj_alanine_pbc_off23_stage13post_prod.dump")
    ap.add_argument("--npbc-mode", choices=("eq", "eq+prod"), default="eq+prod")
    ap.add_argument("--pbc-mode", choices=("eq", "eq+prod"), default="eq")
    ap.add_argument("--atom-ids", default="5,7,9,15,17")
    ap.add_argument("--temp-k", type=float, default=300.0)
    ap.add_argument("--fes-bins-coarse", type=int, default=72)
    ap.add_argument("--fes-bins-fine", type=int, default=144)
    ap.add_argument("--bootstrap-reps", type=int, default=200)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--blocks", type=int, default=10)
    ap.add_argument("--basin-window-deg", type=float, default=20.0)
    ap.add_argument("--top-states", type=int, default=6)
    ap.add_argument("--rolling-window", type=int, default=200)
    ap.add_argument("--outdir", default=None)
    return ap.parse_args()


def dihedral_deg(p0, p1, p2, p3):
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    n1 = np.cross(b0, b1)
    n2 = np.cross(b1, b2)
    n1n = np.linalg.norm(n1)
    n2n = np.linalg.norm(n2)
    b1n = np.linalg.norm(b1)
    if n1n < 1.0e-12 or n2n < 1.0e-12 or b1n < 1.0e-12:
        return float("nan")
    n1 /= n1n
    n2 /= n2n
    m1 = np.cross(n1, b1 / b1n)
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(y, x)))


def wrap_delta_deg(a: np.ndarray, center: float) -> np.ndarray:
    return ((a - center + 180.0) % 360.0) - 180.0


def rolling_mean(x: np.ndarray, window: int) -> np.ndarray:
    y = np.asarray(x, dtype=float)
    n = y.size
    if n == 0:
        return y
    w = max(1, min(int(window), n))
    out = np.full(n, np.nan, dtype=float)
    c = np.cumsum(np.insert(y, 0, 0.0))
    for i in range(n):
        lo = max(0, i - w + 1)
        out[i] = (c[i + 1] - c[lo]) / float(i - lo + 1)
    return out


def write_csv(path: Path, rows: list[tuple], header: list[str]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def quantile_ci(x: np.ndarray, qlo: float = 0.025, qhi: float = 0.975) -> tuple[float, float]:
    a = np.asarray(x, dtype=float)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return float("nan"), float("nan")
    return float(np.quantile(a, qlo)), float(np.quantile(a, qhi))


def block_sem(values: np.ndarray, n_blocks: int) -> float:
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 2:
        return float("nan")
    n_blocks = max(2, min(int(n_blocks), n))
    b = n // n_blocks
    if b <= 0:
        return float("nan")
    means = []
    for i in range(n_blocks):
        lo = i * b
        hi = (i + 1) * b if i < n_blocks - 1 else n
        means.append(float(np.mean(x[lo:hi])))
    if len(means) < 2:
        return float("nan")
    return float(np.std(np.asarray(means), ddof=1) / math.sqrt(len(means)))


def hist2d_prob(phi: np.ndarray, psi: np.ndarray, bins: np.ndarray) -> np.ndarray:
    h, _, _ = np.histogram2d(phi, psi, bins=[bins, bins])
    s = float(np.sum(h))
    if s > 0.0:
        h = h / s
    return h


def hist1d_prob(x: np.ndarray, bins: np.ndarray) -> np.ndarray:
    h, _ = np.histogram(x, bins=bins)
    h = h.astype(float)
    s = float(np.sum(h))
    if s > 0.0:
        h /= s
    return h


def free_energy_from_hist(hist: np.ndarray, temp_k: float) -> np.ndarray:
    p = hist.astype(float)
    total = float(np.sum(p))
    if total > 0.0:
        p = p / total
    eps = 1.0e-15
    f = -KB_KJ_MOL_K * temp_k * np.log(p + eps)
    finite = np.isfinite(f)
    if np.any(finite):
        f = f - np.nanmin(f[finite])
    return f


def js_divergence_from_prob(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float).ravel()
    q = np.asarray(q, dtype=float).ravel()
    p = p / np.sum(p)
    q = q / np.sum(q)
    m = 0.5 * (p + q)
    eps = 1.0e-16
    k1 = np.sum(p * np.log((p + eps) / (m + eps)))
    k2 = np.sum(q * np.log((q + eps) / (m + eps)))
    return float(0.5 * (k1 + k2))


def overlap_coeff_from_prob(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float).ravel()
    q = np.asarray(q, dtype=float).ravel()
    p = p / np.sum(p)
    q = q / np.sum(q)
    return float(np.sum(np.minimum(p, q)))


def ks_distance(a: np.ndarray, b: np.ndarray) -> float:
    aa = np.sort(np.asarray(a, dtype=float))
    bb = np.sort(np.asarray(b, dtype=float))
    if aa.size == 0 or bb.size == 0:
        return float("nan")
    grid = np.unique(np.concatenate([aa, bb]))
    ca = np.searchsorted(aa, grid, side="right") / float(aa.size)
    cb = np.searchsorted(bb, grid, side="right") / float(bb.size)
    return float(np.max(np.abs(ca - cb)))


def basin_population(phi: np.ndarray, psi: np.ndarray, phi_c: float, psi_c: float, halfw: float) -> float:
    dphi = np.abs(wrap_delta_deg(phi, phi_c))
    dpsi = np.abs(wrap_delta_deg(psi, psi_c))
    m = (dphi <= halfw) & (dpsi <= halfw)
    if m.size == 0:
        return float("nan")
    return float(np.mean(m))


def population_in_bin(phi: np.ndarray, psi: np.ndarray, phi_lo: float, phi_hi: float, psi_lo: float, psi_hi: float) -> float:
    m = (phi >= phi_lo) & (phi < phi_hi) & (psi >= psi_lo) & (psi < psi_hi)
    if m.size == 0:
        return float("nan")
    return float(np.mean(m))


def split_blocks(n: int, nb: int) -> list[np.ndarray]:
    if n <= 0:
        return []
    nb = max(1, min(int(nb), n))
    b = n // nb
    out = []
    for i in range(nb):
        lo = i * b
        hi = (i + 1) * b if i < nb - 1 else n
        out.append(np.arange(lo, hi, dtype=int))
    return out


def split_three_blocks(n: int) -> list[tuple[str, np.ndarray]]:
    if n <= 0:
        return []
    idx = np.arange(n, dtype=int)
    cuts = [0, n // 3, (2 * n) // 3, n]
    labels = ["early", "mid", "late"]
    return [(labels[i], idx[cuts[i] : cuts[i + 1]]) for i in range(3)]


def assign_basin(phi: float, psi: float, halfw: float) -> str:
    best = None
    best_d2 = float("inf")
    for (name, pdeg, sdeg) in BASIN_CENTERS:
        dphi = float(abs(((phi - pdeg + 180.0) % 360.0) - 180.0))
        dpsi = float(abs(((psi - sdeg + 180.0) % 360.0) - 180.0))
        if dphi <= halfw and dpsi <= halfw:
            d2 = dphi * dphi + dpsi * dpsi
            if d2 < best_d2:
                best = name
                best_d2 = d2
    return best if best is not None else "other"


def assign_basin_series(phi: np.ndarray, psi: np.ndarray, halfw: float) -> np.ndarray:
    return np.array([assign_basin(float(p), float(s), halfw) for p, s in zip(phi, psi)], dtype=object)


def transition_matrix(states: np.ndarray, state_order: list[str]) -> np.ndarray:
    idx = {s: i for i, s in enumerate(state_order)}
    mat = np.zeros((len(state_order), len(state_order)), dtype=float)
    if states.size < 2:
        return mat
    for a, b in zip(states[:-1], states[1:]):
        ia = idx.get(str(a), idx["other"])
        ib = idx.get(str(b), idx["other"])
        mat[ia, ib] += 1.0
    return mat


def row_normalize(mat: np.ndarray) -> np.ndarray:
    m = np.asarray(mat, dtype=float)
    rs = np.sum(m, axis=1, keepdims=True)
    out = np.divide(m, rs, out=np.zeros_like(m), where=rs > 0.0)
    return out


def dwell_times(states: np.ndarray, target_states: list[str], dt_ps: float) -> dict[str, list[float]]:
    out = {s: [] for s in target_states}
    if states.size == 0:
        return out
    cur = str(states[0])
    run = 1
    for x in states[1:]:
        sx = str(x)
        if sx == cur:
            run += 1
        else:
            if cur in out:
                out[cur].append(run * dt_ps)
            cur = sx
            run = 1
    if cur in out:
        out[cur].append(run * dt_ps)
    return out


def maybe_plot_transition_matrix(out_png: Path, mat_prob: np.ndarray, state_order: list[str], title: str) -> None:
    if not HAVE_PLOT:
        return
    fig, ax = plt.subplots(figsize=(7.2, 6.3), constrained_layout=True)
    im = ax.imshow(mat_prob, origin="lower", aspect="auto", cmap="magma", vmin=0.0, vmax=max(0.05, float(np.nanmax(mat_prob))))
    ax.set_xticks(np.arange(len(state_order)))
    ax.set_yticks(np.arange(len(state_order)))
    ax.set_xticklabels(state_order, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(state_order, fontsize=8)
    ax.set_xlabel("to state")
    ax.set_ylabel("from state")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.9, label="transition probability")
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_dwell_times(out_png: Path, dwell_npbc: dict[str, list[float]], dwell_pbc: dict[str, list[float]], states: list[str]) -> None:
    if not HAVE_PLOT:
        return
    n = len(states)
    fig, axs = plt.subplots(n, 1, figsize=(10.0, max(2.6 * n, 4.5)), constrained_layout=True, sharex=False)
    if n == 1:
        axs = [axs]
    for i, st in enumerate(states):
        ax = axs[i]
        a = np.asarray(dwell_npbc.get(st, []), dtype=float)
        b = np.asarray(dwell_pbc.get(st, []), dtype=float)
        if a.size > 0:
            ax.hist(a, bins=30, density=True, alpha=0.5, label="NPBC")
        if b.size > 0:
            ax.hist(b, bins=30, density=True, alpha=0.5, label="PBC")
        ax.set_title(f"Dwell-time distribution: {st}")
        ax.set_xlabel("dwell time (ps)")
        ax.set_ylabel("density")
        ax.grid(alpha=0.25)
        ax.legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def kabsch_rotation(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    # p, q shape (N,3), centered.
    h = p.T @ q
    u, _s, vt = np.linalg.svd(h)
    r = vt.T @ u.T
    if np.linalg.det(r) < 0.0:
        vt[-1, :] *= -1.0
        r = vt.T @ u.T
    return r


def is_periodic_flag(flag: str) -> bool:
    return str(flag).lower().startswith("p")


def minimum_image_delta(delta: np.ndarray, box_len: np.ndarray, periodic_mask: np.ndarray) -> np.ndarray:
    d = np.asarray(delta, dtype=float).copy()
    for i in range(3):
        if periodic_mask[i] and box_len[i] > 0.0:
            d[i] -= np.round(d[i] / box_len[i]) * box_len[i]
    return d


def reconstruct_backbone_positions_for_dihedrals(
    solute: dict[int, tuple[int, np.ndarray]], box_len: np.ndarray, boundary_flags: tuple[str, str, str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Reconstruct local bonded geometry with minimum-image convention.

    This prevents wrapped PBC coordinates from corrupting phi/psi.
    """
    periodic_mask = np.array([is_periodic_flag(f) for f in boundary_flags], dtype=bool)
    r7 = np.asarray(solute[7][1], dtype=float)
    r5 = r7 + minimum_image_delta(np.asarray(solute[5][1], dtype=float) - r7, box_len, periodic_mask)
    r9 = r7 + minimum_image_delta(np.asarray(solute[9][1], dtype=float) - r7, box_len, periodic_mask)
    r15 = r9 + minimum_image_delta(np.asarray(solute[15][1], dtype=float) - r9, box_len, periodic_mask)
    r17 = r15 + minimum_image_delta(np.asarray(solute[17][1], dtype=float) - r15, box_len, periodic_mask)
    return r5, r7, r9, r15, r17


def parse_dump_lane(path: Path, required_atom_ids: list[int], source_name: str) -> LaneParse:
    required = ["id", "mol", "type", "x", "y", "z"]
    req_set = set(required_atom_ids)
    total_frames = 0
    with_required = 0
    finite = 0
    columns_ok = True

    steps = []
    phis = []
    psis = []
    coords_rows = []
    box_rows = []

    ref_ids = None
    ref_types = None
    boundary_flags = None

    with path.open() as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            total_frames += 1
            step = int(f.readline().strip())
            f.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(f.readline().strip())
            bounds_header = f.readline().strip().split()
            # e.g. ITEM: BOX BOUNDS pp pp pp
            if len(bounds_header) >= 6:
                boundary_flags = (bounds_header[3], bounds_header[4], bounds_header[5])
            else:
                boundary_flags = ("ff", "ff", "ff")
            xlo, xhi = [float(x) for x in f.readline().split()[:2]]
            ylo, yhi = [float(x) for x in f.readline().split()[:2]]
            zlo, zhi = [float(x) for x in f.readline().split()[:2]]
            box_len = np.array([xhi - xlo, yhi - ylo, zhi - zlo], dtype=float)
            hdr = f.readline().strip().split()[2:]
            idx = {k: i for i, k in enumerate(hdr)}
            if any(k not in idx for k in required):
                columns_ok = False
                raise RuntimeError(f"{path} missing required fields {required}")

            solute = {}
            for _ in range(natoms):
                cols = f.readline().split()
                mol = int(cols[idx["mol"]])
                if mol != 1:
                    continue
                aid = int(cols[idx["id"]])
                typ = int(cols[idx["type"]])
                x = float(cols[idx["x"]])
                y = float(cols[idx["y"]])
                z = float(cols[idx["z"]])
                solute[aid] = (typ, np.array([x, y, z], dtype=float))

            if not all(k in solute for k in req_set):
                continue
            with_required += 1

            atom_ids_sorted = np.array(sorted(solute.keys()), dtype=int)
            atom_types_sorted = np.array([solute[aid][0] for aid in atom_ids_sorted], dtype=int)
            atom_coords_sorted = np.array([solute[aid][1] for aid in atom_ids_sorted], dtype=float)

            if ref_ids is None:
                ref_ids = atom_ids_sorted
                ref_types = atom_types_sorted
            else:
                if ref_ids.shape != atom_ids_sorted.shape or np.any(ref_ids != atom_ids_sorted):
                    continue

            d5, d7, d9, d15, d17 = reconstruct_backbone_positions_for_dihedrals(solute, box_len, boundary_flags)
            phi = dihedral_deg(d5, d7, d9, d15)
            psi = dihedral_deg(d7, d9, d15, d17)
            if not (np.isfinite(phi) and np.isfinite(psi)):
                continue
            finite += 1

            steps.append(step)
            phis.append(phi)
            psis.append(psi)
            coords_rows.append(atom_coords_sorted)
            box_rows.append(box_len)

    if ref_ids is None or ref_types is None:
        raise RuntimeError(f"No valid solute frames parsed from {path}")
    if not coords_rows:
        raise RuntimeError(f"No finite phi/psi frames parsed from {path}")
    if boundary_flags is None:
        boundary_flags = ("ff", "ff", "ff")

    return LaneParse(
        source_name=source_name,
        step=np.asarray(steps, dtype=int),
        phi=np.asarray(phis, dtype=float),
        psi=np.asarray(psis, dtype=float),
        solute_coords=np.asarray(coords_rows, dtype=float),
        solute_atom_ids=ref_ids,
        solute_atom_types=ref_types,
        box_lengths=np.asarray(box_rows, dtype=float),
        boundary_flags=boundary_flags,
        total_frames=total_frames,
        frames_with_required_atoms=with_required,
        finite_frames=finite,
        columns_ok=columns_ok,
    )


def unwrap_solute_coords(coords: np.ndarray, box_lengths: np.ndarray, periodic_mask: np.ndarray) -> np.ndarray:
    """Unwrap per-atom coordinates over time using minimum-image continuity."""
    if coords.shape[0] == 0:
        return coords
    out = np.array(coords, copy=True)
    pmask = np.asarray(periodic_mask, dtype=bool)
    for i in range(1, out.shape[0]):
        prev = out[i - 1]
        cur = coords[i]
        delta = cur - coords[i - 1]
        box = np.asarray(box_lengths[i], dtype=float)
        for d in range(3):
            if pmask[d] and box[d] > 0.0:
                delta[:, d] -= np.round(delta[:, d] / box[d]) * box[d]
        out[i] = prev + delta
    return out


def compute_structural_metrics(
    coords: np.ndarray,
    atom_ids: np.ndarray,
    atom_types: np.ndarray,
    rolling_window: int,
) -> dict:
    nframes, natoms, _ = coords.shape
    id_to_idx = {int(a): i for i, a in enumerate(atom_ids.tolist())}
    heavy_idx = np.where(atom_types != 4)[0]
    if heavy_idx.size < 3:
        heavy_idx = np.arange(natoms, dtype=int)

    ref = coords[0]
    ref_center = np.mean(ref[heavy_idx], axis=0)
    ref_fit = ref[heavy_idx] - ref_center
    ref_all = ref - ref_center
    ref_heavy = ref_all[heavy_idx]

    aligned = np.zeros_like(coords)
    rmsd_all = np.zeros(nframes, dtype=float)
    rmsd_heavy = np.zeros(nframes, dtype=float)
    rgyr_all = np.zeros(nframes, dtype=float)
    rgyr_heavy = np.zeros(nframes, dtype=float)
    d_5_17 = np.full(nframes, np.nan, dtype=float)

    idx5 = id_to_idx.get(5, None)
    idx17 = id_to_idx.get(17, None)

    for i in range(nframes):
        cur = coords[i]
        cur_center = np.mean(cur[heavy_idx], axis=0)
        cur_fit = cur[heavy_idx] - cur_center
        rot = kabsch_rotation(cur_fit, ref_fit)
        cur_aligned = (cur - cur_center) @ rot
        aligned[i] = cur_aligned

        d_all = cur_aligned - ref_all
        d_heavy = cur_aligned[heavy_idx] - ref_heavy
        rmsd_all[i] = float(np.sqrt(np.mean(np.sum(d_all * d_all, axis=1))))
        rmsd_heavy[i] = float(np.sqrt(np.mean(np.sum(d_heavy * d_heavy, axis=1))))

        com_all = np.mean(cur, axis=0)
        com_heavy = np.mean(cur[heavy_idx], axis=0)
        rgyr_all[i] = float(np.sqrt(np.mean(np.sum((cur - com_all) ** 2, axis=1))))
        rgyr_heavy[i] = float(np.sqrt(np.mean(np.sum((cur[heavy_idx] - com_heavy) ** 2, axis=1))))

        if idx5 is not None and idx17 is not None:
            d_5_17[i] = float(np.linalg.norm(cur[idx5] - cur[idx17]))

    mean_aligned = np.mean(aligned, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((aligned - mean_aligned) ** 2, axis=2), axis=0))

    return {
        "rmsd_all": rmsd_all,
        "rmsd_heavy": rmsd_heavy,
        "rmsf": rmsf,
        "rgyr_all": rgyr_all,
        "rgyr_heavy": rgyr_heavy,
        "dist_5_17": d_5_17,
        "rmsd_all_roll": rolling_mean(rmsd_all, rolling_window),
        "rmsd_heavy_roll": rolling_mean(rmsd_heavy, rolling_window),
        "rgyr_all_roll": rolling_mean(rgyr_all, rolling_window),
        "rgyr_heavy_roll": rolling_mean(rgyr_heavy, rolling_window),
    }


def build_npbc_time_ns(steps: np.ndarray, sources: np.ndarray) -> np.ndarray:
    # 1 fs timestep => step * 1e-6 ns
    t = np.zeros(steps.size, dtype=float)
    uniq = []
    for s in sources.tolist():
        if not uniq or uniq[-1] != s:
            uniq.append(s)
    offset = 0.0
    for seg in uniq:
        mask = sources == seg
        seg_steps = steps[mask]
        if seg_steps.size == 0:
            continue
        seg_t = (seg_steps - seg_steps.min()) * 1.0e-6
        t[mask] = offset + seg_t
        offset = float(np.max(t[mask]) + 1.0e-6)
    return t


def build_lane_time_ns(steps: np.ndarray, sources: np.ndarray) -> np.ndarray:
    """Build continuous ns timeline across eq/prod segments."""
    return build_npbc_time_ns(steps, sources)


def maybe_plot_fes_triptych(
    out_png: Path,
    f_npbc: np.ndarray,
    f_pbc: np.ndarray,
    title_prefix: str,
) -> None:
    if not HAVE_PLOT:
        return
    delta = f_npbc - f_pbc
    # Force a shared color scale for NPBC/PBC maps so visual comparisons are fair.
    finite_npbc = np.isfinite(f_npbc)
    finite_pbc = np.isfinite(f_pbc)
    pair_vals = np.concatenate((f_npbc[finite_npbc], f_pbc[finite_pbc]))
    if pair_vals.size > 0:
        shared_vmax = float(np.nanmax(pair_vals))
        if not np.isfinite(shared_vmax) or shared_vmax <= 0.0:
            shared_vmax = 1.0
    else:
        shared_vmax = 1.0
    shared_vmin = 0.0

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    ext = [-180, 180, -180, 180]
    im0 = axs[0].imshow(
        f_npbc.T,
        origin="lower",
        extent=ext,
        aspect="auto",
        cmap="viridis",
        vmin=shared_vmin,
        vmax=shared_vmax,
    )
    axs[0].set_title(f"{title_prefix}: NPBC")
    axs[0].set_xlabel("phi (deg)")
    axs[0].set_ylabel("psi (deg)")
    fig.colorbar(im0, ax=axs[0], shrink=0.9)

    im1 = axs[1].imshow(
        f_pbc.T,
        origin="lower",
        extent=ext,
        aspect="auto",
        cmap="viridis",
        vmin=shared_vmin,
        vmax=shared_vmax,
    )
    axs[1].set_title(f"{title_prefix}: PBC")
    axs[1].set_xlabel("phi (deg)")
    axs[1].set_ylabel("psi (deg)")
    fig.colorbar(im1, ax=axs[1], shrink=0.9)

    finite = np.isfinite(delta)
    vmax = float(np.nanmax(np.abs(delta[finite]))) if np.any(finite) else 1.0
    vmax = 1.0 if (not np.isfinite(vmax) or vmax <= 0.0) else vmax
    im2 = axs[2].imshow(
        delta.T, origin="lower", extent=ext, aspect="auto", cmap="coolwarm", vmin=-vmax, vmax=vmax
    )
    axs[2].set_title(f"{title_prefix}: DeltaF NPBC-PBC")
    axs[2].set_xlabel("phi (deg)")
    axs[2].set_ylabel("psi (deg)")
    fig.colorbar(im2, ax=axs[2], shrink=0.9)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_bootstrap_metrics(out_png: Path, jsd_vals: np.ndarray, ov_vals: np.ndarray, jsd_raw: float, ov_raw: float) -> None:
    if not HAVE_PLOT:
        return
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4), constrained_layout=True)
    axs[0].hist(jsd_vals, bins=30, alpha=0.85, color="#1f77b4")
    axs[0].axvline(jsd_raw, color="k", ls="--", lw=1.0, label="raw")
    axs[0].set_title("Bootstrap JSD (balanced)")
    axs[0].set_xlabel("JSD")
    axs[0].set_ylabel("count")
    axs[0].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)

    axs[1].hist(ov_vals, bins=30, alpha=0.85, color="#ff7f0e")
    axs[1].axvline(ov_raw, color="k", ls="--", lw=1.0, label="raw")
    axs[1].set_title("Bootstrap overlap (balanced)")
    axs[1].set_xlabel("overlap coeff")
    axs[1].set_ylabel("count")
    axs[1].grid(alpha=0.25)
    axs[1].legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_basin_populations(out_png: Path, rows: list[dict]) -> None:
    if not HAVE_PLOT or not rows:
        return
    names = [r["basin"] for r in rows]
    x = np.arange(len(names), dtype=float)
    w = 0.27

    pbc = np.array([r["pbc_pop"] for r in rows], dtype=float)
    pbc_sem = np.array([r["pbc_sem"] for r in rows], dtype=float)
    npbc_raw = np.array([r["npbc_raw_pop"] for r in rows], dtype=float)
    npbc_raw_sem = np.array([r["npbc_raw_sem"] for r in rows], dtype=float)
    npbc_bal = np.array([r["npbc_balanced_mean"] for r in rows], dtype=float)
    npbc_bal_lo = np.array([r["npbc_balanced_ci_lo"] for r in rows], dtype=float)
    npbc_bal_hi = np.array([r["npbc_balanced_ci_hi"] for r in rows], dtype=float)
    yerr_bal = np.vstack([npbc_bal - npbc_bal_lo, npbc_bal_hi - npbc_bal])
    yerr_bal = np.clip(np.nan_to_num(yerr_bal, nan=0.0), 0.0, None)

    fig, ax = plt.subplots(figsize=(12, 4.8))
    ax.bar(x - w, pbc, width=w, yerr=pbc_sem, capsize=2, label="PBC raw")
    ax.bar(x, npbc_raw, width=w, yerr=npbc_raw_sem, capsize=2, label="NPBC raw")
    ax.bar(x + w, npbc_bal, width=w, yerr=yerr_bal, capsize=2, label="NPBC balanced")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=22, ha="right")
    ax.set_ylabel("Population")
    ax.set_title("Basin populations")
    ax.grid(alpha=0.25, axis="y")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_timeblock_stability(out_png: Path, rows: list[dict]) -> None:
    if not HAVE_PLOT or not rows:
        return
    lanes = sorted(set(r["lane"] for r in rows))
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.0), constrained_layout=True)
    order = {"early": 0, "mid": 1, "late": 2}
    for lane in lanes:
        rr = sorted([r for r in rows if r["lane"] == lane], key=lambda r: order[r["block"]])
        x = np.arange(len(rr), dtype=float)
        lab = [r["block"] for r in rr]
        jsd_vals = [r["jsd_to_full"] for r in rr]
        ov_vals = [r["overlap_to_full"] for r in rr]
        axs[0].plot(x, jsd_vals, "o-", label=lane)
        axs[1].plot(x, ov_vals, "o-", label=lane)
        axs[0].set_xticks(x)
        axs[1].set_xticks(x)
        axs[0].set_xticklabels(lab)
        axs[1].set_xticklabels(lab)
    axs[0].set_title("JSD to full-lane")
    axs[1].set_title("Overlap to full-lane")
    axs[0].set_ylabel("JSD")
    axs[1].set_ylabel("overlap coeff")
    axs[0].grid(alpha=0.25)
    axs[1].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)
    axs[1].legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_marginals(
    out_png: Path,
    bins: np.ndarray,
    npbc_phi_raw: np.ndarray,
    pbc_phi: np.ndarray,
    npbc_psi_raw: np.ndarray,
    pbc_psi: np.ndarray,
    boot_phi_prob: np.ndarray,
    boot_psi_prob: np.ndarray,
) -> None:
    if not HAVE_PLOT:
        return
    mids = 0.5 * (bins[:-1] + bins[1:])
    phi_raw = hist1d_prob(npbc_phi_raw, bins)
    psi_raw = hist1d_prob(npbc_psi_raw, bins)
    phi_pbc = hist1d_prob(pbc_phi, bins)
    psi_pbc = hist1d_prob(pbc_psi, bins)

    phi_bal_mean = np.mean(boot_phi_prob, axis=0)
    psi_bal_mean = np.mean(boot_psi_prob, axis=0)
    phi_lo = np.quantile(boot_phi_prob, 0.025, axis=0)
    phi_hi = np.quantile(boot_phi_prob, 0.975, axis=0)
    psi_lo = np.quantile(boot_psi_prob, 0.025, axis=0)
    psi_hi = np.quantile(boot_psi_prob, 0.975, axis=0)

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    axs[0].plot(mids, phi_pbc, lw=1.3, label="PBC")
    axs[0].plot(mids, phi_raw, lw=1.3, label="NPBC raw")
    axs[0].plot(mids, phi_bal_mean, lw=1.3, label="NPBC balanced mean")
    axs[0].fill_between(mids, phi_lo, phi_hi, alpha=0.25, label="NPBC balanced 95% CI")
    axs[0].set_title("phi marginal")
    axs[0].set_xlabel("phi (deg)")
    axs[0].set_ylabel("probability")
    axs[0].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)

    axs[1].plot(mids, psi_pbc, lw=1.3, label="PBC")
    axs[1].plot(mids, psi_raw, lw=1.3, label="NPBC raw")
    axs[1].plot(mids, psi_bal_mean, lw=1.3, label="NPBC balanced mean")
    axs[1].fill_between(mids, psi_lo, psi_hi, alpha=0.25, label="NPBC balanced 95% CI")
    axs[1].set_title("psi marginal")
    axs[1].set_xlabel("psi (deg)")
    axs[1].set_ylabel("probability")
    axs[1].grid(alpha=0.25)
    axs[1].legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_fe_profiles_1d(
    out_png: Path,
    out_csv: Path,
    bins: np.ndarray,
    temp_k: float,
    npbc_phi_raw: np.ndarray,
    pbc_phi: np.ndarray,
    npbc_psi_raw: np.ndarray,
    pbc_psi: np.ndarray,
    boot_phi_prob: np.ndarray,
    boot_psi_prob: np.ndarray,
) -> None:
    mids = 0.5 * (bins[:-1] + bins[1:])
    phi_raw = hist1d_prob(npbc_phi_raw, bins)
    psi_raw = hist1d_prob(npbc_psi_raw, bins)
    phi_pbc = hist1d_prob(pbc_phi, bins)
    psi_pbc = hist1d_prob(pbc_psi, bins)

    phi_bal = np.mean(boot_phi_prob, axis=0)
    psi_bal = np.mean(boot_psi_prob, axis=0)

    fe_phi_raw = free_energy_from_hist(phi_raw, temp_k)
    fe_phi_pbc = free_energy_from_hist(phi_pbc, temp_k)
    fe_phi_bal = free_energy_from_hist(phi_bal, temp_k)
    fe_psi_raw = free_energy_from_hist(psi_raw, temp_k)
    fe_psi_pbc = free_energy_from_hist(psi_pbc, temp_k)
    fe_psi_bal = free_energy_from_hist(psi_bal, temp_k)

    fe_phi_boot = np.array([free_energy_from_hist(v, temp_k) for v in boot_phi_prob], dtype=float)
    fe_psi_boot = np.array([free_energy_from_hist(v, temp_k) for v in boot_psi_prob], dtype=float)
    fe_phi_lo = np.quantile(fe_phi_boot, 0.025, axis=0)
    fe_phi_hi = np.quantile(fe_phi_boot, 0.975, axis=0)
    fe_psi_lo = np.quantile(fe_psi_boot, 0.025, axis=0)
    fe_psi_hi = np.quantile(fe_psi_boot, 0.975, axis=0)

    rows = []
    for i in range(mids.size):
        rows.append(
            (
                float(mids[i]),
                float(fe_phi_pbc[i]),
                float(fe_phi_raw[i]),
                float(fe_phi_bal[i]),
                float(fe_phi_lo[i]),
                float(fe_phi_hi[i]),
                float(fe_psi_pbc[i]),
                float(fe_psi_raw[i]),
                float(fe_psi_bal[i]),
                float(fe_psi_lo[i]),
                float(fe_psi_hi[i]),
            )
        )
    write_csv(
        out_csv,
        rows,
        [
            "angle_deg",
            "Fphi_pbc_kJmol",
            "Fphi_npbc_raw_kJmol",
            "Fphi_npbc_bal_mean_kJmol",
            "Fphi_npbc_bal_ci_lo_kJmol",
            "Fphi_npbc_bal_ci_hi_kJmol",
            "Fpsi_pbc_kJmol",
            "Fpsi_npbc_raw_kJmol",
            "Fpsi_npbc_bal_mean_kJmol",
            "Fpsi_npbc_bal_ci_lo_kJmol",
            "Fpsi_npbc_bal_ci_hi_kJmol",
        ],
    )

    if not HAVE_PLOT:
        return

    fig, axs = plt.subplots(1, 2, figsize=(12.2, 4.4), constrained_layout=True)

    axs[0].plot(mids, fe_phi_pbc, lw=1.3, label="PBC")
    axs[0].plot(mids, fe_phi_raw, lw=1.3, label="NPBC raw")
    axs[0].plot(mids, fe_phi_bal, lw=1.3, label="NPBC balanced mean")
    axs[0].fill_between(mids, fe_phi_lo, fe_phi_hi, alpha=0.25, label="NPBC balanced 95% CI")
    axs[0].set_title("1D free-energy profile: phi")
    axs[0].set_xlabel("phi (deg)")
    axs[0].set_ylabel("F (kJ/mol), shifted to min=0")
    axs[0].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)

    axs[1].plot(mids, fe_psi_pbc, lw=1.3, label="PBC")
    axs[1].plot(mids, fe_psi_raw, lw=1.3, label="NPBC raw")
    axs[1].plot(mids, fe_psi_bal, lw=1.3, label="NPBC balanced mean")
    axs[1].fill_between(mids, fe_psi_lo, fe_psi_hi, alpha=0.25, label="NPBC balanced 95% CI")
    axs[1].set_title("1D free-energy profile: psi")
    axs[1].set_xlabel("psi (deg)")
    axs[1].set_ylabel("F (kJ/mol), shifted to min=0")
    axs[1].grid(alpha=0.25)
    axs[1].legend(loc="best", fontsize=8)

    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_bootstrap_zscore_map(out_png: Path, zmap: np.ndarray) -> None:
    if not HAVE_PLOT:
        return
    finite = np.isfinite(zmap)
    vmax = float(np.nanmax(np.abs(zmap[finite]))) if np.any(finite) else 3.0
    vmax = max(3.0, min(vmax, 8.0))
    fig, ax = plt.subplots(figsize=(5.5, 4.6), constrained_layout=True)
    im = ax.imshow(
        zmap.T,
        origin="lower",
        extent=[-180, 180, -180, 180],
        aspect="auto",
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_title("Z-score map: (NPBC_bal - PBC)/sigma_boot")
    ax.set_xlabel("phi (deg)")
    ax.set_ylabel("psi (deg)")
    fig.colorbar(im, ax=ax, shrink=0.9)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_rmsd_timeseries(
    out_png: Path,
    t_npbc: np.ndarray,
    t_pbc: np.ndarray,
    npbc_rmsd_all: np.ndarray,
    npbc_rmsd_heavy: np.ndarray,
    pbc_rmsd_all: np.ndarray,
    pbc_rmsd_heavy: np.ndarray,
    npbc_rmsd_all_roll: np.ndarray,
    npbc_rmsd_heavy_roll: np.ndarray,
    pbc_rmsd_all_roll: np.ndarray,
    pbc_rmsd_heavy_roll: np.ndarray,
) -> None:
    if not HAVE_PLOT:
        return
    fig, axs = plt.subplots(2, 1, figsize=(10.5, 6.5), sharex=False, constrained_layout=True)
    axs[0].plot(t_npbc, npbc_rmsd_all, lw=0.6, alpha=0.35, label="NPBC raw")
    axs[0].plot(t_npbc, npbc_rmsd_all_roll, lw=1.5, label="NPBC rolling")
    axs[0].plot(t_pbc, pbc_rmsd_all, lw=0.6, alpha=0.35, label="PBC raw")
    axs[0].plot(t_pbc, pbc_rmsd_all_roll, lw=1.5, label="PBC rolling")
    axs[0].set_ylabel("RMSD all (A)")
    axs[0].set_title("Solute RMSD vs time")
    axs[0].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)

    axs[1].plot(t_npbc, npbc_rmsd_heavy, lw=0.6, alpha=0.35, label="NPBC raw")
    axs[1].plot(t_npbc, npbc_rmsd_heavy_roll, lw=1.5, label="NPBC rolling")
    axs[1].plot(t_pbc, pbc_rmsd_heavy, lw=0.6, alpha=0.35, label="PBC raw")
    axs[1].plot(t_pbc, pbc_rmsd_heavy_roll, lw=1.5, label="PBC rolling")
    axs[1].set_ylabel("RMSD heavy (A)")
    axs[1].set_xlabel("Time (ns)")
    axs[1].grid(alpha=0.25)
    axs[1].legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_rmsd_distribution(
    out_png: Path,
    npbc_all: np.ndarray,
    pbc_all: np.ndarray,
    npbc_heavy: np.ndarray,
    pbc_heavy: np.ndarray,
) -> None:
    if not HAVE_PLOT:
        return
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    axs[0].hist(npbc_all, bins=50, density=True, alpha=0.5, label="NPBC")
    axs[0].hist(pbc_all, bins=50, density=True, alpha=0.5, label="PBC")
    axs[0].set_title("RMSD all distribution")
    axs[0].set_xlabel("A")
    axs[0].set_ylabel("density")
    axs[0].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)

    axs[1].hist(npbc_heavy, bins=50, density=True, alpha=0.5, label="NPBC")
    axs[1].hist(pbc_heavy, bins=50, density=True, alpha=0.5, label="PBC")
    axs[1].set_title("RMSD heavy distribution")
    axs[1].set_xlabel("A")
    axs[1].set_ylabel("density")
    axs[1].grid(alpha=0.25)
    axs[1].legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_rmsf_per_atom(
    out_png: Path, atom_ids: np.ndarray, atom_types: np.ndarray, npbc_rmsf: np.ndarray, pbc_rmsf: np.ndarray
) -> None:
    if not HAVE_PLOT:
        return
    labels = [f"{aid}:{TYPE_TO_ELEMENT.get(int(t), str(int(t)))}" for aid, t in zip(atom_ids, atom_types)]
    x = np.arange(atom_ids.size, dtype=float)
    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
    axs[0].plot(x, npbc_rmsf, "o-", lw=1.0, ms=3, label="NPBC")
    axs[0].plot(x, pbc_rmsf, "o-", lw=1.0, ms=3, label="PBC")
    axs[0].set_ylabel("RMSF (A)")
    axs[0].set_title("Per-atom RMSF (solute)")
    axs[0].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)

    delta = npbc_rmsf - pbc_rmsf
    axs[1].bar(x, delta, width=0.8)
    axs[1].axhline(0.0, color="k", lw=1.0, ls="--")
    axs[1].set_ylabel("Delta RMSF (NPBC-PBC) (A)")
    axs[1].set_xlabel("Solute atom ID:type")
    axs[1].set_xticks(x)
    axs[1].set_xticklabels(labels, rotation=75, ha="right", fontsize=7)
    axs[1].grid(alpha=0.25, axis="y")
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_rgyr_timeseries(
    out_png: Path,
    t_npbc: np.ndarray,
    t_pbc: np.ndarray,
    npbc_rg_all: np.ndarray,
    npbc_rg_heavy: np.ndarray,
    pbc_rg_all: np.ndarray,
    pbc_rg_heavy: np.ndarray,
    npbc_rg_all_roll: np.ndarray,
    npbc_rg_heavy_roll: np.ndarray,
    pbc_rg_all_roll: np.ndarray,
    pbc_rg_heavy_roll: np.ndarray,
) -> None:
    if not HAVE_PLOT:
        return
    fig, axs = plt.subplots(2, 1, figsize=(10.5, 6.5), sharex=False, constrained_layout=True)
    axs[0].plot(t_npbc, npbc_rg_all, lw=0.6, alpha=0.35, label="NPBC raw")
    axs[0].plot(t_npbc, npbc_rg_all_roll, lw=1.5, label="NPBC rolling")
    axs[0].plot(t_pbc, pbc_rg_all, lw=0.6, alpha=0.35, label="PBC raw")
    axs[0].plot(t_pbc, pbc_rg_all_roll, lw=1.5, label="PBC rolling")
    axs[0].set_ylabel("Rg all (A)")
    axs[0].set_title("Radius of gyration vs time")
    axs[0].grid(alpha=0.25)
    axs[0].legend(loc="best", fontsize=8)

    axs[1].plot(t_npbc, npbc_rg_heavy, lw=0.6, alpha=0.35, label="NPBC raw")
    axs[1].plot(t_npbc, npbc_rg_heavy_roll, lw=1.5, label="NPBC rolling")
    axs[1].plot(t_pbc, pbc_rg_heavy, lw=0.6, alpha=0.35, label="PBC raw")
    axs[1].plot(t_pbc, pbc_rg_heavy_roll, lw=1.5, label="PBC rolling")
    axs[1].set_ylabel("Rg heavy (A)")
    axs[1].set_xlabel("Time (ns)")
    axs[1].grid(alpha=0.25)
    axs[1].legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_distance_timeseries(
    out_png: Path, t_npbc: np.ndarray, t_pbc: np.ndarray, npbc_d: np.ndarray, pbc_d: np.ndarray
) -> None:
    if not HAVE_PLOT:
        return
    fig, ax = plt.subplots(figsize=(10.5, 4.0), constrained_layout=True)
    ax.plot(t_npbc, npbc_d, lw=0.8, alpha=0.8, label="NPBC d(5,17)")
    ax.plot(t_pbc, pbc_d, lw=0.8, alpha=0.8, label="PBC d(5,17)")
    ax.set_title("Backbone distance proxy d(atom5, atom17)")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Distance (A)")
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main() -> None:
    args = parse_args()

    stage13_root = Path(args.stage13_root).resolve()
    npbc_eq_dump = (stage13_root / args.npbc_eq_dump).resolve()
    npbc_prod_dump = (stage13_root / args.npbc_prod_dump).resolve()
    pbc_eq_dump = (stage13_root / args.pbc_eq_dump).resolve()
    pbc_prod_dump = (stage13_root / args.pbc_prod_dump).resolve()
    atom_ids = [int(x.strip()) for x in args.atom_ids.split(",") if x.strip()]
    if atom_ids != [5, 7, 9, 15, 17]:
        raise RuntimeError("This workflow expects atom IDs exactly: 5,7,9,15,17")

    for p in (npbc_eq_dump, pbc_eq_dump):
        if not p.exists():
            raise FileNotFoundError(f"Missing input dump: {p}")
    if args.npbc_mode == "eq+prod" and not npbc_prod_dump.exists():
        raise FileNotFoundError(f"Missing input dump: {npbc_prod_dump}")
    if args.pbc_mode == "eq+prod" and not pbc_prod_dump.exists():
        raise FileNotFoundError(f"Missing input dump: {pbc_prod_dump}")

    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        outdir = (stage13_root / "dihedral_analysis" / stamp).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Parse all lanes
    d_npbc_eq = parse_dump_lane(npbc_eq_dump, atom_ids, "npbc_eq")
    d_npbc_prod = parse_dump_lane(npbc_prod_dump, atom_ids, "npbc_prod") if args.npbc_mode == "eq+prod" else None
    d_pbc_eq = parse_dump_lane(pbc_eq_dump, atom_ids, "pbc_eq")
    d_pbc_prod = parse_dump_lane(pbc_prod_dump, atom_ids, "pbc_prod") if args.pbc_mode == "eq+prod" else None

    # Consistent solute atom ordering/type checks
    if np.any(d_npbc_eq.solute_atom_ids != d_pbc_eq.solute_atom_ids):
        raise RuntimeError("Solute atom ID ordering mismatch across lanes.")
    if np.any(d_npbc_eq.solute_atom_types != d_pbc_eq.solute_atom_types):
        raise RuntimeError("Solute atom type mismatch across lanes.")
    if d_npbc_prod is not None:
        if np.any(d_npbc_eq.solute_atom_ids != d_npbc_prod.solute_atom_ids):
            raise RuntimeError("Solute atom ID ordering mismatch between NPBC eq and NPBC prod.")
        if np.any(d_npbc_eq.solute_atom_types != d_npbc_prod.solute_atom_types):
            raise RuntimeError("Solute atom type mismatch between NPBC eq and NPBC prod.")
    if d_pbc_prod is not None:
        if np.any(d_pbc_eq.solute_atom_ids != d_pbc_prod.solute_atom_ids):
            raise RuntimeError("Solute atom ID ordering mismatch between PBC eq and PBC prod.")
        if np.any(d_pbc_eq.solute_atom_types != d_pbc_prod.solute_atom_types):
            raise RuntimeError("Solute atom type mismatch between PBC eq and PBC prod.")

    # Build combined NPBC series
    if d_npbc_prod is not None:
        npbc_phi = np.concatenate([d_npbc_eq.phi, d_npbc_prod.phi])
        npbc_psi = np.concatenate([d_npbc_eq.psi, d_npbc_prod.psi])
        npbc_step = np.concatenate([d_npbc_eq.step, d_npbc_prod.step])
        npbc_source = np.concatenate(
            [
                np.full(d_npbc_eq.phi.size, "eq", dtype=object),
                np.full(d_npbc_prod.phi.size, "prod", dtype=object),
            ]
        )
        npbc_coords = np.concatenate([d_npbc_eq.solute_coords, d_npbc_prod.solute_coords], axis=0)
        npbc_box = np.concatenate([d_npbc_eq.box_lengths, d_npbc_prod.box_lengths], axis=0)
    else:
        npbc_phi = d_npbc_eq.phi.copy()
        npbc_psi = d_npbc_eq.psi.copy()
        npbc_step = d_npbc_eq.step.copy()
        npbc_source = np.full(d_npbc_eq.phi.size, "eq", dtype=object)
        npbc_coords = d_npbc_eq.solute_coords.copy()
        npbc_box = d_npbc_eq.box_lengths.copy()

    pbc_phi = d_pbc_eq.phi.copy()
    pbc_psi = d_pbc_eq.psi.copy()
    pbc_step = d_pbc_eq.step.copy()
    pbc_source = np.full(d_pbc_eq.phi.size, "eq", dtype=object)
    pbc_coords = d_pbc_eq.solute_coords.copy()
    pbc_box = d_pbc_eq.box_lengths.copy()
    if d_pbc_prod is not None:
        pbc_phi = np.concatenate([d_pbc_eq.phi, d_pbc_prod.phi])
        pbc_psi = np.concatenate([d_pbc_eq.psi, d_pbc_prod.psi])
        pbc_step = np.concatenate([d_pbc_eq.step, d_pbc_prod.step])
        pbc_source = np.concatenate(
            [
                np.full(d_pbc_eq.phi.size, "eq", dtype=object),
                np.full(d_pbc_prod.phi.size, "prod", dtype=object),
            ]
        )
        pbc_coords = np.concatenate([d_pbc_eq.solute_coords, d_pbc_prod.solute_coords], axis=0)
        pbc_box = np.concatenate([d_pbc_eq.box_lengths, d_pbc_prod.box_lengths], axis=0)

    if npbc_phi.size == 0 or pbc_phi.size == 0:
        raise RuntimeError("No finite phi/psi frames parsed.")
    if npbc_phi.size < pbc_phi.size:
        raise RuntimeError(
            f"NPBC finite frames ({npbc_phi.size}) < PBC finite frames ({pbc_phi.size}); cannot balanced-subsample."
        )

    # Input and sanity checks
    finite_frac_npbc = float(np.mean(np.isfinite(npbc_phi) & np.isfinite(npbc_psi)))
    finite_frac_pbc = float(np.mean(np.isfinite(pbc_phi) & np.isfinite(pbc_psi)))
    if finite_frac_npbc < 0.99 or finite_frac_pbc < 0.99:
        raise RuntimeError(
            f"Finite fraction check failed: npbc={finite_frac_npbc:.5f}, pbc={finite_frac_pbc:.5f}"
        )
    for arr, name in ((npbc_phi, "npbc_phi"), (npbc_psi, "npbc_psi"), (pbc_phi, "pbc_phi"), (pbc_psi, "pbc_psi")):
        if np.min(arr) < -180.01 or np.max(arr) > 180.01:
            raise RuntimeError(f"{name} outside [-180,180] bounds.")

    # Save timeseries tables
    write_csv(
        outdir / "npbc_phi_psi_timeseries.csv",
        [(int(i), str(npbc_source[i]), int(npbc_step[i]), float(npbc_phi[i]), float(npbc_psi[i])) for i in range(npbc_phi.size)],
        ["frame_index", "source", "step", "phi_deg", "psi_deg"],
    )
    write_csv(
        outdir / "pbc_phi_psi_timeseries.csv",
        [(int(i), int(pbc_step[i]), float(pbc_phi[i]), float(pbc_psi[i])) for i in range(pbc_phi.size)],
        ["frame_index", "step", "phi_deg", "psi_deg"],
    )

    # Raw FES track
    bins72 = np.linspace(-180.0, 180.0, args.fes_bins_coarse + 1)
    bins144 = np.linspace(-180.0, 180.0, args.fes_bins_fine + 1)
    p_npbc_72 = hist2d_prob(npbc_phi, npbc_psi, bins72)
    p_pbc_72 = hist2d_prob(pbc_phi, pbc_psi, bins72)
    p_npbc_144 = hist2d_prob(npbc_phi, npbc_psi, bins144)
    p_pbc_144 = hist2d_prob(pbc_phi, pbc_psi, bins144)

    f_npbc_72 = free_energy_from_hist(p_npbc_72, args.temp_k)
    f_pbc_72 = free_energy_from_hist(p_pbc_72, args.temp_k)
    f_npbc_144 = free_energy_from_hist(p_npbc_144, args.temp_k)
    f_pbc_144 = free_energy_from_hist(p_pbc_144, args.temp_k)

    raw_metrics = {
        "raw_jsd_72": js_divergence_from_prob(p_npbc_72, p_pbc_72),
        "raw_overlap_72": overlap_coeff_from_prob(p_npbc_72, p_pbc_72),
        "raw_jsd_144": js_divergence_from_prob(p_npbc_144, p_pbc_144),
        "raw_overlap_144": overlap_coeff_from_prob(p_npbc_144, p_pbc_144),
    }

    # Top states (coarse grid)
    comb = p_npbc_72 + p_pbc_72
    order = np.argsort(comb.ravel())[::-1]
    top_states = []
    for idx_flat in order:
        if len(top_states) >= args.top_states:
            break
        i, j = np.unravel_index(idx_flat, comb.shape)
        if comb[i, j] <= 0.0:
            continue
        top_states.append(
            (f"S{len(top_states)+1}", float(bins72[i]), float(bins72[i + 1]), float(bins72[j]), float(bins72[j + 1]))
        )

    # Balanced bootstrap track
    rng = np.random.default_rng(args.seed)
    n_target = pbc_phi.size
    reps = int(args.bootstrap_reps)
    boot_jsd72 = np.zeros(reps, dtype=float)
    boot_ov72 = np.zeros(reps, dtype=float)
    boot_jsd144 = np.zeros(reps, dtype=float)
    boot_ov144 = np.zeros(reps, dtype=float)
    boot_phi_prob = np.zeros((reps, args.fes_bins_coarse), dtype=float)
    boot_psi_prob = np.zeros((reps, args.fes_bins_coarse), dtype=float)
    acc_72 = np.zeros_like(p_npbc_72)
    acc_144 = np.zeros_like(p_npbc_144)
    acc_72_sq = np.zeros_like(p_npbc_72)

    basin_boot = {name: np.zeros(reps, dtype=float) for (name, _, _) in BASIN_CENTERS}
    state_boot = {name: np.zeros(reps, dtype=float) for (name, *_r) in top_states}

    for ib in range(reps):
        idx = rng.choice(npbc_phi.size, size=n_target, replace=False)
        phi_s = npbc_phi[idx]
        psi_s = npbc_psi[idx]

        p_s_72 = hist2d_prob(phi_s, psi_s, bins72)
        p_s_144 = hist2d_prob(phi_s, psi_s, bins144)
        acc_72 += p_s_72
        acc_144 += p_s_144
        acc_72_sq += p_s_72 * p_s_72
        boot_jsd72[ib] = js_divergence_from_prob(p_s_72, p_pbc_72)
        boot_ov72[ib] = overlap_coeff_from_prob(p_s_72, p_pbc_72)
        boot_jsd144[ib] = js_divergence_from_prob(p_s_144, p_pbc_144)
        boot_ov144[ib] = overlap_coeff_from_prob(p_s_144, p_pbc_144)

        boot_phi_prob[ib] = hist1d_prob(phi_s, bins72)
        boot_psi_prob[ib] = hist1d_prob(psi_s, bins72)

        for (name, pdeg, sdeg) in BASIN_CENTERS:
            basin_boot[name][ib] = basin_population(phi_s, psi_s, pdeg, sdeg, args.basin_window_deg)
        for (name, plo, phi, slo, shi) in top_states:
            state_boot[name][ib] = population_in_bin(phi_s, psi_s, plo, phi, slo, shi)

    p_npbc_bal_72 = acc_72 / float(reps)
    p_npbc_bal_144 = acc_144 / float(reps)
    p_npbc_bal_72_std = np.sqrt(np.maximum(acc_72_sq / float(reps) - p_npbc_bal_72 * p_npbc_bal_72, 0.0))
    zmap_72 = (p_npbc_bal_72 - p_pbc_72) / np.maximum(p_npbc_bal_72_std, 1.0e-6)

    f_npbc_bal_72 = free_energy_from_hist(p_npbc_bal_72, args.temp_k)
    f_npbc_bal_144 = free_energy_from_hist(p_npbc_bal_144, args.temp_k)

    balanced_metrics = {
        "balanced_jsd_72_mean": float(np.mean(boot_jsd72)),
        "balanced_jsd_72_std": float(np.std(boot_jsd72, ddof=1)),
        "balanced_jsd_72_ci": quantile_ci(boot_jsd72),
        "balanced_overlap_72_mean": float(np.mean(boot_ov72)),
        "balanced_overlap_72_std": float(np.std(boot_ov72, ddof=1)),
        "balanced_overlap_72_ci": quantile_ci(boot_ov72),
        "balanced_jsd_144_mean": float(np.mean(boot_jsd144)),
        "balanced_jsd_144_std": float(np.std(boot_jsd144, ddof=1)),
        "balanced_jsd_144_ci": quantile_ci(boot_jsd144),
        "balanced_overlap_144_mean": float(np.mean(boot_ov144)),
        "balanced_overlap_144_std": float(np.std(boot_ov144, ddof=1)),
        "balanced_overlap_144_ci": quantile_ci(boot_ov144),
    }

    # Basin and top-state tables
    basin_rows = []
    pbc_blocks = split_blocks(pbc_phi.size, args.blocks)
    npbc_blocks = split_blocks(npbc_phi.size, args.blocks)
    for (name, pdeg, sdeg) in BASIN_CENTERS:
        pbc_pop = basin_population(pbc_phi, pbc_psi, pdeg, sdeg, args.basin_window_deg)
        npbc_raw_pop = basin_population(npbc_phi, npbc_psi, pdeg, sdeg, args.basin_window_deg)
        pbc_sem = block_sem(
            np.array([basin_population(pbc_phi[idx], pbc_psi[idx], pdeg, sdeg, args.basin_window_deg) for idx in pbc_blocks]),
            len(pbc_blocks),
        )
        npbc_raw_sem = block_sem(
            np.array([basin_population(npbc_phi[idx], npbc_psi[idx], pdeg, sdeg, args.basin_window_deg) for idx in npbc_blocks]),
            len(npbc_blocks),
        )
        lo, hi = quantile_ci(basin_boot[name])
        basin_rows.append(
            {
                "basin": name,
                "phi_center_deg": pdeg,
                "psi_center_deg": sdeg,
                "pbc_pop": float(pbc_pop),
                "pbc_sem": float(pbc_sem),
                "npbc_raw_pop": float(npbc_raw_pop),
                "npbc_raw_sem": float(npbc_raw_sem),
                "npbc_balanced_mean": float(np.mean(basin_boot[name])),
                "npbc_balanced_std": float(np.std(basin_boot[name], ddof=1)),
                "npbc_balanced_ci_lo": lo,
                "npbc_balanced_ci_hi": hi,
                "delta_raw_minus_pbc": float(npbc_raw_pop - pbc_pop),
                "delta_balanced_minus_pbc": float(np.mean(basin_boot[name]) - pbc_pop),
            }
        )
    write_csv(
        outdir / "basin_populations.csv",
        [
            (
                r["basin"],
                r["phi_center_deg"],
                r["psi_center_deg"],
                r["pbc_pop"],
                r["pbc_sem"],
                r["npbc_raw_pop"],
                r["npbc_raw_sem"],
                r["npbc_balanced_mean"],
                r["npbc_balanced_std"],
                r["npbc_balanced_ci_lo"],
                r["npbc_balanced_ci_hi"],
                r["delta_raw_minus_pbc"],
                r["delta_balanced_minus_pbc"],
            )
            for r in basin_rows
        ],
        [
            "basin",
            "phi_center_deg",
            "psi_center_deg",
            "pbc_pop",
            "pbc_sem",
            "npbc_raw_pop",
            "npbc_raw_sem",
            "npbc_balanced_mean",
            "npbc_balanced_std",
            "npbc_balanced_ci_lo",
            "npbc_balanced_ci_hi",
            "delta_raw_minus_pbc",
            "delta_balanced_minus_pbc",
        ],
    )

    write_csv(
        outdir / "top_state_populations.csv",
        [
            (
                name,
                plo,
                phi,
                slo,
                shi,
                population_in_bin(pbc_phi, pbc_psi, plo, phi, slo, shi),
                population_in_bin(npbc_phi, npbc_psi, plo, phi, slo, shi),
                float(np.mean(state_boot[name])),
                float(np.std(state_boot[name], ddof=1)),
                quantile_ci(state_boot[name])[0],
                quantile_ci(state_boot[name])[1],
            )
            for (name, plo, phi, slo, shi) in top_states
        ],
        [
            "state",
            "phi_lo",
            "phi_hi",
            "psi_lo",
            "psi_hi",
            "pbc_pop",
            "npbc_raw_pop",
            "npbc_balanced_mean",
            "npbc_balanced_std",
            "npbc_balanced_ci_lo",
            "npbc_balanced_ci_hi",
        ],
    )

    # Time-block stability
    block_rows = []
    for lane, phi_lane, psi_lane, p_full in (
        ("NPBC_raw", npbc_phi, npbc_psi, p_npbc_72),
        ("PBC_raw", pbc_phi, pbc_psi, p_pbc_72),
    ):
        for block_name, idx in split_three_blocks(phi_lane.size):
            p_blk = hist2d_prob(phi_lane[idx], psi_lane[idx], bins72)
            block_rows.append(
                {
                    "lane": lane,
                    "block": block_name,
                    "nframes": int(idx.size),
                    "jsd_to_full": float(js_divergence_from_prob(p_blk, p_full)),
                    "overlap_to_full": float(overlap_coeff_from_prob(p_blk, p_full)),
                }
            )
    write_csv(
        outdir / "timeblock_stability.csv",
        [(r["lane"], r["block"], r["nframes"], r["jsd_to_full"], r["overlap_to_full"]) for r in block_rows],
        ["lane", "block", "nframes", "jsd_to_full", "overlap_to_full"],
    )

    # Kinetics-like basin diagnostics
    basin_state_order = [name for (name, _p, _s) in BASIN_CENTERS] + ["other"]
    states_npbc = assign_basin_series(npbc_phi, npbc_psi, args.basin_window_deg)
    states_pbc = assign_basin_series(pbc_phi, pbc_psi, args.basin_window_deg)
    trans_npbc_counts = transition_matrix(states_npbc, basin_state_order)
    trans_pbc_counts = transition_matrix(states_pbc, basin_state_order)
    trans_npbc_prob = row_normalize(trans_npbc_counts)
    trans_pbc_prob = row_normalize(trans_pbc_counts)
    trans_jsd = js_divergence_from_prob(trans_npbc_prob, trans_pbc_prob)
    trans_overlap = overlap_coeff_from_prob(trans_npbc_prob, trans_pbc_prob)

    # dump stride 500 with dt=1 fs => 0.5 ps per frame
    dt_ps = 0.5
    dwell_states_focus = ["alphaR", "alpha_prime", "C7eq", "PPII"]
    dwell_npbc = dwell_times(states_npbc, dwell_states_focus, dt_ps=dt_ps)
    dwell_pbc = dwell_times(states_pbc, dwell_states_focus, dt_ps=dt_ps)
    dwell_rows = []
    for st in dwell_states_focus:
        a = np.asarray(dwell_npbc.get(st, []), dtype=float)
        b = np.asarray(dwell_pbc.get(st, []), dtype=float)
        dwell_rows.append(
            (
                st,
                int(a.size),
                float(np.mean(a)) if a.size else float("nan"),
                float(np.median(a)) if a.size else float("nan"),
                int(b.size),
                float(np.mean(b)) if b.size else float("nan"),
                float(np.median(b)) if b.size else float("nan"),
            )
        )
    write_csv(
        outdir / "basin_dwell_summary.csv",
        dwell_rows,
        [
            "basin",
            "npbc_n_dwells",
            "npbc_mean_ps",
            "npbc_median_ps",
            "pbc_n_dwells",
            "pbc_mean_ps",
            "pbc_median_ps",
        ],
    )
    write_csv(
        outdir / "transition_matrix_npbc.csv",
        [(state_from, *trans_npbc_prob[i, :].tolist()) for i, state_from in enumerate(basin_state_order)],
        ["from_state", *[f"to_{s}" for s in basin_state_order]],
    )
    write_csv(
        outdir / "transition_matrix_pbc.csv",
        [(state_from, *trans_pbc_prob[i, :].tolist()) for i, state_from in enumerate(basin_state_order)],
        ["from_state", *[f"to_{s}" for s in basin_state_order]],
    )

    # PBC self-baseline (half-vs-half) for equivalence context
    half = pbc_phi.size // 2
    pbc_a_72 = hist2d_prob(pbc_phi[:half], pbc_psi[:half], bins72) if half > 10 else p_pbc_72
    pbc_b_72 = hist2d_prob(pbc_phi[half:], pbc_psi[half:], bins72) if pbc_phi.size - half > 10 else p_pbc_72
    pbc_self_jsd_72 = js_divergence_from_prob(pbc_a_72, pbc_b_72)
    pbc_self_overlap_72 = overlap_coeff_from_prob(pbc_a_72, pbc_b_72)

    # Structural metrics
    npbc_periodic_mask = np.array([is_periodic_flag(f) for f in d_npbc_eq.boundary_flags], dtype=bool)
    pbc_periodic_mask = np.array([is_periodic_flag(f) for f in d_pbc_eq.boundary_flags], dtype=bool)
    npbc_coords_used = (
        unwrap_solute_coords(npbc_coords, npbc_box, npbc_periodic_mask)
        if np.any(npbc_periodic_mask)
        else npbc_coords
    )
    pbc_coords_used = (
        unwrap_solute_coords(pbc_coords, pbc_box, pbc_periodic_mask)
        if np.any(pbc_periodic_mask)
        else pbc_coords
    )

    struct_npbc = compute_structural_metrics(
        npbc_coords_used, d_npbc_eq.solute_atom_ids, d_npbc_eq.solute_atom_types, args.rolling_window
    )
    struct_pbc = compute_structural_metrics(
        pbc_coords_used, d_pbc_eq.solute_atom_ids, d_pbc_eq.solute_atom_types, args.rolling_window
    )
    t_npbc = build_lane_time_ns(npbc_step, npbc_source)
    t_pbc = build_lane_time_ns(pbc_step, pbc_source)

    write_csv(
        outdir / "npbc_structural_timeseries.csv",
        [
            (
                int(i),
                str(npbc_source[i]),
                int(npbc_step[i]),
                float(t_npbc[i]),
                float(struct_npbc["rmsd_all"][i]),
                float(struct_npbc["rmsd_heavy"][i]),
                float(struct_npbc["rgyr_all"][i]),
                float(struct_npbc["rgyr_heavy"][i]),
                float(struct_npbc["dist_5_17"][i]),
            )
            for i in range(npbc_step.size)
        ],
        [
            "frame_index",
            "source",
            "step",
            "time_ns",
            "rmsd_all_A",
            "rmsd_heavy_A",
            "rgyr_all_A",
            "rgyr_heavy_A",
            "dist_5_17_A",
        ],
    )
    write_csv(
        outdir / "pbc_structural_timeseries.csv",
        [
            (
                int(i),
                str(pbc_source[i]),
                int(pbc_step[i]),
                float(t_pbc[i]),
                float(struct_pbc["rmsd_all"][i]),
                float(struct_pbc["rmsd_heavy"][i]),
                float(struct_pbc["rgyr_all"][i]),
                float(struct_pbc["rgyr_heavy"][i]),
                float(struct_pbc["dist_5_17"][i]),
            )
            for i in range(pbc_step.size)
        ],
        [
            "frame_index",
            "source",
            "step",
            "time_ns",
            "rmsd_all_A",
            "rmsd_heavy_A",
            "rgyr_all_A",
            "rgyr_heavy_A",
            "dist_5_17_A",
        ],
    )

    atom_ids = d_npbc_eq.solute_atom_ids
    atom_types = d_npbc_eq.solute_atom_types
    write_csv(
        outdir / "rmsf_per_atom.csv",
        [
            (
                int(atom_ids[i]),
                int(atom_types[i]),
                TYPE_TO_ELEMENT.get(int(atom_types[i]), str(int(atom_types[i]))),
                float(struct_npbc["rmsf"][i]),
                float(struct_pbc["rmsf"][i]),
                float(struct_npbc["rmsf"][i] - struct_pbc["rmsf"][i]),
            )
            for i in range(atom_ids.size)
        ],
        ["atom_id", "type", "element", "npbc_rmsf_A", "pbc_rmsf_A", "delta_npbc_minus_pbc_A"],
    )

    # Summary metrics
    rmsd_metrics = {
        "rmsd_all_npbc_mean": float(np.mean(struct_npbc["rmsd_all"])),
        "rmsd_all_pbc_mean": float(np.mean(struct_pbc["rmsd_all"])),
        "rmsd_heavy_npbc_mean": float(np.mean(struct_npbc["rmsd_heavy"])),
        "rmsd_heavy_pbc_mean": float(np.mean(struct_pbc["rmsd_heavy"])),
        "rmsd_all_ks_distance": ks_distance(struct_npbc["rmsd_all"], struct_pbc["rmsd_all"]),
        "rmsd_heavy_ks_distance": ks_distance(struct_npbc["rmsd_heavy"], struct_pbc["rmsd_heavy"]),
        "rgyr_all_npbc_mean": float(np.mean(struct_npbc["rgyr_all"])),
        "rgyr_all_pbc_mean": float(np.mean(struct_pbc["rgyr_all"])),
        "rgyr_heavy_npbc_mean": float(np.mean(struct_npbc["rgyr_heavy"])),
        "rgyr_heavy_pbc_mean": float(np.mean(struct_pbc["rgyr_heavy"])),
        "dist_5_17_npbc_mean": float(np.nanmean(struct_npbc["dist_5_17"])),
        "dist_5_17_pbc_mean": float(np.nanmean(struct_pbc["dist_5_17"])),
    }
    equivalence_metrics = {
        "pbc_self_jsd_72": float(pbc_self_jsd_72),
        "pbc_self_overlap_72": float(pbc_self_overlap_72),
        "npbc_pbc_jsd_72_raw": float(raw_metrics["raw_jsd_72"]),
        "npbc_pbc_overlap_72_raw": float(raw_metrics["raw_overlap_72"]),
        "jsd_ratio_npbc_pbc_over_pbc_self": float(raw_metrics["raw_jsd_72"] / max(pbc_self_jsd_72, 1.0e-12)),
        "overlap_gap_vs_pbc_self": float(pbc_self_overlap_72 - raw_metrics["raw_overlap_72"]),
        "transition_jsd": float(trans_jsd),
        "transition_overlap": float(trans_overlap),
    }

    jsd72_lo, jsd72_hi = balanced_metrics["balanced_jsd_72_ci"]
    ov72_lo, ov72_hi = balanced_metrics["balanced_overlap_72_ci"]
    robustness = {
        "bootstrap_jsd72_ci_width": float(jsd72_hi - jsd72_lo),
        "bootstrap_overlap72_ci_width": float(ov72_hi - ov72_lo),
        "bootstrap_ci_small_enough_flag": bool((jsd72_hi - jsd72_lo) <= 0.02 and (ov72_hi - ov72_lo) <= 0.05),
        "blockwise_drift_max_jsd_to_full": float(max(r["jsd_to_full"] for r in block_rows)),
        "blockwise_drift_min_overlap_to_full": float(min(r["overlap_to_full"] for r in block_rows)),
    }

    metrics_payload = {
        "raw_metrics": raw_metrics,
        "balanced_metrics": {k: (list(v) if isinstance(v, tuple) else v) for k, v in balanced_metrics.items()},
        "robustness": robustness,
        "structural_metrics": rmsd_metrics,
        "equivalence_metrics": equivalence_metrics,
    }
    (outdir / "metrics_summary.json").write_text(json.dumps(metrics_payload, indent=2) + "\n")

    integrity = {
        "required_columns": ["id", "mol", "type", "x", "y", "z"],
        "files": {
            "npbc_eq": {
                "path": str(npbc_eq_dump),
                "boundary_flags": list(d_npbc_eq.boundary_flags),
                "total_frames": d_npbc_eq.total_frames,
                "frames_with_required_atoms": d_npbc_eq.frames_with_required_atoms,
                "finite_dihedral_frames": d_npbc_eq.finite_frames,
                "columns_ok": d_npbc_eq.columns_ok,
            },
            "npbc_prod": {
                "path": str(npbc_prod_dump) if d_npbc_prod is not None else None,
                "boundary_flags": list(d_npbc_prod.boundary_flags) if d_npbc_prod is not None else None,
                "total_frames": d_npbc_prod.total_frames if d_npbc_prod is not None else 0,
                "frames_with_required_atoms": d_npbc_prod.frames_with_required_atoms if d_npbc_prod is not None else 0,
                "finite_dihedral_frames": d_npbc_prod.finite_frames if d_npbc_prod is not None else 0,
                "columns_ok": d_npbc_prod.columns_ok if d_npbc_prod is not None else None,
            },
            "pbc_eq": {
                "path": str(pbc_eq_dump),
                "boundary_flags": list(d_pbc_eq.boundary_flags),
                "total_frames": d_pbc_eq.total_frames,
                "frames_with_required_atoms": d_pbc_eq.frames_with_required_atoms,
                "finite_dihedral_frames": d_pbc_eq.finite_frames,
                "columns_ok": d_pbc_eq.columns_ok,
            },
            "pbc_prod": {
                "path": str(pbc_prod_dump) if d_pbc_prod is not None else None,
                "boundary_flags": list(d_pbc_prod.boundary_flags) if d_pbc_prod is not None else None,
                "total_frames": d_pbc_prod.total_frames if d_pbc_prod is not None else 0,
                "frames_with_required_atoms": d_pbc_prod.frames_with_required_atoms if d_pbc_prod is not None else 0,
                "finite_dihedral_frames": d_pbc_prod.finite_frames if d_pbc_prod is not None else 0,
                "columns_ok": d_pbc_prod.columns_ok if d_pbc_prod is not None else None,
            },
        },
        "combined_counts": {
            "npbc_total_finite_frames": int(npbc_phi.size),
            "pbc_total_finite_frames": int(pbc_phi.size),
        },
        "dihedral_sanity": {
            "npbc_finite_fraction": finite_frac_npbc,
            "pbc_finite_fraction": finite_frac_pbc,
            "npbc_phi_minmax": [float(np.min(npbc_phi)), float(np.max(npbc_phi))],
            "npbc_psi_minmax": [float(np.min(npbc_psi)), float(np.max(npbc_psi))],
            "pbc_phi_minmax": [float(np.min(pbc_phi)), float(np.max(pbc_phi))],
            "pbc_psi_minmax": [float(np.min(pbc_psi)), float(np.max(pbc_psi))],
            "finite_fraction_pass": bool(finite_frac_npbc > 0.99 and finite_frac_pbc > 0.99),
            "range_pass": True,
        },
    }
    (outdir / "input_integrity.json").write_text(json.dumps(integrity, indent=2) + "\n")

    # Plots
    maybe_plot_fes_triptych(outdir / "01_raw_fes_72.png", f_npbc_72, f_pbc_72, "Raw FES 72x72")
    maybe_plot_fes_triptych(outdir / "02_raw_fes_144.png", f_npbc_144, f_pbc_144, "Raw FES 144x144")
    maybe_plot_fes_triptych(outdir / "03_balanced_fes_72.png", f_npbc_bal_72, f_pbc_72, "Balanced FES 72x72")
    maybe_plot_fes_triptych(outdir / "04_balanced_fes_144.png", f_npbc_bal_144, f_pbc_144, "Balanced FES 144x144")
    maybe_plot_bootstrap_metrics(
        outdir / "05_bootstrap_metric_distributions.png",
        boot_jsd72,
        boot_ov72,
        raw_metrics["raw_jsd_72"],
        raw_metrics["raw_overlap_72"],
    )
    maybe_plot_basin_populations(outdir / "06_basin_populations_raw_vs_balanced.png", basin_rows)
    maybe_plot_timeblock_stability(outdir / "07_timeblock_stability.png", block_rows)
    maybe_plot_marginals(
        outdir / "08_phi_psi_marginals_with_balanced_ci.png",
        bins72,
        npbc_phi,
        pbc_phi,
        npbc_psi,
        pbc_psi,
        boot_phi_prob,
        boot_psi_prob,
    )
    maybe_plot_fe_profiles_1d(
        outdir / "18_free_energy_profiles_1d.png",
        outdir / "free_energy_profiles_1d.csv",
        bins72,
        args.temp_k,
        npbc_phi,
        pbc_phi,
        npbc_psi,
        pbc_psi,
        boot_phi_prob,
        boot_psi_prob,
    )
    maybe_plot_bootstrap_zscore_map(outdir / "09_bootstrap_zscore_map_72.png", zmap_72)
    maybe_plot_rmsd_timeseries(
        outdir / "10_rmsd_timeseries.png",
        t_npbc,
        t_pbc,
        struct_npbc["rmsd_all"],
        struct_npbc["rmsd_heavy"],
        struct_pbc["rmsd_all"],
        struct_pbc["rmsd_heavy"],
        struct_npbc["rmsd_all_roll"],
        struct_npbc["rmsd_heavy_roll"],
        struct_pbc["rmsd_all_roll"],
        struct_pbc["rmsd_heavy_roll"],
    )
    maybe_plot_rmsd_distribution(
        outdir / "11_rmsd_distributions.png",
        struct_npbc["rmsd_all"],
        struct_pbc["rmsd_all"],
        struct_npbc["rmsd_heavy"],
        struct_pbc["rmsd_heavy"],
    )
    maybe_plot_rmsf_per_atom(
        outdir / "12_rmsf_per_atom.png",
        atom_ids,
        atom_types,
        struct_npbc["rmsf"],
        struct_pbc["rmsf"],
    )
    maybe_plot_rgyr_timeseries(
        outdir / "13_rgyr_timeseries.png",
        t_npbc,
        t_pbc,
        struct_npbc["rgyr_all"],
        struct_npbc["rgyr_heavy"],
        struct_pbc["rgyr_all"],
        struct_pbc["rgyr_heavy"],
        struct_npbc["rgyr_all_roll"],
        struct_npbc["rgyr_heavy_roll"],
        struct_pbc["rgyr_all_roll"],
        struct_pbc["rgyr_heavy_roll"],
    )
    maybe_plot_distance_timeseries(
        outdir / "14_distance_5_17_timeseries.png",
        t_npbc,
        t_pbc,
        struct_npbc["dist_5_17"],
        struct_pbc["dist_5_17"],
    )
    maybe_plot_transition_matrix(
        outdir / "15_transition_matrix_npbc.png",
        trans_npbc_prob,
        basin_state_order,
        "NPBC basin transition matrix",
    )
    maybe_plot_transition_matrix(
        outdir / "16_transition_matrix_pbc.png",
        trans_pbc_prob,
        basin_state_order,
        "PBC basin transition matrix",
    )
    maybe_plot_dwell_times(
        outdir / "17_dwell_distributions_focus_basins.png",
        dwell_npbc,
        dwell_pbc,
        dwell_states_focus,
    )

    report = []
    report.append("# Stage13 Alanine Dihedral + Structural Analysis")
    report.append("")
    report.append(f"- Generated: {datetime.now().isoformat(timespec='seconds')}")
    report.append(f"- Stage13 root: `{stage13_root}`")
    report.append(f"- Output folder: `{outdir}`")
    report.append("")
    report.append("## Inputs")
    if d_npbc_prod is not None:
        report.append(f"- NPBC frames (eq+prod): `{npbc_phi.size}`")
    else:
        report.append(f"- NPBC frames (eq only): `{npbc_phi.size}`")
    if d_pbc_prod is not None:
        report.append(f"- PBC frames (eq+prod): `{pbc_phi.size}`")
    else:
        report.append(f"- PBC frames (eq only): `{pbc_phi.size}`")
    report.append(f"- Bootstrap reps: `{reps}`")
    report.append(
        f"- Structural metrics use unwrapped coordinates where periodic boundaries are present: NPBC periodic=`{bool(np.any(npbc_periodic_mask))}`, PBC periodic=`{bool(np.any(pbc_periodic_mask))}`"
    )
    report.append("")
    report.append("## Explicit caution")
    npbc_mode_name = "eq+prod" if d_npbc_prod is not None else "eq"
    pbc_mode_name = "eq+prod" if d_pbc_prod is not None else "eq"
    if npbc_phi.size != pbc_phi.size:
        report.append(f"- Raw NPBC({npbc_mode_name}) vs PBC({pbc_mode_name}) comparisons are length-biased.")
        report.append("- Balanced bootstrap (NPBC subsampled to PBC frame count) is the fair comparison.")
    else:
        report.append(f"- This is a matched comparison (NPBC {npbc_mode_name} vs PBC {pbc_mode_name}).")
        report.append("- Balanced bootstrap is still reported for consistency/uncertainty quantification.")
    report.append("")
    report.append("## Raw FES similarity")
    report.append(f"- JSD 72x72: `{raw_metrics['raw_jsd_72']:.5f}`")
    report.append(f"- Overlap 72x72: `{raw_metrics['raw_overlap_72']:.5f}`")
    report.append(f"- JSD 144x144: `{raw_metrics['raw_jsd_144']:.5f}`")
    report.append(f"- Overlap 144x144: `{raw_metrics['raw_overlap_144']:.5f}`")
    report.append("")
    report.append("## Balanced FES similarity")
    report.append(
        f"- JSD 72x72 mean±std: `{balanced_metrics['balanced_jsd_72_mean']:.5f} ± {balanced_metrics['balanced_jsd_72_std']:.5f}`"
    )
    report.append(
        f"- JSD 72x72 95% CI: `[{balanced_metrics['balanced_jsd_72_ci'][0]:.5f}, {balanced_metrics['balanced_jsd_72_ci'][1]:.5f}]`"
    )
    report.append(
        f"- Overlap 72x72 mean±std: `{balanced_metrics['balanced_overlap_72_mean']:.5f} ± {balanced_metrics['balanced_overlap_72_std']:.5f}`"
    )
    report.append(
        f"- Overlap 72x72 95% CI: `[{balanced_metrics['balanced_overlap_72_ci'][0]:.5f}, {balanced_metrics['balanced_overlap_72_ci'][1]:.5f}]`"
    )
    report.append("")
    report.append("## Structural differences")
    report.append(
        f"- RMSD(all) mean NPBC/PBC: `{rmsd_metrics['rmsd_all_npbc_mean']:.3f} / {rmsd_metrics['rmsd_all_pbc_mean']:.3f}` A"
    )
    report.append(
        f"- RMSD(heavy) mean NPBC/PBC: `{rmsd_metrics['rmsd_heavy_npbc_mean']:.3f} / {rmsd_metrics['rmsd_heavy_pbc_mean']:.3f}` A"
    )
    report.append(
        f"- KS distance RMSD(all/heavy): `{rmsd_metrics['rmsd_all_ks_distance']:.4f}` / `{rmsd_metrics['rmsd_heavy_ks_distance']:.4f}`"
    )
    report.append(
        f"- Rg(all) mean NPBC/PBC: `{rmsd_metrics['rgyr_all_npbc_mean']:.3f} / {rmsd_metrics['rgyr_all_pbc_mean']:.3f}` A"
    )
    report.append(
        f"- d(5,17) mean NPBC/PBC: `{rmsd_metrics['dist_5_17_npbc_mean']:.3f} / {rmsd_metrics['dist_5_17_pbc_mean']:.3f}` A"
    )
    report.append("")
    report.append("## Robustness")
    report.append(f"- Bootstrap CI width JSD72: `{robustness['bootstrap_jsd72_ci_width']:.5f}`")
    report.append(f"- Bootstrap CI width overlap72: `{robustness['bootstrap_overlap72_ci_width']:.5f}`")
    report.append(f"- Small-CI flag: `{robustness['bootstrap_ci_small_enough_flag']}`")
    report.append(
        f"- Max block JSD-to-full: `{robustness['blockwise_drift_max_jsd_to_full']:.5f}`"
    )
    report.append(
        f"- Min block overlap-to-full: `{robustness['blockwise_drift_min_overlap_to_full']:.5f}`"
    )
    report.append("")
    report.append("## Equivalence-context metrics")
    report.append(f"- PBC self JSD(72, half-vs-half): `{equivalence_metrics['pbc_self_jsd_72']:.5f}`")
    report.append(f"- NPBC-vs-PBC raw JSD(72): `{equivalence_metrics['npbc_pbc_jsd_72_raw']:.5f}`")
    report.append(
        f"- JSD ratio (NPBC-PBC / PBC-self): `{equivalence_metrics['jsd_ratio_npbc_pbc_over_pbc_self']:.3f}`"
    )
    report.append(f"- PBC self overlap(72): `{equivalence_metrics['pbc_self_overlap_72']:.5f}`")
    report.append(f"- NPBC-vs-PBC raw overlap(72): `{equivalence_metrics['npbc_pbc_overlap_72_raw']:.5f}`")
    report.append(f"- Transition-matrix JSD: `{equivalence_metrics['transition_jsd']:.5f}`")
    report.append(f"- Transition-matrix overlap: `{equivalence_metrics['transition_overlap']:.5f}`")
    report.append("")
    report.append("## Added outputs")
    report.append("- `transition_matrix_npbc.csv`, `transition_matrix_pbc.csv`")
    report.append("- `basin_dwell_summary.csv`")
    report.append("- `15_transition_matrix_npbc.png`")
    report.append("- `16_transition_matrix_pbc.png`")
    report.append("- `17_dwell_distributions_focus_basins.png`")
    report.append("- `18_free_energy_profiles_1d.png`, `free_energy_profiles_1d.csv`")
    (outdir / "report.md").write_text("\n".join(report) + "\n")

    print(f"Comprehensive stage13 analysis complete: {outdir}")


if __name__ == "__main__":
    main()
