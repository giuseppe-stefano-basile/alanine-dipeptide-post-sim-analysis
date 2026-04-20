#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
import importlib.util
from datetime import datetime
from pathlib import Path

import numpy as np


def load_analysis_module(module_path: Path):
    spec = importlib.util.spec_from_file_location("stage13_dih", module_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage13_dih"] = mod
    spec.loader.exec_module(mod)
    return mod


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Full NPBC vs PBC solute-dynamics comparison for stage13corr (handles unequal trajectory lengths)."
    )
    ap.add_argument(
        "--stage-root",
        default="Alanine_dipeptide/MACE-OFF_2023_small/stage13_alanine_prepared/stage13corr_s991000_20260319_194116",
    )
    ap.add_argument(
        "--analysis-module",
        default="Alanine_dipeptide/MACE-OFF_2023_small/analyze_stage13_dihedrals.py",
    )
    ap.add_argument(
        "--npbc-eq-dump",
        default="npbc/traj_alanine_nbpc_off23_stage13corr_s991000_20260319_194116_eq.dump",
    )
    ap.add_argument(
        "--npbc-prod-dump",
        default="npbc/traj_alanine_nbpc_off23_stage13corr_s991000_20260319_194116_prod.dump",
    )
    ap.add_argument(
        "--pbc-eq-dump",
        default="pbc/traj_alanine_pbc_off23_stage13corr_s991000_20260319_194116_eq.dump",
    )
    ap.add_argument(
        "--pbc-prod-dump",
        default="pbc/traj_alanine_pbc_off23_stage13corr_s991000_20260319_194116_prod.dump",
    )
    ap.add_argument("--bootstrap-reps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--basin-window-deg", type=float, default=20.0)
    ap.add_argument("--rolling-window", type=int, default=200)
    ap.add_argument(
        "--acf-max-lag",
        type=int,
        default=2000,
        help="Maximum lag (frames) for default torsion autocorrelation analysis.",
    )
    ap.add_argument("--temp-k", type=float, default=300.0)
    ap.add_argument(
        "--fes-bins-coarse",
        type=int,
        default=72,
        help="Number of bins per axis for coarse phi/psi FES grid.",
    )
    ap.add_argument(
        "--fes-bins-fine",
        type=int,
        default=144,
        help="Number of bins per axis for fine phi/psi FES grid.",
    )
    ap.add_argument(
        "--density-method",
        choices=("hist", "kde"),
        default="kde",
        help="Probability estimator for FES/marginals.",
    )
    ap.add_argument(
        "--kde-bandwidth-deg",
        type=float,
        default=18.0,
        help="Angular Gaussian bandwidth (deg) for periodic KDE.",
    )
    ap.add_argument(
        "--basin-definition",
        choices=("fixed", "fes"),
        default="fes",
        help="Basin definition: fixed center windows or FES-derived masks.",
    )
    ap.add_argument(
        "--fes-smooth-sigma-bins",
        type=float,
        default=1.0,
        help="Periodic Gaussian smoothing sigma (in grid bins) for FES-based basin masks.",
    )
    ap.add_argument(
        "--bootstrap-mode",
        choices=("simple", "block"),
        default="block",
        help="Bootstrap resampling mode. 'block' preserves short-time correlations.",
    )
    ap.add_argument(
        "--bootstrap-block-len",
        type=int,
        default=0,
        help="Block length (frames) for block bootstrap; 0 => auto from estimated IACT.",
    )
    ap.add_argument(
        "--crop-unsampled-fes",
        action="store_true",
        help="Mask unsampled bins in FES plots (set to NaN so they are visually cropped).",
    )
    ap.add_argument(
        "--crop-min-prob",
        type=float,
        default=0.0,
        help="Bins with probability <= this threshold are considered unsampled for FES cropping.",
    )
    ap.add_argument(
        "--crop-policy",
        choices=("either", "both"),
        default="either",
        help="Cropping rule across lanes: either masks if unsampled in NPBC or PBC, both only if unsampled in both.",
    )
    ap.add_argument("--outdir", default=None)
    return ap.parse_args()


def _must_exist(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")


def circular_delta_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return ((a - b + 180.0) % 360.0) - 180.0


def periodic_kde_1d(samples_deg: np.ndarray, bins: np.ndarray, bandwidth_deg: float) -> np.ndarray:
    x = np.asarray(samples_deg, dtype=float).ravel()
    mids = 0.5 * (bins[:-1] + bins[1:])
    if x.size == 0:
        return np.zeros(mids.size, dtype=float)
    bw = max(float(bandwidth_deg), 1.0e-6)
    d = circular_delta_deg(x[:, None], mids[None, :])
    k = np.exp(-0.5 * (d / bw) ** 2)
    dens = np.sum(k, axis=0) / float(x.size * bw * np.sqrt(2.0 * np.pi))
    s = float(np.sum(dens))
    if s > 0.0:
        dens /= s
    return dens


def periodic_kde_2d(phi_deg: np.ndarray, psi_deg: np.ndarray, bins: np.ndarray, bandwidth_deg: float) -> np.ndarray:
    phi = np.asarray(phi_deg, dtype=float).ravel()
    psi = np.asarray(psi_deg, dtype=float).ravel()
    mids = 0.5 * (bins[:-1] + bins[1:])
    nb = mids.size
    if phi.size == 0 or psi.size == 0:
        return np.zeros((nb, nb), dtype=float)
    bw = max(float(bandwidth_deg), 1.0e-6)
    dphi = circular_delta_deg(phi[:, None], mids[None, :])
    dpsi = circular_delta_deg(psi[:, None], mids[None, :])
    kphi = np.exp(-0.5 * (dphi / bw) ** 2)
    kpsi = np.exp(-0.5 * (dpsi / bw) ** 2)
    dens = (kphi.T @ kpsi) / float(phi.size * (2.0 * np.pi * bw * bw))
    s = float(np.sum(dens))
    if s > 0.0:
        dens /= s
    return dens


def density2d_prob(phi: np.ndarray, psi: np.ndarray, bins: np.ndarray, method: str, kde_bw_deg: float, mod) -> np.ndarray:
    if method == "kde":
        return periodic_kde_2d(phi, psi, bins, kde_bw_deg)
    return mod.hist2d_prob(phi, psi, bins)


def density1d_prob(x: np.ndarray, bins: np.ndarray, method: str, kde_bw_deg: float, mod) -> np.ndarray:
    if method == "kde":
        return periodic_kde_1d(x, bins, kde_bw_deg)
    return mod.hist1d_prob(x, bins)


def _autocorr_1d(y: np.ndarray, max_lag: int) -> np.ndarray:
    x = np.asarray(y, dtype=float).ravel()
    n = x.size
    if n < 2:
        return np.array([1.0], dtype=float)
    x = x - float(np.mean(x))
    c0 = float(np.dot(x, x))
    if c0 <= 0.0:
        return np.array([1.0], dtype=float)
    max_lag = max(1, min(int(max_lag), n - 1))
    ac = np.empty(max_lag + 1, dtype=float)
    ac[0] = 1.0
    for lag in range(1, max_lag + 1):
        ac[lag] = float(np.dot(x[:-lag], x[lag:]) / c0)
    return ac


def estimate_iact_from_series(y: np.ndarray, max_lag: int = 2000) -> float:
    ac = _autocorr_1d(y, max_lag=max_lag)
    if ac.size <= 1:
        return 1.0
    pos = np.where(ac[1:] <= 0.0)[0]
    stop = int(pos[0] + 1) if pos.size > 0 else int(ac.size)
    iact = 1.0 + 2.0 * float(np.sum(ac[1:stop]))
    if not np.isfinite(iact) or iact < 1.0:
        iact = 1.0
    return iact


def estimate_angle_iact_frames(angle_deg: np.ndarray, max_lag: int = 2000) -> float:
    rad = np.deg2rad(np.asarray(angle_deg, dtype=float).ravel())
    if rad.size < 4:
        return 1.0
    iact_cos = estimate_iact_from_series(np.cos(rad), max_lag=max_lag)
    iact_sin = estimate_iact_from_series(np.sin(rad), max_lag=max_lag)
    return max(1.0, iact_cos, iact_sin)


def torsion_autocorr_cosdelta(angle_deg: np.ndarray, max_lag: int) -> np.ndarray:
    x = np.deg2rad(np.asarray(angle_deg, dtype=float).ravel())
    n = int(x.size)
    if n < 2:
        return np.array([1.0], dtype=float)
    mlag = max(1, min(int(max_lag), n - 1))
    ac = np.empty(mlag + 1, dtype=float)
    ac[0] = 1.0
    for lag in range(1, mlag + 1):
        ac[lag] = float(np.mean(np.cos(x[:-lag] - x[lag:])))
    return ac


def torsion_autocorr_centered_components(angle_deg: np.ndarray, max_lag: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rad = np.deg2rad(np.asarray(angle_deg, dtype=float).ravel())
    n = int(rad.size)
    if n < 2:
        z = np.array([1.0], dtype=float)
        return z, z, z
    mlag = max(1, min(int(max_lag), n - 1))
    c = np.cos(rad)
    s = np.sin(rad)
    ac_c = _autocorr_1d(c, max_lag=mlag)
    ac_s = _autocorr_1d(s, max_lag=mlag)
    var_c = float(np.var(c))
    var_s = float(np.var(s))
    w = var_c + var_s
    if w <= 0.0:
        ac_mix = 0.5 * (ac_c + ac_s)
    else:
        ac_mix = (var_c * ac_c + var_s * ac_s) / w
    return ac_mix, ac_c, ac_s


def maybe_plot_torsion_autocorr(
    out_png: Path,
    lag_npbc_ps: np.ndarray,
    lag_pbc_ps: np.ndarray,
    acf_npbc_phi: np.ndarray,
    acf_pbc_phi: np.ndarray,
    acf_npbc_psi: np.ndarray,
    acf_pbc_psi: np.ndarray,
    have_plot: bool,
) -> None:
    if not have_plot:
        return
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    axs[0].plot(lag_npbc_ps, acf_npbc_phi, label="NPBC phi", lw=1.3)
    axs[0].plot(lag_pbc_ps, acf_pbc_phi, label="PBC phi", lw=1.3)
    axs[0].axhline(0.0, color="k", lw=0.8, alpha=0.4)
    axs[0].set_title("phi autocorrelation")
    axs[0].set_xlabel("lag (ps)")
    axs[0].set_ylabel("<cos(delta phi)>")
    axs[0].grid(alpha=0.25)
    axs[0].legend()

    axs[1].plot(lag_npbc_ps, acf_npbc_psi, label="NPBC psi", lw=1.3)
    axs[1].plot(lag_pbc_ps, acf_pbc_psi, label="PBC psi", lw=1.3)
    axs[1].axhline(0.0, color="k", lw=0.8, alpha=0.4)
    axs[1].set_title("psi autocorrelation")
    axs[1].set_xlabel("lag (ps)")
    axs[1].set_ylabel("<cos(delta psi)>")
    axs[1].grid(alpha=0.25)
    axs[1].legend()

    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_torsion_autocorr_centered(
    out_png: Path,
    lag_npbc_ps: np.ndarray,
    lag_pbc_ps: np.ndarray,
    acf_npbc_phi_centered: np.ndarray,
    acf_pbc_phi_centered: np.ndarray,
    acf_npbc_psi_centered: np.ndarray,
    acf_pbc_psi_centered: np.ndarray,
    have_plot: bool,
) -> None:
    if not have_plot:
        return
    import matplotlib.pyplot as plt

    fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
    axs[0].plot(lag_npbc_ps, acf_npbc_phi_centered, label="NPBC phi", lw=1.3)
    axs[0].plot(lag_pbc_ps, acf_pbc_phi_centered, label="PBC phi", lw=1.3)
    axs[0].axhline(0.0, color="k", lw=0.8, alpha=0.4)
    axs[0].set_title("phi autocorrelation (centered sin/cos)")
    axs[0].set_xlabel("lag (ps)")
    axs[0].set_ylabel("ACF")
    axs[0].grid(alpha=0.25)
    axs[0].legend()

    axs[1].plot(lag_npbc_ps, acf_npbc_psi_centered, label="NPBC psi", lw=1.3)
    axs[1].plot(lag_pbc_ps, acf_pbc_psi_centered, label="PBC psi", lw=1.3)
    axs[1].axhline(0.0, color="k", lw=0.8, alpha=0.4)
    axs[1].set_title("psi autocorrelation (centered sin/cos)")
    axs[1].set_xlabel("lag (ps)")
    axs[1].set_ylabel("ACF")
    axs[1].grid(alpha=0.25)
    axs[1].legend()

    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def resolve_block_len_frames(
    user_block_len: int,
    n_target: int,
    npbc_phi: np.ndarray,
    npbc_psi: np.ndarray,
    pbc_phi: np.ndarray,
    pbc_psi: np.ndarray,
) -> tuple[int, float]:
    if user_block_len and user_block_len > 0:
        block = max(2, min(int(user_block_len), max(2, n_target // 2)))
        return block, float(block)
    iacts = [
        estimate_angle_iact_frames(npbc_phi),
        estimate_angle_iact_frames(npbc_psi),
        estimate_angle_iact_frames(pbc_phi),
        estimate_angle_iact_frames(pbc_psi),
    ]
    iact_est = float(np.max(np.asarray(iacts, dtype=float)))
    block = int(max(5, round(iact_est)))
    block = max(2, min(block, max(2, n_target // 8)))
    return block, iact_est


def bootstrap_indices(
    n: int, n_target: int, rng: np.random.Generator, mode: str, block_len: int
) -> np.ndarray:
    if mode == "simple":
        return rng.integers(0, n, size=n_target, dtype=int)
    b = max(2, min(int(block_len), max(2, n_target)))
    out = np.empty(n_target, dtype=int)
    k = 0
    while k < n_target:
        start = int(rng.integers(0, n))
        run = min(b, n_target - k)
        out[k : k + run] = (start + np.arange(run, dtype=int)) % n
        k += run
    return out


def gaussian_kernel_1d(sigma_bins: float) -> np.ndarray:
    s = float(sigma_bins)
    if s <= 0.0:
        return np.array([1.0], dtype=float)
    radius = max(1, int(np.ceil(3.0 * s)))
    x = np.arange(-radius, radius + 1, dtype=float)
    k = np.exp(-0.5 * (x / s) ** 2)
    k /= float(np.sum(k))
    return k


def smooth_periodic_2d(prob: np.ndarray, sigma_bins: float) -> np.ndarray:
    p = np.asarray(prob, dtype=float)
    k = gaussian_kernel_1d(sigma_bins)
    radius = k.size // 2
    if radius <= 0:
        out = p.copy()
    else:
        tmp = np.zeros_like(p)
        for i, w in enumerate(k):
            shift = i - radius
            tmp += w * np.roll(p, shift, axis=0)
        out = np.zeros_like(tmp)
        for i, w in enumerate(k):
            shift = i - radius
            out += w * np.roll(tmp, shift, axis=1)
    s = float(np.sum(out))
    if s > 0.0:
        out /= s
    return out


def _grid_centers_from_bins(bins: np.ndarray) -> np.ndarray:
    return 0.5 * (bins[:-1] + bins[1:])


def _nearest_basin_idx(phi_c: float, psi_c: float, basin_centers: list[tuple[str, float, float]]) -> int:
    best = 0
    best_d2 = float("inf")
    for i, (_name, pc, sc) in enumerate(basin_centers):
        dphi = float(abs(((phi_c - pc + 180.0) % 360.0) - 180.0))
        dpsi = float(abs(((psi_c - sc + 180.0) % 360.0) - 180.0))
        d2 = dphi * dphi + dpsi * dpsi
        if d2 < best_d2:
            best_d2 = d2
            best = i
    return best


def build_fes_basin_masks(
    p_ref: np.ndarray, bins: np.ndarray, basin_centers: list[tuple[str, float, float]], sigma_bins: float
) -> dict[str, np.ndarray]:
    p_s = smooth_periodic_2d(p_ref, sigma_bins=sigma_bins)
    nphi, npsi = p_s.shape
    peak_label = -np.ones((nphi, npsi), dtype=int)
    peak_to_id: dict[tuple[int, int], int] = {}

    def climb(i0: int, j0: int) -> tuple[int, int]:
        i, j = i0, j0
        for _ in range(nphi * npsi):
            cur = p_s[i, j]
            bi, bj = i, j
            bval = cur
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    ni = (i + di) % nphi
                    nj = (j + dj) % npsi
                    v = p_s[ni, nj]
                    if v > bval + 1.0e-15:
                        bval = v
                        bi, bj = ni, nj
            if bi == i and bj == j:
                return i, j
            i, j = bi, bj
        return i, j

    for i in range(nphi):
        for j in range(npsi):
            pi, pj = climb(i, j)
            key = (pi, pj)
            if key not in peak_to_id:
                peak_to_id[key] = len(peak_to_id)
            peak_label[i, j] = peak_to_id[key]

    mids = _grid_centers_from_bins(bins)
    peak_to_basin: dict[int, int] = {}
    for (pi, pj), pid in peak_to_id.items():
        peak_to_basin[pid] = _nearest_basin_idx(float(mids[pi]), float(mids[pj]), basin_centers)

    assign = np.zeros_like(peak_label, dtype=int)
    for pid, bidx in peak_to_basin.items():
        assign[peak_label == pid] = bidx

    masks = {}
    for bidx, (name, _pc, _sc) in enumerate(basin_centers):
        masks[name] = assign == bidx
    return masks


def build_fixed_basin_masks(
    bins: np.ndarray, basin_centers: list[tuple[str, float, float]], halfw_deg: float
) -> dict[str, np.ndarray]:
    mids = _grid_centers_from_bins(bins)
    xx, yy = np.meshgrid(mids, mids, indexing="ij")
    best_idx = -np.ones(xx.shape, dtype=int)
    best_d2 = np.full(xx.shape, np.inf, dtype=float)
    for bidx, (_name, pc, sc) in enumerate(basin_centers):
        dphi = np.abs(circular_delta_deg(xx, pc))
        dpsi = np.abs(circular_delta_deg(yy, sc))
        inside = (dphi <= halfw_deg) & (dpsi <= halfw_deg)
        d2 = dphi * dphi + dpsi * dpsi
        take = inside & (d2 < best_d2)
        best_d2[take] = d2[take]
        best_idx[take] = bidx
    masks = {}
    for bidx, (name, _pc, _sc) in enumerate(basin_centers):
        masks[name] = best_idx == bidx
    return masks


def basin_populations_from_prob(prob: np.ndarray, basin_masks: dict[str, np.ndarray]) -> dict[str, float]:
    p = np.asarray(prob, dtype=float)
    return {name: float(np.sum(p[mask])) for name, mask in basin_masks.items()}


def assign_states_from_basin_masks(
    phi: np.ndarray, psi: np.ndarray, bins: np.ndarray, basin_masks: dict[str, np.ndarray], basin_order: list[str]
) -> np.ndarray:
    bw = float(bins[1] - bins[0])
    nbin = int(bins.size - 1)
    ii = np.floor((phi + 180.0) / bw).astype(int)
    jj = np.floor((psi + 180.0) / bw).astype(int)
    ii = np.clip(ii, 0, nbin - 1)
    jj = np.clip(jj, 0, nbin - 1)
    states = np.full(phi.size, "other", dtype=object)
    for name in basin_order:
        mask = basin_masks[name][ii, jj]
        states[mask] = name
    return states


def maybe_plot_fes_triptych_cropped(
    mod,
    out_png: Path,
    f_npbc: np.ndarray,
    f_pbc: np.ndarray,
    p_npbc: np.ndarray,
    p_pbc: np.ndarray,
    title_prefix: str,
    do_crop: bool,
    min_prob: float,
    policy: str,
) -> dict:
    if not do_crop:
        mod.maybe_plot_fes_triptych(out_png, f_npbc, f_pbc, title_prefix)
        return {"enabled": False, "masked_fraction": 0.0}

    mask = build_crop_mask(p_npbc, p_pbc, True, min_prob, policy)

    if not mod.HAVE_PLOT:
        return {"enabled": True, "masked_fraction": float(np.mean(mask))}

    import matplotlib.pyplot as plt

    f_n = np.array(f_npbc, copy=True)
    f_p = np.array(f_pbc, copy=True)
    f_n[mask] = np.nan
    f_p[mask] = np.nan
    delta = f_n - f_p

    fig, axs = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    ext = [-180, 180, -180, 180]
    cmap_main = plt.get_cmap("viridis").copy()
    cmap_main.set_bad(color="white")
    cmap_delta = plt.get_cmap("coolwarm").copy()
    cmap_delta.set_bad(color="white")

    im0 = axs[0].imshow(f_n.T, origin="lower", extent=ext, aspect="auto", cmap=cmap_main)
    axs[0].set_title(f"{title_prefix}: NPBC")
    axs[0].set_xlabel("phi (deg)")
    axs[0].set_ylabel("psi (deg)")
    fig.colorbar(im0, ax=axs[0], shrink=0.9)

    im1 = axs[1].imshow(f_p.T, origin="lower", extent=ext, aspect="auto", cmap=cmap_main)
    axs[1].set_title(f"{title_prefix}: PBC")
    axs[1].set_xlabel("phi (deg)")
    axs[1].set_ylabel("psi (deg)")
    fig.colorbar(im1, ax=axs[1], shrink=0.9)

    finite = np.isfinite(delta)
    vmax = float(np.nanmax(np.abs(delta[finite]))) if np.any(finite) else 1.0
    vmax = 1.0 if (not np.isfinite(vmax) or vmax <= 0.0) else vmax
    im2 = axs[2].imshow(
        delta.T, origin="lower", extent=ext, aspect="auto", cmap=cmap_delta, vmin=-vmax, vmax=vmax
    )
    axs[2].set_title(f"{title_prefix}: DeltaF NPBC-PBC")
    axs[2].set_xlabel("phi (deg)")
    axs[2].set_ylabel("psi (deg)")
    fig.colorbar(im2, ax=axs[2], shrink=0.9)
    fig.suptitle(f"cropped unsampled bins: policy={policy}, min_prob={min_prob:g}", fontsize=10)
    fig.savefig(out_png, dpi=170)
    plt.close(fig)

    return {
        "enabled": True,
        "policy": policy,
        "min_prob": float(min_prob),
        "masked_fraction": float(np.mean(mask)),
        "masked_bins": int(np.sum(mask)),
        "total_bins": int(mask.size),
    }


def build_crop_mask(
    p_npbc: np.ndarray, p_pbc: np.ndarray, do_crop: bool, min_prob: float, policy: str
) -> np.ndarray:
    if not do_crop:
        return np.zeros_like(p_npbc, dtype=bool)
    mask_npbc = p_npbc <= float(min_prob)
    mask_pbc = p_pbc <= float(min_prob)
    if policy == "both":
        return mask_npbc & mask_pbc
    return mask_npbc | mask_pbc


def masked_renormalized_prob(prob: np.ndarray, mask: np.ndarray) -> np.ndarray:
    p = np.array(prob, copy=True, dtype=float)
    p[mask] = 0.0
    s = float(np.sum(p))
    if s > 0.0:
        p /= s
    return p


def transition_event_rows(
    lane: str,
    states: np.ndarray,
    steps: np.ndarray,
    times_ns: np.ndarray,
    sources: np.ndarray,
    dt_ps: float = 0.5,
) -> list[tuple]:
    out = []
    n = int(states.size)
    if n < 2:
        return out
    run_start = 0
    cur = str(states[0])
    for i in range(1, n):
        nxt = str(states[i])
        if nxt == cur:
            continue
        run_len = i - run_start
        out.append(
            (
                lane,
                int(run_start),
                int(i - 1),
                str(sources[run_start]),
                int(steps[run_start]),
                float(times_ns[run_start]),
                str(sources[i]),
                int(steps[i]),
                float(times_ns[i]),
                cur,
                nxt,
                int(run_len),
                float(run_len * dt_ps),
            )
        )
        run_start = i
        cur = nxt
    return out


def first_entry_rows(
    lane: str,
    states: np.ndarray,
    steps: np.ndarray,
    times_ns: np.ndarray,
    sources: np.ndarray,
    state_order: list[str],
) -> list[tuple]:
    out = []
    seen = set()
    for i, s in enumerate(states):
        st = str(s)
        if st == "other" or st not in state_order or st in seen:
            continue
        seen.add(st)
        out.append((lane, st, int(i), str(sources[i]), int(steps[i]), float(times_ns[i])))
    return out


def maybe_plot_transition_events_timeline(
    out_png: Path,
    transition_rows: list[tuple],
    state_order: list[str],
    first_entry_rows_data: list[tuple],
    have_plot: bool,
) -> None:
    if not have_plot:
        return
    import matplotlib.pyplot as plt

    state_to_i = {s: i for i, s in enumerate(state_order)}
    lanes = ("npbc", "pbc")
    fig, axs = plt.subplots(2, 1, figsize=(12, 7), sharex=True, constrained_layout=True)
    for ax, lane in zip(axs, lanes):
        rows_lane = [r for r in transition_rows if r[0] == lane and str(r[10]) != "other"]
        if len(rows_lane) == 0:
            ax.text(0.5, 0.5, f"{lane.upper()}: no transitions", transform=ax.transAxes, ha="center", va="center")
            continue
        x = np.asarray([float(r[8]) for r in rows_lane], dtype=float)
        y = np.asarray([state_to_i.get(str(r[10]), state_to_i.get("other", 0)) for r in rows_lane], dtype=float)
        c = np.asarray([state_to_i.get(str(r[9]), state_to_i.get("other", 0)) for r in rows_lane], dtype=float)
        ax.scatter(x, y, c=c, cmap="tab20", s=14, alpha=0.85, linewidths=0.0)

        first_lane = [r for r in first_entry_rows_data if r[0] == lane]
        for _ln, st, _fi, _src, _step, tns in first_lane:
            yi = state_to_i.get(str(st), None)
            if yi is None:
                continue
            ax.axvline(float(tns), color="k", ls=":", lw=0.7, alpha=0.25)
            ax.text(float(tns), yi + 0.15, str(st), fontsize=7, rotation=90, alpha=0.7, va="bottom", ha="center")

        ax.set_yticks(np.arange(len(state_order)))
        ax.set_yticklabels(state_order, fontsize=8)
        ax.set_ylabel("to_state")
        ax.set_title(f"{lane.upper()} transition events (dot color = from_state index)")
        ax.grid(alpha=0.25)
    axs[-1].set_xlabel("time (ns)")
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def maybe_plot_transition_event_cumulative(
    out_png: Path, transition_rows: list[tuple], have_plot: bool
) -> None:
    if not have_plot:
        return
    import matplotlib.pyplot as plt

    lanes = ("npbc", "pbc")
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.2), constrained_layout=True)
    for lane in lanes:
        t = np.asarray([float(r[8]) for r in transition_rows if r[0] == lane and str(r[10]) != "other"], dtype=float)
        if t.size == 0:
            continue
        t = np.sort(t)
        c = np.arange(1, t.size + 1, dtype=float)
        ax.step(t, c, where="post", label=f"{lane.upper()} (n={t.size})", lw=1.6)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("cumulative transition events")
    ax.set_title("Cumulative basin-transition events (excluding to_state='other')")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(out_png, dpi=170)
    plt.close(fig)


def main() -> None:
    args = parse_args()

    stage_root = Path(args.stage_root).resolve()
    module_path = Path(args.analysis_module).resolve()
    mod = load_analysis_module(module_path)

    npbc_eq = (stage_root / args.npbc_eq_dump).resolve()
    npbc_prod = (stage_root / args.npbc_prod_dump).resolve()
    pbc_eq = (stage_root / args.pbc_eq_dump).resolve()
    pbc_prod = (stage_root / args.pbc_prod_dump).resolve()
    for p in (npbc_eq, npbc_prod, pbc_eq, pbc_prod):
        _must_exist(p)

    req = [5, 7, 9, 15, 17]
    d_npbc_eq = mod.parse_dump_lane(npbc_eq, req, "npbc_eq")
    d_npbc_prod = mod.parse_dump_lane(npbc_prod, req, "npbc_prod")
    d_pbc_eq = mod.parse_dump_lane(pbc_eq, req, "pbc_eq")
    d_pbc_prod = mod.parse_dump_lane(pbc_prod, req, "pbc_prod")

    if np.any(d_npbc_eq.solute_atom_ids != d_pbc_eq.solute_atom_ids):
        raise RuntimeError("Solute atom-ID ordering mismatch NPBC vs PBC")
    if np.any(d_npbc_eq.solute_atom_types != d_pbc_eq.solute_atom_types):
        raise RuntimeError("Solute atom-type mismatch NPBC vs PBC")

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

    n_npbc = npbc_phi.size
    n_pbc = pbc_phi.size
    n_bal = min(n_npbc, n_pbc)

    if n_bal < 100:
        raise RuntimeError(f"Too few frames for robust comparison (min lane={n_bal})")

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        outdir = (stage_root / "solute_dynamics_comparison" / f"{stamp}_npbc_eqprod_vs_pbc_eqprod").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Structural metrics with periodic-aware unwrapping.
    npbc_periodic_mask = np.array([mod.is_periodic_flag(f) for f in d_npbc_eq.boundary_flags], dtype=bool)
    pbc_periodic_mask = np.array([mod.is_periodic_flag(f) for f in d_pbc_eq.boundary_flags], dtype=bool)
    npbc_coords_used = (
        mod.unwrap_solute_coords(npbc_coords, npbc_box, npbc_periodic_mask) if np.any(npbc_periodic_mask) else npbc_coords
    )
    pbc_coords_used = (
        mod.unwrap_solute_coords(pbc_coords, pbc_box, pbc_periodic_mask) if np.any(pbc_periodic_mask) else pbc_coords
    )

    struct_npbc = mod.compute_structural_metrics(
        npbc_coords_used, d_npbc_eq.solute_atom_ids, d_npbc_eq.solute_atom_types, args.rolling_window
    )
    struct_pbc = mod.compute_structural_metrics(
        pbc_coords_used, d_pbc_eq.solute_atom_ids, d_pbc_eq.solute_atom_types, args.rolling_window
    )

    t_npbc = mod.build_lane_time_ns(npbc_step, npbc_source)
    t_pbc = mod.build_lane_time_ns(pbc_step, pbc_source)

    mod.write_csv(
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
    mod.write_csv(
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
    mod.write_csv(
        outdir / "rmsf_per_atom.csv",
        [
            (
                int(atom_ids[i]),
                int(atom_types[i]),
                mod.TYPE_TO_ELEMENT.get(int(atom_types[i]), str(int(atom_types[i]))),
                float(struct_npbc["rmsf"][i]),
                float(struct_pbc["rmsf"][i]),
                float(struct_npbc["rmsf"][i] - struct_pbc["rmsf"][i]),
            )
            for i in range(atom_ids.size)
        ],
        ["atom_id", "type", "element", "npbc_rmsf_A", "pbc_rmsf_A", "delta_npbc_minus_pbc_A"],
    )

    # Dihedral/FES metrics.
    nbin_coarse = max(8, int(args.fes_bins_coarse))
    nbin_fine = max(8, int(args.fes_bins_fine))
    bins72 = np.linspace(-180.0, 180.0, nbin_coarse + 1)
    bins144 = np.linspace(-180.0, 180.0, nbin_fine + 1)
    basin_order = [name for (name, _pc, _sc) in mod.BASIN_CENTERS]

    p_npbc_72 = density2d_prob(npbc_phi, npbc_psi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
    p_pbc_72 = density2d_prob(pbc_phi, pbc_psi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
    p_npbc_144 = density2d_prob(npbc_phi, npbc_psi, bins144, args.density_method, args.kde_bandwidth_deg, mod)
    p_pbc_144 = density2d_prob(pbc_phi, pbc_psi, bins144, args.density_method, args.kde_bandwidth_deg, mod)

    raw_jsd_72 = mod.js_divergence_from_prob(p_npbc_72, p_pbc_72)
    raw_ov_72 = mod.overlap_coeff_from_prob(p_npbc_72, p_pbc_72)
    raw_jsd_144 = mod.js_divergence_from_prob(p_npbc_144, p_pbc_144)
    raw_ov_144 = mod.overlap_coeff_from_prob(p_npbc_144, p_pbc_144)

    f_npbc_72 = mod.free_energy_from_hist(p_npbc_72, args.temp_k)
    f_pbc_72 = mod.free_energy_from_hist(p_pbc_72, args.temp_k)
    f_npbc_144 = mod.free_energy_from_hist(p_npbc_144, args.temp_k)
    f_pbc_144 = mod.free_energy_from_hist(p_pbc_144, args.temp_k)

    # Shared crop mask on coarse grid for crop-aware basin populations/marginals.
    crop_mask_72 = build_crop_mask(
        p_npbc_72, p_pbc_72, args.crop_unsampled_fes, args.crop_min_prob, args.crop_policy
    )
    if args.crop_unsampled_fes:
        p_npbc_72_basin = masked_renormalized_prob(p_npbc_72, crop_mask_72)
        p_pbc_72_basin = masked_renormalized_prob(p_pbc_72, crop_mask_72)
    else:
        p_npbc_72_basin = p_npbc_72
        p_pbc_72_basin = p_pbc_72

    # Basin definitions from fixed windows or smoothed FES regions.
    if args.basin_definition == "fes":
        pref = 0.5 * (p_npbc_72 + p_pbc_72)
        ps = float(np.sum(pref))
        if ps > 0.0:
            pref /= ps
        basin_masks = build_fes_basin_masks(pref, bins72, mod.BASIN_CENTERS, args.fes_smooth_sigma_bins)
    else:
        basin_masks = build_fixed_basin_masks(bins72, mod.BASIN_CENTERS, args.basin_window_deg)

    raw_basin_npbc = basin_populations_from_prob(p_npbc_72_basin, basin_masks)
    raw_basin_pbc = basin_populations_from_prob(p_pbc_72_basin, basin_masks)

    # Balanced bootstrap: sample both lanes to min(N_npbc, N_pbc).
    rng = np.random.default_rng(args.seed)
    reps = int(args.bootstrap_reps)
    boot_jsd72 = np.zeros(reps, dtype=float)
    boot_ov72 = np.zeros(reps, dtype=float)
    boot_jsd144 = np.zeros(reps, dtype=float)
    boot_ov144 = np.zeros(reps, dtype=float)

    boot_phi_prob_npbc = np.zeros((reps, nbin_coarse), dtype=float)
    boot_psi_prob_npbc = np.zeros((reps, nbin_coarse), dtype=float)
    boot_phi_prob_pbc = np.zeros((reps, nbin_coarse), dtype=float)
    boot_psi_prob_pbc = np.zeros((reps, nbin_coarse), dtype=float)
    boot_p_npbc72 = np.zeros((reps, nbin_coarse, nbin_coarse), dtype=float) if args.crop_unsampled_fes else None
    boot_p_pbc72 = np.zeros((reps, nbin_coarse, nbin_coarse), dtype=float) if args.crop_unsampled_fes else None

    acc_npbc72 = np.zeros_like(p_npbc_72)
    acc_pbc72 = np.zeros_like(p_pbc_72)
    acc_npbc144 = np.zeros_like(p_npbc_144)
    acc_pbc144 = np.zeros_like(p_pbc_144)

    basin_boot = {name: np.zeros((reps, 2), dtype=float) for name in basin_order}
    iact_est_frames = 1.0
    block_len_used = 1
    if args.bootstrap_mode == "block":
        block_len_used, iact_est_frames = resolve_block_len_frames(
            args.bootstrap_block_len, n_bal, npbc_phi, npbc_psi, pbc_phi, pbc_psi
        )
    elif args.bootstrap_block_len > 0:
        block_len_used = max(2, int(args.bootstrap_block_len))

    for ib in range(reps):
        idx_n = bootstrap_indices(n_npbc, n_bal, rng, args.bootstrap_mode, block_len_used)
        idx_p = bootstrap_indices(n_pbc, n_bal, rng, args.bootstrap_mode, block_len_used)
        nphi, npsi = npbc_phi[idx_n], npbc_psi[idx_n]
        pphi, ppsi = pbc_phi[idx_p], pbc_psi[idx_p]

        h_n72 = density2d_prob(nphi, npsi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        h_p72 = density2d_prob(pphi, ppsi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        h_n144 = density2d_prob(nphi, npsi, bins144, args.density_method, args.kde_bandwidth_deg, mod)
        h_p144 = density2d_prob(pphi, ppsi, bins144, args.density_method, args.kde_bandwidth_deg, mod)

        acc_npbc72 += h_n72
        acc_pbc72 += h_p72
        acc_npbc144 += h_n144
        acc_pbc144 += h_p144

        boot_jsd72[ib] = mod.js_divergence_from_prob(h_n72, h_p72)
        boot_ov72[ib] = mod.overlap_coeff_from_prob(h_n72, h_p72)
        boot_jsd144[ib] = mod.js_divergence_from_prob(h_n144, h_p144)
        boot_ov144[ib] = mod.overlap_coeff_from_prob(h_n144, h_p144)

        boot_phi_prob_npbc[ib] = density1d_prob(nphi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        boot_psi_prob_npbc[ib] = density1d_prob(npsi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        boot_phi_prob_pbc[ib] = density1d_prob(pphi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        boot_psi_prob_pbc[ib] = density1d_prob(ppsi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        if boot_p_npbc72 is not None:
            boot_p_npbc72[ib] = h_n72
            boot_p_pbc72[ib] = h_p72

        for name in basin_order:
            mask = basin_masks[name]
            if args.crop_unsampled_fes:
                h_n72_basin = masked_renormalized_prob(h_n72, crop_mask_72)
                h_p72_basin = masked_renormalized_prob(h_p72, crop_mask_72)
                basin_boot[name][ib, 0] = float(np.sum(h_n72_basin[mask]))
                basin_boot[name][ib, 1] = float(np.sum(h_p72_basin[mask]))
            else:
                basin_boot[name][ib, 0] = float(np.sum(h_n72[mask]))
                basin_boot[name][ib, 1] = float(np.sum(h_p72[mask]))

    p_npbc_bal72 = acc_npbc72 / float(reps)
    p_pbc_bal72 = acc_pbc72 / float(reps)
    p_npbc_bal144 = acc_npbc144 / float(reps)
    p_pbc_bal144 = acc_pbc144 / float(reps)

    f_npbc_bal72 = mod.free_energy_from_hist(p_npbc_bal72, args.temp_k)
    f_pbc_bal72 = mod.free_energy_from_hist(p_pbc_bal72, args.temp_k)
    f_npbc_bal144 = mod.free_energy_from_hist(p_npbc_bal144, args.temp_k)
    f_pbc_bal144 = mod.free_energy_from_hist(p_pbc_bal144, args.temp_k)

    # Time-block stability.
    block_rows = []
    for lane, phi, psi, pfull in (
        ("npbc", npbc_phi, npbc_psi, p_npbc_72),
        ("pbc", pbc_phi, pbc_psi, p_pbc_72),
    ):
        for bname, idx in mod.split_three_blocks(phi.size):
            if idx.size < 10:
                continue
            pb = density2d_prob(phi[idx], psi[idx], bins72, args.density_method, args.kde_bandwidth_deg, mod)
            block_rows.append(
                {
                    "lane": lane,
                    "block": bname,
                    "nframes": int(idx.size),
                    "jsd_to_full": mod.js_divergence_from_prob(pb, pfull),
                    "overlap_to_full": mod.overlap_coeff_from_prob(pb, pfull),
                }
            )

    mod.write_csv(
        outdir / "timeblock_stability.csv",
        [(r["lane"], r["block"], r["nframes"], r["jsd_to_full"], r["overlap_to_full"]) for r in block_rows],
        ["lane", "block", "nframes", "jsd_to_full", "overlap_to_full"],
    )

    # Basin transitions / dwell times.
    state_order = basin_order + ["other"]
    states_npbc = assign_states_from_basin_masks(npbc_phi, npbc_psi, bins72, basin_masks, basin_order)
    states_pbc = assign_states_from_basin_masks(pbc_phi, pbc_psi, bins72, basin_masks, basin_order)
    trans_npbc = mod.row_normalize(mod.transition_matrix(states_npbc, state_order))
    trans_pbc = mod.row_normalize(mod.transition_matrix(states_pbc, state_order))
    trans_jsd = mod.js_divergence_from_prob(trans_npbc, trans_pbc)
    trans_ov = mod.overlap_coeff_from_prob(trans_npbc, trans_pbc)

    npbc_dstep = np.diff(npbc_step)
    npbc_dstep = npbc_dstep[npbc_dstep > 0]
    pbc_dstep = np.diff(pbc_step)
    pbc_dstep = pbc_dstep[pbc_dstep > 0]
    npbc_frame_dt_ps = float(np.median(npbc_dstep) * 1.0e-3) if npbc_dstep.size > 0 else 0.5
    pbc_frame_dt_ps = float(np.median(pbc_dstep) * 1.0e-3) if pbc_dstep.size > 0 else 0.5

    # Always-on torsion autocorrelation analysis.
    acf_max_lag_npbc = max(1, min(int(args.acf_max_lag), n_npbc - 1))
    acf_max_lag_pbc = max(1, min(int(args.acf_max_lag), n_pbc - 1))
    acf_npbc_phi = torsion_autocorr_cosdelta(npbc_phi, acf_max_lag_npbc)
    acf_pbc_phi = torsion_autocorr_cosdelta(pbc_phi, acf_max_lag_pbc)
    acf_npbc_psi = torsion_autocorr_cosdelta(npbc_psi, acf_max_lag_npbc)
    acf_pbc_psi = torsion_autocorr_cosdelta(pbc_psi, acf_max_lag_pbc)
    acf_npbc_phi_centered, acf_npbc_phi_c, acf_npbc_phi_s = torsion_autocorr_centered_components(npbc_phi, acf_max_lag_npbc)
    acf_pbc_phi_centered, acf_pbc_phi_c, acf_pbc_phi_s = torsion_autocorr_centered_components(pbc_phi, acf_max_lag_pbc)
    acf_npbc_psi_centered, acf_npbc_psi_c, acf_npbc_psi_s = torsion_autocorr_centered_components(npbc_psi, acf_max_lag_npbc)
    acf_pbc_psi_centered, acf_pbc_psi_c, acf_pbc_psi_s = torsion_autocorr_centered_components(pbc_psi, acf_max_lag_pbc)
    iact_phi_npbc = estimate_angle_iact_frames(npbc_phi, max_lag=acf_max_lag_npbc)
    iact_phi_pbc = estimate_angle_iact_frames(pbc_phi, max_lag=acf_max_lag_pbc)
    iact_psi_npbc = estimate_angle_iact_frames(npbc_psi, max_lag=acf_max_lag_npbc)
    iact_psi_pbc = estimate_angle_iact_frames(pbc_psi, max_lag=acf_max_lag_pbc)

    nlag_common = max(int(acf_npbc_phi.size), int(acf_pbc_phi.size))
    acf_rows = []
    for k in range(nlag_common):
        acf_rows.append(
            (
                int(k),
                float(k * npbc_frame_dt_ps),
                float(k * pbc_frame_dt_ps),
                float(acf_npbc_phi[k]) if k < acf_npbc_phi.size else np.nan,
                float(acf_pbc_phi[k]) if k < acf_pbc_phi.size else np.nan,
                float(acf_npbc_psi[k]) if k < acf_npbc_psi.size else np.nan,
                float(acf_pbc_psi[k]) if k < acf_pbc_psi.size else np.nan,
            )
        )
    mod.write_csv(
        outdir / "torsion_autocorrelation.csv",
        acf_rows,
        [
            "lag_frame",
            "lag_ps_npbc",
            "lag_ps_pbc",
            "phi_acf_npbc",
            "phi_acf_pbc",
            "psi_acf_npbc",
            "psi_acf_pbc",
        ],
    )
    nlag_common_centered = max(int(acf_npbc_phi_centered.size), int(acf_pbc_phi_centered.size))
    acf_rows_centered = []
    for k in range(nlag_common_centered):
        acf_rows_centered.append(
            (
                int(k),
                float(k * npbc_frame_dt_ps),
                float(k * pbc_frame_dt_ps),
                float(acf_npbc_phi_centered[k]) if k < acf_npbc_phi_centered.size else np.nan,
                float(acf_pbc_phi_centered[k]) if k < acf_pbc_phi_centered.size else np.nan,
                float(acf_npbc_psi_centered[k]) if k < acf_npbc_psi_centered.size else np.nan,
                float(acf_pbc_psi_centered[k]) if k < acf_pbc_psi_centered.size else np.nan,
                float(acf_npbc_phi_c[k]) if k < acf_npbc_phi_c.size else np.nan,
                float(acf_npbc_phi_s[k]) if k < acf_npbc_phi_s.size else np.nan,
                float(acf_pbc_phi_c[k]) if k < acf_pbc_phi_c.size else np.nan,
                float(acf_pbc_phi_s[k]) if k < acf_pbc_phi_s.size else np.nan,
                float(acf_npbc_psi_c[k]) if k < acf_npbc_psi_c.size else np.nan,
                float(acf_npbc_psi_s[k]) if k < acf_npbc_psi_s.size else np.nan,
                float(acf_pbc_psi_c[k]) if k < acf_pbc_psi_c.size else np.nan,
                float(acf_pbc_psi_s[k]) if k < acf_pbc_psi_s.size else np.nan,
            )
        )
    mod.write_csv(
        outdir / "torsion_autocorrelation_centered.csv",
        acf_rows_centered,
        [
            "lag_frame",
            "lag_ps_npbc",
            "lag_ps_pbc",
            "phi_acf_centered_npbc",
            "phi_acf_centered_pbc",
            "psi_acf_centered_npbc",
            "psi_acf_centered_pbc",
            "phi_acf_cos_npbc",
            "phi_acf_sin_npbc",
            "phi_acf_cos_pbc",
            "phi_acf_sin_pbc",
            "psi_acf_cos_npbc",
            "psi_acf_sin_npbc",
            "psi_acf_cos_pbc",
            "psi_acf_sin_pbc",
        ],
    )

    transition_rows = []
    transition_rows.extend(transition_event_rows("npbc", states_npbc, npbc_step, t_npbc, npbc_source, npbc_frame_dt_ps))
    transition_rows.extend(transition_event_rows("pbc", states_pbc, pbc_step, t_pbc, pbc_source, pbc_frame_dt_ps))
    mod.write_csv(
        outdir / "transition_events.csv",
        transition_rows,
        [
            "lane",
            "from_frame_start",
            "from_frame_end",
            "from_source",
            "from_step",
            "from_time_ns",
            "to_source",
            "to_step",
            "to_time_ns",
            "from_state",
            "to_state",
            "from_run_len_frames",
            "from_run_len_ps",
        ],
    )

    first_entry = []
    first_entry.extend(first_entry_rows("npbc", states_npbc, npbc_step, t_npbc, npbc_source, state_order))
    first_entry.extend(first_entry_rows("pbc", states_pbc, pbc_step, t_pbc, pbc_source, state_order))
    mod.write_csv(
        outdir / "state_first_entry.csv",
        first_entry,
        ["lane", "state", "frame_index", "source", "step", "time_ns"],
    )

    mod.write_csv(
        outdir / "transition_matrix_npbc.csv",
        [(state_order[i], *trans_npbc[i, :].tolist()) for i in range(len(state_order))],
        ["from_state", *[f"to_{s}" for s in state_order]],
    )
    mod.write_csv(
        outdir / "transition_matrix_pbc.csv",
        [(state_order[i], *trans_pbc[i, :].tolist()) for i in range(len(state_order))],
        ["from_state", *[f"to_{s}" for s in state_order]],
    )

    dwell_states = ["alphaR", "alpha_prime", "C7eq", "PPII"]
    dwell_npbc = mod.dwell_times(states_npbc, dwell_states, dt_ps=0.5)
    dwell_pbc = mod.dwell_times(states_pbc, dwell_states, dt_ps=0.5)
    mod.write_csv(
        outdir / "basin_dwell_summary.csv",
        [
            (
                st,
                int(len(dwell_npbc.get(st, []))),
                float(np.mean(dwell_npbc.get(st, [np.nan])) if len(dwell_npbc.get(st, [])) > 0 else np.nan),
                float(np.median(dwell_npbc.get(st, [np.nan])) if len(dwell_npbc.get(st, [])) > 0 else np.nan),
                int(len(dwell_pbc.get(st, []))),
                float(np.mean(dwell_pbc.get(st, [np.nan])) if len(dwell_pbc.get(st, [])) > 0 else np.nan),
                float(np.median(dwell_pbc.get(st, [np.nan])) if len(dwell_pbc.get(st, [])) > 0 else np.nan),
            )
            for st in dwell_states
        ],
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

    # Basin populations.
    basin_rows = []
    pbc_blocks = mod.split_blocks(pbc_phi.size, 10)
    npbc_blocks = mod.split_blocks(npbc_phi.size, 10)
    for (name, _pc, _sc) in mod.BASIN_CENTERS:
        pbc_pop = raw_basin_pbc[name]
        npbc_pop = raw_basin_npbc[name]
        pbc_sem = mod.block_sem(
            np.array(
                [
                    float(
                        np.sum((
                            masked_renormalized_prob(
                                density2d_prob(
                                    pbc_phi[idx], pbc_psi[idx], bins72, args.density_method, args.kde_bandwidth_deg, mod
                                ),
                                crop_mask_72,
                            )
                            if args.crop_unsampled_fes
                            else density2d_prob(
                                pbc_phi[idx], pbc_psi[idx], bins72, args.density_method, args.kde_bandwidth_deg, mod
                            )
                        )[basin_masks[name]])
                    )
                    for idx in pbc_blocks
                ]
            ),
            len(pbc_blocks),
        )
        npbc_sem = mod.block_sem(
            np.array(
                [
                    float(
                        np.sum((
                            masked_renormalized_prob(
                                density2d_prob(
                                    npbc_phi[idx], npbc_psi[idx], bins72, args.density_method, args.kde_bandwidth_deg, mod
                                ),
                                crop_mask_72,
                            )
                            if args.crop_unsampled_fes
                            else density2d_prob(
                                npbc_phi[idx], npbc_psi[idx], bins72, args.density_method, args.kde_bandwidth_deg, mod
                            )
                        )[basin_masks[name]])
                    )
                    for idx in npbc_blocks
                ]
            ),
            len(npbc_blocks),
        )
        nb = basin_boot[name][:, 0]
        pb = basin_boot[name][:, 1]
        delta_bs = nb - pb
        delta_ci_lo = float(np.quantile(delta_bs, 0.025))
        delta_ci_hi = float(np.quantile(delta_bs, 0.975))
        delta_sig = (delta_ci_lo > 0.0) or (delta_ci_hi < 0.0)
        basin_rows.append(
            {
                "basin": name,
                "pbc_pop": float(pbc_pop),
                "pbc_sem": float(pbc_sem),
                "npbc_raw_pop": float(npbc_pop),
                "npbc_raw_sem": float(npbc_sem),
                "npbc_balanced_mean": float(np.mean(nb)),
                "npbc_balanced_std": float(np.std(nb, ddof=1)),
                "npbc_balanced_ci_lo": float(np.quantile(nb, 0.025)),
                "npbc_balanced_ci_hi": float(np.quantile(nb, 0.975)),
                "pbc_balanced_mean": float(np.mean(pb)),
                "pbc_balanced_ci_lo": float(np.quantile(pb, 0.025)),
                "pbc_balanced_ci_hi": float(np.quantile(pb, 0.975)),
                "delta_balanced_npbc_minus_pbc": float(np.mean(nb) - np.mean(pb)),
                "delta_balanced_ci_lo": delta_ci_lo,
                "delta_balanced_ci_hi": delta_ci_hi,
                "delta_significant": bool(delta_sig),
            }
        )

    mod.write_csv(
        outdir / "basin_populations.csv",
        [
            (
                r["basin"],
                r["pbc_pop"],
                r["pbc_sem"],
                r["npbc_raw_pop"],
                r["npbc_raw_sem"],
                r["npbc_balanced_mean"],
                r["npbc_balanced_std"],
                r["npbc_balanced_ci_lo"],
                r["npbc_balanced_ci_hi"],
                r["pbc_balanced_mean"],
                r["pbc_balanced_ci_lo"],
                r["pbc_balanced_ci_hi"],
                r["delta_balanced_npbc_minus_pbc"],
                r["delta_balanced_ci_lo"],
                r["delta_balanced_ci_hi"],
                int(r["delta_significant"]),
            )
            for r in basin_rows
        ],
        [
            "basin",
            "pbc_raw",
            "pbc_sem",
            "npbc_raw",
            "npbc_sem",
            "npbc_bal_mean",
            "npbc_bal_std",
            "npbc_bal_ci_lo",
            "npbc_bal_ci_hi",
            "pbc_bal_mean",
            "pbc_bal_ci_lo",
            "pbc_bal_ci_hi",
            "delta_bal_npbc_minus_pbc",
            "delta_bal_ci_lo",
            "delta_bal_ci_hi",
            "delta_significant_95",
        ],
    )

    mod.write_csv(
        outdir / "basin_populations_significant_only.csv",
        [
            (
                r["basin"],
                r["npbc_balanced_mean"],
                r["pbc_balanced_mean"],
                r["delta_balanced_npbc_minus_pbc"],
                r["delta_balanced_ci_lo"],
                r["delta_balanced_ci_hi"],
            )
            for r in basin_rows
            if r["delta_significant"]
        ],
        ["basin", "npbc_bal_mean", "pbc_bal_mean", "delta_bal_npbc_minus_pbc", "delta_bal_ci_lo", "delta_bal_ci_hi"],
    )

    summary = {
        "frames": {"npbc": int(n_npbc), "pbc": int(n_pbc), "balanced": int(n_bal)},
        "analysis_config": {
            "fes_bins_coarse": int(nbin_coarse),
            "fes_bins_fine": int(nbin_fine),
            "density_method": str(args.density_method),
            "kde_bandwidth_deg": float(args.kde_bandwidth_deg),
            "basin_definition": str(args.basin_definition),
            "basin_window_deg": float(args.basin_window_deg),
            "fes_smooth_sigma_bins": float(args.fes_smooth_sigma_bins),
            "bootstrap_mode": str(args.bootstrap_mode),
            "bootstrap_block_len_arg": int(args.bootstrap_block_len),
            "bootstrap_block_len_used": int(block_len_used),
            "iact_est_frames": float(iact_est_frames),
            "basin_populations_crop_aware": bool(args.crop_unsampled_fes),
        },
        "dihedral_raw": {
            "jsd_72": float(raw_jsd_72),
            "overlap_72": float(raw_ov_72),
            "jsd_144": float(raw_jsd_144),
            "overlap_144": float(raw_ov_144),
        },
        "dihedral_balanced_bootstrap": {
            "jsd_72_mean": float(np.mean(boot_jsd72)),
            "jsd_72_ci": [float(np.quantile(boot_jsd72, 0.025)), float(np.quantile(boot_jsd72, 0.975))],
            "overlap_72_mean": float(np.mean(boot_ov72)),
            "overlap_72_ci": [float(np.quantile(boot_ov72, 0.025)), float(np.quantile(boot_ov72, 0.975))],
            "jsd_144_mean": float(np.mean(boot_jsd144)),
            "jsd_144_ci": [float(np.quantile(boot_jsd144, 0.025)), float(np.quantile(boot_jsd144, 0.975))],
            "overlap_144_mean": float(np.mean(boot_ov144)),
            "overlap_144_ci": [float(np.quantile(boot_ov144, 0.025)), float(np.quantile(boot_ov144, 0.975))],
        },
        "structural": {
            "rmsd_all_mean_npbc": float(np.mean(struct_npbc["rmsd_all"])),
            "rmsd_all_mean_pbc": float(np.mean(struct_pbc["rmsd_all"])),
            "rmsd_heavy_mean_npbc": float(np.mean(struct_npbc["rmsd_heavy"])),
            "rmsd_heavy_mean_pbc": float(np.mean(struct_pbc["rmsd_heavy"])),
            "rmsd_all_ks": float(mod.ks_distance(struct_npbc["rmsd_all"], struct_pbc["rmsd_all"])),
            "rmsd_heavy_ks": float(mod.ks_distance(struct_npbc["rmsd_heavy"], struct_pbc["rmsd_heavy"])),
            "rgyr_all_mean_npbc": float(np.mean(struct_npbc["rgyr_all"])),
            "rgyr_all_mean_pbc": float(np.mean(struct_pbc["rgyr_all"])),
            "rgyr_heavy_mean_npbc": float(np.mean(struct_npbc["rgyr_heavy"])),
            "rgyr_heavy_mean_pbc": float(np.mean(struct_pbc["rgyr_heavy"])),
            "dist_5_17_mean_npbc": float(np.nanmean(struct_npbc["dist_5_17"])),
            "dist_5_17_mean_pbc": float(np.nanmean(struct_pbc["dist_5_17"])),
            "rmsf_delta_abs_mean": float(np.mean(np.abs(struct_npbc["rmsf"] - struct_pbc["rmsf"]))),
            "rmsf_delta_abs_max": float(np.max(np.abs(struct_npbc["rmsf"] - struct_pbc["rmsf"]))),
        },
        "transition": {"jsd": float(trans_jsd), "overlap": float(trans_ov)},
        "transition_events": {
            "npbc_n": int(sum(1 for r in transition_rows if r[0] == "npbc")),
            "pbc_n": int(sum(1 for r in transition_rows if r[0] == "pbc")),
        },
        "torsion_autocorr": {
            "acf_max_lag_frames_npbc": int(acf_max_lag_npbc),
            "acf_max_lag_frames_pbc": int(acf_max_lag_pbc),
            "frame_dt_ps_npbc": float(npbc_frame_dt_ps),
            "frame_dt_ps_pbc": float(pbc_frame_dt_ps),
            "iact_phi_frames_npbc": float(iact_phi_npbc),
            "iact_phi_frames_pbc": float(iact_phi_pbc),
            "iact_psi_frames_npbc": float(iact_psi_npbc),
            "iact_psi_frames_pbc": float(iact_psi_pbc),
            "iact_phi_ps_npbc": float(iact_phi_npbc * npbc_frame_dt_ps),
            "iact_phi_ps_pbc": float(iact_phi_pbc * pbc_frame_dt_ps),
            "iact_psi_ps_npbc": float(iact_psi_npbc * npbc_frame_dt_ps),
            "iact_psi_ps_pbc": float(iact_psi_pbc * pbc_frame_dt_ps),
            "centered_curve_csv": "torsion_autocorrelation_centered.csv",
            "cosdelta_curve_csv": "torsion_autocorrelation.csv",
        },
        "basin_significance": {
            "n_significant_95": int(sum(1 for r in basin_rows if r["delta_significant"])),
            "n_total": int(len(basin_rows)),
        },
    }
    (outdir / "metrics_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Plot suite.
    crop_stats = {}
    crop_stats["raw_72"] = maybe_plot_fes_triptych_cropped(
        mod,
        outdir / f"01_raw_fes_{nbin_coarse}x{nbin_coarse}.png",
        f_npbc_72,
        f_pbc_72,
        p_npbc_72,
        p_pbc_72,
        f"Raw FES {nbin_coarse}x{nbin_coarse}",
        args.crop_unsampled_fes,
        args.crop_min_prob,
        args.crop_policy,
    )
    crop_stats["raw_144"] = maybe_plot_fes_triptych_cropped(
        mod,
        outdir / f"02_raw_fes_{nbin_fine}x{nbin_fine}.png",
        f_npbc_144,
        f_pbc_144,
        p_npbc_144,
        p_pbc_144,
        f"Raw FES {nbin_fine}x{nbin_fine}",
        args.crop_unsampled_fes,
        args.crop_min_prob,
        args.crop_policy,
    )
    crop_stats["balanced_72"] = maybe_plot_fes_triptych_cropped(
        mod,
        outdir / f"03_balanced_fes_{nbin_coarse}x{nbin_coarse}.png",
        f_npbc_bal72,
        f_pbc_bal72,
        p_npbc_bal72,
        p_pbc_bal72,
        f"Balanced FES {nbin_coarse}x{nbin_coarse}",
        args.crop_unsampled_fes,
        args.crop_min_prob,
        args.crop_policy,
    )
    crop_stats["balanced_144"] = maybe_plot_fes_triptych_cropped(
        mod,
        outdir / f"04_balanced_fes_{nbin_fine}x{nbin_fine}.png",
        f_npbc_bal144,
        f_pbc_bal144,
        p_npbc_bal144,
        p_pbc_bal144,
        f"Balanced FES {nbin_fine}x{nbin_fine}",
        args.crop_unsampled_fes,
        args.crop_min_prob,
        args.crop_policy,
    )
    mod.maybe_plot_bootstrap_metrics(outdir / "05_bootstrap_metric_distributions.png", boot_jsd72, boot_ov72, raw_jsd_72, raw_ov_72)
    mod.maybe_plot_basin_populations(outdir / "06_basin_populations.png", basin_rows)
    mod.maybe_plot_timeblock_stability(outdir / "07_timeblock_stability.png", block_rows)

    # Marginals with CI for BOTH lanes.
    if mod.HAVE_PLOT:
        import matplotlib.pyplot as plt

        mids = 0.5 * (bins72[:-1] + bins72[1:])
        npbc_phi_raw = density1d_prob(npbc_phi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        pbc_phi_raw = density1d_prob(pbc_phi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        npbc_psi_raw = density1d_prob(npbc_psi, bins72, args.density_method, args.kde_bandwidth_deg, mod)
        pbc_psi_raw = density1d_prob(pbc_psi, bins72, args.density_method, args.kde_bandwidth_deg, mod)

        # Crop-aware marginal mode: derive marginals from the same masked 2D map used by cropped FES.
        if args.crop_unsampled_fes:
            crop_mask_2d = build_crop_mask(p_npbc_72, p_pbc_72, True, args.crop_min_prob, args.crop_policy)
            p_npbc_crop = masked_renormalized_prob(p_npbc_72, crop_mask_2d)
            p_pbc_crop = masked_renormalized_prob(p_pbc_72, crop_mask_2d)
            npbc_phi_raw = np.sum(p_npbc_crop, axis=1)
            npbc_psi_raw = np.sum(p_npbc_crop, axis=0)
            pbc_phi_raw = np.sum(p_pbc_crop, axis=1)
            pbc_psi_raw = np.sum(p_pbc_crop, axis=0)

            if boot_p_npbc72 is not None:
                boot_phi_prob_npbc_crop = np.zeros_like(boot_phi_prob_npbc)
                boot_psi_prob_npbc_crop = np.zeros_like(boot_psi_prob_npbc)
                boot_phi_prob_pbc_crop = np.zeros_like(boot_phi_prob_pbc)
                boot_psi_prob_pbc_crop = np.zeros_like(boot_psi_prob_pbc)
                for ib in range(reps):
                    pn = masked_renormalized_prob(boot_p_npbc72[ib], crop_mask_2d)
                    pp = masked_renormalized_prob(boot_p_pbc72[ib], crop_mask_2d)
                    boot_phi_prob_npbc_crop[ib] = np.sum(pn, axis=1)
                    boot_psi_prob_npbc_crop[ib] = np.sum(pn, axis=0)
                    boot_phi_prob_pbc_crop[ib] = np.sum(pp, axis=1)
                    boot_psi_prob_pbc_crop[ib] = np.sum(pp, axis=0)
                boot_phi_prob_npbc = boot_phi_prob_npbc_crop
                boot_psi_prob_npbc = boot_psi_prob_npbc_crop
                boot_phi_prob_pbc = boot_phi_prob_pbc_crop
                boot_psi_prob_pbc = boot_psi_prob_pbc_crop

        npbc_phi_lo = np.quantile(boot_phi_prob_npbc, 0.025, axis=0)
        npbc_phi_hi = np.quantile(boot_phi_prob_npbc, 0.975, axis=0)
        pbc_phi_lo = np.quantile(boot_phi_prob_pbc, 0.025, axis=0)
        pbc_phi_hi = np.quantile(boot_phi_prob_pbc, 0.975, axis=0)
        npbc_psi_lo = np.quantile(boot_psi_prob_npbc, 0.025, axis=0)
        npbc_psi_hi = np.quantile(boot_psi_prob_npbc, 0.975, axis=0)
        pbc_psi_lo = np.quantile(boot_psi_prob_pbc, 0.025, axis=0)
        pbc_psi_hi = np.quantile(boot_psi_prob_pbc, 0.975, axis=0)
        marg_mode = "crop-aware" if args.crop_unsampled_fes else "raw"

        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
        axs[0].plot(mids, npbc_phi_raw, label="NPBC", lw=1.2)
        axs[0].fill_between(mids, npbc_phi_lo, npbc_phi_hi, alpha=0.2)
        axs[0].plot(mids, pbc_phi_raw, label="PBC", lw=1.2)
        axs[0].fill_between(mids, pbc_phi_lo, pbc_phi_hi, alpha=0.2)
        axs[0].set_title(f"phi marginal ({marg_mode} + balanced 95% CI)")
        axs[0].set_xlabel("phi (deg)")
        axs[0].set_ylabel("probability")
        axs[0].grid(alpha=0.25)
        axs[0].legend()

        axs[1].plot(mids, npbc_psi_raw, label="NPBC", lw=1.2)
        axs[1].fill_between(mids, npbc_psi_lo, npbc_psi_hi, alpha=0.2)
        axs[1].plot(mids, pbc_psi_raw, label="PBC", lw=1.2)
        axs[1].fill_between(mids, pbc_psi_lo, pbc_psi_hi, alpha=0.2)
        axs[1].set_title(f"psi marginal ({marg_mode} + balanced 95% CI)")
        axs[1].set_xlabel("psi (deg)")
        axs[1].set_ylabel("probability")
        axs[1].grid(alpha=0.25)
        axs[1].legend()
        fig.savefig(outdir / "08_phi_psi_marginals_dual_ci.png", dpi=170)
        plt.close(fig)

    mod.maybe_plot_rmsd_timeseries(
        outdir / "09_rmsd_timeseries.png",
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
    mod.maybe_plot_rmsd_distribution(
        outdir / "10_rmsd_distributions.png",
        struct_npbc["rmsd_all"],
        struct_pbc["rmsd_all"],
        struct_npbc["rmsd_heavy"],
        struct_pbc["rmsd_heavy"],
    )
    mod.maybe_plot_rmsf_per_atom(outdir / "11_rmsf_per_atom.png", atom_ids, atom_types, struct_npbc["rmsf"], struct_pbc["rmsf"])
    mod.maybe_plot_rgyr_timeseries(
        outdir / "12_rgyr_timeseries.png",
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
    mod.maybe_plot_distance_timeseries(outdir / "13_dist_5_17_timeseries.png", t_npbc, t_pbc, struct_npbc["dist_5_17"], struct_pbc["dist_5_17"])
    mod.maybe_plot_transition_matrix(outdir / "14_transition_matrix_npbc.png", trans_npbc, state_order, "NPBC basin transition matrix")
    mod.maybe_plot_transition_matrix(outdir / "15_transition_matrix_pbc.png", trans_pbc, state_order, "PBC basin transition matrix")
    mod.maybe_plot_dwell_times(outdir / "16_dwell_distributions.png", dwell_npbc, dwell_pbc, dwell_states)
    maybe_plot_transition_events_timeline(
        outdir / "17_transition_events_timeline.png",
        transition_rows,
        state_order,
        first_entry,
        mod.HAVE_PLOT,
    )
    maybe_plot_transition_event_cumulative(
        outdir / "18_transition_events_cumulative.png",
        transition_rows,
        mod.HAVE_PLOT,
    )
    maybe_plot_torsion_autocorr(
        outdir / "19_torsion_autocorrelation.png",
        np.arange(acf_npbc_phi.size, dtype=float) * npbc_frame_dt_ps,
        np.arange(acf_pbc_phi.size, dtype=float) * pbc_frame_dt_ps,
        acf_npbc_phi,
        acf_pbc_phi,
        acf_npbc_psi,
        acf_pbc_psi,
        mod.HAVE_PLOT,
    )
    maybe_plot_torsion_autocorr_centered(
        outdir / "20_torsion_autocorrelation_centered.png",
        np.arange(acf_npbc_phi_centered.size, dtype=float) * npbc_frame_dt_ps,
        np.arange(acf_pbc_phi_centered.size, dtype=float) * pbc_frame_dt_ps,
        acf_npbc_phi_centered,
        acf_pbc_phi_centered,
        acf_npbc_psi_centered,
        acf_pbc_psi_centered,
        mod.HAVE_PLOT,
    )
    (outdir / "fes_crop_stats.json").write_text(json.dumps(crop_stats, indent=2) + "\n")

    report = []
    report.append("# Full Solute Dynamics Comparison (NPBC eq+prod vs PBC eq+prod)")
    report.append("")
    report.append(f"- Frames NPBC/PBC: **{n_npbc} / {n_pbc}** (balanced N={n_bal})")
    report.append(f"- Grid coarse/fine: **{nbin_coarse}x{nbin_coarse} / {nbin_fine}x{nbin_fine}**")
    report.append(f"- Raw JSD(coarse)/Overlap(coarse): **{raw_jsd_72:.5f} / {raw_ov_72:.5f}**")
    report.append(
        f"- Balanced JSD(coarse) mean [95% CI]: **{np.mean(boot_jsd72):.5f} [{np.quantile(boot_jsd72,0.025):.5f}, {np.quantile(boot_jsd72,0.975):.5f}]**"
    )
    report.append(
        f"- Balanced Overlap(coarse) mean [95% CI]: **{np.mean(boot_ov72):.5f} [{np.quantile(boot_ov72,0.025):.5f}, {np.quantile(boot_ov72,0.975):.5f}]**"
    )
    report.append("")
    report.append("## Analysis settings")
    report.append(f"- density_method: `{args.density_method}`")
    report.append(f"- kde_bandwidth_deg: `{args.kde_bandwidth_deg:g}`")
    report.append(f"- basin_definition: `{args.basin_definition}`")
    if args.basin_definition == "fixed":
        report.append(f"- basin_window_deg: `{args.basin_window_deg:g}`")
    else:
        report.append(f"- fes_smooth_sigma_bins: `{args.fes_smooth_sigma_bins:g}`")
    report.append(
        f"- bootstrap_mode: `{args.bootstrap_mode}` (block_len_used={summary['analysis_config']['bootstrap_block_len_used']}, iact_est={summary['analysis_config']['iact_est_frames']:.2f})"
    )
    report.append("")
    report.append("## FES cropping")
    report.append(f"- crop_unsampled_fes: `{args.crop_unsampled_fes}`")
    if args.crop_unsampled_fes:
        report.append(f"- crop_policy: `{args.crop_policy}`")
        report.append(f"- crop_min_prob: `{args.crop_min_prob:g}`")
        report.append(f"- masked_fraction_raw72: `{crop_stats['raw_72'].get('masked_fraction', 0.0):.4f}`")
        report.append(f"- masked_fraction_raw144: `{crop_stats['raw_144'].get('masked_fraction', 0.0):.4f}`")
    report.append("")
    report.append("## Torsion autocorrelation")
    report.append(
        f"- phi IACT NPBC/PBC: {summary['torsion_autocorr']['iact_phi_frames_npbc']:.2f} / {summary['torsion_autocorr']['iact_phi_frames_pbc']:.2f} frames"
    )
    report.append(
        f"- psi IACT NPBC/PBC: {summary['torsion_autocorr']['iact_psi_frames_npbc']:.2f} / {summary['torsion_autocorr']['iact_psi_frames_pbc']:.2f} frames"
    )
    report.append(
        f"- phi IACT NPBC/PBC: {summary['torsion_autocorr']['iact_phi_ps_npbc']:.2f} / {summary['torsion_autocorr']['iact_phi_ps_pbc']:.2f} ps"
    )
    report.append(
        f"- psi IACT NPBC/PBC: {summary['torsion_autocorr']['iact_psi_ps_npbc']:.2f} / {summary['torsion_autocorr']['iact_psi_ps_pbc']:.2f} ps"
    )
    report.append(
        "- files: `torsion_autocorrelation.csv`, `torsion_autocorrelation_centered.csv`, `19_torsion_autocorrelation.png`, `20_torsion_autocorrelation_centered.png`"
    )
    report.append("")
    report.append("## Structural metrics")
    report.append(
        f"- RMSD(all) mean NPBC/PBC: {summary['structural']['rmsd_all_mean_npbc']:.3f} / {summary['structural']['rmsd_all_mean_pbc']:.3f} A"
    )
    report.append(
        f"- RMSD(heavy) mean NPBC/PBC: {summary['structural']['rmsd_heavy_mean_npbc']:.3f} / {summary['structural']['rmsd_heavy_mean_pbc']:.3f} A"
    )
    report.append(
        f"- RMSF |delta| mean/max: {summary['structural']['rmsf_delta_abs_mean']:.3f} / {summary['structural']['rmsf_delta_abs_max']:.3f} A"
    )
    report.append(
        f"- Rg(all) mean NPBC/PBC: {summary['structural']['rgyr_all_mean_npbc']:.3f} / {summary['structural']['rgyr_all_mean_pbc']:.3f} A"
    )
    report.append(f"- Transition JSD/Overlap: {trans_jsd:.5f} / {trans_ov:.5f}")
    report.append(
        f"- Transition events NPBC/PBC: {summary['transition_events']['npbc_n']} / {summary['transition_events']['pbc_n']} (see transition_events.csv, state_first_entry.csv)"
    )
    report.append(
        f"- Significant basin deltas (95% CI excludes 0): {summary['basin_significance']['n_significant_95']} / {summary['basin_significance']['n_total']}"
    )
    (outdir / "report.md").write_text("\n".join(report) + "\n")

    print(outdir)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
