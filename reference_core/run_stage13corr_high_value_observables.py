#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import sys
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


def load_analysis_module(module_path: Path):
    spec = importlib.util.spec_from_file_location("stage13_dih", module_path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["stage13_dih"] = mod
    spec.loader.exec_module(mod)
    return mod


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Extended high-value NPBC/PBC observables for stage13corr (torsion ACF, hydration shells, RDF/CN, H-bonds)."
    )
    ap.add_argument(
        "--stage-root",
        required=True,
        help="Directory containing the selected NPBC/PBC dump files.",
    )
    ap.add_argument(
        "--analysis-module",
        required=True,
        help="Path to analyze_stage13_dihedrals.py.",
    )
    ap.add_argument(
        "--npbc-eq-dump",
        required=True,
    )
    ap.add_argument(
        "--npbc-prod-dump",
        required=True,
    )
    ap.add_argument(
        "--pbc-eq-dump",
        required=True,
    )
    ap.add_argument(
        "--pbc-prod-dump",
        required=True,
    )
    ap.add_argument("--frame-stride", type=int, default=1)
    ap.add_argument("--rmax", type=float, default=8.0)
    ap.add_argument("--dr", type=float, default=0.1)
    ap.add_argument("--coordination-cutoff", type=float, default=3.5)
    ap.add_argument("--shell1-cutoff", type=float, default=3.5)
    ap.add_argument("--shell2-cutoff", type=float, default=5.0)
    ap.add_argument("--hbond-distance-cutoff", type=float, default=2.5)
    ap.add_argument("--hbond-angle-cutoff-deg", type=float, default=150.0)
    ap.add_argument("--donor-bond-cutoff", type=float, default=1.35)
    ap.add_argument("--torsion-max-lag", type=int, default=400)
    ap.add_argument("--outdir", default=None)
    return ap.parse_args()


def _must_exist(path: Path) -> None:
    if not path.exists():
        raise FileNotFoundError(f"Missing file: {path}")


def write_csv(path: Path, rows: list[tuple], header: list[str]) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for r in rows:
            w.writerow(r)


def _is_periodic_flags(boundary_flags: tuple[str, str, str]) -> np.ndarray:
    return np.array([str(f).lower().startswith("p") for f in boundary_flags], dtype=bool)


def _minimum_image(delta: np.ndarray, box_len: np.ndarray, periodic_mask: np.ndarray) -> np.ndarray:
    out = np.asarray(delta, dtype=float).copy()
    for d in range(3):
        if periodic_mask[d] and box_len[d] > 0.0:
            out[..., d] -= np.round(out[..., d] / box_len[d]) * box_len[d]
    return out


def _angle_deg(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    num = np.sum(u * v, axis=-1)
    den = np.linalg.norm(u, axis=-1) * np.linalg.norm(v, axis=-1)
    den = np.where(den > 1.0e-12, den, np.nan)
    c = np.clip(num / den, -1.0, 1.0)
    return np.degrees(np.arccos(c))


def _iter_dump_frames(path: Path):
    with path.open("r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            step = int(f.readline().strip())
            f.readline()  # ITEM: NUMBER OF ATOMS
            natoms = int(f.readline().strip())
            bounds_header = f.readline().strip().split()
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
            for key in ("id", "mol", "type", "x", "y", "z"):
                if key not in idx:
                    raise RuntimeError(f"{path} missing required dump column: {key}")

            arr = np.zeros((natoms, 6), dtype=float)
            for i in range(natoms):
                cols = f.readline().split()
                arr[i, 0] = float(cols[idx["id"]])
                arr[i, 1] = float(cols[idx["mol"]])
                arr[i, 2] = float(cols[idx["type"]])
                arr[i, 3] = float(cols[idx["x"]])
                arr[i, 4] = float(cols[idx["y"]])
                arr[i, 5] = float(cols[idx["z"]])
            yield step, box_len, boundary_flags, arr


def _update_runs(current: dict, finished: list, active_pairs: set) -> None:
    for key in list(current.keys()):
        if key in active_pairs:
            current[key] += 1
        else:
            finished.append(current.pop(key))
    for key in active_pairs:
        if key not in current:
            current[key] = 1


def _circular_autocorr(theta_deg: np.ndarray, max_lag: int) -> np.ndarray:
    theta = np.deg2rad(np.asarray(theta_deg, dtype=float))
    n = theta.size
    if n == 0:
        return np.zeros(1, dtype=float)
    max_lag = int(max(1, min(max_lag, n - 1)))
    out = np.ones(max_lag + 1, dtype=float)
    for lag in range(1, max_lag + 1):
        d = theta[lag:] - theta[:-lag]
        out[lag] = float(np.mean(np.cos(d)))
    return out


def _iact_from_acf(acf: np.ndarray) -> float:
    if acf.size == 0:
        return float("nan")
    cut = acf.size
    for i in range(1, acf.size):
        if acf[i] < 0.0:
            cut = i
            break
    return float(0.5 + np.sum(acf[1:cut]))


@dataclass
class LaneExtended:
    lane: str
    phi: np.ndarray
    psi: np.ndarray
    step: np.ndarray
    acf_phi: np.ndarray
    acf_psi: np.ndarray
    iact_phi_frames: float
    iact_psi_frames: float
    rdf_centers: np.ndarray
    rdf: dict
    shell_timeseries: list
    hbond_timeseries: list
    hbond_lifetimes_ps: dict
    frame_dt_ps: float


def compute_lane_extended(
    lane: str,
    dump_paths: list[Path],
    mod,
    required_ids: list[int],
    frame_stride: int,
    rmax: float,
    dr: float,
    coordination_cutoff: float,
    shell1_cutoff: float,
    shell2_cutoff: float,
    hb_dist_cut: float,
    hb_angle_cut_deg: float,
    donor_bond_cutoff: float,
    torsion_max_lag: int,
) -> LaneExtended:
    # dihedrals from robust parser
    lane_parts = []
    for i, p in enumerate(dump_paths):
        parsed = mod.parse_dump_lane(p, required_ids, f"{lane}_{i}")
        lane_parts.append(parsed)

    phi = np.concatenate([p.phi for p in lane_parts])
    psi = np.concatenate([p.psi for p in lane_parts])
    step = np.concatenate([p.step for p in lane_parts])

    acf_phi = _circular_autocorr(phi, torsion_max_lag)
    acf_psi = _circular_autocorr(psi, torsion_max_lag)
    iact_phi = _iact_from_acf(acf_phi)
    iact_psi = _iact_from_acf(acf_psi)

    # hydration/H-bond metrics
    frame_stride = max(1, int(frame_stride))
    edges = np.arange(0.0, float(rmax) + float(dr) + 1.0e-12, float(dr))
    centers = 0.5 * (edges[:-1] + edges[1:])
    shell_vol = (4.0 / 3.0) * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)

    rdf_counts = {}
    rdf_sites = {}
    shell_rows = []
    hb_rows = []
    steps = []

    donor_map = None
    acc_ids = None
    cat_ids = None
    cat_labels = ("all_heavy", "O_acceptor", "N_site")

    hb_runs_s2w = {}
    hb_runs_w2s = {}
    hb_finished_s2w = []
    hb_finished_w2s = []

    frame_idx = -1
    kept_idx = -1
    for path in dump_paths:
        for st, box_len, boundary_flags, arr in _iter_dump_frames(path):
            frame_idx += 1
            if frame_idx % frame_stride != 0:
                continue
            kept_idx += 1

            ids = arr[:, 0].astype(int)
            mol = arr[:, 1].astype(int)
            typ = arr[:, 2].astype(int)
            xyz = arr[:, 3:6]
            periodic = _is_periodic_flags(boundary_flags)

            sol_mask = mol == 1
            if not np.any(sol_mask):
                continue
            sol_ids = ids[sol_mask]
            sol_types = typ[sol_mask]
            sol_xyz = xyz[sol_mask]
            sol_map = {int(sol_ids[i]): sol_xyz[i] for i in range(sol_ids.size)}

            if donor_map is None:
                heavy_idx = np.where(sol_types != 4)[0]
                h_idx = np.where(sol_types == 4)[0]
                donor_map = {}
                for ih in h_idx:
                    if heavy_idx.size == 0:
                        continue
                    d = sol_xyz[heavy_idx] - sol_xyz[ih]
                    d = _minimum_image(d, box_len, periodic)
                    rr = np.linalg.norm(d, axis=1)
                    j = int(np.argmin(rr))
                    if rr[j] <= donor_bond_cutoff and int(sol_types[heavy_idx[j]]) in (2, 3):
                        donor_map[int(sol_ids[ih])] = int(sol_ids[heavy_idx[j]])

                acc_ids = [int(sol_ids[i]) for i in np.where(sol_types == 3)[0]]
                n_ids = [int(sol_ids[i]) for i in np.where(sol_types == 2)[0]]
                heavy_ids = [int(sol_ids[i]) for i in np.where(sol_types != 4)[0]]
                if len(acc_ids) == 0:
                    acc_ids = heavy_ids[:]
                if len(n_ids) == 0:
                    n_ids = heavy_ids[:]
                cat_ids = {
                    "all_heavy": heavy_ids,
                    "O_acceptor": acc_ids,
                    "N_site": n_ids,
                }
                for c in cat_labels:
                    rdf_counts[c] = np.zeros(edges.size - 1, dtype=float)
                    rdf_sites[c] = 0

            o_mask = (mol > 1) & (typ == 3)
            h_mask = (mol > 1) & (typ == 4)
            o_xyz = xyz[o_mask]
            o_mol = mol[o_mask].astype(int)
            h_xyz = xyz[h_mask]
            h_mol = mol[h_mask].astype(int)
            if o_xyz.shape[0] == 0:
                continue

            steps.append(int(st))

            # radial densities around selected site categories
            for c in cat_labels:
                site_positions = [sol_map[sid] for sid in cat_ids[c] if sid in sol_map]
                if not site_positions:
                    continue
                sp = np.asarray(site_positions, dtype=float)
                for pos in sp:
                    d = _minimum_image(o_xyz - pos, box_len, periodic)
                    r = np.linalg.norm(d, axis=1)
                    h, _ = np.histogram(r, bins=edges)
                    rdf_counts[c] += h.astype(float)
                rdf_sites[c] += sp.shape[0]

            heavy_positions = np.asarray([sol_map[sid] for sid in cat_ids["all_heavy"] if sid in sol_map], dtype=float)
            d_hw = _minimum_image(o_xyz[:, None, :] - heavy_positions[None, :, :], box_len, periodic)
            rmin_hw = np.linalg.norm(d_hw, axis=2).min(axis=1)
            shell1_n = int(np.sum(rmin_hw < shell1_cutoff))
            shell2_n = int(np.sum((rmin_hw >= shell1_cutoff) & (rmin_hw < shell2_cutoff)))
            shell_rows.append((kept_idx, int(st), 0.0, shell1_n, shell2_n))

            water = {}
            for i in range(o_xyz.shape[0]):
                water[int(o_mol[i])] = {"O": o_xyz[i], "H": []}
            for i in range(h_xyz.shape[0]):
                mm = int(h_mol[i])
                if mm in water:
                    water[mm]["H"].append(h_xyz[i])

            acc_positions = np.asarray([sol_map[aid] for aid in acc_ids if aid in sol_map], dtype=float)
            acc_list = [aid for aid in acc_ids if aid in sol_map]

            hb_s2w_pairs = set()
            hb_w2s_pairs = set()

            # solute donor -> water acceptor
            for h_id, d_id in donor_map.items():
                if h_id not in sol_map or d_id not in sol_map:
                    continue
                hpos = sol_map[h_id]
                dpos = sol_map[d_id]
                v_dh = _minimum_image(dpos - hpos, box_len, periodic)
                d_ha = _minimum_image(o_xyz - hpos, box_len, periodic)
                r_ha = np.linalg.norm(d_ha, axis=1)
                cand = np.where(r_ha <= hb_dist_cut)[0]
                if cand.size == 0:
                    continue
                ang = _angle_deg(np.repeat(v_dh[None, :], cand.size, axis=0), d_ha[cand])
                good = cand[np.where(ang >= hb_angle_cut_deg)[0]]
                for j in good:
                    hb_s2w_pairs.add((int(h_id), int(o_mol[j])))

            # water donor -> solute acceptor
            if acc_positions.shape[0] > 0:
                for wm, wd in water.items():
                    if "O" not in wd or len(wd["H"]) == 0:
                        continue
                    opos = wd["O"]
                    for hpos in wd["H"]:
                        v_oh = _minimum_image(opos - hpos, box_len, periodic)
                        d_ha = _minimum_image(acc_positions - hpos, box_len, periodic)
                        r_ha = np.linalg.norm(d_ha, axis=1)
                        cand = np.where(r_ha <= hb_dist_cut)[0]
                        if cand.size == 0:
                            continue
                        ang = _angle_deg(np.repeat(v_oh[None, :], cand.size, axis=0), d_ha[cand])
                        valid = np.where(ang >= hb_angle_cut_deg)[0]
                        for k in valid:
                            acc_id = int(acc_list[int(cand[k])])
                            hb_w2s_pairs.add((int(wm), acc_id))

            _update_runs(hb_runs_s2w, hb_finished_s2w, hb_s2w_pairs)
            _update_runs(hb_runs_w2s, hb_finished_w2s, hb_w2s_pairs)

            hb_rows.append(
                (
                    kept_idx,
                    int(st),
                    0.0,
                    int(len(hb_s2w_pairs)),
                    int(len(hb_w2s_pairs)),
                    int(len(hb_s2w_pairs) + len(hb_w2s_pairs)),
                )
            )

    for k in list(hb_runs_s2w.keys()):
        hb_finished_s2w.append(hb_runs_s2w[k])
    for k in list(hb_runs_w2s.keys()):
        hb_finished_w2s.append(hb_runs_w2s[k])

    if len(steps) > 1:
        frame_dt_ps = float(np.median(np.diff(np.asarray(steps, dtype=float))) * 0.001)
    else:
        frame_dt_ps = 0.5
    if not np.isfinite(frame_dt_ps) or frame_dt_ps <= 0.0:
        frame_dt_ps = 0.5

    shell_rows = [(r[0], r[1], float(r[0] * frame_dt_ps * 1.0e-3), r[3], r[4]) for r in shell_rows]
    hb_rows = [(r[0], r[1], float(r[0] * frame_dt_ps * 1.0e-3), r[3], r[4], r[5]) for r in hb_rows]

    rdf = {}
    for c in cat_labels:
        if rdf_sites[c] <= 0:
            rho = np.full_like(centers, np.nan, dtype=float)
        else:
            rho = rdf_counts[c] / (float(rdf_sites[c]) * shell_vol)
        cn = float(np.sum(rho[centers <= coordination_cutoff] * shell_vol[centers <= coordination_cutoff]))
        rdf[c] = {"rho": rho, "cn": cn}

    hb_s2w_ps = np.asarray(hb_finished_s2w, dtype=float) * frame_dt_ps
    hb_w2s_ps = np.asarray(hb_finished_w2s, dtype=float) * frame_dt_ps
    hb_all_ps = np.concatenate([hb_s2w_ps, hb_w2s_ps]) if hb_s2w_ps.size + hb_w2s_ps.size > 0 else np.array([], dtype=float)

    return LaneExtended(
        lane=lane,
        phi=phi,
        psi=psi,
        step=step,
        acf_phi=acf_phi,
        acf_psi=acf_psi,
        iact_phi_frames=iact_phi,
        iact_psi_frames=iact_psi,
        rdf_centers=centers,
        rdf=rdf,
        shell_timeseries=shell_rows,
        hbond_timeseries=hb_rows,
        hbond_lifetimes_ps={"solute_to_water": hb_s2w_ps, "water_to_solute": hb_w2s_ps, "all": hb_all_ps},
        frame_dt_ps=frame_dt_ps,
    )


def _rolling(x: np.ndarray, w: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    if x.size == 0:
        return x
    w = max(1, min(int(w), x.size))
    out = np.empty_like(x)
    c = np.cumsum(np.insert(x, 0, 0.0))
    for i in range(x.size):
        lo = max(0, i - w + 1)
        out[i] = (c[i + 1] - c[lo]) / float(i - lo + 1)
    return out


def main() -> None:
    args = parse_args()

    stage_root = Path(args.stage_root).resolve()
    mod = load_analysis_module(Path(args.analysis_module).resolve())

    npbc_paths = [(stage_root / args.npbc_eq_dump).resolve(), (stage_root / args.npbc_prod_dump).resolve()]
    pbc_paths = [(stage_root / args.pbc_eq_dump).resolve(), (stage_root / args.pbc_prod_dump).resolve()]
    for p in [*npbc_paths, *pbc_paths]:
        _must_exist(p)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    if args.outdir:
        outdir = Path(args.outdir).resolve()
    else:
        outdir = (stage_root / "solute_dynamics_comparison" / f"{stamp}_high_value_observables").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    req = [5, 7, 9, 15, 17]

    npbc = compute_lane_extended(
        lane="NPBC",
        dump_paths=npbc_paths,
        mod=mod,
        required_ids=req,
        frame_stride=args.frame_stride,
        rmax=args.rmax,
        dr=args.dr,
        coordination_cutoff=args.coordination_cutoff,
        shell1_cutoff=args.shell1_cutoff,
        shell2_cutoff=args.shell2_cutoff,
        hb_dist_cut=args.hbond_distance_cutoff,
        hb_angle_cut_deg=args.hbond_angle_cutoff_deg,
        donor_bond_cutoff=args.donor_bond_cutoff,
        torsion_max_lag=args.torsion_max_lag,
    )
    pbc = compute_lane_extended(
        lane="PBC",
        dump_paths=pbc_paths,
        mod=mod,
        required_ids=req,
        frame_stride=args.frame_stride,
        rmax=args.rmax,
        dr=args.dr,
        coordination_cutoff=args.coordination_cutoff,
        shell1_cutoff=args.shell1_cutoff,
        shell2_cutoff=args.shell2_cutoff,
        hb_dist_cut=args.hbond_distance_cutoff,
        hb_angle_cut_deg=args.hbond_angle_cutoff_deg,
        donor_bond_cutoff=args.donor_bond_cutoff,
        torsion_max_lag=args.torsion_max_lag,
    )

    # Torsion ACF CSV
    max_lag = min(npbc.acf_phi.size, pbc.acf_phi.size, npbc.acf_psi.size, pbc.acf_psi.size)
    acf_rows = []
    for lag in range(max_lag):
        acf_rows.append(
            (
                lag,
                lag * npbc.frame_dt_ps,
                float(npbc.acf_phi[lag]),
                float(pbc.acf_phi[lag]),
                float(npbc.acf_psi[lag]),
                float(pbc.acf_psi[lag]),
            )
        )
    write_csv(
        outdir / "torsion_autocorrelation.csv",
        acf_rows,
        ["lag_frame", "lag_ps", "acf_phi_npbc", "acf_phi_pbc", "acf_psi_npbc", "acf_psi_pbc"],
    )

    # Shell occupancy timeseries
    write_csv(
        outdir / "hydration_shell_occupancy_npbc.csv",
        npbc.shell_timeseries,
        ["frame", "step", "time_ns", "n_shell1", "n_shell2"],
    )
    write_csv(
        outdir / "hydration_shell_occupancy_pbc.csv",
        pbc.shell_timeseries,
        ["frame", "step", "time_ns", "n_shell1", "n_shell2"],
    )

    # HBond timeseries and lifetimes
    write_csv(
        outdir / "hbond_timeseries_npbc.csv",
        npbc.hbond_timeseries,
        ["frame", "step", "time_ns", "n_solute_to_water", "n_water_to_solute", "n_total"],
    )
    write_csv(
        outdir / "hbond_timeseries_pbc.csv",
        pbc.hbond_timeseries,
        ["frame", "step", "time_ns", "n_solute_to_water", "n_water_to_solute", "n_total"],
    )

    def _hb_stats(arr: np.ndarray) -> tuple[float, float, float]:
        if arr.size == 0:
            return float("nan"), float("nan"), float("nan")
        return float(np.mean(arr)), float(np.median(arr)), float(np.std(arr, ddof=1) if arr.size > 1 else 0.0)

    hb_rows = []
    for label in ["solute_to_water", "water_to_solute", "all"]:
        a = npbc.hbond_lifetimes_ps[label]
        b = pbc.hbond_lifetimes_ps[label]
        am, amd, asd = _hb_stats(a)
        bm, bmd, bsd = _hb_stats(b)
        hb_rows.append((label, am, amd, asd, bm, bmd, bsd, am - bm))
    write_csv(
        outdir / "hbond_lifetime_summary.csv",
        hb_rows,
        [
            "class",
            "npbc_mean_ps",
            "npbc_median_ps",
            "npbc_std_ps",
            "pbc_mean_ps",
            "pbc_median_ps",
            "pbc_std_ps",
            "delta_mean_npbc_minus_pbc_ps",
        ],
    )

    # Radial density / coordination
    rdf_rows = []
    coord_rows = []
    for cat in ["all_heavy", "O_acceptor", "N_site"]:
        rn = npbc.rdf[cat]["rho"]
        rp = pbc.rdf[cat]["rho"]
        coord_rows.append((cat, float(npbc.rdf[cat]["cn"]), float(pbc.rdf[cat]["cn"]), float(npbc.rdf[cat]["cn"] - pbc.rdf[cat]["cn"])))
        for i, r in enumerate(npbc.rdf_centers):
            delta = float(rn[i] - rp[i])
            rel = float(delta / rp[i]) if np.isfinite(rp[i]) and abs(rp[i]) > 1.0e-12 else float("nan")
            rdf_rows.append((cat, float(r), float(rn[i]), float(rp[i]), delta, rel))

    write_csv(
        outdir / "site_water_radial_density_compare.csv",
        rdf_rows,
        ["site_category", "r_A", "rho_npbc_molA3", "rho_pbc_molA3", "delta_npbc_minus_pbc", "rel_delta_vs_pbc"],
    )
    write_csv(
        outdir / "coordination_summary.csv",
        coord_rows,
        ["site_category", "CN_npbc", "CN_pbc", "delta_npbc_minus_pbc"],
    )

    # Plots
    if HAVE_PLOT:
        # 1) torsion ACF
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
        lags_ps = np.asarray([r[1] for r in acf_rows], dtype=float)
        axs[0].plot(lags_ps, [r[2] for r in acf_rows], label="NPBC")
        axs[0].plot(lags_ps, [r[3] for r in acf_rows], label="PBC")
        axs[0].set_title("phi autocorrelation")
        axs[0].set_xlabel("lag (ps)")
        axs[0].set_ylabel("Cphi(lag)")
        axs[0].grid(alpha=0.25)
        axs[0].legend(loc="best", fontsize=8)

        axs[1].plot(lags_ps, [r[4] for r in acf_rows], label="NPBC")
        axs[1].plot(lags_ps, [r[5] for r in acf_rows], label="PBC")
        axs[1].set_title("psi autocorrelation")
        axs[1].set_xlabel("lag (ps)")
        axs[1].set_ylabel("Cpsi(lag)")
        axs[1].grid(alpha=0.25)
        axs[1].legend(loc="best", fontsize=8)
        fig.savefig(outdir / "17_torsion_autocorrelation.png", dpi=170)
        plt.close(fig)

        # 2) hydration shell occupancy
        n_shell_npbc = np.asarray([[r[2], r[3], r[4]] for r in npbc.shell_timeseries], dtype=float)
        n_shell_pbc = np.asarray([[r[2], r[3], r[4]] for r in pbc.shell_timeseries], dtype=float)
        fig, axs = plt.subplots(1, 2, figsize=(12, 4.2), constrained_layout=True)
        axs[0].plot(n_shell_npbc[:, 0], n_shell_npbc[:, 1], alpha=0.3, lw=0.7, label="NPBC s1 raw")
        axs[0].plot(n_shell_pbc[:, 0], n_shell_pbc[:, 1], alpha=0.3, lw=0.7, label="PBC s1 raw")
        axs[0].plot(n_shell_npbc[:, 0], _rolling(n_shell_npbc[:, 1], 100), lw=1.8, label="NPBC s1 roll")
        axs[0].plot(n_shell_pbc[:, 0], _rolling(n_shell_pbc[:, 1], 100), lw=1.8, label="PBC s1 roll")
        axs[0].set_title("Hydration shell-1 occupancy")
        axs[0].set_xlabel("time (ns)")
        axs[0].set_ylabel("N waters")
        axs[0].grid(alpha=0.25)
        axs[0].legend(loc="best", fontsize=8)

        axs[1].plot(n_shell_npbc[:, 0], n_shell_npbc[:, 2], alpha=0.3, lw=0.7, label="NPBC s2 raw")
        axs[1].plot(n_shell_pbc[:, 0], n_shell_pbc[:, 2], alpha=0.3, lw=0.7, label="PBC s2 raw")
        axs[1].plot(n_shell_npbc[:, 0], _rolling(n_shell_npbc[:, 2], 100), lw=1.8, label="NPBC s2 roll")
        axs[1].plot(n_shell_pbc[:, 0], _rolling(n_shell_pbc[:, 2], 100), lw=1.8, label="PBC s2 roll")
        axs[1].set_title("Hydration shell-2 occupancy")
        axs[1].set_xlabel("time (ns)")
        axs[1].set_ylabel("N waters")
        axs[1].grid(alpha=0.25)
        axs[1].legend(loc="best", fontsize=8)
        fig.savefig(outdir / "18_hydration_shell_occupancy.png", dpi=170)
        plt.close(fig)

        # 3) HBond counts
        hb_npbc = np.asarray([[r[2], r[3], r[4], r[5]] for r in npbc.hbond_timeseries], dtype=float)
        hb_pbc = np.asarray([[r[2], r[3], r[4], r[5]] for r in pbc.hbond_timeseries], dtype=float)
        fig, axs = plt.subplots(1, 3, figsize=(15, 4.2), constrained_layout=True)
        for k, title in [(1, "solute->water"), (2, "water->solute"), (3, "total")]:
            axs[k - 1].plot(hb_npbc[:, 0], hb_npbc[:, k], alpha=0.3, lw=0.7, label="NPBC raw")
            axs[k - 1].plot(hb_pbc[:, 0], hb_pbc[:, k], alpha=0.3, lw=0.7, label="PBC raw")
            axs[k - 1].plot(hb_npbc[:, 0], _rolling(hb_npbc[:, k], 100), lw=1.8, label="NPBC roll")
            axs[k - 1].plot(hb_pbc[:, 0], _rolling(hb_pbc[:, k], 100), lw=1.8, label="PBC roll")
            axs[k - 1].set_title(f"HBonds: {title}")
            axs[k - 1].set_xlabel("time (ns)")
            axs[k - 1].set_ylabel("count")
            axs[k - 1].grid(alpha=0.25)
            axs[k - 1].legend(loc="best", fontsize=8)
        fig.savefig(outdir / "19_hbond_counts_timeseries.png", dpi=170)
        plt.close(fig)

        # 4) site-water radial densities
        fig, axs = plt.subplots(1, 3, figsize=(16, 4.2), constrained_layout=True)
        for i, cat in enumerate(["all_heavy", "O_acceptor", "N_site"]):
            axs[i].plot(npbc.rdf_centers, npbc.rdf[cat]["rho"], label="NPBC")
            axs[i].plot(pbc.rdf_centers, pbc.rdf[cat]["rho"], label="PBC")
            axs[i].set_title(f"site-water density: {cat}")
            axs[i].set_xlabel("r (A)")
            axs[i].set_ylabel("rho_O(r) [mol/A^3]")
            axs[i].grid(alpha=0.25)
            axs[i].legend(loc="best", fontsize=8)
        fig.savefig(outdir / "20_site_water_radial_density.png", dpi=170)
        plt.close(fig)

    metrics = {
        "frames": {
            "npbc": int(npbc.phi.size),
            "pbc": int(pbc.phi.size),
        },
        "torsion_autocorrelation": {
            "iact_phi_frames_npbc": float(npbc.iact_phi_frames),
            "iact_phi_frames_pbc": float(pbc.iact_phi_frames),
            "iact_psi_frames_npbc": float(npbc.iact_psi_frames),
            "iact_psi_frames_pbc": float(pbc.iact_psi_frames),
            "frame_dt_ps_npbc": float(npbc.frame_dt_ps),
            "frame_dt_ps_pbc": float(pbc.frame_dt_ps),
        },
        "coordination": {
            cat: {
                "npbc": float(npbc.rdf[cat]["cn"]),
                "pbc": float(pbc.rdf[cat]["cn"]),
                "delta_npbc_minus_pbc": float(npbc.rdf[cat]["cn"] - pbc.rdf[cat]["cn"]),
            }
            for cat in ["all_heavy", "O_acceptor", "N_site"]
        },
        "hbond_lifetime_ps": {
            "solute_to_water_mean_npbc": float(np.mean(npbc.hbond_lifetimes_ps["solute_to_water"])) if npbc.hbond_lifetimes_ps["solute_to_water"].size else float("nan"),
            "solute_to_water_mean_pbc": float(np.mean(pbc.hbond_lifetimes_ps["solute_to_water"])) if pbc.hbond_lifetimes_ps["solute_to_water"].size else float("nan"),
            "water_to_solute_mean_npbc": float(np.mean(npbc.hbond_lifetimes_ps["water_to_solute"])) if npbc.hbond_lifetimes_ps["water_to_solute"].size else float("nan"),
            "water_to_solute_mean_pbc": float(np.mean(pbc.hbond_lifetimes_ps["water_to_solute"])) if pbc.hbond_lifetimes_ps["water_to_solute"].size else float("nan"),
            "all_mean_npbc": float(np.mean(npbc.hbond_lifetimes_ps["all"])) if npbc.hbond_lifetimes_ps["all"].size else float("nan"),
            "all_mean_pbc": float(np.mean(pbc.hbond_lifetimes_ps["all"])) if pbc.hbond_lifetimes_ps["all"].size else float("nan"),
        },
    }
    (outdir / "high_value_metrics_summary.json").write_text(json.dumps(metrics, indent=2) + "\n")

    report = []
    report.append("# Stage13corr high-value observables (NPBC vs PBC)")
    report.append("")
    report.append(f"- Frames NPBC/PBC: {npbc.phi.size} / {pbc.phi.size}")
    report.append(f"- IACT phi (frames) NPBC/PBC: {npbc.iact_phi_frames:.2f} / {pbc.iact_phi_frames:.2f}")
    report.append(f"- IACT psi (frames) NPBC/PBC: {npbc.iact_psi_frames:.2f} / {pbc.iact_psi_frames:.2f}")
    report.append("")
    report.append("## Coordination numbers (waters around solute sites)")
    for cat in ["all_heavy", "O_acceptor", "N_site"]:
        report.append(
            f"- {cat}: NPBC={npbc.rdf[cat]['cn']:.3f}, PBC={pbc.rdf[cat]['cn']:.3f}, delta={npbc.rdf[cat]['cn'] - pbc.rdf[cat]['cn']:+.3f}"
        )
    report.append("")
    report.append("## H-bond lifetime means (ps)")
    for lbl in ["solute_to_water", "water_to_solute", "all"]:
        a = npbc.hbond_lifetimes_ps[lbl]
        b = pbc.hbond_lifetimes_ps[lbl]
        am = float(np.mean(a)) if a.size else float("nan")
        bm = float(np.mean(b)) if b.size else float("nan")
        report.append(f"- {lbl}: NPBC={am:.3f}, PBC={bm:.3f}, delta={am - bm:+.3f}")
    (outdir / "high_value_report.md").write_text("\n".join(report) + "\n")

    print(outdir)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
