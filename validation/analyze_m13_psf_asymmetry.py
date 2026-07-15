"""Diagnose asymmetric PSF-minus-aperture scatter in the M13 validation run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.ndimage import shift as image_shift
from scipy.spatial import cKDTree
from scipy.stats import skew, spearmanr


def _numeric(frame: pd.DataFrame, columns: tuple[str, ...]) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def _quantile_stats(values: pd.Series | np.ndarray) -> dict[str, float | int]:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return {
            "n": 0,
            "median": np.nan,
            "q05": np.nan,
            "q16": np.nan,
            "q84": np.nan,
            "q95": np.nan,
            "lower_width": np.nan,
            "upper_width": np.nan,
            "quantile_asymmetry": np.nan,
            "skewness": np.nan,
        }
    q05, q16, median, q84, q95 = np.percentile(data, (5, 16, 50, 84, 95))
    lower = float(median - q16)
    upper = float(q84 - median)
    denominator = lower + upper
    return {
        "n": int(data.size),
        "median": float(median),
        "q05": float(q05),
        "q16": float(q16),
        "q84": float(q84),
        "q95": float(q95),
        "lower_width": lower,
        "upper_width": upper,
        "quantile_asymmetry": float((upper - lower) / denominator) if denominator > 0 else 0.0,
        "skewness": float(skew(data, bias=False)) if data.size >= 3 else np.nan,
    }


def _match_catalogues(result_dir: Path, filename: str) -> pd.DataFrame:
    psf = pd.read_csv(result_dir / "cmd_psf" / f"photometry_{filename}.tsv", sep="\t")
    aperture = pd.read_csv(
        result_dir / "step7_forced_phot" / f"photometry_{filename}.tsv",
        sep="\t",
        low_memory=False,
    )
    _numeric(
        psf,
        (
            "x_fit",
            "y_fit",
            "flux_psf_e",
            "snr_psf",
            "qfit",
            "cfit",
            "reduced_chi2",
            "iter_found",
            "flags_psf",
        ),
    )
    _numeric(
        aperture,
        ("x_fit", "y_fit", "flux_e", "snr", "mag_inst", "sky", "sky_std"),
    )

    psf_ok = (
        np.isfinite(psf["x_fit"])
        & np.isfinite(psf["y_fit"])
        & np.isfinite(psf["flux_psf_e"])
        & (psf["flux_psf_e"] > 0)
        & (psf["flags_psf"] == 0)
    )
    aperture_ok = (
        np.isfinite(aperture["x_fit"])
        & np.isfinite(aperture["y_fit"])
        & np.isfinite(aperture["flux_e"])
        & (aperture["flux_e"] > 0)
        & np.isfinite(aperture["snr"])
    )
    if "bad_phot_flag" in aperture:
        aperture_ok &= ~_as_bool(aperture["bad_phot_flag"])
    if "off_frame_flag" in aperture:
        aperture_ok &= ~_as_bool(aperture["off_frame_flag"])

    psf = psf.loc[psf_ok].copy().reset_index(drop=True)
    aperture = aperture.loc[aperture_ok].copy().reset_index(drop=True)
    aperture_xy = aperture[["x_fit", "y_fit"]].to_numpy()
    tree = cKDTree(aperture_xy)
    distance, aperture_index = tree.query(psf[["x_fit", "y_fit"]].to_numpy(), k=1)
    pairs = pd.DataFrame(
        {
            "psf_index": np.arange(len(psf)),
            "aperture_index": aperture_index,
            "match_distance_px": distance,
        }
    )
    pairs = pairs.loc[pairs["match_distance_px"] <= 1.5]
    pairs = pairs.sort_values("match_distance_px").drop_duplicates("aperture_index")
    p = psf.iloc[pairs["psf_index"].to_numpy()].reset_index(drop=True)
    a = aperture.iloc[pairs["aperture_index"].to_numpy()].reset_index(drop=True)

    if len(aperture_xy) >= 2:
        nearest = tree.query(aperture_xy, k=2)[0][:, 1]
    else:
        nearest = np.full(len(aperture), np.nan)
    selected_nearest = nearest[pairs["aperture_index"].to_numpy()]
    output = pd.DataFrame(
        {
            "x": p["x_fit"].to_numpy(),
            "y": p["y_fit"].to_numpy(),
            "snr_aperture": a["snr"].to_numpy(),
            "snr_psf": p["snr_psf"].to_numpy(),
            "flux_psf_e": p["flux_psf_e"].to_numpy(),
            "flux_aperture_e": a["flux_e"].to_numpy(),
            "qfit": p["qfit"].to_numpy(),
            "reduced_chi2": p["reduced_chi2"].to_numpy(),
            "iter_found": p["iter_found"].to_numpy(),
            "sky": a["sky"].to_numpy(),
            "sky_std": a["sky_std"].to_numpy(),
            "nearest_neighbor_px": selected_nearest,
            "match_distance_px": pairs["match_distance_px"].to_numpy(),
        }
    )
    if "detected_flag" in a:
        output["detected"] = _as_bool(a["detected_flag"]).to_numpy()
    else:
        output["detected"] = True
    output["delta_mag"] = -2.5 * np.log10(
        output["flux_psf_e"] / output["flux_aperture_e"]
    )
    return output.replace([np.inf, -np.inf], np.nan).dropna(subset=["delta_mag"])


def _epsf_metrics(path: Path) -> dict[str, float | str]:
    model = np.asarray(fits.getdata(path), dtype=float)
    weights = np.clip(model, 0.0, None)
    yy, xx = np.mgrid[: model.shape[0], : model.shape[1]]
    total = float(weights.sum())
    cx = float(np.sum(weights * xx) / total)
    cy = float(np.sum(weights * yy) / total)
    dx = xx - cx
    dy = yy - cy
    covariance = np.array(
        [
            [np.sum(weights * dx * dx), np.sum(weights * dx * dy)],
            [np.sum(weights * dx * dy), np.sum(weights * dy * dy)],
        ],
        dtype=float,
    ) / total
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    order = np.argsort(eigenvalues)[::-1]
    major, minor = eigenvalues[order]
    vector = eigenvectors[:, order[0]]
    ellipticity = 1.0 - np.sqrt(max(minor, 0.0) / max(major, 1e-12))
    angle = float(np.degrees(np.arctan2(vector[1], vector[0])))

    geometric_x = 0.5 * (model.shape[1] - 1)
    geometric_y = 0.5 * (model.shape[0] - 1)
    centered = image_shift(
        model,
        shift=(geometric_y - cy, geometric_x - cx),
        order=3,
        mode="nearest",
        prefilter=True,
    )
    rotated = np.rot90(centered, 2)
    denominator = 2.0 * np.sum(np.abs(centered))
    asymmetry = float(np.sum(np.abs(centered - rotated)) / denominator)
    return {
        "file": path.name,
        "centroid_x_oversampled_px": cx,
        "centroid_y_oversampled_px": cy,
        "centroid_offset_native_px": float(np.hypot(cx - geometric_x, cy - geometric_y) / 2.0),
        "ellipticity": float(ellipticity),
        "position_angle_deg": angle,
        "rotation_asymmetry": asymmetry,
        "positive_sum": total,
    }


def _stack_residuals(
    residual: np.ndarray,
    stars: pd.DataFrame,
    *,
    half_size: int = 10,
) -> tuple[np.ndarray, int, float]:
    stamps = []
    height, width = residual.shape
    for row in stars.itertuples(index=False):
        x_center = int(round(float(row.x)))
        y_center = int(round(float(row.y)))
        if (
            x_center - half_size < 0
            or y_center - half_size < 0
            or x_center + half_size >= width
            or y_center + half_size >= height
        ):
            continue
        stamp = residual[
            y_center - half_size : y_center + half_size + 1,
            x_center - half_size : x_center + half_size + 1,
        ].astype(float, copy=True)
        stamp = image_shift(
            stamp,
            shift=(y_center - float(row.y), x_center - float(row.x)),
            order=3,
            mode="nearest",
            prefilter=True,
        )
        flux_adu = float(row.flux_psf_e) / 0.689
        if np.isfinite(flux_adu) and flux_adu > 0:
            stamps.append(stamp / flux_adu)
    if not stamps:
        empty = np.full((2 * half_size + 1, 2 * half_size + 1), np.nan)
        return empty, 0, np.nan
    stack = np.nanmedian(np.stack(stamps), axis=0)
    denominator = 2.0 * np.nansum(np.abs(stack))
    asymmetry = float(np.nansum(np.abs(stack - np.rot90(stack, 2))) / denominator)
    return stack, len(stamps), asymmetry


def _group_rows(data: pd.DataFrame, group_type: str, group_column: str) -> list[dict]:
    rows = []
    for name, group in data.groupby(group_column, observed=True, sort=True):
        row = {"group_type": group_type, "group": str(name)}
        row.update(_quantile_stats(group["delta_centered_mag"]))
        rows.append(row)
    return rows


def _plot_group_quantiles(ax, data: pd.DataFrame, column: str, title: str) -> None:
    groups = list(data.groupby(column, observed=True, sort=True))
    positions = np.arange(len(groups))
    medians = []
    lower = []
    upper = []
    labels = []
    for name, group in groups:
        stats = _quantile_stats(group["delta_centered_mag"])
        medians.append(stats["median"])
        lower.append(stats["lower_width"])
        upper.append(stats["upper_width"])
        labels.append(f"{name}\nN={stats['n']}")
    ax.errorbar(
        positions,
        medians,
        yerr=np.vstack([lower, upper]),
        fmt="o",
        color="#0072B2",
        ecolor="#56B4E9",
        capsize=3,
    )
    ax.axhline(0.0, color="0.45", linewidth=0.8, linestyle="--")
    ax.set_xticks(positions, labels, fontsize=8)
    ax.set_ylabel("centered PSF - aperture (mag)")
    ax.set_title(title, loc="left")
    ax.grid(axis="y", alpha=0.18)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project",
        type=Path,
        default=Path(r"E:\APEX_validation\reprocess\M13"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation/real_gui_run"),
    )
    args = parser.parse_args()

    result_dir = args.project / "result"
    psf_dir = result_dir / "cmd_psf"
    meta_path = next(psf_dir.glob("residual_meta_*.json"))
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    filename = meta["file"]
    core = meta["core_cut"]
    index = pd.read_csv(psf_dir / "photometry_index.csv").iloc[0]
    fwhm_px = 6.942195606701463
    pixel_scale_arcsec = 2.7281805081475565 / fwhm_px

    data = _match_catalogues(result_dir, filename)
    data["radius_px"] = np.hypot(
        data["x"] - float(core["center_x"]),
        data["y"] - float(core["center_y"]),
    )
    data["radius_arcmin"] = data["radius_px"] * pixel_scale_arcsec / 60.0
    data["nearest_neighbor_fwhm"] = data["nearest_neighbor_px"] / fwhm_px
    reference = data.loc[data["snr_aperture"] >= 100, "delta_mag"]
    if len(reference) < 20:
        reference = data.loc[data["snr_aperture"] >= 50, "delta_mag"]
    reference_median = float(np.median(reference))
    data["delta_centered_mag"] = data["delta_mag"] - reference_median

    data["snr_bin"] = pd.cut(
        data["snr_aperture"],
        bins=[0, 5, 10, 20, 50, 100, np.inf],
        labels=["<5", "5-10", "10-20", "20-50", "50-100", ">=100"],
    )
    core_arcmin = float(core["radius_px"]) * pixel_scale_arcsec / 60.0
    data["radius_bin"] = pd.cut(
        data["radius_arcmin"],
        bins=[0, core_arcmin, 1.5, 3.0, 5.0, np.inf],
        labels=["core", "0.9-1.5'", "1.5-3'", "3-5'", ">=5'"],
    )
    data["neighbor_bin"] = pd.cut(
        data["nearest_neighbor_fwhm"],
        bins=[0, 1, 2, 4, np.inf],
        labels=["<1 FWHM", "1-2 FWHM", "2-4 FWHM", ">=4 FWHM"],
    )
    data["qfit_bin"] = pd.cut(
        data["qfit"],
        bins=[-np.inf, 0.1, 0.2, 0.4, np.inf],
        labels=["<=0.1", "0.1-0.2", "0.2-0.4", ">0.4"],
    )
    data["detection_class"] = np.where(data["detected"], "Step4 detected", "forced only")
    data["quadrant"] = np.select(
        [
            (data["x"] >= core["center_x"]) & (data["y"] >= core["center_y"]),
            (data["x"] < core["center_x"]) & (data["y"] >= core["center_y"]),
            (data["x"] < core["center_x"]) & (data["y"] < core["center_y"]),
        ],
        ["NE", "NW", "SW"],
        default="SE",
    )

    group_rows = []
    for group_type, column in (
        ("snr", "snr_bin"),
        ("radius", "radius_bin"),
        ("neighbor", "neighbor_bin"),
        ("qfit", "qfit_bin"),
        ("detection", "detection_class"),
        ("quadrant", "quadrant"),
    ):
        group_rows.extend(_group_rows(data, group_type, column))
    group_frame = pd.DataFrame(group_rows)

    correlations = {}
    for name, values in (
        ("log10_snr", np.log10(data["snr_aperture"])),
        ("radius_arcmin", data["radius_arcmin"]),
        ("neighbor_fwhm", data["nearest_neighbor_fwhm"]),
        ("qfit", data["qfit"]),
        ("reduced_chi2", data["reduced_chi2"]),
        ("sky", data["sky"]),
        ("sky_std", data["sky_std"]),
        ("x", data["x"]),
        ("y", data["y"]),
    ):
        finite = np.isfinite(values) & np.isfinite(data["delta_centered_mag"])
        coefficient, p_value = spearmanr(
            np.asarray(values)[finite],
            data.loc[finite, "delta_centered_mag"],
        )
        correlations[name] = {"rho": float(coefficient), "p_value": float(p_value)}

    epsf_rows = []
    validation_root = args.project.parent
    for label in ("NGC6811", "M67", "M13"):
        model_path = next((validation_root / label / "result" / "cmd_psf").glob("epsf_model_*.fits"))
        row = {"project": label}
        row.update(_epsf_metrics(model_path))
        epsf_rows.append(row)
    epsf_frame = pd.DataFrame(epsf_rows)
    m13_epsf = epsf_frame.loc[epsf_frame["project"] == "M13"].iloc[0]

    final_residual = np.asarray(
        fits.getdata(psf_dir / meta["iters"][-1]["residual_path"]),
        dtype=float,
    )
    stack_sample = data.loc[
        (data["snr_aperture"] >= 50)
        & (data["qfit"] <= 0.3)
        & (data["nearest_neighbor_fwhm"] >= 2.0)
    ].copy()
    stacks = []
    all_stack, all_n, all_asymmetry = _stack_residuals(final_residual, stack_sample)
    stacks.append(("all", all_stack, all_n, all_asymmetry))
    for quadrant in ("NE", "NW", "SW", "SE"):
        stack, count, asymmetry = _stack_residuals(
            final_residual,
            stack_sample.loc[stack_sample["quadrant"] == quadrant],
        )
        stacks.append((quadrant, stack, count, asymmetry))

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 9.2), constrained_layout=True)
    ax = axes[0, 0]
    ax.hexbin(
        data["snr_aperture"],
        data["delta_centered_mag"],
        xscale="log",
        gridsize=42,
        mincnt=1,
        cmap="Blues",
        bins="log",
    )
    ax.axhline(0.0, color="0.35", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Step7 aperture SNR")
    ax.set_ylabel("centered PSF - aperture (mag)")
    ax.set_title("(a) asymmetric low-SNR tail", loc="left")

    ax = axes[0, 1]
    for label, selection, color in (
        ("SNR 5-10", (data["snr_aperture"] >= 5) & (data["snr_aperture"] < 10), "#D55E00"),
        ("SNR >= 50", data["snr_aperture"] >= 50, "#0072B2"),
    ):
        values = data.loc[selection, "delta_centered_mag"]
        ax.hist(values, bins=40, histtype="step", density=True, linewidth=1.6, label=label, color=color)
    ax.axvline(0.0, color="0.35", linewidth=0.8, linestyle="--")
    ax.set_xlabel("centered PSF - aperture (mag)")
    ax.set_ylabel("density")
    ax.set_title("(b) tail collapses at high SNR", loc="left")
    ax.legend()

    _plot_group_quantiles(axes[0, 2], data, "radius_bin", "(c) radial dependence")
    _plot_group_quantiles(axes[1, 0], data, "neighbor_bin", "(d) crowding dependence")

    ax = axes[1, 1]
    color_limit = float(np.percentile(np.abs(data["delta_centered_mag"]), 90))
    scatter_plot = ax.scatter(
        data["x"],
        data["y"],
        c=data["delta_centered_mag"],
        s=8,
        alpha=0.65,
        cmap="coolwarm",
        vmin=-color_limit,
        vmax=color_limit,
        linewidths=0,
    )
    ax.add_patch(
        Circle(
            (float(core["center_x"]), float(core["center_y"])),
            float(core["radius_px"]),
            fill=False,
            color="black",
            linewidth=1.1,
        )
    )
    ax.set_aspect("equal")
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    ax.set_title("(e) detector-position pattern", loc="left")
    fig.colorbar(scatter_plot, ax=ax, label="centered delta mag", shrink=0.78)

    ax = axes[1, 2]
    epsf_path = psf_dir / str(m13_epsf["file"])
    epsf = np.asarray(fits.getdata(epsf_path), dtype=float)
    ax.imshow(epsf, origin="lower", cmap="magma", interpolation="nearest")
    ax.set_title(
        "(f) M13 ePSF\n"
        f"ellipticity={float(m13_epsf['ellipticity']):.3f}, "
        f"A180={float(m13_epsf['rotation_asymmetry']):.3f}",
        loc="left",
    )
    ax.set_xlabel("oversampled x")
    ax.set_ylabel("oversampled y")
    fig.suptitle(
        f"M13 PSF asymmetry diagnostics: N={len(data)}, high-SNR offset={reference_median:+.3f} mag",
        fontsize=15,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    diagnostic_path = args.output_dir / "m13_psf_asymmetry_diagnostics.png"
    fig.savefig(diagnostic_path, dpi=200, facecolor="white")
    plt.close(fig)

    residual_values = np.concatenate(
        [stack[np.isfinite(stack)] for _, stack, _, _ in stacks if np.any(np.isfinite(stack))]
    )
    residual_limit = float(np.percentile(np.abs(residual_values), 99.0))
    fig, axes = plt.subplots(1, 5, figsize=(14.0, 3.1), constrained_layout=True)
    residual_stack_rows = []
    for ax, (label, stack, count, asymmetry) in zip(axes, stacks):
        ax.imshow(
            stack * 1e4,
            origin="lower",
            cmap="coolwarm",
            vmin=-residual_limit * 1e4,
            vmax=residual_limit * 1e4,
            interpolation="nearest",
        )
        ax.set_title(f"{label}: N={count}, A180={asymmetry:.3f}", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        residual_stack_rows.append(
            {
                "region": label,
                "n_stars": count,
                "rotation_asymmetry": asymmetry,
                "peak_abs_flux_fraction": float(np.nanmax(np.abs(stack))),
                "sum_abs_flux_fraction": float(np.nansum(np.abs(stack))),
                "sum_signed_flux_fraction": float(np.nansum(stack)),
            }
        )
    fig.suptitle("Median normalized residual stacks (SNR >= 50, qfit <= 0.3, separation >= 2 FWHM)")
    residual_path = args.output_dir / "m13_psf_residual_stacks.png"
    fig.savefig(residual_path, dpi=210, facecolor="white")
    plt.close(fig)

    group_path = args.output_dir / "m13_psf_asymmetry_groups.csv"
    epsf_path_out = args.output_dir / "m13_epsf_shape_metrics.csv"
    matched_path = args.output_dir / "m13_psf_asymmetry_matched.csv"
    summary_path = args.output_dir / "m13_psf_asymmetry_summary.json"
    group_frame.to_csv(group_path, index=False)
    epsf_frame.to_csv(epsf_path_out, index=False)
    data.to_csv(matched_path, index=False)
    summary = {
        "file": filename,
        "n_matched": int(len(data)),
        "reference_snr_min": 100 if int(np.sum(data["snr_aperture"] >= 100)) >= 20 else 50,
        "reference_offset_mag": reference_median,
        "overall_centered": _quantile_stats(data["delta_centered_mag"]),
        "core_radius_px": float(core["radius_px"]),
        "core_radius_arcmin": core_arcmin,
        "n_psf": int(index["n"]),
        "n_goodmag": int(index["n_goodmag"]),
        "correlations": correlations,
        "epsf_metrics": epsf_rows,
        "residual_stacks": residual_stack_rows,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(diagnostic_path.resolve())
    print(residual_path.resolve())
    print(group_path.resolve())
    print(summary_path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
