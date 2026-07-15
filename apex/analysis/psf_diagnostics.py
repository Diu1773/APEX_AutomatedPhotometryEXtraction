"""Final PSF diagnostics shared by the Step 8 GUI and validation tools."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Circle
from scipy.ndimage import shift as image_shift
from scipy.spatial import cKDTree


_TRUE_VALUES = {"true", "1", "yes", "y", "t"}


def _numeric(frame: pd.DataFrame, columns: tuple[str, ...]) -> None:
    for column in columns:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(_TRUE_VALUES)


def _column(frame: pd.DataFrame, name: str, default: float = np.nan) -> np.ndarray:
    if name not in frame.columns:
        return np.full(len(frame), default, dtype=float)
    return pd.to_numeric(frame[name], errors="coerce").to_numpy(dtype=float)


def load_psf_final_diagnostic_data(
    result_dir: str | Path,
    filename: str,
    *,
    match_radius_px: float = 1.5,
) -> pd.DataFrame:
    """Position-match clean PSF measurements to Step 7 aperture measurements."""
    result_dir = Path(result_dir)
    psf_path = result_dir / "cmd_psf" / f"photometry_{filename}.tsv"
    aperture_path = result_dir / "step7_forced_phot" / f"photometry_{filename}.tsv"
    if not psf_path.exists():
        raise FileNotFoundError(f"PSF photometry not found: {psf_path}")
    if not aperture_path.exists():
        raise FileNotFoundError(f"Aperture photometry not found: {aperture_path}")

    psf = pd.read_csv(psf_path, sep="\t", low_memory=False)
    aperture = pd.read_csv(aperture_path, sep="\t", low_memory=False)
    _numeric(
        psf,
        (
            "x_fit",
            "y_fit",
            "flux_psf_e",
            "snr_psf",
            "qfit",
            "qfit_noise_ratio",
            "reduced_chi2",
            "iter_found",
            "flags_psf",
            "fit_window_px",
            "fit_window_energy",
            "psf_nea_px",
        ),
    )
    _numeric(
        aperture,
        ("x_fit", "y_fit", "flux_e", "snr", "sky", "sky_std"),
    )

    required_psf = {"x_fit", "y_fit", "flux_psf_e"}
    required_aperture = {"x_fit", "y_fit", "flux_e", "snr"}
    if missing := required_psf.difference(psf.columns):
        raise ValueError(f"PSF table is missing columns: {sorted(missing)}")
    if missing := required_aperture.difference(aperture.columns):
        raise ValueError(f"Aperture table is missing columns: {sorted(missing)}")

    aperture_geometry_ok = np.isfinite(aperture["x_fit"]) & np.isfinite(aperture["y_fit"])
    if "off_frame_flag" in aperture.columns:
        aperture_geometry_ok &= ~_as_bool(aperture["off_frame_flag"])
    aperture_geometry = aperture.loc[aperture_geometry_ok, ["x_fit", "y_fit"]]
    geometry_xy = aperture_geometry.to_numpy(dtype=float)
    nearest_by_index: dict[object, float] = {}
    if len(geometry_xy) >= 2:
        geometry_tree = cKDTree(geometry_xy)
        nearest = geometry_tree.query(geometry_xy, k=2)[0][:, 1]
        nearest_by_index = dict(zip(aperture_geometry.index, nearest))
    elif len(geometry_xy) == 1:
        nearest_by_index = {aperture_geometry.index[0]: np.nan}

    psf_ok = (
        np.isfinite(psf["x_fit"])
        & np.isfinite(psf["y_fit"])
        & np.isfinite(psf["flux_psf_e"])
        & (psf["flux_psf_e"] > 0)
    )
    if "flags_psf" in psf.columns:
        psf_ok &= psf["flags_psf"].fillna(1).eq(0)

    aperture_ok = (
        aperture_geometry_ok
        & np.isfinite(aperture["flux_e"])
        & (aperture["flux_e"] > 0)
        & np.isfinite(aperture["snr"])
        & (aperture["snr"] > 0)
    )
    if "bad_phot_flag" in aperture.columns:
        aperture_ok &= ~_as_bool(aperture["bad_phot_flag"])

    psf_clean = psf.loc[psf_ok].copy().reset_index(drop=True)
    aperture_clean = aperture.loc[aperture_ok].copy()
    if psf_clean.empty or aperture_clean.empty:
        return pd.DataFrame()

    aperture_xy = aperture_clean[["x_fit", "y_fit"]].to_numpy(dtype=float)
    tree = cKDTree(aperture_xy)
    distance, aperture_position = tree.query(
        psf_clean[["x_fit", "y_fit"]].to_numpy(dtype=float),
        k=1,
    )
    pairs = pd.DataFrame(
        {
            "psf_index": np.arange(len(psf_clean), dtype=int),
            "aperture_position": aperture_position.astype(int),
            "match_distance_px": distance,
        }
    )
    pairs = pairs.loc[pairs["match_distance_px"] <= float(match_radius_px)]
    pairs = pairs.sort_values("match_distance_px").drop_duplicates("aperture_position")
    if pairs.empty:
        return pd.DataFrame()

    p = psf_clean.iloc[pairs["psf_index"].to_numpy()].reset_index(drop=True)
    selected_aperture_index = aperture_clean.index[pairs["aperture_position"].to_numpy()]
    a = aperture_clean.loc[selected_aperture_index].reset_index(drop=True)
    nearest = np.asarray(
        [nearest_by_index.get(index, np.nan) for index in selected_aperture_index],
        dtype=float,
    )

    output = pd.DataFrame(
        {
            "x": p["x_fit"].to_numpy(dtype=float),
            "y": p["y_fit"].to_numpy(dtype=float),
            "snr_aperture": a["snr"].to_numpy(dtype=float),
            "snr_psf": _column(p, "snr_psf"),
            "flux_psf_e": p["flux_psf_e"].to_numpy(dtype=float),
            "flux_aperture_e": a["flux_e"].to_numpy(dtype=float),
            "qfit": _column(p, "qfit"),
            "qfit_noise_ratio": _column(p, "qfit_noise_ratio"),
            "reduced_chi2": _column(p, "reduced_chi2"),
            "iter_found": _column(p, "iter_found"),
            "fit_window_px": _column(p, "fit_window_px"),
            "fit_window_energy": _column(p, "fit_window_energy"),
            "psf_nea_px": _column(p, "psf_nea_px"),
            "sky": _column(a, "sky"),
            "sky_std": _column(a, "sky_std"),
            "nearest_neighbor_px": nearest,
            "match_distance_px": pairs["match_distance_px"].to_numpy(dtype=float),
        }
    )
    if "detected_flag" in a.columns:
        output["detected"] = _as_bool(a["detected_flag"]).to_numpy(dtype=bool)
    else:
        output["detected"] = True
    output["delta_mag"] = -2.5 * np.log10(
        output["flux_psf_e"] / output["flux_aperture_e"]
    )
    return output.replace([np.inf, -np.inf], np.nan).dropna(subset=["delta_mag"])


def epsf_shape_metrics(model: np.ndarray, *, oversampling: float = 2.0) -> dict[str, float]:
    """Measure second-moment ellipticity and 180-degree rotational asymmetry."""
    data = np.asarray(model, dtype=float)
    weights = np.clip(np.nan_to_num(data, nan=0.0), 0.0, None)
    total = float(weights.sum())
    empty = {
        "centroid_x_oversampled_px": np.nan,
        "centroid_y_oversampled_px": np.nan,
        "centroid_offset_native_px": np.nan,
        "ellipticity": np.nan,
        "position_angle_deg": np.nan,
        "rotation_asymmetry": np.nan,
        "positive_sum": total,
    }
    if data.ndim != 2 or total <= 0:
        return empty

    yy, xx = np.mgrid[: data.shape[0], : data.shape[1]]
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

    geometric_x = 0.5 * (data.shape[1] - 1)
    geometric_y = 0.5 * (data.shape[0] - 1)
    centered = image_shift(
        data,
        shift=(geometric_y - cy, geometric_x - cx),
        order=3,
        mode="nearest",
        prefilter=True,
    )
    denominator = 2.0 * np.sum(np.abs(centered))
    asymmetry = (
        float(np.sum(np.abs(centered - np.rot90(centered, 2))) / denominator)
        if denominator > 0
        else np.nan
    )
    return {
        "centroid_x_oversampled_px": cx,
        "centroid_y_oversampled_px": cy,
        "centroid_offset_native_px": float(
            np.hypot(cx - geometric_x, cy - geometric_y) / max(float(oversampling), 1e-12)
        ),
        "ellipticity": float(ellipticity),
        "position_angle_deg": float(np.degrees(np.arctan2(vector[1], vector[0]))),
        "rotation_asymmetry": asymmetry,
        "positive_sum": total,
    }


def _quantile_stats(values: pd.Series | np.ndarray) -> dict[str, float | int]:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return {
            "n": 0,
            "median": np.nan,
            "q16": np.nan,
            "q84": np.nan,
            "robust_scatter": np.nan,
            "rmse": np.nan,
            "outlier_fraction_0p2mag": np.nan,
        }
    q16, median, q84 = np.percentile(data, (16, 50, 84))
    return {
        "n": int(data.size),
        "median": float(median),
        "q16": float(q16),
        "q84": float(q84),
        "robust_scatter": float(1.4826 * np.median(np.abs(data - median))),
        "rmse": float(np.sqrt(np.mean(np.square(data)))),
        "outlier_fraction_0p2mag": float(np.mean(np.abs(data) > 0.2)),
    }


def _plot_group_quantiles(ax, data: pd.DataFrame, column: str, title: str) -> None:
    groups = list(data.groupby(column, observed=True, sort=True))
    if not groups:
        ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(title, loc="left")
        return
    positions = np.arange(len(groups))
    stats = [_quantile_stats(group["delta_centered_mag"]) for _, group in groups]
    medians = np.asarray([row["median"] for row in stats], dtype=float)
    lower = medians - np.asarray([row["q16"] for row in stats], dtype=float)
    upper = np.asarray([row["q84"] for row in stats], dtype=float) - medians
    labels = [f"{name}\nN={row['n']}" for (name, _), row in zip(groups, stats)]
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
    ax.set_xticks(positions, labels, fontsize=7)
    ax.set_ylabel("centered PSF - aperture (mag)")
    ax.set_title(title, loc="left")
    ax.grid(axis="y", alpha=0.18)


def _reference_offset(data: pd.DataFrame) -> tuple[float, float, int]:
    for threshold in (100.0, 50.0):
        reference = data.loc[data["snr_aperture"] >= threshold, "delta_mag"]
        reference = reference[np.isfinite(reference)]
        if len(reference) >= 20:
            return float(np.median(reference)), threshold, int(len(reference))
    finite = data.loc[np.isfinite(data["delta_mag"]), ["snr_aperture", "delta_mag"]]
    if finite.empty:
        return np.nan, np.nan, 0
    threshold = float(finite["snr_aperture"].quantile(0.8))
    reference = finite.loc[finite["snr_aperture"] >= threshold, "delta_mag"]
    return float(np.median(reference)), threshold, int(len(reference))


def _radial_bins(
    radius: pd.Series,
    *,
    core_radius: float,
    unit: str,
) -> pd.Series:
    finite = radius[np.isfinite(radius)]
    if finite.empty:
        return pd.Series(pd.Categorical([np.nan] * len(radius)), index=radius.index)

    if np.isfinite(core_radius) and core_radius > 0:
        edges = np.asarray([0.0, core_radius, 2.0 * core_radius, 4.0 * core_radius, np.inf])
        labels = [
            f"<{core_radius:.2g}{unit}",
            f"{core_radius:.2g}-{2 * core_radius:.2g}{unit}",
            f"{2 * core_radius:.2g}-{4 * core_radius:.2g}{unit}",
            f">={4 * core_radius:.2g}{unit}",
        ]
    else:
        quantiles = np.unique(np.nanpercentile(finite, (0, 25, 50, 75, 100)))
        if len(quantiles) < 3:
            return pd.Series(pd.Categorical(["all"] * len(radius)), index=radius.index)
        quantiles[0] = 0.0
        quantiles[-1] = np.inf
        edges = quantiles
        labels = [
            f"{left:.2g}-{right:.2g}{unit}" if np.isfinite(right) else f">={left:.2g}{unit}"
            for left, right in zip(edges[:-1], edges[1:])
        ]
    return pd.cut(radius, bins=edges, labels=labels, include_lowest=True)


def draw_psf_final_diagnostics(
    figure,
    data: pd.DataFrame,
    epsf_model: np.ndarray | None,
    *,
    filename: str,
    fwhm_px: float,
    pixel_scale_arcsec: float = np.nan,
    core_center: tuple[float, float] = (np.nan, np.nan),
    core_radius_px: float = np.nan,
    epsf_reference: pd.DataFrame | None = None,
) -> dict[str, object]:
    """Draw the six-panel final diagnostic and return JSON-safe summary values."""
    figure.clear()
    axes = figure.subplots(2, 3)
    if data is None or data.empty:
        for ax in axes.ravel():
            ax.set_axis_off()
        axes[0, 0].text(
            0.5,
            0.5,
            "No clean PSF/aperture matches for this frame.",
            ha="center",
            va="center",
            transform=axes[0, 0].transAxes,
        )
        figure.suptitle(f"{filename} PSF final diagnostics")
        return {"file": filename, "status": "CHECK", "warnings": ["no matched data"], "n_matched": 0}

    frame = data.copy()
    frame = frame.loc[
        np.isfinite(frame["snr_aperture"])
        & (frame["snr_aperture"] > 0)
        & np.isfinite(frame["delta_mag"])
    ].copy()
    reference_offset, reference_snr, reference_n = _reference_offset(frame)
    frame["delta_centered_mag"] = frame["delta_mag"] - reference_offset

    center_x, center_y = (float(core_center[0]), float(core_center[1]))
    if not np.isfinite(center_x):
        center_x = float(np.nanmedian(frame["x"]))
    if not np.isfinite(center_y):
        center_y = float(np.nanmedian(frame["y"]))
    frame["radius_px"] = np.hypot(frame["x"] - center_x, frame["y"] - center_y)

    finite_scale = np.isfinite(pixel_scale_arcsec) and pixel_scale_arcsec > 0
    finite_fwhm = np.isfinite(fwhm_px) and fwhm_px > 0
    if finite_scale:
        frame["radius_plot"] = frame["radius_px"] * float(pixel_scale_arcsec) / 60.0
        radial_unit = "'"
        radial_core = float(core_radius_px) * float(pixel_scale_arcsec) / 60.0
    elif finite_fwhm:
        frame["radius_plot"] = frame["radius_px"] / float(fwhm_px)
        radial_unit = " FWHM"
        radial_core = float(core_radius_px) / float(fwhm_px)
    else:
        frame["radius_plot"] = frame["radius_px"]
        radial_unit = " px"
        radial_core = float(core_radius_px)
    frame["radius_bin"] = _radial_bins(
        frame["radius_plot"],
        core_radius=radial_core,
        unit=radial_unit,
    )
    if finite_fwhm:
        frame["nearest_neighbor_fwhm"] = frame["nearest_neighbor_px"] / float(fwhm_px)
        frame["neighbor_bin"] = pd.cut(
            frame["nearest_neighbor_fwhm"],
            bins=[0, 1, 2, 4, np.inf],
            labels=["<1 FWHM", "1-2 FWHM", "2-4 FWHM", ">=4 FWHM"],
            include_lowest=True,
        )
    else:
        frame["neighbor_bin"] = _radial_bins(
            frame["nearest_neighbor_px"],
            core_radius=np.nan,
            unit=" px",
        )

    ax = axes[0, 0]
    ax.hexbin(
        frame["snr_aperture"],
        frame["delta_centered_mag"],
        xscale="log",
        gridsize=42,
        mincnt=1,
        cmap="Blues",
        bins="log",
    )
    ax.axhline(0.0, color="0.35", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Step 7 aperture SNR")
    ax.set_ylabel("centered PSF - aperture (mag)")
    ax.set_title("(a) low-SNR tail", loc="left")

    ax = axes[0, 1]
    low_snr = frame.loc[frame["snr_aperture"].between(5, 10, inclusive="left"), "delta_centered_mag"]
    high_snr = frame.loc[frame["snr_aperture"] >= 50, "delta_centered_mag"]
    hist_values = pd.concat([low_snr, high_snr]).replace([np.inf, -np.inf], np.nan).dropna()
    if hist_values.empty:
        bins = np.linspace(-0.5, 0.5, 30)
    else:
        lo, hi = np.nanpercentile(hist_values, (1, 99))
        margin = max(0.05, 0.08 * max(hi - lo, 0.1))
        bins = np.linspace(lo - margin, hi + margin, 34)
    if len(low_snr):
        ax.hist(low_snr, bins=bins, density=True, histtype="step", linewidth=1.5, label="SNR 5-10")
    if len(high_snr):
        ax.hist(high_snr, bins=bins, density=True, histtype="step", linewidth=1.5, label="SNR >= 50")
    ax.axvline(0.0, color="0.35", linewidth=0.8, linestyle="--")
    ax.set_xlabel("centered PSF - aperture (mag)")
    ax.set_ylabel("density")
    ax.set_title("(b) high-SNR convergence", loc="left")
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)

    _plot_group_quantiles(axes[0, 2], frame, "radius_bin", "(c) radial dependence")
    _plot_group_quantiles(axes[1, 0], frame, "neighbor_bin", "(d) crowding dependence")

    ax = axes[1, 1]
    centered = frame["delta_centered_mag"].to_numpy(dtype=float)
    finite_centered = centered[np.isfinite(centered)]
    color_limit = float(np.nanpercentile(np.abs(finite_centered), 90)) if finite_centered.size else 0.1
    color_limit = max(color_limit, 0.05)
    scatter = ax.scatter(
        frame["x"],
        frame["y"],
        c=centered,
        s=8,
        alpha=0.65,
        cmap="coolwarm",
        vmin=-color_limit,
        vmax=color_limit,
        linewidths=0,
    )
    if np.isfinite(core_radius_px) and core_radius_px > 0:
        ax.add_patch(Circle((center_x, center_y), float(core_radius_px), fill=False, color="black", linewidth=1.0))
    n_epsf_reference = 0
    if isinstance(epsf_reference, pd.DataFrame) and not epsf_reference.empty:
        reference_frame = epsf_reference.copy()
        if "selected" in reference_frame.columns:
            selected_reference = _as_bool(reference_frame["selected"])
            reference_frame = reference_frame.loc[selected_reference]
        if {"x", "y"} <= set(reference_frame.columns):
            ref_x = pd.to_numeric(reference_frame["x"], errors="coerce")
            ref_y = pd.to_numeric(reference_frame["y"], errors="coerce")
            finite_ref = np.isfinite(ref_x) & np.isfinite(ref_y)
            n_epsf_reference = int(np.count_nonzero(finite_ref))
            if n_epsf_reference:
                ax.scatter(
                    ref_x[finite_ref],
                    ref_y[finite_ref],
                    marker="s",
                    s=34,
                    facecolors="none",
                    edgecolors="black",
                    linewidths=0.8,
                    label="ePSF references",
                )
                ax.legend(loc="upper right", fontsize=7)
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")
    ax.set_title("(e) detector pattern + ePSF references", loc="left")
    figure.colorbar(scatter, ax=ax, label="centered delta mag", fraction=0.046, pad=0.04)

    metrics = epsf_shape_metrics(epsf_model) if epsf_model is not None else epsf_shape_metrics(np.asarray([]))
    fit_window_values = _column(frame, "fit_window_px")
    fit_window_values = fit_window_values[np.isfinite(fit_window_values)]
    fit_window_px = (
        float(np.median(fit_window_values)) if fit_window_values.size else np.nan
    )
    fit_energy_values = _column(frame, "fit_window_energy")
    fit_energy_values = fit_energy_values[np.isfinite(fit_energy_values)]
    fit_window_energy = (
        float(np.median(fit_energy_values)) if fit_energy_values.size else np.nan
    )
    qfit_noise = _column(frame, "qfit_noise_ratio")
    qfit_noise = qfit_noise[np.isfinite(qfit_noise)]
    median_qfit_noise = float(np.median(qfit_noise)) if qfit_noise.size else np.nan
    qfit_noise_gt3_fraction = (
        float(np.mean(qfit_noise > 3.0)) if qfit_noise.size else np.nan
    )
    ax = axes[1, 2]
    if epsf_model is None or np.asarray(epsf_model).ndim != 2:
        ax.text(0.5, 0.5, "ePSF model not found", ha="center", va="center", transform=ax.transAxes)
        ax.set_axis_off()
    else:
        image = ax.imshow(np.asarray(epsf_model, dtype=float), origin="lower", cmap="magma", interpolation="nearest")
        figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        ax.set_xlabel("oversampled x")
        ax.set_ylabel("oversampled y")
    ax.set_title(
        "(f) ePSF shape\n"
        f"ellipticity={metrics['ellipticity']:.3f}, A180={metrics['rotation_asymmetry']:.3f}\n"
        f"fit={fit_window_px:.0f}px, energy={fit_window_energy:.3f}",
        loc="left",
    )

    low_stats = _quantile_stats(low_snr)
    high_stats = _quantile_stats(high_snr)
    high_scatter = float(high_stats["robust_scatter"])
    warnings: list[str] = []
    if len(frame) < 30:
        warnings.append("fewer than 30 clean matches")
    if low_stats["n"] >= 5 and abs(float(low_stats["median"])) > 0.10:
        warnings.append("low-SNR median bias exceeds 0.10 mag")
    if low_stats["n"] >= 20 and float(low_stats["robust_scatter"]) > 0.30:
        warnings.append("low-SNR robust scatter exceeds 0.30 mag")
    if np.isfinite(high_scatter) and high_scatter > 0.10:
        warnings.append("high-SNR robust scatter exceeds 0.10 mag")
    if high_stats["n"] >= 20 and (
        float(high_stats["rmse"]) > 0.20
        or float(high_stats["outlier_fraction_0p2mag"]) > 0.05
    ):
        warnings.append("high-SNR non-Gaussian tail exceeds RMSE/outlier limits")
    if np.isfinite(metrics["ellipticity"]) and metrics["ellipticity"] > 0.20:
        warnings.append("ePSF ellipticity exceeds 0.20")
    if np.isfinite(metrics["rotation_asymmetry"]) and metrics["rotation_asymmetry"] > 0.15:
        warnings.append("ePSF A180 exceeds 0.15")
    if np.isfinite(median_qfit_noise) and median_qfit_noise > 1.5:
        warnings.append("median qfit/expected-noise exceeds 1.5")
    if np.isfinite(qfit_noise_gt3_fraction) and qfit_noise_gt3_fraction > 0.10:
        warnings.append("more than 10% of fits exceed qfit/expected-noise 3")

    figure.suptitle(
        f"{Path(filename).stem} PSF final diagnostics: N={len(frame)}, "
        f"high-SNR offset={reference_offset:+.3f} mag",
        fontsize=12,
    )
    figure.subplots_adjust(left=0.07, right=0.97, bottom=0.08, top=0.90, wspace=0.35, hspace=0.38)
    return {
        "file": filename,
        "status": "CHECK" if warnings else "OK",
        "warnings": warnings,
        "n_matched": int(len(frame)),
        "high_snr_reference_offset_mag": float(reference_offset),
        "high_snr_reference_threshold": float(reference_snr),
        "high_snr_reference_n": int(reference_n),
        "low_snr_5_10_n": int(low_stats["n"]),
        "low_snr_5_10_median_centered_mag": float(low_stats["median"]),
        "low_snr_5_10_robust_scatter_mag": float(low_stats["robust_scatter"]),
        "low_snr_5_10_rmse_mag": float(low_stats["rmse"]),
        "low_snr_5_10_outlier_fraction_0p2mag": float(
            low_stats["outlier_fraction_0p2mag"]
        ),
        "high_snr_n": int(high_stats["n"]),
        "high_snr_robust_scatter_mag": high_scatter,
        "high_snr_rmse_mag": float(high_stats["rmse"]),
        "high_snr_outlier_fraction_0p2mag": float(
            high_stats["outlier_fraction_0p2mag"]
        ),
        "epsf_ellipticity": float(metrics["ellipticity"]),
        "epsf_rotation_asymmetry": float(metrics["rotation_asymmetry"]),
        "epsf_reference_n": n_epsf_reference,
        "fit_window_px": fit_window_px,
        "fit_window_energy": fit_window_energy,
        "median_qfit_noise_ratio": median_qfit_noise,
        "qfit_noise_gt3_fraction": qfit_noise_gt3_fraction,
        "core_center_x_px": center_x,
        "core_center_y_px": center_y,
        "core_radius_px": float(core_radius_px),
        "fwhm_px": float(fwhm_px),
        "pixel_scale_arcsec": float(pixel_scale_arcsec),
    }
