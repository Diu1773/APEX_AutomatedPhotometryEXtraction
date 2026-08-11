"""Figure 14: synthetic crowded-field validation of APEX PSF fitting.

The experiment uses known source positions with small centroid/flux seed
errors. It tests fitting and deblending, not blind detection completeness.
"""

from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

from apex.gui.workflow.cmd.step8_psf_photometry import (  # noqa: E402
    _allstar_build_model,
    _allstar_fit,
)
from apex_paper_style import C, DOUBLE_COL, PALETTE, apply_paper_style, save_fig  # noqa: E402


SEED = 20260713
IMAGE_SHAPE = (256, 256)
SOURCE_COUNTS = (80, 160, 320, 640)
SIGMA_PX = 1.4
FWHM_PX = 2.3548 * SIGMA_PX
BACKGROUND_RMS = 5.0
GAIN = 1.0
FIT_SHAPE = 9
STAMP_SIZE = 21

OUTDIR = REPO / "validation" / "paper" / "figures"
DATADIR = REPO / "validation" / "paper" / "data" / "psf_iteration"


def _psf(dx, dy):
    values = np.exp(
        -0.5 * (np.asarray(dx, dtype=float) ** 2 + np.asarray(dy, dtype=float) ** 2)
        / SIGMA_PX**2
    )
    return values / (2.0 * np.pi * SIGMA_PX**2)


def _nearest_distance(xy: np.ndarray) -> np.ndarray:
    distances, _ = cKDTree(xy).query(xy, k=2, workers=1)
    return np.asarray(distances[:, 1], dtype=float)


def _aperture_flux(image: np.ndarray, xy: np.ndarray) -> np.ndarray:
    output = np.full(len(xy), np.nan, dtype=float)
    yy_full, xx_full = np.indices(image.shape)
    for index, (x_center, y_center) in enumerate(xy):
        radius = np.hypot(xx_full - x_center, yy_full - y_center)
        aperture = radius <= 2.0 * FWHM_PX
        annulus = (radius >= 2.5 * FWHM_PX) & (radius <= 3.5 * FWHM_PX)
        sky = float(np.median(image[annulus])) if np.any(annulus) else 0.0
        output[index] = float(np.sum(image[aperture] - sky))
    return output


def _run_field(n_sources: int, rng: np.random.Generator) -> tuple[pd.DataFrame, dict]:
    height, width = IMAGE_SHAPE
    xy_true = np.column_stack([
        rng.uniform(12.0, width - 12.0, n_sources),
        rng.uniform(12.0, height - 12.0, n_sources),
    ])
    flux_true = 10.0 ** rng.uniform(2.7, 4.5, n_sources)
    noiseless = _allstar_build_model(
        IMAGE_SHAPE,
        xy_true[:, 0],
        xy_true[:, 1],
        flux_true,
        _psf,
        STAMP_SIZE,
    )
    image = noiseless + rng.normal(0.0, BACKGROUND_RMS, IMAGE_SHAPE)
    seed_xy = xy_true + rng.normal(0.0, 0.15, xy_true.shape)
    seed_flux = flux_true * rng.uniform(0.85, 1.15, n_sources)

    started = time.perf_counter()
    free = _allstar_fit(
        image,
        seed_xy,
        seed_flux,
        _psf,
        fit_shape=FIT_SHAPE,
        stamp_size=STAMP_SIZE,
        max_iter=4,
        flux_conv=0.01,
        max_shift=3.0,
        background_rms=BACKGROUND_RMS,
        gain=GAIN,
        initial_positions=seed_xy,
        position_bound=3.0,
    )
    free_xy = np.column_stack([free["x_fit"], free["y_fit"]]).astype(float)
    fixed = _allstar_fit(
        image,
        free_xy,
        np.asarray(free["flux_fit"], dtype=float),
        _psf,
        fit_shape=FIT_SHAPE,
        stamp_size=STAMP_SIZE,
        max_iter=2,
        flux_conv=0.01,
        max_shift=3.0,
        background_rms=BACKGROUND_RMS,
        gain=GAIN,
        initial_positions=free_xy,
        position_fixed=True,
    )
    elapsed = time.perf_counter() - started

    fit_flux = np.asarray(fixed["flux_fit"], dtype=float)
    fit_flags = (
        np.asarray(free["flags"], dtype=np.int32)
        | np.asarray(fixed["flags"], dtype=np.int32)
    )
    aperture_flux = _aperture_flux(image, xy_true)
    mag_error_psf = -2.5 * np.log10(np.clip(fit_flux / flux_true, 1e-12, None))
    mag_error_aperture = -2.5 * np.log10(
        np.clip(aperture_flux / flux_true, 1e-12, None)
    )
    nearest = _nearest_distance(xy_true)
    fitted_model = _allstar_build_model(
        IMAGE_SHAPE,
        np.asarray(fixed["x_fit"], dtype=float),
        np.asarray(fixed["y_fit"], dtype=float),
        fit_flux,
        _psf,
        STAMP_SIZE,
    )
    residual_rms = float(np.std(image - fitted_model))

    rows = pd.DataFrame({
        "n_sources": n_sources,
        "source_density_per_1e4px": n_sources / (height * width) * 1e4,
        "nearest_px": nearest,
        "flux_true": flux_true,
        "mag_error_psf": mag_error_psf,
        "mag_error_aperture": mag_error_aperture,
        "qfit": np.asarray(fixed["qfit"], dtype=float),
        "cfit": np.asarray(fixed["cfit"], dtype=float),
        "reduced_chi2": np.asarray(fixed["reduced_chi2"], dtype=float),
        "flags": fit_flags,
    })
    clean = (rows["flags"] == 0) & np.isfinite(rows["mag_error_psf"])
    summary = {
        "n_sources": n_sources,
        "density_per_1e4px": float(rows["source_density_per_1e4px"].iloc[0]),
        "elapsed_s": elapsed,
        "clean_fraction": float(np.mean(clean)),
        "recovery_fraction_0p1mag": float(np.mean(clean & (np.abs(rows["mag_error_psf"]) < 0.1))),
        "median_abs_error_psf": float(np.median(np.abs(rows.loc[clean, "mag_error_psf"]))),
        "median_abs_error_aperture": float(np.median(np.abs(rows["mag_error_aperture"]))),
        "residual_rms": residual_rms,
    }
    return rows, summary


def _binned_median(x: np.ndarray, y: np.ndarray, edges: np.ndarray):
    x_mid, median, scatter = [], [], []
    for low, high in zip(edges[:-1], edges[1:]):
        keep = np.isfinite(x) & np.isfinite(y) & (x >= low) & (x < high)
        if np.sum(keep) < 8:
            continue
        values = y[keep]
        center = float(np.median(values))
        mad = 1.4826 * float(np.median(np.abs(values - center)))
        x_mid.append(float(np.median(x[keep])))
        median.append(center)
        scatter.append(mad / np.sqrt(np.sum(keep)))
    return np.asarray(x_mid), np.asarray(median), np.asarray(scatter)


def main() -> int:
    apply_paper_style()
    rng = np.random.default_rng(SEED)
    row_tables = []
    summaries = []
    for count in SOURCE_COUNTS:
        rows, summary = _run_field(count, rng)
        row_tables.append(rows)
        summaries.append(summary)
        print(json.dumps(summary, sort_keys=True))

    all_rows = pd.concat(row_tables, ignore_index=True)
    summary_df = pd.DataFrame(summaries)
    DATADIR.mkdir(parents=True, exist_ok=True)
    all_rows.to_csv(DATADIR / "source_metrics.csv", index=False)
    summary_df.to_csv(DATADIR / "summary.csv", index=False)
    (DATADIR / "run.json").write_text(
        json.dumps(
            {
                "seed": SEED,
                "image_shape": IMAGE_SHAPE,
                "source_counts": SOURCE_COUNTS,
                "fwhm_px": FWHM_PX,
                "background_rms": BACKGROUND_RMS,
                "one_core": True,
            },
            indent=2,
        ),
        encoding="utf-8",
    )

    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 6.1))
    ax_bias, ax_qfit, ax_recovery, ax_time = axes.ravel()

    edges = np.array([0, 3, 4, 5, 7, 10, 15, 30, np.inf])
    for column, label, color, marker in (
        ("mag_error_psf", "APEX PSF", C["data"], "o"),
        ("mag_error_aperture", "aperture", C["accent"], "s"),
    ):
        x_mid, median, error = _binned_median(
            all_rows["nearest_px"].to_numpy(float),
            all_rows[column].to_numpy(float),
            edges,
        )
        ax_bias.errorbar(x_mid, median, yerr=error, marker=marker, color=color, label=label)
    ax_bias.axhline(0, color=PALETTE["grey"], lw=0.8, ls=":")
    ax_bias.set_xscale("log")
    ax_bias.set_xlabel("nearest-neighbour separation (px)")
    ax_bias.set_ylabel("recovered - true magnitude")
    ax_bias.set_title("(a) blending bias", loc="left")
    ax_bias.legend()

    clean_rows = all_rows[all_rows["flags"] == 0]
    ax_qfit.scatter(
        clean_rows["qfit"],
        np.abs(clean_rows["mag_error_psf"]),
        s=7,
        alpha=0.25,
        color=C["data"],
        edgecolors="none",
    )
    ax_qfit.set_xscale("log")
    ax_qfit.set_yscale("log")
    ax_qfit.set_xlabel("qfit")
    ax_qfit.set_ylabel("absolute magnitude error")
    ax_qfit.set_title("(b) fit metric vs truth", loc="left")

    ax_recovery.plot(
        summary_df["density_per_1e4px"],
        summary_df["recovery_fraction_0p1mag"],
        "o-",
        color=C["data"],
        label=r"clean and $|\Delta m|<0.1$",
    )
    ax_recovery.plot(
        summary_df["density_per_1e4px"],
        summary_df["clean_fraction"],
        "s--",
        color=C["accent"],
        label="flags = 0",
    )
    ax_recovery.set_ylim(0, 1.03)
    ax_recovery.set_xlabel(r"source density ($10^{-4}$ px$^{-2}$)")
    ax_recovery.set_ylabel("fraction")
    ax_recovery.set_title("(c) clean recovery", loc="left")
    ax_recovery.legend()

    ax_time.plot(
        summary_df["n_sources"],
        summary_df["elapsed_s"],
        "o-",
        color=PALETTE["green"],
    )
    ax_time.set_xlabel("fitted sources")
    ax_time.set_ylabel("wall time (s, one CPU core)")
    ax_time.set_title("(d) CPU scaling", loc="left")

    fig.suptitle("APEX iterative PSF photometry: synthetic crowded fields", y=1.01)
    fig.tight_layout()
    paths = save_fig(fig, "fig14_psf_iteration", OUTDIR)
    plt.close(fig)
    print("saved:", ", ".join(str(path) for path in paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
