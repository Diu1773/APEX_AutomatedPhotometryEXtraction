"""Figure 15: visual demonstration of residual-source PSF iteration.

Bright sources are supplied as the initial catalogue. Fainter sources are
deliberately omitted, recovered from the first residual image, and included in
a joint second PSF fit. This isolates the residual-detection workflow from WCS
and Step 4 catalogue completeness.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import matplotlib.pyplot as plt
import numpy as np
from photutils.detection import DAOStarFinder
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

from apex.gui.workflow.cmd.step8_psf_photometry import (  # noqa: E402
    _allstar_build_model,
    _allstar_fit,
)
from apex_paper_style import DOUBLE_COL, apply_paper_style, save_fig  # noqa: E402


SEED = 20260714
IMAGE_SHAPE = (256, 256)
SIGMA_PX = 1.4
FWHM_PX = 2.3548 * SIGMA_PX
BACKGROUND_RMS = 5.0
FIT_SHAPE = 9
STAMP_SIZE = 21
OUTDIR = REPO / "validation" / "paper" / "figures"
DATADIR = REPO / "validation" / "paper" / "data" / "psf_residual_sequence"


def _psf(dx, dy):
    values = np.exp(
        -0.5
        * (np.asarray(dx, dtype=float) ** 2 + np.asarray(dy, dtype=float) ** 2)
        / SIGMA_PX**2
    )
    return values / (2.0 * np.pi * SIGMA_PX**2)


def _fit(image, xy, flux, *, position_fixed=False):
    result = _allstar_fit(
        image,
        xy,
        flux,
        _psf,
        fit_shape=FIT_SHAPE,
        stamp_size=STAMP_SIZE,
        max_iter=2 if position_fixed else 4,
        flux_conv=0.01,
        max_shift=3.0,
        background_rms=BACKGROUND_RMS,
        gain=1.0,
        initial_positions=xy,
        position_bound=3.0,
        position_fixed=position_fixed,
    )
    result_xy = np.column_stack([result["x_fit"], result["y_fit"]]).astype(float)
    result_flux = np.asarray(result["flux_fit"], dtype=float)
    return result, result_xy, result_flux


def _display_limits(image):
    finite = np.asarray(image, dtype=float)[np.isfinite(image)]
    return float(np.percentile(finite, 1.0)), float(np.percentile(finite, 99.8))


def _show(ax, image, title, *, residual=False, extent=None):
    if residual:
        limit = float(np.percentile(np.abs(image[np.isfinite(image)]), 99.5))
        ax.imshow(
            image,
            origin="lower",
            cmap="coolwarm",
            vmin=-limit,
            vmax=limit,
            interpolation="nearest",
            extent=extent,
        )
    else:
        low, high = _display_limits(image)
        ax.imshow(
            np.arcsinh(np.clip(image - low, 0.0, None)),
            origin="lower",
            cmap="gray",
            vmin=0.0,
            vmax=np.arcsinh(max(high - low, 1.0)),
            interpolation="nearest",
            extent=extent,
        )
    ax.set_title(title, loc="left")
    ax.set_xlabel("x (px)")
    ax.set_ylabel("y (px)")


def main() -> int:
    apply_paper_style()
    rng = np.random.default_rng(SEED)
    height, width = IMAGE_SHAPE

    n_initial = 55
    initial_true_xy = np.column_stack(
        [
            rng.uniform(14.0, width - 14.0, n_initial),
            rng.uniform(14.0, height - 14.0, n_initial),
        ]
    )
    initial_flux = 10.0 ** rng.uniform(3.6, 4.3, n_initial)

    companion_indices = rng.choice(n_initial, size=18, replace=False)
    angle = rng.uniform(0.0, 2.0 * np.pi, len(companion_indices))
    separation = rng.uniform(3.8, 6.0, len(companion_indices))
    companion_xy = initial_true_xy[companion_indices] + np.column_stack(
        [np.cos(angle) * separation, np.sin(angle) * separation]
    )
    isolated_xy = np.column_stack(
        [rng.uniform(14.0, width - 14.0, 12), rng.uniform(14.0, height - 14.0, 12)]
    )
    omitted_xy = np.vstack([companion_xy, isolated_xy])
    omitted_flux = 10.0 ** rng.uniform(2.65, 3.25, len(omitted_xy))

    truth_xy = np.vstack([initial_true_xy, omitted_xy])
    truth_flux = np.concatenate([initial_flux, omitted_flux])
    noiseless = _allstar_build_model(
        IMAGE_SHAPE,
        truth_xy[:, 0],
        truth_xy[:, 1],
        truth_flux,
        _psf,
        STAMP_SIZE,
    )
    image = noiseless + rng.normal(0.0, BACKGROUND_RMS, IMAGE_SHAPE)

    seed_xy = initial_true_xy + rng.normal(0.0, 0.12, initial_true_xy.shape)
    seed_flux = initial_flux * rng.uniform(0.9, 1.1, n_initial)
    _, first_xy, first_flux = _fit(image, seed_xy, seed_flux)
    first_model = _allstar_build_model(
        IMAGE_SHAPE,
        first_xy[:, 0],
        first_xy[:, 1],
        first_flux,
        _psf,
        STAMP_SIZE,
    )
    first_residual = image - first_model

    detections = DAOStarFinder(
        fwhm=FWHM_PX,
        threshold=4.0 * BACKGROUND_RMS,
        sharplo=0.1,
        sharphi=1.0,
        roundlo=-0.9,
        roundhi=0.9,
        exclude_border=True,
    )(first_residual)
    if detections is None:
        raise RuntimeError("Residual detector returned no sources")
    detected_xy = np.column_stack([detections["xcentroid"], detections["ycentroid"]]).astype(float)
    detected_flux = np.asarray(detections["flux"], dtype=float)

    distance_to_initial, _ = cKDTree(first_xy).query(detected_xy, k=1, workers=1)
    unique = distance_to_initial >= 2.0
    detected_xy = detected_xy[unique]
    detected_flux = detected_flux[unique]

    distance_to_omitted, omitted_match = cKDTree(omitted_xy).query(
        detected_xy, k=1, workers=1
    )
    matched = distance_to_omitted <= 1.5
    recovered_truth = np.unique(omitted_match[matched])

    combined_xy = np.vstack([first_xy, detected_xy])
    combined_flux = np.concatenate([first_flux, np.clip(detected_flux, 1.0, None)])
    _, second_xy, second_flux = _fit(image, combined_xy, combined_flux)
    _, final_xy, final_flux = _fit(
        image, second_xy, second_flux, position_fixed=True
    )
    final_model = _allstar_build_model(
        IMAGE_SHAPE,
        final_xy[:, 0],
        final_xy[:, 1],
        final_flux,
        _psf,
        STAMP_SIZE,
    )
    final_residual = image - final_model

    DATADIR.mkdir(parents=True, exist_ok=True)
    summary = {
        "seed": SEED,
        "initial_sources": int(n_initial),
        "deliberately_omitted": int(len(omitted_xy)),
        "residual_candidates": int(len(detected_xy)),
        "matched_candidates": int(np.sum(matched)),
        "unique_omitted_recovered": int(len(recovered_truth)),
        "recovery_fraction": float(len(recovered_truth) / len(omitted_xy)),
        "false_candidate_fraction": float(np.mean(~matched)) if len(matched) else 0.0,
        "first_residual_rms": float(np.std(first_residual)),
        "final_residual_rms": float(np.std(final_residual)),
    }
    (DATADIR / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))

    fig, axes = plt.subplots(2, 3, figsize=(DOUBLE_COL, 5.4))
    _show(axes[0, 0], image, "(a) original: 85 true sources")
    _show(axes[0, 1], first_model, "(b) model: 55 initial sources")
    _show(axes[0, 2], first_residual, "(c) subtraction residual", residual=True)
    axes[0, 2].scatter(
        detected_xy[:, 0], detected_xy[:, 1], s=26, facecolors="none",
        edgecolors="#009E73", linewidths=0.8, label="residual detections"
    )
    axes[0, 2].legend(loc="upper right", fontsize=7)

    best = int(np.argmax(omitted_flux[: len(companion_xy)]))
    zoom_center = companion_xy[best]
    half = 13
    x0 = max(0, int(np.floor(zoom_center[0] - half)))
    x1 = min(width, int(np.ceil(zoom_center[0] + half)))
    y0 = max(0, int(np.floor(zoom_center[1] - half)))
    y1 = min(height, int(np.ceil(zoom_center[1] + half)))
    extent = (x0, x1, y0, y1)
    zoom = np.s_[y0:y1, x0:x1]

    _show(axes[1, 0], image[zoom], "(d) blended pair: original", extent=extent)
    axes[1, 0].scatter(
        omitted_xy[best, 0], omitted_xy[best, 1], marker="+", s=70,
        color="#D55E00", linewidths=1.2, label="omitted companion"
    )
    axes[1, 0].legend(loc="upper right", fontsize=7)
    _show(
        axes[1, 1], first_residual[zoom], "(e) companion after subtraction",
        residual=True, extent=extent
    )
    in_zoom = (
        (detected_xy[:, 0] >= x0) & (detected_xy[:, 0] < x1)
        & (detected_xy[:, 1] >= y0) & (detected_xy[:, 1] < y1)
    )
    axes[1, 1].scatter(
        detected_xy[in_zoom, 0], detected_xy[in_zoom, 1], s=50,
        facecolors="none", edgecolors="#009E73", linewidths=1.2
    )
    _show(
        axes[1, 2], final_residual[zoom], "(f) residual after joint refit",
        residual=True, extent=extent
    )

    fig.suptitle("Residual detection and iterative PSF refitting")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    saved = save_fig(fig, "fig15_psf_residual_sequence", OUTDIR)
    print("saved:", ", ".join(sorted(saved)))
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
