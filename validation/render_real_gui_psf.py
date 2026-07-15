"""Render actual Step 8 FITS snapshots from a GUI-worker validation run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


def _load(path: Path) -> np.ndarray:
    return np.asarray(fits.getdata(path), dtype=float)


def _limits(image: np.ndarray, high: float = 99.8) -> tuple[float, float]:
    finite = image[np.isfinite(image)]
    return float(np.percentile(finite, 1.0)), float(np.percentile(finite, high))


def _show_signal(ax, image, title, *, extent=None):
    low, high = _limits(image)
    display = np.arcsinh(np.clip(image - low, 0.0, None))
    ax.imshow(
        display,
        origin="lower",
        cmap="gray",
        vmin=0.0,
        vmax=np.arcsinh(max(high - low, 1.0)),
        interpolation="nearest",
        extent=extent,
    )
    ax.set_title(title, loc="left")


def _show_residual(ax, image, title, *, extent=None):
    finite = image[np.isfinite(image)]
    limit = float(np.percentile(np.abs(finite), 99.5))
    ax.imshow(
        image,
        origin="lower",
        cmap="coolwarm",
        vmin=-limit,
        vmax=limit,
        interpolation="nearest",
        extent=extent,
    )
    ax.set_title(title, loc="left")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    psf_dir = args.run_dir / "cmd_psf"
    meta_path = next(psf_dir.glob("residual_meta_*.json"))
    metadata = json.loads(meta_path.read_text(encoding="utf-8"))
    snapshots = metadata["iters"]
    first = snapshots[0]
    final = snapshots[-1]

    model_first = _load(psf_dir / first["model_path"])
    residual_first = _load(psf_dir / first["residual_path"])
    residual_final = _load(psf_dir / final["residual_path"])
    sky_subtracted = model_first + residual_first
    candidate_xy = np.load(psf_dir / first["candidatexy_path"])
    candidate_pruned = any(int(row.get("n_pruned", 0)) > 0 for row in snapshots)
    if len(candidate_xy):
        zoom_x, zoom_y = map(float, candidate_xy[0])
    else:
        maximum = np.unravel_index(np.nanargmax(residual_first), residual_first.shape)
        zoom_y, zoom_x = map(float, maximum)

    half = 35
    height, width = sky_subtracted.shape
    x0 = max(0, int(np.floor(zoom_x - half)))
    x1 = min(width, int(np.ceil(zoom_x + half)))
    y0 = max(0, int(np.floor(zoom_y - half)))
    y1 = min(height, int(np.ceil(zoom_y + half)))
    zoom = np.s_[y0:y1, x0:x1]
    extent = (x0, x1, y0, y1)

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 8.2))
    fig.patch.set_facecolor("white")
    for ax in axes.ravel():
        ax.set_facecolor("white")
    _show_signal(axes[0, 0], sky_subtracted, "(a) real M60 frame: sky subtracted")
    _show_signal(axes[0, 1], model_first, "(b) model: 99 Step 4 sources")
    _show_residual(axes[0, 2], residual_first, "(c) residual after first PSF fit")
    if len(candidate_xy):
        axes[0, 2].scatter(
            candidate_xy[:, 0], candidate_xy[:, 1], s=70, facecolors="none",
            edgecolors="#009E73", linewidths=1.4,
            label=("trial candidate, later rejected" if candidate_pruned else "accepted residual source")
        )
        axes[0, 2].legend(loc="upper right", fontsize=8)

    _show_signal(
        axes[1, 0], sky_subtracted[zoom], "(d) trial candidate: original",
        extent=extent
    )
    _show_residual(
        axes[1, 1], residual_first[zoom], "(e) candidate in first residual",
        extent=extent
    )
    _show_residual(
        axes[1, 2], residual_final[zoom],
        ("(f) final residual after QC rejection" if candidate_pruned else "(f) residual after joint refit"),
        extent=extent
    )
    if len(candidate_xy):
        for ax in axes[1, :2]:
            ax.scatter(
                candidate_xy[:, 0], candidate_xy[:, 1], s=90,
                facecolors="none", edgecolors="#009E73", linewidths=1.5
            )

    for ax in axes.ravel():
        ax.set_xlabel("x (px)")
        ax.set_ylabel("y (px)")
    fig.suptitle(
        "Actual GUI-worker run: residual detection and iterative PSF refit",
        fontsize=17,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    output = args.output or (args.run_dir / "real_psf_residual_sequence.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output,
        dpi=220,
        bbox_inches="tight",
        facecolor="white",
        transparent=False,
    )
    plt.close(fig)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
