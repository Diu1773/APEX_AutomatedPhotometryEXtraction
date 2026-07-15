"""Render PSF validation figures for the completed reprocessing projects."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.spatial import cKDTree


PROJECTS = ("NGC6811", "M67", "M13")
SNR_CUTS = (5, 10, 20, 50, 100)


def _load_fits(path: Path) -> np.ndarray:
    return np.asarray(fits.getdata(path), dtype=np.float32)


def _finite_percentile(values: np.ndarray, percentile: float) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 1.0
    return float(np.percentile(finite, percentile))


def _show_signal(ax, image: np.ndarray, title: str) -> None:
    low = _finite_percentile(image, 1.0)
    high = _finite_percentile(image, 99.8)
    scaled = np.arcsinh(np.clip(image - low, 0.0, None))
    ax.imshow(
        scaled,
        origin="lower",
        cmap="gray",
        vmin=0.0,
        vmax=np.arcsinh(max(high - low, 1.0)),
        interpolation="nearest",
    )
    ax.set_title(title, loc="left", fontsize=11)


def _show_residual(ax, image: np.ndarray, title: str) -> None:
    limit = _finite_percentile(np.abs(image), 99.5)
    ax.imshow(
        image,
        origin="lower",
        cmap="coolwarm",
        vmin=-max(limit, 1.0),
        vmax=max(limit, 1.0),
        interpolation="nearest",
    )
    ax.set_title(title, loc="left", fontsize=11)


def _read_project(root: Path, label: str) -> dict:
    result_dir = root / label / "result"
    psf_dir = result_dir / "cmd_psf"
    meta_path = next(psf_dir.glob("residual_meta_*.json"))
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    index = pd.read_csv(psf_dir / "photometry_index.csv").iloc[0]
    first = meta["iters"][0]
    final = meta["iters"][-1]
    model = _load_fits(psf_dir / first["model_path"])
    first_residual = _load_fits(psf_dir / first["residual_path"])
    return {
        "label": label,
        "result_dir": result_dir,
        "psf_dir": psf_dir,
        "meta": meta,
        "index": index,
        "science": model + first_residual,
        "model": model,
        "residual": _load_fits(psf_dir / final["residual_path"]),
    }


def _overlay_core(ax, core: dict) -> None:
    if not core.get("enabled"):
        return
    circle = Circle(
        (float(core["center_x"]), float(core["center_y"])),
        float(core["radius_px"]),
        fill=False,
        edgecolor="#F0E442",
        linewidth=1.4,
    )
    ax.add_patch(circle)


def _render_overview(projects: list[dict], output: Path) -> None:
    fig, axes = plt.subplots(3, 3, figsize=(14.5, 10.0), constrained_layout=True)
    for row, project in enumerate(projects):
        label = project["label"]
        index = project["index"]
        _show_signal(axes[row, 0], project["science"], f"{label}: sky-subtracted frame")
        _show_signal(axes[row, 1], project["model"], "first-pass PSF model")
        _show_residual(axes[row, 2], project["residual"], "final residual")
        core = project["meta"].get("core_cut", {})
        for ax in axes[row]:
            _overlay_core(ax, core)
            ax.set_xticks([])
            ax.set_yticks([])
        annotation = (
            f"Step4={int(index['epsf_n_detected'])}  "
            f"fit={int(index['n'])}  good={int(index['n_goodmag'])}  "
            f"new={int(index['n_new_iter'])}\n"
            f"median qfit={float(index['median_qfit']):.3f}  "
            f"reduced chi2={float(index['median_reduced_chi2']):.2f}"
        )
        if core.get("enabled"):
            annotation += (
                f"  core r={float(core['radius_px']):.0f}px"
                f"  excluded={int(core['n_excluded_init'])}"
            )
        axes[row, 0].text(
            0.012,
            0.018,
            annotation,
            transform=axes[row, 0].transAxes,
            color="white",
            fontsize=8.5,
            va="bottom",
            bbox={"facecolor": "black", "alpha": 0.62, "edgecolor": "none", "pad": 3},
        )
    fig.suptitle("APEX iterative PSF photometry on completed reprocessing projects", fontsize=16)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=190, facecolor="white")
    plt.close(fig)


def _as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "yes"})


def _match_aperture(project: dict) -> pd.DataFrame:
    meta = project["meta"]
    filename = meta["file"]
    psf_path = project["psf_dir"] / f"photometry_{filename}.tsv"
    aperture_path = (
        project["result_dir"]
        / "step7_forced_phot"
        / f"photometry_{filename}.tsv"
    )
    psf = pd.read_csv(psf_path, sep="\t")
    aperture = pd.read_csv(aperture_path, sep="\t", low_memory=False)

    for frame, columns in (
        (psf, ("x_fit", "y_fit", "flux_psf_e", "snr_psf", "flags_psf")),
        (aperture, ("x_fit", "y_fit", "flux_e", "snr")),
    ):
        for column in columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

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
    tree = cKDTree(aperture[["x_fit", "y_fit"]].to_numpy())
    distance, aperture_index = tree.query(psf[["x_fit", "y_fit"]].to_numpy(), k=1)
    matched = pd.DataFrame(
        {
            "psf_index": np.arange(len(psf)),
            "aperture_index": aperture_index,
            "distance_px": distance,
        }
    )
    matched = matched.loc[matched["distance_px"] <= 1.5]
    matched = matched.sort_values("distance_px").drop_duplicates("aperture_index")
    p = psf.iloc[matched["psf_index"].to_numpy()].reset_index(drop=True)
    a = aperture.iloc[matched["aperture_index"].to_numpy()].reset_index(drop=True)
    flux_ratio = p["flux_psf_e"].to_numpy() / a["flux_e"].to_numpy()
    return pd.DataFrame(
        {
            "snr_aperture": a["snr"].to_numpy(),
            "mag_aperture": pd.to_numeric(a["mag_inst"], errors="coerce").to_numpy(),
            "delta_mag": -2.5 * np.log10(flux_ratio),
            "distance_px": matched["distance_px"].to_numpy(),
        }
    )


def _robust_stats(values: np.ndarray) -> tuple[float, float, float]:
    median = float(np.median(values))
    scatter = float(1.4826 * np.median(np.abs(values - median)))
    p90 = float(np.percentile(np.abs(values - median), 90.0))
    return median, scatter, p90


def _render_comparison(projects: list[dict], output: Path, stats_path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.7), sharey=False, constrained_layout=True)
    stats_rows = []
    colors = plt.colormaps["viridis"](np.linspace(0.15, 0.9, len(SNR_CUTS)))
    rng = np.random.default_rng(20260714)

    for ax, project in zip(axes, projects):
        matched = _match_aperture(project)
        matched = matched.replace([np.inf, -np.inf], np.nan).dropna(subset=["snr_aperture", "delta_mag"])
        draw = matched
        if len(draw) > 3500:
            draw = draw.iloc[rng.choice(len(draw), 3500, replace=False)]
        ax.scatter(
            draw["snr_aperture"],
            draw["delta_mag"],
            s=7,
            color="#4C78A8",
            alpha=0.22,
            linewidths=0,
        )
        ax.set_xscale("log")
        ax.axhline(0.0, color="0.45", linewidth=0.8, linestyle="--")

        y_values = matched["delta_mag"].to_numpy()
        low, high = np.percentile(y_values, (1.0, 99.0))
        pad = max((high - low) * 0.12, 0.03)
        ax.set_ylim(low - pad, high + pad)
        for color, cut in zip(colors, SNR_CUTS):
            selected = matched.loc[matched["snr_aperture"] >= cut, "delta_mag"].to_numpy()
            if len(selected) < 5:
                continue
            median, scatter, p90 = _robust_stats(selected)
            stats_rows.append(
                {
                    "project": project["label"],
                    "snr_cut": cut,
                    "n_match": len(selected),
                    "median_psf_minus_aperture_mag": median,
                    "robust_scatter_mag": scatter,
                    "p90_abs_centered_mag": p90,
                }
            )
            ax.scatter(cut, median, s=48, color=color, edgecolor="black", linewidth=0.45, zorder=3)
            ax.errorbar(cut, median, yerr=scatter, color=color, linewidth=1.2, capsize=2, zorder=2)

        all_median, all_scatter, _ = _robust_stats(y_values)
        ax.set_title(
            f"{project['label']}  matched={len(matched)}\n"
            f"all: median={all_median:+.3f}, robust scatter={all_scatter:.3f} mag",
            fontsize=10.5,
        )
        ax.set_xlabel("Step7 aperture SNR")
        ax.set_ylabel("PSF - aperture instrumental magnitude")
        ax.grid(alpha=0.18, linewidth=0.6)

    fig.suptitle("PSF-to-aperture consistency versus SNR", fontsize=15)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=200, facecolor="white")
    plt.close(fig)
    pd.DataFrame(stats_rows).to_csv(stats_path, index=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(r"E:\APEX_validation\reprocess"))
    parser.add_argument("--output-dir", type=Path, default=Path("validation/real_gui_run"))
    args = parser.parse_args()

    projects = [_read_project(args.root, label) for label in PROJECTS]
    overview = args.output_dir / "reprocess_psf_overview.png"
    comparison = args.output_dir / "reprocess_psf_aperture_comparison.png"
    stats = args.output_dir / "reprocess_psf_aperture_stats.csv"
    _render_overview(projects, overview)
    _render_comparison(projects, comparison, stats)
    print(overview.resolve())
    print(comparison.resolve())
    print(stats.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
