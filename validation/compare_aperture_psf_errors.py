"""Compare aperture and current APEX PSF photometry with reported errors.

The comparison is an agreement and error-calibration diagnostic. Aperture
photometry is not treated as truth. A constant bright, isolated-star offset is
removed before evaluating the PSF-minus-aperture residual distribution.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm


REPO = Path(__file__).absolute().parents[1]
DEFAULT_OUTPUT = REPO / "validation" / "psf_archive" / "aperture_psf"

SNR_EDGES = np.asarray([5.0, 10.0, 20.0, 50.0, 100.0, np.inf])
SNR_LABELS = ("5-10", "10-20", "20-50", "50-100", ">=100")
CROWD_EDGES = np.asarray([0.0, 1.5, 3.0, 6.0, np.inf])
CROWD_LABELS = ("<1.5", "1.5-3", "3-6", ">=6")


@dataclass(frozen=True)
class ClusterSpec:
    slug: str
    title: str
    cluster_type: str
    filename: str
    psf_result: Path
    aperture_result: Path


def _numeric(frame: pd.DataFrame, columns: tuple[str, ...]) -> None:
    for column in columns:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")


def _as_bool(series: pd.Series) -> np.ndarray:
    return series.astype(str).str.strip().str.lower().isin(
        {"true", "1", "yes", "y"}
    ).to_numpy(dtype=bool)


def _robust_scatter(values) -> float:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size < 2:
        return np.nan
    median = float(np.median(data))
    return float(1.4826 * np.median(np.abs(data - median)))


def _read_fwhm(result_dir: Path, filename: str) -> float:
    candidates = (
        result_dir / "cache" / f"detect_{filename}.json",
        result_dir / "step4_detection" / f"detect_{filename}.json",
    )
    for path in candidates:
        if not path.exists():
            continue
        metadata = json.loads(path.read_text(encoding="utf-8"))
        value = float(metadata.get("fwhm_px", np.nan))
        if np.isfinite(value) and value > 0:
            return value
    raise FileNotFoundError(f"No valid Step 4 FWHM metadata for {filename}")


def _load_matches(spec: ClusterSpec) -> tuple[pd.DataFrame, dict[str, float]]:
    psf_path = spec.psf_result / "cmd_psf" / f"photometry_{spec.filename}.tsv"
    aperture_path = (
        spec.aperture_result
        / "step7_forced_phot"
        / f"photometry_{spec.filename}.tsv"
    )
    if not psf_path.exists():
        raise FileNotFoundError(psf_path)
    if not aperture_path.exists():
        raise FileNotFoundError(aperture_path)

    psf = pd.read_csv(psf_path, sep="\t")
    aperture = pd.read_csv(aperture_path, sep="\t", low_memory=False)
    _numeric(
        psf,
        (
            "x_fit",
            "y_fit",
            "mag_psf",
            "mag_psf_err",
            "snr_psf",
            "qfit_noise_ratio",
            "reduced_chi2",
            "flags_psf",
        ),
    )
    _numeric(
        aperture,
        ("x_fit", "y_fit", "mag_inst", "mag_err", "snr", "master_id"),
    )

    psf_good = (
        np.isfinite(psf["x_fit"])
        & np.isfinite(psf["y_fit"])
        & np.isfinite(psf["mag_psf"])
        & np.isfinite(psf["mag_psf_err"])
        & (psf["mag_psf_err"] > 0)
        & np.isfinite(psf["snr_psf"])
        & (psf["snr_psf"] > 0)
        & (psf["flags_psf"].fillna(0) == 0)
    )
    aperture_good = (
        np.isfinite(aperture["x_fit"])
        & np.isfinite(aperture["y_fit"])
        & np.isfinite(aperture["mag_inst"])
        & np.isfinite(aperture["mag_err"])
        & (aperture["mag_err"] > 0)
        & np.isfinite(aperture["snr"])
        & (aperture["snr"] > 0)
    )
    if "bad_phot_flag" in aperture:
        aperture_good &= ~_as_bool(aperture["bad_phot_flag"])
    if "off_frame_flag" in aperture:
        aperture_good &= ~_as_bool(aperture["off_frame_flag"])

    psf = psf.loc[psf_good].copy().reset_index(drop=True)
    aperture = aperture.loc[aperture_good].copy().reset_index(drop=True)
    if psf.empty or aperture.empty:
        raise RuntimeError(f"No valid photometry rows for {spec.title}")

    aperture_xy = aperture[["x_fit", "y_fit"]].to_numpy(dtype=float)
    tree = cKDTree(aperture_xy)
    distance, aperture_index = tree.query(
        psf[["x_fit", "y_fit"]].to_numpy(dtype=float), k=1, workers=1
    )
    pairs = pd.DataFrame(
        {
            "psf_index": np.arange(len(psf), dtype=int),
            "aperture_index": aperture_index.astype(int),
            "match_distance_px": distance,
        }
    )
    pairs = pairs[pairs["match_distance_px"] <= 1.0]
    pairs = pairs.sort_values("match_distance_px").drop_duplicates(
        "aperture_index"
    )
    p = psf.iloc[pairs["psf_index"].to_numpy(dtype=int)].reset_index(drop=True)
    a = aperture.iloc[pairs["aperture_index"].to_numpy(dtype=int)].reset_index(
        drop=True
    )

    fwhm_px = _read_fwhm(spec.psf_result, spec.filename)
    if len(aperture_xy) >= 2:
        nearest_px = tree.query(aperture_xy, k=2, workers=1)[0][:, 1]
        matched_neighbor = nearest_px[pairs["aperture_index"].to_numpy(dtype=int)]
    else:
        matched_neighbor = np.full(len(pairs), np.nan)

    matched = pd.DataFrame(
        {
            "cluster": spec.slug,
            "cluster_type": spec.cluster_type,
            "file": spec.filename,
            "x": p["x_fit"].to_numpy(dtype=float),
            "y": p["y_fit"].to_numpy(dtype=float),
            "master_id": pd.to_numeric(
                a.get("master_id", pd.Series(np.arange(len(a)))), errors="coerce"
            ).to_numpy(),
            "mag_psf": p["mag_psf"].to_numpy(dtype=float),
            "mag_aperture": a["mag_inst"].to_numpy(dtype=float),
            "err_psf": p["mag_psf_err"].to_numpy(dtype=float),
            "err_aperture": a["mag_err"].to_numpy(dtype=float),
            "snr_psf": p["snr_psf"].to_numpy(dtype=float),
            "snr_aperture": a["snr"].to_numpy(dtype=float),
            "nearest_neighbor_px": matched_neighbor,
            "nearest_neighbor_fwhm": matched_neighbor / fwhm_px,
            "match_distance_px": pairs["match_distance_px"].to_numpy(dtype=float),
        }
    )
    for optional in ("qfit_noise_ratio", "reduced_chi2"):
        if optional in p:
            matched[optional] = p[optional].to_numpy(dtype=float)

    matched["snr_min"] = np.minimum(matched["snr_psf"], matched["snr_aperture"])
    matched["err_combined"] = np.hypot(
        matched["err_psf"], matched["err_aperture"]
    )
    matched["delta_raw"] = matched["mag_psf"] - matched["mag_aperture"]

    reference = matched[
        (matched["snr_min"] >= 50.0)
        & (matched["nearest_neighbor_fwhm"] >= 6.0)
    ]
    if len(reference) < 20:
        reference = matched[
            (matched["snr_min"] >= 20.0)
            & (matched["nearest_neighbor_fwhm"] >= 3.0)
        ]
    if len(reference) < 10:
        reference = matched[matched["snr_min"] >= 20.0]
    if reference.empty:
        reference = matched
    offset = float(np.median(reference["delta_raw"]))
    matched["delta_centered"] = matched["delta_raw"] - offset
    matched["z_combined"] = matched["delta_centered"] / matched["err_combined"]
    metadata = {
        "fwhm_px": fwhm_px,
        "offset_mag": offset,
        "n_offset_reference": int(len(reference)),
        "n_psf_valid": int(len(psf)),
        "n_aperture_valid": int(len(aperture)),
        "n_matched": int(len(matched)),
    }
    return matched, metadata


def _binned_summary(data: pd.DataFrame, column: str, edges, labels) -> pd.DataFrame:
    work = data.copy()
    work["bin"] = pd.cut(
        work[column], bins=edges, labels=labels, right=False, include_lowest=True
    )
    rows = []
    for label in labels:
        subset = work[work["bin"] == label]
        if subset.empty:
            continue
        scatter = _robust_scatter(subset["delta_centered"])
        reported = float(np.median(subset["err_combined"]))
        rows.append(
            {
                "bin": str(label),
                "n": int(len(subset)),
                "median_x": float(np.median(subset[column])),
                "median_delta_mag": float(np.median(subset["delta_centered"])),
                "robust_scatter_mag": scatter,
                "median_err_psf_mag": float(np.median(subset["err_psf"])),
                "median_err_aperture_mag": float(
                    np.median(subset["err_aperture"])
                ),
                "median_err_combined_mag": reported,
                "scatter_to_reported": scatter / reported if reported > 0 else np.nan,
                "coverage_1sigma": float(
                    np.mean(np.abs(subset["delta_centered"]) <= subset["err_combined"])
                ),
                "robust_z_scatter": _robust_scatter(subset["z_combined"]),
            }
        )
    return pd.DataFrame(rows)


def _method_error_curve(
    data: pd.DataFrame,
    snr_column: str,
    error_column: str,
) -> tuple[np.ndarray, np.ndarray]:
    work = data[[snr_column, error_column]].replace([np.inf, -np.inf], np.nan)
    work = work.dropna()
    work = work[(work[snr_column] >= SNR_EDGES[0]) & (work[error_column] > 0)]
    categories = pd.cut(
        work[snr_column],
        bins=SNR_EDGES,
        labels=SNR_LABELS,
        right=False,
        include_lowest=True,
    )
    x_values, y_values = [], []
    for label in SNR_LABELS:
        subset = work[categories == label]
        if subset.empty:
            continue
        x_values.append(float(np.median(subset[snr_column])))
        y_values.append(float(np.median(subset[error_column])))
    return np.asarray(x_values), np.asarray(y_values)


def _sample_summary(
    spec: ClusterSpec,
    matched: pd.DataFrame,
    metadata: dict[str, float],
    snr_bins: pd.DataFrame,
) -> dict[str, float | int | str]:
    science = matched[matched["snr_min"] >= 5.0]
    high_isolated = science[
        (science["snr_min"] >= 50.0)
        & (science["nearest_neighbor_fwhm"] >= 6.0)
    ]
    scatter = _robust_scatter(science["delta_centered"])
    reported = float(np.median(science["err_combined"]))
    return {
        "cluster": spec.slug,
        "title": spec.title,
        "cluster_type": spec.cluster_type,
        "file": spec.filename,
        **metadata,
        "n_snr_ge_5": int(len(science)),
        "median_delta_mag": float(np.median(science["delta_centered"])),
        "robust_scatter_mag": scatter,
        "median_err_psf_mag": float(np.median(science["err_psf"])),
        "median_err_aperture_mag": float(np.median(science["err_aperture"])),
        "median_err_combined_mag": reported,
        "scatter_to_reported": scatter / reported if reported > 0 else np.nan,
        "coverage_1sigma": float(
            np.mean(np.abs(science["delta_centered"]) <= science["err_combined"])
        ),
        "n_high_snr_isolated": int(len(high_isolated)),
        "high_snr_isolated_centered_bias_mag": (
            float(np.median(high_isolated["delta_centered"]))
            if len(high_isolated)
            else np.nan
        ),
        "high_snr_isolated_scatter_mag": _robust_scatter(
            high_isolated["delta_centered"]
        ),
        "snr_bin_count": int(len(snr_bins)),
    }


def _plot_cluster(
    spec: ClusterSpec,
    matched: pd.DataFrame,
    snr_bins: pd.DataFrame,
    crowd_bins: pd.DataFrame,
    summary: dict[str, float | int | str],
    output_dir: Path,
) -> tuple[Path, Path]:
    data = matched[matched["snr_min"] >= 5.0].copy()
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 8.2))
    colors = {"psf": "#0072B2", "aperture": "#D55E00", "combined": "#009E73"}

    ax = axes[0, 0]
    corrected_psf = data["mag_psf"] - float(summary["offset_mag"])
    ax.hexbin(
        data["mag_aperture"], corrected_psf, gridsize=44, mincnt=1,
        bins="log", cmap="Blues", linewidths=0
    )
    low = float(np.nanpercentile(np.r_[data["mag_aperture"], corrected_psf], 1))
    high = float(np.nanpercentile(np.r_[data["mag_aperture"], corrected_psf], 99))
    ax.plot([low, high], [low, high], color="0.25", ls="--", lw=1)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    ax.set_xlabel("aperture instrumental magnitude")
    ax.set_ylabel("PSF magnitude - bright-star offset")
    ax.set_title("(a) same-frame agreement", loc="left")

    ax = axes[0, 1]
    scatter = ax.scatter(
        data["snr_min"], data["delta_centered"],
        c=np.clip(data["nearest_neighbor_fwhm"], 0, 8), s=9, alpha=0.42,
        cmap="viridis", vmin=0, vmax=8, edgecolors="none"
    )
    ax.axhline(0, color="0.3", ls="--", lw=0.9)
    if not snr_bins.empty:
        x = snr_bins["median_x"].to_numpy(float)
        y = snr_bins["median_delta_mag"].to_numpy(float)
        e = snr_bins["robust_scatter_mag"].to_numpy(float)
        ax.errorbar(x, y, yerr=e, fmt="o-", color="#CC79A7", lw=1.4, capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel("minimum of aperture and PSF SNR")
    ax.set_ylabel("centered PSF - aperture (mag)")
    ax.set_title("(b) residual and crowding", loc="left")
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.02)
    colorbar.set_label("nearest neighbor (FWHM)")

    ax = axes[0, 2]
    show = data[(data["err_psf"] > 0) & (data["err_aperture"] > 0)]
    take = np.linspace(0, max(len(show) - 1, 0), min(len(show), 800), dtype=int)
    sampled = show.iloc[take] if len(show) else show
    ax.scatter(sampled["snr_psf"], sampled["err_psf"], s=4, alpha=0.08, color=colors["psf"])
    ax.scatter(sampled["snr_aperture"], sampled["err_aperture"], s=4, alpha=0.08, color=colors["aperture"])
    psf_x, psf_y = _method_error_curve(data, "snr_psf", "err_psf")
    aperture_x, aperture_y = _method_error_curve(
        data, "snr_aperture", "err_aperture"
    )
    ax.plot(psf_x, psf_y, "o-", color=colors["psf"], label="PSF reported")
    ax.plot(aperture_x, aperture_y, "s-", color=colors["aperture"], label="aperture reported")
    theory_x = np.geomspace(
        max(5.0, float(data[["snr_psf", "snr_aperture"]].min().min())),
        float(data[["snr_psf", "snr_aperture"]].max().max()),
        100,
    )
    ax.plot(theory_x, 1.0857 / theory_x, "--", color="0.25", label="1.0857 / SNR")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("minimum SNR")
    ax.set_ylabel("reported magnitude error")
    ax.set_title("(c) individual error models", loc="left")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 0]
    positions = np.arange(len(snr_bins))
    width = 0.36
    ax.bar(
        positions - width / 2, snr_bins["robust_scatter_mag"], width,
        color="#56B4E9", label="observed robust scatter"
    )
    ax.bar(
        positions + width / 2, snr_bins["median_err_combined_mag"], width,
        color=colors["combined"], label="reported combined error"
    )
    for pos, row in enumerate(snr_bins.itertuples(index=False)):
        ax.text(
            pos, max(row.robust_scatter_mag, row.median_err_combined_mag) * 1.08,
            f"{row.scatter_to_reported:.1f}x\nN={row.n}", ha="center", va="bottom", fontsize=8
        )
    positive_scales = np.r_[
        snr_bins["robust_scatter_mag"].to_numpy(float),
        snr_bins["median_err_combined_mag"].to_numpy(float),
    ]
    positive_scales = positive_scales[np.isfinite(positive_scales) & (positive_scales > 0)]
    if positive_scales.size:
        ax.set_ylim(float(positive_scales.min()) * 0.70, float(positive_scales.max()) * 1.9)
    ax.set_yscale("log")
    ax.set_xticks(positions, snr_bins["bin"])
    ax.set_xlabel("SNR bin")
    ax.set_ylabel("magnitude scale")
    ax.set_title("(d) observed versus reported", loc="left")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 1]
    z_all = data.loc[data["snr_min"] >= 10, "z_combined"].to_numpy(float)
    z_all = z_all[np.isfinite(z_all) & (np.abs(z_all) <= 8)]
    z_iso = data.loc[
        (data["snr_min"] >= 50) & (data["nearest_neighbor_fwhm"] >= 6),
        "z_combined",
    ].to_numpy(float)
    z_iso = z_iso[np.isfinite(z_iso) & (np.abs(z_iso) <= 8)]
    bins = np.linspace(-6, 6, 49)
    if len(z_all):
        ax.hist(z_all, bins=bins, density=True, histtype="step", lw=1.6, color="#0072B2", label=f"SNR>=10, N={len(z_all)}")
    if len(z_iso):
        ax.hist(z_iso, bins=bins, density=True, histtype="step", lw=1.6, color="#D55E00", label=f"SNR>=50, >=6 FWHM, N={len(z_iso)}")
    xx = np.linspace(-6, 6, 400)
    ax.plot(xx, norm.pdf(xx), "--", color="0.25", label="N(0,1)")
    ax.set_xlabel("(centered PSF - aperture) / combined error")
    ax.set_ylabel("density")
    ax.set_title("(e) normalized residual", loc="left")
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1, 2]
    positions = np.arange(len(crowd_bins))
    ax.errorbar(
        positions, crowd_bins["median_delta_mag"],
        yerr=crowd_bins["robust_scatter_mag"], fmt="o-", color="#0072B2",
        ecolor="#56B4E9", capsize=4, label="median +/- robust scatter"
    )
    for pos, row in enumerate(crowd_bins.itertuples(index=False)):
        ax.text(
            pos, row.median_delta_mag + row.robust_scatter_mag * 1.08,
            f"{row.scatter_to_reported:.1f}x\nN={row.n}", ha="center", va="bottom", fontsize=8
        )
    ax.axhline(0, color="0.3", ls="--", lw=0.9)
    ax.set_xticks(positions, crowd_bins["bin"])
    ax.set_xlabel("nearest-neighbor separation (FWHM)")
    ax.set_ylabel("centered PSF - aperture (mag)")
    ax.set_title("(f) crowding and error ratio", loc="left")
    lower = crowd_bins["median_delta_mag"] - crowd_bins["robust_scatter_mag"]
    upper = crowd_bins["median_delta_mag"] + crowd_bins["robust_scatter_mag"]
    if len(crowd_bins):
        span = max(float(upper.max() - lower.min()), 0.02)
        ax.set_ylim(float(lower.min()) - 0.18 * span, float(upper.max()) + 0.42 * span)

    for ax in axes.ravel():
        ax.grid(alpha=0.18)
    fig.suptitle(
        f"{spec.title} ({spec.cluster_type}): aperture vs APEX PSF error diagnostic\n"
        f"N={summary['n_snr_ge_5']}, offset={summary['offset_mag']:+.3f} mag, "
        f"scatter/reported={summary['scatter_to_reported']:.2f}; neither method is treated as truth",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.008,
        "Reported combined error is the quadrature sum; shared pixels make the independence assumption approximate.",
        ha="center",
        va="bottom",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.025, 1, 0.94))
    output_dir.mkdir(parents=True, exist_ok=True)
    png_path = output_dir / f"{spec.slug}_aperture_vs_psf_error.png"
    pdf_path = output_dir / f"{spec.slug}_aperture_vs_psf_error.pdf"
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    return png_path, pdf_path


def _default_specs(args) -> list[ClusterSpec]:
    return [
        ClusterSpec(
            slug="m13",
            title="M13",
            cluster_type="globular cluster",
            filename="pp_messier13-0003-B.fit",
            psf_result=args.m13_psf_result,
            aperture_result=args.m13_aperture_result,
        ),
        ClusterSpec(
            slug="ngc6811",
            title="NGC 6811",
            cluster_type="open cluster",
            filename="pp_NGC6811-0006-B.fit",
            psf_result=args.ngc6811_psf_result,
            aperture_result=args.ngc6811_aperture_result,
        ),
    ]


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--m13-psf-result",
        type=Path,
        default=REPO / "validation" / "real_gui_run" / "psf_auto_v3" / "m13_local_radius_q3" / "result",
    )
    parser.add_argument(
        "--m13-aperture-result",
        type=Path,
        default=Path(r"E:\APEX_validation\reprocess\M13\result"),
    )
    parser.add_argument(
        "--ngc6811-psf-result",
        type=Path,
        default=REPO / "validation" / "real_gui_run" / "psf_auto_v3" / "ngc6811_current_qc" / "result",
    )
    parser.add_argument(
        "--ngc6811-aperture-result",
        type=Path,
        default=Path(r"E:\APEX_validation\reprocess\NGC6811\result"),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries = []
    manifest = []
    for spec in _default_specs(args):
        matched, metadata = _load_matches(spec)
        science = matched[matched["snr_min"] >= 5.0]
        snr_bins = _binned_summary(science, "snr_min", SNR_EDGES, SNR_LABELS)
        crowd_bins = _binned_summary(
            science, "nearest_neighbor_fwhm", CROWD_EDGES, CROWD_LABELS
        )
        summary = _sample_summary(spec, matched, metadata, snr_bins)
        png_path, pdf_path = _plot_cluster(
            spec, matched, snr_bins, crowd_bins, summary, args.output_dir
        )
        matched.to_csv(args.output_dir / f"{spec.slug}_matched.csv", index=False)
        snr_bins.to_csv(args.output_dir / f"{spec.slug}_snr_bins.csv", index=False)
        crowd_bins.to_csv(args.output_dir / f"{spec.slug}_crowding_bins.csv", index=False)
        summaries.append(summary)
        manifest.append(
            {
                "cluster": spec.slug,
                "psf_result": str(spec.psf_result),
                "aperture_result": str(spec.aperture_result),
                "file": spec.filename,
                "figure_png": str(png_path),
                "figure_pdf": str(pdf_path),
            }
        )

    summary_frame = pd.DataFrame(summaries)
    summary_frame.to_csv(args.output_dir / "aperture_psf_error_summary.csv", index=False)
    (args.output_dir / "aperture_psf_error_manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    print(summary_frame.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
