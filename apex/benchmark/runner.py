"""End-to-end artificial-star benchmark runner."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
from typing import Any

from astropy.io import fits
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import tomllib
except ImportError:  # pragma: no cover - Python 3.10
    import tomli as tomllib

from apex.benchmark.artificial_stars import (
    extract_empirical_psf,
    inject_catalog,
    sample_injection_catalog,
)
from apex.benchmark.metrics import (
    count_new_false_detections,
    fit_completeness_logistic,
    magnitude_bin_summary,
    magnitude_point_summary,
    mark_baseline_confounding,
    match_truth_to_detections,
    summarize_results,
    summarize_by_placement,
)
from apex.config.parameters_cmd import Parameters
from apex.utils.noise_params import resolve_effective_noise_params
from apex.utils.photometry_utils import phot_vectorized


@dataclass
class BenchmarkConfig:
    input_fits: str = ""
    parameter_file: str = "parameters.toml"
    output_root: str = "benchmark/runs"
    seed: int = 20260607
    trials: int = 3
    stars_per_trial: int = 40
    magnitude_min: float = 15.0
    magnitude_max: float = 22.0
    magnitude_grid: list[float] = field(default_factory=list)
    stars_per_magnitude_per_trial: int = 0
    magnitude_bin_width: float = 1.0
    completeness_bootstrap_samples: int = 1000
    save_injected_fits: bool = True
    zeropoint_mag: float | None = None
    allow_initial_zeropoint_fallback: bool = False
    match_radius_px: float = 2.0
    baseline_confound_radius_px: float = 2.0
    isolated_fraction: float = 0.5
    isolated_min_real_sep_fwhm: float = 3.0
    min_injected_sep_fwhm: float = 5.0
    psf_radius_fwhm: float = 4.0
    psf_isolation_fwhm: float = 4.0
    psf_max_stars: int = 40
    psf_min_stars: int = 3
    psf_allow_nonisolated_fallback: bool = False
    aperture_scale_fwhm: float | None = None
    annulus_scale_fwhm: float | None = None
    annulus_width_fwhm: float | None = None
    detection_overrides: dict[str, Any] = field(default_factory=dict)


def load_benchmark_config(path: str | Path) -> BenchmarkConfig:
    config_path = Path(path)
    with config_path.open("rb") as handle:
        raw = tomllib.load(handle)
    section = dict(raw.get("benchmark", raw))
    detection_overrides = dict(raw.get("detection_overrides", {}))
    known = set(BenchmarkConfig.__dataclass_fields__)
    unknown = sorted(set(section) - known)
    if unknown:
        raise ValueError(f"Unknown benchmark config keys: {', '.join(unknown)}")
    section["detection_overrides"] = detection_overrides
    return BenchmarkConfig(**section)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_dump(path: Path, payload: Any) -> None:
    def clean(value: Any) -> Any:
        if isinstance(value, dict):
            return {str(key): clean(item) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [clean(item) for item in value]
        if isinstance(value, np.ndarray):
            return clean(value.tolist())
        if isinstance(value, (np.integer,)):
            return int(value)
        if isinstance(value, (np.floating, float)):
            number = float(value)
            return number if np.isfinite(number) else None
        if isinstance(value, Path):
            return str(value)
        return value

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(clean(payload), indent=2, sort_keys=True, allow_nan=False, default=str),
        encoding="utf-8",
    )


class _BenchmarkParameters:
    def __init__(self, base: Parameters, input_path: Path, result_dir: Path, overrides: dict):
        import copy

        self.P = copy.deepcopy(base.P)
        self.param_file = base.param_file
        self.param_hash = base.param_hash
        self.P.data_dir = input_path.parent
        self.P.result_dir = result_dir
        self.P.cache_dir = result_dir / "cache"
        self.P.file_path_map = {input_path.name: str(input_path)}
        self.P.max_workers = 1
        self.P.parallel_max_workers = 1
        for key, value in overrides.items():
            if not hasattr(self.P, key):
                raise ValueError(f"Unknown detection override: {key}")
            setattr(self.P, key, value)

    def get_file_path(self, filename: str) -> Path:
        mapped = self.P.file_path_map.get(filename)
        return Path(mapped) if mapped else Path(self.P.data_dir) / filename


def run_apex_detection(
    input_path: Path,
    base_params: Parameters,
    output_dir: Path,
    detection_overrides: dict[str, Any],
) -> tuple[pd.DataFrame, dict, float]:
    """Run the production Step 4 worker synchronously for one FITS file."""
    from apex.gui.workflow.step4_source_detection import DetectionWorker

    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    params = _BenchmarkParameters(base_params, input_path, output_dir, detection_overrides)
    p = params.P
    # Per-filter sigma may be present-but-None (e.g. empty sigma_by_filter), in
    # which case getattr returns None rather than the default — fall back to the
    # base detect_sigma, then to 3.2. (A sigma is always > 0, so `or` is safe.)
    sigma_map = {
        key: float(
            getattr(p, f"detect_sigma_{key}", None)
            or getattr(p, "detect_sigma", None)
            or 3.2
        )
        for key in ("g", "r", "i")
    }
    worker = DetectionWorker(
        [input_path.name],
        params,
        input_path.parent,
        cache_dir,
        use_cropped=False,
        filter_sigma_map=sigma_map,
    )
    started = time.perf_counter()
    worker.run()
    elapsed = time.perf_counter() - started

    step_dir = output_dir / "step4_detection"
    csv_path = step_dir / f"detect_{input_path.name}.csv"
    json_path = step_dir / f"detect_{input_path.name}.json"
    if not csv_path.exists() or not json_path.exists():
        raise RuntimeError(f"Step 4 did not produce detection outputs for {input_path.name}")
    detections = pd.read_csv(csv_path)
    metadata = json.loads(json_path.read_text(encoding="utf-8"))
    return detections, metadata, elapsed


def _photometry_radii(config: BenchmarkConfig, p: Any, fwhm_px: float) -> tuple[float, float, float]:
    aperture_scale = (
        float(config.aperture_scale_fwhm)
        if config.aperture_scale_fwhm is not None
        else float(getattr(p, "forced_r_ap_scale", 0.8))
    )
    annulus_scale = (
        float(config.annulus_scale_fwhm)
        if config.annulus_scale_fwhm is not None
        else float(getattr(p, "fitsky_annulus_scale", 4.0))
    )
    annulus_width = (
        float(config.annulus_width_fwhm)
        if config.annulus_width_fwhm is not None
        else float(getattr(p, "fitsky_dannulus_scale", 2.0))
    )
    r_ap = max(float(getattr(p, "min_r_ap_px", 4.0)), aperture_scale * fwhm_px)
    ann_gap = float(getattr(p, "annulus_min_gap_px", 6.0))
    r_in = max(r_ap + ann_gap, annulus_scale * fwhm_px)
    r_out = r_in + max(ann_gap, annulus_width * fwhm_px)
    return r_ap, r_in, r_out


def _measure_positions(
    image: np.ndarray,
    positions: np.ndarray,
    header: fits.Header,
    p: Any,
    radii: tuple[float, float, float],
) -> np.ndarray:
    noise = resolve_effective_noise_params(p, header)
    r_ap, r_in, r_out = radii
    flux_e, *_ = phot_vectorized(
        image,
        positions,
        r_ap,
        r_in,
        r_out,
        gain=noise.gain_e_per_adu,
        rn_param_e=noise.rdnoise_e,
        sky_frame_e=np.nan,
        sky_sigma_mode=str(getattr(p, "sky_sigma_mode", "local")),
        sky_sigma_includes_rn=bool(getattr(p, "sky_sigma_includes_rn", True)),
        min_n_sky_for_local=int(getattr(p, "sky_sigma_min_n_sky", 50)),
        sat_adu=float(getattr(p, "saturation_adu", 60000.0)),
        datamax_adu=float(getattr(p, "datamax_adu", 55000.0)),
        sigma_clip_val=float(getattr(p, "phot_sigma_clip", 3.0)),
        maxiters=int(getattr(p, "phot_max_iter", 5)),
    )
    return np.asarray(flux_e, dtype=float)


def _version_manifest() -> dict[str, str]:
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for package in ("numpy", "pandas", "scipy", "astropy", "photutils", "sep", "PyQt5"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def resolve_benchmark_zeropoint(
    input_path: Path,
    header: fits.Header,
    p: Any,
    explicit_zeropoint: float | None,
    *,
    allow_initial_fallback: bool,
) -> tuple[float, str, dict[str, Any]]:
    """Resolve a zero point matching the total-electron Step 7 convention."""
    if explicit_zeropoint is not None:
        value = float(explicit_zeropoint)
        if not np.isfinite(value):
            raise ValueError("zeropoint_mag must be finite")
        return value, "config:zeropoint_mag", {"color_term_assumed": 0.0}

    result_dir = Path(getattr(p, "result_dir", ""))
    filter_name = str(header.get("FILTER", "")).strip().lower()
    frame_path = result_dir / "cmd_zeropoint" / "frame_zeropoint.csv"
    if frame_path.exists():
        frame_df = pd.read_csv(frame_path)
        if {"file", "zp_frame"} <= set(frame_df.columns):
            rows = frame_df[frame_df["file"].astype(str) == input_path.name].copy()
            if filter_name and "filter" in rows.columns:
                filtered = rows[
                    rows["filter"].astype(str).str.strip().str.lower() == filter_name
                ]
                if not filtered.empty:
                    rows = filtered
            values = pd.to_numeric(rows["zp_frame"], errors="coerce")
            values = values[np.isfinite(values)]
            if not values.empty:
                row = rows.loc[values.index[0]]
                return (
                    float(values.iloc[0]),
                    f"step10:{frame_path.name}:frame",
                    {
                        "frame_zeropoint_path": str(frame_path),
                        "zp_scatter": row.get("zp_scatter"),
                        "n_ref": row.get("n_ref"),
                        "color_term_assumed": 0.0,
                    },
                )

    coeff_path = result_dir / "cmd_zeropoint" / "zp_fit_coefficients.csv"
    if coeff_path.exists() and filter_name:
        coeff_df = pd.read_csv(coeff_path)
        if {"filter", "zp"} <= set(coeff_df.columns):
            rows = coeff_df[
                coeff_df["filter"].astype(str).str.strip().str.lower() == filter_name
            ]
            values = pd.to_numeric(rows["zp"], errors="coerce")
            values = values[np.isfinite(values)]
            if not values.empty:
                row = rows.loc[values.index[0]]
                return (
                    float(values.iloc[0]),
                    f"step10:{coeff_path.name}:filter",
                    {
                        "zeropoint_coefficients_path": str(coeff_path),
                        "zp_scatter": row.get("scatter_rms"),
                        "n_ref": row.get("N"),
                        "color_term_assumed": 0.0,
                    },
                )

    if allow_initial_fallback:
        value = float(getattr(p, "zp_initial", 25.0))
        return value, "parameters:zp_initial", {"color_term_assumed": 0.0}

    raise RuntimeError(
        f"No calibrated Step 10 zero point found for {input_path.name}. "
        "Run CMD Step 10, set benchmark.zeropoint_mag explicitly, or enable "
        "allow_initial_zeropoint_fallback for engineering-only tests."
    )


def _git_manifest(repo_root: Path) -> dict[str, Any]:
    def run(*args: str) -> str:
        try:
            return subprocess.check_output(
                ["git", *args],
                cwd=repo_root,
                text=True,
                stderr=subprocess.DEVNULL,
                timeout=5,
            ).strip()
        except Exception:
            return ""

    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
        "dirty": bool(run("status", "--porcelain")),
    }


def _write_plots(stars: pd.DataFrame, bins: pd.DataFrame, report_dir: Path) -> None:
    report_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    if len(bins):
        axes[0].plot(bins["magnitude_center"], bins["completeness"], marker="o")
    axes[0].set(xlabel="Injected magnitude", ylabel="Completeness", ylim=(-0.03, 1.03))
    axes[0].grid(alpha=0.25)

    eligible = ~stars["baseline_confounded"].astype(bool)
    mag_error = pd.to_numeric(stars["forced_mag_error"], errors="coerce")
    good = eligible & np.isfinite(mag_error)
    axes[1].scatter(
        stars.loc[good, "magnitude_true"],
        mag_error[good],
        s=9,
        alpha=0.55,
    )
    axes[1].axhline(0.0, color="black", linewidth=1)
    axes[1].set(xlabel="Injected magnitude", ylabel="Recovered - expected mag")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(report_dir / "benchmark_summary.png", dpi=160)
    plt.close(fig)


def _write_precision_plot(
    points: pd.DataFrame,
    fit: dict,
    report_dir: Path,
) -> None:
    if points.empty or not fit:
        return
    x = points["magnitude"].to_numpy(float)
    y = points["completeness"].to_numpy(float)
    lower = np.maximum(0.0, y - points["ci95_low"].to_numpy(float))
    upper = np.maximum(0.0, points["ci95_high"].to_numpy(float) - y)
    grid = np.linspace(float(x.min()) - 0.1, float(x.max()) + 0.1, 400)
    model = 1.0 / (
        1.0 + np.exp((grid - float(fit["m50"])) / float(fit["width_mag"]))
    )
    fig, ax = plt.subplots(figsize=(7.5, 5.0))
    ax.errorbar(x, y, yerr=np.vstack([lower, upper]), fmt="o", capsize=3, label="Injection points")
    ax.plot(grid, model, linewidth=2, label="Binomial logistic fit")
    ax.axhline(0.5, color="black", linewidth=1, linestyle="--", alpha=0.7)
    ax.axvline(float(fit["m50"]), color="#c43d3d", linewidth=1.2, linestyle="--")
    ax.set(
        xlabel="Injected calibrated magnitude",
        ylabel="Completeness",
        ylim=(-0.04, 1.04),
    )
    ax.grid(alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(report_dir / "completeness_logistic_fit.png", dpi=180)
    plt.close(fig)


def _write_markdown_report(
    path: Path,
    summary: dict,
    config: BenchmarkConfig,
    fwhm_px: float,
    psf_meta: dict,
) -> None:
    fit = summary.get("completeness_fit") or {}
    fit_lines = ""
    if fit:
        fit_lines = (
            f"- 90% completeness magnitude: {fit['m90']:.4f} mag\n"
            f"- 50% completeness magnitude: {fit['m50']:.4f} mag "
            f"(95% CI {fit.get('m50_ci95_low', float('nan')):.4f} to "
            f"{fit.get('m50_ci95_high', float('nan')):.4f})\n"
            f"- 10% completeness magnitude: {fit['m10']:.4f} mag\n"
        )
    if config.magnitude_grid:
        magnitude_setup = (
            f"- Magnitude grid: {', '.join(f'{value:.3f}' for value in config.magnitude_grid)}\n"
            f"- Stars per magnitude per trial: {config.stars_per_magnitude_per_trial}"
        )
        stars_per_trial = len(config.magnitude_grid) * config.stars_per_magnitude_per_trial
    else:
        magnitude_setup = f"- Magnitude range: {config.magnitude_min} to {config.magnitude_max}"
        stars_per_trial = config.stars_per_trial
    text = f"""# APEX Artificial-Star Benchmark

Generated: {datetime.now(timezone.utc).isoformat()}

## Scope

This run injects an empirical PSF into an existing FITS observation, reruns the
production APEX Step 4 detector, and measures differential aperture flux with
the same `phot_vectorized` function used by Step 7. It estimates selection and
measurement behavior for this image and configuration; it is not an absolute
calibration test.

## Summary

- Injected stars: {summary['n_injected']}
- Completeness-eligible stars: {summary['n_completeness_eligible']}
- Baseline-confounded stars: {summary['n_baseline_confounded']}
- Recovered eligible stars: {summary['n_recovered']}
- Completeness: {summary['completeness']:.4f}
- Position RMSE: {summary['position_rmse_px']:.4f} px
- Forced magnitude bias: {summary['forced_mag_bias_median']:.4f} mag
- Forced magnitude MAD scatter: {summary['forced_mag_scatter_mad']:.4f} mag
- New detections: {summary['new_detections']}
- New unmatched detections: {summary['new_false_detections']}
{fit_lines}

## Effective Setup

- Trials: {config.trials}
- Stars per trial: {stars_per_trial}
{magnitude_setup}
- Zeropoint: {summary['zeropoint_mag']:.6f}
- Zeropoint source: {summary['zeropoint_source']}
- Magnitude definition: total source electrons, zero injected color term
- FWHM: {fwhm_px:.3f} px
- Empirical PSF stars used: {psf_meta['n_used']}
- Empirical PSF size: {psf_meta['stamp_size']} px
- Random seed: {config.seed}
"""
    path.write_text(text, encoding="utf-8")


def run_benchmark(
    config: BenchmarkConfig,
    *,
    input_override: str | Path | None = None,
    output_override: str | Path | None = None,
) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    input_path = Path(input_override or config.input_fits).expanduser().resolve()
    if not input_path.exists():
        raise FileNotFoundError(f"Input FITS does not exist: {input_path}")
    parameter_path = Path(config.parameter_file)
    if not parameter_path.is_absolute():
        parameter_path = (repo_root / parameter_path).resolve()
    if not parameter_path.exists():
        raise FileNotFoundError(f"Parameter file does not exist: {parameter_path}")

    if output_override:
        run_dir = Path(output_override).expanduser().resolve()
    else:
        output_root = Path(config.output_root)
        if not output_root.is_absolute():
            output_root = repo_root / output_root
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        run_dir = output_root / f"{input_path.stem}_{stamp}_seed{config.seed}"
    if run_dir.exists() and any(run_dir.iterdir()):
        raise FileExistsError(f"Benchmark output directory is not empty: {run_dir}")
    run_dir.mkdir(parents=True, exist_ok=True)

    base_params = Parameters(parameter_path)
    p = base_params.P
    with fits.open(input_path, memmap=False) as hdul:
        image = np.asarray(hdul[0].data, dtype=float)
        header = hdul[0].header.copy()
    if image.ndim != 2:
        raise ValueError(f"Expected a 2D FITS image, got shape {image.shape}")

    baseline_dir = run_dir / "baseline"
    baseline_detections, baseline_meta, baseline_seconds = run_apex_detection(
        input_path,
        base_params,
        baseline_dir,
        config.detection_overrides,
    )
    fwhm_px = float(baseline_meta.get("fwhm_px") or getattr(p, "fwhm_seed_px", 6.0))
    if not np.isfinite(fwhm_px) or fwhm_px <= 0:
        fwhm_px = float(getattr(p, "fwhm_seed_px", 6.0))
    noise = resolve_effective_noise_params(p, header)
    psf_result = extract_empirical_psf(
        image,
        baseline_detections,
        fwhm_px,
        saturation_adu=float(getattr(p, "saturation_adu", 60000.0)),
        radius_fwhm=config.psf_radius_fwhm,
        isolation_fwhm=config.psf_isolation_fwhm,
        max_stars=config.psf_max_stars,
        min_stars=config.psf_min_stars,
        allow_nonisolated_fallback=config.psf_allow_nonisolated_fallback,
    )
    fits.writeto(run_dir / "empirical_psf.fits", psf_result.data.astype(np.float32), overwrite=True)

    zeropoint, zeropoint_source, zeropoint_meta = resolve_benchmark_zeropoint(
        input_path,
        header,
        p,
        config.zeropoint_mag,
        allow_initial_fallback=config.allow_initial_zeropoint_fallback,
    )
    radii = _photometry_radii(config, p, fwhm_px)
    real_xy = baseline_detections[["x", "y"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    all_rows: list[pd.DataFrame] = []
    trial_rows = []
    magnitude_grid = np.asarray(config.magnitude_grid, dtype=float)
    if magnitude_grid.size:
        if not np.isfinite(magnitude_grid).all():
            raise ValueError("magnitude_grid must contain only finite values")
        if len(np.unique(magnitude_grid)) != len(magnitude_grid):
            raise ValueError("magnitude_grid must not contain duplicate values")
        if config.stars_per_magnitude_per_trial <= 0:
            raise ValueError(
                "stars_per_magnitude_per_trial must be positive when magnitude_grid is used"
            )
        stars_per_trial = int(
            len(magnitude_grid) * config.stars_per_magnitude_per_trial
        )
    else:
        if config.stars_per_trial <= 0:
            raise ValueError("stars_per_trial must be positive")
        stars_per_trial = int(config.stars_per_trial)

    for trial in range(config.trials):
        trial_seed = int(config.seed + trial)
        rng = np.random.default_rng(trial_seed)
        truth = sample_injection_catalog(
            image.shape,
            real_xy,
            count=stars_per_trial,
            magnitude_min=(
                float(magnitude_grid.min()) if magnitude_grid.size else config.magnitude_min
            ),
            magnitude_max=(
                float(magnitude_grid.max()) if magnitude_grid.size else config.magnitude_max
            ),
            fwhm_px=fwhm_px,
            psf_size=psf_result.stamp_size,
            rng=rng,
            isolated_fraction=config.isolated_fraction,
            isolated_min_real_sep_fwhm=config.isolated_min_real_sep_fwhm,
            min_injected_sep_fwhm=max(
                config.min_injected_sep_fwhm,
                radii[2] / fwhm_px + 0.5 * psf_result.stamp_size / fwhm_px,
            ),
            extra_margin_px=radii[2],
        )
        if magnitude_grid.size:
            assigned = np.repeat(
                magnitude_grid,
                int(config.stars_per_magnitude_per_trial),
            )
            rng.shuffle(assigned)
            truth["magnitude_true"] = assigned
        injected, expected_signal, _, truth = inject_catalog(
            image,
            psf_result.data,
            truth,
            gain_e_per_adu=noise.gain_e_per_adu,
            zeropoint_mag=zeropoint,
            rng=rng,
        )
        truth.insert(0, "trial", trial)
        truth["trial_seed"] = trial_seed

        trial_dir = run_dir / f"trial_{trial:03d}"
        trial_dir.mkdir(parents=True, exist_ok=True)
        injected_path = trial_dir / f"{input_path.stem}_injected_{trial:03d}.fits"
        trial_header = header.copy()
        trial_header["APEXBEN"] = (True, "APEX artificial-star benchmark")
        trial_header["BENCHTR"] = (trial, "Benchmark trial index")
        trial_header["BENCHSED"] = (trial_seed, "Benchmark random seed")
        trial_header["BENCHZP"] = (zeropoint, "Injected-star zero point")
        trial_header.add_history("Artificial stars injected with an empirical PSF")
        fits.writeto(injected_path, injected.astype(np.float32), trial_header, overwrite=True)
        truth.to_csv(trial_dir / "truth.csv", index=False)

        detections, detection_meta, detection_seconds = run_apex_detection(
            injected_path,
            base_params,
            trial_dir / "apex",
            config.detection_overrides,
        )
        detections.to_csv(trial_dir / "detections.csv", index=False)
        matched = match_truth_to_detections(
            truth,
            detections,
            radius_px=config.match_radius_px,
        )
        matched = mark_baseline_confounding(
            matched,
            baseline_detections,
            radius_px=config.baseline_confound_radius_px,
        )

        positions = truth[["x_true", "y_true"]].to_numpy(float)
        base_flux = _measure_positions(image, positions, header, p, radii)
        injected_flux = _measure_positions(injected, positions, trial_header, p, radii)
        expected_aperture_flux = _measure_positions(
            expected_signal, positions, fits.Header(), p, radii
        )
        delta_flux = injected_flux - base_flux
        matched["baseline_aperture_flux_e"] = base_flux
        matched["injected_aperture_flux_e"] = injected_flux
        matched["expected_aperture_flux_e"] = expected_aperture_flux
        matched["recovered_delta_flux_e"] = delta_flux
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = injected_flux / expected_aperture_flux
            differential_ratio = delta_flux / expected_aperture_flux
            matched["forced_flux_ratio"] = ratio
            matched["forced_mag_error"] = np.where(
                np.isfinite(ratio) & (ratio > 0),
                -2.5 * np.log10(ratio),
                np.nan,
            )
            matched["differential_flux_ratio"] = differential_ratio
            matched["differential_mag_error"] = np.where(
                np.isfinite(differential_ratio) & (differential_ratio > 0),
                -2.5 * np.log10(differential_ratio),
                np.nan,
            )

        new_count, false_count = count_new_false_detections(
            matched,
            detections,
            baseline_detections,
            baseline_radius_px=config.baseline_confound_radius_px,
        )
        matched["new_detections_trial"] = new_count
        matched["new_false_detections_trial"] = false_count
        matched.to_csv(trial_dir / "matched_stars.csv", index=False)
        if not config.save_injected_fits:
            injected_path.unlink(missing_ok=True)
        all_rows.append(matched)
        trial_rows.append(
            {
                "trial": trial,
                "seed": trial_seed,
                "n_injected": len(truth),
                "n_detected_total": len(detections),
                "new_detections": new_count,
                "new_false_detections": false_count,
                "detection_seconds": detection_seconds,
                "fwhm_px": detection_meta.get("fwhm_px"),
            }
        )

    stars = pd.concat(all_rows, ignore_index=True)
    stars.to_csv(run_dir / "stars.csv", index=False)
    trials_df = pd.DataFrame(trial_rows)
    trials_df.to_csv(run_dir / "trials.csv", index=False)
    bins = magnitude_bin_summary(stars, bin_width=float(config.magnitude_bin_width))
    bins.to_csv(run_dir / "magnitude_bins.csv", index=False)
    points = magnitude_point_summary(stars)
    points.to_csv(run_dir / "magnitude_points.csv", index=False)
    summary = summarize_results(stars)
    summary["by_placement"] = summarize_by_placement(stars)
    summary["new_detections"] = int(trials_df["new_detections"].sum())
    summary["new_false_detections"] = int(trials_df["new_false_detections"].sum())
    summary["baseline_detection_seconds"] = baseline_seconds
    summary["trial_detection_seconds_total"] = float(trials_df["detection_seconds"].sum())
    summary["fwhm_px"] = fwhm_px
    summary["zeropoint_mag"] = zeropoint
    summary["zeropoint_source"] = zeropoint_source
    summary["aperture_radii_px"] = {
        "aperture": radii[0],
        "annulus_inner": radii[1],
        "annulus_outer": radii[2],
    }
    completeness_fit = {}
    try:
        completeness_fit = fit_completeness_logistic(
            stars,
            bootstrap_samples=config.completeness_bootstrap_samples,
            seed=config.seed,
        )
        summary["completeness_fit"] = completeness_fit
        _json_dump(run_dir / "completeness_fit.json", completeness_fit)
    except RuntimeError as exc:
        summary["completeness_fit_error"] = str(exc)
    fit_by_placement = {}
    if "placement_stratum" in stars.columns:
        for placement, group in stars.groupby("placement_stratum", dropna=False):
            try:
                fit_by_placement[str(placement)] = fit_completeness_logistic(
                    group.reset_index(drop=True),
                    bootstrap_samples=config.completeness_bootstrap_samples,
                    seed=config.seed + len(fit_by_placement) + 1,
                )
            except RuntimeError as exc:
                fit_by_placement[str(placement)] = {"error": str(exc)}
    if fit_by_placement:
        summary["completeness_fit_by_placement"] = fit_by_placement
        _json_dump(run_dir / "completeness_fit_by_placement.json", fit_by_placement)
    _json_dump(run_dir / "summary.json", summary)

    psf_meta = {
        "n_candidates": psf_result.n_candidates,
        "n_used": psf_result.n_used,
        "stamp_size": psf_result.stamp_size,
        "selection_mode": psf_result.selection_mode,
    }
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_fits": str(input_path),
        "input_sha256": _sha256(input_path),
        "parameter_file": str(parameter_path),
        "parameter_sha256": _sha256(parameter_path),
        "effective_config": asdict(config),
        "versions": _version_manifest(),
        "git": _git_manifest(repo_root),
        "instrument": {
            "gain_e_per_adu": noise.gain_e_per_adu,
            "rdnoise_e": noise.rdnoise_e,
            "gain_source": noise.gain_source,
            "rdnoise_source": noise.rdnoise_source,
            "zeropoint_mag": zeropoint,
            "zeropoint_source": zeropoint_source,
            "magnitude_definition": (
                "m_cal=-2.5*log10(total_source_electrons)+zeropoint; "
                "color_term=0 for injected stars"
            ),
            **zeropoint_meta,
        },
        "empirical_psf": psf_meta,
        "baseline_detection": baseline_meta,
        "command": sys.argv,
    }
    _json_dump(run_dir / "manifest.json", manifest)
    _write_plots(stars, bins, run_dir / "reports")
    _write_precision_plot(points, completeness_fit, run_dir / "reports")
    _write_markdown_report(
        run_dir / "reports" / "README.md",
        summary,
        config,
        fwhm_px,
        psf_meta,
    )
    return run_dir
