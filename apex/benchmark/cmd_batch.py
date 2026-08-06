"""Batch orchestration for CMD artificial-star benchmarks."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

try:
    import tomllib
except ImportError:  # pragma: no cover - Python 3.10
    import tomli as tomllib

from apex.benchmark.runner import load_benchmark_config, run_benchmark


@dataclass
class CmdBatchConfig:
    project_root: str = ""
    dataset_name: str = ""
    benchmark_config: str = "benchmark/configs/baseline.toml"
    output_root: str = "benchmark/runs/cmd_batch"
    filters: list[str] = field(default_factory=lambda: ["g", "r", "i"])
    quality_quantiles: dict[str, float] = field(
        default_factory=lambda: {"best": 0.0, "median": 0.5, "worst": 1.0}
    )
    coarse_offsets_from_zeropoint: list[float] = field(
        default_factory=lambda: [
            -12.0,
            -11.5,
            -11.0,
            -10.5,
            -10.0,
            -9.5,
            -9.0,
            -8.5,
            -8.0,
            -7.5,
            -7.0,
            -6.5,
            -6.0,
        ]
    )
    coarse_trials: int = 5
    coarse_stars_per_magnitude_per_trial: int = 5
    coarse_bootstrap_samples: int = 300
    precision_half_width_mag: float = 0.6
    precision_step_mag: float = 0.2
    precision_trials: int = 20
    precision_stars_per_magnitude_per_trial: int = 10
    precision_bootstrap_samples: int = 1000
    save_injected_fits: bool = False


def load_cmd_batch_config(path: str | Path) -> CmdBatchConfig:
    config_path = Path(path)
    from apex.config.config_io import load_config_data
    raw, config_path = load_config_data(config_path)  # JSON authority
    section = dict(raw.get("cmd_batch", raw))
    known = set(CmdBatchConfig.__dataclass_fields__)
    unknown = sorted(set(section) - known)
    if unknown:
        raise ValueError(f"Unknown CMD batch config keys: {', '.join(unknown)}")
    return CmdBatchConfig(**section)


def _normal_filter(value: Any) -> str:
    return str(value).strip().lower().rstrip("'")


def _resolve_frame_path(project_root: Path, filename: str) -> Path:
    direct = project_root / filename
    if direct.is_file():
        return direct.resolve()
    matches = [path for path in project_root.rglob(filename) if path.is_file()]
    if len(matches) == 1:
        return matches[0].resolve()
    if not matches:
        raise FileNotFoundError(f"CMD frame is missing: {project_root / filename}")
    raise RuntimeError(f"CMD frame name is ambiguous under {project_root}: {filename}")


def select_cmd_frames(
    project_root: str | Path,
    *,
    filters: list[str] | tuple[str, ...] = ("g", "r", "i"),
    quality_quantiles: dict[str, float] | None = None,
) -> pd.DataFrame:
    """Select unique frames nearest requested within-filter FWHM quantiles."""
    root = Path(project_root).expanduser().resolve()
    result_dir = root / "result"
    index_path = result_dir / "step7_forced_phot" / "photometry_index.csv"
    zeropoint_path = result_dir / "cmd_zeropoint" / "frame_zeropoint.csv"
    if not index_path.is_file():
        raise FileNotFoundError(f"Step 7 photometry index is missing: {index_path}")
    if not zeropoint_path.is_file():
        raise FileNotFoundError(f"Step 10 frame zeropoints are missing: {zeropoint_path}")

    index = pd.read_csv(index_path)
    zeropoints = pd.read_csv(zeropoint_path)
    required_index = {"file", "filter", "fwhm_px"}
    required_zp = {"file", "filter", "zp_frame"}
    if not required_index <= set(index.columns):
        missing = sorted(required_index - set(index.columns))
        raise ValueError(f"Step 7 index is missing columns: {', '.join(missing)}")
    if not required_zp <= set(zeropoints.columns):
        missing = sorted(required_zp - set(zeropoints.columns))
        raise ValueError(f"Step 10 table is missing columns: {', '.join(missing)}")

    index = index.copy()
    zeropoints = zeropoints.copy()
    index["filter_key"] = index["filter"].map(_normal_filter)
    zeropoints["filter_key"] = zeropoints["filter"].map(_normal_filter)
    index["fwhm_px"] = pd.to_numeric(index["fwhm_px"], errors="coerce")
    zeropoints["zp_frame"] = pd.to_numeric(zeropoints["zp_frame"], errors="coerce")
    if "status" in index.columns:
        index = index[index["status"].astype(str).str.strip().str.lower() == "ok"]
    if "wcs_ok" in index.columns:
        wcs_text = index["wcs_ok"].astype(str).str.strip().str.lower()
        index = index[wcs_text.isin({"true", "1", "yes"})]

    zp_columns = ["file", "filter_key", "zp_frame"]
    for optional in ("zp_scatter", "n_ref", "snr_med", "airmass"):
        if optional in zeropoints.columns:
            zp_columns.append(optional)
    merged = index.merge(
        zeropoints[zp_columns],
        on=["file", "filter_key"],
        how="inner",
        validate="one_to_one",
    )
    merged = merged[
        np.isfinite(merged["fwhm_px"]) & np.isfinite(merged["zp_frame"])
    ].copy()

    quantiles = quality_quantiles or {"best": 0.0, "median": 0.5, "worst": 1.0}
    if not quantiles:
        raise ValueError("At least one quality quantile is required")
    for label, quantile in quantiles.items():
        if not label or not 0.0 <= float(quantile) <= 1.0:
            raise ValueError(f"Invalid quality quantile {label!r}: {quantile}")

    rows: list[dict[str, Any]] = []
    for filter_name in filters:
        filter_key = _normal_filter(filter_name)
        group = merged[merged["filter_key"] == filter_key].copy()
        if len(group) < len(quantiles):
            raise RuntimeError(
                f"Filter {filter_key} has {len(group)} calibrated frames; "
                f"{len(quantiles)} unique selections were requested"
            )
        available = set(group.index)
        for label, quantile in quantiles.items():
            target = float(group["fwhm_px"].quantile(float(quantile)))
            candidates = group.loc[list(available)].copy()
            candidates["distance_to_target"] = (candidates["fwhm_px"] - target).abs()
            selected = candidates.sort_values(
                ["distance_to_target", "fwhm_px", "file"],
                kind="stable",
            ).iloc[0]
            available.remove(selected.name)
            row = selected.to_dict()
            row.update(
                {
                    "condition": str(label),
                    "quality_quantile": float(quantile),
                    "target_fwhm_px": target,
                    "input_fits": str(_resolve_frame_path(root, str(selected["file"]))),
                }
            )
            rows.append(row)

    selected = pd.DataFrame(rows)
    leading = [
        "filter_key",
        "condition",
        "quality_quantile",
        "file",
        "input_fits",
        "fwhm_px",
        "target_fwhm_px",
        "zp_frame",
    ]
    return selected[leading + [column for column in selected.columns if column not in leading]]


def _precision_grid(m50: float, half_width: float, step: float) -> list[float]:
    if not np.isfinite(m50):
        raise ValueError("Coarse m50 must be finite")
    if half_width <= 0 or step <= 0:
        raise ValueError("Precision half-width and step must be positive")
    count = int(round((2.0 * half_width) / step)) + 1
    if count < 3:
        raise ValueError("Precision magnitude grid must contain at least three points")
    return np.linspace(m50 - half_width, m50 + half_width, count).round(6).tolist()


def _read_summary(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "summary.json"
    if not path.is_file():
        raise RuntimeError(f"Benchmark summary is missing: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def _read_magnitude_points(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "magnitude_points.csv"
    if not path.is_file():
        raise RuntimeError(f"Benchmark magnitude points are missing: {path}")
    points = pd.read_csv(path)
    required = {"magnitude", "completeness"}
    if not required <= set(points.columns):
        missing = sorted(required - set(points.columns))
        raise RuntimeError(f"Magnitude points are missing columns: {', '.join(missing)}")
    points = points.copy()
    points["magnitude"] = pd.to_numeric(points["magnitude"], errors="coerce")
    points["completeness"] = pd.to_numeric(points["completeness"], errors="coerce")
    points = points[
        np.isfinite(points["magnitude"]) & np.isfinite(points["completeness"])
    ].sort_values("magnitude")
    if points.empty:
        raise RuntimeError("Magnitude points do not contain finite completeness values")
    return points


def _coarse_m50_seed(points: pd.DataFrame, fit: dict[str, Any], frame_name: str) -> float:
    """Choose a precision-grid center only when the coarse points bracket 50%."""
    min_completeness = float(points["completeness"].min())
    max_completeness = float(points["completeness"].max())
    min_mag = float(points["magnitude"].min())
    max_mag = float(points["magnitude"].max())
    if min_completeness > 0.5:
        raise RuntimeError(
            f"Coarse magnitude grid is too bright for {frame_name}: "
            f"completeness stays above 50% through {max_mag:.3f} mag"
        )
    if max_completeness < 0.5:
        raise RuntimeError(
            f"Coarse magnitude grid is too faint for {frame_name}: "
            f"completeness stays below 50% from {min_mag:.3f} mag"
        )

    fit_m50 = fit.get("m50")
    if fit_m50 is not None and np.isfinite(float(fit_m50)):
        fit_m50 = float(fit_m50)
        if min_mag <= fit_m50 <= max_mag:
            return fit_m50

    mags = points["magnitude"].to_numpy(float)
    completeness = points["completeness"].to_numpy(float)
    for index in range(len(points) - 1):
        c0 = completeness[index]
        c1 = completeness[index + 1]
        if (c0 - 0.5) == 0:
            return float(mags[index])
        if (c0 - 0.5) * (c1 - 0.5) <= 0:
            if c0 == c1:
                return float(0.5 * (mags[index] + mags[index + 1]))
            fraction = (0.5 - c0) / (c1 - c0)
            return float(mags[index] + fraction * (mags[index + 1] - mags[index]))
    nearest = int(np.argmin(np.abs(completeness - 0.5)))
    return float(mags[nearest])


def _summary_row(
    selection: pd.Series,
    *,
    stage: str,
    run_dir: Path,
    summary: dict[str, Any],
) -> dict[str, Any]:
    fit = summary.get("completeness_fit") or {}
    return {
        "filter": selection["filter_key"],
        "condition": selection["condition"],
        "file": selection["file"],
        "input_fits": selection["input_fits"],
        "fwhm_px_step7": float(selection["fwhm_px"]),
        "fwhm_px_benchmark": summary.get("fwhm_px"),
        "zp_frame": float(selection["zp_frame"]),
        "stage": stage,
        "m90": fit.get("m90"),
        "m50": fit.get("m50"),
        "m50_ci95_low": fit.get("m50_ci95_low"),
        "m50_ci95_high": fit.get("m50_ci95_high"),
        "m10": fit.get("m10"),
        "completeness": summary.get("completeness"),
        "forced_mag_bias_median": summary.get("forced_mag_bias_median"),
        "forced_mag_scatter_mad": summary.get("forced_mag_scatter_mad"),
        "run_dir": str(run_dir),
    }


def _write_selection_plot(selected: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for filter_name, group in selected.groupby("filter_key", sort=True):
        ordered = group.sort_values("quality_quantile")
        ax.plot(
            ordered["condition"],
            ordered["fwhm_px"],
            marker="o",
            linewidth=1.5,
            label=str(filter_name),
        )
    ax.set(xlabel="Selected observing condition", ylabel="Step 7 FWHM (px)")
    ax.grid(alpha=0.25)
    ax.legend(title="Filter")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def _write_result_plot(results: pd.DataFrame, path: Path) -> None:
    precision = results[results["stage"] == "precision"].copy()
    if precision.empty:
        return
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    markers = {"g": "o", "r": "s", "i": "^"}
    for filter_name, group in precision.groupby("filter", sort=True):
        group = group.sort_values("fwhm_px_benchmark")
        y = group["m50"].astype(float)
        yerr = np.vstack(
            [
                np.maximum(0.0, y - group["m50_ci95_low"].astype(float)),
                np.maximum(0.0, group["m50_ci95_high"].astype(float) - y),
            ]
        )
        ax.errorbar(
            group["fwhm_px_benchmark"],
            y,
            yerr=yerr,
            marker=markers.get(str(filter_name), "o"),
            linewidth=1.7,
            capsize=3,
            label=str(filter_name),
        )
        for _, row in group.iterrows():
            ax.annotate(
                str(row["condition"]),
                (float(row["fwhm_px_benchmark"]), float(row["m50"])),
                textcoords="offset points",
                xytext=(5, 5),
                fontsize=8,
            )
    ax.set(xlabel="Benchmark FWHM (px)", ylabel="50% completeness magnitude")
    ax.grid(alpha=0.25)
    ax.legend(title="Filter")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run_cmd_batch(
    config: CmdBatchConfig,
    *,
    project_override: str | Path | None = None,
    output_override: str | Path | None = None,
    select_only: bool = False,
) -> Path:
    """Select CMD frames and optionally run coarse plus precision benchmarks."""
    repo_root = Path(__file__).resolve().parents[2]
    project_root = Path(project_override or config.project_root).expanduser().resolve()
    if not project_root.is_dir():
        raise FileNotFoundError(f"CMD project root does not exist: {project_root}")
    output_root = Path(output_override or config.output_root).expanduser()
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    selected = select_cmd_frames(
        project_root,
        filters=config.filters,
        quality_quantiles=config.quality_quantiles,
    )
    selected.to_csv(output_root / "selected_frames.csv", index=False)
    _write_selection_plot(selected, output_root / "selected_fwhm.png")
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_name": config.dataset_name or project_root.name,
        "project_root": str(project_root),
        "effective_config": asdict(config),
        "select_only": bool(select_only),
    }
    (output_root / "cmd_batch_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    if select_only:
        return output_root

    benchmark_config_path = Path(config.benchmark_config)
    if not benchmark_config_path.is_absolute():
        benchmark_config_path = repo_root / benchmark_config_path
    base = load_benchmark_config(benchmark_config_path)
    result_rows: list[dict[str, Any]] = []
    frames_root = output_root / "frames"
    for index, selection in selected.reset_index(drop=True).iterrows():
        frame_root = (
            frames_root
            / str(selection["filter_key"])
            / f"{selection['condition']}_{Path(str(selection['file'])).stem}"
        )
        if frame_root.exists() and any(frame_root.iterdir()):
            raise FileExistsError(
                f"CMD batch frame output already exists; choose a new output directory: {frame_root}"
            )
        zeropoint = float(selection["zp_frame"])
        coarse_grid = sorted(
            zeropoint + float(offset)
            for offset in config.coarse_offsets_from_zeropoint
        )
        seed = int(base.seed + index * 100_000)
        coarse_config = replace(
            base,
            seed=seed,
            trials=int(config.coarse_trials),
            magnitude_grid=coarse_grid,
            stars_per_magnitude_per_trial=int(
                config.coarse_stars_per_magnitude_per_trial
            ),
            completeness_bootstrap_samples=int(config.coarse_bootstrap_samples),
            save_injected_fits=bool(config.save_injected_fits),
            zeropoint_mag=zeropoint,
            allow_initial_zeropoint_fallback=False,
        )
        coarse_dir = run_benchmark(
            coarse_config,
            input_override=selection["input_fits"],
            output_override=frame_root / "coarse",
        )
        coarse_summary = _read_summary(coarse_dir)
        coarse_points = _read_magnitude_points(coarse_dir)
        result_rows.append(
            _summary_row(
                selection,
                stage="coarse",
                run_dir=coarse_dir,
                summary=coarse_summary,
            )
        )
        fit = coarse_summary.get("completeness_fit") or {}
        coarse_m50 = _coarse_m50_seed(coarse_points, fit, str(selection["file"]))

        precision_config = replace(
            base,
            seed=seed + 50_000,
            trials=int(config.precision_trials),
            magnitude_grid=_precision_grid(
                float(coarse_m50),
                float(config.precision_half_width_mag),
                float(config.precision_step_mag),
            ),
            stars_per_magnitude_per_trial=int(
                config.precision_stars_per_magnitude_per_trial
            ),
            completeness_bootstrap_samples=int(config.precision_bootstrap_samples),
            save_injected_fits=bool(config.save_injected_fits),
            zeropoint_mag=zeropoint,
            allow_initial_zeropoint_fallback=False,
        )
        precision_dir = run_benchmark(
            precision_config,
            input_override=selection["input_fits"],
            output_override=frame_root / "precision",
        )
        precision_summary = _read_summary(precision_dir)
        result_rows.append(
            _summary_row(
                selection,
                stage="precision",
                run_dir=precision_dir,
                summary=precision_summary,
            )
        )
        results = pd.DataFrame(result_rows)
        results.to_csv(output_root / "cmd_batch_summary.csv", index=False)
        results[results["stage"] == "precision"].to_csv(
            output_root / "cmd_precision_summary.csv",
            index=False,
        )
        _write_result_plot(results, output_root / "m50_vs_fwhm.png")

    return output_root
