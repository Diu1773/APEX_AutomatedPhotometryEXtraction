"""Combine multiple CMD validation reports into paper-facing tables and plots."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class ValidationInput:
    label: str
    path: str
    role: str = ""


@dataclass
class CmdCombinedValidationConfig:
    inputs: list[ValidationInput] = field(default_factory=list)
    output_root: str = "benchmark/runs/cmd_validation_combined"


def parse_validation_input(value: str) -> ValidationInput:
    """Parse label=path or label:role=path syntax."""
    if "=" not in value:
        path = Path(value)
        return ValidationInput(label=path.name, path=str(path), role="")
    left, path = value.split("=", 1)
    if ":" in left:
        label, role = left.split(":", 1)
    else:
        label, role = left, ""
    label = label.strip()
    if not label:
        raise ValueError(f"Validation input has an empty label: {value!r}")
    return ValidationInput(label=label, role=role.strip(), path=path.strip())


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _read_table(root: Path, name: str, dataset: ValidationInput) -> pd.DataFrame:
    path = root / name
    if not path.is_file():
        raise FileNotFoundError(f"Required validation table is missing: {path}")
    table = pd.read_csv(path)
    table.insert(0, "dataset", dataset.label)
    table.insert(1, "dataset_role", dataset.role)
    table.insert(2, "validation_root", str(root))
    return table


def collect_combined_tables(
    inputs: list[ValidationInput],
) -> dict[str, pd.DataFrame]:
    if not inputs:
        raise ValueError("At least one validation input is required")
    table_names = {
        "precision": "precision.csv",
        "m50_loss": "ast_m50_loss.csv",
        "false_positive": "false_positive_summary.csv",
        "repeatability_bins": "repeatability_by_magnitude.csv",
        "zeropoint_residuals": "zeropoint_residual_summary.csv",
        "frame_zeropoints": "frame_zeropoint_summary.csv",
    }
    collected: dict[str, list[pd.DataFrame]] = {key: [] for key in table_names}
    for dataset in inputs:
        root = Path(dataset.path).expanduser().resolve()
        for key, filename in table_names.items():
            collected[key].append(_read_table(root, filename, dataset))
    return {key: pd.concat(parts, ignore_index=True) for key, parts in collected.items()}


def build_dataset_summary(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    precision = tables["precision"]
    loss = tables["m50_loss"]
    false_positive = tables["false_positive"]
    repeatability = tables["repeatability_bins"]
    zp = tables["zeropoint_residuals"]
    for dataset, group in precision.groupby("dataset", sort=False):
        role = str(group["dataset_role"].iloc[0])
        fp = false_positive[false_positive["dataset"] == dataset]
        rep = repeatability[repeatability["dataset"] == dataset]
        zpg = zp[zp["dataset"] == dataset]
        lossg = loss[loss["dataset"] == dataset]
        rows.append(
            {
                "dataset": dataset,
                "dataset_role": role,
                "n_precision_frames": int(len(group)),
                "m50_loss_min": float(_num(lossg["m50_loss_best_minus_worst"]).min()),
                "m50_loss_median": float(_num(lossg["m50_loss_best_minus_worst"]).median()),
                "m50_loss_max": float(_num(lossg["m50_loss_best_minus_worst"]).max()),
                "false_per_1000_injected": float(
                    1000.0
                    * _num(fp["new_false_detections"]).sum()
                    / max(1.0, _num(fp["n_injected"]).sum())
                ),
                "repeatability_median_binned_mad": float(
                    _num(rep["repeatability_mad_median"]).median()
                ),
                "zp_clipped_rms_median": float(_num(zpg["residual_rms_clipped"]).median()),
                "zp_fit_scatter_median": float(_num(zpg["fit_scatter_rms"]).median()),
            }
        )
    return pd.DataFrame(rows)


def build_filter_summary(tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    loss = tables["m50_loss"].copy()
    false_positive = tables["false_positive"].copy()
    repeatability = tables["repeatability_bins"].copy()
    zp = tables["zeropoint_residuals"].copy()
    rows = []
    for (dataset, filter_name), group in loss.groupby(["dataset", "filter"], sort=True):
        fp = false_positive[
            (false_positive["dataset"] == dataset) & (false_positive["filter"] == filter_name)
        ]
        rep = repeatability[
            (repeatability["dataset"] == dataset) & (repeatability["filter"] == filter_name)
        ]
        zpg = zp[(zp["dataset"] == dataset) & (zp["filter"] == filter_name)]
        rows.append(
            {
                "dataset": dataset,
                "dataset_role": group["dataset_role"].iloc[0],
                "filter": filter_name,
                "m50_best": float(group["m50_best"].iloc[0]),
                "m50_worst": float(group["m50_worst"].iloc[0]),
                "m50_loss": float(group["m50_loss_best_minus_worst"].iloc[0]),
                "false_per_1000_injected": float(
                    1000.0
                    * _num(fp["new_false_detections"]).sum()
                    / max(1.0, _num(fp["n_injected"]).sum())
                ),
                "repeatability_binned_mad_median": float(
                    _num(rep["repeatability_mad_median"]).median()
                ),
                "zp_clipped_rms": float(_num(zpg["residual_rms_clipped"]).iloc[0])
                if len(zpg)
                else np.nan,
                "zp_fit_scatter": float(_num(zpg["fit_scatter_rms"]).iloc[0])
                if len(zpg)
                else np.nan,
            }
        )
    return pd.DataFrame(rows).sort_values(["dataset", "filter"])


def _write_m50_loss_plot(filter_summary: pd.DataFrame, output_dir: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    datasets = list(dict.fromkeys(filter_summary["dataset"]))
    filters = sorted(filter_summary["filter"].unique())
    x = np.arange(len(filters))
    width = 0.8 / max(1, len(datasets))
    for index, dataset in enumerate(datasets):
        group = filter_summary[filter_summary["dataset"] == dataset].set_index("filter")
        values = [group.loc[f, "m50_loss"] if f in group.index else np.nan for f in filters]
        ax.bar(x + (index - (len(datasets) - 1) / 2) * width, values, width, label=dataset)
    ax.set_xticks(x)
    ax.set_xticklabels(filters)
    ax.set_ylabel("m50 loss best-to-worst seeing (mag)")
    ax.set_xlabel("Filter")
    ax.grid(alpha=0.25, axis="y")
    ax.legend(title="Dataset")
    fig.tight_layout()
    fig.savefig(output_dir / "combined_m50_loss_by_filter.png", dpi=180)
    plt.close(fig)


def _write_m50_vs_fwhm_plot(precision: pd.DataFrame, output_dir: Path) -> None:
    datasets = list(dict.fromkeys(precision["dataset"]))
    fig, axes = plt.subplots(1, len(datasets), figsize=(7.0 * len(datasets), 5.0), sharey=True)
    axes = np.atleast_1d(axes)
    markers = {"g": "o", "r": "s", "i": "^"}
    for ax, dataset in zip(axes, datasets):
        sub = precision[precision["dataset"] == dataset].copy()
        for filter_name, group in sub.groupby("filter", sort=True):
            group = group.sort_values("fwhm_px_benchmark")
            y = _num(group["m50"])
            yerr = np.vstack(
                [
                    np.maximum(0.0, y - _num(group["m50_ci95_low"])),
                    np.maximum(0.0, _num(group["m50_ci95_high"]) - y),
                ]
            )
            ax.errorbar(
                _num(group["fwhm_px_benchmark"]),
                y,
                yerr=yerr,
                marker=markers.get(str(filter_name), "o"),
                linewidth=1.6,
                capsize=3,
                label=str(filter_name),
            )
        role = str(sub["dataset_role"].iloc[0])
        title = dataset if not role else f"{dataset}\n{role}"
        ax.set_title(title)
        ax.set_xlabel("Benchmark FWHM (px)")
        ax.grid(alpha=0.25)
        ax.legend(title="Filter")
    axes[0].set_ylabel("50% completeness magnitude")
    fig.tight_layout()
    fig.savefig(output_dir / "combined_m50_vs_fwhm.png", dpi=180)
    plt.close(fig)


def _write_metric_radar_plot(dataset_summary: pd.DataFrame, output_dir: Path) -> None:
    # A compact normalized comparison. Values are not a score; each axis is scaled to max.
    metrics = [
        "m50_loss_median",
        "false_per_1000_injected",
        "repeatability_median_binned_mad",
        "zp_clipped_rms_median",
    ]
    labels = ["m50 loss", "false+", "repeatability", "ZP RMS"]
    values = dataset_summary[metrics].astype(float).copy()
    maxima = values.max(axis=0).replace(0, np.nan)
    normalized = values / maxima
    angles = np.linspace(0, 2 * np.pi, len(metrics), endpoint=False).tolist()
    angles += angles[:1]
    fig, ax = plt.subplots(figsize=(6.2, 6.2), subplot_kw={"projection": "polar"})
    for idx, row in dataset_summary.reset_index(drop=True).iterrows():
        vals = normalized.iloc[idx].fillna(0).tolist()
        vals += vals[:1]
        ax.plot(angles, vals, linewidth=1.8, label=row["dataset"])
        ax.fill(angles, vals, alpha=0.08)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels)
    ax.set_yticklabels([])
    ax.set_title("Normalized validation stress metrics")
    ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.1))
    fig.tight_layout()
    fig.savefig(output_dir / "combined_validation_stress_metrics.png", dpi=180)
    plt.close(fig)


def _write_repeatability_plot(repeatability: pd.DataFrame, output_dir: Path) -> None:
    datasets = list(dict.fromkeys(repeatability["dataset"]))
    fig, axes = plt.subplots(1, len(datasets), figsize=(7.0 * len(datasets), 4.8), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, dataset in zip(axes, datasets):
        sub = repeatability[repeatability["dataset"] == dataset]
        for filter_name, group in sub.groupby("filter", sort=True):
            group = group.sort_values("mag_median")
            ax.plot(
                _num(group["mag_median"]),
                _num(group["repeatability_mad_median"]),
                marker="o",
                linewidth=1.5,
                label=str(filter_name),
            )
        ax.set_yscale("log")
        ax.set_title(str(dataset))
        ax.set_xlabel("Median calibrated magnitude")
        ax.grid(alpha=0.25, which="both")
        ax.legend(title="Filter")
    axes[0].set_ylabel("Repeated-frame robust scatter (mag)")
    fig.tight_layout()
    fig.savefig(output_dir / "combined_repeatability.png", dpi=180)
    plt.close(fig)


def _write_report(
    output_dir: Path,
    *,
    config: CmdCombinedValidationConfig,
    dataset_summary: pd.DataFrame,
    filter_summary: pd.DataFrame,
) -> None:
    lines = [
        "# Combined CMD Validation Report",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Datasets",
        "",
    ]
    for item in config.inputs:
        role = f" - {item.role}" if item.role else ""
        lines.append(f"- {item.label}{role}: {item.path}")
    lines.extend(["", "## Dataset Summary", ""])
    for _, row in dataset_summary.iterrows():
        role = f" ({row['dataset_role']})" if row["dataset_role"] else ""
        lines.append(
            f"- {row['dataset']}{role}: median m50 loss "
            f"{row['m50_loss_median']:.3f} mag, false+ "
            f"{row['false_per_1000_injected']:.2f}/1000, repeatability "
            f"{row['repeatability_median_binned_mad']:.4f} mag, ZP RMS "
            f"{row['zp_clipped_rms_median']:.4f} mag"
        )

    lines.extend(["", "## Filter-Level Summary", ""])
    for _, row in filter_summary.iterrows():
        lines.append(
            f"- {row['dataset']} {row['filter']}: m50 loss {row['m50_loss']:.3f} mag, "
            f"false+ {row['false_per_1000_injected']:.2f}/1000, "
            f"repeatability {row['repeatability_binned_mad_median']:.4f} mag"
        )

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- dataset_summary.csv",
            "- filter_summary.csv",
            "- precision_combined.csv",
            "- combined_m50_loss_by_filter.png",
            "- combined_m50_vs_fwhm.png",
            "- combined_repeatability.png",
            "- combined_validation_stress_metrics.png",
            "",
            "## Configuration",
            "",
            "```json",
            json.dumps(asdict(config), indent=2, sort_keys=True),
            "```",
            "",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def run_combined_validation(
    config: CmdCombinedValidationConfig,
    *,
    output_override: str | Path | None = None,
) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    output_root = Path(output_override or config.output_root).expanduser()
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    tables = collect_combined_tables(config.inputs)
    for name, table in tables.items():
        table.to_csv(output_root / f"{name}_combined.csv", index=False)

    dataset_summary = build_dataset_summary(tables)
    filter_summary = build_filter_summary(tables)
    dataset_summary.to_csv(output_root / "dataset_summary.csv", index=False)
    filter_summary.to_csv(output_root / "filter_summary.csv", index=False)

    _write_m50_loss_plot(filter_summary, output_root)
    _write_m50_vs_fwhm_plot(tables["precision"], output_root)
    _write_repeatability_plot(tables["repeatability_bins"], output_root)
    _write_metric_radar_plot(dataset_summary, output_root)
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "effective_config": asdict(config),
    }
    (output_root / "combined_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_report(
        output_root,
        config=config,
        dataset_summary=dataset_summary,
        filter_summary=filter_summary,
    )
    return output_root
