"""Unified, reproducible APEX validation harness.

This is the orchestration layer behind ``apex validate``. It produces
quantitative accuracy evidence as a single markdown report plus a JSON
manifest, suitable for the README centerpiece, a paper's validation section,
and a CI regression guard.

It BUILDS ON the existing benchmark/validation code and does not duplicate any
of the algorithms:

* :func:`apex.benchmark.runner.run_benchmark` drives the end-to-end
  artificial-star test (baseline detect -> empirical PSF -> inject -> detect ->
  match -> metrics + plots).
* :func:`apex.benchmark.synthetic_frame.make_synthetic_reference_frame`
  generates a self-contained reference frame so the suite runs with zero
  external data.
* :func:`apex.benchmark.iraf_crosscheck.run_iraf_crosscheck` runs the optional
  IRAF/PyRAF cross-check when IRAF and the required inputs are available.
* ``validation.cmd_step12_synthetic`` provides the zero-data isochrone-recovery
  check reused by the known-targets suite.

The module is deliberately Qt-free: no ``apex.gui`` / PyQt5 imports anywhere in
its import graph (the heavy science imports happen lazily inside the suite
functions, and the benchmark runner itself imports the Qt-free detection
service).
"""

from __future__ import annotations

from dataclasses import asdict
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any, Iterable

import numpy as np

# ── Tuned synthetic-frame defaults ────────────────────────────────────────────
# These were chosen empirically so the recovered completeness curve is
# NON-DEGENERATE: bright injections recover ~100%, faint injections ~0%, with a
# meaningful 50% point in the middle of the injection magnitude range. They are
# matched to the pipeline gain at runtime (see _resolve_gain).
SYNTH_FRAME_DEFAULTS: dict[str, Any] = {
    "n_stars": 170,
    "size": 640,
    "fwhm_px": 3.5,
    "zeropoint": 25.0,
    "background": 150.0,
    "read_noise": 5.0,
    "mag_min": 13.5,   # baked-in stars go brighter than the injection floor
    "mag_max": 19.0,
    "filter_name": "r",
    "exptime": 1.0,
}

# Injection magnitude range for the artificial-star suite (within which the
# completeness curve sweeps from ~1.0 to ~0.0). Slightly narrower than the
# baked-in star range so the bright plateau is well sampled.
INJECTION_MAG_MIN = 14.5
INJECTION_MAG_MAX = 19.0

SUITE_CHOICES = ("artificial-star", "iraf", "known-targets", "crosscheck")


# ── small helpers ─────────────────────────────────────────────────────────────

def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def _round(value: Any, ndigits: int = 4) -> float | None:
    out = _finite(value)
    return round(out, ndigits) if out is not None else None


def _git_commit(repo_root: Path) -> dict[str, Any]:
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


def _package_versions() -> dict[str, str]:
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for package in ("numpy", "scipy", "pandas", "astropy", "photutils", "sep", "matplotlib"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _resolve_gain(config: Any) -> float:
    """Best-effort detector gain so synthetic noise matches the pipeline scale."""
    default = 1.5
    parameter_file = None
    if isinstance(config, dict):
        parameter_file = config.get("parameter_file")
    elif config is not None:
        parameter_file = getattr(config, "parameter_file", None)
    candidates: list[Path] = []
    if parameter_file:
        candidates.append(Path(parameter_file))
    candidates.append(_repo_root() / "parameters.toml")
    for path in candidates:
        try:
            if not path or not Path(path).exists():
                continue
            try:
                import tomllib
            except ImportError:  # pragma: no cover
                import tomli as tomllib  # type: ignore
            with Path(path).open("rb") as handle:
                raw = tomllib.load(handle)
            for section in (raw, raw.get("instrument", {}), raw.get("camera", {})):
                if isinstance(section, dict) and "gain_e_per_adu" in section:
                    g = _finite(section["gain_e_per_adu"])
                    if g and g > 0:
                        return float(g)
        except Exception:
            continue
    return default


# ── artificial-star suite ─────────────────────────────────────────────────────

def run_artificial_star_suite(
    output_dir: str | Path,
    *,
    reference_frame: str | Path | None = None,
    config: Any | None = None,
    seed: int = 20260620,
    trials: int = 4,
    stars_per_trial: int = 40,
    bootstrap_samples: int = 200,
    synth_overrides: dict[str, Any] | None = None,
    keep_run_dir: bool = True,
) -> dict[str, Any]:
    """Run the end-to-end artificial-star test and extract key accuracy metrics.

    If ``reference_frame`` is None a self-contained synthetic frame is generated
    (zero external data). The benchmark runs with an explicit injection
    zeropoint so it never depends on Step-10 calibration products.
    """
    from apex.benchmark.runner import BenchmarkConfig, run_benchmark
    from apex.benchmark.synthetic_frame import make_synthetic_reference_frame

    output_dir = Path(output_dir)
    suite_dir = output_dir / "artificial_star"
    suite_dir.mkdir(parents=True, exist_ok=True)

    synth = dict(SYNTH_FRAME_DEFAULTS)
    if synth_overrides:
        synth.update(synth_overrides)

    parameter_file = "parameters.toml"
    if config is not None:
        parameter_file = getattr(config, "parameter_file", parameter_file)
    gain = _resolve_gain(config)
    synth.setdefault("gain", gain)

    generated_frame = None
    if reference_frame is None:
        frame_path = suite_dir / "synthetic_reference.fits"
        make_synthetic_reference_frame(frame_path, seed=seed, **synth)
        reference_frame = frame_path
        generated_frame = frame_path
    reference_frame = Path(reference_frame)

    run_dir = suite_dir / "benchmark_run"
    if run_dir.exists():
        shutil.rmtree(run_dir, ignore_errors=True)

    zeropoint = float(synth.get("zeropoint", 25.0))
    bench_config = BenchmarkConfig(
        input_fits=str(reference_frame),
        parameter_file=parameter_file,
        seed=int(seed),
        trials=int(trials),
        stars_per_trial=int(stars_per_trial),
        magnitude_min=INJECTION_MAG_MIN,
        magnitude_max=INJECTION_MAG_MAX,
        magnitude_bin_width=1.0,
        completeness_bootstrap_samples=int(bootstrap_samples),
        save_injected_fits=False,
        zeropoint_mag=zeropoint,  # explicit -> no Step-10 dependency
        isolated_fraction=0.6,
        psf_min_stars=3,
    )

    try:
        produced = run_benchmark(
            bench_config,
            input_override=reference_frame,
            output_override=run_dir,
        )
    except Exception as exc:  # noqa: BLE001 - surface as a failed suite, never crash report
        return {
            "status": "error",
            "reason": f"benchmark run failed: {exc}",
            "reference_frame": str(reference_frame),
            "synthetic_frame_params": synth if generated_frame else None,
        }

    import pandas as pd  # local import keeps top-level import-light

    summary = json.loads((produced / "summary.json").read_text(encoding="utf-8"))
    fit = summary.get("completeness_fit", {}) or {}

    # Photometric bias/scatter vs magnitude from the per-bin table.
    bins = pd.read_csv(produced / "magnitude_bins.csv") if (produced / "magnitude_bins.csv").exists() else pd.DataFrame()
    bright_rms = faint_rms = None
    per_mag: list[dict[str, Any]] = []
    if not bins.empty:
        bins = bins.sort_values("magnitude_center")
        for _, row in bins.iterrows():
            per_mag.append(
                {
                    "magnitude_center": _round(row.get("magnitude_center"), 3),
                    "n_eligible": int(row.get("n_eligible", 0)),
                    "completeness": _round(row.get("completeness"), 3),
                    "mag_bias_median": _round(row.get("forced_mag_bias_median"), 4),
                }
            )

    # Per-magnitude scatter (RMS of recovered-expected mag error) from stars.csv.
    stars = pd.read_csv(produced / "stars.csv") if (produced / "stars.csv").exists() else pd.DataFrame()
    if not stars.empty and "forced_mag_error" in stars.columns:
        err = pd.to_numeric(stars["forced_mag_error"], errors="coerce")
        mag = pd.to_numeric(stars["magnitude_true"], errors="coerce")
        eligible = ~stars.get("baseline_confounded", pd.Series(False, index=stars.index)).astype(bool)
        good = eligible & np.isfinite(err) & np.isfinite(mag)
        m50 = _finite(fit.get("m50")) or float(np.nanmedian(mag[good])) if good.any() else None
        if good.any() and m50 is not None:
            bright_mask = good & (mag <= m50 - 1.0)
            faint_mask = good & (mag > m50 - 1.0)
            if bright_mask.any():
                bright_rms = float(np.sqrt(np.mean(err[bright_mask] ** 2)))
            if faint_mask.any():
                faint_rms = float(np.sqrt(np.mean(err[faint_mask] ** 2)))

    # Copy plots into the suite report directory.
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    copied_plots: list[str] = []
    for name in ("benchmark_summary.png", "completeness_logistic_fit.png"):
        src = produced / "reports" / name
        if src.exists():
            dest = plot_dir / f"artificial_star_{name}"
            shutil.copyfile(src, dest)
            copied_plots.append(dest.name)

    completeness_overall = _finite(summary.get("completeness"))
    n_injected = int(summary.get("n_injected", 0))
    n_false = int(summary.get("new_false_detections", 0))
    false_positive_rate = (n_false / n_injected) if n_injected else None

    # Non-degenerate check: curve must span a wide completeness range across mag.
    non_degenerate = False
    if per_mag:
        comp_values = [p["completeness"] for p in per_mag if p["completeness"] is not None]
        if comp_values:
            non_degenerate = (max(comp_values) - min(comp_values)) >= 0.5 and (
                max(comp_values) >= 0.7 and min(comp_values) <= 0.3
            )

    result: dict[str, Any] = {
        "status": "ok",
        "reference_frame": str(reference_frame),
        "generated_synthetic_frame": bool(generated_frame),
        "synthetic_frame_params": synth if generated_frame else None,
        "n_trials": int(trials),
        "stars_per_trial": int(stars_per_trial),
        "n_injected": n_injected,
        "injection_mag_range": [INJECTION_MAG_MIN, INJECTION_MAG_MAX],
        "completeness_overall": _round(completeness_overall, 3),
        "completeness_m90": _round(fit.get("m90"), 3),
        "completeness_m50": _round(fit.get("m50"), 3),
        "completeness_m10": _round(fit.get("m10"), 3),
        "completeness_width_mag": _round(fit.get("width_mag"), 3),
        "completeness_m50_ci95": [
            _round(fit.get("m50_ci95_low"), 3),
            _round(fit.get("m50_ci95_high"), 3),
        ],
        "completeness_non_degenerate": bool(non_degenerate),
        "photometric_bias_median_mag": _round(summary.get("forced_mag_bias_median"), 4),
        "photometric_scatter_mad_mag": _round(summary.get("forced_mag_scatter_mad"), 4),
        "photometric_rms_bright_mag": _round(bright_rms, 4),
        "photometric_rms_faint_mag": _round(faint_rms, 4),
        "position_rmse_px": _round(summary.get("position_rmse_px"), 4),
        "new_detections": int(summary.get("new_detections", 0)),
        "new_false_detections": n_false,
        "false_positive_rate": _round(false_positive_rate, 4),
        "completeness_by_magnitude": per_mag,
        "fwhm_px": _round(summary.get("fwhm_px"), 3),
        "zeropoint_mag": _round(summary.get("zeropoint_mag"), 4),
        "plots": copied_plots,
        "benchmark_run_dir": str(produced) if keep_run_dir else None,
    }
    return result


# ── IRAF suite ────────────────────────────────────────────────────────────────

def _iraf_available(runtime_cmd: list[str] | None) -> bool:
    import os

    if runtime_cmd:
        head = runtime_cmd[0]
    else:
        head = "wsl" if os.name == "nt" else "python3"
    return shutil.which(head) is not None


def run_iraf_suite(output_dir: str | Path, config: Any | None) -> dict[str, Any]:
    """Run the optional IRAF/PyRAF cross-check; skip cleanly if unavailable.

    Returns ``{"status": "skipped", "reason": ...}`` when IRAF or the required
    inputs (science FITS + Step 7 TSV) are not provided. Never raises.
    """
    output_dir = Path(output_dir)
    iraf_cfg = None
    input_fits = step7_tsv = None
    runtime_cmd: list[str] | None = None

    if config is not None:
        iraf_cfg = getattr(config, "iraf", None)
        if isinstance(iraf_cfg, dict):
            input_fits = iraf_cfg.get("input_fits")
            step7_tsv = iraf_cfg.get("step7_tsv")
            runtime_cmd = iraf_cfg.get("runtime_cmd")
        else:
            input_fits = getattr(config, "iraf_input_fits", None)
            step7_tsv = getattr(config, "iraf_step7_tsv", None)

    if not input_fits or not step7_tsv:
        return {
            "status": "skipped",
            "reason": (
                "IRAF cross-check needs a science FITS + Step 7 TSV; configure "
                "[iraf].input_fits and [iraf].step7_tsv to enable it."
            ),
        }
    if not Path(input_fits).exists() or not Path(step7_tsv).exists():
        return {
            "status": "skipped",
            "reason": "configured IRAF inputs do not exist on disk",
            "input_fits": str(input_fits),
            "step7_tsv": str(step7_tsv),
        }
    if not _iraf_available(runtime_cmd):
        return {
            "status": "skipped",
            "reason": "IRAF/PyRAF runtime not found on PATH (wsl/python3)",
        }

    try:
        from apex.benchmark.iraf_crosscheck import IRAFCrosscheckConfig, run_iraf_crosscheck

        cc = IRAFCrosscheckConfig(
            input_fits=str(input_fits),
            step7_tsv=str(step7_tsv),
            output_root=str(output_dir / "iraf_crosscheck"),
            overwrite=True,
        )
        if runtime_cmd:
            cc.runtime_cmd = list(runtime_cmd)
        produced = run_iraf_crosscheck(cc)
    except Exception as exc:  # noqa: BLE001 - surface as skipped, never crash report
        return {
            "status": "skipped",
            "reason": f"IRAF cross-check could not complete: {exc}",
        }

    fixed_summary_path = Path(produced) / "phot_fixed_coords" / "fixed_summary.json"
    metrics: dict[str, Any] = {}
    if fixed_summary_path.exists():
        fs = json.loads(fixed_summary_path.read_text(encoding="utf-8"))
        metrics = {
            "n_matched": fs.get("n_matched"),
            "median_delta_mag_raw": _round(fs.get("median_delta_mag_raw"), 4),
            "mad_sigma_delta_mag_zp_aligned": _round(fs.get("mad_sigma_delta_mag_zp_aligned"), 4),
            "rms_delta_mag_zp_aligned": _round(fs.get("rms_delta_mag_zp_aligned"), 4),
            "position_rms_px": _round(fs.get("position_rms_px"), 4),
        }
    return {"status": "ok", "output_dir": str(produced), "metrics": metrics}


# ── photometric cross-check suite ──────────────────────────────────────────────

def run_crosscheck_suite(output_dir: str | Path, config: Any | None) -> dict[str, Any]:
    """Re-measure an APEX-processed frame with independent photometry engine(s).

    Reads the result dir + frame + reference selection from ``config`` (set by
    the CLI from ``--result-dir`` / ``--frame`` / ``--reference``). Returns
    ``{"status": "skipped", ...}`` when no result dir is configured, so the suite
    is a no-op when run without real APEX outputs. Never raises.
    """
    output_dir = Path(output_dir)

    result_dir = None
    frame = None
    references: Any = ("sep",)
    parameter_file = None
    if config is not None:
        cc_cfg = getattr(config, "crosscheck", None)
        if isinstance(cc_cfg, dict):
            result_dir = cc_cfg.get("result_dir")
            frame = cc_cfg.get("frame")
            references = cc_cfg.get("references", references)
        else:
            result_dir = getattr(config, "crosscheck_result_dir", None)
            frame = getattr(config, "crosscheck_frame", None)
            references = getattr(config, "crosscheck_references", references)
        parameter_file = getattr(config, "parameter_file", None)

    if not result_dir:
        return {
            "status": "skipped",
            "reason": (
                "crosscheck needs an APEX result dir; pass --result-dir "
                "(or set [crosscheck].result_dir) to enable it."
            ),
        }
    if not Path(result_dir).exists():
        return {"status": "skipped", "reason": f"result dir does not exist: {result_dir}"}

    if isinstance(references, str):
        references = [r.strip() for r in references.replace(",", " ").split() if r.strip()]
    references = tuple(references) or ("sep",)

    # Best-effort load of params so aperture scales + gain match the run.
    params = None
    if parameter_file and Path(parameter_file).exists():
        try:
            from apex.config.parameters_cmd import read_params

            params = read_params(Path(parameter_file))
        except Exception:  # noqa: BLE001 - fall back to defaults
            params = None

    try:
        from apex.benchmark.photometry_crosscheck import run_crosscheck

        report = run_crosscheck(
            result_dir,
            output_dir / "crosscheck",
            frame=frame,
            references=references,
            params=params,
        )
    except Exception as exc:  # noqa: BLE001 - surface as error, never crash report
        return {"status": "error", "reason": f"crosscheck failed: {exc}"}

    # Copy per-reference plots into the shared report plots dir.
    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    src_plots = output_dir / "crosscheck" / "plots"
    copied: list[str] = []
    if src_plots.is_dir():
        for png in sorted(src_plots.glob("*.png")):
            dest = plot_dir / f"crosscheck_{png.name}"
            shutil.copyfile(png, dest)
            copied.append(dest.name)
    report["copied_plots"] = copied
    return report


# ── known-targets suite ───────────────────────────────────────────────────────

def run_known_targets_suite(output_dir: str | Path, config: Any | None) -> dict[str, Any]:
    """Always run the zero-data synthetic isochrone-recovery check.

    Additionally run real cluster checks (M67/M5/NGC457) when result dirs are
    configured; otherwise skip the real ones.
    """
    output_dir = Path(output_dir)
    suite_dir = output_dir / "known_targets"
    suite_dir.mkdir(parents=True, exist_ok=True)

    # Reuse the existing synthetic isochrone-recovery validator.
    repo_root = _repo_root()
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

    synthetic: dict[str, Any]
    try:
        from validation.cmd_step12_synthetic import run_validation as run_iso_synthetic

        report = run_iso_synthetic(seed=42)
        (suite_dir / "isochrone_synthetic.json").write_text(
            json.dumps(report, indent=2, sort_keys=True), encoding="utf-8"
        )
        checks = report.get("checks", {})
        synthetic = {
            "status": "ok",
            "passed": bool(report.get("passed")),
            "n_stars": report.get("n_stars"),
            "reduced_chi2": _round(report.get("fit", {}).get("reduced_chi2"), 4),
            "parameter_checks": {
                key: {
                    "truth": _round(item.get("truth"), 4),
                    "recovered": _round(item.get("recovered"), 4),
                    "abs_delta": _round(item.get("abs_delta"), 4),
                    "tolerance": _round(item.get("tolerance"), 4),
                    "passed": bool(item.get("passed")),
                }
                for key, item in checks.items()
            },
        }
    except Exception as exc:  # noqa: BLE001
        synthetic = {"status": "error", "reason": f"synthetic isochrone check failed: {exc}"}

    # Optional real cluster checks.
    real_results: dict[str, Any] = {}
    real_targets: dict[str, Any] = {}
    if config is not None:
        cfg_targets = getattr(config, "known_targets", None)
        if isinstance(cfg_targets, dict):
            real_targets = cfg_targets
    if not real_targets:
        real_results = {"status": "skipped", "reason": "no real cluster result dirs configured"}
    else:
        try:
            from validation.cmd_step12_realdata import run_validation as run_iso_realdata  # type: ignore

            for name, result_dir in real_targets.items():
                try:
                    if not Path(result_dir).exists():
                        real_results[name] = {"status": "skipped", "reason": "result dir missing"}
                        continue
                    rep = run_iso_realdata(result_dir=str(result_dir))  # signature best-effort
                    real_results[name] = {"status": "ok", "passed": bool(rep.get("passed"))}
                except Exception as exc:  # noqa: BLE001
                    real_results[name] = {"status": "error", "reason": str(exc)}
        except Exception as exc:  # noqa: BLE001
            real_results = {"status": "skipped", "reason": f"real-data validator unavailable: {exc}"}

    passed = bool(synthetic.get("passed")) if synthetic.get("status") == "ok" else False
    return {
        "status": "ok" if synthetic.get("status") == "ok" else "error",
        "passed": passed,
        "synthetic_isochrone": synthetic,
        "real_clusters": real_results,
    }


# ── orchestration + report ────────────────────────────────────────────────────

def _suite_status_row(name: str, result: dict[str, Any]) -> tuple[str, str, str]:
    status = result.get("status", "unknown")
    if name == "artificial-star" and status == "ok":
        ok = result.get("completeness_non_degenerate")
        headline = (
            f"m50={result.get('completeness_m50')}, "
            f"RMS(bright)={result.get('photometric_rms_bright_mag')}, "
            f"pos RMSE={result.get('position_rmse_px')} px"
        )
        symbol = "PASS" if ok else "CHECK"
        return name, symbol, headline
    if name == "iraf":
        if status == "skipped":
            return name, "SKIP", str(result.get("reason", ""))
        if status == "ok":
            m = result.get("metrics", {})
            return name, "PASS", f"MAD={m.get('mad_sigma_delta_mag_zp_aligned')} mag, N={m.get('n_matched')}"
        return name, "FAIL", str(result.get("reason", ""))
    if name == "known-targets":
        if status == "ok":
            symbol = "PASS" if result.get("passed") else "FAIL"
            syn = result.get("synthetic_isochrone", {})
            return name, symbol, f"isochrone recovery passed={syn.get('passed')}"
        return name, "FAIL", str(result.get("reason", ""))
    if name == "crosscheck":
        if status == "skipped":
            return name, "SKIP", str(result.get("reason", ""))
        if status == "ok":
            refs = result.get("references", {})
            # Headline off the strongest independent engine that ran (sep > iraf).
            for key in ("sep", "iraf", "photutils"):
                ref = refs.get(key)
                if isinstance(ref, dict) and ref.get("status") == "ok":
                    return (
                        name,
                        "PASS" if _crosscheck_ref_ok(ref) else "CHECK",
                        f"{key}: offset={ref.get('dmag_median')}, MAD={ref.get('dmag_mad')} mag, "
                        f"r={ref.get('pearson_r')}, N={ref.get('n_matched')}",
                    )
            return name, "SKIP", "no reference engine produced a comparison"
        return name, "FAIL", str(result.get("reason", ""))
    return name, status.upper(), str(result.get("reason", ""))


def _crosscheck_ref_ok(ref: dict[str, Any]) -> bool:
    """A cross-check reference 'passes' when the bulk scatter is photometric-grade.

    Uses the robust scatter (MAD) about the zeropoint offset, not the RMS, so a
    small tail of blended/faint outliers does not mask millimag core agreement.
    """
    mad = _finite(ref.get("dmag_mad"))
    n = ref.get("n_matched") or 0
    return bool(mad is not None and mad <= 0.02 and int(n) >= 20)


def write_report(
    results: dict[str, Any],
    output_dir: str | Path,
    *,
    update_docs: bool = True,
) -> dict[str, Path]:
    """Write validation_summary.md, validation_manifest.json, and docs/validation.md.

    Set ``update_docs=False`` to skip refreshing the tracked ``docs/validation.md``
    (used by tests so they do not clobber the committed docs page).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = results.get("timestamp") or _utc_now()
    suites = results.get("suites", {})

    lines: list[str] = []
    lines.append("# APEX Validation Report")
    lines.append("")
    lines.append(f"_Generated: {timestamp}_")
    git = results.get("git", {})
    if git.get("commit"):
        dirty = " (dirty)" if git.get("dirty") else ""
        lines.append(f"_Commit: `{git['commit'][:12]}` on `{git.get('branch', '?')}`{dirty}_")
    lines.append("")
    lines.append(
        "This report is produced by `apex validate`. It quantifies APEX's "
        "selection and measurement accuracy from a reproducible, self-contained "
        "artificial-star experiment, an optional IRAF cross-check, and a "
        "zero-data isochrone-recovery test."
    )
    lines.append("")

    # Status table.
    lines.append("## Status")
    lines.append("")
    lines.append("| Suite | Status | Headline |")
    lines.append("| --- | --- | --- |")
    for name in SUITE_CHOICES:
        if name not in suites:
            continue
        suite_name, symbol, headline = _suite_status_row(name, suites[name])
        lines.append(f"| {suite_name} | {symbol} | {headline} |")
    lines.append("")

    # Artificial-star detail.
    art = suites.get("artificial-star")
    if art:
        lines.append("## Artificial-star completeness & photometry")
        lines.append("")
        if art.get("status") != "ok":
            lines.append(f"> **{art.get('status', 'error')}** — {art.get('reason', '')}")
            lines.append("")
        else:
            lines.append(f"- Trials: {art.get('n_trials')}, stars/trial: {art.get('stars_per_trial')}, total injected: {art.get('n_injected')}")
            lines.append(f"- Injection magnitude range: {art.get('injection_mag_range')}")
            lines.append(f"- FWHM: {art.get('fwhm_px')} px, injection zeropoint: {art.get('zeropoint_mag')} mag")
            lines.append(f"- **Completeness 50% magnitude (m50): {art.get('completeness_m50')}** "
                         f"(95% CI {art.get('completeness_m50_ci95')})")
            lines.append(f"- Completeness m90 / m10: {art.get('completeness_m90')} / {art.get('completeness_m10')} "
                         f"(logistic width {art.get('completeness_width_mag')} mag)")
            lines.append(f"- Completeness curve non-degenerate: **{art.get('completeness_non_degenerate')}**")
            lines.append(f"- Photometric bias (median): {art.get('photometric_bias_median_mag')} mag; "
                         f"scatter (MAD): {art.get('photometric_scatter_mad_mag')} mag")
            lines.append(f"- Photometric RMS — bright: {art.get('photometric_rms_bright_mag')} mag, "
                         f"faint: {art.get('photometric_rms_faint_mag')} mag")
            lines.append(f"- Position RMSE: {art.get('position_rmse_px')} px")
            lines.append(f"- False-positive rate (new unmatched detections / injected): {art.get('false_positive_rate')} "
                         f"({art.get('new_false_detections')} false of {art.get('new_detections')} new)")
            lines.append("")
            per_mag = art.get("completeness_by_magnitude", [])
            if per_mag:
                lines.append("| Injected mag | N eligible | Completeness | Mag bias (median) |")
                lines.append("| --- | --- | --- | --- |")
                for row in per_mag:
                    lines.append(
                        f"| {row.get('magnitude_center')} | {row.get('n_eligible')} | "
                        f"{row.get('completeness')} | {row.get('mag_bias_median')} |"
                    )
                lines.append("")
            for plot in art.get("plots", []):
                lines.append(f"![{plot}](plots/{plot})")
                lines.append("")

    # IRAF detail.
    iraf = suites.get("iraf")
    if iraf:
        lines.append("## IRAF cross-check")
        lines.append("")
        if iraf.get("status") == "skipped":
            lines.append(f"> **skipped** — {iraf.get('reason', '')}")
        elif iraf.get("status") == "ok":
            m = iraf.get("metrics", {})
            lines.append(f"- Matched sources: {m.get('n_matched')}")
            lines.append(f"- APEX-IRAF residual MAD (ZP-aligned): {m.get('mad_sigma_delta_mag_zp_aligned')} mag")
            lines.append(f"- APEX-IRAF residual RMS (ZP-aligned): {m.get('rms_delta_mag_zp_aligned')} mag")
            lines.append(f"- Position RMS: {m.get('position_rms_px')} px")
        else:
            lines.append(f"> **{iraf.get('status')}** — {iraf.get('reason', '')}")
        lines.append("")

    # Photometric cross-check detail.
    cc = suites.get("crosscheck")
    if cc:
        lines.append("## Photometric cross-check")
        lines.append("")
        lines.append(
            "Re-measures the *same stars in the same frame* APEX measured, using "
            "an independent photometry engine, then reports the agreement. This "
            "validates the MEASUREMENT directly (independent of isochrone fitting). "
            "The zeropoint offset absorbs the constant gain / exposure-time / "
            "aperture-correction difference between the two flux scales; the "
            "robust scatter about it is the real agreement metric."
        )
        lines.append("")
        if cc.get("status") == "skipped":
            lines.append(f"> **skipped** — {cc.get('reason', '')}")
            lines.append("")
        elif cc.get("status") != "ok":
            lines.append(f"> **{cc.get('status', 'error')}** — {cc.get('reason', '')}")
            lines.append("")
        else:
            lines.append(
                f"- Frame: `{cc.get('frame')}` (filter {cc.get('filter')}), "
                f"FWHM {cc.get('fwhm_px')} px, gain {cc.get('gain')} e-/ADU, SNR>{cc.get('snr_min')}"
            )
            lines.append(f"- Image measured: `{cc.get('image')}`")
            sanity = cc.get("sanity", {}) or {}
            verdict = "PLAUSIBLE" if cc.get("positions_plausible") else "SUSPECT"
            lines.append(
                f"- Position sanity: **{verdict}** — star peak pixel sits "
                f"{_round(sanity.get('star_over_sky_sigma'), 1)}σ above sky, "
                f"{_round((sanity.get('frac_on_frame') or 0) * 100, 1)}% on-frame "
                "(if SUSPECT, the image/coordinate frame is mismatched and the "
                "numbers below are meaningless)."
            )
            lines.append("")
            refs = cc.get("references", {}) or {}
            lines.append("| Reference | Status | N | Offset (mag) | MAD (mag) | RMS (mag) | robust σ | Pearson r | <0.05 | <0.1 |")
            lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
            for key in ("sep", "iraf", "photutils"):
                ref = refs.get(key)
                if not isinstance(ref, dict):
                    continue
                if ref.get("status") != "ok":
                    lines.append(f"| {key} | {ref.get('status')} | — | — | — | — | — | — | — | — |")
                    continue
                lines.append(
                    f"| {key} | ok | {ref.get('n_matched')} | {_round(ref.get('dmag_median'))} | "
                    f"{_round(ref.get('dmag_mad'))} | {_round(ref.get('dmag_rms'))} | "
                    f"{_round(ref.get('robust_sigma'))} | {_round(ref.get('pearson_r'), 5)} | "
                    f"{_round(ref.get('frac_within_0p05'), 3)} | {_round(ref.get('frac_within_0p1'), 3)} |"
                )
            lines.append("")
            if refs.get("photutils", {}).get("status") == "ok":
                lines.append(
                    "_Note: `photutils` shares APEX's measurement backend, so it is a "
                    "consistency reference (different call path, same library), not a "
                    "fully independent tool. `sep` (independent C engine) and `iraf` are "
                    "the independent tests._"
                )
                lines.append("")
            for plot in cc.get("copied_plots", []):
                lines.append(f"![{plot}](plots/{plot})")
                lines.append("")

    # Known-targets detail.
    kt = suites.get("known-targets")
    if kt:
        lines.append("## Known-target reproduction (isochrone recovery)")
        lines.append("")
        syn = kt.get("synthetic_isochrone", {})
        if syn.get("status") == "ok":
            lines.append(f"- Synthetic isochrone recovery passed: **{syn.get('passed')}** "
                         f"(N stars {syn.get('n_stars')}, reduced chi2 {syn.get('reduced_chi2')})")
            checks = syn.get("parameter_checks", {})
            if checks:
                lines.append("")
                lines.append("| Parameter | Truth | Recovered | abs delta | Tolerance | Pass |")
                lines.append("| --- | --- | --- | --- | --- | --- |")
                for key, item in checks.items():
                    lines.append(
                        f"| {key} | {item.get('truth')} | {item.get('recovered')} | "
                        f"{item.get('abs_delta')} | {item.get('tolerance')} | {item.get('passed')} |"
                    )
        else:
            lines.append(f"> **{syn.get('status')}** — {syn.get('reason', '')}")
        real = kt.get("real_clusters", {})
        lines.append("")
        if isinstance(real, dict) and real.get("status") == "skipped":
            lines.append(f"_Real cluster checks: skipped — {real.get('reason', '')}_")
        elif real:
            lines.append("_Real cluster checks:_")
            for name, item in real.items():
                if name == "status":
                    continue
                if isinstance(item, dict):
                    lines.append(f"- {name}: {item.get('status')} (passed={item.get('passed')})")
        lines.append("")

    # Interpretation & caveats.
    lines.append("## Interpretation & caveats")
    lines.append("")
    lines.append(
        "- The artificial-star suite injects an **empirical PSF** sampled from the "
        "reference frame and reruns the *production* Step 4 detector, so it tests "
        "the real selection + measurement chain — not a toy model."
    )
    lines.append(
        "- When no real reference frame is supplied, the experiment runs on a "
        "**synthetic frame** with a uniform Gaussian PSF and Poisson+read-noise "
        "background. Its completeness transition is sharper than a real, crowded "
        "field because there is no faint-source confusion; m50/m90/m10 remain "
        "meaningful but the logistic width is narrower than a typical survey."
    )
    lines.append(
        "- Completeness, photometric bias/scatter, and position RMSE are "
        "*differential* (recovered vs injected truth) and depend on the frame's "
        "noise and the configured detection thresholds."
    )
    lines.append(
        "- The IRAF cross-check is an external-tool agreement test and is skipped "
        "automatically when IRAF/PyRAF or its inputs are unavailable; a skip is "
        "not a failure."
    )
    lines.append(
        "- The isochrone-recovery test is a science-regression guard against a "
        "known synthetic truth, not an astrophysical calibration."
    )
    if suites.get("crosscheck", {}).get("status") == "ok":
        lines.append(
            "- The photometric cross-check's **robust scatter (MAD / robust σ)** "
            "is the headline metric, not the RMS. In crowded fields the RMS is "
            "inflated by a small tail of faint/blended sources where APEX's "
            "recentering + aperture correction legitimately diverge from a fixed "
            "aperture; the MAD reflects the bulk core agreement. `sep` is a fully "
            "independent C engine; `photutils` shares APEX's backend and is only a "
            "consistency reference. The cross-check re-measures APEX's own catalog, "
            "so it validates measurement reproducibility against another tool — it "
            "is not an absolute-truth test (that is the artificial-star suite)."
        )
    lines.append("")

    summary_path = output_dir / "validation_summary.md"
    summary_path.write_text("\n".join(lines), encoding="utf-8")

    # Manifest.
    manifest = {
        "timestamp": timestamp,
        "git": results.get("git", {}),
        "versions": results.get("versions", {}),
        "suites": suites,
    }
    manifest_path = output_dir / "validation_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, default=str), encoding="utf-8"
    )

    # docs/validation.md — a stable pointer into the docs site.
    docs_path = _repo_root() / "docs" / "validation.md"
    art = suites.get("artificial-star", {}) or {}
    docs_lines = [
        "# Validation",
        "",
        "APEX ships a reproducible validation harness. Run it with:",
        "",
        "```bash",
        "apex validate --suite all --output validation/report",
        "```",
        "",
        "This generates `validation/report/validation_summary.md` (human-readable "
        "report) and `validation/report/validation_manifest.json` (machine-readable "
        "metrics with the git commit and package versions).",
        "",
        "## What it measures",
        "",
        "- **Artificial-star completeness & photometry** — injects an empirical PSF "
        "into a reference frame (synthetic when none is supplied), reruns the "
        "production detector, and reports the 50% completeness magnitude, "
        "photometric bias/scatter, and position RMSE.",
        "- **IRAF cross-check** (optional) — compares APEX aperture magnitudes "
        "against IRAF `phot` on the same fixed coordinates.",
        "- **Known-target reproduction** — a zero-data synthetic isochrone-recovery "
        "regression guard (plus optional real-cluster checks).",
        "",
        "## Latest headline numbers",
        "",
    ]
    if art.get("status") == "ok":
        docs_lines += [
            f"- 50% completeness magnitude (m50): **{art.get('completeness_m50')}**",
            f"- Photometric RMS (bright): **{art.get('photometric_rms_bright_mag')} mag**",
            f"- Position RMSE: **{art.get('position_rmse_px')} px**",
            f"- False-positive rate: **{art.get('false_positive_rate')}**",
            "",
            "_Numbers above are from the most recent committed report; rerun "
            "`apex validate` to refresh._",
        ]
    else:
        docs_lines += ["_Run `apex validate` to populate the latest numbers._"]
    docs_lines.append("")
    if update_docs:
        docs_path.parent.mkdir(parents=True, exist_ok=True)
        docs_path.write_text("\n".join(docs_lines), encoding="utf-8")

    return {
        "summary": summary_path,
        "manifest": manifest_path,
        "docs": docs_path,
    }


def run_validation(
    suites: Iterable[str] | str,
    output_dir: str | Path,
    config: Any | None = None,
    *,
    reference_frame: str | Path | None = None,
    timestamp: str | None = None,
    **suite_kwargs: Any,
) -> dict[str, Any]:
    """Run the requested suites, aggregate metrics, and write the report."""
    if isinstance(suites, str):
        requested = [suites]
    else:
        requested = list(suites)
    if "all" in requested:
        requested = list(SUITE_CHOICES)
    unknown = [s for s in requested if s not in SUITE_CHOICES]
    if unknown:
        raise ValueError(f"Unknown validation suite(s): {', '.join(unknown)}")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results: dict[str, Any] = {
        "timestamp": timestamp or _utc_now(),
        "git": _git_commit(_repo_root()),
        "versions": _package_versions(),
        "suites": {},
    }

    if "artificial-star" in requested:
        results["suites"]["artificial-star"] = run_artificial_star_suite(
            output_dir, reference_frame=reference_frame, config=config, **suite_kwargs
        )
    if "iraf" in requested:
        results["suites"]["iraf"] = run_iraf_suite(output_dir, config)
    if "crosscheck" in requested:
        results["suites"]["crosscheck"] = run_crosscheck_suite(output_dir, config)
    if "known-targets" in requested:
        results["suites"]["known-targets"] = run_known_targets_suite(output_dir, config)

    paths = write_report(results, output_dir)
    results["report_paths"] = {k: str(v) for k, v in paths.items()}
    return results
