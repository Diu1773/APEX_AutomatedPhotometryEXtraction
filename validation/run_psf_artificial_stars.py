"""Run a resumable real-image Step 4 + Step 8 artificial-star benchmark."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from apex.benchmark.psf_artificial_stars import (  # noqa: E402
    DEFAULT_TARGET_SNRS,
    add_forced_truth_to_step7,
    aggregate_recovery_metrics,
    apply_recovery_quality_policy,
    inject_flux_catalog,
    match_injections_to_products,
    measure_preinjection_psf_residual,
    optimal_psf_flux_for_snr,
    oversampled_epsf_to_native_kernel,
    psf_noise_equivalent_area,
    sample_stratified_injections,
    write_benchmark_outputs,
)


def _source_name(path: Path) -> str:
    return path.name


def _detect_path(result_dir: Path, filename: str) -> Path | None:
    for path in (result_dir / "step4_detection" / f"detect_{filename}.csv", result_dir / "cache" / f"detect_{filename}.csv"):
        if path.exists():
            return path
    return next(iter(sorted((result_dir / "step4_detection").glob("detect_*.csv"))), None) if (result_dir / "step4_detection").exists() else None


def _psf_path(result_dir: Path, filename: str) -> Path | None:
    exact = result_dir / "cmd_psf" / f"photometry_{filename}.tsv"
    if exact.exists():
        return exact
    paths = sorted((result_dir / "cmd_psf").glob("photometry_*.tsv")) if (result_dir / "cmd_psf").exists() else []
    return paths[0] if len(paths) == 1 else None


def _baseline_epsf(result_dir: Path, explicit: Path | None) -> Path:
    if explicit is not None:
        if not explicit.exists():
            raise FileNotFoundError(explicit)
        return explicit
    paths = sorted((result_dir / "cmd_psf").glob("epsf_model_*.fits"))
    if not paths:
        raise FileNotFoundError("no baseline cmd_psf/epsf_model_*.fits; pass --epsf-model")
    return paths[0]


def _baseline_step7(baseline_dir: Path, filename: str, explicit: Path | None) -> Path | None:
    if explicit is not None:
        return explicit if explicit.exists() else None
    step7 = baseline_dir / "step7_forced_phot"
    exact = step7 / f"photometry_{filename}.tsv"
    if exact.exists():
        return exact
    paths = sorted(step7.glob("photometry_*.tsv")) if step7.exists() else []
    return paths[0] if len(paths) == 1 else None


def fit_frame_moffat(
    image: np.ndarray,
    step7_table: Path,
    fwhm_px: float,
    *,
    n_stars: int = 40,
    min_snr: float = 50.0,
) -> tuple[float, float, int]:
    """Fit one Moffat profile to the frame's bright isolated real stars.

    The injection experiment needs a truth shape that neither engine built.
    The baseline ePSF is APEX's own model, so implanting it hands APEX a truth
    drawn from its model family — the completeness numbers inherit that
    favour. Fitting an analytic Moffat directly to the frame's real stars
    (astropy fitter, no Step 8 or DAOPHOT machinery involved) gives a truth
    that approximates the same real PSF while belonging to neither engine.

    Returns (gamma, alpha, n_used) — the median of per-star fits. Stars are
    the brightest isolated ones: SNR >= min_snr, nearest catalogue neighbour
    beyond 6 FWHM, whole stamp inside the frame, no saturation flag.
    """
    from astropy.modeling import fitting, functional_models
    from scipy.spatial import cKDTree

    table = pd.read_csv(step7_table, sep="\t")
    for column in ("x", "y", "snr"):
        table[column] = pd.to_numeric(table.get(column), errors="coerce")
    ok = np.isfinite(table["x"]) & np.isfinite(table["y"]) & np.isfinite(table["snr"])
    if "is_saturated" in table.columns:
        ok &= ~table["is_saturated"].astype(str).str.lower().isin({"true", "1"})
    table = table[ok]

    positions = table[["x", "y"]].to_numpy(float)
    nn, _ = cKDTree(positions).query(positions, k=2)
    isolated = nn[:, 1] >= 6.0 * fwhm_px
    candidates = table[isolated & (table["snr"] >= min_snr)]
    candidates = candidates.sort_values("snr", ascending=False).head(int(n_stars))

    half = max(3, int(round(2.0 * fwhm_px)))
    gamma_init = fwhm_px / (2.0 * float(np.sqrt(2.0 ** (1.0 / 2.5) - 1.0)))
    fitter = fitting.TRFLSQFitter()
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1].astype(float)

    gammas: list[float] = []
    alphas: list[float] = []
    for _, star in candidates.iterrows():
        cx, cy = int(round(star["x"])), int(round(star["y"]))
        if not (half <= cx < image.shape[1] - half and half <= cy < image.shape[0] - half):
            continue
        stamp = np.asarray(
            image[cy - half:cy + half + 1, cx - half:cx + half + 1], dtype=float)
        edge = np.concatenate([stamp[0], stamp[-1], stamp[1:-1, 0], stamp[1:-1, -1]])
        stamp = stamp - float(np.median(edge))
        peak = float(stamp.max())
        if not np.isfinite(peak) or peak <= 0:
            continue
        model = functional_models.Moffat2D(
            amplitude=peak, x_0=0.0, y_0=0.0, gamma=gamma_init, alpha=2.5,
            bounds={"gamma": (0.5, 12.0 * fwhm_px), "alpha": (1.2, 60.0),
                    "x_0": (-2.0, 2.0), "y_0": (-2.0, 2.0)},
        )
        try:
            fit = fitter(model, xx, yy, stamp, maxiter=300)
        except Exception:
            continue
        gamma, alpha = float(fit.gamma.value), float(fit.alpha.value)
        # Judge the profile, not the parameters. gamma and alpha are degenerate
        # for a PSF without measurable wings: as alpha grows the Moffat becomes
        # a Gaussian and gamma grows with it, so both run to whatever bound is
        # set while the *shape* stays put. Rejecting on `alpha < 9.9` therefore
        # threw away every star on the LCO 0.4 m frame — 39 of 40 — even though
        # each fit reproduced the measured FWHM to 2 % (2026-08-14). What the
        # injection needs is the shape, so that is what is checked; a fit that
        # wandered off-centre or that does not reproduce the frame's FWHM is a
        # neighbour or a defect and still goes.
        implied_fwhm = 2.0 * gamma * float(np.sqrt(2.0 ** (1.0 / alpha) - 1.0))
        if (np.hypot(float(fit.x_0.value), float(fit.y_0.value)) > 1.5
                or not np.isfinite(implied_fwhm)
                or not (0.70 * fwhm_px < implied_fwhm < 1.40 * fwhm_px)):
            continue
        gammas.append(gamma)
        alphas.append(alpha)

    if len(gammas) < 8:
        raise RuntimeError(
            f"Moffat 적합에 쓸 별이 부족하다: {len(gammas)}개 (최소 8)")
    return float(np.median(gammas)), float(np.median(alphas)), len(gammas)


def make_moffat_sampler(gamma: float, alpha: float, size: int):
    """Return a kernel sampler evaluating the fitted Moffat analytically.

    Same contract as the ePSF sampler: odd native grid, sub-pixel phases are
    the source offset from the nearest pixel centre, kernel sums to one so the
    requested flux is the flux actually implanted within the window.
    """
    from astropy.modeling import functional_models

    if size % 2 == 0:
        size += 1
    half = size // 2
    offsets = np.arange(-half, half + 1, dtype=float)
    yy, xx = np.meshgrid(offsets, offsets, indexing="ij")
    model = functional_models.Moffat2D(
        amplitude=1.0, x_0=0.0, y_0=0.0, gamma=float(gamma), alpha=float(alpha))

    def sampler(phase_x: float, phase_y: float) -> np.ndarray:
        values = model(xx - float(phase_x), yy - float(phase_y))
        values = np.where(np.isfinite(values), np.maximum(values, 0.0), 0.0)
        total = float(values.sum())
        if not np.isfinite(total) or total <= 0:
            raise ValueError("Moffat kernel has no positive finite flux")
        return (values / total).astype(np.float64)

    return sampler


def _baseline_residual(result_dir: Path, filename: str) -> Path | None:
    psf_dir = result_dir / "cmd_psf"
    for prefix in ("starsub_", "residual_"):
        exact = psf_dir / f"{prefix}{filename}"
        if exact.exists():
            return exact
    return None


def _copy_baseline_seed(
    baseline_dir: Path,
    trial_result: Path,
    filename: str,
    explicit: Path | None,
) -> Path | None:
    source = _baseline_step7(baseline_dir, filename, explicit)
    if source is None:
        return None
    target_dir = trial_result / "step7_forced_phot"
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / f"photometry_{filename}.tsv"
    shutil.copy2(source, target)
    index = source.parent / "photometry_index.csv"
    if index.exists():
        shutil.copy2(index, target_dir / index.name)
    return target


def _copy_baseline_detection(
    baseline_dir: Path,
    trial_result: Path,
    filename: str,
) -> Path:
    source = _detect_path(baseline_dir, filename)
    if source is None:
        raise FileNotFoundError("baseline Step4 detection table is unavailable")
    target_dir = trial_result / "step4_detection"
    target_dir.mkdir(parents=True, exist_ok=True)
    target = target_dir / f"detect_{filename}.csv"
    shutil.copy2(source, target)
    source_summary = source.with_suffix(".json")
    if source_summary.exists():
        shutil.copy2(source_summary, target.with_suffix(".json"))
    return target


def _estimate_background_rms(image: np.ndarray) -> float:
    array = np.asarray(image)
    stride = max(1, int(np.ceil(np.sqrt(array.size / 1_000_000))))
    sample = np.asarray(array[::stride, ::stride], dtype=float)
    _, _, std = sigma_clipped_stats(sample, sigma=3.0, maxiters=5)
    value = float(std)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("could not estimate a positive background RMS")
    return value


def _read_background_rms(
    detection_path: Path | None,
    override: float | None,
    image: np.ndarray,
) -> float:
    if override is not None:
        value = float(override)
        if np.isfinite(value) and value > 0:
            return value
        raise ValueError("--background-rms-adu must be positive")
    if detection_path is not None:
        summary_path = detection_path.with_suffix(".json")
        if summary_path.exists():
            try:
                payload = json.loads(summary_path.read_text(encoding="utf-8"))
                value = float(payload.get("bkg_rms", np.nan))
                if np.isfinite(value) and value > 0:
                    return value
            except (OSError, TypeError, ValueError, json.JSONDecodeError):
                pass
    return _estimate_background_rms(image)


def _read_fwhm(detection_path: Path | None, override: float | None) -> float:
    if override is not None:
        value = float(override)
        if np.isfinite(value) and value > 0:
            return value
        raise ValueError("--fwhm-px must be positive")
    if detection_path is not None:
        summary_path = detection_path.with_suffix(".json")
        if summary_path.exists():
            try:
                payload = json.loads(summary_path.read_text(encoding="utf-8"))
                value = float(payload.get("fwhm_px", np.nan))
                if np.isfinite(value) and value > 0:
                    return value
            except (OSError, TypeError, ValueError, json.JSONDecodeError):
                pass
        df = pd.read_csv(
            detection_path,
            usecols=lambda name: name in {
                "fwhm_px",
                "fwhm_status",
                "fwhm_ratio_to_frame",
            },
        )
        if "fwhm_px" in df:
            values_series = pd.to_numeric(df["fwhm_px"], errors="coerce")
            if "fwhm_status" in df:
                status = df["fwhm_status"].fillna("").astype(str).str.lower()
                bounded = values_series[status == "ok"].to_numpy(float)
                bounded = bounded[np.isfinite(bounded) & (bounded > 0)]
                if len(bounded):
                    return float(np.median(bounded))
            if "fwhm_ratio_to_frame" in df:
                ratio = pd.to_numeric(
                    df["fwhm_ratio_to_frame"], errors="coerce"
                ).to_numpy(float)
                values = values_series.to_numpy(float)
                inferred = values / ratio
                inferred = inferred[np.isfinite(inferred) & (inferred > 0)]
                if len(inferred):
                    return float(np.median(inferred))
            values = values_series.to_numpy(float)
            values = values[np.isfinite(values) & (values > 0)]
            if len(values):
                return float(np.median(values))
    raise ValueError("FWHM is unavailable; pass --fwhm-px")


def _parameter_gain(parameter_file: Path) -> float:
    try:
        from apex.config.parameters_cmd import Parameters
        params = Parameters(parameter_file)
        value = float(getattr(params.P, "gain_e_per_adu"))
        if value > 0:
            return value
    except Exception:
        pass
    return 1.0


def _build_truth(
    positions: pd.DataFrame,
    target_snrs: tuple[float, ...],
    gain: float,
    background_rms: float,
    nea: float,
    rng: np.random.Generator,
    *,
    kernel_sampler=None,
) -> pd.DataFrame:
    truth = positions.copy()
    values = np.resize(np.asarray(target_snrs, dtype=float), len(truth))
    rng.shuffle(values)
    if kernel_sampler is None:
        neas = np.full(len(truth), float(nea), dtype=float)
    else:
        neas = np.asarray([
            psf_noise_equivalent_area(
                kernel_sampler(
                    float(row.x_true) - round(float(row.x_true)),
                    float(row.y_true) - round(float(row.y_true)),
                )
            )
            for row in truth.itertuples(index=False)
        ], dtype=float)
    fluxes = [
        optimal_psf_flux_for_snr(
            target_snr,
            gain_e_per_adu=gain,
            background_rms_adu=background_rms,
            psf_nea_px=phase_nea,
        )
        for target_snr, phase_nea in zip(values, neas)
    ]
    truth["target_snr"] = values
    truth["true_flux_e"] = [item[0] for item in fluxes]
    truth["true_flux_adu"] = [item[1] for item in fluxes]
    truth["gain_e_per_adu"] = gain
    truth["background_rms_adu"] = background_rms
    truth["psf_nea_px"] = neas
    return truth


def _trial_command(args: argparse.Namespace, trial_data: Path, trial_result: Path, filename: str) -> list[str]:
    command = [
        sys.executable,
        str(REPO / "validation" / "run_real_gui_psf.py"),
        "--data-dir", str(trial_data),
        "--result-dir", str(trial_result),
        "--parameter-file", str(args.parameter_file),
        "--fit-engine", "apex_iterative",
        "--fitter-max-iter", "8",
        "--fit-shape-fwhm-mult", str(args.fit_shape_fwhm_mult),
        "--fit-window-mode", args.fit_window_mode,
        "--fit-encircled-energy", str(args.fit_encircled_energy),
        "--postfit-qfit-noise-max", str(args.postfit_qfit_noise_max),
        "--residual-passes", "2",
        "--epsf-max-stars", "0",
        "--epsf-contamination", "on",
        "--flux-scale", args.flux_scale,
        "--use-grouper", args.use_grouper,
        "--grouper-max-size", str(args.grouper_max_size),
        "--grouper-radius-fwhm", str(args.grouper_radius_fwhm),
        "--forced-match-radius-fwhm", str(args.forced_match_radius_fwhm),
        "--core-cut", "off",
    ]
    if args.reuse_baseline_detection:
        command.append("--skip-step4")
    if args.center_x is not None and args.center_y is not None:
        command.extend([
            "--core-center-mode", "manual",
            "--core-x", str(args.center_x),
            "--core-y", str(args.center_y),
        ])
    else:
        command.extend(["--core-center-mode", "auto"])
    command.append(filename)
    return command


def _remove_trial_dir(trial_dir: Path, output_dir: Path) -> None:
    target = trial_dir.resolve()
    root = output_dir.resolve()
    if target.parent != root or not target.name.startswith("trial_"):
        raise RuntimeError(f"refusing to remove unsafe trial path: {target}")
    shutil.rmtree(target)


def _run_one_trial(args: argparse.Namespace, trial_number: int, image: np.ndarray, header: fits.Header, source_path: Path, kernel: np.ndarray, kernel_sampler, baseline_residual: np.ndarray | None, real_positions: np.ndarray, fwhm_px: float, gain: float, background_rms: float, nea: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    trial_dir = args.output_dir / f"trial_{trial_number:04d}"
    marker = trial_dir / "complete.json"
    if marker.exists() and not args.force and (trial_dir / "truth.csv").exists() and (trial_dir / "recovery.csv").exists():
        truth = pd.read_csv(trial_dir / "truth.csv")
        if baseline_residual is not None:
            truth = measure_preinjection_psf_residual(
                baseline_residual,
                truth,
                kernel_sampler=kernel_sampler,
            )
            truth.to_csv(trial_dir / "truth.csv", index=False)
        filename = _source_name(source_path)
        detection_path = _detect_path(trial_dir / "result", filename)
        psf_path = _psf_path(trial_dir / "result", filename)
        if detection_path is not None and psf_path is not None:
            recovery = match_injections_to_products(
                truth,
                pd.read_csv(detection_path),
                pd.read_csv(psf_path, sep="\t"),
                radius_px=args.match_radius_px,
            )
            recovery["trial"] = trial_number
        else:
            recovery = apply_recovery_quality_policy(
                pd.read_csv(trial_dir / "recovery.csv")
            )
        recovery.to_csv(trial_dir / "recovery.csv", index=False)
        return truth, recovery
    if args.force and trial_dir.exists():
        _remove_trial_dir(trial_dir, args.output_dir)
    trial_data = trial_dir / "data"
    trial_result = trial_dir / "result"
    trial_data.mkdir(parents=True, exist_ok=True)
    trial_result.mkdir(parents=True, exist_ok=True)
    filename = _source_name(source_path)
    if args.reuse_baseline_detection:
        _copy_baseline_detection(
            args.baseline_result_dir,
            trial_result,
            filename,
        )
    rng = np.random.default_rng(args.seed + trial_number - 1)
    positions = sample_stratified_injections(
        image.shape, real_positions, count=args.injections, fwhm_px=fwhm_px, rng=rng,
        center_xy=(args.center_x, args.center_y) if args.center_x is not None and args.center_y is not None else None,
        pixel_scale_arcsec=args.pixel_scale_arcsec, psf_size=kernel.shape[0],
        min_real_sep_fwhm=args.min_real_sep_fwhm, min_injected_sep_fwhm=args.min_injected_sep_fwhm,
        pair_fraction=args.pair_fraction,
        pair_separations_fwhm=tuple(args.pair_separations_fwhm),
    )
    truth = _build_truth(
        positions,
        tuple(args.target_snr),
        gain,
        background_rms,
        nea,
        rng,
        kernel_sampler=kernel_sampler,
    )
    if baseline_residual is not None:
        truth = measure_preinjection_psf_residual(
            baseline_residual,
            truth,
            kernel_sampler=kernel_sampler,
        )
    injected, _, _, truth = inject_flux_catalog(
        image,
        kernel,
        truth,
        gain_e_per_adu=gain,
        rng=rng,
        kernel_sampler=kernel_sampler,
        return_layers=False,
    )
    trial_source = trial_data / filename
    fits.PrimaryHDU(data=np.asarray(injected, dtype=np.float32), header=header).writeto(trial_source, overwrite=True)
    del injected
    step7_target = _copy_baseline_seed(
        args.baseline_result_dir,
        trial_result,
        filename,
        args.baseline_step7,
    )
    copied_step7_seed = step7_target is not None
    forced_truth_seeded = False
    if not args.blind_only:
        if step7_target is not None:
            step7 = pd.read_csv(step7_target, sep="\t")
        else:
            step7_dir = trial_result / "step7_forced_phot"
            step7_dir.mkdir(parents=True, exist_ok=True)
            step7_target = step7_dir / f"photometry_{filename}.tsv"
            step7 = pd.DataFrame()
        step7 = add_forced_truth_to_step7(step7, truth, filename=filename)
        step7.to_csv(step7_target, sep="\t", index=False)
        forced_truth_seeded = True
    env = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[name] = "1"
    with (trial_dir / "run.log").open("w", encoding="utf-8") as log:
        completed = subprocess.run(_trial_command(args, trial_data, trial_result, filename), cwd=REPO, env=env, stdout=log, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"trial {trial_number} production validation failed; see {trial_dir / 'run.log'}")
    detection_path = _detect_path(trial_result, filename)
    psf_path = _psf_path(trial_result, filename)
    if detection_path is None or psf_path is None:
        raise FileNotFoundError(f"trial {trial_number} did not produce Step4 and Step8 tables")
    detections = pd.read_csv(detection_path)
    psf = pd.read_csv(psf_path, sep="\t")
    recovery = match_injections_to_products(truth, detections, psf, radius_px=args.match_radius_px)
    truth["trial"] = trial_number
    recovery["trial"] = trial_number
    summary = aggregate_recovery_metrics(recovery)
    write_benchmark_outputs(
        truth,
        recovery,
        summary,
        trial_dir,
        metadata={
            "trial": trial_number,
            "seed": args.seed + trial_number - 1,
            "copied_step7_seed": copied_step7_seed,
            "forced_truth_seeded": forced_truth_seeded,
            "reused_baseline_detection": bool(args.reuse_baseline_detection),
        },
    )
    marker.write_text(json.dumps({"status": "complete", "trial": trial_number, "returncode": completed.returncode}, indent=2), encoding="utf-8")
    return truth, recovery


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-fits", type=Path, required=True)
    parser.add_argument("--baseline-result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path, default=REPO / "parameters.toml")
    parser.add_argument("--epsf-model", type=Path)
    parser.add_argument("--baseline-step7", type=Path)
    parser.add_argument("--trials", type=int, default=2)
    parser.add_argument("--injections", type=int, default=25)
    parser.add_argument("--target-snr", type=float, nargs="+", default=list(DEFAULT_TARGET_SNRS))
    parser.add_argument("--seed", type=int, default=20260714)
    parser.add_argument("--gain-e-per-adu", type=float)
    parser.add_argument("--background-rms-adu", type=float)
    parser.add_argument("--fwhm-px", type=float)
    parser.add_argument("--fit-shape-fwhm-mult", type=float, default=2.4)
    parser.add_argument(
        "--fit-window-mode", choices=("auto", "manual"), default="manual"
    )
    parser.add_argument("--fit-encircled-energy", type=float, default=0.90)
    parser.add_argument("--postfit-qfit-noise-max", type=float, default=3.0)
    parser.add_argument("--pixel-scale-arcsec", type=float)
    parser.add_argument("--center-x", type=float)
    parser.add_argument("--center-y", type=float)
    parser.add_argument("--min-real-sep-fwhm", type=float, default=0.75)
    parser.add_argument("--min-injected-sep-fwhm", type=float, default=3.0)
    # A sparse field cannot populate the tight-crowding bins on its own, so
    # without deliberate companions the benchmark quietly becomes an
    # isolated-star test on anything but a globular cluster.
    parser.add_argument("--pair-fraction", type=float, default=0.0,
                        help="fraction of injections placed as close pairs "
                             "(0 keeps the previous field-driven behaviour)")
    parser.add_argument("--pair-separations-fwhm", type=float, nargs="+",
                        default=[0.8, 1.2, 2.0, 3.0],
                        help="companion separations sampled for each pair")
    parser.add_argument("--match-radius-px", type=float, default=1.5)
    parser.add_argument("--use-grouper", choices=("on", "off"), default="off")
    parser.add_argument("--grouper-max-size", type=int, default=3)
    parser.add_argument("--grouper-radius-fwhm", type=float, default=1.5)
    parser.add_argument("--forced-match-radius-fwhm", type=float, default=1.25)
    parser.add_argument(
        "--flux-scale",
        choices=("on", "off"),
        default="off",
        help="Apply the optional aperture-based PSF flux scale in production Step 8.",
    )
    parser.add_argument(
        "--blind-only",
        action="store_true",
        help="Do not append injected truth positions to the copied Step 7 seed table.",
    )
    parser.add_argument(
        "--reuse-baseline-detection",
        action="store_true",
        help=(
            "Copy the baseline Step4 table and run only production Step8. "
            "Use this for PSF-fitting A/B tests, not Step4 completeness tests."
        ),
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--inject-kernel", choices=("epsf", "moffat"), default="epsf",
        help=(
            "인공별의 참 모양. epsf = 기준선 ePSF (APEX 모형 계열이라 APEX 에 "
            "유리한 상한). moffat = 프레임의 밝고 고립된 실제 별에 독립 적합한 "
            "해석 Moffat — 두 엔진 어느 쪽 모형도 아니다."
        ),
    )
    parser.add_argument("--moffat-stars", type=int, default=40)
    args = parser.parse_args()
    if args.trials <= 0 or args.injections <= 0:
        parser.error("--trials and --injections must be positive")
    if not args.source_fits.exists():
        raise FileNotFoundError(args.source_fits)
    args.baseline_result_dir = args.baseline_result_dir.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with fits.open(args.source_fits, memmap=True) as hdul:
        image = np.asarray(hdul[0].data)
        header = hdul[0].header.copy()
    detection_path = _detect_path(args.baseline_result_dir, args.source_fits.name)
    if detection_path is not None:
        detection = pd.read_csv(
            detection_path,
            usecols=lambda name: name in {"x", "y"},
        )
        real_positions = detection[["x", "y"]].apply(pd.to_numeric, errors="coerce").to_numpy(float) if {"x", "y"} <= set(detection.columns) else np.empty((0, 2))
    else:
        real_positions = np.empty((0, 2))
    fwhm_px = _read_fwhm(detection_path, args.fwhm_px)
    gain = float(args.gain_e_per_adu) if args.gain_e_per_adu is not None else _parameter_gain(args.parameter_file)
    background_rms = _read_background_rms(
        detection_path,
        args.background_rms_adu,
        image,
    )
    moffat_meta: dict = {}
    if args.inject_kernel == "moffat":
        step7_path = _baseline_step7(
            args.baseline_result_dir, args.source_fits.name, args.baseline_step7)
        if step7_path is None:
            raise FileNotFoundError("moffat 주입에는 기준선 step7 표가 필요하다")
        gamma, alpha, n_fit = fit_frame_moffat(
            image, step7_path, fwhm_px, n_stars=args.moffat_stars)
        implied_fwhm = 2.0 * gamma * float(np.sqrt(2.0 ** (1.0 / alpha) - 1.0))
        size = int(round(4.0 * fwhm_px)) | 1
        kernel_sampler = make_moffat_sampler(gamma, alpha, size)
        moffat_meta = {
            "inject_kernel": "moffat", "moffat_gamma": gamma,
            "moffat_alpha": alpha, "moffat_fit_stars": n_fit,
            "moffat_implied_fwhm_px": implied_fwhm,
        }
        print(
            f"주입 커널 = 독립 Moffat: gamma={gamma:.3f} alpha={alpha:.3f} "
            f"(별 {n_fit}개, 함의 FWHM {implied_fwhm:.2f}px vs 측정 {fwhm_px:.2f}px)")
    else:
        epsf_path = _baseline_epsf(args.baseline_result_dir, args.epsf_model)
        with fits.open(epsf_path, memmap=False) as hdul:
            epsf_data = np.asarray(hdul[0].data, dtype=float).copy()
            epsf_header = hdul[0].header.copy()
        def kernel_sampler(phase_x: float, phase_y: float) -> np.ndarray:
            return oversampled_epsf_to_native_kernel(
                epsf_data,
                header=epsf_header,
                phase_x=phase_x,
                phase_y=phase_y,
            )
        moffat_meta = {"inject_kernel": "epsf"}

    kernel = kernel_sampler(0.0, 0.0)
    nea = psf_noise_equivalent_area(kernel)
    residual_path = _baseline_residual(args.baseline_result_dir, args.source_fits.name)
    if residual_path is not None:
        with fits.open(residual_path, memmap=True) as hdul:
            baseline_residual = np.asarray(hdul[0].data)
    else:
        baseline_residual = None
    truths: list[pd.DataFrame] = []
    recoveries: list[pd.DataFrame] = []
    for trial in range(1, args.trials + 1):
        truth, recovery = _run_one_trial(args, trial, image, header, args.source_fits, kernel, kernel_sampler, baseline_residual, real_positions, fwhm_px, gain, background_rms, nea)
        truths.append(truth)
        recoveries.append(recovery)
        print(f"trial {trial}/{args.trials}: {len(recovery)} injections")
    truth_all = pd.concat(truths, ignore_index=True)
    recovery_all = pd.concat(recoveries, ignore_index=True)
    summary = aggregate_recovery_metrics(recovery_all)
    write_benchmark_outputs(truth_all, recovery_all, summary, args.output_dir, metadata={"source_fits": str(args.source_fits), "trials": args.trials, "injections_per_trial": args.injections, "seed": args.seed, "fwhm_px": fwhm_px, "gain_e_per_adu": gain, "background_rms_adu": background_rms, "psf_nea_px": nea, "flux_scale": args.flux_scale, "fit_shape_fwhm_mult": args.fit_shape_fwhm_mult, "fit_window_mode": args.fit_window_mode, "fit_encircled_energy": args.fit_encircled_energy, "postfit_qfit_noise_max": args.postfit_qfit_noise_max, "use_grouper": args.use_grouper, "grouper_max_size": args.grouper_max_size, "grouper_radius_fwhm": args.grouper_radius_fwhm, "forced_match_radius_fwhm": args.forced_match_radius_fwhm, "reused_baseline_detection": bool(args.reuse_baseline_detection), "baseline_residual": str(residual_path) if residual_path is not None else None, **moffat_meta})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
