"""Run the GUI Step 4 and Step 8 workers headlessly on real FITS frames.

This is an execution harness, not a separate photometry implementation. It
uses the same DetectionWorker and Step6PSFWorker classes launched by the GUI,
while keeping validation outputs outside the user's project result directory.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from PyQt5.QtCore import QCoreApplication  # noqa: E402

from apex.config.parameters_cmd import Parameters  # noqa: E402
from apex.gui.workflow.cmd.step8_psf_photometry import Step6PSFWorker  # noqa: E402
from apex.gui.workflow.step4_source_detection import DetectionWorker  # noqa: E402


def _print_signal(prefix: str):
    def emit(*values):
        rendered = " | ".join(str(value) for value in values)
        message = f"[{prefix}] {rendered}"
        encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
        safe_message = message.encode(encoding, errors="backslashreplace").decode(encoding)
        print(safe_message, flush=True)

    return emit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path)
    parser.add_argument("--fit-engine", choices=("apex_iterative", "allstar"), default="apex_iterative")
    parser.add_argument("--fitter-max-iter", type=int)
    parser.add_argument("--fit-shape-fwhm-mult", type=float)
    parser.add_argument("--fit-window-mode", choices=("auto", "manual"))
    parser.add_argument("--fit-encircled-energy", type=float)
    parser.add_argument("--postfit-qfit-noise-max", type=float)
    parser.add_argument("--residual-passes", type=int)
    parser.add_argument("--epsf-max-stars", type=int)
    parser.add_argument("--epsf-contamination", choices=("on", "off"))
    parser.add_argument("--flux-scale", choices=("on", "off"))
    parser.add_argument("--use-grouper", choices=("on", "off"))
    parser.add_argument("--grouper-max-size", type=int)
    parser.add_argument("--grouper-radius-fwhm", type=float)
    parser.add_argument("--forced-match-radius-fwhm", type=float)
    parser.add_argument("--core-cut", choices=("on", "off"))
    parser.add_argument("--core-center-mode", choices=("auto", "image", "manual"))
    parser.add_argument("--core-x", type=float)
    parser.add_argument("--core-y", type=float)
    parser.add_argument("--core-radius", type=float)
    parser.add_argument(
        "--compact-output",
        action="store_true",
        help="Keep GUI-default residual products instead of every iteration image.",
    )
    parser.add_argument(
        "--skip-step4",
        action="store_true",
        help="Reuse an existing Step4 table in result-dir and run Step8 only.",
    )
    parser.add_argument("files", nargs="+")
    args = parser.parse_args()

    missing = [name for name in args.files if not (args.data_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing FITS files: {missing}")

    parameter_file = args.parameter_file or (REPO / "parameters.toml")
    if not parameter_file.exists():
        raise FileNotFoundError(f"Missing parameter file: {parameter_file}")
    if args.core_center_mode == "manual" and (args.core_x is None or args.core_y is None):
        raise ValueError("--core-center-mode manual requires --core-x and --core-y")

    app = QCoreApplication.instance() or QCoreApplication([])
    params = Parameters(parameter_file)
    params.P.data_dir = str(args.data_dir)
    params.P.result_dir = str(args.result_dir)
    params.P.cache_dir = str(args.result_dir / "cache")
    params.P.max_workers = 1
    params.P.psf_parallel_workers = 1
    params.P.force_redetect = not args.skip_step4
    params.P.force_rephot = True
    params.P.psf_save_residuals = True
    params.P.psf_save_model_image = True
    params.P.psf_save_all_iter_residuals = not args.compact_output
    params.P.psf_fit_engine = args.fit_engine
    if args.fitter_max_iter is not None:
        params.P.psf_fitter_max_iter = max(1, int(args.fitter_max_iter))
    if args.fit_shape_fwhm_mult is not None:
        params.P.psf_fit_shape_fwhm_mult = min(
            4.0, max(1.0, float(args.fit_shape_fwhm_mult))
        )
        if args.fit_window_mode is None:
            params.P.psf_fit_window_mode = "manual"
    if args.fit_window_mode is not None:
        params.P.psf_fit_window_mode = args.fit_window_mode
    if args.fit_encircled_energy is not None:
        params.P.psf_fit_encircled_energy = min(
            0.995, max(0.50, float(args.fit_encircled_energy))
        )
    if args.postfit_qfit_noise_max is not None:
        params.P.psf_postfit_qfit_max = max(
            0.0, float(args.postfit_qfit_noise_max)
        )
    if args.residual_passes is not None:
        params.P.psf_max_iter = max(1, int(args.residual_passes))
    if args.epsf_max_stars is not None:
        params.P.psf_n_stars_max = max(0, int(args.epsf_max_stars))
    if args.epsf_contamination is not None:
        params.P.psf_epsf_contamination_filter = args.epsf_contamination == "on"
    if args.flux_scale is not None:
        params.P.psf_flux_scale_correction = args.flux_scale == "on"
    if args.use_grouper is not None:
        params.P.psf_use_grouper = args.use_grouper == "on"
    if args.grouper_max_size is not None:
        params.P.psf_grouper_max_size = min(3, max(1, int(args.grouper_max_size)))
    if args.grouper_radius_fwhm is not None:
        params.P.psf_grouper_radius_fwhm = min(
            5.0, max(0.5, float(args.grouper_radius_fwhm))
        )
    if args.forced_match_radius_fwhm is not None:
        params.P.psf_forced_match_radius_fwhm = min(
            3.0, max(0.1, float(args.forced_match_radius_fwhm))
        )
    if args.core_cut is not None:
        params.P.psf_core_cut_enable = args.core_cut == "on"
    if args.core_center_mode is not None:
        params.P.psf_core_cut_center_mode = args.core_center_mode
    if args.core_x is not None:
        params.P.psf_core_cut_x_px = float(args.core_x)
    if args.core_y is not None:
        params.P.psf_core_cut_y_px = float(args.core_y)
    if args.core_radius is not None:
        params.P.psf_core_cut_radius_px = max(0.0, float(args.core_radius))

    args.result_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = args.result_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    run_metadata = {
        "data_dir": str(args.data_dir),
        "result_dir": str(args.result_dir),
        "files": args.files,
        "one_cpu_worker": True,
        "compact_output": bool(args.compact_output),
        "step4_skipped": bool(args.skip_step4),
        "parameter_file": str(parameter_file),
        "psf_fit_engine": getattr(params.P, "psf_fit_engine", None),
        "psf_model_mode": getattr(params.P, "psf_model_mode", None),
        "psf_max_iter": getattr(params.P, "psf_max_iter", None),
        "psf_fitter_max_iter": getattr(params.P, "psf_fitter_max_iter", None),
        "psf_fit_shape_fwhm_mult": getattr(
            params.P, "psf_fit_shape_fwhm_mult", None
        ),
        "psf_fit_window_mode": getattr(params.P, "psf_fit_window_mode", None),
        "psf_fit_encircled_energy": getattr(
            params.P, "psf_fit_encircled_energy", None
        ),
        "psf_postfit_qfit_noise_max": getattr(
            params.P, "psf_postfit_qfit_max", None
        ),
        "psf_n_stars_max": getattr(params.P, "psf_n_stars_max", None),
        "psf_epsf_contamination_filter": getattr(
            params.P,
            "psf_epsf_contamination_filter",
            None,
        ),
        "psf_flux_scale_correction": getattr(
            params.P,
            "psf_flux_scale_correction",
            None,
        ),
        "psf_use_grouper": getattr(params.P, "psf_use_grouper", None),
        "psf_grouper_max_size": getattr(params.P, "psf_grouper_max_size", None),
        "psf_grouper_radius_fwhm": getattr(
            params.P, "psf_grouper_radius_fwhm", None
        ),
        "psf_forced_match_radius_fwhm": getattr(
            params.P, "psf_forced_match_radius_fwhm", None
        ),
        "psf_field_mode": getattr(params.P, "psf_field_mode", None),
        "psf_core_cut_enable": getattr(params.P, "psf_core_cut_enable", None),
        "psf_core_cut_center_mode": getattr(params.P, "psf_core_cut_center_mode", None),
        "psf_core_cut_x_px": getattr(params.P, "psf_core_cut_x_px", None),
        "psf_core_cut_y_px": getattr(params.P, "psf_core_cut_y_px", None),
        "psf_core_cut_radius_px": getattr(params.P, "psf_core_cut_radius_px", None),
    }

    if args.skip_step4:
        run_metadata["step4_elapsed_s"] = 0.0
        run_metadata["step4_summary"] = {"status": "reused"}
    else:
        detection_summary: list[dict] = []
        detector = DetectionWorker(
            args.files, params, args.data_dir, cache_dir, use_cropped=False
        )
        detector.progress.connect(_print_signal("STEP4"))
        detector.worker_status.connect(_print_signal("STEP4-WORKER"))
        detector.file_done.connect(_print_signal("STEP4-DONE"))
        detector.error.connect(_print_signal("STEP4-ERROR"))

        def detection_finished(value):
            detection_summary.append(value)
            app.quit()

        detector.finished.connect(detection_finished)
        started = time.perf_counter()
        detector.start()
        app.exec_()
        detector.wait()
        run_metadata["step4_elapsed_s"] = time.perf_counter() - started
        run_metadata["step4_summary"] = (
            detection_summary[-1] if detection_summary else {}
        )

    psf_summary: list[dict] = []
    worker = Step6PSFWorker(
        args.files,
        params,
        args.data_dir,
        args.result_dir,
        cache_dir,
        use_cropped=False,
    )
    worker.progress.connect(_print_signal("STEP8"))
    worker.worker_status.connect(_print_signal("STEP8-WORKER"))
    worker.frame_done.connect(_print_signal("STEP8-DONE"))
    worker.error.connect(_print_signal("STEP8-ERROR"))
    worker.log.connect(_print_signal("STEP8-LOG"))
    def psf_finished(value):
        psf_summary.append(value)
        app.quit()

    worker.finished.connect(psf_finished)
    started = time.perf_counter()
    worker.start()
    app.exec_()
    worker.wait()
    run_metadata["step8_elapsed_s"] = time.perf_counter() - started
    run_metadata["step8_summary"] = psf_summary[-1] if psf_summary else {}

    (args.result_dir / "real_gui_run.json").write_text(
        json.dumps(run_metadata, indent=2, default=str), encoding="utf-8"
    )
    print(json.dumps(run_metadata, indent=2, default=str), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
