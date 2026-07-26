from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("PyQt5")

from apex.analysis.cmd.isochrone_data import (
    extract_multicolor_arrays,
    infer_color_mag,
    resolve_wide_csv,
)
from apex.analysis.cmd.isochrone_fitter_v2 import (
    FitBounds,
    IsochroneFitterV2,
)
from apex.gui.workflow.cmd.step10_zeropoint_calibration import (
    ZeropointCalibrationWorker,
    resolve_cmd_photometry_input,
)
from apex.gui.workflow.cmd.step8_psf_photometry import (
    build_psf_output_signature,
    write_psf_output_signature,
)
from apex.gui.workflow.cmd.step12_isochrone_model import (
    _bands_from_df,
    _preferred_err_col,
    _preferred_mag_col,
)
from apex.utils.gaia_transforms import GAIA_TO_BAND
from apex.utils.photometry_provenance import summarize_photometry_table
from apex.utils.step_paths import step6_refbuild_dir, step7_forced_phot_dir
from apex.utils.step_paths_cmd import step8_psf_dir


def _poly_ascending(x: np.ndarray, coeffs: list[float]) -> np.ndarray:
    return sum(float(value) * x**power for power, value in enumerate(coeffs))


def _write_psf_project(root: Path, n_stars: int = 36) -> dict[str, np.ndarray]:
    ids = np.arange(1, n_stars + 1, dtype=int)
    source_ids = 10_000 + ids
    bp_rp = np.linspace(0.45, 1.45, n_stars)
    gaia_g = np.linspace(12.0, 17.0, n_stars)
    gaia_bp = gaia_g + 0.55 * bp_rp
    gaia_rp = gaia_bp - bp_rp

    reference: dict[str, np.ndarray] = {}
    psf_magnitudes: dict[str, np.ndarray] = {}
    frames = {"B": "synthetic_B.fits", "V": "synthetic_V.fits"}
    zeropoints = {"B": 4.25, "V": 3.80}
    for band in frames:
        coeffs = GAIA_TO_BAND[band][0]
        reference[band] = gaia_g - _poly_ascending(bp_rp, coeffs)
        psf_magnitudes[band] = reference[band] - zeropoints[band]

    step7_dir = step7_forced_phot_dir(root)
    step7_dir.mkdir(parents=True)
    step7_index = step7_dir / "photometry_index.csv"
    pd.DataFrame(
        {"file": list(frames.values()), "filter": list(frames.keys())}
    ).to_csv(step7_index, index=False)

    psf_dir = step8_psf_dir(root)
    psf_dir.mkdir(parents=True)
    pd.DataFrame(
        {"file": list(frames.values()), "filter": list(frames.keys())}
    ).to_csv(psf_dir / "photometry_index.csv", index=False)

    for band, frame in frames.items():
        aperture_table = step7_dir / f"photometry_{frame}.tsv"
        pd.DataFrame(
            {
                "ID": ids,
                "source_id": source_ids,
                "det_uid": ids,
                "x_fit": np.linspace(100.0, 900.0, n_stars),
                "y_fit": np.linspace(150.0, 950.0, n_stars),
                "FILTER": band,
                "mag": psf_magnitudes[band] + 0.15,
                "mag_err": np.full(n_stars, 0.03),
                "snr": np.full(n_stars, 50.0),
            }
        ).to_csv(aperture_table, sep="\t", index=False)
        pd.DataFrame(
            {
                "det_uid": ids,
                "seed_uid": ids,
                "x_fit": np.linspace(100.05, 900.05, n_stars),
                "y_fit": np.linspace(150.05, 950.05, n_stars),
                "FILTER": band,
                "mag_psf": psf_magnitudes[band],
                "mag_psf_err": np.full(n_stars, 0.01),
                "snr_psf": np.full(n_stars, 100.0),
                "flags_psf": np.zeros(n_stars, dtype=int),
            }
        ).to_csv(psf_dir / f"photometry_{frame}.tsv", sep="\t", index=False)

    master_dir = step6_refbuild_dir(root)
    master_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "ID": ids,
            "source_id": source_ids,
            "x_ref": np.linspace(100.0, 900.0, n_stars),
            "y_ref": np.linspace(150.0, 950.0, n_stars),
            "gaia_G": gaia_g,
            "gaia_BP": gaia_bp,
            "gaia_RP": gaia_rp,
        }
    ).to_csv(master_dir / "ref_catalog.tsv", sep="\t", index=False)

    (root / "project_state.json").write_text(
        json.dumps({"step_data": {"psf_photometry": {"skip_psf": False}}}),
        encoding="utf-8",
    )
    signature_params = SimpleNamespace(
        P=SimpleNamespace(result_dir=root, data_dir=root, cache_dir=root / "cache")
    )
    signature = build_psf_output_signature(
        signature_params,
        list(frames.values()),
    )
    write_psf_output_signature(root, signature)
    return {**psf_magnitudes, "ref_B": reference["B"], "ref_V": reference["V"]}


def _run_step10(root: Path) -> dict:
    p = SimpleNamespace(
        result_dir=root,
        data_dir=root,
        cache_dir=root / "cache",
        min_master_gaia_matches=10,
        min_snr_for_mag=3.0,
        gaia_snr_calib_min=5.0,
        cmd_snr_calib_min=5.0,
        frame_zp_min_n=5,
        phot_ref_require_apcorr_candidate=False,
        phot_ref_apcorr_min_keep=5,
        phot_use_qc_pass_only=False,
        cmd_apply_extinction=False,
        cmd_extinction_mode="absorb",
        zp_clip_sigma=3.0,
        zp_fit_iters=4,
        zp_slope_absmax=1.0,
        gaia_gi_min=-0.5,
        gaia_gi_max=3.5,
    )
    worker = ZeropointCalibrationWorker(
        SimpleNamespace(P=p), p.data_dir, p.result_dir, p.cache_dir
    )
    summaries: list[dict] = []
    errors: list[str] = []
    worker.finished.connect(summaries.append)
    worker.error.connect(errors.append)
    worker.run()
    assert not errors, errors[0] if errors else ""
    assert summaries and summaries[-1].get("ok") is True
    return summaries[-1]


def _isochrone_grid_from_cmd(df: pd.DataFrame) -> np.ndarray:
    b = pd.to_numeric(df["mag_std_B"], errors="coerce").to_numpy(float) - 10.0
    v = pd.to_numeric(df["mag_std_V"], errors="coerce").to_numpy(float) - 10.0
    mass = np.linspace(0.4, 1.4, len(df))
    rows = []
    for age in (8.9, 9.0, 9.1):
        for metallicity in (-0.1, 0.0, 0.1):
            rows.append(
                np.column_stack(
                    [
                        np.zeros(len(df)),
                        np.full(len(df), metallicity),
                        np.full(len(df), age),
                        b,
                        v,
                        mass,
                    ]
                )
            )
    return np.vstack(rows)


def test_complete_psf_reaches_step12_isochrone_grid_scan(tmp_path):
    truth = _write_psf_project(tmp_path)
    source = resolve_cmd_photometry_input(tmp_path)
    assert source["source"] == "psf"

    _run_step10(tmp_path)
    wide_path = resolve_wide_csv(tmp_path)
    cmd = pd.read_csv(wide_path)

    assert summarize_photometry_table(cmd) == {
        "source": "psf",
        "mag_column": "mag_psf",
        "mag_error_column": "mag_psf_err",
    }
    assert np.allclose(cmd["mag_inst_B"], truth["B"], atol=1e-10)
    assert np.allclose(cmd["mag_inst_V"], truth["V"], atol=1e-10)
    assert np.nanmedian(np.abs(cmd["mag_std_B"] - truth["ref_B"])) < 0.02
    assert np.nanmedian(np.abs(cmd["mag_std_V"] - truth["ref_V"])) < 0.02

    color_pairs, magnitude_bands = _bands_from_df(cmd)
    assert ("B", "V") in color_pairs
    assert {"B", "V"}.issubset(magnitude_bands)
    b_col, b_system = _preferred_mag_col(cmd.columns, "B")
    v_col, v_system = _preferred_mag_col(cmd.columns, "V")
    assert (b_col, b_system) == ("mag_cal_B", "Calibrated")
    assert (v_col, v_system) == ("mag_cal_V", "Calibrated")
    assert _preferred_err_col(cmd.columns, "B") == "mag_cal_err_B"
    assert _preferred_err_col(cmd.columns, "V") == "mag_cal_err_V"
    assert np.nanmedian(np.abs(cmd[b_col] - truth["ref_B"])) < 0.02
    assert np.nanmedian(np.abs(cmd[v_col] - truth["ref_V"])) < 0.02

    color, mag_band = infer_color_mag(cmd)
    assert color == ("B", "V")
    assert mag_band == "V"
    obs, err, labels, meta = extract_multicolor_arrays(
        cmd,
        colors=[color],
        mag_band=mag_band,
        data_snr_min=5.0,
        max_stars=0,
        seed=42,
    )
    assert labels == ["B-V", "V"]
    assert meta["columns"] == {"B": "mag_std_B", "V": "mag_std_V"}
    assert len(obs) == len(cmd) == 36

    fitter = IsochroneFitterV2(
        tmp_path / "synthetic_isochrone.dat",
        col_mh=1,
        col_age=2,
        col_g=3,
        col_r=4,
        col_mag=4,
        col_mass=5,
        iso_data=_isochrone_grid_from_cmd(cmd),
    )
    result = fitter.fit_grid_scan(
        obs[:, 0],
        obs[:, 1],
        err[:, 0],
        err[:, 1],
        bounds=FitBounds(
            log_age=(8.9, 9.1),
            metallicity=(-0.1, 0.1),
            distance_mod=(9.5, 10.5),
            extinction_gr=(0.0, 0.05),
        ),
        local_maxiter=10,
    )
    assert result.n_evaluated == 9
    assert result.best_fit.n_stars == 36
    assert result.best_fit.converged is True
    assert np.isfinite(result.best_fit.chi2)
