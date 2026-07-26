from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from apex.analysis.light_curve.photometry_source_service import (
    load_filter_photometry_timeseries,
    load_lightcurve_frame_photometry,
    resolve_lightcurve_photometry_source,
)
import apex.analysis.light_curve.photometry_source_service as source_service
from apex.core.project_state import ProjectState
from apex.gui.main_window import _migrate_lc_optional_psf_state
from apex.gui.workflow.lc.step8_target_selection import (
    TargetComparisonSelectionWindow,
    _simbad_type_is_variable,
)
from apex.gui.workflow.lc.step9_lightcurve_builder import LightCurveBuilderWindow
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_cmd import step8_psf_dir
from apex.utils.step_paths_lc import step8_selection_dir


def _write_complete_sources(root, frame: str = "frame001.fits"):
    aperture_dir = step7_forced_phot_dir(root)
    aperture_dir.mkdir(parents=True)
    aperture_index = aperture_dir / "photometry_index.csv"
    pd.DataFrame({"file": [frame], "filter": ["g"]}).to_csv(
        aperture_index, index=False
    )
    pd.DataFrame(
        {
            "ID": [101, 102, 103],
            "source_id": [1001, 1002, 1003],
            "det_uid": [1, 2, -1],
            "x": [9.5, 19.5, 30.0],
            "y": [10.5, 20.5, 31.0],
            "mag": [12.3, 13.4, 14.5],
            "mag_err": [0.02, 0.04, 0.05],
        }
    ).to_csv(aperture_dir / f"photometry_{frame}.tsv", sep="\t", index=False)

    psf_dir = step8_psf_dir(root)
    psf_dir.mkdir(parents=True)
    pd.DataFrame({"file": [frame], "filter": ["g"]}).to_csv(
        psf_dir / "photometry_index.csv", index=False
    )
    pd.DataFrame(
        {
            "det_uid": [1, 2, -100],
            "x_fit": [10.0, 20.0, 40.0],
            "y_fit": [11.0, 21.0, 41.0],
            "mag_psf": [12.1, 13.2, 15.0],
            "mag_psf_err": [0.01, 0.03, 0.06],
            "flags_psf": [0, 2, 0],
            "qfit": [0.05, 0.8, 0.1],
            "cfit": [0.01, -0.2, 0.02],
            "reduced_chi2": [1.1, 5.0, 1.2],
            "snr_psf": [100.0, 20.0, 10.0],
            "iter_found": [1, 1, 2],
        }
    ).to_csv(psf_dir / f"photometry_{frame}.tsv", sep="\t", index=False)
    stat = aperture_index.stat()
    aperture_table = aperture_dir / f"photometry_{frame}.tsv"
    table_stat = aperture_table.stat()
    signature = {
        "frames": [frame],
        "inputs": {
            "step7_index": {
                "path": str(aperture_index.resolve()),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            },
            "frames": [
                {
                    "file": frame,
                    "step7_tsv": {
                        "path": str(aperture_table.resolve()),
                        "size": table_stat.st_size,
                        "mtime_ns": table_stat.st_mtime_ns,
                    },
                }
            ],
        },
    }
    (psf_dir / "psf_output_signature.json").write_text(
        json.dumps(signature), encoding="utf-8"
    )
    return frame


def test_complete_psf_output_is_selected_and_normalized(tmp_path):
    frame = _write_complete_sources(tmp_path)
    source = resolve_lightcurve_photometry_source(tmp_path)
    assert source["source"] == "psf"
    loaded = load_lightcurve_frame_photometry(tmp_path, frame, source)
    assert loaded["ID"].tolist() == [101, 102, 103]
    assert loaded["source_id"].tolist() == [1001, 1002, 1003]
    assert loaded.loc[0, "mag"] == 12.1
    assert np.isnan(loaded.loc[1, "mag"])
    assert np.isnan(loaded.loc[2, "mag"])
    assert loaded.loc[0, "x"] == 10.0
    assert loaded.loc[0, "qfit"] == 0.05
    assert loaded.loc[0, "reduced_chi2"] == 1.1
    assert bool(loaded.loc[0, "psf_qc_clean"])
    assert not bool(loaded.loc[1, "psf_qc_clean"])
    assert set(loaded["photometry_source"]) == {"psf"}
    assert set(loaded["mag_input_column"]) == {"mag_psf"}
    assert set(loaded["mag_error_input_column"]) == {"mag_psf_err"}


def test_psf_resolution_trusts_tables_older_than_completion_signature(
    tmp_path, monkeypatch
):
    frame = _write_complete_sources(tmp_path)
    output_table = step8_psf_dir(tmp_path) / f"photometry_{frame}.tsv"
    real_read_csv = source_service.pd.read_csv
    output_reads = []

    def tracked_read_csv(path, *args, **kwargs):
        if Path(path) == output_table:
            output_reads.append(path)
        return real_read_csv(path, *args, **kwargs)

    monkeypatch.setattr(source_service.pd, "read_csv", tracked_read_csv)
    source = resolve_lightcurve_photometry_source(tmp_path)

    assert source["source"] == "psf"
    assert output_reads == []


def test_psf_resolution_revalidates_table_changed_after_signature(tmp_path):
    frame = _write_complete_sources(tmp_path)
    psf_dir = step8_psf_dir(tmp_path)
    output_table = psf_dir / f"photometry_{frame}.tsv"
    signature_path = psf_dir / "psf_output_signature.json"
    broken = pd.read_csv(output_table, sep="\t").drop(columns=["mag_psf"])
    broken.to_csv(output_table, sep="\t", index=False)
    newer = signature_path.stat().st_mtime_ns + 1_000_000_000
    os.utime(output_table, ns=(newer, newer))

    source = resolve_lightcurve_photometry_source(tmp_path)

    assert source["source"] == "aperture"
    assert "Incomplete PSF table" in source["reason"]


def test_skip_state_forces_aperture_even_when_psf_exists(tmp_path):
    frame = _write_complete_sources(tmp_path)
    state = ProjectState(tmp_path / ".state")
    state.store_step_data("psf_photometry", {"skip_psf": True})
    source = resolve_lightcurve_photometry_source(tmp_path, state)
    assert source["source"] == "aperture"
    loaded = load_lightcurve_frame_photometry(tmp_path, frame, source)
    assert loaded.loc[0, "mag"] == 12.3
    assert set(loaded["photometry_source"]) == {"aperture"}
    assert set(loaded["mag_input_column"]) == {"mag"}
    assert set(loaded["mag_error_input_column"]) == {"mag_err"}


def test_mirrored_skip_state_is_used_for_an_external_dataset(tmp_path):
    _write_complete_sources(tmp_path)
    (tmp_path / "project_state.json").write_text(
        json.dumps({"step_data": {"psf_photometry": {"skip_psf": True}}}),
        encoding="utf-8",
    )

    source = resolve_lightcurve_photometry_source(tmp_path)

    assert source["source"] == "aperture"
    assert source["psf_skipped"] is True


def test_changed_step7_signature_falls_back_to_aperture(tmp_path):
    _write_complete_sources(tmp_path)
    index = step7_forced_phot_dir(tmp_path) / "photometry_index.csv"
    index.write_text(index.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    source = resolve_lightcurve_photometry_source(tmp_path)
    assert source["source"] == "aperture"
    assert "Step 7 changed" in source["reason"]


def test_legacy_lc_progress_is_shifted_and_psf_marked_skipped(tmp_path):
    state = ProjectState(tmp_path)
    state.state["completed_steps"] = list(range(9))
    state.state["current_step"] = 8
    state.save()
    assert _migrate_lc_optional_psf_state(state) is True
    assert state.state["completed_steps"] == list(range(10))
    assert state.state["current_step"] == 9
    assert state.get_step_data("psf_photometry")["skip_psf"] is True
    assert _migrate_lc_optional_psf_state(state) is False


def test_instrumental_magnitude_map_includes_non_gaia_sources(tmp_path):
    forced_dir = step7_forced_phot_dir(tmp_path)
    forced_dir.mkdir(parents=True)
    for frame, mags in (("a.fits", [15.0, 16.0]), ("b.fits", [15.2, 16.4])):
        pd.DataFrame(
            {
                "source_id": [-1, 2001],
                "ID": [1, 2],
                "det_uid": [10, 11],
                "mag_inst": mags,
                "bad_phot_flag": [False, False],
            }
        ).to_csv(forced_dir / f"photometry_{frame}.tsv", sep="\t", index=False)

    window = TargetComparisonSelectionWindow.__new__(TargetComparisonSelectionWindow)
    window.params = SimpleNamespace(P=SimpleNamespace(result_dir=tmp_path))
    window.filter_frames = {"g": ["a.fits", "b.fits"]}
    window._instrumental_mag_cache = {}

    medians = window._get_instrumental_mag_map("g")

    assert medians[-1] == 15.1
    assert medians[2001] == 16.2


def test_filter_timeseries_loader_filters_frames_and_bad_measurements(tmp_path):
    forced_dir = step7_forced_phot_dir(tmp_path)
    forced_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "file": ["g1.fits", "r1.fits"],
            "filter": ["g", "r"],
            "JD": [2460000.1, 2460000.2],
            "airmass": [1.2, 1.4],
        }
    ).to_csv(forced_dir / "photometry_index.csv", index=False)
    for frame, offset in (("g1.fits", 0.0), ("r1.fits", 1.0)):
        pd.DataFrame(
            {
                "source_id": [1001, 1002],
                "ID": [1, 2],
                "mag": [13.0 + offset, 14.0 + offset],
                "mag_err": [0.01, 0.02],
                "bad_phot_flag": [False, True],
                "x": [10.0, 20.0],
                "y": [11.0, 21.0],
            }
        ).to_csv(forced_dir / f"photometry_{frame}.tsv", sep="\t", index=False)

    output, source = load_filter_photometry_timeseries(tmp_path, "g")

    assert source["source"] == "aperture"
    assert output["frame"].tolist() == ["g1.fits"]
    assert output["star_id"].tolist() == [1001]
    assert output.iloc[0]["time"] == 2460000.1
    assert output.iloc[0]["airmass"] == 1.2
    assert output.iloc[0]["photometry_source"] == "aperture"
    assert output.iloc[0]["mag_input_column"] == "mag"
    assert output.iloc[0]["mag_error_input_column"] == "mag_err"


def test_filter_timeseries_loader_stops_between_frame_reads(tmp_path, monkeypatch):
    forced_dir = step7_forced_phot_dir(tmp_path)
    forced_dir.mkdir(parents=True)
    pd.DataFrame(
        {"file": ["g1.fits", "g2.fits"], "filter": ["g", "g"]}
    ).to_csv(forced_dir / "photometry_index.csv", index=False)
    reads = []

    def fake_frame_loader(root, filename, source, filter_name=None):
        reads.append(filename)
        return pd.DataFrame(
            {
                "source_id": [1001],
                "mag": [13.0],
                "mag_err": [0.01],
            }
        )

    monkeypatch.setattr(
        source_service, "load_lightcurve_frame_photometry", fake_frame_loader
    )
    output, _ = load_filter_photometry_timeseries(
        tmp_path,
        "g",
        should_stop=lambda: bool(reads),
    )

    assert reads == ["g1.fits"]
    assert output.empty


def test_legacy_forced_photometry_directory_is_loaded(tmp_path):
    forced_dir = tmp_path / "step5_photometry"
    forced_dir.mkdir()
    frame = "legacy001.fits"
    pd.DataFrame(
        {"file": [frame], "filter": ["g"], "JD": [2460001.25], "night_id": [2]}
    ).to_csv(forced_dir / "photometry_index.csv", index=False)
    pd.DataFrame(
        {
            "det_uid": [7],
            "mag": [13.2],
            "mag_err": [0.01],
            "xcenter": [10.0],
            "ycenter": [11.0],
        }
    ).to_csv(forced_dir / f"{frame}_photometry.tsv", sep="\t", index=False)
    identity_dir = tmp_path / "step8_idmatch" / "20260328"
    identity_dir.mkdir(parents=True)
    gaia_source_id = 1_387_320_379_874_985_728
    pd.DataFrame(
        {"det_idx": [7], "source_id": [gaia_source_id], "x": [10.0], "y": [11.0]}
    ).to_csv(identity_dir / f"idmatch_{frame}.csv", index=False)

    source = resolve_lightcurve_photometry_source(tmp_path)
    output, _ = load_filter_photometry_timeseries(tmp_path, "g")

    assert source["source"] == "aperture"
    assert source["directory"] == forced_dir
    assert output["frame"].tolist() == [frame]
    assert output["star_id"].dtype == "int64"
    assert output["star_id"].tolist() == [gaia_source_id]
    assert output.iloc[0]["time"] == 2460001.25
    assert output.iloc[0]["night_id"] == 2


def test_candidate_pool_rescreens_check_and_respects_catalog_exclusions():
    rows = []
    for frame_index in range(10):
        for star_id in range(1, 8):
            rows.append(
                {
                    "frame": f"f{frame_index}",
                    "star_id": star_id,
                    "mag": 12.0 + 0.2 * star_id,
                    "mag_err": 0.01,
                }
            )
    window = TargetComparisonSelectionWindow.__new__(TargetComparisonSelectionWindow)
    window.params = SimpleNamespace(P=SimpleNamespace(comparison_auto_pool_max=30))
    window.filter_check_stars = {"g": 3}
    window.filter_rejected_sources = {"g": {2}}
    window.filter_comparisons = {"g": set()}
    window._gaia_variable_flag_map = {4: "VARIABLE"}

    pool, report, _ = window._comparison_candidate_pool(
        "g",
        pd.DataFrame(rows),
        target_id=1,
        target_count=3,
        catalog_reject_ids={5},
    )

    assert pool == [3, 6, 7]
    reasons = report.set_index("star_id")["basic_reason"].to_dict()
    assert reasons[1] == "target"
    assert reasons[2] == "manual_reject"
    assert reasons[3] == ""
    assert reasons[4] == "gaia_variable"
    assert reasons[5] == "simbad_variable"


def test_target_cannot_be_assigned_as_check_star(monkeypatch):
    window = TargetComparisonSelectionWindow.__new__(TargetComparisonSelectionWindow)
    window.current_filter = "g"
    window.selected_source_id = 101
    window.target_source_id = 101
    window.filter_check_stars = {"g": None}
    messages = []
    monkeypatch.setattr(
        "apex.gui.workflow.lc.step8_target_selection.QMessageBox.information",
        lambda *args: messages.append(args[-1]),
    )

    window.set_check_selected()

    assert window.filter_check_stars["g"] is None
    assert messages == ["Target cannot be the independent check star."]


def test_assigning_target_removes_existing_check_role():
    window = TargetComparisonSelectionWindow.__new__(TargetComparisonSelectionWindow)
    window.current_filter = "g"
    window.selected_source_id = 101
    window.target_source_id = None
    window.master_ids = {101}
    window.filter_targets = {"g": None}
    window.filter_check_stars = {"g": 101}
    window.filter_rejected_sources = {"g": set()}
    window.comparison_ids = set()
    window.filter_comparisons = {"g": set()}
    window._filter_stability_cache = {}
    window.save_selection = lambda: None
    window._queue_selection_save = lambda: None
    window._refresh_role_ui = lambda **_kwargs: None
    window.log = lambda _message: None

    window.set_target_selected()

    assert window.target_source_id == 101
    assert window.filter_check_stars["g"] is None


def test_selection_validation_requires_three_comparisons_and_check(tmp_path):
    output_dir = step8_selection_dir(tmp_path)
    output_dir.mkdir(parents=True)
    window = TargetComparisonSelectionWindow.__new__(TargetComparisonSelectionWindow)
    window.params = SimpleNamespace(P=SimpleNamespace(result_dir=tmp_path))
    window.filter_frames = {"g": ["g1.fits"]}
    selection_path = output_dir / "selection_g.json"
    selection_path.write_text(
        json.dumps(
            {
                "target_source_id": 1,
                "comparison_source_ids": [2, 3, 4],
                "check_source_id": None,
            }
        ),
        encoding="utf-8",
    )

    assert window.validate_step() is False

    selection_path.write_text(
        json.dumps(
            {
                "target_source_id": 1,
                "comparison_source_ids": [2, 3, 4],
                "check_source_id": 5,
            }
        ),
        encoding="utf-8",
    )
    assert window.validate_step() is True


def test_simbad_variable_types_are_screened_conservatively():
    assert _simbad_type_is_variable("V*") is True
    assert _simbad_type_is_variable("RR*") is True
    assert _simbad_type_is_variable("PulsV*") is True
    assert _simbad_type_is_variable("Star") is False


def test_multi_dataset_builder_forces_one_aperture_source_when_mixed(tmp_path):
    psf_root = tmp_path / "psf"
    aperture_root = tmp_path / "aperture"
    _write_complete_sources(psf_root)
    forced_dir = step7_forced_phot_dir(aperture_root)
    forced_dir.mkdir(parents=True)
    pd.DataFrame({"file": ["frame001.fits"], "filter": ["g"]}).to_csv(
        forced_dir / "photometry_index.csv", index=False
    )

    window = LightCurveBuilderWindow.__new__(LightCurveBuilderWindow)
    window.datasets = [("psf", psf_root), ("aperture", aperture_root)]
    window.params = SimpleNamespace(P=SimpleNamespace(result_dir=psf_root))
    window.project_state = None
    window._photometry_source_cache = {}
    window._force_aperture_for_datasets = False
    window._photometry_cache = {}
    window._photometry_cache_dir = None
    window._diff_series_cache = {}
    window._check_series_cache = {}

    window._refresh_photometry_source_policy()

    assert window._force_aperture_for_datasets is True
    source = window._photometry_source_for_dir(psf_root)
    assert source["source"] == "aperture"
    assert source["mag_column"] == "mag"
    assert source["mag_error_column"] == "mag_err"
