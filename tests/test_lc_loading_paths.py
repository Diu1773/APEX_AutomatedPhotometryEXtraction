from __future__ import annotations

from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import EarthLocation
from astropy.time import Time

import apex.gui.workflow.lc.step10_detrend_merge as detrend_module
from apex.gui.workflow.lc.step8_target_selection import TargetComparisonSelectionWindow
from apex.gui.workflow.lc.step10_detrend_merge import DetrendNightMergeWindow
from apex.utils.step_paths_lc import step1_dir


def test_step9_preassign_uses_loaded_aggregate_ids_without_frame_scan():
    window = TargetComparisonSelectionWindow.__new__(TargetComparisonSelectionWindow)
    window.filter_catalogs = {
        "V": pd.DataFrame({"source_id": pd.Series([11, 12], dtype="int64")})
    }
    window.filter_master_ids = {"V": {12, 13}}
    window._id_registry = {"V": {10: 110}}
    window._load_id_registry = lambda _flt: None
    window._load_global_id_map = lambda _flt: None
    window._get_or_assign_stable_id = lambda _flt, sid: sid + 100
    window._save_id_registry = lambda _flt: None
    window.log = lambda _message: None

    def fail_frame_scan(_flt):
        raise AssertionError("per-frame idmatch tables should not be scanned")

    window._collect_filter_source_ids = fail_frame_scan

    window._preassign_all_ids("V")

    assert window.sid_to_id == {10: 110, 11: 111, 12: 112, 13: 113}
    assert window.id_to_sid == {110: 10, 111: 11, 112: 12, 113: 13}


def test_step11_fills_airmass_from_header_times_without_opening_fits(
    tmp_path, monkeypatch
):
    location = EarthLocation(lat=37.0 * u.deg, lon=127.0 * u.deg, height=100.0 * u.m)
    jd_values = np.array([2460000.0, 2460000.001, 2460000.002])
    target_ra = Time(jd_values[1], format="jd", scale="utc").sidereal_time(
        "apparent", longitude=location.lon
    ).deg

    headers_dir = step1_dir(tmp_path)
    headers_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "Filename": ["a.fits", "b.fits", "c.fits"],
            "JD": jd_values,
            "AIRMASS": [np.nan, np.nan, np.nan],
        }
    ).to_csv(headers_dir / "headers.csv", index=False)

    window = DetrendNightMergeWindow.__new__(DetrendNightMergeWindow)
    window.datasets = [("night", tmp_path)]
    window.raw_df = pd.DataFrame(
        {
            "dataset": ["night", "night", "night"],
            "file": ["a.fits", "b.fits", "c.fits"],
            "airmass": [np.nan, np.nan, np.nan],
        }
    )
    window.params = SimpleNamespace(
        P=SimpleNamespace(
            result_dir=tmp_path,
            site_lat_deg=37.0,
            site_lon_deg=127.0,
            site_alt_m=100.0,
            site_tz_offset_hours=9.0,
            target_ra_deg=target_ra,
            target_dec_deg=37.0,
            airmass_formula="Kasten & Young (1989)",
        ),
        get_file_path=lambda _filename: tmp_path / "missing.fits",
    )
    messages = []
    window.log = messages.append

    def fail_fits_open(*_args, **_kwargs):
        raise AssertionError("FITS fallback should not run")

    monkeypatch.setattr(detrend_module.fits, "open", fail_fits_open)

    window._fill_airmass_from_headers()

    values = pd.to_numeric(window.raw_df["airmass"], errors="coerce").to_numpy()
    assert np.all(np.isfinite(values))
    assert np.all((values >= 0.95) & (values < 1.1))
    assert any("header_times=3" in message for message in messages)
