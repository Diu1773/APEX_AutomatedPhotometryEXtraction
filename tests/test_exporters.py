"""Tests for community-standard photometry exporters (apex.io.exporters).

Qt-free by construction: these import only apex.io.exporters + pandas/numpy, so
they run in the headless CI job. Do NOT import apex.gui here.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from apex.io.exporters import (
    AAVSO_COLUMNS,
    EXOCLOCK_COLUMNS,
    EXOFOP_COLUMNS,
    export_aavso_extended,
    export_exoclock,
    export_exofop,
)


def _synthetic_lc(n: int = 6) -> pd.DataFrame:
    """A small APEX-shaped differential light curve (raw + detrended columns)."""
    rng = np.random.default_rng(42)
    jd = 2459000.5 + np.arange(n) * 0.01
    diff_mag_raw = 12.5 + rng.normal(0, 0.01, n)
    return pd.DataFrame(
        {
            "file": [f"frame_{i:03d}.fits" for i in range(n)],
            "filter": ["V"] * n,
            "date": ["2024-06-20"] * n,
            "night_id": [1] * n,
            "JD": jd,
            "BJD_TDB": jd + 0.0003,
            "rel_time_hr": (jd - jd[0]) * 24.0,
            "mag": diff_mag_raw + 3.0,
            "mag_err": np.full(n, 0.012),
            "comp_avg": np.full(n, 3.0),
            "comp_err": np.full(n, 0.008),
            "diff_mag_raw": diff_mag_raw,
            "diff_mag_corr": diff_mag_raw - 0.002,
            "diff_err": np.full(n, 0.014),
            "airmass": 1.2 + np.arange(n) * 0.01,
        }
    )


def _meta() -> dict:
    return {
        "obscode": "ABC",
        "target": "TEST-VAR-01",
        "filter": "V",
        "observer": "Tester",
        "telescope": "0.4m SCT",
    }


# ── AAVSO ─────────────────────────────────────────────────────────────────────

def test_aavso_header_and_columns(tmp_path):
    out = export_aavso_extended(_synthetic_lc(), _meta(), tmp_path / "aavso.txt")
    assert out.exists()
    lines = out.read_text(encoding="utf-8").splitlines()

    # Required AAVSO Extended File Format header keys.
    assert lines[0] == "#TYPE=EXTENDED"
    header_block = lines[:6]
    assert "#OBSCODE=ABC" in header_block
    assert any(l.startswith("#SOFTWARE=APEX") for l in header_block)
    assert "#DELIM=," in header_block
    assert any(l.startswith("#DATE=") for l in header_block)
    assert any(l.startswith("#OBSTYPE=") for l in header_block)

    # Column header line is exactly the 15 AAVSO columns, comma-delimited.
    col_line = lines[6]
    assert col_line.split(",") == AAVSO_COLUMNS
    assert len(AAVSO_COLUMNS) == 15


def test_aavso_data_row_parseable(tmp_path):
    df = _synthetic_lc()
    out = export_aavso_extended(df, _meta(), tmp_path / "aavso.txt")

    parsed = pd.read_csv(out, comment="#")
    assert list(parsed.columns) == AAVSO_COLUMNS
    assert len(parsed) == len(df)

    first = parsed.iloc[0]
    assert first["NAME"] == "TEST-VAR-01"
    assert first["FILT"] == "V"
    # MAG should match the corrected differential mag (preferred over raw).
    assert math.isclose(float(first["MAG"]), float(df["diff_mag_corr"].iloc[0]), abs_tol=1e-4)
    assert math.isclose(float(first["MERR"]), 0.014, abs_tol=1e-4)
    assert math.isclose(float(first["AMASS"]), float(df["airmass"].iloc[0]), abs_tol=1e-4)


def test_aavso_unknown_fields_use_na(tmp_path):
    df = _synthetic_lc()
    out = export_aavso_extended(df, {"obscode": "ABC", "target": "X"}, tmp_path / "aavso.txt")
    parsed = pd.read_csv(out, comment="#", keep_default_na=False)
    # No check star supplied -> KNAME/KMAG are the AAVSO 'na' sentinel.
    assert (parsed["KNAME"] == "na").all()
    assert (parsed["KMAG"] == "na").all()
    assert (parsed["CNAME"] == "na").all()


def test_aavso_requires_target(tmp_path):
    with pytest.raises(ValueError):
        export_aavso_extended(_synthetic_lc(), {"obscode": "ABC"}, tmp_path / "x.txt")


# ── ExoClock ──────────────────────────────────────────────────────────────────

def test_exoclock_header_columns_and_flux_normalization(tmp_path):
    df = _synthetic_lc()
    out = export_exoclock(df, _meta(), tmp_path / "exoclock.csv")
    assert out.exists()
    text = out.read_text(encoding="utf-8")
    assert "# target=TEST-VAR-01" in text
    assert any(l.startswith("# software=APEX") for l in text.splitlines())

    parsed = pd.read_csv(out, comment="#")
    assert list(parsed.columns) == EXOCLOCK_COLUMNS == ["BJD_TDB", "flux", "flux_err"]
    assert len(parsed) == len(df)
    # Time written from BJD_TDB and sorted ascending.
    assert parsed["BJD_TDB"].is_monotonic_increasing
    # Normalized flux baseline ~ 1.0 (median magnitude maps to flux 1).
    assert abs(float(np.median(parsed["flux"])) - 1.0) < 0.05
    assert (parsed["flux_err"] > 0).all()


def test_exoclock_flux_from_magnitude_roundtrip(tmp_path):
    df = _synthetic_lc()
    out = export_exoclock(df, _meta(), tmp_path / "exoclock.csv")
    parsed = pd.read_csv(out, comment="#")

    # Invert flux -> delta-mag and compare against (mag - median mag).
    mags = df["diff_mag_corr"].to_numpy(float)
    mref = float(np.median(mags))
    expected_flux_sorted = np.sort(np.power(10.0, -0.4 * (mags - mref)))[::-1]
    got_flux_sorted = np.sort(parsed["flux"].to_numpy(float))[::-1]
    assert np.allclose(got_flux_sorted, expected_flux_sorted, atol=1e-5)


# ── ExoFOP ────────────────────────────────────────────────────────────────────

def test_exofop_header_columns_and_airmass(tmp_path):
    df = _synthetic_lc()
    out = export_exofop(df, _meta(), tmp_path / "exofop.csv")
    assert out.exists()
    parsed = pd.read_csv(out, comment="#")
    assert list(parsed.columns) == EXOFOP_COLUMNS == [
        "BJD_TDB",
        "rel_flux",
        "rel_flux_err",
        "airmass",
    ]
    assert len(parsed) == len(df)
    assert parsed["BJD_TDB"].is_monotonic_increasing
    assert math.isclose(
        float(parsed["airmass"].iloc[0]), float(df["airmass"].iloc[0]), abs_tol=1e-4
    )


# ── shared input validation ───────────────────────────────────────────────────

@pytest.mark.parametrize("fn", [export_aavso_extended, export_exoclock, export_exofop])
def test_empty_table_raises(tmp_path, fn):
    with pytest.raises((ValueError, TypeError)):
        fn(pd.DataFrame(), _meta(), tmp_path / "out.txt")


@pytest.mark.parametrize("fn", [export_aavso_extended, export_exoclock, export_exofop])
def test_missing_time_column_raises(tmp_path, fn):
    bad = pd.DataFrame({"mag": [12.0, 12.1], "filter": ["V", "V"]})
    with pytest.raises(ValueError):
        fn(bad, _meta(), tmp_path / "out.txt")


def test_falls_back_to_raw_when_no_corrected_mag(tmp_path):
    df = _synthetic_lc().drop(columns=["diff_mag_corr"])
    out = export_aavso_extended(df, _meta(), tmp_path / "aavso.txt")
    parsed = pd.read_csv(out, comment="#")
    assert math.isclose(
        float(parsed["MAG"].iloc[0]), float(df["diff_mag_raw"].iloc[0]), abs_tol=1e-4
    )
