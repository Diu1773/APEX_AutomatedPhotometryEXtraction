import sys
from pathlib import Path

import pandas as pd
import pytest

from apex.benchmark.iraf_crosscheck import (
    _match_iraf_to_apex,
    _parse_threshold_grid,
    add_iraf_calibrated_equivalent_columns,
    build_error_model_by_magnitude,
    build_fixed_comparison,
    load_apex_reference_table,
    parse_txdump,
    select_fixed_coordinate_sources,
    summarize_mag_comparison,
    windows_to_wsl_path,
    write_iraf_coords,
)


@pytest.mark.skipif(
    sys.platform != "win32",
    reason="windows_to_wsl_path maps Windows drive paths for WSL; only used on Windows",
)
def test_windows_to_wsl_path_converts_drive_path():
    assert windows_to_wsl_path(r"E:\data\frame.fit") == "/mnt/e/data/frame.fit"


def test_load_apex_reference_table_filters_bad_rows(tmp_path):
    path = tmp_path / "photometry_frame.tsv"
    pd.DataFrame(
        {
            "x": [1.0, 2.0, 3.0, 4.0],
            "y": [1.0, 2.0, 3.0, 4.0],
            "mag_inst": [-10.0, -11.0, -12.0, -13.0],
            "snr": [30.0, 5.0, 40.0, 50.0],
            "bad_phot_flag": [False, False, True, False],
            "is_saturated": [False, False, False, True],
        }
    ).to_csv(path, sep="\t", index=False)

    ref = load_apex_reference_table(path, min_snr=20.0)

    assert len(ref) == 1
    assert ref.loc[0, "x"] == pytest.approx(1.0)


def test_select_fixed_coordinate_sources_is_magnitude_stratified():
    ref = pd.DataFrame(
        {
            "x": range(100),
            "y": range(100),
            "mag_inst": [-14.0 + i * 0.04 for i in range(100)],
        }
    )

    selected = select_fixed_coordinate_sources(ref, max_sources=16, seed=1, bins=4)

    assert len(selected) == 16
    assert selected["iraf_id"].tolist() == list(range(1, 17))
    assert selected["mag_inst"].min() < -13.0
    assert selected["mag_inst"].max() > -11.0


def test_write_iraf_coords_preserves_iraf_id_order(tmp_path):
    ref = pd.DataFrame({"iraf_id": [2, 1], "x": [20.0, 10.0], "y": [21.0, 11.0]})
    path = write_iraf_coords(ref, tmp_path / "coords.coo")

    # Written values carry IRAF's 1-based pixel origin (+1 on both axes);
    # ordering by iraf_id is what this test pins. The origin conversion itself
    # is covered by tests/test_iraf_pixel_origin.py.
    assert path.read_text(encoding="ascii").splitlines() == [
        "11.000000 12.000000",
        "21.000000 22.000000",
    ]


def test_parse_txdump_reads_indef_as_nan(tmp_path):
    path = tmp_path / "phot.txt"
    path.write_text("1 10.1 11.2 15.0 0.02 100.0 5.0 40\n2 20 21 INDEF INDEF 90 4 35\n")

    table = parse_txdump(path)

    assert len(table) == 2
    assert table.loc[0, "iraf_mag"] == pytest.approx(15.0)
    assert pd.isna(table.loc[1, "iraf_mag"])


def test_fixed_comparison_reports_median_removed_scatter():
    reference = pd.DataFrame(
        {
            "iraf_id": [1, 2, 3],
            "x": [10.0, 20.0, 30.0],
            "y": [10.0, 20.0, 30.0],
            # IRAF-style count-rate convention: mag_inst = INSTRUMENTAL_ZMAG(25)
            # - 2.5*log10(flux_e/exptime), i.e. old raw values shifted by +25.
            "mag_inst": [15.0, 14.0, 13.0],
            "apcorr": [1.0, 1.0, 1.0],
            "gain_e_per_adu": [1.0, 1.0, 1.0],
        }
    )
    iraf = pd.DataFrame(
        {
            "iraf_id": [1, 2, 3],
            "iraf_x": [10.0, 20.0, 30.0],
            "iraf_y": [10.0, 20.1, 30.0],
            "iraf_mag": [15.1, 14.0, 13.0],
            "iraf_merr": [0.01, 0.01, 0.01],
            "iraf_msky": [0.0, 0.0, 0.0],
            "iraf_stdev": [1.0, 1.0, 1.0],
            "iraf_nsky": [10, 10, 10],
        }
    )

    comparison = build_fixed_comparison(reference, iraf, zmag=25.0)
    summary = summarize_mag_comparison(comparison)

    assert comparison["delta_mag_raw"].tolist() == pytest.approx([-0.1, 0.0, 0.0])
    assert summary["n_matched"] == 3
    assert summary["mad_sigma_delta_mag_centered"] == pytest.approx(0.0)


def test_fixed_comparison_aligns_iraf_to_frame_zeropoint():
    reference = pd.DataFrame(
        {
            "iraf_id": [1, 2],
            "x": [10.0, 20.0],
            "y": [10.0, 20.0],
            # Count-rate convention (+25 vs old raw); frame ZP shifts by -25 to
            # keep the calibrated magnitude (mag_inst + frame_zp) invariant.
            "mag_inst": [15.0, 14.0],
        }
    )
    iraf = pd.DataFrame(
        {
            "iraf_id": [1, 2],
            "iraf_x": [10.0, 20.0],
            "iraf_y": [10.0, 20.0],
            "iraf_mag": [14.8, 13.8],
        }
    )

    comparison = build_fixed_comparison(reference, iraf, zmag=25.0, frame_zeropoint=-2.5)
    summary = summarize_mag_comparison(comparison)

    assert comparison["apex_mag_cal_frame"].tolist() == pytest.approx([12.5, 11.5])
    assert summary["frame_zeropoint_mag"] == pytest.approx(-2.5)
    assert summary["iraf_to_apex_zp_offset"] == pytest.approx(-2.3)
    assert comparison["iraf_mag_aligned_to_apex"].tolist() == pytest.approx([12.5, 11.5])
    assert summary["mad_sigma_delta_mag_zp_aligned"] == pytest.approx(0.0)


def test_iraf_calibrated_equivalent_applies_gain_apcorr_and_zp():
    table = pd.DataFrame(
        {
            "iraf_mag": [20.0],
            # Count-rate convention: mag_inst already carries INSTRUMENTAL_ZMAG.
            "mag_inst": [14.247425011],
            "gain_e_per_adu": [2.0],
            "apcorr": [1.5],
            "apex_mag_cal_frame": [12.247425011],
        }
    )

    out = add_iraf_calibrated_equivalent_columns(table, zmag=25.0, frame_zeropoint=-2.0)

    # IRAF MAG=20 (zmag=25, itime=exptime) -> count rate 100 e-/s. Applying
    # gain=2 and apcorr=1.5 gives 300 e-/s, so the APEX-equivalent instrumental
    # magnitude is INSTRUMENTAL_ZMAG - 2.5*log10(300) = 25 - 6.192803 = 18.807197.
    assert out.loc[0, "iraf_mag_inst_apcorr_equiv"] == pytest.approx(18.807196863, abs=1e-6)
    assert out.loc[0, "iraf_mag_cal_apcorr_zp"] == pytest.approx(16.807196863, abs=1e-6)
    assert out.loc[0, "delta_mag_apcorr_zp"] == pytest.approx(-4.559771852, abs=1e-6)
    assert out.loc[0, "delta_mag_apcorr_zp_centered"] == pytest.approx(0.0)


def test_error_model_by_magnitude_summarizes_formal_errors_and_residuals():
    table = pd.DataFrame(
        {
            "apex_mag_cal_frame": [12.0, 12.1, 12.2, 12.3, 12.4, 13.0, 13.1, 13.2, 13.3, 13.4],
            "delta_mag_zp_aligned": [-0.01, 0.0, 0.01, 0.0, 0.02, -0.02, 0.0, 0.02, 0.0, 0.01],
            "mag_err": [0.03] * 5 + [0.05] * 5,
            "iraf_merr": [0.02] * 5 + [0.04] * 5,
        }
    )

    bins = build_error_model_by_magnitude(table, bin_width=1.0, min_count=5)

    assert len(bins) == 2
    assert bins.loc[0, "apex_mag_err_median"] == pytest.approx(0.03)
    assert bins.loc[0, "iraf_merr_median"] == pytest.approx(0.02)
    assert bins.loc[0, "residual_mad_sigma"] > 0


def test_match_iraf_to_apex_is_unique_and_radius_limited():
    apex = pd.DataFrame({"x": [10.0, 20.0], "y": [10.0, 20.0], "mag_inst": [-10.0, -11.0]})
    iraf = pd.DataFrame(
        {
            "iraf_x": [10.1, 10.2, 30.0],
            "iraf_y": [10.1, 10.2, 30.0],
            "iraf_mag": [15.0, 15.1, 16.0],
        }
    )

    matched = _match_iraf_to_apex(iraf, apex, radius_px=1.0)

    assert len(matched) == 1
    assert matched.loc[0, "match_distance_px"] == pytest.approx(2**0.5 * 0.1)


def test_parse_threshold_grid_sorts_descending_and_rejects_empty():
    assert _parse_threshold_grid([5, 12, 9, 9]) == [12.0, 9.0, 5.0]
    with pytest.raises(ValueError):
        _parse_threshold_grid([0, -1])
