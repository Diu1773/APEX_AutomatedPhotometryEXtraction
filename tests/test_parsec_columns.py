from __future__ import annotations

from apex.utils.parsec_columns import read_parsec_header


def test_read_parsec_header_detects_sdss_columns(tmp_path):
    path = tmp_path / "sdss.dat"
    path.write_text(
        "# Photometric system: SDSS <i>ugriz</i>\n"
        "# Zini MH logAge Mini mbolmag umag gmag rmag imag zmag\n"
        "0.01 0.0 9.0 1.0 4.0 5.0 4.0 3.0 2.0 1.0\n",
        encoding="utf-8",
    )

    info = read_parsec_header(path)

    assert info.photometric_system == "SDSS ugriz"
    assert info.band_columns == {"u": 5, "g": 6, "r": 7, "i": 8, "z": 9}
    assert "B" not in info.band_columns
    assert "V" not in info.band_columns


def test_read_parsec_header_preserves_bessell_band_case(tmp_path):
    path = tmp_path / "bessell.dat"
    path.write_text(
        "# Photometric system: UBVRIJHK\n"
        "# Zini MH logAge Mini mbolmag UXmag BXmag Bmag Vmag Rmag Imag\n"
        "0.01 0.0 9.0 1.0 4.0 5.0 5.0 4.0 3.0 2.0 1.0\n",
        encoding="utf-8",
    )

    info = read_parsec_header(path)

    assert info.band_columns == {"U": 5, "B": 7, "V": 8, "R": 9, "I": 10}
    assert "r" not in info.band_columns


def test_read_parsec_header_accepts_apex_normalized_basti_table(tmp_path):
    path = tmp_path / "basti_johnson.dat"
    path.write_text(
        "# Isochrone model: BaSTI-IAC updated models\n"
        "# Photometric system: Johnson-Cousins / Bessell (BaSTI-IAC)\n"
        "# Zini MH logAge Mini int_IMF Mass logL logTe logg label "
        "Umag BXmag Bmag Vmag Rmag Imag Jmag Hmag Kmag Lpmag Lmag Mmag\n"
        "0.01258 -0.08 9.0 1.0 0 1.0 0.0 3.8 0 1 "
        "4 3 3.1 2.6 2.3 2.0 1.8 1.7 1.6 1.5 1.5 1.5\n",
        encoding="ascii",
    )

    info = read_parsec_header(path)

    assert info.photometric_system == "Johnson-Cousins / Bessell (BaSTI-IAC)"
    assert info.band_columns == {"B": 12, "V": 13, "R": 14, "I": 15}
