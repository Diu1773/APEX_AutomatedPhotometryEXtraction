from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

from apex.utils.gaia_catalog_service import GaiaCatalogService
from apex.utils.step_paths import step5_wcs_dir


def _params():
    return SimpleNamespace(
        P=SimpleNamespace(
            gaia_mag_max=20.0,
            gaia_wcs_mag_max=18.0,
            gaia_retry=1,
            gaia_backoff_s=0.0,
            gaia_allow_no_cache=True,
        )
    )


def test_gaia_catalog_service_uses_valid_cache_and_filters_mag(tmp_path):
    result_dir = tmp_path / "result"
    s5 = step5_wcs_dir(result_dir)
    s5.mkdir(parents=True)
    Table.from_pandas(
        pd.DataFrame(
            {
                "source_id": [1, 2, 3],
                "ra": [10.0, 10.001, 9.999],
                "dec": [20.0, 20.001, 19.999],
                "phot_g_mean_mag": [15.0, 18.5, 17.0],
                "phot_variable_flag": ["", "", ""],
            }
        )
    ).write(s5 / "gaia_fov.ecsv", format="ascii.ecsv")
    (s5 / "gaia_fov_meta.json").write_text(
        (
            '{'
            '"center_ra_deg": 10.0, '
            '"center_dec_deg": 20.0, '
            '"radius_deg": 0.5, '
            '"mag_max": 20.0, '
            '"n_stars": 3, '
            '"gaia_source": "test"'
            '}'
        ),
        encoding="utf-8",
    )

    service = GaiaCatalogService(_params(), result_dir)
    df, source = service.load_or_query(SkyCoord(10.0 * u.deg, 20.0 * u.deg), 0.3)

    assert source == "cache"
    assert list(df["source_id"]) == [1, 3]
    assert set(["ra", "dec", "phot_g_mean_mag"]) <= set(df.columns)


def test_gaia_catalog_service_normalizes_legacy_ra_dec_cache(tmp_path):
    result_dir = tmp_path / "result"
    s5 = step5_wcs_dir(result_dir)
    s5.mkdir(parents=True)
    Table.from_pandas(
        pd.DataFrame(
            {
                "source_id": [11, 12],
                "ra_deg": [250.0, 250.002],
                "dec_deg": [36.0, 35.998],
                "phot_g_mean_mag": [13.0, 14.0],
                "phot_variable_flag": ["", ""],
            }
        )
    ).write(s5 / "gaia_fov.ecsv", format="ascii.ecsv")
    (s5 / "gaia_fov_meta.json").write_text(
        (
            '{'
            '"center_ra_deg": 250.0, '
            '"center_dec_deg": 36.0, '
            '"radius_deg": 0.5, '
            '"mag_max": 18.0, '
            '"n_stars": 2, '
            '"gaia_source": "legacy"'
            '}'
        ),
        encoding="utf-8",
    )

    service = GaiaCatalogService(_params(), result_dir)
    df, source = service.load_or_query(SkyCoord(250.0 * u.deg, 36.0 * u.deg), 0.3)

    assert source == "cache"
    assert list(df["source_id"]) == [11, 12]
    assert "ra" in df.columns
    assert "dec" in df.columns
