from __future__ import annotations

from types import SimpleNamespace

import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

from apex.utils.gaia_catalog_service import GaiaCatalogService
from apex.utils.step_paths import step5_wcs_dir


def _contract_complete(df: pd.DataFrame) -> pd.DataFrame:
    """Fill in every contract column the fixture does not care about.

    A cached catalogue missing a contract column is deliberately treated as
    stale (it would silently disable a step10 quality cut), so a fixture that
    hard-codes a column list goes out of date the moment the contract grows.
    Deriving it here means these tests exercise cache *reuse*, which is what
    they are about, instead of failing on an unrelated schema change.
    """
    from apex.utils.gaia_columns import GAIA_COLUMNS

    out = df.copy()
    have = {str(c).strip().lower() for c in out.columns}
    for col in GAIA_COLUMNS:
        # Alias-aware: a fixture that deliberately uses a legacy spelling
        # (ra_deg for ra) already supplies that column, and adding the
        # canonical name alongside it would defeat the very test that checks
        # the legacy spelling is understood.
        spellings = {col.name.lower(), *(a.lower() for a in col.aliases)}
        if col.vizier:
            spellings.add(col.vizier.rsplit(".", 1)[-1].strip('"').lower())
        if spellings & have:
            continue
        out[col.name] = "" if col.name == "phot_variable_flag" else 1.0
    return out


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
        _contract_complete(
            pd.DataFrame(
                {
                    "source_id": [1, 2, 3],
                    "ra": [10.0, 10.001, 9.999],
                    "dec": [20.0, 20.001, 19.999],
                    "phot_g_mean_mag": [15.0, 18.5, 17.0],
                    "ruwe": [1.0, 1.1, 0.9],
                }
            )
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
        _contract_complete(
            pd.DataFrame(
                {
                    "source_id": [11, 12],
                    "ra_deg": [250.0, 250.002],
                    "dec_deg": [36.0, 35.998],
                    "phot_g_mean_mag": [13.0, 14.0],
                    "ruwe": [1.0, 1.2],
                }
            )
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
