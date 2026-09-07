import pandas as pd
from astropy.io import fits

from apex.analysis.merge.workspace_scan import (
    load_master_catalogs_by_filter,
    normalize_filter_key,
    scan_merge_input_workspace,
)
from apex.utils.astro_utils import get_filter_from_fits
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_lc import step8_selection_dir, step9_lc_dir


def _write_catalog(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_canonical_filter_key_preserves_johnson_and_sdss_case():
    assert normalize_filter_key("V") == "V"
    assert normalize_filter_key("v") == "V"
    assert normalize_filter_key("R") == "R"
    assert normalize_filter_key("r") == "r"


def test_workspace_scan_keeps_uppercase_johnson_filter_names(tmp_path):
    result_dir = tmp_path
    selection_dir = step8_selection_dir(result_dir)
    selection_dir.mkdir(parents=True)
    step7_forced_phot_dir(result_dir).mkdir(parents=True)
    step9_lc_dir(result_dir).mkdir(parents=True)

    _write_catalog(selection_dir / "master_catalog_V.tsv", [{"ID": 1, "source_id": 1001}])
    (selection_dir / "selection_V.json").write_text(
        '{"filter": "V", "target_id": 1, "target_source_id": 1001}',
        encoding="utf-8",
    )
    (step7_forced_phot_dir(result_dir) / "photometry_index.csv").write_text(
        "file,filter\nframe_V.fits,V\n",
        encoding="utf-8",
    )
    (step9_lc_dir(result_dir) / "lightcurve_ID1_raw.csv").write_text("JD,filter\n1,V\n", encoding="utf-8")

    catalogs = load_master_catalogs_by_filter(result_dir)
    scan = scan_merge_input_workspace(result_dir)

    assert "V" in catalogs
    assert "v" not in catalogs
    assert scan["filters"] == ["V"]


def test_get_filter_from_fits_returns_canonical_filter_case(tmp_path):
    path = tmp_path / "frame.fits"
    fits.PrimaryHDU(header=fits.Header({"FILTER": "v"})).writeto(path)

    assert get_filter_from_fits(path) == "V"


# --- LCO·MuSCAT 계열의 프라임 표기 -------------------------------------------
#
# 이 관측소들은 g' 를 `gp` 로 적는다. 별칭표에 `g'` 와 `gprime` 은 있었지만 이
# 철자가 없어서, MuSCAT3 자료를 넣으면 필터가 `gp` 로 그대로 흘러 SDSS 로
# 인식되지 않았다 (2026-09-07, 목표 둘 착수 중 발견).

import pytest

from apex.utils.common_helpers import (
    normalize_filter_key,
    photometric_system_label,
)


@pytest.mark.parametrize("raw,want", [
    ("gp", "g"), ("rp", "r"), ("ip", "i"), ("up", "u"), ("zp", "z"),
    ("GP", "g"), ("Rp", "r"),
])
def test_prime_p_spelling_maps_to_sdss(raw, want):
    assert normalize_filter_key(raw) == want


def test_z_short_is_not_folded_into_z():
    """z-short 는 z' 와 다른 대역이다. 같은 키로 묶으면 없는 등가를 주장한다."""
    assert normalize_filter_key("zs") == "zs"


def test_muscat3_bands_are_recognised_as_sdss():
    keys = [normalize_filter_key(f) for f in ("gp", "rp", "ip")]
    assert keys == ["g", "r", "i"]
    assert photometric_system_label(*keys) == "SDSS"


def test_the_existing_spellings_still_work():
    for raw in ("g'", "gprime", "SDSS-G", "sdss_g"):
        assert normalize_filter_key(raw) == "g", raw
