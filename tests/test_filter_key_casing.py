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


# --- 사용자가 자기 규칙을 설정으로 줄 수 있어야 한다 --------------------------
#
# 기본 별칭표에 이름을 더하는 것은 임시방편이다. 관측소마다 필터를 자기 방식으로
# 적으므로, 새 관측소를 만날 때마다 코드를 고쳐야 한다면 그 설계가 틀린 것이다
# (사용자 지적, 2026-09-08). 설정 `[filters].aliases` 로 받는다.

from apex.utils.astro_utils import (
    filter_aliases_in_effect,
    normalize_filter_name,
    register_filter_aliases,
)


@pytest.fixture(autouse=True)
def _clear_user_aliases():
    """한 시험의 별칭이 다음 시험으로 새지 않게 한다."""
    register_filter_aliases({})
    yield
    register_filter_aliases({})


def test_a_site_can_name_its_own_filters():
    register_filter_aliases({"Bessell-V": "V", "F01": "B", "무명필터": "R"})
    assert normalize_filter_name("Bessell-V") == "V"
    assert normalize_filter_name("F01") == "B"
    assert normalize_filter_name("무명필터") == "R"


def test_user_aliases_beat_the_built_in_table():
    """관측소가 우리 표와 다른 뜻으로 쓰면 그쪽이 옳다."""
    assert normalize_filter_name("gp") == "g"          # 기본표
    register_filter_aliases({"gp": "V"})
    assert normalize_filter_name("gp") == "V"


def test_the_built_in_table_still_works_for_names_the_user_did_not_give():
    register_filter_aliases({"Bessell-V": "V"})
    assert normalize_filter_name("sdss_r") == "r"
    assert normalize_filter_name("halpha") == "Ha"


def test_case_does_not_matter_on_the_left():
    register_filter_aliases({"Bessell-V": "V"})
    for spelling in ("bessell-v", "BESSELL-V", "Bessell-V"):
        assert normalize_filter_name(spelling) == "V", spelling


def test_empty_values_are_dropped_not_applied():
    """설정에서 지우다 만 줄이 필터 이름을 빈칸으로 만들면 안 된다."""
    register_filter_aliases({"zs": "", "F01": "  ", "ok": "z"})
    assert normalize_filter_name("zs") == "zs"
    assert normalize_filter_name("F01") == "F01"
    assert normalize_filter_name("ok") == "z"


def test_what_is_in_effect_is_inspectable():
    register_filter_aliases({"Bessell-V": "V"})
    eff = filter_aliases_in_effect()
    assert eff["bessell-v"] == "V"     # 사용자 것
    assert eff["gp"] == "g"            # 기본표


def test_the_config_section_reaches_the_normalizer(tmp_path):
    """설정 파일 -> read_params -> normalize 까지 실제로 이어지는가."""
    import json as _json
    from pathlib import Path

    from apex.config.parameters_cmd import read_params

    base = _json.loads(Path("parameters.example.json").read_text(encoding="utf-8"))
    base["filters"] = {"aliases": {"Bessell-V": "V", "zs": "z"}}
    cfg = tmp_path / "apex_config.json"
    cfg.write_text(_json.dumps(base, ensure_ascii=False), encoding="utf-8")

    params = read_params(str(cfg))
    assert params.P.filter_aliases == {"Bessell-V": "V", "zs": "z"}
    assert normalize_filter_name("Bessell-V") == "V"
    assert normalize_filter_name("zs") == "z"
