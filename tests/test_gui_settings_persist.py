"""A setting a window offers has to survive pressing Save.

The chain is: a widget shows `getattr(P, name, literal)`, the user edits it,
the dialog writes it back onto `P`, and `save_toml` walks the key map writing
`getattr(P, attr)` into the file. Break any link and the window still says
"Parameters saved." while nothing is written — the value lives until the app
closes and then goes.

Reproduced on 2026-08-16 with the Extinction (Airmass Fit) tool:

    loaded, P has extfit/extinction knob : False
    file value                           : 3.0
    save_toml returned                   : True   <- dialog shows "Parameters saved."
    file value after save                : 3.0    <- the 9.9 never landed
    value after reload                   : attribute does not exist

Thirty settings were in that state, seventeen of them in that one tool. They
are wired now. This test is the guard: a name a window reads off `P` or writes
to it must either arrive (so the map persists it) or be listed below as
something that is not a setting.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from apex.config.parameter_map import CMD_TOML_KEY_MAP, LC_TOML_KEY_MAP
from apex.config.parameters_cmd import Parameters, read_params

REPO = Path(__file__).absolute().parents[1]
GUI = [p for p in (REPO / "apex/gui").rglob("*.py") if "__pycache__" not in p.parts]

READS = re.compile(r'getattr\(\s*(?:self\.)?(?:params\.)?P\s*,\s*"([a-z_][a-z0-9_]*)"')
WRITES = re.compile(
    r'(?:setattr\(\s*self\.params\.P\s*,\s*"([a-z_][a-z0-9_]*)"'
    r'|self\.params\.P\.([a-z_][a-z0-9_]*)\s*=)'
)

# Names a window touches on `P` that are deliberately not settings. Anything
# else must be reachable, or the window is offering a value it cannot keep.
NOT_SETTINGS = {
    # Runtime state the window fills in, never read from a file.
    "file_path_map": "airmass window builds this at run time",
    "readnoise": "Step 0 reads the measured value off the frames",
    # Legacy duplicates. `target.ra_deg` maps to `target_ra_deg`, which is what
    # the code actually reads; these bare names are the old spelling.
    "ra_deg": "superseded by target_ra_deg",
    "dec_deg": "superseded by target_dec_deg",
    # Per-window image cache sizes: tuning knobs with no dialog, sized from the
    # step's own working set rather than from a preference.
    "step3_fits_cache_size": "internal cache size, no dialog",
    "step4_fits_cache_size": "internal cache size, no dialog",
    "step8_fits_cache_size": "internal cache size, no dialog",
    # A real setting, deliberately left unreachable in CMD: LC's map gives it 5
    # (the config file's value) where the CMD code falls back to 10, so wiring it
    # would change which extinction fits CMD accepts. Needs a decision, not a
    # quiet move — see docs/audit/CONFIG_REACHABILITY.md.
    "extfit_min_points": "CMD 코드 10 vs 파일 5 — 배선하면 적합 결과가 바뀐다",
}


def _mode_of(path: Path) -> str:
    """Which modes open this window. Steps 1-7 and the tools open in both."""
    parts = str(path).replace("\\", "/")
    if "/workflow/cmd/" in parts:
        return "cmd"
    if "/workflow/lc/" in parts:
        return "lc"
    return "both"


def _gui_setting_names() -> dict[str, tuple[str, str]]:
    """attr -> (mode that opens the window, file it was found in)."""
    names: dict[str, tuple[str, str]] = {}
    for path in GUI:
        text = path.read_text(encoding="utf-8", errors="replace")
        for match in READS.finditer(text):
            names.setdefault(match.group(1), (_mode_of(path), path.name))
        for match in WRITES.finditer(text):
            names.setdefault(match.group(1) or match.group(2), (_mode_of(path), path.name))
    return names


@pytest.fixture(scope="module")
def blank(tmp_path_factory):
    config = tmp_path_factory.mktemp("ws") / "apex_config.json"
    config.write_text(json.dumps({"io": {"result_dir": ".", "data_dir": "."}}),
                      encoding="utf-8")
    return config


@pytest.mark.parametrize("mode", ["cmd", "lc"])
def test_every_setting_a_window_offers_can_be_saved(mode, tmp_path):
    """Per mode, not the union of both.

    Checking `CMD map | LC map` was too lenient and hid the biggest instance:
    Step 6 builds the master catalogue and is shared by both modes, but all
    seventeen of its settings sat in the LC-only map. Opened from CMD the window
    showed them, saved none of them, and reverted on the next launch.
    """
    from apex.config.parameters_lc import read_params as read_lc

    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({"io": {"result_dir": ".", "data_dir": "."}}),
                      encoding="utf-8")
    key_map = LC_TOML_KEY_MAP if mode == "lc" else CMD_TOML_KEY_MAP
    P = (read_lc if mode == "lc" else read_params)(config).P

    names = _gui_setting_names()
    assert len(names) > 250, "GUI 를 못 훑었다 — 정규식을 확인할 것"

    mapped = {row[1] for row in key_map}
    broken = sorted(
        name for name, (window_mode, _file) in names.items()
        if name not in NOT_SETTINGS
        and window_mode in (mode, "both")
        and not hasattr(P, name) and name not in mapped
    )
    assert not broken, (
        f"{mode.upper()} 모드에서 창이 저장할 수 없는 설정을 내놓는다 — 맵에 행을 "
        f"넣거나, 설정이 아니면 NOT_SETTINGS 에 이유와 함께 적을 것: "
        + ", ".join(f"{n} ({names[n][1]})" for n in broken)
    )


def test_the_not_settings_list_has_no_stale_entries(blank):
    """Something that became reachable should leave the exemption list."""
    names = _gui_setting_names()
    stale = sorted(n for n in NOT_SETTINGS if n not in names)
    assert not stale, f"GUI 가 더는 안 쓰는 이름이 목록에 남아 있다: {stale}"


def test_a_window_edit_survives_save_and_reload(tmp_path):
    """The round trip the Extinction tool used to fail."""
    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({
        "io": {"result_dir": str(tmp_path / "result"), "data_dir": str(tmp_path / "data")},
    }), encoding="utf-8")

    params = Parameters(config)
    assert params.P.extinction_snr_min == pytest.approx(10.0)

    params.P.extinction_snr_min = 9.9          # what the dialog does on Save
    assert params.save_toml() is True

    written = json.loads(config.read_text(encoding="utf-8"))
    assert written["extinction_fit"]["snr_min"] == pytest.approx(9.9)
    assert Parameters(config).P.extinction_snr_min == pytest.approx(9.9)
