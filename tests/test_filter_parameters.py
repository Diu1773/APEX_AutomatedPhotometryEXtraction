from __future__ import annotations

from apex.utils.filter_parameters import (
    build_active_filter_float_map,
    normalize_filter_float_map,
)


def test_normalize_filter_float_map_preserves_photometric_system_case():
    result = normalize_filter_float_map({
        "R": 3.4,
        "r": 4.1,
        "Rc": 3.6,
        "invalid": "not-a-number",
    })

    assert result == {"R": 3.6, "r": 4.1}


def test_build_active_filter_map_replaces_stale_sloan_keys_for_bvr_data():
    result = build_active_filter_float_map(
        {"g": 3.0, "r": 3.1, "i": 3.2},
        ["B", "V", "R"],
        3.2,
    )

    assert result == {"B": 3.2, "R": 3.2, "V": 3.2}


def test_build_active_filter_map_keeps_exact_canonical_overrides():
    result = build_active_filter_float_map(
        {"B": 3.3, "V": 3.4, "R": 3.5, "r": 4.2},
        ["B", "V", "R"],
        3.2,
    )

    assert result == {"B": 3.3, "R": 3.5, "V": 3.4}


def test_build_active_filter_map_keeps_config_when_scan_is_empty():
    result = build_active_filter_float_map({"Ha": 4.0}, [], 3.2)

    assert result == {"Ha": 4.0}
