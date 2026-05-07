from __future__ import annotations

from apex.utils.cache_utils import (
    DETECTION_CACHE_SCHEMA_VERSION,
    build_detection_cache_signature,
    detection_cache_signature_matches,
    normalize_detect_engine,
)


def test_normalize_detect_engine_aliases():
    assert normalize_detect_engine("sextractor") == "sep"
    assert normalize_detect_engine("photutils") == "segm"
    assert normalize_detect_engine("unknown") == "segm"


def test_detection_cache_signature_matches_current_file(tmp_path):
    source = tmp_path / "frame.fit"
    source.write_text("fits", encoding="utf-8")
    payload = build_detection_cache_signature(
        source,
        use_cropped=False,
        detect_engine="sourceextractor",
    )
    current = build_detection_cache_signature(
        source,
        use_cropped=False,
        detect_engine="sep",
    )

    assert payload["cache_schema"] == DETECTION_CACHE_SCHEMA_VERSION
    assert detection_cache_signature_matches(
        payload,
        current,
        min_schema=DETECTION_CACHE_SCHEMA_VERSION,
        current_engine="sep",
    )


def test_detection_cache_rejects_engine_and_schema_mismatch(tmp_path):
    source = tmp_path / "frame.fit"
    source.write_text("fits", encoding="utf-8")
    payload = build_detection_cache_signature(source, use_cropped=False, detect_engine="sep")
    current = build_detection_cache_signature(source, use_cropped=False, detect_engine="segm")

    assert not detection_cache_signature_matches(
        {**payload, "cache_schema": DETECTION_CACHE_SCHEMA_VERSION - 1},
        current,
        min_schema=DETECTION_CACHE_SCHEMA_VERSION,
        current_engine="sep",
    )
    assert not detection_cache_signature_matches(
        payload,
        current,
        min_schema=DETECTION_CACHE_SCHEMA_VERSION,
        current_engine="segm",
    )
