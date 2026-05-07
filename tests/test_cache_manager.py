from __future__ import annotations

import json

from apex.core.cache_manager import (
    CACHE_MANIFEST_VERSION,
    CacheInputSignature,
    CacheManifest,
    StepCacheManager,
    safe_cache_key,
)


def test_safe_cache_key_handles_paths_and_empty_values():
    assert safe_cache_key("frame 01/blue.fit").startswith("frame_01_blue.fit")
    assert safe_cache_key("   ")


def test_cache_input_signature_detects_size_and_mtime_changes(tmp_path):
    source = tmp_path / "frame.fit"
    source.write_text("abc", encoding="utf-8")
    saved = CacheInputSignature.from_path(source)

    source.write_text("abcd", encoding="utf-8")
    current = CacheInputSignature.from_path(source)

    assert not saved.matches(current)


def test_cache_manifest_roundtrip():
    manifest = CacheManifest(
        step_id="step4_detection",
        cache_schema_version=4,
        parameter_hash="abc123",
        input_files=[CacheInputSignature(path="x.fit", size=10, mtime_ns=20)],
        payload_paths={"json": "detect_x.json"},
        dependency_versions={"sep": "1.4"},
        extra={"engine": "sep"},
    )

    data = manifest.to_dict()
    loaded = CacheManifest.from_dict(data)

    assert loaded.manifest_version == CACHE_MANIFEST_VERSION
    assert loaded.step_id == "step4_detection"
    assert loaded.input_files[0].size == 10
    assert loaded.extra["engine"] == "sep"


def test_step_cache_manager_writes_and_validates_manifest(tmp_path):
    source = tmp_path / "frame.fit"
    payload = tmp_path / "detect_frame.json"
    source.write_text("fits", encoding="utf-8")
    payload.write_text("{}", encoding="utf-8")

    manager = StepCacheManager(
        tmp_path / "cache",
        step_id="step4_detection",
        cache_schema_version=4,
        parameter_hash="p1",
    )
    manifest = manager.build_manifest(
        input_paths=[source],
        payload_paths={"json": payload},
        dependency_versions={"test": "1"},
        extra={"engine": "sep"},
    )
    manifest_path = manager.write_manifest("frame.fit", manifest)

    assert manifest_path.exists()
    result = manager.validate_key(
        "frame.fit",
        input_paths=[source],
        required_payloads=["json"],
    )
    assert result.valid
    assert result.reason == "ok"


def test_step_cache_manager_reports_invalid_reasons(tmp_path):
    source = tmp_path / "frame.fit"
    payload = tmp_path / "detect_frame.json"
    source.write_text("fits", encoding="utf-8")
    payload.write_text("{}", encoding="utf-8")

    manager = StepCacheManager(tmp_path / "cache", "step4_detection", 4, parameter_hash="p1")
    manifest = manager.build_manifest(input_paths=[source], payload_paths={"json": payload})

    wrong_step = StepCacheManager(tmp_path / "cache", "step5_photometry", 4, parameter_hash="p1")
    assert wrong_step.validate_manifest(manifest, input_paths=[source]).reason == "step_id_mismatch"

    wrong_schema = StepCacheManager(tmp_path / "cache", "step4_detection", 5, parameter_hash="p1")
    assert wrong_schema.validate_manifest(manifest, input_paths=[source]).reason == "cache_schema_mismatch"

    wrong_hash = StepCacheManager(tmp_path / "cache", "step4_detection", 4, parameter_hash="p2")
    assert wrong_hash.validate_manifest(manifest, input_paths=[source]).reason == "parameter_hash_mismatch"

    payload.unlink()
    assert manager.validate_manifest(manifest, input_paths=[source], required_payloads=["json"]).reason == "payload_not_found:json"


def test_manifest_file_is_plain_json(tmp_path):
    source = tmp_path / "frame.fit"
    source.write_text("fits", encoding="utf-8")
    manager = StepCacheManager(tmp_path / "cache", "step4_detection", 4)
    path = manager.write_manifest("frame.fit", manager.build_manifest(input_paths=[source]))

    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["manifest_version"] == CACHE_MANIFEST_VERSION
    assert data["step_id"] == "step4_detection"
