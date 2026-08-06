"""Unit tests for apex.config.config_io — JSON authority + TOML migration."""

from __future__ import annotations

import json
import os
import time
import warnings

import pytest

from apex.config.config_io import (
    CONFIG_BASENAME,
    load_config_data,
    migrate_config_path,
    resolve_config_path,
    save_config_data,
)

TOML_BODY = """
schema_version = 6

[io]
data_dir = "E:\\\\somewhere\\\\sci"
result_dir = "E:\\\\somewhere\\\\result"

[target]
name = "TESTOBJ"
ra_deg = 132.8
"""


def _write_toml(path, body=TOML_BODY):
    path.write_text(body, encoding="utf-8")
    return path


def test_resolve_mapping_rules(tmp_path):
    assert resolve_config_path(tmp_path) == tmp_path / CONFIG_BASENAME
    assert resolve_config_path(tmp_path / "apex_config.json").name == "apex_config.json"
    assert resolve_config_path(tmp_path / "parameters.toml").name == "apex_config.json"
    assert (resolve_config_path(tmp_path / "parameters_result_psf.toml").name
            == "apex_config_result_psf.json")
    assert resolve_config_path(tmp_path / "custom.toml").name == "custom.json"


def test_migrates_legacy_toml_once(tmp_path):
    toml = _write_toml(tmp_path / "parameters.toml")
    data, json_path = load_config_data(toml)
    assert json_path == tmp_path / CONFIG_BASENAME
    assert json_path.exists()
    assert data["target"]["name"] == "TESTOBJ"
    assert data["_meta"]["migrated_from"] == "parameters.toml"
    # 두 번째 로드는 JSON 에서 (TOML 을 다시 읽지 않음)
    on_disk = json.loads(json_path.read_text(encoding="utf-8"))
    on_disk["target"]["name"] = "EDITED_IN_JSON"
    save_config_data(json_path, on_disk)
    data2, _ = load_config_data(toml)
    assert data2["target"]["name"] == "EDITED_IN_JSON"


def test_newer_toml_is_ignored_with_warning(tmp_path):
    toml = _write_toml(tmp_path / "parameters.toml")
    _, json_path = load_config_data(toml)
    past = time.time() - 100
    os.utime(json_path, (past, past))     # JSON 을 과거로
    _write_toml(toml)                     # TOML 을 지금으로
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        data, _ = load_config_data(toml)
    assert any("IGNORED" in str(w.message) for w in caught)
    assert data["target"]["name"] == "TESTOBJ"  # JSON 내용 유지


def test_missing_everything_yields_empty(tmp_path):
    data, json_path = load_config_data(tmp_path / "parameters.toml")
    assert data == {}
    assert not json_path.exists()   # 빈 설정은 파일을 만들지 않는다


def test_directory_request_uses_canonical_name(tmp_path):
    _write_toml(tmp_path / "parameters.toml")
    json_path = migrate_config_path(tmp_path)
    assert json_path == tmp_path / CONFIG_BASENAME
    assert json_path.exists()


def test_save_is_atomic_and_utf8(tmp_path):
    target = tmp_path / CONFIG_BASENAME
    assert save_config_data(target, {"한글": "값", "n": 1})
    loaded = json.loads(target.read_text(encoding="utf-8"))
    assert loaded == {"한글": "값", "n": 1}
    leftovers = [p for p in tmp_path.iterdir() if p.suffix == ".tmp"]
    assert not leftovers


def test_variant_tomls_do_not_collide(tmp_path):
    a = _write_toml(tmp_path / "parameters.toml")
    b = _write_toml(tmp_path / "parameters_variant.toml",
                    TOML_BODY.replace("TESTOBJ", "VARIANT"))
    _, ja = load_config_data(a)
    _, jb = load_config_data(b)
    assert ja != jb
    assert json.loads(ja.read_text(encoding="utf-8"))["target"]["name"] == "TESTOBJ"
    assert json.loads(jb.read_text(encoding="utf-8"))["target"]["name"] == "VARIANT"
