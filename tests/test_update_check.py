"""Tests for the Qt-free update-check version logic (no network)."""

from __future__ import annotations

import pytest

from apex.utils.update_check import is_newer, parse_version


@pytest.mark.parametrize("value,expected", [
    ("0.1.0", (0, 1, 0)),
    ("v1.2.3", (1, 2, 3)),
    ("V2.0", (2, 0)),
    ("1.2.3-rc1", (1, 2, 3)),
    ("1.4.0+build5", (1, 4, 0)),
    ("garbage", (0,)),
    ("", (0,)),
])
def test_parse_version(value, expected):
    assert parse_version(value) == expected


@pytest.mark.parametrize("remote,local,expected", [
    ("0.2.0", "0.1.0", True),
    ("0.1.1", "0.1.0", True),
    ("1.0.0", "0.9.9", True),
    ("0.1.0", "0.1.0", False),     # equal -> not newer
    ("0.1.0", "0.2.0", False),     # older -> not newer
    ("0.1", "0.1.0", False),       # 0.1 == 0.1.0 after padding
    ("v1.2.0", "1.1.9", True),     # tolerant of leading v
])
def test_is_newer(remote, local, expected):
    assert is_newer(remote, local) is expected
