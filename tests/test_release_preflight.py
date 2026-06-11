from __future__ import annotations

import subprocess

from deploy.verify_release import _check_no_untracked_package_sources


def test_release_preflight_rejects_untracked_package_source(tmp_path):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    source = tmp_path / "apex" / "new_module.py"
    source.parent.mkdir()
    source.write_text("VALUE = 1\n", encoding="utf-8")

    errors: list[str] = []
    _check_no_untracked_package_sources(tmp_path, errors)

    assert errors
    assert "apex/new_module.py" in errors[0]


def test_release_preflight_accepts_staged_package_source(tmp_path):
    subprocess.run(["git", "init", "-q"], cwd=tmp_path, check=True)
    source = tmp_path / "apex" / "new_module.py"
    source.parent.mkdir()
    source.write_text("VALUE = 1\n", encoding="utf-8")
    subprocess.run(["git", "add", "apex/new_module.py"], cwd=tmp_path, check=True)

    errors: list[str] = []
    _check_no_untracked_package_sources(tmp_path, errors)

    assert errors == []
