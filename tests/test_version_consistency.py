"""Guard: the version must stay single-sourced from apex.__version__.

deploy/version.txt drives the Windows build stamp and main.py shows it in the
GUI footer. This test fails if any of them drift from apex.__version__.
"""

from __future__ import annotations

from pathlib import Path

import apex

_ROOT = Path(__file__).resolve().parent.parent


def test_deploy_version_txt_matches_package():
    version_txt = _ROOT / "deploy" / "version.txt"
    assert version_txt.exists(), "deploy/version.txt is missing"
    assert version_txt.read_text(encoding="utf-8").strip() == apex.__version__


def test_pyproject_version_is_dynamic_from_package():
    try:
        import tomllib  # py3.11+
    except ModuleNotFoundError:  # pragma: no cover - py3.10
        import tomli as tomllib  # type: ignore

    data = tomllib.loads((_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    # Version is declared dynamic and sourced from apex.__version__.
    assert "version" in data["project"].get("dynamic", [])
    attr = data["tool"]["setuptools"]["dynamic"]["version"]["attr"]
    assert attr == "apex.__version__"
