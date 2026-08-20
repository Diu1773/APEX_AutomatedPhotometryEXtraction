"""Ranges to install with, a lock to reproduce with — and neither may rot.

`pyproject.toml` and `requirements.txt` keep version ranges on purpose: an
install should not fail because a dependency moved a patch release.
`requirements-lock.txt` is the other half — the exact versions the published
measurements were made on.

The reason both exist is one measurement. On 2026-08-17 five PSF magnitudes out
of 22,305 moved by 4.8e-05 between two runs of identical code, and the cause was
scipy 1.18.0 against 1.17.1. All five were already flagged bad and none reached
a figure, but a dependency that can move a number quietly is one that has to be
written down (D-011).

What can go wrong with this arrangement is that the two drift: a range narrows
past what the lock pins, or a core package is added to `pyproject` and never
locked. Both are silent — the install still works, and only a reproduction
attempt months later finds out.
"""

from __future__ import annotations

import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).absolute().parents[1]
LOCK = ROOT / "requirements-lock.txt"
PYPROJECT = ROOT / "pyproject.toml"


def _locked() -> dict[str, str]:
    out = {}
    for line in LOCK.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "==" not in line:
            continue
        name, _, version = line.partition("==")
        out[name.strip().lower().replace("_", "-")] = version.strip()
    return out


def _core_requirements() -> list[str]:
    text = PYPROJECT.read_text(encoding="utf-8")
    block = re.search(r"^dependencies\s*=\s*\[(.*?)^\]", text, re.S | re.M)
    assert block, "pyproject has no core dependencies block"
    return [m.group(1) for m in re.finditer(r'"([^"]+)"', block.group(1))]


def test_the_lock_exists_and_pins_exactly():
    assert LOCK.exists(), "requirements-lock.txt is missing"
    locked = _locked()
    assert len(locked) > 30, f"only {len(locked)} pins — did the freeze run?"
    for name, version in locked.items():
        assert re.match(r"^[\w.+!-]+$", version), f"{name} is pinned to {version!r}"


def test_the_lock_says_which_python_it_was_taken_on():
    """A pin set is only reproducible against the interpreter it came from."""
    head = LOCK.read_text(encoding="utf-8")[:2000]
    assert re.search(r"python\s+3\.\d+", head), (
        "the lock does not record the Python version it was captured on")
    assert "commit" in head, "the lock does not record the commit it was captured at"


def test_every_core_dependency_is_locked():
    """A package added to `pyproject` and never locked reproduces as whatever
    pip happens to pick — which is the situation the lock exists to end."""
    locked = _locked()
    missing = []
    for spec in _core_requirements():
        if ";" in spec:                     # environment marker: may not install here
            continue
        name = re.split(r"[<>=!~\[]", spec, 1)[0].strip().lower().replace("_", "-")
        if name and name not in locked:
            missing.append(name)
    assert not missing, (
        f"core dependencies with no pin in requirements-lock.txt: {missing} — "
        "regenerate it with `pip freeze` from the deployment venv")


def test_the_ranges_still_admit_the_locked_versions():
    """The two halves must agree. A range that narrowed past its own lock means
    the recorded environment can no longer be installed by the recorded rules.
    """
    packaging = pytest.importorskip("packaging.specifiers")
    from packaging.version import Version

    locked = _locked()
    violations = []
    for spec in _core_requirements():
        if ";" in spec:
            continue
        m = re.match(r"^([\w.\-]+)\s*(.*)$", spec.strip())
        if not m:
            continue
        name = m.group(1).lower().replace("_", "-")
        rng = m.group(2).strip()
        if name not in locked or not rng:
            continue
        try:
            spec_set = packaging.SpecifierSet(rng)
        except Exception:                               # noqa: BLE001
            continue
        if Version(locked[name]) not in spec_set:
            violations.append(f"{name}: locked {locked[name]}, range {rng}")
    assert not violations, (
        "pyproject ranges exclude the locked versions: " + "; ".join(violations))


def test_the_ranges_are_still_ranges():
    """The lock is what pins. If `pyproject` starts pinning too, a normal
    install breaks the moment one dependency ships a patch release."""
    pinned = [s for s in _core_requirements()
              if "==" in s and ";" not in s]
    assert not pinned, (
        f"pyproject pins exactly: {pinned} — pinning belongs in "
        "requirements-lock.txt so that installing APEX stays possible")
