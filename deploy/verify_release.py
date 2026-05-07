"""Validate a local APEX Windows release bundle.

The checks are filesystem-only so they can run immediately after PyInstaller
without launching the PyQt GUI.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def _read_version(project_root: Path) -> str:
    version_path = project_root / "deploy" / "version.txt"
    try:
        version = version_path.read_text(encoding="utf-8").strip()
    except OSError:
        version = ""
    return version or "0.1.0"


def _resource_root(bundle_dir: Path) -> Path:
    internal = bundle_dir / "_internal"
    return internal if internal.is_dir() else bundle_dir


def _check_file(path: Path, errors: list[str]) -> None:
    if not path.is_file():
        errors.append(f"Missing file: {path}")


def _check_dir(path: Path, errors: list[str]) -> None:
    if not path.is_dir():
        errors.append(f"Missing directory: {path}")


def verify(project_root: Path) -> list[str]:
    project_root = project_root.resolve()
    version = _read_version(project_root)
    bundle_dir = project_root / "dist" / "APEX"
    resource_root = _resource_root(bundle_dir)
    setup_dir = project_root / "release" / "Setup"

    errors: list[str] = []
    _check_dir(bundle_dir, errors)
    _check_file(bundle_dir / "APEX.exe", errors)
    _check_file(resource_root / "parameters.example.toml", errors)
    _check_file(resource_root / "README.md", errors)
    _check_file(resource_root / "apex" / "resources" / "logo_base.svg", errors)
    _check_file(resource_root / "apex" / "resources" / "logo_cmd.svg", errors)
    _check_file(resource_root / "apex" / "resources" / "logo_lc.svg", errors)
    _check_file(resource_root / "apex" / "resources" / "apex.ico", errors)

    _check_dir(setup_dir, errors)
    _check_file(setup_dir / "setup.exe", errors)
    _check_file(setup_dir / f"APEX-Portable-{version}-x64.zip", errors)
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="APEX repository root.",
    )
    args = parser.parse_args(argv)

    errors = verify(args.project_root)
    if errors:
        for err in errors:
            print(f"[ERROR] {err}")
        return 1
    print("[OK] Release bundle looks complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
