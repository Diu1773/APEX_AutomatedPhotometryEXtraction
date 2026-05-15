"""Validate APEX release inputs and a local Windows release bundle."""

from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
import traceback
from pathlib import Path


SOURCE_IMPORTS = (
    "apex.gui.workflow.ui_helpers",
    "apex.gui.main_window",
    "apex.gui.workflow.step6_ref_build",
    "apex.gui.workflow.step7_forced_aperture_phot",
    "apex.gui.workflow.cmd.step8_psf_photometry",
    "astroquery.gaia",
    "astroquery.simbad",
    "astroquery.utils.tap.core",
    "certifi",
)

REQUIRED_SOURCE_FILES = (
    "main.py",
    "parameters.example.toml",
    "README.md",
    "apex/gui/workflow/ui_helpers.py",
    "apex/resources/logo_base.svg",
    "apex/resources/logo_cmd.svg",
    "apex/resources/logo_lc.svg",
)


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


def verify_source(project_root: Path) -> list[str]:
    project_root = project_root.resolve()
    errors: list[str] = []
    for rel_path in REQUIRED_SOURCE_FILES:
        _check_file(project_root / rel_path, errors)

    root_text = str(project_root)
    if root_text not in sys.path:
        sys.path.insert(0, root_text)

    for module_name in SOURCE_IMPORTS:
        try:
            importlib.import_module(module_name)
        except Exception:
            errors.append(
                f"Source import failed: {module_name}\n{traceback.format_exc()}"
            )
    return errors


def _check_exe_smoke(exe_path: Path, errors: list[str]) -> None:
    if not exe_path.is_file():
        return
    try:
        completed = subprocess.run(
            [str(exe_path), "--smoke"],
            cwd=str(exe_path.parent),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=90,
        )
    except Exception as exc:
        errors.append(f"Smoke test could not start: {exe_path} ({exc})")
        return
    if completed.returncode != 0:
        detail = (completed.stdout or "") + (completed.stderr or "")
        detail = detail.strip()
        if len(detail) > 4000:
            detail = detail[:4000] + "\n... output truncated ..."
        errors.append(
            f"Smoke test failed: {exe_path} returned {completed.returncode}"
            + (f"\n{detail}" if detail else "")
        )


def verify(project_root: Path, *, smoke: bool = True) -> list[str]:
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
    if smoke:
        _check_exe_smoke(bundle_dir / "APEX.exe", errors)
    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent,
        help="APEX repository root.",
    )
    parser.add_argument(
        "--source-only",
        action="store_true",
        help="Validate build-critical source files and imports only.",
    )
    parser.add_argument(
        "--skip-smoke",
        action="store_true",
        help="Skip APEX.exe --smoke when validating a built bundle.",
    )
    args = parser.parse_args(argv)

    errors = (
        verify_source(args.project_root)
        if args.source_only
        else verify(args.project_root, smoke=not args.skip_smoke)
    )
    if errors:
        for err in errors:
            print(f"[ERROR] {err}")
        return 1
    if args.source_only:
        print("[OK] Release source preflight passed.")
    else:
        print("[OK] Release bundle looks complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
