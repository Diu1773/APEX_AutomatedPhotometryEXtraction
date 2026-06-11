"""Validate APEX release inputs and a local Windows release bundle."""

from __future__ import annotations

import argparse
import importlib
import subprocess
import sys
import traceback
from pathlib import Path


CORE_SOURCE_IMPORTS = (
    "apex.gui.workflow.ui_helpers",
    "apex.gui.main_window",
    "apex.gui.workflow.target_resolver",
    "apex.gui.workflow.step6_ref_build",
    "apex.gui.workflow.step7_forced_aperture_phot",
    "apex.gui.workflow.cmd.step8_psf_photometry",
    "astroquery.gaia",
    "astroquery.simbad",
    "astroquery.utils.tap.core",
    "certifi",
)

REQUIRED_SOURCE_FILES = (
    ".gitignore",
    "main.py",
    "parameters.example.toml",
    "README.md",
    "apex/gui/workflow/ui_helpers.py",
    "apex/resources/logo_base.svg",
    "apex/resources/logo_cmd.svg",
    "apex/resources/logo_lc.svg",
)

REQUIRED_GITIGNORE_PATTERNS = (
    "parameters.toml",
    "apex/.state/",
    "isochrone/",
    ".pytest_cache/",
    "*.log",
)

RUNTIME_TRACKED_EXACT = {
    "parameters.toml",
}
RUNTIME_TRACKED_PREFIXES = (
    "apex/.state/",
    "isochrone/",
    ".pytest_cache/",
    "build/",
    "dist/",
    "release/Setup/",
)
RUNTIME_TRACKED_SUFFIXES = (
    ".log",
    ".pyc",
    ".pyo",
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


def _check_gitignore(project_root: Path, errors: list[str]) -> None:
    path = project_root / ".gitignore"
    try:
        patterns = {
            line.strip()
            for line in path.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        }
    except OSError as exc:
        errors.append(f"Could not read .gitignore: {exc}")
        return
    missing = [pattern for pattern in REQUIRED_GITIGNORE_PATTERNS if pattern not in patterns]
    if missing:
        errors.append("Missing runtime ignore rules in .gitignore: " + ", ".join(missing))


def _check_no_tracked_runtime_files(project_root: Path, errors: list[str]) -> None:
    try:
        completed = subprocess.run(
            ["git", "ls-files"],
            cwd=str(project_root),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
        )
    except Exception:
        return
    if completed.returncode != 0:
        return

    offenders: list[str] = []
    for rel in completed.stdout.splitlines():
        rel = rel.strip().replace("\\", "/")
        if not rel:
            continue
        if rel in RUNTIME_TRACKED_EXACT:
            offenders.append(rel)
        elif rel.startswith(RUNTIME_TRACKED_PREFIXES):
            offenders.append(rel)
        elif rel.endswith(RUNTIME_TRACKED_SUFFIXES):
            offenders.append(rel)
    if offenders:
        preview = ", ".join(offenders[:12])
        if len(offenders) > 12:
            preview += f", ... (+{len(offenders) - 12} more)"
        errors.append("Runtime/build artifacts are tracked by git: " + preview)


def _check_no_untracked_package_sources(project_root: Path, errors: list[str]) -> None:
    """Reject releases that depend on Python package files absent from Git."""
    try:
        completed = subprocess.run(
            ["git", "ls-files", "--others", "--exclude-standard", "--", "apex"],
            cwd=str(project_root),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=20,
        )
    except Exception:
        return
    if completed.returncode != 0:
        return

    offenders = sorted(
        rel.strip().replace("\\", "/")
        for rel in completed.stdout.splitlines()
        if rel.strip().lower().endswith((".py", ".pyi"))
    )
    if offenders:
        preview = ", ".join(offenders[:12])
        if len(offenders) > 12:
            preview += f", ... (+{len(offenders) - 12} more)"
        errors.append(
            "Untracked Python package sources would be missing from a clean checkout: "
            + preview
        )


def verify_source(project_root: Path) -> list[str]:
    project_root = project_root.resolve()
    errors: list[str] = []
    for rel_path in REQUIRED_SOURCE_FILES:
        _check_file(project_root / rel_path, errors)
    _check_gitignore(project_root, errors)
    _check_no_tracked_runtime_files(project_root, errors)
    _check_no_untracked_package_sources(project_root, errors)

    root_text = str(project_root)
    if root_text not in sys.path:
        sys.path.insert(0, root_text)

    source_imports = list(CORE_SOURCE_IMPORTS)
    try:
        from apex.gui.tools.registry import iter_tool_modules
        source_imports.extend(iter_tool_modules())
    except Exception:
        errors.append(
            "Source import failed: apex.gui.tools.registry\n"
            + traceback.format_exc()
        )

    for module_name in source_imports:
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
            # Cold first launch of the frozen onedir bundle imports the whole
            # scientific + GUI stack (numpy/scipy/astropy/photutils/matplotlib/
            # PyQt5/astroquery) and builds the matplotlib font cache and astropy
            # IERS data on a CI runner with no warm cache. 90 s was too tight and
            # tripped a spurious timeout even though the import itself succeeds;
            # give cold runners real headroom.
            timeout=240,
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
    _check_file(setup_dir / f"setup-APEX-{version}.exe", errors)
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
