"""APEX headless command-line interface.

This is the scriptable front end for APEX. It deliberately avoids importing
PyQt5 (or any GUI module) at import time so that ``apex doctor`` / ``apex
config`` / ``apex version`` run on headless servers, in CI, and inside the
PyInstaller smoke test without a display.

Commands
--------
    apex version                 Print the APEX version and exit.
    apex doctor [--network]      Diagnose the runtime: Python, dependencies,
                                 external WCS solvers, and (optionally) network
                                 reachability to Gaia/SIMBAD.
    apex config init             Create parameters.toml from the bundled example.
    apex config path             Print the resolved parameters.toml path.
    apex config show             Print the active parameters.toml.
    apex gui [--mode cmd|lc]     Launch the desktop GUI (imports PyQt5 lazily).

The intent is for this module to grow an ``apex run`` subcommand once the
per-step science logic is decoupled from the Qt workers (see the deployment
roadmap). Today it owns the install-friction surface: the single biggest
adoption barrier for a new photometry tool is "it would not start".
"""

from __future__ import annotations

import argparse
import importlib
import os
import shutil
import sys
from pathlib import Path
from typing import Optional

from apex import __version__

# ── Result formatting ────────────────────────────────────────────────────────

_OK = "[ OK ]"
_WARN = "[WARN]"
_FAIL = "[FAIL]"

# Independent photometry engines for `apex validate --suite crosscheck`. Kept in
# sync with apex.benchmark.photometry_crosscheck.REFERENCE_CHOICES but declared
# here so the CLI stays import-light (no science imports at arg-parse time).
REFERENCE_CHOICES = ("sep", "iraf", "photutils")


def _repo_root() -> Path:
    """Best-effort location of the working/install directory.

    Frozen builds run beside the executable; source checkouts run from the
    repository root (the parent of the ``apex`` package).
    """
    if getattr(sys, "frozen", False):
        return Path(sys.executable).parent
    return Path(__file__).resolve().parent.parent


def _example_toml() -> Optional[Path]:
    if getattr(sys, "frozen", False):
        candidate = Path(getattr(sys, "_MEIPASS", _repo_root())) / "parameters.example.toml"
    else:
        candidate = _repo_root() / "parameters.example.toml"
    return candidate if candidate.exists() else None


def _params_path() -> Path:
    return _repo_root() / "parameters.toml"


# ── version ──────────────────────────────────────────────────────────────────

def _cmd_version(_args: argparse.Namespace) -> int:
    print(f"APEX {__version__}")
    print(f"Python {sys.version.split()[0]} ({sys.platform})")
    return 0


# ── doctor ───────────────────────────────────────────────────────────────────

# (import name, friendly label, required?)
_CORE_DEPS = [
    ("numpy", "numpy", True),
    ("scipy", "scipy", True),
    ("pandas", "pandas", True),
    ("matplotlib", "matplotlib", True),
    ("astropy", "astropy", True),
    ("photutils", "photutils", True),
    ("astroquery", "astroquery", True),
    ("sep", "sep", True),
    ("bottleneck", "bottleneck", True),
    ("pydantic", "pydantic", True),
]

# PyQt5/PyOpenGL and the GUI-tool science deps live in the `gui` extra, so they
# are OPTIONAL for the headless core (`pip install apex-photometry`). The
# desktop app needs them: `pip install apex-photometry[gui]`.
_OPTIONAL_DEPS = [
    ("PyQt5.QtCore", "PyQt5 (desktop GUI; install '[gui]')", False),
    ("OpenGL", "PyOpenGL (3D viewers; '[gui]')", False),
    ("batman", "batman-package (transit models; '[gui]')", False),
    ("emcee", "emcee (MCMC; '[gui]')", False),
    ("PIL", "Pillow", False),
]


def _module_version(import_name: str) -> str:
    top = import_name.split(".")[0]
    try:
        from importlib.metadata import version as _v  # py3.8+

        return _v(top)
    except Exception:
        pass
    try:
        mod = importlib.import_module(import_name)
        return str(getattr(mod, "__version__", "?"))
    except Exception:
        return "?"


def _check_imports(deps, results: list) -> None:
    for import_name, label, required in deps:
        try:
            importlib.import_module(import_name)
        except Exception as exc:  # noqa: BLE001 - report any failure
            tag = _FAIL if required else _WARN
            note = "required" if required else "optional"
            results.append((tag, f"{label}: not importable ({note}) - {exc}"))
            continue
        results.append((_OK, f"{label} {_module_version(import_name)}"))


def _check_python(results: list) -> None:
    major, minor = sys.version_info[:2]
    if (major, minor) >= (3, 10):
        results.append((_OK, f"Python {major}.{minor} (>= 3.10)"))
    else:
        results.append((_FAIL, f"Python {major}.{minor} is too old; APEX needs >= 3.10"))


def _check_external_solvers(results: list) -> None:
    # solve-field (astrometry.net)
    solve_field = shutil.which("solve-field")
    if solve_field:
        results.append((_OK, f"astrometry.net solve-field on PATH: {solve_field}"))
    else:
        results.append((_WARN, "astrometry.net solve-field not on PATH (optional WCS solver)"))

    # ASTAP - read the configured exe from parameters.toml when present.
    astap_exe = None
    pp = _params_path()
    if pp.exists():
        try:
            from apex.config.parameters_cmd import read_params

            astap_exe = getattr(read_params(pp), "astap_exe", None)
        except Exception:
            astap_exe = None
    candidate = astap_exe or "astap_cli.exe"
    resolved = shutil.which(candidate) or (candidate if Path(candidate).exists() else None)
    if resolved:
        results.append((_OK, f"ASTAP solver found: {resolved}"))
    else:
        results.append((_WARN, f"ASTAP not found ('{candidate}'); set [wcs].astap_exe (optional)"))


def _check_config(results: list) -> None:
    pp = _params_path()
    if pp.exists():
        results.append((_OK, f"parameters.toml present: {pp}"))
    elif _example_toml():
        results.append((_WARN, "parameters.toml missing - run 'apex config init' to create it"))
    else:
        results.append((_FAIL, "Neither parameters.toml nor parameters.example.toml found"))


def _check_network(results: list) -> None:
    try:
        import requests
    except Exception:
        results.append((_WARN, "network check skipped: requests not importable"))
        return
    endpoints = [
        ("Gaia TAP", "https://gea.esac.esa.int/tap-server/tap/availability"),
        ("SIMBAD", "https://simbad.cds.unistra.fr/simbad/"),
    ]
    for label, url in endpoints:
        try:
            resp = requests.get(url, timeout=8)
            if resp.status_code < 500:
                results.append((_OK, f"{label} reachable (HTTP {resp.status_code})"))
            else:
                results.append((_WARN, f"{label} returned HTTP {resp.status_code}"))
        except Exception as exc:  # noqa: BLE001
            results.append((_WARN, f"{label} unreachable - {exc}"))


def _cmd_doctor(args: argparse.Namespace) -> int:
    results: list = []
    print("APEX doctor - runtime diagnostics\n")
    _check_python(results)
    print("Core dependencies:")
    core: list = []
    _check_imports(_CORE_DEPS, core)
    print("Optional dependencies:")
    optional: list = []
    _check_imports(_OPTIONAL_DEPS, optional)
    print("External WCS solvers:")
    solvers: list = []
    _check_external_solvers(solvers)
    print("Configuration:")
    config: list = []
    _check_config(config)
    network: list = []
    if args.network:
        print("Network:")
        _check_network(network)

    # Single ordered render so failures are easy to scan.
    all_results = results + core + optional + solvers + config + network
    for tag, msg in all_results:
        print(f"  {tag} {msg}")

    n_fail = sum(1 for tag, _ in all_results if tag == _FAIL)
    n_warn = sum(1 for tag, _ in all_results if tag == _WARN)
    print()
    if n_fail:
        print(f"{_FAIL} {n_fail} blocking problem(s), {n_warn} warning(s).")
        return 1
    if n_warn:
        print(f"{_WARN} 0 blocking problems, {n_warn} warning(s) - APEX should run.")
        return 0
    print(f"{_OK} All checks passed.")
    return 0


# ── config ───────────────────────────────────────────────────────────────────

def _cmd_config(args: argparse.Namespace) -> int:
    action = args.config_action
    pp = _params_path()
    if action == "path":
        print(pp)
        return 0
    if action == "init":
        if pp.exists():
            print(f"{_WARN} parameters.toml already exists: {pp}")
            return 0
        example = _example_toml()
        if not example:
            print(f"{_FAIL} parameters.example.toml not found; cannot initialize.")
            return 1
        shutil.copy(example, pp)
        print(f"{_OK} Created {pp} from {example.name}")
        return 0
    if action == "show":
        if not pp.exists():
            print(f"{_FAIL} parameters.toml not found; run 'apex config init' first.")
            return 1
        print(pp.read_text(encoding="utf-8"))
        return 0
    return 2


# ── gui ──────────────────────────────────────────────────────────────────────

def _setup_pipeline_logger():
    import logging

    logger = logging.getLogger("apex.pipeline")
    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter("%(message)s"))
        logger.addHandler(handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    return logger


def _cmd_run(args: argparse.Namespace) -> int:
    """Run the shared headless pipeline (Steps 1-7)."""
    try:
        from apex.pipeline import RunContext, PipelineRunner, get_steps, parse_step_range
    except Exception as exc:  # noqa: BLE001
        print(f"{_FAIL} Could not load the pipeline: {exc}")
        return 1

    config_path = args.config or str(_params_path())
    logger = _setup_pipeline_logger()
    try:
        ctx = RunContext.build(
            args.mode,
            config_path,
            result_dir=args.result_dir,
            data_dir=args.data_dir,
            force=args.force,
            dry_run=args.dry_run,
            logger=logger,
        )
    except (FileNotFoundError, ValueError) as exc:
        print(f"{_FAIL} {exc}")
        return 1

    only = parse_step_range(args.steps) if args.steps else None
    runner = PipelineRunner(get_steps(args.mode))
    report = runner.run(ctx, only=only)

    print("\nPipeline summary:")
    for r in report.results:
        print(f"  [{r.status:>16}] step {r.index} {r.key:<11} {r.message}")
    if report.success:
        print(f"\n{_OK} Pipeline run completed.")
        return 0
    print(f"\n{_WARN} Pipeline stopped before completing all requested steps.")
    return 1


# ── export ─────────────────────────────────────────────────────────────────

def _cmd_export(args: argparse.Namespace) -> int:
    """Export a light-curve table to a community photometry format.

    Qt-free: imports pandas and the (GUI-free) ``apex.io.exporters`` lazily so
    the rest of the CLI stays import-light.
    """
    try:
        import pandas as pd
    except Exception as exc:  # noqa: BLE001
        print(f"{_FAIL} pandas is required for 'apex export': {exc}")
        return 1

    in_path = Path(args.input)
    if not in_path.exists():
        print(f"{_FAIL} input file not found: {in_path}")
        return 1

    # Sniff the delimiter so both APEX CSVs and TSVs load.
    sep = "\t" if in_path.suffix.lower() in (".tsv", ".tab") else None
    try:
        lc_table = pd.read_csv(in_path, sep=sep, engine="python", comment="#")
    except Exception as exc:  # noqa: BLE001
        print(f"{_FAIL} could not read {in_path}: {exc}")
        return 1

    meta = {
        "obscode": args.obscode,
        "target": args.target,
        "filter": args.filter,
        "observer": args.observer,
        "telescope": args.telescope,
    }
    meta = {k: v for k, v in meta.items() if v}

    try:
        from apex.io.exporters import (
            export_aavso_extended,
            export_exoclock,
            export_exofop,
        )
    except Exception as exc:  # noqa: BLE001
        print(f"{_FAIL} could not load exporters: {exc}")
        return 1

    exporters = {
        "aavso": export_aavso_extended,
        "exoclock": export_exoclock,
        "exofop": export_exofop,
    }
    try:
        written = exporters[args.format](lc_table, meta, Path(args.output))
    except (ValueError, TypeError, KeyError) as exc:
        print(f"{_FAIL} export failed: {exc}")
        return 1

    print(f"{_OK} Wrote {args.format} file: {written}")
    return 0


def _load_validate_config(config_path: Optional[str]):
    """Load an optional TOML config for the validation harness.

    Returns a lightweight namespace with the fields the suites look for
    (``parameter_file``, ``iraf``, ``known_targets``). Qt-free.
    """
    if not config_path:
        return None
    path = Path(config_path)
    if not path.exists():
        print(f"{_WARN} validate config not found, ignoring: {path}")
        return None
    from apex.utils.io_utils import load_toml
    raw = load_toml(path)  # BOM-tolerant (PowerShell-written files)
    section = raw.get("validate", raw)
    import types

    ns = types.SimpleNamespace()
    ns.parameter_file = section.get("parameter_file", "parameters.toml")
    ns.iraf = section.get("iraf", {})
    ns.known_targets = section.get("known_targets", {})
    ns.crosscheck = section.get("crosscheck", {})
    return ns


def _cmd_validate(args: argparse.Namespace) -> int:
    """Run the reproducible validation harness and write a report."""
    try:
        from apex.benchmark.validate import run_validation
    except Exception as exc:  # noqa: BLE001
        print(f"{_FAIL} Could not load the validation harness: {exc}")
        return 1

    suite = getattr(args, "suite", "all") or "all"
    output_dir = Path(getattr(args, "output", None) or "validation/report")
    config = _load_validate_config(getattr(args, "config", None))

    # CLI flags for the crosscheck suite override (or supply) the TOML config.
    result_dir = getattr(args, "result_dir", None)
    frame = getattr(args, "frame", None)
    reference = getattr(args, "reference", None)
    if result_dir or frame or reference:
        import types

        if config is None:
            config = types.SimpleNamespace(
                parameter_file="parameters.toml", iraf={}, known_targets={}, crosscheck={}
            )
        cc = dict(getattr(config, "crosscheck", {}) or {})
        if result_dir:
            cc["result_dir"] = result_dir
        if frame:
            cc["frame"] = frame
        if reference:
            cc["references"] = list(REFERENCE_CHOICES) if reference == "all" else [reference]
        config.crosscheck = cc

    print(f"APEX validate - suite={suite}, output={output_dir}")
    try:
        results = run_validation(
            suite,
            output_dir,
            config=config,
            reference_frame=getattr(args, "reference_frame", None),
        )
    except ValueError as exc:
        print(f"{_FAIL} {exc}")
        return 1
    except Exception as exc:  # noqa: BLE001
        print(f"{_FAIL} validation failed: {exc}")
        return 1

    suites = results.get("suites", {})
    print()
    n_fail = 0
    for name, result in suites.items():
        status = result.get("status", "unknown")
        if name in ("iraf", "crosscheck") and status == "skipped":
            tag = _WARN
        elif status == "ok" and result.get("status") == "ok":
            tag = _OK
            if name == "artificial-star" and not result.get("completeness_non_degenerate"):
                tag = _WARN
            if name == "known-targets" and not result.get("passed"):
                tag = _FAIL
                n_fail += 1
        else:
            tag = _FAIL
            n_fail += 1
        print(f"  {tag} {name}: {status}")

    paths = results.get("report_paths", {})
    print()
    print(f"{_OK} Report: {paths.get('summary')}")
    print(f"{_OK} Manifest: {paths.get('manifest')}")
    print(f"{_OK} Docs: {paths.get('docs')}")
    return 1 if n_fail else 0


def _cmd_gui(args: argparse.Namespace) -> int:
    """Launch the desktop GUI. PyQt5 is imported lazily here only."""
    try:
        from apex.gui_entry import main as gui_main
    except Exception as exc:  # noqa: BLE001
        print(f"{_FAIL} Could not load the GUI launcher: {exc}")
        return 1
    return gui_main(mode=getattr(args, "mode", None))


# ── parser ───────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="apex",
        description="APEX - Automated Photometry EXtraction (headless CLI).",
    )
    parser.add_argument("--version", action="version", version=f"APEX {__version__}")
    sub = parser.add_subparsers(dest="command")

    sub.add_parser("version", help="Print the APEX version.")

    p_doctor = sub.add_parser("doctor", help="Diagnose the runtime environment.")
    p_doctor.add_argument(
        "--network",
        action="store_true",
        help="Also probe Gaia/SIMBAD reachability (adds a few seconds).",
    )

    p_config = sub.add_parser("config", help="Manage parameters.toml.")
    p_config.add_argument(
        "config_action",
        choices=["init", "path", "show"],
        help="init: copy the example; path: print location; show: print contents.",
    )

    p_run = sub.add_parser("run", help="Run the shared headless pipeline (Steps 1-7).")
    p_run.add_argument("--mode", choices=["cmd", "lc"], required=True,
                       help="Pipeline mode.")
    p_run.add_argument("--config", default=None,
                       help="Path to parameters.toml (default: ./parameters.toml).")
    p_run.add_argument("--steps", default=None,
                       help="Step subset, e.g. '1-7', '4', '2,4,6' (default: all).")
    p_run.add_argument("--result-dir", default=None, help="Override [io].result_dir.")
    p_run.add_argument("--data-dir", default=None, help="Override [io].data_dir.")
    p_run.add_argument("--force", action="store_true",
                       help="Re-run steps even if their outputs already exist.")
    p_run.add_argument("--dry-run", action="store_true",
                       help="Show the resolved plan without executing.")

    p_export = sub.add_parser(
        "export",
        help="Export an APEX light curve to a community format (AAVSO/ExoClock/ExoFOP).",
    )
    p_export.add_argument("--format", choices=["aavso", "exoclock", "exofop"], required=True,
                          help="Target output format.")
    p_export.add_argument("--input", required=True,
                          help="Input light-curve table (APEX .csv or .tsv).")
    p_export.add_argument("--output", required=True, help="Output file path.")
    p_export.add_argument("--obscode", default=None,
                          help="AAVSO/observatory observer code.")
    p_export.add_argument("--target", default=None,
                          help="Target/star designation (required by AAVSO/ExoClock/ExoFOP).")
    p_export.add_argument("--filter", default=None,
                          help="Band code override (e.g. V); falls back to the table's filter column.")
    p_export.add_argument("--observer", default=None, help="Observer name (ExoClock/ExoFOP).")
    p_export.add_argument("--telescope", default=None, help="Telescope name (ExoClock/ExoFOP).")

    p_validate = sub.add_parser(
        "validate",
        help="Run the reproducible validation harness and write a report.",
    )
    p_validate.add_argument(
        "--suite",
        choices=["artificial-star", "iraf", "known-targets", "crosscheck", "all"],
        default="all",
        help="Which validation suite(s) to run (default: all).",
    )
    p_validate.add_argument(
        "--reference-frame",
        default=None,
        help="Optional reference FITS for the artificial-star suite "
        "(default: generate a self-contained synthetic frame).",
    )
    p_validate.add_argument(
        "--result-dir",
        default=None,
        help="APEX result dir for the crosscheck suite (the dir containing "
        "step7_forced_phot/).",
    )
    p_validate.add_argument(
        "--frame",
        default=None,
        help="Specific frame name to cross-check (default: auto-pick a "
        "representative status==ok frame near the median FWHM).",
    )
    p_validate.add_argument(
        "--reference",
        choices=["sep", "iraf", "photutils", "all"],
        default="sep",
        help="Independent photometry engine for the crosscheck suite "
        "(default: sep).",
    )
    p_validate.add_argument(
        "--output",
        default="validation/report",
        help="Output directory for the report (default: validation/report).",
    )
    p_validate.add_argument(
        "--config",
        default=None,
        help="Optional TOML config (parameter_file, [iraf], [known_targets], [crosscheck]).",
    )

    p_gui = sub.add_parser("gui", help="Launch the desktop GUI.")
    p_gui.add_argument("--mode", choices=["cmd", "lc"], default=None,
                       help="Skip the launcher and open a mode directly.")
    return parser


def main(argv: Optional[list] = None) -> int:
    # Windows consoles default to a legacy codepage (e.g. cp949) that cannot
    # encode characters that show up in dependency versions or file paths.
    for stream in (sys.stdout, sys.stderr):
        try:
            stream.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
        except Exception:
            pass
    parser = build_parser()
    args = parser.parse_args(argv)
    command = args.command or "version"
    handlers = {
        "version": _cmd_version,
        "doctor": _cmd_doctor,
        "config": _cmd_config,
        "run": _cmd_run,
        "export": _cmd_export,
        "validate": _cmd_validate,
        "gui": _cmd_gui,
    }
    handler = handlers.get(command)
    if handler is None:
        parser.print_help()
        return 2
    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
