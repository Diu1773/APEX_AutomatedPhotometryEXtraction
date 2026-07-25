#!/usr/bin/env python3
"""Audit every APEX window against a target screen size.

Why this exists
---------------
``apex/gui/layout_rules.py`` documents four root causes of the recurring
clipping complaints. Three of them are invisible until someone opens the
window on a small screen and notices — which is how they kept coming back.
This walks every window on a real project and reports them mechanically.

The checks, and the bug each one was written for:

``window_over_budget``
    The window's ``minimumSizeHint`` is larger than the screen. The OS clamps
    the window down, Qt refuses to lay out below the minimum, and everything
    past the cut (nav row, sliders, Run buttons) is unreachable — there is no
    scrollbar, because the window itself is the thing that overflowed.

``hclip_in_scroll``
    A ``QScrollArea`` with the horizontal scrollbar off, whose content is wider
    than its viewport. The overflow is not merely hidden, it is *unreachable*.
    This is what ``scroll_wrap`` caused when it hid a fixed-width control
    column's width minimum along with its height minimum.

``nested_scroll``
    A scroll area inside another scroll area: the user cannot tell which
    surface they are scrolling.

``clipped_text``
    A button/label/checkbox laid out narrower than its ``sizeHint`` — the label
    is cut, usually by a hardcoded pixel width measured on another font.

``starved_scroll``
    A scroll viewport so small relative to its content that scrolling it is
    useless (a 150 px window onto an 800 px plot).

Usage
-----
    .venv-deploy\\Scripts\\python.exe scripts/layout_audit.py PROJECT_DIR
    .venv-deploy\\Scripts\\python.exe scripts/layout_audit.py PROJECT_DIR --width 1280 --height 704

PROJECT_DIR must contain a ``parameters.toml`` for a project whose pipeline has
been run, so the windows populate with real content. Windows that cannot be
constructed are listed as skipped — never silently counted as passing.
"""
from __future__ import annotations

import argparse
import inspect
import importlib
import sys
import tempfile
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
if str(_REPO) not in sys.path:
    sys.path.insert(0, str(_REPO))

# Default target: the smallest screen APEX is expected to be usable on — a
# 16" 2560x1600 laptop panel at 200% scaling is 1280x800 logical, minus the
# margin layout_rules keeps for the title bar.
DEFAULT_W, DEFAULT_H = 1280, 704

STEP_NAMES = [
    "File Selection", "Image Crop", "Sky Preview & QC", "Source Detection",
    "WCS Plate Solving", "Master Catalog Build", "Forced Aperture Phot",
    "PSF Photometry", "Step 9", "Step 10", "Step 11", "Step 12",
]

# (module, class or None to auto-detect, mode, step index)
TARGETS = [
    ("apex.gui.workflow.step1_file_selection_common", "CommonFileSelectionWindow", "cmd", 0),
    ("apex.gui.workflow.step2_crop_selector", "CropSelectorWindow", "cmd", 1),
    ("apex.gui.workflow.step3_sky_preview", "SkyPreviewWindow", "cmd", 2),
    ("apex.gui.workflow.step4_source_detection", "SourceDetectionWindow", "cmd", 3),
    ("apex.gui.workflow.step6_ref_build", "RefBuildWindow", "cmd", 5),
    ("apex.gui.workflow.step7_forced_aperture_phot", "ForcedPhotWindow", "cmd", 6),
    ("apex.gui.workflow.cmd.step8_psf_photometry", "PSFPhotometryWindow", "cmd", 7),
    ("apex.gui.workflow.cmd.step9_master_id_editor", "MasterIdEditorWindow", "cmd", 8),
    ("apex.gui.workflow.cmd.step10_zeropoint_calibration", "ZeropointCalibrationWindow", "cmd", 9),
    ("apex.gui.workflow.cmd.step11_cmd_plot", "CmdPlotWindow", "cmd", 10),
    ("apex.gui.workflow.cmd.step12_isochrone_model", "IsochroneModelWindow", "cmd", 11),
    ("apex.gui.workflow.lc.step8_target_selection", "TargetComparisonSelectionWindow", "lc", 8),
    ("apex.gui.workflow.lc.step9_lightcurve_builder", "LightCurveBuilderWindow", "lc", 9),
    ("apex.gui.workflow.lc.step10_detrend_merge", "DetrendNightMergeWindow", "lc", 10),
    ("apex.gui.workflow.lc.step11_period_analysis", "PeriodAnalysisWindow", "lc", 11),
]


def _tool_targets():
    from apex.gui.tools.registry import TOOL_SPECS
    return [(spec.module, None, spec.modes[0], 11) for spec in TOOL_SPECS]


def _make_env(project_toml: Path, mode: str):
    if mode == "cmd":
        from apex.config.parameters_cmd import read_params
    else:
        from apex.config.parameters_lc import read_params
    from apex.core import FileManager, ProjectState
    from PyQt5.QtWidgets import QMainWindow

    params = read_params(project_toml)
    fm = FileManager(params)
    ps = ProjectState(Path(tempfile.mkdtemp()))
    ps.assign_steps(list(STEP_NAMES))

    class StubMain(QMainWindow):
        """Enough of MainWindowWorkflow for a step window to build itself."""

        def __init__(self):
            super().__init__()
            self.params, self.project_state, self.file_manager = params, ps, fm
            self.mode, self.step_names = mode, list(STEP_NAMES)

        def __getattr__(self, name):
            if name.startswith(("on_", "open_", "update_", "append_", "show_")):
                return lambda *a, **k: None
            raise AttributeError(name)

    return params, fm, ps, StubMain()


def _auto_class(module):
    from PyQt5.QtWidgets import QWidget
    best = None
    for name, obj in vars(module).items():
        if (inspect.isclass(obj) and obj.__module__ == module.__name__
                and issubclass(obj, QWidget)
                and any(k in name for k in ("Window", "Dialog", "Tool"))):
            if best is None or len(name) > len(best.__name__):
                best = obj
    return best


def _build(cls, env, mode, idx, params):
    _, fm, ps, main = env
    by_name = {
        "params": params, "file_manager": fm, "project_state": ps,
        "main_window": main, "parent": None, "step_index": idx,
        "step_name": STEP_NAMES[min(idx, len(STEP_NAMES) - 1)], "mode": mode,
        "data_dir": Path(params.P.data_dir), "result_dir": Path(params.P.result_dir),
    }
    args, missing = {}, []
    for name, p in list(inspect.signature(cls.__init__).parameters.items())[1:]:
        if name in by_name:
            args[name] = by_name[name]
        elif p.default is not inspect.Parameter.empty or p.kind in (
                p.VAR_POSITIONAL, p.VAR_KEYWORD):
            continue
        else:
            missing.append(name)
    if missing:
        raise TypeError("unknown required args: " + ", ".join(missing))
    return cls(**args)


def audit_window(win, budget_w: int, budget_h: int) -> list[tuple[str, str, str]]:
    """Return [(check, where, detail)] for one already-shown window."""
    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import (
        QAbstractScrollArea, QCheckBox, QLabel, QPushButton, QRadioButton,
        QScrollArea, QWidget,
    )

    found: list[tuple[str, str, str]] = []

    mh = win.minimumSizeHint()
    if mh.height() > budget_h or mh.width() > budget_w:
        found.append((
            "window_over_budget", win.__class__.__name__,
            f"minimum {mh.width()}x{mh.height()} vs screen {budget_w}x{budget_h} "
            f"(over by {max(0, mh.width()-budget_w)}x{max(0, mh.height()-budget_h)})",
        ))

    # Qt-internal plumbing is legitimately zero-width when a header is hidden.
    INTERNAL = {"QHeaderView", "QAbstractButton", "QWidget"}

    for w in win.findChildren(QWidget):
        if not w.isVisible():
            continue
        cls_name = w.__class__.__name__
        gw, gh = w.width(), w.height()
        internal = cls_name in INTERNAL or str(w.objectName()).startswith("qt_")

        if (gw <= 0 or gh <= 0) and not internal:
            found.append(("zero_size", cls_name, f"{gw}x{gh} px"))
            continue

        if isinstance(w, (QPushButton, QLabel, QCheckBox, QRadioButton)):
            text = (w.text() or "").strip()
            wraps = getattr(w, "wordWrap", None)
            if text and not (callable(wraps) and wraps()):
                need = w.sizeHint().width()
                if need > gw + 1:
                    found.append((
                        "clipped_text", f"{cls_name} {text[:32]!r}",
                        f"needs {need} px, has {gw}",
                    ))

        if isinstance(w, QScrollArea):
            inner = w.widget()
            if inner is None:
                continue
            vp = w.viewport()

            # unreachable horizontal overflow
            if (w.horizontalScrollBarPolicy() == Qt.ScrollBarAlwaysOff
                    and inner.width() > vp.width() + 1):
                found.append((
                    "hclip_in_scroll", cls_name,
                    f"content {inner.width()} px in a {vp.width()} px viewport with "
                    "the horizontal scrollbar off — the overflow cannot be reached",
                ))

            if inner.findChildren(QScrollArea):
                found.append((
                    "nested_scroll", cls_name,
                    f"{len(inner.findChildren(QScrollArea))} scroll area(s) inside "
                    "another — ambiguous which one scrolls",
                ))

            need_h = inner.minimumSizeHint().height()
            if need_h > 0 and 0 < vp.height() < min(160, need_h * 0.45):
                found.append((
                    "starved_scroll", cls_name,
                    f"viewport {vp.height()} px onto {need_h} px of content "
                    f"({100 * vp.height() / need_h:.0f}%)",
                ))

    return list(dict.fromkeys(found))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("project", type=Path, help="project dir holding parameters.toml")
    ap.add_argument("--width", type=int, default=DEFAULT_W)
    ap.add_argument("--height", type=int, default=DEFAULT_H)
    ap.add_argument("--shots", type=Path, default=None,
                    help="directory to save a capture of every window with findings")
    args = ap.parse_args()

    toml = args.project if args.project.suffix == ".toml" else args.project / "parameters.toml"
    if not toml.exists():
        print(f"no parameters.toml at {toml}", file=sys.stderr)
        return 2

    from PyQt5.QtCore import Qt
    from PyQt5.QtWidgets import QApplication, QTabWidget

    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    app = QApplication(sys.argv[:1])
    from apex.gui.theme import apply_theme
    apply_theme(app)

    envs = {m: _make_env(toml, m) for m in ("cmd", "lc")}
    if args.shots:
        args.shots.mkdir(parents=True, exist_ok=True)

    results, skipped = [], []
    for module, clsname, mode, idx in TARGETS + _tool_targets():
        label = clsname or module.rsplit(".", 1)[-1]
        try:
            mod = importlib.import_module(module)
            cls = getattr(mod, clsname) if clsname else _auto_class(mod)
            if cls is None:
                skipped.append((label, "no window class found"))
                continue
            params = envs[mode][0]
            win = _build(cls, envs[mode], mode, idx, params)
            win.show()
            for _ in range(8):
                app.processEvents()
            win.setMinimumSize(0, 0)
            win.setGeometry(30, 20, args.width, args.height)
            for _ in range(10):
                app.processEvents()

            found = audit_window(win, args.width, args.height)
            for tabs in win.findChildren(QTabWidget):
                for i in range(tabs.count()):
                    tabs.setCurrentIndex(i)
                    for _ in range(6):
                        app.processEvents()
                    found += audit_window(win, args.width, args.height)
            found = list(dict.fromkeys(found))
            results.append((cls.__name__, found))
            if found and args.shots:
                win.grab().save(str(args.shots / f"{cls.__name__}.png"))
            win.close()
            for _ in range(3):
                app.processEvents()
        except Exception as exc:
            skipped.append((label, f"{type(exc).__name__}: {exc}"))

    bad = [r for r in results if r[1]]
    print("=" * 78)
    print(f"layout audit @ {args.width}x{args.height}")
    print("=" * 78)
    print(f"{len(results)} windows checked - {len(bad)} with findings - "
          f"{len(skipped)} skipped\n")

    by_check: dict[str, int] = {}
    for _, found in results:
        for check, _, _ in found:
            by_check[check] = by_check.get(check, 0) + 1
    if by_check:
        print("findings by check:")
        for check, n in sorted(by_check.items(), key=lambda kv: -kv[1]):
            print(f"   {check:22s} {n}")
        print()

    for name, found in sorted(bad, key=lambda r: -len(r[1])):
        print(f"* {name}  ({len(found)})")
        groups: dict[str, list] = {}
        for check, where, detail in found:
            groups.setdefault(check, []).append((where, detail))
        for check, items in groups.items():
            print(f"    {check} ({len(items)})")
            for where, detail in items[:5]:
                print(f"       - {where}: {detail}")
            if len(items) > 5:
                print(f"       ... {len(items)-5} more")
        print()

    clean = [n for n, f in results if not f]
    if clean:
        print("clean: " + ", ".join(clean) + "\n")
    if skipped:
        print("skipped (could not construct - check manually):")
        for label, why in skipped:
            print(f"   - {label}: {why[:100]}")

    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
