from __future__ import annotations

import ast
from pathlib import Path

from apex.gui.tools.registry import TOOL_SPECS, iter_tool_modules, iter_tools_for_mode


ROOT = Path(__file__).resolve().parents[1]


def _main_window_methods() -> set[str]:
    tree = ast.parse((ROOT / "apex/gui/main_window.py").read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.ClassDef) and node.name == "MainWindowWorkflow":
            return {item.name for item in node.body if isinstance(item, ast.FunctionDef)}
    raise AssertionError("MainWindowWorkflow class not found")


def test_tool_specs_have_unique_ids_and_shortcuts_per_mode():
    ids = [spec.id for spec in TOOL_SPECS]
    assert len(ids) == len(set(ids))

    for mode in ("cmd", "lc"):
        shortcuts = [spec.shortcut for spec in iter_tools_for_mode(mode) if spec.shortcut]
        assert len(shortcuts) == len(set(shortcuts))


def test_tool_spec_launchers_exist_on_main_window():
    methods = _main_window_methods()
    missing = sorted({spec.launcher for spec in TOOL_SPECS} - methods)
    assert missing == []


def test_tool_specs_have_release_smoke_modules():
    modules = list(iter_tool_modules())
    assert len(modules) == len(set(modules))
    assert len(modules) == len(TOOL_SPECS)
    assert all(module.startswith("apex.gui.tools.") for module in modules)


def test_common_and_mode_specific_tools_are_exposed_once():
    cmd_ids = [spec.id for spec in iter_tools_for_mode("cmd")]
    lc_ids = [spec.id for spec in iter_tools_for_mode("lc")]

    for common_id in ("qa_report", "iraf_photometry", "extinction_fit"):
        assert cmd_ids.count(common_id) == 1
        assert lc_ids.count(common_id) == 1

    assert "cmd_isochrone" not in cmd_ids
    assert "cmd_isochrone" not in lc_ids
    assert "variable_star" in lc_ids
    assert "variable_star" not in cmd_ids
