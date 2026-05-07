"""Central registry for Tools-menu entries.

Tool implementations can stay in their current modules while menu ownership
is centralized here. Add new GUI tools by adding one ToolSpec instead of editing
the menu-building logic in main_window.py.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable


Mode = str


@dataclass(frozen=True)
class ToolSpec:
    id: str
    label: str
    launcher: str
    modes: tuple[Mode, ...] = ("cmd", "lc")
    shortcut: str | None = None
    separator_before: bool = False


TOOL_SPECS: tuple[ToolSpec, ...] = (
    ToolSpec(
        id="qa_report",
        label="QA Reports",
        launcher="open_qa_report",
        shortcut="Ctrl+R",
    ),
    ToolSpec(
        id="iraf_photometry",
        label="IRAF/DAOPHOT Tool",
        launcher="open_iraf_tool",
        shortcut="Ctrl+I",
    ),
    ToolSpec(
        id="extinction_fit",
        label="Extinction Fit Tool",
        launcher="open_extinction_tool",
        modes=("lc",),
        shortcut="Ctrl+E",
    ),
    ToolSpec(
        id="airmass_debug",
        label="Airmass Header Debug",
        launcher="open_airmass_debug_tool",
        modes=("lc",),
    ),
    ToolSpec(
        id="multi_night_merger",
        label="Multi-Night Light Curve Merger",
        launcher="open_multi_night_merger",
        modes=("lc",),
        shortcut="Ctrl+M",
        separator_before=True,
    ),
    ToolSpec(
        id="variable_star",
        label="Variable Star Analysis",
        launcher="open_variable_star_tool",
        modes=("lc",),
        shortcut="Ctrl+Shift+V",
        separator_before=True,
    ),
    ToolSpec(
        id="transit",
        label="Exoplanet Transit Analysis",
        launcher="open_transit_tool",
        modes=("lc",),
        shortcut="Ctrl+Shift+T",
    ),
    ToolSpec(
        id="eclipsing_binary",
        label="Eclipsing Binary Analysis",
        launcher="open_eb_tool",
        modes=("lc",),
        shortcut="Ctrl+Shift+B",
    ),
    ToolSpec(
        id="extinction_fit",
        label="Extinction (Airmass Fit)",
        launcher="open_extinction_fit_cmd",
        modes=("cmd",),
    ),
    ToolSpec(
        id="cmd_isochrone",
        label="CMD + Isochrone (From Results)...",
        launcher="open_cmd_iso_tool",
        modes=("cmd",),
    ),
    ToolSpec(
        id="gaia_3d_viewer",
        label="Gaia 3D Cluster Viewer",
        launcher="open_gaia_3d_viewer",
        modes=("cmd",),
    ),
    ToolSpec(
        id="cluster_structure",
        label="Analyze Cluster Structure",
        launcher="open_cluster_structure_tool",
        modes=("cmd",),
    ),
)


def iter_tools_for_mode(mode: Mode) -> Iterable[ToolSpec]:
    """Yield Tools-menu specs visible in the requested application mode."""
    mode_key = str(mode or "").lower()
    for spec in TOOL_SPECS:
        if mode_key in spec.modes:
            yield spec

