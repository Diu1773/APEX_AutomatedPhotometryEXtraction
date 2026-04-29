#!/usr/bin/env python3
"""
APEX flat restructure migration script.

Copies files from 3-package layout (common, cmd, lightcurve) into flat
single-software layout under apex/, rewrites imports, creates new files,
then syntax-checks all new Python files.

Does NOT delete the old structure until after validation.
"""

import shutil
import sys
import ast
import re
from pathlib import Path

ROOT = Path("/mnt/c/Users/bmffr/Desktop/Result/Automated_Photometry_EXtraction")
APEX = ROOT / "apex"

errors: list[str] = []
created: list[str] = []
skipped: list[str] = []


# ── helpers ──────────────────────────────────────────────────────────────────

def ensure(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_init(path: Path, content: str = "") -> None:
    path.write_text(content, encoding="utf-8")
    created.append(str(path))


def copy_and_rewrite(src: Path, dst: Path, rewrites: list[tuple[str, str]]) -> bool:
    """Copy src → dst applying string rewrites in order (longest-key-first already sorted by caller)."""
    try:
        text = src.read_text(encoding="utf-8")
        for old, new in rewrites:
            text = text.replace(old, new)
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(text, encoding="utf-8")
        created.append(str(dst))
        return True
    except Exception as e:
        errors.append(f"COPY FAILED {src} → {dst}: {e}")
        return False


# ── Import rewrite rules (longer patterns FIRST to avoid partial replacements) ──
# Applied in this order to every file we copy.
# Note: relative imports inside the moved files are also handled here via
# pattern-specific rewrites injected per file group.

COMMON_REWRITES: list[tuple[str, str]] = [
    # Rule 8: from apex.lightcurve.config → from apex.config
    ("from apex.lightcurve.config", "from apex.config"),
    # Rule 9: from apex.lightcurve.utils → from apex.utils
    ("from apex.lightcurve.utils", "from apex.utils"),
    # Rule 10: from apex.lightcurve.analysis → from apex.analysis
    ("from apex.lightcurve.analysis", "from apex.analysis"),
    # Rule 6: from apex.cmd.config → from apex.config
    ("from apex.cmd.config", "from apex.config"),
    # Rule 7: from apex.cmd.utils → from apex.utils
    ("from apex.cmd.utils", "from apex.utils"),
    # Rule 3: from apex.common.gui.workflow.shared → from apex.gui.workflow
    ("from apex.common.gui.workflow.shared", "from apex.gui.workflow"),
    # Rule 4: from apex.common.gui.workflow.step_window_base → from apex.gui.workflow.step_window_base
    ("from apex.common.gui.workflow.step_window_base", "from apex.gui.workflow.step_window_base"),
    # Rule 5: from apex.common.gui.widgets → from apex.gui.widgets
    ("from apex.common.gui.widgets", "from apex.gui.widgets"),
    # Rule 1: from apex.common.core → from apex.core
    ("from apex.common.core", "from apex.core"),
    # Rule 2: from apex.common.utils → from apex.utils
    ("from apex.common.utils", "from apex.utils"),
    # Rule 12: AperturePhotometryWindow import path
    (
        "from apex.cmd.gui.workflow.step5_aperture_photometry import AperturePhotometryWindow",
        "from apex.gui.workflow.step5_aperture_photometry import AperturePhotometryWindow",
    ),
    # Rule 11: from apex.cmd.gui.workflow → from apex.gui.workflow.cmd
    ("from apex.cmd.gui.workflow", "from apex.gui.workflow.cmd"),
    # step_paths_base rename
    ("from apex.utils.step_paths_base", "from apex.utils.step_paths"),
    ("apex.utils.step_paths_base", "apex.utils.step_paths"),
    ("step_paths_base", "step_paths"),  # bare module ref in *-imports
]

# Additional per-group rewrites for relative imports
# (applied after COMMON_REWRITES, before writing)

def apply_rewrites(text: str, extra: list[tuple[str, str]]) -> str:
    for old, new in COMMON_REWRITES + extra:
        text = text.replace(old, new)
    return text


def copy_rewrite(src: Path, dst: Path, extra: list[tuple[str, str]] = None) -> bool:
    try:
        text = src.read_text(encoding="utf-8")
        text = apply_rewrites(text, extra or [])
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(text, encoding="utf-8")
        created.append(str(dst))
        return True
    except Exception as e:
        errors.append(f"COPY FAILED {src} → {dst}: {e}")
        return False


# ── relative-import rewrite helpers per package location ────────────────────

# For files that were in apex/common/core/ → apex/core/
# Relative imports: ..utils.X → apex.utils.X,  ..utils → apex.utils
CORE_EXTRA = [
    ("from ..utils.step_paths_base", "from apex.utils.step_paths"),
    ("from ..utils.step_paths", "from apex.utils.step_paths"),
    ("from ..utils import", "from apex.utils import"),
    ("from ..utils.", "from apex.utils."),
]

# For files in apex/common/utils/ → apex/utils/
# These have single-dot relative imports: from .something
# No changes needed for within-utils relative imports — they stay as .X
# But step_paths_base rename:
UTILS_EXTRA = [
    # Within apex/utils/, .step_paths_base → .step_paths (for header_cache, qc_utils)
    ("from .step_paths_base import", "from .step_paths import"),
    ("from .step_paths_base ", "from .step_paths "),
]

# For shared workflow files (apex/common/gui/workflow/shared/) → apex/gui/workflow/
# They have 4-dot relative imports: from ....utils.X → from apex.utils.X
# and from ..step_window_base → from .step_window_base (stays same level now)
SHARED_WORKFLOW_EXTRA = [
    # 4-dot up to common → now absolute
    ("from ....utils.step_paths_base", "from apex.utils.step_paths"),
    ("from ....utils.step_paths", "from apex.utils.step_paths"),
    ("from ....utils.constants", "from apex.utils.constants"),
    ("from ....utils.photometry_utils", "from apex.utils.photometry_utils"),
    ("from ....utils.qc_utils", "from apex.utils.qc_utils"),
    ("from ....utils.cache_utils", "from apex.utils.cache_utils"),
    ("from ....utils.io_utils", "from apex.utils.io_utils"),
    ("from ....utils.common_helpers", "from apex.utils.common_helpers"),
    ("from ....utils.astro_utils", "from apex.utils.astro_utils"),
    ("from ....utils import", "from apex.utils import"),
    ("from ....utils.", "from apex.utils."),
    # ..step_window_base stays as .step_window_base
    ("from ..step_window_base", "from .step_window_base"),
]

# For files in apex/cmd/gui/workflow/ → apex/gui/workflow/cmd/
# 3-dot = ..utils → apex.utils, ..config → apex.config
# 4-dot = ....common.utils → apex.utils
CMD_WORKFLOW_EXTRA = [
    ("from ....common.utils.step_paths_base", "from apex.utils.step_paths"),
    ("from ....common.utils.astro_utils", "from apex.utils.astro_utils"),
    ("from ....common.utils.photometry_utils", "from apex.utils.photometry_utils"),
    ("from ....common.utils.constants", "from apex.utils.constants"),
    ("from ....common.utils.io_utils", "from apex.utils.io_utils"),
    ("from ....common.utils.common_helpers", "from apex.utils.common_helpers"),
    ("from ....common.utils.qc_utils", "from apex.utils.qc_utils"),
    ("from ....common.utils.cache_utils", "from apex.utils.cache_utils"),
    ("from ....common.utils", "from apex.utils"),
    ("from ....common.core", "from apex.core"),
    ("from ...utils.step_paths", "from apex.utils.step_paths_cmd"),
    ("from ...utils import normalize_filter_name", "from apex.utils.astro_utils import normalize_filter_name"),
    ("from ...utils import", "from apex.utils import"),
    ("from ...utils.", "from apex.utils."),
    ("from ...config", "from apex.config"),
    ("from ...analysis.isochrone_fitter_v2", "from apex.analysis.cmd.isochrone_fitter_v2"),
    ("from ...analysis.isochrone_fitter", "from apex.analysis.cmd.isochrone_fitter"),
    # step_window_base: ..step_window_base or .step_window_base
    ("from ..step_window_base", "from apex.gui.workflow.step_window_base"),
    # Intra-workflow: from .workflow.step → from apex.gui.workflow.cmd.step
    ("from ..workflow.step11_zeropoint_calibration", "from apex.gui.workflow.cmd.step11_zeropoint_calibration"),
    ("from ..workflow.step13_isochrone_model", "from apex.gui.workflow.cmd.step13_isochrone_model"),
]

# For files in apex/cmd/gui/tools/ → apex/gui/tools/
CMD_TOOLS_EXTRA = [
    ("from ....common.utils.step_paths_base", "from apex.utils.step_paths"),
    ("from ....common.utils.astro_utils", "from apex.utils.astro_utils"),
    ("from ....common.utils.photometry_utils", "from apex.utils.photometry_utils"),
    ("from ....common.utils.constants", "from apex.utils.constants"),
    ("from ....common.utils.io_utils", "from apex.utils.io_utils"),
    ("from ....common.utils.common_helpers", "from apex.utils.common_helpers"),
    ("from ....common.utils.qc_utils", "from apex.utils.qc_utils"),
    ("from ....common.utils.cache_utils", "from apex.utils.cache_utils"),
    ("from ....common.utils", "from apex.utils"),
    ("from ....common.core", "from apex.core"),
    ("from ...utils.step_paths", "from apex.utils.step_paths_cmd"),
    ("from ...utils import", "from apex.utils import"),
    ("from ...utils.", "from apex.utils."),
    ("from ...config", "from apex.config"),
    # Intra-tool: relative to workflow
    ("from ..workflow.step11_zeropoint_calibration", "from apex.gui.workflow.cmd.step11_zeropoint_calibration"),
    ("from ..workflow.step13_isochrone_model", "from apex.gui.workflow.cmd.step13_isochrone_model"),
    # step_window_base
    ("from ..workflow.step_window_base", "from apex.gui.workflow.step_window_base"),
    # cluster_structure __init__
    ("from .window import", "from .window import"),  # no change needed
]

# For CMD cluster_structure sub-package (apex/cmd/gui/tools/cluster_structure/) → apex/gui/tools/cluster_structure/
CMD_CLUSTER_EXTRA = [
    ("from ....utils.step_paths", "from apex.utils.step_paths_cmd"),
    ("from ....utils.", "from apex.utils."),
    ("from ....utils import", "from apex.utils import"),
    ("from .....common.utils", "from apex.utils"),
]

# For LC gui/workflow files → apex/gui/workflow/lc/ (lc-specific) or apex/gui/workflow/ (shared)
LC_WORKFLOW_EXTRA = [
    ("from ....common.utils.step_paths_base", "from apex.utils.step_paths"),
    ("from ....common.utils.astro_utils", "from apex.utils.astro_utils"),
    ("from ....common.utils.photometry_utils", "from apex.utils.photometry_utils"),
    ("from ....common.utils.constants", "from apex.utils.constants"),
    ("from ....common.utils.io_utils", "from apex.utils.io_utils"),
    ("from ....common.utils.common_helpers", "from apex.utils.common_helpers"),
    ("from ....common.utils.qc_utils", "from apex.utils.qc_utils"),
    ("from ....common.utils.cache_utils", "from apex.utils.cache_utils"),
    ("from ....common.utils", "from apex.utils"),
    ("from ....common.core", "from apex.core"),
    ("from ...utils.step_paths", "from apex.utils.step_paths_lc"),
    ("from ...utils.photometry_loader", "from apex.utils.photometry_loader"),
    ("from ...utils.run_workspace", "from apex.utils.run_workspace"),
    ("from ...utils import", "from apex.utils import"),
    ("from ...utils.", "from apex.utils."),
    ("from ...config", "from apex.config"),
    ("from ...analysis.", "from apex.analysis."),
    ("from ...core.", "from apex.core."),
    # step_window_base (shared)
    ("from ..step_window_base", "from apex.gui.workflow.step_window_base"),
    # Tools → moved to apex.gui.tools
    ("from .workflow.qa_report_window", "from apex.gui.tools.qa_report"),
    ("from .workflow.extinction_fit_window", "from apex.gui.tools.extinction_fit"),
    ("from .workflow.iraf_photometry_window", "from apex.gui.tools.iraf_photometry"),
    ("from .workflow.iraf_comparison_window", "from apex.gui.tools.iraf_comparison"),
    ("from .workflow.airmass_header_debug_tool", "from apex.gui.tools.airmass_debug"),
    # LC workflow cross-references
    ("from .workflow.step12_period_analysis import PeriodAnalysisWorker",
     "from apex.gui.workflow.lc.step12_period_analysis import PeriodAnalysisWorker"),
]

# For LC gui/tools files → apex/gui/tools/
LC_TOOLS_EXTRA = [
    ("from ....common.utils.step_paths_base", "from apex.utils.step_paths"),
    ("from ....common.utils.astro_utils", "from apex.utils.astro_utils"),
    ("from ....common.utils.photometry_utils", "from apex.utils.photometry_utils"),
    ("from ....common.utils.constants", "from apex.utils.constants"),
    ("from ....common.utils.io_utils", "from apex.utils.io_utils"),
    ("from ....common.utils.common_helpers", "from apex.utils.common_helpers"),
    ("from ....common.utils.qc_utils", "from apex.utils.qc_utils"),
    ("from ....common.utils.cache_utils", "from apex.utils.cache_utils"),
    ("from ....common.utils", "from apex.utils"),
    ("from ....common.core", "from apex.core"),
    ("from ...utils.step_paths", "from apex.utils.step_paths_lc"),
    ("from ...utils.photometry_loader", "from apex.utils.photometry_loader"),
    ("from ...utils.run_workspace", "from apex.utils.run_workspace"),
    ("from ...utils import", "from apex.utils import"),
    ("from ...utils.", "from apex.utils."),
    ("from ...config", "from apex.config"),
    ("from ...analysis.", "from apex.analysis."),
    ("from ...core.", "from apex.core."),
    # Intra: workflow references
    ("from ..workflow.step12_period_analysis import PeriodAnalysisWorker",
     "from apex.gui.workflow.lc.step12_period_analysis import PeriodAnalysisWorker"),
    ("from ...analysis.light_curve.period_analysis_service", "from apex.analysis.light_curve.period_analysis_service"),
    ("from ...analysis.merge.", "from apex.analysis.merge."),
    ("from ...analysis.merge", "from apex.analysis.merge"),
]

# For LC analysis files → apex/analysis/
LC_ANALYSIS_EXTRA = [
    ("from ...utils.step_paths", "from apex.utils.step_paths_lc"),
    ("from ...utils.photometry_loader", "from apex.utils.photometry_loader"),
    ("from ...utils.run_workspace", "from apex.utils.run_workspace"),
    ("from ...utils.io_utils", "from apex.utils.io_utils"),
    ("from ...utils.", "from apex.utils."),
    ("from ...utils import", "from apex.utils import"),
    ("from ...config", "from apex.config"),
    # For period_io_service: from ...utils.step_paths
    ("from apex.utils.step_paths_lc import step12_period_dir", "from apex.utils.step_paths_lc import step12_period_dir"),
]


# ═══════════════════════════════════════════════════════════════════════════════
# 1. CREATE DIRECTORY STRUCTURE
# ═══════════════════════════════════════════════════════════════════════════════

print("Creating directory structure...")

dirs = [
    APEX / "core",
    APEX / "config",
    APEX / "utils",
    APEX / "analysis" / "light_curve",
    APEX / "analysis" / "merge",
    APEX / "analysis" / "cmd",
    APEX / "gui",
    APEX / "gui" / "widgets",
    APEX / "gui" / "workflow",
    APEX / "gui" / "workflow" / "cmd",
    APEX / "gui" / "workflow" / "lc",
    APEX / "gui" / "tools",
    APEX / "gui" / "tools" / "cluster_structure",
]

for d in dirs:
    ensure(d)

print(f"Created {len(dirs)} directories")


# ═══════════════════════════════════════════════════════════════════════════════
# 2. COPY CORE FILES (from apex/common/core/)
# ═══════════════════════════════════════════════════════════════════════════════

print("\nCopying core files...")

common_core = APEX / "common" / "core"
new_core = APEX / "core"

for fname in ["project_state.py", "file_manager.py", "instrument.py"]:
    copy_rewrite(common_core / fname, new_core / fname, CORE_EXTRA)

# core/__init__.py
write_init(new_core / "__init__.py",
    "from .instrument import InstrumentConfig  # noqa: F401\n"
    "from .file_manager import FileManager  # noqa: F401\n"
    "from .project_state import ProjectState  # noqa: F401\n"
)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. COPY UTILS FILES
# ═══════════════════════════════════════════════════════════════════════════════

print("Copying utils files...")

common_utils = APEX / "common" / "utils"
lc_utils = APEX / "lightcurve" / "utils"
new_utils = APEX / "utils"

# From common/utils/ (canonical source)
for fname in [
    "logging_utils.py", "constants.py", "astro_utils.py", "io_utils.py",
    "common_helpers.py", "cache_utils.py", "photometry_utils.py", "qc_utils.py",
    "param_loader.py",
]:
    copy_rewrite(common_utils / fname, new_utils / fname, UTILS_EXTRA)

# header_cache — only in common/utils
copy_rewrite(common_utils / "header_cache.py", new_utils / "header_cache.py", UTILS_EXTRA)

# step_paths.py: rename of step_paths_base.py (content identical, just renamed)
copy_rewrite(common_utils / "step_paths_base.py", new_utils / "step_paths.py",
    [("step_paths_base", "step_paths")])

# step_paths_cmd.py: from cmd/utils/step_paths.py, update base import
copy_rewrite(
    APEX / "cmd" / "utils" / "step_paths.py",
    new_utils / "step_paths_cmd.py",
    [
        ("from apex.common.utils.step_paths_base import *", "from apex.utils.step_paths import *"),
        ("from apex.common.utils.step_paths_base import step_dir, _as_path",
         "from apex.utils.step_paths import step_dir, _as_path"),
        ("from apex.common.utils.step_paths_base", "from apex.utils.step_paths"),
        ("apex.common.utils.step_paths_base", "apex.utils.step_paths"),
    ]
)

# step_paths_lc.py: from lightcurve/utils/step_paths.py, update base import
copy_rewrite(
    lc_utils / "step_paths.py",
    new_utils / "step_paths_lc.py",
    [
        ("from apex.common.utils.step_paths_base import *", "from apex.utils.step_paths import *"),
        ("from apex.common.utils.step_paths_base import step_dir, _as_path",
         "from apex.utils.step_paths import step_dir, _as_path"),
        ("from apex.common.utils.step_paths_base", "from apex.utils.step_paths"),
        ("apex.common.utils.step_paths_base", "apex.utils.step_paths"),
    ]
)

# LC-specific utils: photometry_loader, run_workspace
for fname in ["photometry_loader.py", "run_workspace.py"]:
    copy_rewrite(lc_utils / fname, new_utils / fname, LC_ANALYSIS_EXTRA)

# utils/__init__.py
write_init(new_utils / "__init__.py",
    "from .logging_utils import get_logger, setup_logger  # noqa: F401\n"
)


# ═══════════════════════════════════════════════════════════════════════════════
# 4. COPY CONFIG FILES
# ═══════════════════════════════════════════════════════════════════════════════

print("Copying config files...")

new_config = APEX / "config"

# parameters_cmd.py: from cmd/config/parameters.py
copy_rewrite(
    APEX / "cmd" / "config" / "parameters.py",
    new_config / "parameters_cmd.py",
    []
)

# parameters_lc.py: from lightcurve/config/parameters.py
copy_rewrite(
    APEX / "lightcurve" / "config" / "parameters.py",
    new_config / "parameters_lc.py",
    []
)

# schema.py: from lightcurve/config/schema.py (more complete)
copy_rewrite(
    APEX / "lightcurve" / "config" / "schema.py",
    new_config / "schema.py",
    [
        ("from apex.lightcurve.config", "from apex.config"),
        ("from apex.lightcurve.utils", "from apex.utils"),
    ]
)

# config/__init__.py
write_init(new_config / "__init__.py",
    "from .parameters_cmd import Parameters as CmdParameters  # noqa: F401\n"
    "from .parameters_lc import Parameters as LcParameters  # noqa: F401\n"
)


# ═══════════════════════════════════════════════════════════════════════════════
# 5. COPY ANALYSIS FILES (from lightcurve/analysis/)
# ═══════════════════════════════════════════════════════════════════════════════

print("Copying analysis files...")

lc_analysis = APEX / "lightcurve" / "analysis"
new_analysis = APEX / "analysis"

# light_curve/
for fname in [
    "asteroid.py", "detrend_output_service.py", "eclipse.py", "global_ensemble.py",
    "lightcurve_output_service.py", "loader.py", "period_analysis_service.py",
    "period_io_service.py", "variable_star.py",
]:
    src = lc_analysis / "light_curve" / fname
    if src.exists():
        copy_rewrite(src, new_analysis / "light_curve" / fname, LC_ANALYSIS_EXTRA)

write_init(new_analysis / "light_curve" / "__init__.py")

# merge/
for fname in ["id_match.py", "workspace_build.py", "workspace_scan.py"]:
    src = lc_analysis / "merge" / fname
    if src.exists():
        copy_rewrite(src, new_analysis / "merge" / fname, LC_ANALYSIS_EXTRA)

write_init(new_analysis / "merge" / "__init__.py")

# analysis/cmd/ — CMD isochrone fitters
cmd_analysis = APEX / "cmd" / "analysis"
for fname in ["isochrone_fitter.py", "isochrone_fitter_v2.py"]:
    src = cmd_analysis / fname
    if src.exists():
        copy_rewrite(src, new_analysis / "cmd" / fname,
            [
                ("from apex.cmd.config", "from apex.config"),
                ("from apex.cmd.utils", "from apex.utils"),
                ("from ..config", "from apex.config"),
                ("from ..utils", "from apex.utils"),
            ]
        )

write_init(new_analysis / "cmd" / "__init__.py")
write_init(new_analysis / "__init__.py")


# ═══════════════════════════════════════════════════════════════════════════════
# 6. COPY GUI WIDGETS
# ═══════════════════════════════════════════════════════════════════════════════

print("Copying GUI widgets...")

copy_rewrite(
    APEX / "common" / "gui" / "widgets" / "image_viewer.py",
    APEX / "gui" / "widgets" / "image_viewer.py",
    []
)
write_init(APEX / "gui" / "widgets" / "__init__.py")
write_init(APEX / "gui" / "__init__.py")


# ═══════════════════════════════════════════════════════════════════════════════
# 7. COPY WORKFLOW FILES
# ═══════════════════════════════════════════════════════════════════════════════

print("Copying workflow shared files...")

# step_window_base.py: from common/gui/workflow/
copy_rewrite(
    APEX / "common" / "gui" / "workflow" / "step_window_base.py",
    APEX / "gui" / "workflow" / "step_window_base.py",
    []
)

# Shared steps 2-8: from common/gui/workflow/shared/
shared_src = APEX / "common" / "gui" / "workflow" / "shared"
shared_dst = APEX / "gui" / "workflow"

for fname in [
    "step2_crop_selector.py",
    "step3_sky_preview.py",
    "step4_source_detection.py",
    "step6_wcs_plate_solving.py",
    "step7_ref_build.py",
    "step8_star_id_matching.py",
]:
    copy_rewrite(shared_src / fname, shared_dst / fname, SHARED_WORKFLOW_EXTRA)

# step5_aperture_worker.py: standalone worker from shared/
copy_rewrite(
    shared_src / "step5_aperture_worker.py",
    shared_dst / "step5_aperture_worker.py",
    SHARED_WORKFLOW_EXTRA
)

# step5_aperture_photometry.py: FULL file from cmd/gui/workflow/ (NOT shared worker)
copy_rewrite(
    APEX / "cmd" / "gui" / "workflow" / "step5_aperture_photometry.py",
    shared_dst / "step5_aperture_photometry.py",
    CMD_WORKFLOW_EXTRA + [
        # The cmd version imports from apex.common.gui.workflow.shared.step5_aperture_worker
        ("from apex.common.gui.workflow.shared.step5_aperture_worker",
         "from apex.gui.workflow.step5_aperture_worker"),
    ]
)

write_init(APEX / "gui" / "workflow" / "__init__.py")

print("Copying CMD workflow files...")

# CMD-specific workflow files
cmd_wf_src = APEX / "cmd" / "gui" / "workflow"
cmd_wf_dst = APEX / "gui" / "workflow" / "cmd"

for src_name, dst_name in [
    ("step1_file_selection_window.py", "step1_file_selection.py"),
    ("step6_psf_photometry.py",        "step6_psf_photometry.py"),
    ("step10_master_id_editor.py",     "step10_master_id_editor.py"),
    ("step11_zeropoint_calibration.py","step11_zeropoint_calibration.py"),
    ("step12_cmd_plot.py",             "step12_cmd_plot.py"),
    ("step13_isochrone_model.py",      "step13_isochrone_model.py"),
]:
    copy_rewrite(cmd_wf_src / src_name, cmd_wf_dst / dst_name, CMD_WORKFLOW_EXTRA)

write_init(cmd_wf_dst / "__init__.py")

print("Copying LC workflow files...")

# LC-specific workflow files
lc_wf_src = APEX / "lightcurve" / "gui" / "workflow"
lc_wf_dst = APEX / "gui" / "workflow" / "lc"

for src_name, dst_name in [
    ("step1_file_selection_window.py",         "step1_file_selection.py"),
    ("step9_target_comparison_selection.py",   "step9_target_selection.py"),
    ("step10_light_curve_builder.py",          "step10_lightcurve_builder.py"),
    ("step11_detrend_merge.py",                "step11_detrend_merge.py"),
    ("step12_period_analysis.py",              "step12_period_analysis.py"),
]:
    copy_rewrite(lc_wf_src / src_name, lc_wf_dst / dst_name, LC_WORKFLOW_EXTRA)

write_init(lc_wf_dst / "__init__.py")


# ═══════════════════════════════════════════════════════════════════════════════
# 8. COPY TOOLS FILES
# ═══════════════════════════════════════════════════════════════════════════════

print("Copying tools files...")

lc_tools = APEX / "lightcurve" / "gui" / "tools"
lc_wf = APEX / "lightcurve" / "gui" / "workflow"
cmd_tools = APEX / "cmd" / "gui" / "tools"
new_tools = APEX / "gui" / "tools"

# LC-sourced tools (renamed: drop _window suffix)
lc_tool_map = [
    (lc_wf / "extinction_fit_window.py",       new_tools / "extinction_fit.py"),
    (lc_wf / "iraf_photometry_window.py",      new_tools / "iraf_photometry.py"),
    (lc_wf / "iraf_comparison_window.py",      new_tools / "iraf_comparison.py"),
    (lc_wf / "qa_report_window.py",            new_tools / "qa_report.py"),
    (lc_wf / "airmass_header_debug_tool.py",   new_tools / "airmass_debug.py"),
    (lc_wf / "aperture_overlay_tab.py",        new_tools / "aperture_overlay.py"),
    (lc_tools / "variable_star_tool.py",       new_tools / "variable_star.py"),
    (lc_tools / "multi_night_merger_tool.py",  new_tools / "multi_night_merger.py"),
    (lc_tools / "transit_tool.py",             new_tools / "transit_tool.py"),
    (lc_tools / "eb_tool.py",                  new_tools / "eb_tool.py"),
]

lc_tool_extra = LC_TOOLS_EXTRA + [
    # LC tools import from ..workflow.step12... fix
    ("from ..workflow.step12_period_analysis import PeriodAnalysisWorker",
     "from apex.gui.workflow.lc.step12_period_analysis import PeriodAnalysisWorker"),
    ("from ..workflow.step12_period_analysis", "from apex.gui.workflow.lc.step12_period_analysis"),
]

for src, dst in lc_tool_map:
    if src.exists():
        copy_rewrite(src, dst, lc_tool_extra)
    else:
        errors.append(f"MISSING SOURCE: {src}")

# CMD-specific tools
cmd_tool_map = [
    (cmd_tools / "gaia_3d_viewer_window.py",   new_tools / "gaia_3d_viewer.py"),
    (cmd_tools / "cmd_iso_tool_window.py",      new_tools / "cmd_iso_tool.py"),
]

for src, dst in cmd_tool_map:
    if src.exists():
        copy_rewrite(src, dst, CMD_TOOLS_EXTRA)
    else:
        errors.append(f"MISSING SOURCE: {src}")

# cluster_structure sub-package
cs_src = cmd_tools / "cluster_structure"
cs_dst = new_tools / "cluster_structure"
for fname in ["__init__.py", "analysis.py", "io.py", "window.py"]:
    src = cs_src / fname
    if src.exists():
        copy_rewrite(src, cs_dst / fname, CMD_CLUSTER_EXTRA + CMD_TOOLS_EXTRA)
    else:
        errors.append(f"MISSING SOURCE cluster_structure/{fname}: {src}")

write_init(new_tools / "__init__.py")


# ═══════════════════════════════════════════════════════════════════════════════
# 9. CREATE UNIFIED main_window.py
# ═══════════════════════════════════════════════════════════════════════════════

print("Creating unified main_window.py...")

main_window_content = '''"""
APEX Main Window — unified CMD + LC workflow.

Usage:
    MainWindowWorkflow(mode="cmd")  ← cluster photometry
    MainWindowWorkflow(mode="lc")   ← light curve analysis
"""

from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QGroupBox, QMessageBox, QFileDialog,
    QAction, QApplication, QLineEdit, QTextEdit,
    QPlainTextEdit, QComboBox, QSpinBox, QDoubleSpinBox
)
from PyQt5.QtCore import Qt, pyqtSignal, QObject, QEvent
from PyQt5.QtGui import QFont
from pathlib import Path
from typing import Optional, List


class StepButton(QPushButton):
    """Step button with completion/accessibility status indication."""

    def __init__(self, step_number: int, step_name: str, parent=None):
        super().__init__(parent)
        self.step_number = step_number
        self.step_name = step_name
        self.completed = False
        self.accessible = False
        self.setText(f"Step {step_number}: {step_name}")
        self.setMinimumHeight(50)
        self.setMinimumWidth(300)
        self.update_appearance()

    def set_completed(self, completed: bool):
        self.completed = completed
        self.update_appearance()

    def set_accessible(self, accessible: bool):
        self.accessible = accessible
        self.update_appearance()

    def update_appearance(self):
        if self.completed:
            self.setStyleSheet("""
                QPushButton {
                    background-color: #4CAF50; color: white;
                    font-size: 14px; font-weight: bold;
                    border: 2px solid #45a049; border-radius: 5px;
                    text-align: left; padding: 10px;
                }
                QPushButton:hover { background-color: #45a049; }
            """)
            self.setText(f"\\u2713 Step {self.step_number}: {self.step_name}")
            self.setEnabled(True)
        elif self.accessible:
            self.setStyleSheet("""
                QPushButton {
                    background-color: #2196F3; color: white;
                    font-size: 14px; font-weight: bold;
                    border: 2px solid #1976D2; border-radius: 5px;
                    text-align: left; padding: 10px;
                }
                QPushButton:hover { background-color: #1976D2; }
            """)
            self.setText(f"\\u25cb Step {self.step_number}: {self.step_name}")
            self.setEnabled(True)
        else:
            self.setStyleSheet("""
                QPushButton {
                    background-color: #E0E0E0; color: #999999;
                    font-size: 14px;
                    border: 2px solid #CCCCCC; border-radius: 5px;
                    text-align: left; padding: 10px;
                }
                QPushButton:disabled { background-color: #E0E0E0; color: #999999; }
            """)
            self.setText(f"\\U0001f512 Step {self.step_number}: {self.step_name} (Locked)")
            self.setEnabled(False)


class ShortcutRouter(QObject):
    def __init__(self, main_window: "MainWindowWorkflow"):
        super().__init__(main_window)
        self.main_window = main_window

    @staticmethod
    def _is_text_input(widget) -> bool:
        if isinstance(widget, (QLineEdit, QTextEdit, QPlainTextEdit)):
            return True
        if isinstance(widget, QComboBox) and widget.isEditable():
            return True
        return False

    def eventFilter(self, obj, event):
        if event.type() != QEvent.KeyPress:
            return False
        if self._is_text_input(QApplication.focusWidget()):
            return False
        key = event.key()
        if key not in (Qt.Key_Period, Qt.Key_BracketLeft, Qt.Key_BracketRight):
            return False
        target = self.main_window.current_step_window
        if target is None:
            target = QApplication.activeWindow()
        if target is None:
            return False
        if key == Qt.Key_Period and hasattr(target, "cycle_filter"):
            target.cycle_filter()
            return True
        if key == Qt.Key_BracketLeft:
            if hasattr(target, "navigate_frame"):
                target.navigate_frame(-1)
                return True
            if hasattr(target, "step_frame"):
                target.step_frame(-1)
                return True
        if key == Qt.Key_BracketRight:
            if hasattr(target, "navigate_frame"):
                target.navigate_frame(1)
                return True
            if hasattr(target, "step_frame"):
                target.step_frame(1)
                return True
        return False


class MainWindowWorkflow(QMainWindow):
    """Unified APEX main window for CMD (cluster) and LC (light curve) modes."""

    step_requested = pyqtSignal(int)
    log_message = pyqtSignal(str)

    def __init__(self, mode: str = "lc", param_file: Optional[str] = None):
        super().__init__()
        if mode not in ("cmd", "lc"):
            raise ValueError(f"mode must be \\'cmd\\' or \\'lc\\', got {mode!r}")
        self.mode = mode

        try:
            from apex.utils.param_loader import resolve_param_file
            if mode == "cmd":
                from apex.config.parameters_cmd import Parameters
            else:
                from apex.config.parameters_lc import Parameters
            self.params = Parameters(resolve_param_file(self, param_file))

            from apex.core import InstrumentConfig, FileManager, ProjectState
            self.instrument = InstrumentConfig(self.params)
            self.file_manager = FileManager(self.params)

            project_root = Path(__file__).parent.parent
            state_dir = project_root / ".state" / mode
            state_dir.mkdir(parents=True, exist_ok=True)
            self.project_state = ProjectState(state_dir)

            if mode == "lc":
                self._bootstrap_file_selection_state()

        except Exception as e:
            QMessageBox.critical(self, "Initialization Error",
                                 f"Failed to load parameters:\\n{e}")
            raise

        if mode == "cmd":
            self.step_names = [
                "File Selection",
                "Image Crop",
                "Sky Preview & QC",
                "Source Detection",
                "Aperture Photometry",
                "PSF Photometry",
                "WCS Plate Solving",
                "Reference Catalog Build",
                "Star ID Matching",
                "Master ID Editor",
                "Zeropoint Calibration",
                "CMD Plot",
                "Isochrone Model",
            ]
        else:
            self.step_names = [
                "File Selection",
                "Image Crop",
                "Sky Preview & QC",
                "Source Detection",
                "Aperture Photometry",
                "WCS Plate Solving",
                "Reference Build",
                "Star ID Matching",
                "Target/Comparison Selection",
                "Light Curve Builder",
                "Detrend & Night Merge",
                "Period Analysis",
            ]

        self.project_state.assign_steps(self.step_names)
        self.step_buttons: List[StepButton] = []
        self.current_step_window = None
        self.varstar_window = None

        self.setup_ui()
        self.setup_menu()
        self.update_step_buttons()

        mode_label = "CMD Cluster Photometry" if mode == "cmd" else "LC Light Curve Analysis"
        self.append_log(f"APEX {mode_label} initialized")
        self.append_log(f"Project: {self.project_state.state[\'project_name\']}")

        self._shortcut_router = ShortcutRouter(self)
        QApplication.instance().installEventFilter(self._shortcut_router)

        if mode == "lc" and hasattr(self, "_offline_data_dir"):
            QMessageBox.warning(
                self, "Previous Data Path Unavailable",
                f"The last-used data path is inaccessible:\\n\\n"
                f"  {self._offline_data_dir}\\n\\n"
                "External drive may be disconnected.\\n"
                "Please set a data path in Step 1."
            )

    # ── LC file-selection bootstrap ──────────────────────────────────────────

    def _bootstrap_file_selection_state(self) -> None:
        state_data = self.project_state.get_step_data("file_selection")
        if not state_data:
            return
        data_dir = state_data.get("data_dir")
        if data_dir:
            self.params.P.data_dir = Path(data_dir)
            saved_result_dir = state_data.get("result_dir")
            self.params.P.result_dir = (
                Path(saved_result_dir) if saved_result_dir
                else self.params.P.data_dir / "result"
            )
            self.params.P.cache_dir = self.params.P.result_dir / "cache"
            try:
                self.params.P.result_dir.mkdir(parents=True, exist_ok=True)
                self.params.P.cache_dir.mkdir(parents=True, exist_ok=True)
            except OSError:
                self._offline_data_dir = str(data_dir)
        prefix = state_data.get("filename_prefix")
        if prefix:
            self.params.P.filename_prefix = prefix
        ref_frame = state_data.get("reference_frame")
        if ref_frame:
            self.file_manager.ref_filename = ref_frame
        multi_night = bool(state_data.get("multi_night"))
        root_dir = state_data.get("root_dir") or data_dir
        night_dirs = [Path(p) for p in state_data.get("night_dirs", []) if p]
        if multi_night and night_dirs:
            root_path = Path(root_dir) if root_dir else self.params.P.data_dir
            self.file_manager.set_multi_night_dirs(root_path, night_dirs)
        else:
            self.file_manager.clear_multi_night_dirs()

    # ── UI setup ─────────────────────────────────────────────────────────────

    def setup_ui(self):
        mode_title = "CMD Cluster Photometry" if self.mode == "cmd" else "Light Curve Analysis"
        self.setWindowTitle(f"APEX — {mode_title}")
        self.setMinimumSize(800, 700)

        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        title = QLabel("Aperture Photometry Toolkit")
        title.setFont(QFont("Arial", 18, QFont.Bold))
        title.setAlignment(Qt.AlignCenter)
        layout.addWidget(title)

        subtitle = QLabel(f"KNUEMAO Observatory — {mode_title}")
        subtitle.setFont(QFont("Arial", 10))
        subtitle.setAlignment(Qt.AlignCenter)
        layout.addWidget(subtitle)

        settings_layout = QHBoxLayout()
        settings_layout.addStretch()
        btn_settings = QPushButton("\\u2699 Instrument Settings")
        btn_settings.setFont(QFont("Arial", 11, QFont.Bold))
        btn_settings.setMinimumHeight(40)
        btn_settings.setStyleSheet("""
            QPushButton {
                background-color: #9C27B0; color: white;
                border: 2px solid #7B1FA2; border-radius: 5px; padding: 5px 15px;
            }
            QPushButton:hover { background-color: #7B1FA2; }
        """)
        btn_settings.clicked.connect(self.open_settings)
        settings_layout.addWidget(btn_settings)
        settings_layout.addStretch()
        layout.addLayout(settings_layout)

        progress_group = QGroupBox("Workflow Progress")
        progress_layout = QVBoxLayout(progress_group)
        self.progress_label = QLabel(f"Progress: 0/{len(self.step_names)} steps finished")
        self.progress_label.setFont(QFont("Arial", 10, QFont.Bold))
        progress_layout.addWidget(self.progress_label)
        layout.addWidget(progress_group)

        steps_group = QGroupBox("Processing Steps")
        steps_layout = QVBoxLayout(steps_group)
        for i, step_name in enumerate(self.step_names):
            btn = StepButton(i + 1, step_name)
            btn.clicked.connect(lambda checked, idx=i: self.open_step(idx))
            self.step_buttons.append(btn)
            steps_layout.addWidget(btn)
        layout.addWidget(steps_group)

        action_layout = QHBoxLayout()
        btn_resume = QPushButton("Resume Next Step")
        btn_resume.setFont(QFont("Arial", 11, QFont.Bold))
        btn_resume.setMinimumHeight(40)
        btn_resume.clicked.connect(self.resume_next_step)
        action_layout.addWidget(btn_resume)
        btn_reset = QPushButton("Reset Progress")
        btn_reset.setMinimumHeight(40)
        btn_reset.clicked.connect(self.reset_progress)
        action_layout.addWidget(btn_reset)
        layout.addLayout(action_layout)

        log_group = QGroupBox("Activity Log")
        log_layout = QVBoxLayout(log_group)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setMaximumHeight(120)
        self.log_text.setFont(QFont("Courier", 8))
        log_layout.addWidget(self.log_text)
        layout.addWidget(log_group)

        self.status_bar = self.statusBar()
        self.status_bar.showMessage("Ready")

    # ── Menu setup ───────────────────────────────────────────────────────────

    def setup_menu(self):
        menubar = self.menuBar()

        file_menu = menubar.addMenu("&File")
        action_save = QAction("&Save Project State", self)
        action_save.setShortcut("Ctrl+S")
        action_save.triggered.connect(self.save_project_state)
        file_menu.addAction(action_save)
        action_export = QAction("&Export Summary...", self)
        action_export.triggered.connect(self.export_summary)
        file_menu.addAction(action_export)
        file_menu.addSeparator()
        action_exit = QAction("E&xit", self)
        action_exit.setShortcut("Ctrl+Q")
        action_exit.triggered.connect(self.close)
        file_menu.addAction(action_exit)

        tools_menu = menubar.addMenu("&Tools")
        action_params = QAction("View &Parameters", self)
        action_params.triggered.connect(self.show_parameters)
        tools_menu.addAction(action_params)
        tools_menu.addSeparator()

        action_qa = QAction("QA Reports", self)
        action_qa.setShortcut("Ctrl+R")
        action_qa.triggered.connect(self.open_qa_report)
        tools_menu.addAction(action_qa)

        action_iraf = QAction("IRAF/DAOPHOT Tool", self)
        action_iraf.setShortcut("Ctrl+I")
        action_iraf.triggered.connect(self.open_iraf_tool)
        tools_menu.addAction(action_iraf)

        if self.mode == "lc":
            action_ext = QAction("Extinction Fit Tool", self)
            action_ext.setShortcut("Ctrl+E")
            action_ext.triggered.connect(self.open_extinction_tool)
            tools_menu.addAction(action_ext)

            action_airmass = QAction("Airmass Header Debug", self)
            action_airmass.triggered.connect(self.open_airmass_debug_tool)
            tools_menu.addAction(action_airmass)

            tools_menu.addSeparator()

            action_merger = QAction("Multi-Night Light Curve Merger", self)
            action_merger.setShortcut("Ctrl+M")
            action_merger.triggered.connect(self.open_multi_night_merger)
            tools_menu.addAction(action_merger)

            tools_menu.addSeparator()

            action_varstar = QAction("Variable Star Analysis", self)
            action_varstar.setShortcut("Ctrl+Shift+V")
            action_varstar.triggered.connect(self.open_variable_star_tool)
            tools_menu.addAction(action_varstar)

            action_transit = QAction("Exoplanet Transit Analysis", self)
            action_transit.setShortcut("Ctrl+Shift+T")
            action_transit.triggered.connect(self.open_transit_tool)
            tools_menu.addAction(action_transit)

            action_eb = QAction("Eclipsing Binary Analysis", self)
            action_eb.setShortcut("Ctrl+Shift+B")
            action_eb.triggered.connect(self.open_eb_tool)
            tools_menu.addAction(action_eb)

        elif self.mode == "cmd":
            action_ext_fit = QAction("Extinction (Airmass Fit)", self)
            action_ext_fit.triggered.connect(self.open_extinction_fit_cmd)
            tools_menu.addAction(action_ext_fit)

            action_cmd_prev = QAction("CMD + Isochrone (From Results)...", self)
            action_cmd_prev.triggered.connect(self.open_cmd_iso_tool)
            tools_menu.addAction(action_cmd_prev)

            action_gaia_3d = QAction("Gaia 3D Cluster Viewer", self)
            action_gaia_3d.triggered.connect(self.open_gaia_3d_viewer)
            tools_menu.addAction(action_gaia_3d)

            action_cluster = QAction("Analyze Cluster Structure", self)
            action_cluster.triggered.connect(self.open_cluster_structure_tool)
            tools_menu.addAction(action_cluster)

        help_menu = menubar.addMenu("&Help")
        action_about = QAction("&About", self)
        action_about.triggered.connect(self.show_about)
        help_menu.addAction(action_about)

    # ── Step button state ────────────────────────────────────────────────────

    def update_step_buttons(self):
        completed_count = len(self.project_state.state["completed_steps"])
        self.progress_label.setText(
            f"Progress: {completed_count}/{len(self.step_names)} steps finished"
        )
        for i, btn in enumerate(self.step_buttons):
            completed = self.project_state.is_step_completed(i)
            accessible = self.project_state.is_step_accessible(i)
            btn.set_completed(completed)
            btn.set_accessible(accessible)
            btn.setEnabled(accessible)

    # ── Step dispatch ────────────────────────────────────────────────────────

    def open_step(self, step_index: int):
        if not self.project_state.is_step_accessible(step_index):
            prev_idx = step_index - 1
            prev_name = (self.step_names[prev_idx]
                         if 0 <= prev_idx < len(self.step_names) else "previous step")
            QMessageBox.warning(self, "Step Not Accessible",
                                f"Please finish Step {step_index}: {prev_name} first.")
            return

        self.project_state.set_current_step(step_index)
        self.append_log(f"Opening Step {step_index + 1}: {self.step_names[step_index]}")

        if self.current_step_window:
            self.current_step_window.close()

        win = self._open_step_window(step_index)
        if win is None:
            QMessageBox.information(self, "Step Not Implemented",
                                    f"Step {step_index + 1} is not yet implemented.")
            return

        self.current_step_window = win
        win.show()

    def _open_step_window(self, step_index: int):  # noqa: C901 (complexity ok)
        p, fm, ps = self.params, self.file_manager, self.project_state

        # ── Step 0: File selection (mode-specific) ──
        if step_index == 0:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step1_file_selection import FileSelectionWindow
            else:
                from apex.gui.workflow.lc.step1_file_selection import FileSelectionWindow
            return FileSelectionWindow(p, fm, ps, self)

        # ── Steps 1-4: shared ──
        elif step_index == 1:
            from apex.gui.workflow.step2_crop_selector import CropSelectorWindow
            return CropSelectorWindow(p, fm, ps, self)
        elif step_index == 2:
            from apex.gui.workflow.step3_sky_preview import SkyPreviewWindow
            return SkyPreviewWindow(p, fm, ps, self)
        elif step_index == 3:
            from apex.gui.workflow.step4_source_detection import SourceDetectionWindow
            return SourceDetectionWindow(p, fm, ps, self)
        elif step_index == 4:
            from apex.gui.workflow.step5_aperture_photometry import AperturePhotometryWindow
            return AperturePhotometryWindow(p, fm, ps, self)

        # ── Step 5: PSF (CMD only) / WCS (LC) ──
        elif step_index == 5:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step6_psf_photometry import PSFPhotometryWindow
                return PSFPhotometryWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.step6_wcs_plate_solving import WcsPlateSolvingWindow
                return WcsPlateSolvingWindow(p, fm, ps, self)

        # ── Steps 6-7: WCS/Ref in CMD mode; Ref/IDMatch in LC ──
        elif step_index == 6:
            if self.mode == "cmd":
                from apex.gui.workflow.step6_wcs_plate_solving import WcsPlateSolvingWindow
                return WcsPlateSolvingWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.step7_ref_build import RefBuildWindow
                return RefBuildWindow(p, fm, ps, self)

        elif step_index == 7:
            if self.mode == "cmd":
                from apex.gui.workflow.step7_ref_build import RefBuildWindow
                return RefBuildWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.step8_star_id_matching import StarIdMatchingWindow
                return StarIdMatchingWindow(p, fm, ps, self)

        # ── Step 8: IDMatch (CMD) / Target selection (LC) ──
        elif step_index == 8:
            if self.mode == "cmd":
                from apex.gui.workflow.step8_star_id_matching import StarIdMatchingWindow
                return StarIdMatchingWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step9_target_selection import TargetComparisonSelectionWindow
                return TargetComparisonSelectionWindow(p, fm, ps, self)

        # ── CMD steps 9-12 ──
        elif step_index == 9:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step10_master_id_editor import MasterIdEditorWindow
                return MasterIdEditorWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step10_lightcurve_builder import LightCurveBuilderWindow
                return LightCurveBuilderWindow(p, fm, ps, self)

        elif step_index == 10:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step11_zeropoint_calibration import ZeropointCalibrationWindow
                return ZeropointCalibrationWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step11_detrend_merge import DetrendNightMergeWindow
                return DetrendNightMergeWindow(p, fm, ps, self)

        elif step_index == 11:
            if self.mode == "cmd":
                from apex.gui.workflow.cmd.step12_cmd_plot import CmdPlotWindow
                return CmdPlotWindow(p, fm, ps, self)
            else:
                from apex.gui.workflow.lc.step12_period_analysis import PeriodAnalysisWindow
                return PeriodAnalysisWindow(p, fm, ps, self)

        elif step_index == 12 and self.mode == "cmd":
            from apex.gui.workflow.cmd.step13_isochrone_model import IsochroneModelWindow
            return IsochroneModelWindow(p, fm, ps, self)

        return None

    def on_step_completed(self, step_index: int):
        self.project_state.mark_step_completed(step_index)
        self.update_step_buttons()
        self.append_log(f"\\u2713 Step {step_index + 1} finished: {self.step_names[step_index]}")

    def resume_next_step(self):
        next_step = self.project_state.get_next_incomplete_step()
        if next_step is not None:
            if not self.project_state.is_step_accessible(next_step):
                prev_idx = next_step - 1
                prev_name = (self.step_names[prev_idx]
                             if 0 <= prev_idx < len(self.step_names) else "previous step")
                QMessageBox.warning(self, "Step Not Accessible",
                                    f"Please finish earlier steps first.\\n"
                                    f"Required now: Step {next_step}: {prev_name}")
                return
            self.open_step(next_step)
        else:
            QMessageBox.information(self, "Workflow Finished",
                                    "All workflow steps are finished.")

    def reset_progress(self):
        reply = QMessageBox.question(
            self, "Reset Progress",
            "Are you sure you want to reset all progress?\\n"
            "This will clear completion status but keep your data files.",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No
        )
        if reply == QMessageBox.Yes:
            self.project_state.reset()
            self.update_step_buttons()
            self.append_log("Progress reset")

    def save_project_state(self):
        self.project_state.save()
        self.append_log("Project state saved")
        QMessageBox.information(self, "Saved", "Project state saved successfully.")

    def export_summary(self):
        summary = self.project_state.export_summary()
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Summary",
            str(self.params.P.result_dir / "project_summary.txt"),
            "Text Files (*.txt)"
        )
        if file_path:
            Path(file_path).write_text(summary, encoding="utf-8")
            self.append_log(f"Summary exported to {file_path}")

    def show_parameters(self):
        self.params.print_summary()
        QMessageBox.information(self, "Parameters",
                                "Parameter summary printed to console.")

    def append_log(self, message: str):
        from datetime import datetime
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.append(f"[{timestamp}] {message}")

    def show_about(self):
        mode_label = "CMD Cluster Photometry" if self.mode == "cmd" else "LC Light Curve Analysis"
        QMessageBox.about(
            self, "About APEX",
            f"<h2>APEX — Automated Photometry EXtraction</h2>"
            f"<p><b>Mode: {mode_label}</b></p>"
            "<p>KNUEMAO Observatory — CDK500 + Moravian C3-61000</p>"
            "<p>Version 2.0.0</p>"
        )

    # ── Tool launchers ───────────────────────────────────────────────────────

    def open_qa_report(self, tab: int = 0):
        from apex.gui.tools.qa_report import QAReportWindow
        self.qa_window = QAReportWindow(self.params, self.params.P.result_dir, parent=None)
        self.qa_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        if hasattr(self.qa_window, "tabs") and tab >= 0:
            self.qa_window.tabs.setCurrentIndex(
                min(tab, self.qa_window.tabs.count() - 1))
        self.qa_window.show()
        self.qa_window.raise_()
        self.qa_window.activateWindow()
        self.append_log("Opened QA Report window")

    def open_iraf_tool(self):
        from apex.gui.tools.iraf_photometry import IRAFPhotometryWindow
        self.iraf_window = IRAFPhotometryWindow(
            self.params, self.params.P.data_dir, self.params.P.result_dir,
            self.project_state, parent=None
        )
        self.iraf_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.iraf_window.show()
        self.iraf_window.raise_()
        self.iraf_window.activateWindow()
        self.append_log("Opened IRAF/DAOPHOT Tool")

    # ── LC tools ─────────────────────────────────────────────────────────────

    def open_extinction_tool(self):
        from apex.gui.tools.extinction_fit import ExtinctionFitWindow
        self.extinction_window = ExtinctionFitWindow(
            self.params, self.params.P.data_dir, self.params.P.result_dir, parent=None
        )
        self.extinction_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.extinction_window.show()
        self.extinction_window.raise_()
        self.extinction_window.activateWindow()
        self.append_log("Opened Extinction (Airmass Fit) Tool")

    def open_airmass_debug_tool(self):
        from apex.gui.tools.airmass_debug import AirmassHeaderDebugToolWindow
        self.airmass_debug_window = AirmassHeaderDebugToolWindow(
            self.params, self.project_state, parent=None, file_manager=self.file_manager
        )
        self.airmass_debug_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.airmass_debug_window.show()
        self.airmass_debug_window.raise_()
        self.airmass_debug_window.activateWindow()
        self.append_log("Opened Airmass Header Debug Tool")

    def open_multi_night_merger(self):
        from apex.gui.tools.multi_night_merger import MultiNightMergerWindow
        self.merger_window = MultiNightMergerWindow(
            self.params, self.project_state, main_window=self
        )
        self.merger_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.hide()
        self.merger_window.show()
        self.merger_window.raise_()
        self.merger_window.activateWindow()
        self.append_log("Opened Multi-Night Light Curve Merger")

    def open_variable_star_tool(self):
        from apex.gui.tools.variable_star import VariableStarToolWindow
        if self.varstar_window is None:
            self.varstar_window = VariableStarToolWindow(
                self.params, self.project_state, parent=None
            )
            self.varstar_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.varstar_window.show()
        self.varstar_window.raise_()
        self.varstar_window.activateWindow()
        self.append_log("Opened Variable Star Analysis Tool")

    def open_transit_tool(self):
        from apex.gui.tools.transit_tool import TransitToolWindow
        self.transit_window = TransitToolWindow(
            self.params, self.project_state, parent=None
        )
        self.transit_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.transit_window.show()
        self.transit_window.raise_()
        self.transit_window.activateWindow()
        self.append_log("Opened Exoplanet Transit Analysis Tool")

    def open_eb_tool(self):
        from apex.gui.tools.eb_tool import EclipsingBinaryToolWindow
        self.eb_window = EclipsingBinaryToolWindow(
            self.params, self.project_state, parent=None
        )
        self.eb_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.eb_window.show()
        self.eb_window.raise_()
        self.eb_window.activateWindow()
        self.append_log("Opened Eclipsing Binary Analysis Tool")

    # ── CMD tools ────────────────────────────────────────────────────────────

    def open_extinction_fit_cmd(self):
        from apex.gui.tools.extinction_fit import ExtinctionFitWindow
        self.ext_window = ExtinctionFitWindow(
            self.params, self.params.P.data_dir, self.params.P.result_dir, parent=None
        )
        self.ext_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.ext_window.show()
        self.ext_window.raise_()
        self.ext_window.activateWindow()
        self.append_log("Opened Extinction Fit window")

    def open_cmd_iso_tool(self):
        start_dir = str(getattr(self.params.P, "result_dir", Path.cwd()))
        selected = QFileDialog.getExistingDirectory(self, "Select Result Folder", start_dir)
        if not selected:
            return
        from apex.gui.tools.cmd_iso_tool import CmdIsoToolWindow
        self.cmd_tool_window = CmdIsoToolWindow(
            self.params, self.file_manager, self.project_state, self,
            initial_result_dir=Path(selected), parent=None
        )
        self.cmd_tool_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.cmd_tool_window.show()
        self.cmd_tool_window.raise_()
        self.cmd_tool_window.activateWindow()
        self.append_log(f"Opened CMD + Isochrone tool: {selected}")

    def open_gaia_3d_viewer(self):
        from apex.gui.tools.gaia_3d_viewer import Gaia3DViewerWindow
        self.gaia_3d_window = Gaia3DViewerWindow(
            self.params, self.params.P.result_dir, parent=None
        )
        self.gaia_3d_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.gaia_3d_window.show()
        self.gaia_3d_window.raise_()
        self.gaia_3d_window.activateWindow()
        self.append_log(f"Opened Gaia 3D Viewer: {self.params.P.result_dir}")

    def open_cluster_structure_tool(self):
        from apex.gui.tools.cluster_structure import ClusterStructureWindow
        self.cluster_structure_window = ClusterStructureWindow(
            self.params, self.params.P.result_dir, parent=None
        )
        self.cluster_structure_window.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
        self.cluster_structure_window.show()
        self.cluster_structure_window.raise_()
        self.cluster_structure_window.activateWindow()
        self.append_log(f"Opened Cluster Structure Tool: {self.params.P.result_dir}")

    # ── Instrument settings ──────────────────────────────────────────────────

    def open_settings(self):
        from PyQt5.QtWidgets import (
            QDialog, QVBoxLayout, QLabel, QLineEdit,
            QDialogButtonBox, QGroupBox, QFormLayout
        )
        dialog = QDialog(self)
        dialog.setWindowTitle("Instrument Settings")
        dialog.setMinimumWidth(650)
        layout = QVBoxLayout(dialog)
        title = QLabel("Instrument Configuration")
        title.setFont(QFont("Arial", 12, QFont.Bold))
        layout.addWidget(title)

        tel_group = QGroupBox("Telescope")
        tel_layout = QFormLayout(tel_group)
        tel_name_edit = QLineEdit(self.instrument.telescope_name)
        tel_layout.addRow("Name:", tel_name_edit)
        tel_aperture_edit = QLineEdit(str(self.instrument.aperture_mm))
        tel_aperture_edit.setMaximumWidth(150)
        tel_layout.addRow("Aperture (mm):", tel_aperture_edit)
        tel_focal_edit = QLineEdit(str(self.instrument.focal_length_mm))
        tel_focal_edit.setMaximumWidth(150)
        tel_layout.addRow("Focal Length (mm):", tel_focal_edit)
        layout.addWidget(tel_group)

        cam_group = QGroupBox("Camera")
        cam_layout = QFormLayout(cam_group)
        cam_name_edit = QLineEdit(self.instrument.camera_name)
        cam_layout.addRow("Name:", cam_name_edit)
        cam_pixsize_edit = QLineEdit(str(self.instrument.pix_size_um))
        cam_pixsize_edit.setMaximumWidth(150)
        cam_layout.addRow("Pixel Size (\\u03bcm):", cam_pixsize_edit)
        cam_nx_edit = QLineEdit(str(self.instrument.sensor_nx_1x))
        cam_nx_edit.setMaximumWidth(150)
        cam_layout.addRow("Sensor Width (px):", cam_nx_edit)
        cam_ny_edit = QLineEdit(str(self.instrument.sensor_ny_1x))
        cam_ny_edit.setMaximumWidth(150)
        cam_layout.addRow("Sensor Height (px):", cam_ny_edit)
        cam_binning_edit = QLineEdit(str(self.instrument.binning))
        cam_binning_edit.setMaximumWidth(150)
        cam_layout.addRow("Binning:", cam_binning_edit)
        layout.addWidget(cam_group)

        params_group = QGroupBox("Camera Parameters")
        params_layout = QFormLayout(params_group)
        gain_edit = QLineEdit(str(self.params.P.gain_e_per_adu))
        gain_edit.setMaximumWidth(150)
        params_layout.addRow("Gain (e-/ADU):", gain_edit)
        rdnoise_edit = QLineEdit(str(self.params.P.rdnoise_e))
        rdnoise_edit.setMaximumWidth(150)
        params_layout.addRow("Read Noise (e-):", rdnoise_edit)
        saturation_edit = QLineEdit(str(self.params.P.saturation_adu))
        saturation_edit.setMaximumWidth(150)
        params_layout.addRow("Saturation (ADU):", saturation_edit)
        layout.addWidget(params_group)

        parallel_group = QGroupBox("Parallel Processing")
        parallel_layout = QFormLayout(parallel_group)
        parallel_workers_spin = QSpinBox()
        parallel_workers_spin.setRange(0, 16)
        parallel_workers_spin.setValue(int(getattr(self.params.P, "max_workers", 0)))
        parallel_workers_spin.setToolTip("0 = auto (use ~75% of CPU cores)")
        parallel_layout.addRow("Max Workers (0=auto):", parallel_workers_spin)
        layout.addWidget(parallel_group)

        site_group = QGroupBox("Observatory Location")
        site_layout = QFormLayout(site_group)
        site_lat_edit = QLineEdit(str(getattr(self.params.P, "site_lat_deg", 0.0)))
        site_lat_edit.setMaximumWidth(150)
        site_layout.addRow("Latitude (deg):", site_lat_edit)
        site_lon_edit = QLineEdit(str(getattr(self.params.P, "site_lon_deg", 0.0)))
        site_lon_edit.setMaximumWidth(150)
        site_layout.addRow("Longitude (deg):", site_lon_edit)
        site_alt_edit = QLineEdit(str(getattr(self.params.P, "site_alt_m", 0.0)))
        site_alt_edit.setMaximumWidth(150)
        site_layout.addRow("Altitude (m):", site_alt_edit)
        site_tz_edit = QLineEdit(str(getattr(self.params.P, "site_tz_offset_hours", 0.0)))
        site_tz_edit.setMaximumWidth(150)
        site_layout.addRow("UTC Offset (hours):", site_tz_edit)
        layout.addWidget(site_group)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dialog.accept)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        if dialog.exec_() == QDialog.Accepted:
            try:
                self.instrument.telescope_name = tel_name_edit.text().strip()
                self.instrument.aperture_mm = float(tel_aperture_edit.text())
                self.instrument.focal_length_mm = float(tel_focal_edit.text())
                self.instrument.focal_ratio = self.instrument.focal_length_mm / self.instrument.aperture_mm
                self.instrument.camera_name = cam_name_edit.text().strip()
                self.instrument.pix_size_um = float(cam_pixsize_edit.text())
                self.instrument.sensor_nx_1x = int(cam_nx_edit.text())
                self.instrument.sensor_ny_1x = int(cam_ny_edit.text())
                self.instrument.binning = int(cam_binning_edit.text())
                self.params.P.gain_e_per_adu = float(gain_edit.text())
                self.params.P.rdnoise_e = float(rdnoise_edit.text())
                self.params.P.saturation_adu = float(saturation_edit.text())
                self.params.P.site_lat_deg = float(site_lat_edit.text())
                self.params.P.site_lon_deg = float(site_lon_edit.text())
                self.params.P.site_alt_m = float(site_alt_edit.text())
                self.params.P.site_tz_offset_hours = float(site_tz_edit.text())
                self.params.P.max_workers = int(parallel_workers_spin.value())
                if not self.params.save_toml():
                    self.append_log("Warning: could not save settings to TOML")
                QMessageBox.information(self, "Settings Saved",
                                        "Instrument settings updated.")
                self.append_log("Instrument settings updated")
            except ValueError as e:
                QMessageBox.warning(self, "Invalid Input",
                                    f"Please enter valid numeric values.\\n{e}")

    # ── Window close ─────────────────────────────────────────────────────────

    def closeEvent(self, event):
        self.project_state.save()
        if self.current_step_window:
            self.current_step_window.close()
        event.accept()
'''

(APEX / "gui" / "main_window.py").write_text(main_window_content, encoding="utf-8")
created.append(str(APEX / "gui" / "main_window.py"))


# ═══════════════════════════════════════════════════════════════════════════════
# 10. UPDATE ENTRY POINTS
# ═══════════════════════════════════════════════════════════════════════════════

print("Updating entry points...")

cmd_main_content = '''#!/usr/bin/env python3
"""APEX CMD — Cluster Photometry entry point."""

import sys
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message=".*datfix.*MJD-OBS.*", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning, module="astropy")
try:
    from astropy.wcs import FITSFixedWarning
    warnings.filterwarnings("ignore", category=FITSFixedWarning)
except ImportError:
    pass
warnings.filterwarnings("ignore", message=".*tight_layout.*", category=UserWarning)

_HERE = Path(__file__).resolve().parent          # apex/cmd/
_ROOT = _HERE.parent.parent                       # Automated_Photometry_EXtraction/
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from PyQt5.QtWidgets import QApplication, QMessageBox
from PyQt5.QtGui import QFont, QFontDatabase


def _configure_fonts(app: QApplication) -> None:
    candidate_font_files = [
        Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / "malgun.ttf",
        Path("/mnt/c/Windows/Fonts/malgun.ttf"),
    ]
    candidate_families = [
        "Malgun Gothic", "맑은 고딕", "NanumGothic", "Nanum Gothic",
        "Noto Sans CJK KR", "Noto Sans KR", "AppleGothic", "Arial Unicode MS",
    ]
    resolved_family = None
    try:
        import matplotlib
        from matplotlib import font_manager
        for fp in candidate_font_files:
            if not fp.exists():
                continue
            try:
                font_manager.fontManager.addfont(str(fp))
                resolved_family = font_manager.FontProperties(fname=str(fp)).get_name()
                break
            except Exception:
                continue
        if resolved_family is None:
            available = {f.name.lower(): f.name for f in font_manager.fontManager.ttflist}
            for fam in candidate_families:
                if fam.lower() in available:
                    resolved_family = available[fam.lower()]
                    break
        if resolved_family:
            matplotlib.rcParams.update({
                "font.family": [resolved_family, "DejaVu Sans"],
                "axes.unicode_minus": False,
            })
    except Exception:
        pass
    try:
        db = QFontDatabase()
        qt_fams = {f.lower(): f for f in db.families()}
        qt_fam = resolved_family and qt_fams.get(resolved_family.lower())
        if qt_fam is None:
            for fam in candidate_families:
                qt_fam = qt_fams.get(fam.lower())
                if qt_fam:
                    break
        if qt_fam:
            app.setFont(QFont(qt_fam, 9))
    except Exception:
        pass


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("APEX CMD")
    app.setOrganizationName("APEX Project")
    _configure_fonts(app)
    os.chdir(_HERE)

    try:
        from apex.gui.main_window import MainWindowWorkflow
        window = MainWindowWorkflow(mode="cmd")
        window.show()
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"Failed to start APEX CMD: {e}")
        print(tb)
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle("APEX CMD — Startup Error")
        msg.setText(str(e))
        msg.setDetailedText(tb)
        msg.exec_()
        return 1

    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
'''

(APEX / "cmd" / "main.py").write_text(cmd_main_content, encoding="utf-8")
created.append(str(APEX / "cmd" / "main.py"))

lc_main_content = '''#!/usr/bin/env python3
"""APEX LC — Light Curve Analysis entry point."""

import sys
import os
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message=".*datfix.*MJD-OBS.*", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning, module="astropy")
try:
    from astropy.wcs import FITSFixedWarning
    warnings.filterwarnings("ignore", category=FITSFixedWarning)
except ImportError:
    pass
warnings.filterwarnings("ignore", message=".*tight_layout.*", category=UserWarning)

_HERE = Path(__file__).resolve().parent          # apex/lightcurve/
_ROOT = _HERE.parent.parent                       # Automated_Photometry_EXtraction/
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from PyQt5.QtWidgets import QApplication, QMessageBox
from PyQt5.QtGui import QFont, QFontDatabase


def _configure_fonts(app: QApplication) -> None:
    candidate_font_files = [
        Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / "malgun.ttf",
        Path("/mnt/c/Windows/Fonts/malgun.ttf"),
    ]
    candidate_families = [
        "Malgun Gothic", "맑은 고딕", "NanumGothic", "Nanum Gothic",
        "Noto Sans CJK KR", "Noto Sans KR", "AppleGothic", "Arial Unicode MS",
    ]
    resolved_family = None
    try:
        import matplotlib
        from matplotlib import font_manager
        for fp in candidate_font_files:
            if not fp.exists():
                continue
            try:
                font_manager.fontManager.addfont(str(fp))
                resolved_family = font_manager.FontProperties(fname=str(fp)).get_name()
                break
            except Exception:
                continue
        if resolved_family is None:
            available = {f.name.lower(): f.name for f in font_manager.fontManager.ttflist}
            for fam in candidate_families:
                if fam.lower() in available:
                    resolved_family = available[fam.lower()]
                    break
        if resolved_family:
            matplotlib.rcParams.update({
                "font.family": [resolved_family, "DejaVu Sans"],
                "axes.unicode_minus": False,
            })
    except Exception:
        pass
    try:
        db = QFontDatabase()
        qt_fams = {f.lower(): f for f in db.families()}
        qt_fam = resolved_family and qt_fams.get(resolved_family.lower())
        if qt_fam is None:
            for fam in candidate_families:
                qt_fam = qt_fams.get(fam.lower())
                if qt_fam:
                    break
        if qt_fam:
            app.setFont(QFont(qt_fam, 9))
    except Exception:
        pass


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("APEX LC")
    app.setOrganizationName("APEX Project")
    _configure_fonts(app)
    os.chdir(_HERE)

    try:
        from apex.gui.main_window import MainWindowWorkflow
        window = MainWindowWorkflow(mode="lc")
        window.show()
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"Failed to start APEX LC: {e}")
        print(tb)
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle("APEX LC — Startup Error")
        msg.setText(str(e))
        msg.setDetailedText(tb)
        msg.exec_()
        return 1

    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
'''

(APEX / "lightcurve" / "main.py").write_text(lc_main_content, encoding="utf-8")
created.append(str(APEX / "lightcurve" / "main.py"))


# ═══════════════════════════════════════════════════════════════════════════════
# 11. CREATE DOCUMENTATION
# ═══════════════════════════════════════════════════════════════════════════════

print("Creating documentation...")

readme_content = '''# APEX — Automated Photometry EXtraction

APEX is a PyQt5-based GUI toolkit for aperture and PSF photometry of astronomical images. It supports two operational modes:

- **CMD mode** (`apex/cmd/`): Cluster photometry pipeline — detection through CMD diagram and isochrone fitting (13 steps).
- **LC mode** (`apex/lightcurve/`): Light curve analysis pipeline — multi-night photometry, detrending, and period analysis (12 steps).

Both modes share a common core (steps 1–8: file selection, crop, sky preview, source detection, aperture photometry, WCS plate solving, reference catalog build, star ID matching) and diverge at step 9.

## Requirements

- Python 3.10+
- PyQt5
- astropy >= 5.0
- photutils >= 1.5
- numpy
- pandas
- scipy
- matplotlib
- tomli / tomllib (Python 3.11+ built-in)
- astroquery (Gaia access)
- astrometry.net (local installation for WCS solving)

## Installation

```bash
cd /path/to/Automated_Photometry_EXtraction
pip install -r requirements.txt
```

## Quickstart

```bash
# CMD mode (cluster photometry)
python apex/cmd/main.py

# LC mode (light curve analysis)
python apex/lightcurve/main.py
```

Or from the project root launcher:
```bash
python main.py
```

## Screenshot

<!-- TODO: add screenshot -->

## License

<!-- TODO: add license -->
'''

(ROOT / "README.md").write_text(readme_content, encoding="utf-8")
created.append(str(ROOT / "README.md"))

arch_content = '''# APEX Architecture

## Package Overview

```
apex/
  core/                    Shared core: ProjectState, FileManager, InstrumentConfig
  config/
    parameters_cmd.py      Parameters class for CMD mode (TOML-backed)
    parameters_lc.py       Parameters class for LC mode (TOML-backed)
    schema.py              TOML schema validators
  utils/
    step_paths.py          Shared step dir helpers (steps 1-8)
    step_paths_cmd.py      CMD-specific step dir helpers (steps 9-13)
    step_paths_lc.py       LC-specific step dir helpers (steps 9-12)
    photometry_utils.py    Aperture photometry utilities
    astro_utils.py         Airmass, BJD, filter normalization
    io_utils.py            CSV/ECSV int64 safe readers
    cache_utils.py         Detection/WCS cache management
    param_loader.py        TOML parameter file resolution
    photometry_loader.py   Frame photometry CSV loader (LC)
    run_workspace.py       Multi-night workspace helpers (LC)
    ... (logging, constants, qc, header_cache, common_helpers)
  analysis/
    light_curve/           LC science: lightcurve, detrend, period, eclipse, asteroid
    merge/                 Multi-night workspace scan/build/id-match (LC)
    cmd/                   CMD science: isochrone_fitter, isochrone_fitter_v2
  gui/
    main_window.py         UNIFIED main window (mode="cmd" or "lc")
    widgets/
      image_viewer.py      Zoomable FITS image viewer
    workflow/
      step_window_base.py  Base class for all step windows
      step2_crop_selector.py   }
      step3_sky_preview.py     }  Shared steps 2-8
      step4_source_detection.py}  (same window for both modes)
      step5_aperture_worker.py }
      step5_aperture_photometry.py}
      step6_wcs_plate_solving.py }
      step7_ref_build.py        }
      step8_star_id_matching.py }
      cmd/
        step1_file_selection.py   CMD file selection
        step6_psf_photometry.py   PSF photometry (CMD step 5, index 5)
        step10_master_id_editor.py
        step11_zeropoint_calibration.py
        step12_cmd_plot.py
        step13_isochrone_model.py
      lc/
        step1_file_selection.py   LC file selection (multi-night aware)
        step9_target_selection.py
        step10_lightcurve_builder.py
        step11_detrend_merge.py
        step12_period_analysis.py
    tools/
      extinction_fit.py     Bouguer extinction / zeropoint fitting
      iraf_photometry.py    IRAF/DAOPHOT integration
      iraf_comparison.py    IRAF comparison photometry
      qa_report.py          QA / publication validation report
      airmass_debug.py      Airmass header diagnostics
      aperture_overlay.py   Aperture overlay visualizer
      variable_star.py      Variable star classification (LC)
      multi_night_merger.py Multi-night LC merge tool (LC)
      transit_tool.py       Exoplanet transit fitting (LC)
      eb_tool.py            Eclipsing binary fitting (LC)
      gaia_3d_viewer.py     Gaia 3D cluster viewer (CMD)
      cmd_iso_tool.py       CMD + isochrone from results (CMD)
      cluster_structure/    Cluster structure analysis (CMD)
```

## Mode Concept

Both modes launch from `MainWindowWorkflow(mode=...)`:

- **Shared steps (1-8)** use the same window classes regardless of mode.
- **CMD-only steps**: PSF photometry (step 5), master ID editor, zeropoint, CMD plot, isochrone.
- **LC-only steps**: target/comparison selection, light curve builder, detrend/merge, period analysis.

Step index → file dispatch is in `_open_step_window(step_index)` of `main_window.py`.

## Adding a New Shared Step

1. Create `apex/gui/workflow/step_N_xxx.py` inheriting `StepWindowBase`.
2. Add a `stepN_xxx_dir()` helper to `apex/utils/step_paths.py`.
3. Wire up in `main_window._open_step_window()` for both modes.

## Adding a New Mode-Specific Step

1. Create `apex/gui/workflow/cmd/stepN_xxx.py` or `apex/gui/workflow/lc/stepN_xxx.py`.
2. Add path helpers to `step_paths_cmd.py` or `step_paths_lc.py`.
3. Add the step name to the appropriate `step_names` list and wire dispatch in `main_window.py`.

## Data Flow

```
params.P.data_dir/          Raw FITS files
params.P.result_dir/
  step1_file_selection/     FITS scan manifest
  step2_crop/               Crop region + cropped images
  step3_sky_preview/        Sky QC metadata
  step4_detection/          Source catalogs + frame QC
  step5_aperture/           Aperture photometry CSVs
  step6_wcs/                WCS-solved FITS headers
  step7_refbuild/           Master star catalog, Gaia IDs
  step8_idmatch/            Per-frame star ID matches
  [cmd_*/ or lc_*/]         Mode-specific outputs
  cache/                    Intermediate caches (header scan, detect, WCS)
```

## Step Directory Conventions

Each step writes to a named subdirectory of `result_dir`, defined in `step_paths*.py`.
- Shared names (step1-8): defined in `step_paths.py`.
- CMD names (cmd_psf, cmd_selection, cmd_zeropoint, cmd_plot, cmd_isochrone): in `step_paths_cmd.py`.
- LC names (lc_selection, lc_lightcurve, lc_detrend, lc_period): in `step_paths_lc.py`.
'''

(ROOT / "ARCHITECTURE.md").write_text(arch_content, encoding="utf-8")
created.append(str(ROOT / "ARCHITECTURE.md"))


# ═══════════════════════════════════════════════════════════════════════════════
# 12. SYNTAX CHECK ALL NEW FILES
# ═══════════════════════════════════════════════════════════════════════════════

print("\nRunning syntax check on all new files...")

import glob

new_py_files = glob.glob(
    str(APEX / "**" / "*.py"), recursive=True
)

# Only check files we just created (in new directories)
new_dirs_str = [
    str(APEX / "core"),
    str(APEX / "config"),
    str(APEX / "utils"),
    str(APEX / "analysis"),
    str(APEX / "gui"),
]

syntax_errors = []
syntax_ok = 0

for fpath in sorted(new_py_files):
    if "__pycache__" in fpath:
        continue
    # Only check files in the new flat directories
    is_new = any(fpath.startswith(nd) for nd in new_dirs_str)
    # Also check the updated entry points
    is_entrypoint = fpath in [
        str(APEX / "cmd" / "main.py"),
        str(APEX / "lightcurve" / "main.py"),
    ]
    if not is_new and not is_entrypoint:
        continue
    try:
        with open(fpath, encoding="utf-8") as f:
            source = f.read()
        ast.parse(source, filename=fpath)
        syntax_ok += 1
    except SyntaxError as e:
        syntax_errors.append(f"SYNTAX ERROR in {fpath}: {e}")
    except Exception as e:
        syntax_errors.append(f"READ ERROR in {fpath}: {e}")

print(f"Syntax OK: {syntax_ok} files")
if syntax_errors:
    print(f"Syntax ERRORS: {len(syntax_errors)} files")
    for e in syntax_errors:
        print(f"  {e}")


# ═══════════════════════════════════════════════════════════════════════════════
# 13. DELETION OF OLD STRUCTURE
# ═══════════════════════════════════════════════════════════════════════════════

if syntax_errors:
    print("\nSkipping deletion: fix syntax errors first.")
else:
    print("\nDeleting old structure...")

    # Delete apex/common/ entirely
    common_dir = APEX / "common"
    if common_dir.exists():
        shutil.rmtree(common_dir)
        print(f"  Deleted {common_dir}")

    # From apex/cmd/ — keep __init__.py and main.py (already updated), delete the rest
    for subdir in ["config", "core", "utils", "analysis"]:
        d = APEX / "cmd" / subdir
        if d.exists():
            shutil.rmtree(d)
            print(f"  Deleted {d}")
    # Delete cmd/gui/ entirely (all workflow+tools moved to apex/gui/)
    cmd_gui = APEX / "cmd" / "gui"
    if cmd_gui.exists():
        shutil.rmtree(cmd_gui)
        print(f"  Deleted {cmd_gui}")

    # From apex/lightcurve/ — keep __init__.py and main.py (already updated), delete the rest
    for subdir in ["config", "core", "utils", "analysis", "gui"]:
        d = APEX / "lightcurve" / subdir
        if d.exists():
            shutil.rmtree(d)
            print(f"  Deleted {d}")

    print("Deletion complete.")


# ═══════════════════════════════════════════════════════════════════════════════
# 14. FINAL REPORT
# ═══════════════════════════════════════════════════════════════════════════════

print(f"\n{'='*60}")
print("MIGRATION SUMMARY")
print(f"{'='*60}")
print(f"Files created/updated : {len(created)}")
print(f"Syntax errors         : {len(syntax_errors)}")
print(f"Other errors          : {len(errors)}")

if errors:
    print("\nERRORS:")
    for e in errors:
        print(f"  {e}")

if syntax_errors:
    print("\nSYNTAX ERRORS (fix these before re-running):")
    for e in syntax_errors:
        print(f"  {e}")

print("\nDIVERGED TOOL FILES (LC version used for tools/ — CMD versions differ):")
diverged = [
    ("extinction_fit_window", "lc/gui/workflow/ vs cmd/gui/tools/ — LC has more features (LC USED)"),
    ("iraf_photometry_window", "lc/gui/workflow/ vs cmd/gui/tools/ — LC has more subprocess handling (LC USED)"),
    ("iraf_comparison_window", "lc/gui/workflow/ vs cmd/gui/tools/ — CMD has extra TOML/photutils code (LC USED)"),
    ("qa_report_window", "lc/gui/workflow/ vs cmd/gui/tools/ — LC has photometry_loader (LC USED)"),
]
for name, note in diverged:
    print(f"  {name}: {note}")

print("\nORPHANED FILES (not in target spec, not copied):")
orphaned = [
    "apex/cmd/gui/tools/aperture_overlay_panel.py — no importers found, dropped",
    "apex/cmd/gui/tools/aperture_photometry_worker.py — no importers found, dropped",
    "apex/lightcurve/gui/main_window_lightcurve.py — not referenced, dropped",
    "apex/lightcurve/gui/light_curve/tabs.py — not referenced, dropped",
    "apex/cmd/gui/workflow/step8_ref_build.py — shim to common/shared, superseded",
    "apex/cmd/gui/workflow/step9_star_id_matching.py — shim to common/shared, superseded",
    "apex/cmd/gui/workflow/step2_crop_selector.py — duplicate of shared, dropped",
    "apex/cmd/gui/workflow/step3_sky_preview.py — duplicate of shared, dropped",
    "apex/cmd/gui/workflow/step4_source_detection.py — duplicate of shared, dropped",
    "apex/cmd/gui/workflow/step7_wcs_plate_solving.py — duplicate of shared, dropped",
    "apex/lightcurve/gui/workflow/step2_crop_selector.py — duplicate of shared, dropped",
    "apex/lightcurve/gui/workflow/step3_sky_preview.py — duplicate of shared, dropped",
    "apex/lightcurve/gui/workflow/step4_source_detection.py — duplicate of shared, dropped",
    "apex/lightcurve/gui/workflow/step5_aperture_photometry.py — superseded by cmd version, dropped",
    "apex/lightcurve/gui/workflow/step6_wcs_plate_solving.py — duplicate of shared, dropped",
    "apex/lightcurve/gui/workflow/step7_ref_build.py — duplicate of shared, dropped",
    "apex/lightcurve/gui/workflow/step8_star_id_matching.py — duplicate of shared, dropped",
    "apex/lightcurve/gui/workflow/aperture_photometry_worker.py — superseded by shared worker, dropped",
    "apex/lightcurve/gui/workflow/step_window_base.py — superseded by common version, dropped",
    "apex/cmd/gui/workflow/step_window_base.py — superseded by common version, dropped",
]
for o in orphaned:
    print(f"  {o}")

print("\nPRE-EXISTING IMPORT ISSUES (preserved as-is, not introduced by migration):")
preexisting = [
    "step10_dir — used in lc/utils/step_paths but not defined there (was pre-existing)",
    "step5_dir, step6_dir, step7_dir, step8_dir — used in cmd/gui/tools but not defined in cmd step_paths (pre-existing)",
    "step11_dir, step12_dir, step13_dir — used in cmd workflow files but not defined (pre-existing)",
]
for p in preexisting:
    print(f"  {p}")

print(f"\nDone.")
