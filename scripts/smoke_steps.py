"""Import-smoke check for APEX workflow steps and tool windows.

This catches migration misses such as renamed modules or missing shared helper
functions before a GUI button closes silently.
"""

from __future__ import annotations

import importlib
from pathlib import Path
import sys


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


MODULES = [
    "apex.gui.main_window",
    "apex.gui.workflow.step2_crop_selector",
    "apex.gui.workflow.step3_sky_preview",
    "apex.gui.workflow.step4_source_detection",
    "apex.gui.workflow.step5_wcs_plate_solving",
    "apex.gui.workflow.step6_ref_build",
    "apex.gui.workflow.step7_forced_aperture_phot",
    "apex.gui.workflow.cmd.step1_file_selection",
    "apex.gui.workflow.cmd.step8_psf_photometry",
    "apex.gui.workflow.cmd.step9_master_id_editor",
    "apex.gui.workflow.cmd.step10_zeropoint_calibration",
    "apex.gui.workflow.cmd.step11_cmd_plot",
    "apex.gui.workflow.cmd.step12_isochrone_model",
    "apex.gui.workflow.lc.step1_file_selection",
    "apex.gui.workflow.lc.step8_target_selection",
    "apex.gui.workflow.lc.step9_lightcurve_builder",
    "apex.gui.workflow.lc.step10_detrend_merge",
    "apex.gui.workflow.lc.step11_period_analysis",
    "apex.gui.tools.gaia_3d_viewer",
    "apex.gui.tools.extinction_fit",
    "apex.gui.tools.iraf_photometry",
    "apex.gui.tools.multi_night_merger",
    "apex.gui.tools.variable_star",
    "apex.gui.tools.transit_tool",
    "apex.gui.tools.eb_tool",
]


def main() -> int:
    failed: list[tuple[str, str, str]] = []
    for module_name in MODULES:
        try:
            importlib.import_module(module_name)
        except Exception as exc:
            failed.append((module_name, type(exc).__name__, str(exc)))
            print(f"FAIL {module_name}: {type(exc).__name__}: {exc}")
        else:
            print(f"OK   {module_name}")

    print(f"TOTAL_FAIL {len(failed)}")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
