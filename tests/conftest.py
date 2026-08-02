"""Pytest configuration.

Lets the test suite run in a HEADLESS environment (no PyQt5 installed), which
is how the cross-platform `test` CI job validates the headless core install.
When PyQt5 is unavailable, the handful of tests that import real GUI worker /
window classes at module load are skipped from collection; everything else
(the Qt-free analysis, pipeline, utils, and config tests) still runs.

With PyQt5 present (local dev, the `test-gui` CI job), nothing is ignored and
the full suite runs.
"""

from __future__ import annotations

import importlib.util

_HAS_PYQT5 = importlib.util.find_spec("PyQt5") is not None

# Test modules whose MODULE-LEVEL imports pull in real GUI classes (which import
# PyQt5). Verified individually: these fail to import without PyQt5; all other
# test modules — including those importing the Qt-free apex.gui.tools.registry
# and apex.gui.workflow.step1_header_target helpers — import cleanly headless.
_GUI_DEPENDENT_TESTS = [
    "test_iso_cache_worker.py",
    "test_isochrone_fitter_v2.py",
    "test_lc_night_classification.py",
    "test_step10_nonlinearity.py",
    "test_step6_union_master.py",
    "test_step8_psf_detection_loader.py",
    "test_variable_star_phase_plot.py",
]

collect_ignore: list = [] if _HAS_PYQT5 else list(_GUI_DEPENDENT_TESTS)
