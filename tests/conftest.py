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
# test modules import cleanly headless.
#
# Two left this list on 2026-08-16 when Steps 8 and 10 moved their calculation
# to apex.analysis: those tests exercise the photometry itself, which is the
# part that should never have needed a widget toolkit.
_GUI_DEPENDENT_TESTS = [
    "test_iso_cache_worker.py",
    "test_isochrone_fitter_v2.py",
    "test_lc_night_classification.py",
    "test_step6_union_master.py",
    "test_variable_star_phase_plot.py",
]

collect_ignore: list = [] if _HAS_PYQT5 else list(_GUI_DEPENDENT_TESTS)
