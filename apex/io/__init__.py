"""Community-standard photometry output exporters for APEX.

This package converts APEX light-curve tables into the file formats accepted by
the major amateur/pro-am photometry archives:

- :func:`export_aavso_extended` — AAVSO Extended File Format (variable-star
  observations submitted to the AAVSO International Database).
- :func:`export_exoclock` — ExoClock-style transit photometry (BJD_TDB time,
  normalized flux + error).
- :func:`export_exofop` — ExoFOP/TFOP-SG1 seeing-limited photometry table.

LAYER RULE: this package is Qt-free. It must NOT import ``apex.gui`` or
``PyQt5`` (so the exporters run in headless CI and behind ``apex export``).
"""

from __future__ import annotations

from apex.io.exporters import (
    export_aavso_extended,
    export_exoclock,
    export_exofop,
)

__all__ = [
    "export_aavso_extended",
    "export_exoclock",
    "export_exofop",
]
