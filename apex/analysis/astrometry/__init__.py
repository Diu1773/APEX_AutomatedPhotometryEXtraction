"""APEX internal plate-solving package.

Astrometry.net-style geometric quad matching with optional Gaia-projected
prior (target center + pixel scale from Step 1 + Instrument config).

Entry point: ``apex.analysis.astrometry.solver.solve``.
"""

from .solver import solve, WCSSolveResult

__all__ = ["solve", "WCSSolveResult"]
