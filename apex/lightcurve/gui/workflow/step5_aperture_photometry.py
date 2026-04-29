"""
LC Step 5 — Aperture Photometry.

Re-uses the CMD implementation: shared ApertureWorker (growth curve + apcorr)
combined with Step5PhotWorker (per-star photometry at detection positions).

Output written to step5_aperture/:
    aperture_by_frame.csv, apcorr_summary.csv, growth_curve.csv
    photometry_{fname}.tsv, photometry_index.csv
"""
from apex.cmd.gui.workflow.step5_aperture_photometry import AperturePhotometryWindow  # noqa: F401

__all__ = ["AperturePhotometryWindow"]
