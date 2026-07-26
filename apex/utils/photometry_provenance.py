"""Canonical provenance labels for aperture and PSF magnitude products."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any


APERTURE_SOURCE = "aperture"
PSF_SOURCE = "psf"
MIXED_SOURCE = "mixed"
UNKNOWN_SOURCE = "unknown"


def infer_photometry_source(mag_column: str | None, source: str | None = None) -> str:
    value = str(source or "").strip().lower()
    if value in {APERTURE_SOURCE, PSF_SOURCE, MIXED_SOURCE}:
        return value
    column = str(mag_column or "").strip().lower()
    if column.startswith("mag_psf"):
        return PSF_SOURCE
    if column:
        return APERTURE_SOURCE
    return UNKNOWN_SOURCE


def build_photometry_provenance(
    source: str | None = None,
    mag_column: str | None = None,
    mag_error_column: str | None = None,
) -> dict[str, str]:
    normalized = infer_photometry_source(mag_column, source)
    if not mag_column:
        mag_column = {
            PSF_SOURCE: "mag_psf",
            APERTURE_SOURCE: "mag",
            MIXED_SOURCE: MIXED_SOURCE,
        }.get(normalized, UNKNOWN_SOURCE)
    if not mag_error_column:
        mag_error_column = {
            PSF_SOURCE: "mag_psf_err",
            APERTURE_SOURCE: "mag_err",
            MIXED_SOURCE: MIXED_SOURCE,
        }.get(normalized, UNKNOWN_SOURCE)
    return {
        "source": normalized,
        "mag_column": str(mag_column),
        "mag_error_column": str(mag_error_column),
    }


def format_photometry_provenance(info: dict[str, Any] | None) -> str:
    payload = info or {}
    provenance = build_photometry_provenance(
        payload.get("source"),
        payload.get("mag_column") or payload.get("mag_input_column"),
        payload.get("mag_error_column") or payload.get("mag_error_input_column"),
    )
    labels = {
        APERTURE_SOURCE: "Aperture",
        PSF_SOURCE: "PSF",
        MIXED_SOURCE: "Mixed",
        UNKNOWN_SOURCE: "Unknown",
    }
    return f"Photometry: {labels[provenance['source']]} | MAG: {provenance['mag_column']}"


def collapse_provenance_values(values: Iterable[Any]) -> str:
    unique = {
        str(value).strip()
        for value in values
        if value is not None and str(value).strip() and str(value).strip().lower() != "nan"
    }
    if not unique:
        return ""
    if len(unique) == 1:
        return next(iter(unique))
    return MIXED_SOURCE


def summarize_photometry_table(table: Any) -> dict[str, str]:
    columns = list(getattr(table, "columns", []))

    def _collect(prefix: str) -> str:
        names = [name for name in columns if name == prefix or name.startswith(prefix + "_")]
        values: list[Any] = []
        for name in names:
            try:
                values.extend(table[name].dropna().tolist())
            except Exception:
                continue
        return collapse_provenance_values(values)

    source = _collect("photometry_source")
    mag_column = _collect("mag_input_column")
    error_column = _collect("mag_error_input_column")
    return build_photometry_provenance(source, mag_column, error_column)
