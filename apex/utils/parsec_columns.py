"""Read photometric band columns from APEX-compatible isochrone headers."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path


_PARSEC_MAG_TOKEN_TO_BAND = {
    # SDSS
    "umag": "u",
    "gmag": "g",
    "rmag": "r",
    "imag": "i",
    "zmag": "z",
    # Johnson-Cousins / Bessell
    "UXmag": "U",
    "Bmag": "B",
    "Vmag": "V",
    "Rmag": "R",
    "Imag": "I",
}


@dataclass(frozen=True)
class ParsecHeaderInfo:
    photometric_system: str
    band_columns: dict[str, int]
    column_names: tuple[str, ...]

    @property
    def available_bands(self) -> tuple[str, ...]:
        return tuple(self.band_columns)


def read_parsec_header(path: str | Path, max_lines: int = 512) -> ParsecHeaderInfo:
    """Parse column names and the photometric system from an isochrone file.

    The original implementation targeted PARSEC CMD output. It also accepts
    APEX-normalized tables from other model families, including BaSTI, when
    they expose the same metadata names and ``*mag`` band tokens.
    """
    source = Path(path)
    photometric_system = ""
    column_names: tuple[str, ...] = ()

    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle):
            if line_number >= max_lines:
                break
            stripped = line.strip()
            if not stripped:
                continue
            if not stripped.startswith("#"):
                break

            comment = stripped.lstrip("#").strip()
            if "Photometric system:" in comment:
                value = comment.split("Photometric system:", 1)[1].strip()
                photometric_system = re.sub(r"<[^>]+>", "", value).strip()

            tokens = tuple(comment.split())
            has_grid_metadata = {"Zini", "MH", "logAge", "Mini"}.issubset(tokens)
            has_supported_band = any(
                token in tokens for token in _PARSEC_MAG_TOKEN_TO_BAND
            )
            if has_grid_metadata and has_supported_band:
                column_names = tokens

    band_columns = {
        band: column_names.index(token)
        for token, band in _PARSEC_MAG_TOKEN_TO_BAND.items()
        if token in column_names
    }

    if not photometric_system:
        bands = set(band_columns)
        if bands & {"B", "V", "R", "I"}:
            photometric_system = "Johnson-Cousins / Bessell"
        elif bands & {"u", "g", "r", "i", "z"}:
            photometric_system = "SDSS ugriz"
        else:
            photometric_system = "Unknown"

    return ParsecHeaderInfo(
        photometric_system=photometric_system,
        band_columns=band_columns,
        column_names=column_names,
    )
