"""
I/O utility functions
APEX I/O helpers.
"""

from __future__ import annotations
from pathlib import Path
import json
import time
from collections import deque
from typing import Union, Optional
from decimal import Decimal, InvalidOperation

import pandas as pd


def frame_bytes_from_header(path, dtype_bytes: int = 4) -> Optional[int]:
    """In-memory size of one frame, from ``NAXIS1``/``NAXIS2`` alone.

    The worker pool has to be sized against free RAM *before* any frame is
    loaded, so the size has to come from the header: a few kB to read rather
    than the 61 MB the frame itself costs.  ``dtype_bytes`` defaults to 4
    because the pipeline loads frames as float32.

    Returns None when the header cannot be read — callers treat that as
    "unknown", not as zero.
    """
    from astropy.io import fits

    try:
        with fits.open(path, memmap=False) as hdul:
            for hdu in hdul:
                nx = hdu.header.get("NAXIS1")
                ny = hdu.header.get("NAXIS2")
                if nx and ny:
                    return int(nx) * int(ny) * int(dtype_bytes)
    except Exception:
        return None
    return None


def science_hdu_index(hdul) -> int:
    """Index of the HDU that carries the image and its keywords.

    A ``.fits.fz`` file — how public archives serve frames — keeps an **empty**
    ``PrimaryHDU`` at index 0 and the image in a ``CompImageHDU`` at index 1.
    A multi-extension file does the same. Reading index 0 there returns a
    header with no ``OBSTYPE``, ``EXPTIME``, ``FILTER``, ``DATE-OBS`` or
    ``GAIN`` at all, so every keyword lookup silently returns the default.
    Measured on LCO MuSCAT3 frames, 2026-09-07 (``Main/FAILURES.md`` F-302).

    The choice is made from ``NAXIS`` alone, so no pixel data is read.
    """
    for i, hdu in enumerate(hdul):
        try:
            if int(hdu.header.get("NAXIS") or 0) > 0:
                return i
        except (TypeError, ValueError):
            continue
    return 0


def read_fits_header(path, ext=None):
    """Header of the image-bearing HDU. Use instead of ``fits.getheader(path)``.

    ``fits.getheader`` defaults to HDU 0, which is empty in a compressed file.
    Pass ``ext`` to force a particular index.
    """
    from astropy.io import fits

    with fits.open(path, memmap=False) as hdul:
        i = science_hdu_index(hdul) if ext is None else int(ext)
        return hdul[i].header.copy()


def read_fits_image(path, dtype=None, ext=None):
    """``(data, header)`` of the image-bearing HDU.

    Use instead of ``hdul[0].data`` / ``hdul[0].header`` so a compressed or
    multi-extension frame reads the same as a plain one.
    """
    from astropy.io import fits

    with fits.open(path, memmap=False) as hdul:
        i = science_hdu_index(hdul) if ext is None else int(ext)
        data = hdul[i].data
        if data is None:
            raise ValueError(f"no image data in {path}")
        if dtype is not None:
            data = data.astype(dtype, copy=False)
        else:
            data = data.copy()
        return data, hdul[i].header.copy()


def load_toml(path: Union[str, Path]) -> dict:
    """Parse a TOML file, tolerating a UTF-8 BOM.

    PowerShell 5.1's ``-Encoding utf8`` prefixes files with a BOM, which
    ``tomllib.load`` rejects as "Invalid statement (at line 1, column 1)".
    Decoding with ``utf-8-sig`` accepts both BOM'd and clean files, so an
    externally rewritten parameters.toml can't brick every entry point.
    """
    try:
        import tomllib
    except ImportError:  # pragma: no cover - py3.10
        import tomli as tomllib  # type: ignore
    return tomllib.loads(Path(path).read_bytes().decode("utf-8-sig"))


def _parse_int64_col(series: pd.Series) -> pd.array:
    """Convert a string/object series of source_ids to pandas Int64.

    - Exact integer strings: "2823345641878527872" → preserved with full precision
    - Float strings (old format): "2823345641878528000.0" → int (already-rounded)
    - Blank / nan / NA → pd.NA

    Fast paths (avoid the per-element Decimal loop, which cost ~4 s on a
    40 k-row Gaia ECSV cache):
      1. Already an integer dtype  → just view as Int64.
      2. Clean integer strings     → vectorised ``pd.to_numeric``.
    The slow element-wise Decimal parse only runs for genuinely mixed /
    float-formatted data.
    """
    # Fast path 1: column is already integer (typical for an ECSV cache
    # where astropy preserved source_id as int64) — no parsing needed.
    try:
        if pd.api.types.is_integer_dtype(series.dtype):
            return series.astype("Int64").array
    except Exception:
        pass

    # Fast path 2: object/string column that is purely integer strings.
    # 19-digit Gaia ids fit in int64 (max 9.2e18), so to_numeric is exact
    # as long as no value needs float (decimal point / exponent / NaN).
    try:
        s = series.astype("string").str.strip()
        s = s.replace({"": pd.NA, "nan": pd.NA, "NaN": pd.NA,
                       "<NA>": pd.NA, "None": pd.NA})
        non_null = s.dropna()
        if len(non_null) and not non_null.str.contains(r"[.eE]", regex=True).any():
            return pd.array(pd.to_numeric(s, errors="raise"), dtype="Int64")
    except Exception:
        pass

    def _parse(s):
        if pd.isna(s):
            return pd.NA
        sval = str(s).strip()
        if sval in ("", "nan", "NaN", "<NA>", "None"):
            return pd.NA
        try:
            # Keep exact precision for integer-like decimal/exponent strings.
            if "." in sval or "e" in sval.lower():
                d = Decimal(sval)
                if not d.is_finite():
                    return pd.NA
                if d != d.to_integral_value():
                    return pd.NA
                return int(d)
            return int(sval)
        except (ValueError, OverflowError, InvalidOperation):
            return pd.NA
    return pd.array([_parse(v) for v in series], dtype="Int64")


def coerce_int64_source_id(series: pd.Series) -> pd.Series:
    """Coerce arbitrary source_id values to nullable Int64 without float round-trip."""
    if series is None:
        return pd.Series(pd.array([], dtype="Int64"))
    parsed = _parse_int64_col(series)
    return pd.Series(parsed, index=series.index, dtype="Int64")


# Alias used by workflow code that expects parse_int64_series.
def parse_int64_series(series: pd.Series) -> pd.Series:
    """Convert a series to pandas nullable Int64 without float precision loss.

    Alias for coerce_int64_source_id for cross-project compatibility.
    """
    return coerce_int64_source_id(series)


def parse_int64_scalar(value):
    """Convert one Gaia source_id-like value to int, or pd.NA when invalid."""
    parsed = coerce_int64_source_id(pd.Series([value]))
    out = parsed.iloc[0] if len(parsed) else pd.NA
    return pd.NA if pd.isna(out) else int(out)


#: Every column that can hold a 19-digit Gaia DR3 identifier. float64 runs out
#: of mantissa at ~9e15, so any of these landing in a float column loses its
#: last three digits — that is how ``8.14859719594028e+17`` was written into
#: master_catalog_*.tsv instead of ``814859719594028032``.
ID_LIKE_COLUMNS = ("source_id", "gaia_source_id", "gaia_id", "star_id")


def normalize_id_columns(df: pd.DataFrame, columns=ID_LIKE_COLUMNS) -> pd.DataFrame:
    """Cast identifier columns back to nullable Int64, in place.

    Call it right before writing a table so ``to_csv`` prints whole digits, and
    after building one from row dicts (a dict column holding ``np.nan`` becomes
    float64 at construction time and the ids are already rounded by then; use
    ``pd.NA`` there instead).

    This stops further loss — it cannot recover digits a float column has
    already dropped.
    """
    if df is None or len(df) == 0:
        return df
    for col in columns:
        if col in df.columns:
            df[col] = coerce_int64_source_id(df[col]).astype("Int64")
    return df


def read_csv_int64_source_id(path: Union[str, Path], sep: str = ",", **kwargs) -> pd.DataFrame:
    """Read a CSV/TSV file preserving 19-digit Gaia source_id precision.

    pandas default read_csv promotes a column with mixed integer/NaN to float64,
    silently rounding the last 3-4 digits of 19-digit Gaia source_ids.
    Every identifier column is read as string and converted to Int64.
    """
    dtype = dict(kwargs.pop("dtype", None) or {})
    for col in ID_LIKE_COLUMNS:
        dtype.setdefault(col, str)
    try:
        df = pd.read_csv(path, sep=sep, dtype=dtype, **kwargs)
    except ValueError:
        # A dtype for a column the file does not have is fine on modern pandas;
        # fall back for any reader that rejects the mapping outright.
        df = pd.read_csv(path, sep=sep, dtype={"source_id": str}, **kwargs)
    return normalize_id_columns(df)


def _fast_read_ecsv(path: Union[str, Path]) -> Optional[pd.DataFrame]:
    """Read an ECSV body with pandas, bypassing astropy.

    An ECSV file is a whitespace-delimited table whose schema lives in
    ``#``-comment lines; the first non-comment line is the column header.
    For the all-numeric Gaia catalogue this lets pandas parse it in
    ~0.3 s vs ~3.7 s for ``Table.read`` + ``to_pandas`` on a 40 k-row
    cache.  source_id columns are read as strings to preserve 19-digit
    precision.  Returns None if the structure looks unexpected so the
    caller can fall back to the robust astropy reader.
    """
    try:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            comment="#",
            dtype={"source_id": str, "gaia_source_id": str},
        )
        if df.empty or df.shape[1] < 2:
            return None
        # ECSV header line itself starts with no '#', so the column names
        # must look like real identifiers, not numeric data.
        if any(str(c).replace(".", "", 1).lstrip("-").isdigit() for c in df.columns):
            return None
        return df
    except Exception:
        return None


def read_ecsv_int64_source_id(path: Union[str, Path], **kwargs) -> pd.DataFrame:
    """Read an Astropy ECSV table preserving Gaia source_id precision."""
    if not kwargs:
        df = _fast_read_ecsv(path)
        if df is not None:
            for col in ("source_id", "gaia_source_id"):
                if col in df.columns:
                    df[col] = coerce_int64_source_id(df[col])
            return df

    from astropy.table import Table

    table = Table.read(str(path), format="ascii.ecsv", **kwargs)
    df = table.to_pandas()
    for col in ("source_id", "gaia_source_id"):
        if col in df.columns:
            df[col] = coerce_int64_source_id(df[col])
    return df


def load_night_assignments(result_dir: Path) -> dict[str, int]:
    """Load filename → night_id from step1/night_assignments.json."""
    from .step_paths import step1_dir
    na_path = step1_dir(result_dir) / "night_assignments.json"
    if not na_path.exists():
        na_path = result_dir / "night_assignments.json"
    if not na_path.exists():
        return {}
    try:
        data = json.loads(na_path.read_text(encoding="utf-8"))
        return {k: int(v) for k, v in data.get("night_assignments", {}).items()}
    except Exception:
        return {}


def load_headers_table(result_dir: Path) -> pd.DataFrame:
    """Read headers.csv from step1 dir (fallback: result_dir)."""
    from .step_paths import step1_dir
    headers_path = step1_dir(result_dir) / "headers.csv"
    if not headers_path.exists():
        headers_path = result_dir / "headers.csv"
    if not headers_path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(headers_path)
    except Exception:
        return pd.DataFrame()


def load_file_path_map(result_dir: Path) -> dict[str, str]:
    """Load filename → original FITS path mapping from step1/file_path_map.json."""
    from .step_paths import step1_dir
    path = step1_dir(result_dir) / "file_path_map.json"
    if not path.exists():
        path = result_dir / "file_path_map.json"
    if not path.exists():
        project_state_path = Path(result_dir) / "project_state.json"
        if not project_state_path.exists():
            return {}
        try:
            state = json.loads(project_state_path.read_text(encoding="utf-8"))
            step_data = state.get("step_data", {})
            file_sel = step_data.get("file_selection", {}) if isinstance(step_data, dict) else {}
            data = file_sel.get("file_path_map", {}) if isinstance(file_sel, dict) else {}
        except Exception:
            return {}
    try:
        if path.exists():
            data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}
    out: dict[str, str] = {}
    for key, value in data.items():
        if not key or not value:
            continue
        out[str(key)] = str(value)
    return out


class TailLogger:
    """
    Logger that maintains a tail buffer of recent messages
    Useful for displaying recent activity in GUI
    """

    def __init__(self, log_path: Path, tail: int = 5, enable_console: bool = True):
        """
        Initialize tail logger

        Args:
            log_path: Path to log file
            tail: Number of recent messages to keep in buffer
            enable_console: Whether to print to console
        """
        self.log_path = Path(log_path)
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(self.log_path, "a", encoding="utf-8")
        self.buf = deque(maxlen=max(1, tail))
        self.enable_console = enable_console

        # Try to import IPython clear_output for Jupyter support
        try:
            from IPython.display import clear_output
            self._clear = lambda: clear_output(wait=True)
        except Exception:
            self._clear = lambda: None

    def write(self, msg: str):
        """Write message to log file and buffer"""
        ts = time.strftime("%H:%M:%S")
        line = f"[{ts}] {msg}"
        self.fh.write(line + "\n")
        self.fh.flush()

        if self.enable_console:
            self.buf.append(line)
            self._clear()
            print("\n".join(self.buf))

    def get_recent(self) -> list[str]:
        """Get recent messages from buffer"""
        return list(self.buf)

    def close(self):
        """Close log file"""
        try:
            self.fh.close()
        except:
            pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
