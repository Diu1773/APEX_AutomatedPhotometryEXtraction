"""Re-anchor Step-10 standard magnitudes on an external standard-star catalog.

Why this exists (validated on real data, 2026-08-04 — see
``validation/psf_crossinstrument/REPORT_UB_DEGENERACY.md``): Step 10 derives
zero-points against *Gaia-transformed* reference magnitudes.  For bands the
Gaia BP/RP spectra barely cover (Johnson U: the ``approx`` relation carries
sigma ~ 0.20 mag) the reference itself is systematically wrong — on M67 the
U zero-point was off by −0.13 mag and B by +0.05, which drove the isochrone
MCMC to a confident, wrong [M/H] (−0.83 vs literature ~0.0).  Re-anchoring
every band on one homogeneous standard-star system (e.g. Montgomery+1993,
tied to Landolt) restored all four cluster parameters.

This module is the Qt-free core: fetch a VizieR standard catalog (cached),
cross-match against the Step-10 wide CMD table by sky position, measure the
per-band median offset ``standard − mag_std``, and shift ``mag_std_<band>``.
Instrumental photometry and ``mag_cal_*`` are never touched — only the
zero-point anchor moves.

The offsets are only trustworthy when the external catalog is *independent*
of the default reference (VizieR standard photometry is; another
Gaia-derived catalog would not be).
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Optional

import numpy as np
import pandas as pd

MAD_TO_SIGMA = 1.4826

#: |offset| above this is almost certainly a wrong catalog / wrong field,
#: not a zero-point error — refuse to apply and surface a warning instead.
MAX_SANE_OFFSET_MAG = 1.0

#: Johnson magnitudes reconstructable from Vmag + colour columns, in the
#: order the chain must be evaluated. Direct "<band>mag" columns win.
_COLOR_CHAIN: dict[str, tuple[tuple[str, str, float], ...]] = {
    # band: ((base_band, colour_column, sign), ...) applied left to right
    "V": (),
    "B": (("V", "B-V", +1.0),),
    "U": (("V", "B-V", +1.0), ("B", "U-B", +1.0)),
    "R": (("V", "V-R", -1.0),),
    "I": (("V", "V-I", -1.0),),
}


@dataclass
class AnchorBandResult:
    band: str
    offset: float          # median(standard − mag_std); add to mag_std
    sigma: float           # robust (MAD-based) scatter of the residuals
    n: int                 # matched stars used
    color_slope: float     # d(offset)/d(colour) diagnostic (0 if no colour)
    applied: bool


@dataclass
class AnchorResult:
    catalog: str
    n_catalog: int
    n_matched: int
    sep_median_arcsec: float
    bands: list[AnchorBandResult] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    residuals: Optional[pd.DataFrame] = None

    @property
    def applied_bands(self) -> list[str]:
        return [b.band for b in self.bands if b.applied]


def _sanitize_catalog_id(catalog_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", catalog_id.strip())


def fetch_standard_catalog(
    catalog_id: str,
    cache_dir: Path,
    table_index: int = 0,
    log: Optional[Callable[[str], None]] = None,
) -> pd.DataFrame:
    """Fetch a VizieR catalog (all columns, all rows); cache as CSV.

    The cache makes re-runs offline and reproducible: delete the file under
    ``<cache_dir>`` to force a re-download.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f"{_sanitize_catalog_id(catalog_id)}_t{table_index}.csv"
    if cache.exists():
        if log:
            log(f"[ANCHOR] catalog cache hit: {cache.name}")
        return pd.read_csv(cache)

    from astroquery.vizier import Vizier

    v = Vizier(columns=["**"], row_limit=-1)
    tables = v.get_catalogs(catalog_id)
    if not tables:
        raise ValueError(f"VizieR returned no tables for {catalog_id!r}")
    df = tables[table_index].to_pandas()
    df.to_csv(cache, index=False, encoding="utf-8")
    if log:
        log(f"[ANCHOR] fetched {catalog_id} table#{table_index}: "
            f"{len(df)} rows -> cached {cache.name}")
    return df


def standard_positions_deg(std: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Extract ICRS degrees from a VizieR table (numeric or sexagesimal)."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    for ra_col, de_col in (("_RA.icrs", "_DE.icrs"), ("RAJ2000", "DEJ2000"),
                           ("RA_ICRS", "DE_ICRS")):
        if ra_col in std.columns and de_col in std.columns:
            ra_num = pd.to_numeric(std[ra_col], errors="coerce")
            de_num = pd.to_numeric(std[de_col], errors="coerce")
            if ra_num.notna().all() and de_num.notna().all():
                return ra_num.to_numpy(float), de_num.to_numpy(float)
            coords = SkyCoord(std[ra_col].astype(str).values,
                              std[de_col].astype(str).values,
                              unit=(u.hourangle, u.deg))
            return coords.ra.deg, coords.dec.deg
    raise ValueError("no usable position columns in the standard catalog "
                     "(_RA.icrs/_DE.icrs, RAJ2000/DEJ2000, RA_ICRS/DE_ICRS)")


def resolve_band_magnitude(std: pd.DataFrame, band: str) -> Optional[pd.Series]:
    """Standard magnitude for ``band``: direct ``<band>mag`` column, else the
    Johnson Vmag + colour chain (U = V + (B−V) + (U−B), R = V − (V−R), …)."""
    direct = f"{band}mag"
    if direct in std.columns:
        return pd.to_numeric(std[direct], errors="coerce")
    chain = _COLOR_CHAIN.get(band)
    if chain is None or "Vmag" not in std.columns:
        return None
    mag = pd.to_numeric(std["Vmag"], errors="coerce")
    for _base, color_col, sign in chain:
        if color_col not in std.columns:
            return None
        mag = mag + sign * pd.to_numeric(std[color_col], errors="coerce")
    return mag


def compute_anchor(
    wide: pd.DataFrame,
    std: pd.DataFrame,
    catalog_id: str,
    bands: Optional[Iterable[str]] = None,
    match_radius_arcsec: float = 1.5,
    min_stars: int = 20,
    log: Optional[Callable[[str], None]] = None,
) -> AnchorResult:
    """Cross-match and measure per-band offsets. Does not modify ``wide``."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    _log = log or (lambda _msg: None)
    result = AnchorResult(catalog=catalog_id, n_catalog=len(std),
                          n_matched=0, sep_median_arcsec=float("nan"))

    if bands is None:
        bands = [c[len("mag_std_"):] for c in wide.columns
                 if c.startswith("mag_std_") and "err" not in c]
    bands = list(bands)

    pos_ok = (pd.to_numeric(wide.get("ra_deg"), errors="coerce").notna()
              & pd.to_numeric(wide.get("dec_deg"), errors="coerce").notna())
    w = wide[pos_ok]
    if w.empty:
        result.warnings.append("wide table has no rows with ra_deg/dec_deg")
        return result

    try:
        std_ra, std_de = standard_positions_deg(std)
    except ValueError as exc:
        result.warnings.append(str(exc))
        return result

    c_apex = SkyCoord(pd.to_numeric(w["ra_deg"]).values * u.deg,
                      pd.to_numeric(w["dec_deg"]).values * u.deg)
    c_std = SkyCoord(std_ra * u.deg, std_de * u.deg)
    idx, sep, _ = c_apex.match_to_catalog_sky(c_std)
    m = sep.arcsec < float(match_radius_arcsec)
    result.n_matched = int(m.sum())
    if result.n_matched == 0:
        result.warnings.append(
            f"no matches within {match_radius_arcsec}\" — wrong field or catalog?")
        return result
    result.sep_median_arcsec = float(np.median(sep.arcsec[m]))

    ww = w[m].reset_index(drop=True)
    ss = std.iloc[idx[m]].reset_index(drop=True)

    # colour axis for the slope diagnostic (first two requested bands)
    color_vals: Optional[np.ndarray] = None
    if len(bands) >= 2:
        c1 = resolve_band_magnitude(ss, bands[0])
        c2 = resolve_band_magnitude(ss, bands[1])
        if c1 is not None and c2 is not None:
            color_vals = (c1 - c2).to_numpy(float)

    res_rows: dict[str, np.ndarray] = {"sep_arcsec": sep.arcsec[m]}
    if "ID" in ww.columns:
        res_rows["ID"] = ww["ID"].to_numpy()

    for band in bands:
        col = f"mag_std_{band}"
        if col not in ww.columns:
            continue
        std_mag = resolve_band_magnitude(ss, band)
        if std_mag is None:
            result.warnings.append(
                f"[{band}] standard catalog has no {band}mag nor a Vmag+colour chain")
            continue
        d = std_mag.to_numpy(float) - pd.to_numeric(ww[col], errors="coerce").to_numpy(float)
        ok = np.isfinite(d)
        n = int(ok.sum())
        res_rows[f"d{band}"] = d
        if n < int(min_stars):
            result.bands.append(AnchorBandResult(band, float("nan"), float("nan"),
                                                 n, 0.0, applied=False))
            result.warnings.append(f"[{band}] only {n} matches (< {min_stars}) — not applied")
            continue
        med = float(np.median(d[ok]))
        sig = float(MAD_TO_SIGMA * np.median(np.abs(d[ok] - med)))
        slope = 0.0
        if color_vals is not None:
            cok = ok & np.isfinite(color_vals)
            if cok.sum() >= 10:
                slope = float(np.polyfit(color_vals[cok], d[cok], 1)[0])
        applied = abs(med) <= MAX_SANE_OFFSET_MAG
        if not applied:
            result.warnings.append(
                f"[{band}] offset {med:+.3f} exceeds {MAX_SANE_OFFSET_MAG} mag — "
                f"wrong catalog/field suspected, not applied")
        result.bands.append(AnchorBandResult(band, med, sig, n, slope, applied))
        _log(f"[ANCHOR][{band}] standard - APEX = {med:+.4f} "
             f"(robust sigma {sig:.4f}, N={n}, colour slope {slope:+.3f})")

    result.residuals = pd.DataFrame(res_rows)
    return result


def apply_anchor(wide: pd.DataFrame, result: AnchorResult) -> pd.DataFrame:
    """Return a copy of ``wide`` with ``mag_std_<band>`` shifted by each
    applied offset. ``mag_inst_*`` / ``mag_cal_*`` are untouched."""
    out = wide.copy()
    for band in result.bands:
        if band.applied and np.isfinite(band.offset):
            out[f"mag_std_{band.band}"] = (
                pd.to_numeric(out[f"mag_std_{band.band}"], errors="coerce")
                + band.offset
            )
    return out


@dataclass
class CatalogCandidate:
    catalog_id: str
    description: str
    bands: list[str]          # requested bands this catalog can supply
    n_sample: int             # rows seen in the probe query
    in_field: bool            # sample positions fall within the search radius


def _probe_candidate(
    catalog_id: str,
    ra_deg: float,
    dec_deg: float,
    radius_deg: float,
    bands: list[str],
    fetch_sample,
) -> Optional[CatalogCandidate]:
    """Check one candidate: does it have usable positions + >=2 of the
    requested bands, with rows inside the field?  ``fetch_sample`` returns a
    small DataFrame for the catalog (injected so tests need no network)."""
    try:
        sample = fetch_sample(catalog_id)
    except Exception:
        return None
    if sample is None or len(sample) == 0:
        return None
    try:
        ra, de = standard_positions_deg(sample)
    except Exception:
        return None
    d_ra = (ra - ra_deg) * np.cos(np.radians(dec_deg))
    sep = np.hypot(d_ra, de - dec_deg)
    in_field = bool(np.nanmin(sep) < radius_deg)
    have = [b for b in bands if resolve_band_magnitude(sample, b) is not None
            and resolve_band_magnitude(sample, b).notna().any()]
    if len(have) < 2:
        return None
    return CatalogCandidate(catalog_id, "", have, len(sample), in_field)


def discover_standard_catalogs(
    target_name: str,
    ra_deg: float,
    dec_deg: float,
    bands: Iterable[str],
    radius_deg: float = 0.5,
    max_candidates: int = 10,
    log: Optional[Callable[[str], None]] = None,
) -> list[CatalogCandidate]:
    """Find VizieR photometric catalogs that could anchor this field.

    Strategy (what a human does, automated): keyword-search VizieR for the
    target ('<name> photometry' / '<name> CCD photometry' / '<name> UBV'),
    then probe each hit — usable positions, >=2 requested bands, sample rows
    inside ``radius_deg`` of the field centre.  Candidates that pass are
    ranked by (in_field, #bands, sample size).  Network-bound; call it from
    a worker, never the GUI thread.
    """
    from astroquery.vizier import Vizier

    _log = log or (lambda _msg: None)
    bands = list(bands)
    name = str(target_name or "").strip()
    # VizieR keyword search treats "NGC 6811" and "NGC6811" differently —
    # try both spellings.
    variants = {name, name.replace(" ", "")} - {""}
    queries = [f"{v} {suffix}" for v in variants
               for suffix in ("photometry", "CCD photometry", "UBV")]
    seen: dict[str, str] = {}
    for query in queries:
        try:
            cats = Vizier.find_catalogs(query) or {}
        except Exception as exc:
            _log(f"[ANCHOR][discover] query {query!r} failed: {exc}")
            continue
        for cid, meta in cats.items():
            if cid not in seen:
                seen[cid] = str(getattr(meta, "description", "") or "")
    _log(f"[ANCHOR][discover] {len(seen)} keyword candidates for {name!r}")

    def fetch_sample(cid: str) -> Optional[pd.DataFrame]:
        v = Vizier(columns=["**"], row_limit=30)
        tables = v.get_catalogs(cid)
        return tables[0].to_pandas() if tables else None

    # Pre-rank by description before probing: keyword hits are NOT
    # relevance-ordered (M67 returned 301 candidates), so probing the first
    # N dict entries can miss the right catalog entirely (false negative).
    # Photometry-looking descriptions mentioning the target go first.
    def _desc_score(item: tuple[str, str]) -> int:
        desc = item[1].lower()
        score = 0
        for kw in ("photometr", "ubv", "ugriz", "ubvri", "standard", "ccd"):
            if kw in desc:
                score += 2
        for v in variants:
            if v.lower() in desc:
                score += 3
        return score

    ranked = sorted(seen.items(), key=_desc_score, reverse=True)

    out: list[CatalogCandidate] = []
    for cid, desc in ranked[:max_candidates]:
        cand = _probe_candidate(cid, ra_deg, dec_deg, radius_deg, bands, fetch_sample)
        if cand is None:
            continue
        cand.description = desc
        out.append(cand)
        _log(f"[ANCHOR][discover] {cid}: bands={cand.bands} "
             f"in_field={cand.in_field} | {desc[:60]}")
    out.sort(key=lambda c: (c.in_field, len(c.bands), c.n_sample), reverse=True)
    return out


def anchor_qc_frame(result: AnchorResult) -> pd.DataFrame:
    """Per-band QC table for ``standard_anchor_offsets.csv``."""
    return pd.DataFrame([
        {
            "band": b.band,
            "offset_mag": b.offset,
            "robust_sigma": b.sigma,
            "n_matched": b.n,
            "color_slope": b.color_slope,
            "applied": b.applied,
            "catalog": result.catalog,
            "sep_median_arcsec": result.sep_median_arcsec,
        }
        for b in result.bands
    ])
