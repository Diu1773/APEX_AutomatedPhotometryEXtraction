"""Config rows to `IsochroneFitConfig` — one translation, two callers.

The fit is decided by things that are not visible in its output: which colours
were fitted, how wide the age window was, and — the one that bites — where the
walls of the [M/H], distance-modulus and reddening boxes were. A posterior
median sitting on a wall is the wall, not a measurement.

This module exists because those settings were being decided in two places. The
headless step read them from the config; the desktop dialog hardcoded its own,
and the two disagreed. The GUI's box was `[M/H] in (-1.0, +0.5)`, which cannot
contain a globular cluster — M13 is near -1.5 — and the service narrows the box
around an [M/H] prior with `max(lo, ...)`/`min(hi, ...)`, so asking for a
metal-poor prior inside that box produced an *inverted* range: lo = -1.000,
hi = -1.360. The grid mask `(mh >= lo) & (mh <= hi)` is then empty and the fit
cannot start. Measured 2026-08-17.

So both callers build the configuration here. The dialog still owns the fields
its widgets own — colours, age window, walkers, the prior checkboxes — and
replaces only those.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

__all__ = [
    "colors_from_params",
    "prior_from_params",
    "build_fit_config",
    "missing_decisive_settings",
    "check_bounds",
]


def colors_from_params(params) -> list[tuple[str, str]]:
    """`isochrone.colors` as pairs, e.g. "B-V,V-R" -> [("B","V"), ("V","R")]."""
    raw = str(getattr(params.P, "iso_colors", "") or "").strip()
    pairs = []
    for token in raw.replace(";", ",").split(","):
        token = token.strip()
        if not token or "-" not in token:
            continue
        left, _, right = token.partition("-")
        if left.strip() and right.strip():
            pairs.append((left.strip(), right.strip()))
    return pairs


def prior_from_params(params, name: str) -> Optional[Tuple[float, float]]:
    """A prior written as "value,sigma"; absent means no prior."""
    raw = str(getattr(params.P, name, "") or "").strip()
    if not raw:
        return None
    parts = [p.strip() for p in raw.replace(";", ",").split(",")]
    if len(parts) != 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None


def build_fit_config(params, iso_file: Path | str | None = None):
    """Turn the config rows into an `IsochroneFitConfig`.

    Kept separate from any step so a test can check the translation without
    paying for an MCMC.
    """
    from apex.analysis.cmd.isochrone_fit_service import IsochroneFitConfig

    P = params.P
    path = iso_file or (str(getattr(P, "iso_file_path", "") or "").strip() or None)
    return IsochroneFitConfig(
        colors=colors_from_params(params),
        mag_band=str(getattr(P, "iso_mag_band", "g") or "g"),
        iso_file=str(path) if path else None,
        age_bounds=(float(P.iso_age_min), float(P.iso_age_max)),
        mh_bounds=(float(P.iso_mh_min), float(P.iso_mh_max)),
        dm_bounds=(float(P.iso_dm_min), float(P.iso_dm_max)),
        ecolor_bounds=(float(P.iso_ecolor_min), float(P.iso_ecolor_max)),
        mh_prior=prior_from_params(params, "iso_mh_prior"),
        ecolor_prior=prior_from_params(params, "iso_ecolor_prior"),
        dm_prior=prior_from_params(params, "iso_dm_prior"),
        parallax_distance_prior=bool(P.iso_parallax_distance_prior),
        parallax_dm_sigma=float(P.iso_parallax_dm_sigma),
        parallax_dm_window=float(P.iso_parallax_dm_window),
        use_membership=bool(P.iso_use_membership),
        data_snr_min=float(P.iso_data_snr_min),
        fit_snr_min=float(P.iso_fit_snr_min),
        max_stars=int(P.iso_max_stars),
        n_walkers=int(P.iso_n_walkers),
        n_burn=int(P.iso_n_burn),
        n_steps=int(P.iso_n_steps),
        f_bin=float(P.iso_f_bin),
        f_field=float(P.iso_f_field),
        err_floor=float(P.iso_err_floor),
        seed=int(P.iso_seed),
    )


def missing_decisive_settings(params) -> list[str]:
    """What must be written down before a batch fit means anything.

    Only two: the colours, because the fit is undefined without them, and the
    isochrone grid, because there is no sensible default file.
    """
    missing = []
    if not colors_from_params(params):
        missing.append('isochrone.colors (예: "B-V,V-R")')
    if not str(getattr(params.P, "iso_file_path", "") or "").strip():
        missing.append("isochrone.file_path (PARSEC 격자 파일)")
    return missing


def check_bounds(config) -> list[str]:
    """Boxes a prior would turn inside out, named before the fit starts.

    The service tightens each box around its prior and never widens it. That is
    the right behaviour, but it means a prior outside the box leaves lo > hi,
    and the only downstream symptom is an empty grid — a failure that reads as
    "the data are bad" rather than "the walls are in the wrong place".
    """
    problems = []
    checks = (
        ("[M/H]", config.mh_bounds, config.mh_prior, 0.04),
        ("E(colour)", config.ecolor_bounds, config.ecolor_prior, 0.008),
        ("(m-M)", config.dm_bounds, config.dm_prior, 0.0),
    )
    for label, box, prior, floor in checks:
        if prior is None or box is None:
            continue
        lo, hi = float(box[0]), float(box[1])
        centre, sigma = float(prior[0]), float(prior[1])
        width = max(2.0 * sigma, floor)
        lo, hi = max(lo, centre - width), min(hi, centre + width)
        if lo >= hi:
            problems.append(
                f"{label}: 사전값 {centre:+.3f}±{sigma:.3f} 가 탐색 상자 "
                f"[{float(box[0]):+.3f}, {float(box[1]):+.3f}] 밖이라 "
                f"범위가 [{lo:+.3f}, {hi:+.3f}] 로 뒤집힌다"
            )
    return problems
