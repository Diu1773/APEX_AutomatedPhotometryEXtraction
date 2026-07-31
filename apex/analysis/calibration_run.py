"""Run detector calibration over a scanned frame set — Qt-free (Step 0 core).

This is the one implementation of "build the masters, match each light, write
``pp_*.fits`` and ``calibration.json``". The GUI worker, the headless pipeline
step and the reprocess scripts all call :func:`run_calibration`, so what the
user sees in the window and what a batch run produces cannot drift apart.

Progress and log lines go out through plain callables; cancellation is a plain
predicate. Nothing here imports Qt.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
from astropy.io import fits

from apex.analysis import calibration as cal
from apex.analysis import calibration_scan as scan
from apex.analysis.calibration import CalibrationOptions
from apex.analysis.calibration_scan import FrameInfo

ProgressFn = Callable[[int, int, str], None]
LogFn = Callable[[str], None]
StopFn = Callable[[], bool]

ALL_NIGHTS = "__all__"


def _noop(*_args, **_kwargs) -> None:
    return None


def _never() -> bool:
    return False


def accept_dark(match: Optional[scan.DarkMatch], light: FrameInfo,
                opts: CalibrationOptions, log: LogFn, counters: Dict[str, int]):
    """Apply the temperature/exposure policy to a dark match.

    Returns the dark group key to use, or None to calibrate without a dark.
    A mismatch is never silent: it is logged, and ``strict_temp`` turns the
    warning into a refusal.
    """
    if match is None:
        return None
    if not match.within_temp_tol:
        counters["temp_mismatch"] += 1
        verdict = "REFUSED (strict)" if opts.strict_temp else "used anyway"
        log(f"  [warn] {light.name}: nearest dark is ΔT="
            f"{match.delta_temp_c:.2f}°C away (tolerance "
            f"{opts.temp_match_tol_c:g}°C) — {verdict}")
        if opts.strict_temp:
            counters["dark_refused"] += 1
            return None
    # A dark of a different exposure is only defensible if it gets scaled.
    if match.delta_exp_s > 1e-3 and not opts.dark_scale:
        counters["exp_mismatch"] += 1
        log(f"  [warn] {light.name}: dark exposure {match.exp:g}s vs light "
            f"{light.exp:g}s and dark scaling is OFF — dark is not "
            f"exposure-matched")
    return match.key


def run_calibration(frames: List[FrameInfo], night: str, out_dir: Path,
                    opts: CalibrationOptions,
                    progress: Optional[ProgressFn] = None,
                    log: Optional[LogFn] = None,
                    stop: Optional[StopFn] = None) -> Dict:
    """Calibrate every light frame of ``night`` (or :data:`ALL_NIGHTS`).

    Writes ``<out_dir>/masters/``, ``<out_dir>/calibrated/<night>/pp_*.fits``
    and ``<out_dir>/calibration.json``, and returns the summary dict that was
    written.
    """
    progress = progress or _noop
    log = log or _noop
    stop = stop or _never

    out_dir = Path(out_dir)
    masters_dir = out_dir / "masters"
    masters_dir.mkdir(parents=True, exist_ok=True)

    nights = scan.nights(frames) if night == ALL_NIGHTS else [night]
    nights = [n for n in nights if scan.group_for_night(frames, n)["light"]]
    if not nights:
        raise ValueError("No light frames found to calibrate")

    # Master cache keyed by the exact source files, so a master built from a
    # shared bias/dark pool (or a single pre-built master frame) is built ONCE
    # and reused across every night; per-night flats build per night.
    cache: Dict = {}
    written: set = set()
    counters = {"temp_mismatch": 0, "exp_mismatch": 0, "dark_refused": 0}
    summary: Dict = {"nights": {}, "masters": [], "options": opts.to_mapping()}

    def _write_master(fname, arr, prov, label, extra=None):
        if fname in written:
            return
        # Stamp a MASTER header so APEX's own masters are auto-detected and
        # reused verbatim on a later run (parity with AIPPI).
        hdr = fits.Header()
        kind = str(prov.get("type", "")).capitalize()
        if kind:
            hdr["IMAGETYP"] = f"Master {kind}"
        if prov.get("exptime") is not None:
            hdr["EXPTIME"] = float(prov["exptime"])
        hdr["NCOMBINE"] = int(prov.get("n_frames", 1))
        for key, value in (extra or {}).items():
            hdr[key] = value
        fits.writeto(masters_dir / fname, arr.astype(np.float32), header=hdr,
                     overwrite=True)
        written.add(fname)
        summary["masters"].append({"file": fname, **prov})
        tag = " [pre-built master, used directly]" if prov.get("master_input") else ""
        log(f"{label}: {prov['n_frames']} frame(s){tag} → {fname}")

    def _bias(paths):
        if not paths:
            return None
        key = ("bias", frozenset(paths))
        if key not in cache:
            master, prov = cal.build_master_bias(paths, opts)
            cache[key] = master
            _write_master("master_bias.fits", master, prov, "master bias")
        return cache[key]

    def _dark(frames_in, dark_key, master_bias, night_tag):
        paths = [f.path for f in frames_in]
        if not paths:
            return None, 1.0
        key = ("dark", frozenset(paths))
        if key not in cache:
            master, dark_exp, prov = cal.build_master_dark(
                paths, opts, master_bias=master_bias)
            cache[key] = (master, dark_exp)
            exp, bucket = dark_key
            tag = night_tag or "global"   # keep per-night darks distinct on disk
            extra = {"CCD-TEMP": float(bucket)} if bucket is not None else None
            _write_master(f"master_dark_{exp:g}s_{bucket}C_{tag}.fits", master,
                          prov, f"master dark {exp:g}s/{bucket}°C [{tag}]",
                          extra=extra)
        return cache[key]

    def _flat(flat_frames, filt, master_bias):
        if not flat_frames:
            return None
        paths = [f.path for f in flat_frames]
        key = ("flat", frozenset(paths))
        if key not in cache:
            master, prov = cal.build_master_flat(paths, opts, master_bias=master_bias)
            cache[key] = master
            tag = flat_frames[0].night or "global"    # the flats' own night
            _write_master(f"master_flat_{filt}_{tag}.fits", master, prov,
                          f"master flat {filt} [{tag}]", extra={"FILTER": filt})
        return cache[key]

    total = sum(len(scan.group_for_night(frames, n)["light"]) for n in nights)
    done = 0
    no_flat = 0
    log(f"Calibrating {total} light frames across {len(nights)} night(s)…")

    for night_key in nights:
        group = scan.group_for_night(frames, night_key)
        calibrated_dir = out_dir / "calibrated" / (night_key or "undated")
        calibrated_dir.mkdir(parents=True, exist_ok=True)
        # The bias pool is resolved per light (below), not once per night: a
        # pool of mixed geometry cannot be stacked into one master.
        night_recs: List[Dict] = []
        for light in group["light"]:
            if stop():
                log("Stopped by user.")
                break
            progress(done, total, f"[{night_key}] {light.name}")
            # Every calibration frame must share the light's geometry;
            # mismatched binning/size is filtered out at match time rather than
            # raising a broadcast error mid-subtraction.
            master_bias = _bias([f.path for f in group["bias"]
                                 if scan.compatible_geometry(light, f)])
            match = scan.match_dark_detail(
                group["dark"], light.exp, light.temp, opts.temp_match_tol_c,
                light=light, fallback=group.get("dark_library"))
            dark_key = accept_dark(match, light, opts, log, counters)
            master_dark, dark_exp = (
                _dark(match.frames, dark_key, master_bias, match.night)
                if dark_key else (None, 1.0))
            flat_filt = scan.match_flat(group["flat"], light.filt, light=light)
            master_flat = _flat(group["flat"][flat_filt], flat_filt, master_bias) \
                if flat_filt else None
            if master_flat is None:
                no_flat += 1
                log(f"  [warn] {light.name}: no flat for '{light.filt}' — flat skipped")

            data, header, qc = cal.calibrate_light_file(
                light.path, opts, master_bias=master_bias,
                master_dark=master_dark, dark_exp=dark_exp,
                master_flat=master_flat)
            fits.writeto(calibrated_dir / f"pp_{light.name}", data,
                         header=header, overwrite=True)
            done += 1
            night_recs.append({
                "input": light.name, "output": f"pp_{light.name}",
                "filter": light.filt,
                "dark": f"{dark_key[0]:g}s/{dark_key[1]}C" if dark_key else None,
                "flat": flat_filt,
                # Provenance for the dark match: how far the chosen dark sat
                # from this light in temperature and exposure, and whether it
                # came from this night or from a shared dark library.
                "dark_delta_temp_c": (match.delta_temp_c if match is not None else None),
                "dark_delta_exp_s": (match.delta_exp_s if match is not None else None),
                "dark_night": (match.night if match is not None else None),
                "dark_source": (match.source if match is not None else None),
                "median": qc.get("median"), "neg_pct": qc.get("neg_pct"),
            })
        summary["nights"][night_key] = {
            "calibrated_dir": str(calibrated_dir),
            "n_calibrated": len(night_recs), "frames": night_recs,
        }
    progress(total, total, "Done")

    summary["n_calibrated"] = done
    summary["n_missing_flat"] = no_flat
    summary["n_temp_mismatch"] = counters["temp_mismatch"]
    summary["n_exp_mismatch"] = counters["exp_mismatch"]
    summary["n_dark_refused"] = counters["dark_refused"]
    summary["n_nights"] = len(nights)
    summary["calibrated_root"] = str(out_dir / "calibrated")
    (out_dir / "calibration.json").write_text(
        json.dumps(summary, indent=2, default=float), encoding="utf-8")
    return summary
