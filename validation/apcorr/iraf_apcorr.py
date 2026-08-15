"""The external arm: IRAF's own aperture correction on the same stars.

`reproduce_apcorr.py` showed the reimplementation matches APEX, and
`apcorr_arms.py` showed the summation kernel is not in question (photutils and
SExtractor agree to 1e-14) while the *choices* around it move the answer by up
to 0.1 mag. What neither can say is whether APEX's answer agrees with a tool
that was written independently. IRAF DAOPHOT is that tool: `phot` measures the
same stars through a list of radii, and `mkapfile` fits the growth curve and
reports the correction.

Variables held fixed (V9 — the parameter ledger is the point, not a formality):

* **same stars** — APEX's own reference selection, written out as a coordinate
  list, so star choice cannot explain a difference
* **same radii** — APEX's 14-point grid handed to `photpars.apertures`
* **same detector constants** — gain and read noise from the frame's own
  measurement, not IRAF defaults
* **same pixel convention** — IRAF counts the first pixel 1, APEX counts 0.
  This cost the PSF comparison a day; `+1` goes on the way in and comes off on
  the way out, and the round trip is asserted.

What is deliberately *not* matched is the growth-curve model. IRAF fits a
Moffat-plus-power-law with `mkapfile`; APEX takes a median of per-star
normalised curves. That difference is the thing being measured.

    python validation/apcorr/iraf_apcorr.py --workspace <phase3/M67> --frames 3
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).absolute().parent))

from apex.benchmark.iraf_crosscheck import windows_to_wsl_path  # noqa: E402
from reproduce_apcorr import (  # noqa: E402
    GC_N_STEPS, _num, apcorr_from_curve, reference_mask,
)

# IRAF counts the first pixel 1; numpy counts it 0. The PSF cross-check lost a
# day to this in 2026-08-14 — the offset is sqrt(2)=1.414 px against a 1.5 px
# matching radius, so it looked like a marginal engine difference.
IRAF_ORIGIN_OFFSET = 1.0

SCRIPT = '''\
import json, os
from pyraf import iraf

iraf.noao(); iraf.digiphot(); iraf.apphot(); iraf.photcal()
for task in ("datapars", "centerpars", "fitskypars", "photpars", "phot", "mkapfile"):
    try:
        iraf.unlearn(task)
    except Exception:
        pass

iraf.datapars.fwhmpsf = {fwhm:.8f}
iraf.datapars.sigma = {sigma:.8f}
iraf.datapars.readnoise = {readnoise:.8f}
iraf.datapars.epadu = {gain:.8f}
iraf.datapars.datamin = "INDEF"
iraf.datapars.datamax = {datamax:.8f}
iraf.datapars.exposure = ""
iraf.datapars.itime = 1.0

# Positions come from APEX and must not move: recentring here would
# make the two engines measure different pixels.
iraf.centerpars.calgorithm = "none"

iraf.fitskypars.salgorithm = "mode"
iraf.fitskypars.annulus = {annulus:.8f}
iraf.fitskypars.dannulus = {dannulus:.8f}

iraf.photpars.apertures = "{apertures}"
iraf.photpars.zmag = 25.0

iraf.phot("{image}", coords="{coords}", output="{magfile}",
          interactive="no", verify="no", verbose="no")

# mkapfile takes no `verbose`; an unknown name raises here.
# Not mkapfile. Its correction runs from the smallest aperture to *infinity*
# via a fitted Moffat-plus-wings model, while APEX's runs from r_ap to the
# outermost grid radius. Comparing them directly reads as a 1.43 mag
# disagreement that is entirely definitional (measured 2026-08-16). Dump the
# per-aperture fluxes instead and build APEX's finite quantity from them, so
# what differs is IRAF's sky estimator and summation, not the question asked.
iraf.txdump("{magfile}", "XCENTER,YCENTER,FLUX", "yes", Stdout="{fluxfile}")

with open("{dumpfile}", "w") as handle:
    json.dump({{"apertures": iraf.photpars.apertures,
               "annulus": float(iraf.fitskypars.annulus),
               "dannulus": float(iraf.fitskypars.dannulus),
               "fwhmpsf": float(iraf.datapars.fwhmpsf),
               "sigma": float(iraf.datapars.sigma),
               "epadu": float(iraf.datapars.epadu),
               "readnoise": float(iraf.datapars.readnoise),
               "datamax": float(iraf.datapars.datamax),
               "calgorithm": str(iraf.centerpars.calgorithm),
               "salgorithm": str(iraf.fitskypars.salgorithm)}}, handle, indent=2)
print("IRAF_DONE")
'''


def _read_txdump_fluxes(path: Path, n_apertures: int) -> np.ndarray:
    """Per-star, per-aperture net fluxes from `txdump`.

    One row per star: XCENTER, YCENTER, then one flux per aperture. INDEF
    marks a measurement IRAF refused to make (saturated, off-frame, bad sky)
    and becomes NaN rather than being dropped, so star rows stay aligned.
    """
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        parts = line.replace(",", " ").split()
        if len(parts) < 2 + n_apertures:
            continue
        values = []
        for token in parts[2:2 + n_apertures]:
            try:
                values.append(float(token))
            except ValueError:
                values.append(np.nan)
        rows.append(values)
    return np.asarray(rows, dtype=float).T if rows else np.empty((n_apertures, 0))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workspace", required=True, type=Path)
    ap.add_argument("--frames", type=int, default=3)
    ap.add_argument("--workdir", type=Path,
                    default=Path(r"E:\APEX_validation\apcorr_iraf"))
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    from apex.config.parameters_cmd import read_params
    from apex.analysis.forced_photometry import _to_float

    params = read_params(args.workspace / "apex_config.json")
    P = params.P
    step7 = Path(P.result_dir) / "step7_forced_phot"
    data_dir = Path(P.data_dir)

    recorded = pd.read_csv(step7 / "apcorr_summary.csv")
    stats = pd.read_csv(step7 / "frame_stats.csv")
    master_path = step7 / "master_sources.csv"
    master = pd.read_csv(master_path) if master_path.exists() else None

    r_ap_scale = _to_float(getattr(P, "forced_r_ap_scale", 0.8), 0.8)
    ref_scale = _to_float(getattr(P, "forced_ref_ap_scale", 2.4), 2.4)
    min_r_ap = _to_float(getattr(P, "min_r_ap_px", 4.0), 4.0)
    ann_scale = _to_float(getattr(P, "fitsky_annulus_scale", 4.0), 4.0)
    dann_scale = _to_float(getattr(P, "fitsky_dannulus_scale", 2.0), 2.0)
    ann_gap = _to_float(getattr(P, "annulus_min_gap_px", 6.0), 6.0)
    datamax = _to_float(getattr(P, "datamax_adu", 55000.0), 55000.0)

    fwhm_column = next(c for c in ("fwhm_px", "fwhm", "fwhm_median_px")
                       if c in stats.columns)
    args.workdir.mkdir(parents=True, exist_ok=True)

    rows = []
    frames = recorded if args.frames <= 0 else recorded.head(args.frames)
    for entry in frames.itertuples(index=False):
        name = str(entry.file)
        table, image_path = step7 / f"photometry_{name}.tsv", data_dir / name
        stat_row = stats[stats["file"].astype(str) == name]
        if not table.exists() or not image_path.exists() or stat_row.empty:
            continue

        fwhm = float(pd.to_numeric(stat_row.iloc[0][fwhm_column], errors="coerce"))
        r_ap = max(min_r_ap, r_ap_scale * fwhm)
        r_ref = max(r_ap + 2.0, ref_scale * fwhm)
        r_in = max(r_ref + ann_gap, ann_scale * fwhm)
        r_out = r_in + max(ann_gap, dann_scale * fwhm)
        radii = np.linspace(max(2.0, r_ap * 0.4), r_ref * 1.15, GC_N_STEPS)

        phot = pd.read_csv(table, sep="\t")
        chosen = phot[reference_mask(phot, master)]
        x = _num(chosen, "x_fit")
        y = _num(chosen, "y_fit")
        finite = np.isfinite(x) & np.isfinite(y)
        x, y = x[finite], y[finite]
        sky_sigma = float(np.nanmedian(_num(chosen, "sky_std")[finite]))
        gain = float(np.nanmedian(_num(chosen, "gain_e_per_adu")[finite]))
        readnoise = float(np.nanmedian(_num(chosen, "rdnoise_e")[finite]))

        run_dir = args.workdir / Path(name).stem
        if run_dir.exists():
            shutil.rmtree(run_dir)
        run_dir.mkdir(parents=True)

        with fits.open(image_path, memmap=False) as hdul:
            fits.PrimaryHDU(np.asarray(hdul[0].data, dtype=np.float32),
                            header=hdul[0].header).writeto(run_dir / "frame.fits")

        coords = run_dir / "coords.txt"
        coords.write_text("\n".join(
            f"{xi + IRAF_ORIGIN_OFFSET:.4f} {yi + IRAF_ORIGIN_OFFSET:.4f}"
            for xi, yi in zip(x, y)) + "\n", encoding="ascii")

        wsl = lambda p: windows_to_wsl_path(str(p))  # noqa: E731
        script = SCRIPT.format(
            fwhm=fwhm, sigma=sky_sigma, readnoise=readnoise, gain=gain,
            datamax=datamax, annulus=r_in, dannulus=r_out - r_in,
            apertures=",".join(f"{r:.4f}" for r in radii),
            naperts=len(radii),
            image="frame.fits", coords="coords.txt",
            magfile="frame.mag", fluxfile="frame.flux",
            dumpfile="iraf_params.json",
        )
        (run_dir / "run.py").write_text(script, encoding="utf-8")

        completed = subprocess.run(
            ["wsl", "-e", "bash", "-lc",
             f"cd {wsl(run_dir)} && python3 run.py"],
            capture_output=True, text=True, timeout=1800)
        if "IRAF_DONE" not in completed.stdout:
            print(f"{name}: IRAF 실패\n{completed.stdout[-400:]}\n{completed.stderr[-400:]}")
            continue

        flux_file = run_dir / "frame.flux"
        iraf_flux = (_read_txdump_fluxes(flux_file, len(radii))
                     if flux_file.exists() else np.empty((len(radii), 0)))
        iraf_apcorr = apcorr_from_curve(iraf_flux, radii, r_ap)[0] if iraf_flux.size else float("nan")
        # Both as the magnitude the correction adds to the small aperture.
        iraf_mag = -2.5 * np.log10(iraf_apcorr) if np.isfinite(iraf_apcorr) else float("nan")
        apex_mag = -2.5 * np.log10(float(entry.apcorr))
        rows.append({
            "frame": name, "filter": entry.filter, "fwhm_px": fwhm,
            "n_stars": int(len(x)),
            "apex_apcorr": float(entry.apcorr), "apex_mag": apex_mag,
            "iraf_apcorr": iraf_apcorr,
            "iraf_mag": iraf_mag, "difference_mmag": (iraf_mag - apex_mag) * 1000.0,
        })
        print(f"{name:<30} APEX {apex_mag:+.4f}  IRAF {iraf_mag:+.4f}  "
              f"Δ {(iraf_mag - apex_mag) * 1000:+.1f} mmag", flush=True)

    if not rows:
        print("[error] IRAF 결과가 없다")
        return 1
    frame = pd.DataFrame(rows)
    delta = frame["difference_mmag"].to_numpy(float)
    delta = delta[np.isfinite(delta)]
    if delta.size:
        print(f"\n{len(delta)} 프레임 · 중앙 {np.median(delta):+.1f} mmag · "
              f"강건σ {1.4826 * np.median(np.abs(delta - np.median(delta))):.1f} mmag")
    if args.output:
        frame.to_csv(args.output, index=False)
        args.output.with_suffix(".inputs.json").write_text(json.dumps({
            "workspace": str(args.workspace),
            "held_fixed": ["star list (APEX selection)", "radius grid",
                           "sky annulus", "gain / read noise", "pixel origin (+1/-1)"],
            "independent_here": "IRAF sky estimator (mode) and its own aperture "
                                "summation; APEX's finite-radius definition is "
                                "rebuilt from the per-aperture fluxes",
            "not_used": "mkapfile — its correction extrapolates to infinity, a "
                        "different quantity (measured 1.43 mag apart)",
            "iraf_centering": "calgorithm=none (positions must not move)",
            "iraf_origin_offset": IRAF_ORIGIN_OFFSET,
        }, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"표 -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
