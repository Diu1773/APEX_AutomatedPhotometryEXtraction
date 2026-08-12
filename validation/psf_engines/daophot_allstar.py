"""Run IRAF DAOPHOT PSF + ALLSTAR on a frame, for comparison against Step 8.

The engine scorecard scores APEX's PSF photometry ✕ on its A axis — "was it
measured against an existing engine?" — and the reason is narrow: APEX already
has a 2,591-line IRAF cross-check harness, but that harness only ever calls
`phot`, IRAF's *aperture* task. It loads the daophot package and never uses
`psf` or `allstar`. Step 7 has been compared to IRAF; step 8 never has.

This module closes that. It drives the classical DAOPHOT chain
(`phot` → `pstselect` → `psf` → `allstar`) through PyRAF under WSL, on the same
pixels APEX measured, and writes a flat table that
`compare_recovery.py` can put beside APEX's.

Fairness
--------
Two kinds of number are set here, and they are set differently on purpose.

*Detector constants* — gain, read noise, the good-data window, exposure time —
describe the data, not the method. They are read from APEX's own configuration
so both engines are told the same thing about the same pixels. Letting these
differ would make the comparison meaningless.

*Method parameters* — the PSF radius, the fitting radius, the analytic function
— are choices each engine's authors made. Matching them to APEX would be
tuning DAOPHOT toward APEX's answer. They follow the standard DAOPHOT guidance
instead (Massey & Davis 1992, *A User's Guide to Stellar CCD Photometry with
IRAF*): `psfrad ≈ 4·FWHM + 1`, `fitrad ≈ FWHM`, `function = auto`. Because
APEX's fit window is larger (an automatic window holding 90 % of the encircled
energy, about 1.7·FWHM in these frames), `--fitrad-fwhm` runs the same chain at
APEX's radius as a sensitivity check, so the reader can see whether any
difference is the engine or the window.

Star positions are supplied rather than detected. APEX's step 8 receives
positions from steps 4 and 7 and does not detect independently, so feeding
DAOPHOT the same list compares PSF photometry to PSF photometry instead of
folding two detection stages into the result. Step 8 can still add sources
during its residual passes; that difference is reported, not hidden.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

REPO = Path(__file__).absolute().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from apex.benchmark.iraf_crosscheck import windows_to_wsl_path  # noqa: E402

# DAOPHOT writes into the working directory and refuses to overwrite, so every
# run gets a clean directory and short ASCII names.
IMAGE_STEM = "frame"

PYRAF_TEMPLATE = '''\
import os, time, json
from pyraf import iraf

iraf.noao(); iraf.digiphot(); iraf.daophot()
for task in ("datapars", "daopars", "centerpars", "fitskypars", "photpars",
             "findpars", "phot", "pstselect", "psf", "allstar"):
    try:
        iraf.unlearn(task)
    except Exception:
        pass

# --- detector constants: these describe the data, taken from APEX's config ---
iraf.datapars.noise = "poisson"
iraf.datapars.fwhmpsf = {fwhm:.8f}
iraf.datapars.sigma = {sigma:.8f}
iraf.datapars.readnoise = {readnoise:.8f}
iraf.datapars.epadu = {gain:.8f}
iraf.datapars.itime = {exptime:.8f}
iraf.datapars.datamin = {datamin:.8f}
iraf.datapars.datamax = {datamax:.8f}

# --- method parameters: standard DAOPHOT guidance, not matched to APEX -------
# `phot` only supplies ALLSTAR's initial sky and magnitude at positions that
# are already known, so it must not move them. With recentering on, a faint
# star sitting a fraction of a FWHM from a bright neighbour has its centroid
# dragged onto the neighbour and is then lost -- measured here as 7 of 25
# implanted stars surviving `phot` while 1575 of 1599 real ones did. ALLSTAR
# still refines positions during the fit (daopars.recenter), which is the
# stage where APEX refines them too.
iraf.centerpars.calgorithm = "{calgorithm}"
iraf.centerpars.cbox = {cbox:.8f}
iraf.fitskypars.salgorithm = "mode"
iraf.fitskypars.annulus = {annulus:.8f}
iraf.fitskypars.dannulus = {dannulus:.8f}
iraf.photpars.apertures = "{aperture:.4f}"
iraf.photpars.zmag = {zmag:.4f}

iraf.daopars.function = "{function}"
iraf.daopars.varorder = {varorder}
iraf.daopars.psfrad = {psfrad:.8f}
iraf.daopars.fitrad = {fitrad:.8f}
iraf.daopars.recenter = "yes"
iraf.daopars.fitsky = "yes"
iraf.daopars.sannulus = {annulus:.8f}
iraf.daopars.wsannulus = {dannulus:.8f}
iraf.daopars.maxiter = {maxiter}
iraf.daopars.maxnstar = 20000
iraf.daopars.nclean = {nclean}
iraf.daopars.mergerad = "INDEF"

# Every parameter the five tasks will actually use, dumped from IRAF itself
# rather than from the assignments above. A default that `unlearn` restored, or
# a value IRAF coerced, would otherwise never reach the published table — and a
# comparison against another package is not usable in a paper without it.
import re as _re
_pars = {{}}
_ROW = _re.compile(r"^\\s*\\(?([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.*?)\\)?\\s{{2,}}(.*)$")
for _task in ("datapars", "daopars", "centerpars", "fitskypars", "photpars"):
    _dump = "%s.par" % _task
    iraf.lpar(_task, Stdout=_dump)
    _pars[_task] = {{}}
    for _line in open(_dump, encoding="utf-8", errors="replace"):
        _m = _ROW.match(_line.rstrip("\\n"))
        if _m:
            _pars[_task][_m.group(1)] = {{"value": _m.group(2).strip(),
                                         "description": _m.group(3).strip()}}
json.dump(_pars, open("iraf_parameters.json", "w"), indent=1, sort_keys=True)

IMAGE = "{image}"
COORDS = "{coords}"
stages = {{}}

def stage(name, fn):
    start = time.time()
    fn()
    stages[name] = time.time() - start

stage("phot", lambda: iraf.phot(
    IMAGE, coords=COORDS, output="phot.mag", verify="no", interactive="no",
    verbose="no", Stdout=os.devnull))

stage("pstselect", lambda: iraf.pstselect(
    IMAGE, photfile="phot.mag", pstfile="psf.pst", maxnpsf={maxnpsf},
    verify="no", interactive="no", verbose="no", Stdout=os.devnull))

stage("psf", lambda: iraf.psf(
    IMAGE, photfile="phot.mag", pstfile="psf.pst", psfimage="psf.fits",
    opstfile="psf.opst", groupfile="psf.grp",
    verify="no", interactive="no", verbose="no", Stdout=os.devnull))

stage("allstar", lambda: iraf.allstar(
    IMAGE, photfile="phot.mag", psfimage="psf.fits", allstarfile="all.als",
    rejfile="all.arj", subimage="sub.fits",
    verify="no", verbose="no", Stdout=os.devnull))

iraf.txdump("all.als", fields="ID,XCENTER,YCENTER,MAG,MERR,MSKY,NITER,SHARPNESS,CHI",
            expr="yes", headers="no", Stdout="allstar.txt")
iraf.txdump("phot.mag", fields="ID,XCENTER,YCENTER,MAG,MERR,MSKY",
            expr="yes", headers="no", Stdout="phot.txt")

n_psf_stars = sum(1 for line in open("psf.opst")
                  if line.strip() and not line.startswith("#"))
json.dump({{"stages": stages, "n_psf_stars": n_psf_stars}},
          open("timing.json", "w"))
print("DAOPHOT_OK")
'''

ALLSTAR_COLUMNS = ("id", "x", "y", "mag", "merr", "msky", "niter",
                   "sharpness", "chi")


def frame_statistics(path: Path) -> dict:
    """Sky level and noise, measured the same way APEX measures them."""
    data = fits.getdata(path).astype(float)
    header = fits.getheader(path)
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    return {
        "sky_mean": float(mean), "sky_median": float(median),
        "sky_sigma": float(std),
        "exptime": float(header.get("EXPTIME", header.get("EXPOSURE", 1.0))),
    }


def read_positions(path: Path) -> pd.DataFrame:
    """Star positions to measure. Accepts an APEX step-7 TSV or a plain CSV."""
    sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    table = pd.read_csv(path, sep=sep)
    names = {c.lower(): c for c in table.columns}
    for x_key, y_key in (("x", "y"), ("x_fit", "y_fit"),
                         ("xcentroid", "ycentroid"), ("x_pix", "y_pix")):
        if x_key in names and y_key in names:
            out = pd.DataFrame({
                "x": pd.to_numeric(table[names[x_key]], errors="coerce"),
                "y": pd.to_numeric(table[names[y_key]], errors="coerce"),
            })
            return out[np.isfinite(out["x"]) & np.isfinite(out["y"])]
    raise SystemExit(f"좌표 열을 못 찾았다: {sorted(table.columns)[:12]}")


def run(args: argparse.Namespace) -> dict:
    work = Path(args.workdir)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    # IRAF is unhappy with a `.fit` extension and with long Windows paths, so
    # the frame is copied in under a short name.
    image = work / f"{IMAGE_STEM}.fits"
    shutil.copy2(args.frame, image)

    stats = frame_statistics(image)
    positions = read_positions(Path(args.positions))
    coords = work / "stars.coo"
    positions.to_csv(coords, sep=" ", header=False, index=False,
                     float_format="%.4f")
    print(f"프레임 {image.name} · 좌표 {len(positions)}개 · "
          f"sky {stats['sky_median']:.1f} ± {stats['sky_sigma']:.2f} ADU")

    fwhm = float(args.fwhm)
    fitrad = args.fitrad_fwhm * fwhm if args.fitrad_fwhm else fwhm
    script = work / "run_daophot.py"
    script.write_text(PYRAF_TEMPLATE.format(
        fwhm=fwhm,
        sigma=stats["sky_sigma"],
        readnoise=args.readnoise,
        gain=args.gain,
        exptime=stats["exptime"],
        datamin=args.datamin,
        datamax=args.datamax,
        calgorithm=args.calgorithm,
        cbox=2.0 * fwhm,
        annulus=args.annulus_fwhm * fwhm,
        dannulus=args.dannulus_fwhm * fwhm,
        aperture=args.aperture_fwhm * fwhm,
        zmag=args.zmag,
        function=args.function,
        varorder=args.varorder,
        psfrad=4.0 * fwhm + 1.0,
        fitrad=fitrad,
        maxiter=args.maxiter, nclean=args.nclean,
        maxnpsf=args.maxnpsf,
        image=f"{IMAGE_STEM}.fits",
        coords="stars.coo",
    ), encoding="utf-8")

    runtime = ["wsl", "python3"] if sys.platform == "win32" else ["python3"]
    started = time.perf_counter()
    completed = subprocess.run(
        [*runtime, windows_to_wsl_path(script) if runtime[0] == "wsl"
         else str(script)],
        cwd=str(work), text=True, capture_output=True, check=False,
    )
    elapsed = time.perf_counter() - started
    (work / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (work / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if "DAOPHOT_OK" not in completed.stdout:
        print(completed.stdout[-2000:])
        print(completed.stderr[-3000:])
        raise SystemExit(f"DAOPHOT 실패 (returncode={completed.returncode})")

    table = pd.read_csv(work / "allstar.txt", sep=r"\s+", header=None,
                        names=ALLSTAR_COLUMNS, na_values=["INDEF"])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)

    timing = json.loads((work / "timing.json").read_text(encoding="utf-8"))
    # The parameter table belongs next to the numbers it produced, not only in
    # the scratch work directory that a later run overwrites.
    iraf_parameters = json.loads(
        (work / "iraf_parameters.json").read_text(encoding="utf-8"))
    (output.parent / f"{output.stem}_iraf_parameters.json").write_text(
        json.dumps(iraf_parameters, indent=1, sort_keys=True), encoding="utf-8")
    fitted = int(np.isfinite(pd.to_numeric(table["mag"], errors="coerce")).sum())
    summary = {
        "frame": str(args.frame),
        "n_input_positions": int(len(positions)),
        "n_allstar_rows": int(len(table)),
        "n_valid_mag": fitted,
        "n_psf_stars": timing.get("n_psf_stars"),
        "elapsed_s": elapsed,
        "stage_seconds": timing.get("stages", {}),
        "parameters": {
            "fwhm_px": fwhm, "psfrad_px": 4.0 * fwhm + 1.0,
            "fitrad_px": fitrad, "fitrad_fwhm": fitrad / fwhm,
            "function": args.function, "varorder": args.varorder,
            "nclean": args.nclean,
            "gain_e_per_adu": args.gain, "readnoise_e": args.readnoise,
            "datamin_adu": args.datamin, "datamax_adu": args.datamax,
            "sky_sigma_adu": stats["sky_sigma"], "zmag": args.zmag,
            "exptime": stats["exptime"],
        },
        "output": str(output),
        "iraf_parameters": iraf_parameters,
    }
    (output.parent / f"{output.stem}_summary.json").write_text(
        json.dumps(summary, indent=1), encoding="utf-8")

    print(f"ALLSTAR {fitted}/{len(positions)} 적합 · {elapsed:.0f}s "
          f"(PSF 기준성 {timing.get('n_psf_stars')}개, "
          f"fitrad {fitrad / fwhm:.2f}xFWHM)")
    for name, seconds in timing.get("stages", {}).items():
        print(f"  {name:>10}: {seconds:6.1f}s")
    print(f"saved -> {output}")
    return summary


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--frame", required=True)
    ap.add_argument("--positions", required=True,
                    help="APEX step7 TSV 또는 x,y CSV")
    ap.add_argument("--output", required=True)
    ap.add_argument("--workdir", default=r"E:\APEX_validation\psf_engines\daophot_work")
    ap.add_argument("--fwhm", type=float, required=True)
    # Detector constants default to the M13 workspace values; pass explicitly
    # for any other instrument.
    ap.add_argument("--gain", type=float, default=0.68)
    ap.add_argument("--readnoise", type=float, default=2.35)
    ap.add_argument("--datamin", type=float, default=0.1)
    ap.add_argument("--datamax", type=float, default=55000.0)
    ap.add_argument("--zmag", type=float, default=25.0)
    # Method parameters — DAOPHOT guidance unless deliberately overridden.
    ap.add_argument("--aperture-fwhm", type=float, default=1.0)
    ap.add_argument("--annulus-fwhm", type=float, default=4.0)
    ap.add_argument("--dannulus-fwhm", type=float, default=2.0)
    ap.add_argument("--fitrad-fwhm", type=float, default=None,
                    help="기본은 DAOPHOT 표준(1xFWHM). APEX 창(약 1.7)으로 "
                         "맞춰 보려면 명시할 것")
    ap.add_argument("--calgorithm", default="none",
                    choices=("none", "centroid", "gauss", "ofilter"),
                    help="`phot` 단계의 재중심. 좌표를 주는 강제측광이므로 "
                         "기본은 none — 켜면 희미한 별의 중심이 이웃으로 끌려간다")
    ap.add_argument("--function", default="auto")
    ap.add_argument("--varorder", type=int, default=0,
                    help="0 = 시야 내 일정. APEX 는 프레임당 ePSF 하나이므로 "
                         "0 이 대응된다")
    ap.add_argument("--maxiter", type=int, default=50)
    ap.add_argument("--nclean", type=int, default=0,
                    help="psf 가 PSF 별을 정제하는 반복 횟수. 0 은 정제 없음이라 "
                         "DAOPHOT 에 불리하다 — 사람이 눈으로 걸러내던 단계의 대체물")
    ap.add_argument("--maxnpsf", type=int, default=60)
    args = ap.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
