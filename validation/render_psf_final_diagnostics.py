"""Render the same final PSF diagnostic figure exposed by the Step 8 GUI."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from astropy.io import fits
from matplotlib.figure import Figure
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from apex.analysis.psf_diagnostics import (
    draw_psf_final_diagnostics,
    load_psf_final_diagnostic_data,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path)
    parser.add_argument(
        "--frame",
        help="FITS filename to render when the result directory contains multiple frames",
    )
    parser.add_argument("--fwhm-px", type=float)
    parser.add_argument("--pixel-scale", type=float, default=np.nan)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    psf_dir = args.result_dir / "cmd_psf"
    if args.frame:
        meta_path = psf_dir / f"residual_meta_{args.frame}.json"
        if not meta_path.exists():
            raise FileNotFoundError(f"PSF metadata not found: {meta_path}")
    else:
        meta_path = next(psf_dir.glob("residual_meta_*.json"))
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    filename = str(meta["file"])
    data = load_psf_final_diagnostic_data(args.result_dir, filename)
    frame_stem = Path(filename).stem
    epsf_path = next(psf_dir.glob(f"epsf_model_*_{frame_stem}.fits"))
    epsf_model = np.asarray(fits.getdata(epsf_path), dtype=float)

    detect_path = args.result_dir / "cache" / f"detect_{filename}.json"
    detect = {}
    if detect_path.exists():
        detect = json.loads(detect_path.read_text(encoding="utf-8"))
    fwhm_px = args.fwhm_px
    if fwhm_px is None:
        fwhm_px = float(detect.get("fwhm_px", np.nan))
    pixel_scale = float(args.pixel_scale)
    if not np.isfinite(pixel_scale):
        fwhm_arcsec = float(detect.get("fwhm_arcsec", np.nan))
        if np.isfinite(fwhm_arcsec) and np.isfinite(fwhm_px) and fwhm_px > 0:
            pixel_scale = fwhm_arcsec / fwhm_px

    reference = None
    reference_meta = meta.get("epsf_reference", {})
    reference_name = reference_meta.get("catalog_path")
    if reference_name and (psf_dir / reference_name).exists():
        reference = pd.read_csv(psf_dir / reference_name)

    core = meta.get("core_cut", {})
    figure = Figure(figsize=(16.0, 9.3))
    summary = draw_psf_final_diagnostics(
        figure,
        data,
        epsf_model,
        filename=filename,
        fwhm_px=fwhm_px,
        pixel_scale_arcsec=pixel_scale,
        core_center=(
            float(core.get("center_x", np.nan)),
            float(core.get("center_y", np.nan)),
        ),
        core_radius_px=float(core.get("radius_px", np.nan)),
        epsf_reference=reference,
    )

    output = args.output or (args.result_dir / "psf_final_diagnostics.png")
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=200, facecolor="white")
    summary_path = output.with_suffix(".json")
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(output.resolve())
    print(summary_path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
