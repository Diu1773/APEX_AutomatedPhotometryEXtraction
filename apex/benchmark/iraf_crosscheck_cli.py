"""Command-line entry point for IRAF/PyRAF cross-check benchmarks."""

from __future__ import annotations

import argparse

from apex.benchmark.iraf_crosscheck import IRAFCrosscheckConfig, run_iraf_crosscheck


def _float_list(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare APEX Step 7 photometry against IRAF/DAOPHOT via PyRAF."
    )
    parser.add_argument("--input", required=True, help="Input science FITS frame.")
    parser.add_argument("--step7", required=True, help="APEX Step 7 photometry TSV for the frame.")
    parser.add_argument("--output", default="benchmark/runs/iraf_crosscheck", help="Output directory.")
    parser.add_argument(
        "--mode",
        choices=["phot_fixed_coords", "daofind_phot", "both"],
        default="both",
        help="Cross-check mode.",
    )
    parser.add_argument("--filter", default="", help="Filter label override.")
    parser.add_argument("--fwhm-px", type=float, help="FWHM override in pixels.")
    parser.add_argument("--max-fixed-coords", type=int, default=500)
    parser.add_argument("--min-snr", type=float, default=20.0)
    parser.add_argument("--zmag", type=float, default=25.0)
    parser.add_argument("--gain", type=float, help="IRAF datapars.epadu override.")
    parser.add_argument("--rdnoise", type=float, help="IRAF datapars.readnoise override.")
    parser.add_argument("--sigma-adu", type=float, help="IRAF datapars.sigma override.")
    parser.add_argument("--recenter", action="store_true", help="Allow IRAF phot centroid recentering.")
    parser.add_argument("--aperture-scale", type=float, default=0.8)
    parser.add_argument("--annulus-scale", type=float, default=4.0)
    parser.add_argument("--dannulus-scale", type=float, default=2.0)
    parser.add_argument("--threshold-grid", default="12,9,7,5")
    parser.add_argument("--daofind-max-sources", type=int, default=2500)
    parser.add_argument("--daofind-max-ratio-to-apex", type=float, default=2.0)
    parser.add_argument("--match-radius-px", type=float, default=2.0)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--runtime-cmd",
        nargs="+",
        default=None,
        help="Runtime command before the generated PyRAF script, e.g. wsl python3.",
    )
    args = parser.parse_args()

    config = IRAFCrosscheckConfig(
        input_fits=args.input,
        step7_tsv=args.step7,
        output_root=args.output,
        mode=args.mode,
        filter_name=args.filter,
        fwhm_px=args.fwhm_px,
        aperture_scale_fwhm=args.aperture_scale,
        annulus_scale_fwhm=args.annulus_scale,
        dannulus_scale_fwhm=args.dannulus_scale,
        max_fixed_coords=args.max_fixed_coords,
        min_snr=args.min_snr,
        zmag=args.zmag,
        gain_e_per_adu=args.gain,
        rdnoise_e=args.rdnoise,
        sigma_adu=args.sigma_adu,
        recenter=bool(args.recenter),
        daofind_threshold_grid=_float_list(args.threshold_grid),
        daofind_max_sources=args.daofind_max_sources,
        daofind_max_ratio_to_apex=args.daofind_max_ratio_to_apex,
        daofind_match_radius_px=args.match_radius_px,
        overwrite=bool(args.overwrite),
        runtime_cmd=list(args.runtime_cmd or []),
    )
    output_dir = run_iraf_crosscheck(config)
    print(f"IRAF cross-check complete: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

