"""Command-line entry point for IRAF daofind batch cross-checks."""

from __future__ import annotations

import argparse

from apex.benchmark.iraf_crosscheck_batch import (
    IRAFCrosscheckBatchConfig,
    run_iraf_crosscheck_batch,
)


def _float_list(text: str) -> list[float]:
    return [float(part.strip()) for part in text.split(",") if part.strip()]


def _str_list(text: str) -> list[str]:
    return [part.strip() for part in text.split(",") if part.strip()]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run IRAF/DAOPHOT daofind+phot against all selected APEX CMD frames."
    )
    parser.add_argument("--project-root", required=True, help="Completed APEX CMD project root.")
    parser.add_argument("--output", default="benchmark/runs/iraf_daofind_batch", help="Output directory.")
    parser.add_argument("--filters", default="g,r,i", help="Comma-separated filters.")
    parser.add_argument("--threshold-grid", default="12,9,7,5", help="Comma-separated IRAF thresholds.")
    parser.add_argument("--min-snr", type=float, default=20.0)
    parser.add_argument("--zmag", type=float, default=25.0)
    parser.add_argument("--daofind-max-sources", type=int, default=2500)
    parser.add_argument("--daofind-max-ratio-to-apex", type=float, default=2.0)
    parser.add_argument("--match-radius-px", type=float, default=2.0)
    parser.add_argument("--limit", type=int, help="Limit selected frames for smoke runs.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--no-resume", action="store_true", help="Do not skip completed frame outputs.")
    parser.add_argument(
        "--runtime-cmd",
        nargs="+",
        default=None,
        help="Runtime command before the generated PyRAF script, e.g. wsl python3.",
    )
    args = parser.parse_args()

    config = IRAFCrosscheckBatchConfig(
        project_root=args.project_root,
        output_root=args.output,
        filters=_str_list(args.filters),
        threshold_grid=_float_list(args.threshold_grid),
        min_snr=float(args.min_snr),
        zmag=float(args.zmag),
        daofind_max_sources=int(args.daofind_max_sources),
        daofind_max_ratio_to_apex=float(args.daofind_max_ratio_to_apex),
        match_radius_px=float(args.match_radius_px),
        limit=args.limit,
        overwrite=bool(args.overwrite),
        resume=not bool(args.no_resume),
        runtime_cmd=list(args.runtime_cmd or []),
    )
    output_dir = run_iraf_crosscheck_batch(config)
    print(f"IRAF daofind batch complete: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

