"""Command-line entry point for CMD validation reports."""

from __future__ import annotations

import argparse

from apex.benchmark.cmd_validation import CmdValidationConfig, run_cmd_validation


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build a CMD validation report from APEX benchmark outputs."
    )
    parser.add_argument("--batch-root", required=True, help="Completed CMD batch output directory.")
    parser.add_argument("--project-root", required=True, help="Completed CMD project root.")
    parser.add_argument("--output", default="benchmark/runs/cmd_validation", help="Output directory.")
    parser.add_argument("--repeatability-min-frames", type=int, default=3)
    parser.add_argument("--repeatability-min-snr", type=float, default=5.0)
    parser.add_argument("--magnitude-bin-width", type=float, default=0.5)
    args = parser.parse_args()

    config = CmdValidationConfig(
        batch_root=args.batch_root,
        project_root=args.project_root,
        output_root=args.output,
        repeatability_min_frames=args.repeatability_min_frames,
        repeatability_min_snr=args.repeatability_min_snr,
        magnitude_bin_width=args.magnitude_bin_width,
    )
    output_dir = run_cmd_validation(config)
    print(f"CMD validation complete: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
