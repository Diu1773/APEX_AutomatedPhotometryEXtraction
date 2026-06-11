"""Command-line entry point for combined CMD validation reports."""

from __future__ import annotations

import argparse

from apex.benchmark.cmd_combined_validation import (
    CmdCombinedValidationConfig,
    parse_validation_input,
    run_combined_validation,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Combine multiple CMD validation reports into one summary."
    )
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="Validation input in label=path or label:role=path form. Repeatable.",
    )
    parser.add_argument(
        "--output",
        default="benchmark/runs/cmd_validation_combined",
        help="Output directory.",
    )
    args = parser.parse_args()

    config = CmdCombinedValidationConfig(
        inputs=[parse_validation_input(value) for value in args.input],
        output_root=args.output,
    )
    output_dir = run_combined_validation(config)
    print(f"Combined CMD validation complete: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
