"""Command-line entry point for CMD batch benchmarks."""

from __future__ import annotations

import argparse

from apex.benchmark.cmd_batch import load_cmd_batch_config, run_cmd_batch


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Select CMD quality frames and run APEX artificial-star benchmarks."
    )
    parser.add_argument(
        "--config",
        default="benchmark/configs/cmd_batch.toml",
        help="CMD batch TOML configuration.",
    )
    parser.add_argument("--project-root", help="Override the CMD project root.")
    parser.add_argument("--output", help="Override the batch output directory.")
    parser.add_argument(
        "--select-only",
        action="store_true",
        help="Write the selected best/median/worst frames without running injections.",
    )
    args = parser.parse_args()

    config = load_cmd_batch_config(args.config)
    output_dir = run_cmd_batch(
        config,
        project_override=args.project_root,
        output_override=args.output,
        select_only=args.select_only,
    )
    print(f"CMD batch complete: {output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
