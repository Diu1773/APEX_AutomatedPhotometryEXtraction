"""Command-line entry point for APEX benchmarks."""

from __future__ import annotations

import argparse

from apex.benchmark.runner import load_benchmark_config, run_benchmark


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run an automated artificial-star benchmark using APEX production code."
    )
    parser.add_argument(
        "--config",
        default="benchmark/configs/baseline.toml",
        help="Benchmark TOML configuration.",
    )
    parser.add_argument("--input", help="Override the input FITS path.")
    parser.add_argument("--output", help="Override the output directory.")
    args = parser.parse_args()

    config = load_benchmark_config(args.config)
    run_dir = run_benchmark(
        config,
        input_override=args.input,
        output_override=args.output,
    )
    print(f"Benchmark complete: {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
