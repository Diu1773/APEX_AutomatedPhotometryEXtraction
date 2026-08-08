"""B3: what the light-curve build actually reads (docs/audit/APEX_PERF_DEV_PLAN.md).

O3 proposes a shared cache, a frame x star matrix and a bounded preload. Before
any of that, the numbers it is supposed to improve have to exist: how many
frame tables get parsed, how many rows that is, and what the working set costs.

Read *counts* are deterministic, so unlike wall time they are meaningful even
on a machine that is busy with something else — which is why this can run
before the idle-machine measurements.

Two paths are profiled because they are shaped differently:

* ``service`` — ``load_filter_photometry_timeseries``, the Qt-free path. One
  pass over the frames, producing a long frame/star table.
* ``per_star`` — the shape the GUI builder uses: ``_build_star_mag_series`` is
  called once per star and walks every frame each time. With a warm cache that
  is free; the point of measuring is to show what it costs when the cache is
  dropped, which is what happens in the builder whenever ``result_dir``
  changes (``_get_photometry_df`` clears the whole cache, keyed on one
  directory at a time).

Usage:
    python -X utf8 benchmark/b3_lc_load.py [--result-dir DIR] [--stars 10]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from apex.benchmark import resources  # noqa: E402
from apex.utils import photometry_loader as PL  # noqa: E402

DEFAULT_RESULT = Path(r"E:\APEX_validation\reprocess\YZBoo_2n\result")
OUT = REPO / "benchmark" / "perf" / "20260809" / "b3_lc_load.json"


def profile_service(result_dir: Path, filt: str) -> dict:
    """One pass over every frame of a filter (the Qt-free service path)."""
    from apex.analysis.light_curve.photometry_source_service import (
        load_filter_photometry_timeseries,
    )

    PL.reset_load_counters()
    with resources.measure(f"service_{filt}") as m:
        table, source = load_filter_photometry_timeseries(result_dir, filt)
    m.update(PL.get_load_counters())
    m["filter"] = filt
    m["path"] = "service"
    m["source"] = source.get("source")
    m["table_rows"] = int(len(table))
    m["stars"] = int(table["star_id"].nunique()) if "star_id" in table else 0
    m["frames_in_table"] = (int(table["frame_order"].nunique())
                            if "frame_order" in table else 0)
    return m


def profile_per_star(result_dir: Path, filt: str, n_stars: int) -> dict:
    """The builder's shape: walk every frame once per star, cache dropped.

    This is not a strawman — it is what the GUI builder degrades to whenever
    the frame it wants belongs to a different ``result_dir`` than the one the
    cache currently holds, because the cache keeps a single directory.
    """
    from apex.analysis.light_curve.photometry_source_service import (
        load_lightcurve_frame_photometry,
        resolve_lightcurve_photometry_source,
    )
    import pandas as pd

    source = resolve_lightcurve_photometry_source(result_dir, None)
    index_path = result_dir / "step7_forced_phot" / "photometry_index.csv"
    index = pd.read_csv(index_path)
    column = "filter" if "filter" in index.columns else "FILTER"
    from apex.utils.common_helpers import normalize_filter_key

    wanted = normalize_filter_key(filt)
    files = [Path(str(v)).name for v, f in zip(index["file"], index[column])
             if normalize_filter_key(f) == wanted]

    PL.reset_load_counters()
    with resources.measure(f"per_star_{filt}") as m:
        for _ in range(n_stars):
            for fname in files:
                load_lightcurve_frame_photometry(result_dir, fname, source)
    m.update(PL.get_load_counters())
    m["filter"] = filt
    m["path"] = "per_star_uncached"
    m["stars"] = n_stars
    m["frames"] = len(files)
    return m


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result-dir", default=str(DEFAULT_RESULT))
    ap.add_argument("--filters", default="g,r,i")
    ap.add_argument("--stars", type=int, default=5,
                    help="stars for the uncached per-star profile")
    args = ap.parse_args()

    resources.warn_if_busy()
    result_dir = Path(args.result_dir)
    records = []

    print(f"{'path':18s}{'filt':>5}{'frames':>8}{'rows':>10}"
          f"{'stars':>7}{'wall':>9}{'USS':>9}")
    for filt in args.filters.split(","):
        m = profile_service(result_dir, filt)
        records.append(m)
        print(f"{m['path']:18s}{filt:>5}{m['frames_loaded']:>8}"
              f"{m['rows_loaded']:>10,}{m['stars']:>7}"
              f"{m['wall_s']:>9.1f}{m['peak_uss_mb'] or 0:>9,.0f}")

    for filt in args.filters.split(","):
        m = profile_per_star(result_dir, filt, args.stars)
        records.append(m)
        print(f"{m['path']:18s}{filt:>5}{m['frames_loaded']:>8}"
              f"{m['rows_loaded']:>10,}{m['stars']:>7}"
              f"{m['wall_s']:>9.1f}{m['peak_uss_mb'] or 0:>9,.0f}")

    resources.save_metrics({"label": "b3_lc_load", "result_dir": str(result_dir),
                            "runs": records}, OUT)
    print(f"\nsaved -> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
