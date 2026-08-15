"""What the seed-merge fix changes in the real cluster products.

The fix (`6b382d9`) stopped a forced-catalogue star from claiming a detected
star's seed and deleting it, which cost one member of every close pair. It was
found and measured on implanted stars; the published cluster photometry was
measured before it and has never been remeasured. This compares the two runs
where it matters: how many stars come back, whether their magnitudes moved, and
whether the move is concentrated in the close pairs the fix was about.

`det_uid` is the same identity in both runs — the detection is shared, only
Step 8's fit differs — so stars pair exactly instead of by position.

    seedfix_product_delta.py --target M13 [--target M5 ...]

Reads `<result>/cmd_psf_PRE_SEEDFIX` against `<result>/cmd_psf`, and the two
zero-point directories if both are present.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PHASE3 = Path(r"E:\APEX_validation\phase3")
CLOSE_FWHM = 2.0          # the separation the fix was about


def _load(folder: Path, frame: str) -> pd.DataFrame | None:
    path = folder / f"photometry_{frame}"
    if not path.exists():
        return None
    table = pd.read_csv(path, sep="\t")
    for column in ("mag_psf", "neighbor_dist_fwhm", "flags_psf", "snr_psf"):
        if column in table:
            table[column] = pd.to_numeric(table[column], errors="coerce")
    return table


def _frame_delta(before: pd.DataFrame, after: pd.DataFrame) -> dict:
    """Paired on det_uid; production keeps only flags_psf == 0."""
    def clean(table):
        keep = np.isfinite(table["mag_psf"]) & (table["flags_psf"].fillna(-1) == 0)
        return table[keep].set_index("det_uid")

    a, b = clean(before), clean(after)
    shared = a.index.intersection(b.index)
    delta = (b.loc[shared, "mag_psf"] - a.loc[shared, "mag_psf"]).to_numpy(float)
    separation = b.loc[shared, "neighbor_dist_fwhm"].to_numpy(float)
    close = np.isfinite(separation) & (separation < CLOSE_FWHM)
    return {
        "n_before": len(a), "n_after": len(b),
        "gained": int(len(b.index.difference(a.index))),
        "lost": int(len(a.index.difference(b.index))),
        "n_paired": int(len(shared)),
        "median_delta_mmag": float(np.nanmedian(delta)) * 1000 if len(delta) else np.nan,
        "rms_delta_mmag": float(np.sqrt(np.nanmean(delta ** 2))) * 1000 if len(delta) else np.nan,
        "n_moved_over_10mmag": int(np.sum(np.abs(delta) > 0.010)),
        "n_close": int(close.sum()),
        "close_median_delta_mmag": (float(np.nanmedian(delta[close])) * 1000
                                    if close.sum() else np.nan),
        "close_rms_delta_mmag": (float(np.sqrt(np.nanmean(delta[close] ** 2))) * 1000
                                 if close.sum() else np.nan),
    }


def _zeropoint_delta(result: Path) -> dict | None:
    """Frame zero-points are the product the paper quotes; did they move?"""
    pairs = []
    for name in ("frame_zeropoint.csv", "median_by_ID_filter.csv"):
        old = result / "cmd_zeropoint_PRE_SEEDFIX" / name
        new = result / "cmd_zeropoint" / name
        if old.exists() and new.exists():
            pairs.append((name, pd.read_csv(old), pd.read_csv(new)))
    if not pairs:
        return None
    # Pair on identity, not row order: the star list itself changes slightly,
    # so a positional subtraction would compare two different stars.
    keys = {"frame_zeropoint.csv": (["file", "filter"], "zp_frame"),
            "median_by_ID_filter.csv": (["ID", "FILTER"], "mag_cal_med")}
    out = {}
    for name, before, after in pairs:
        index, column = keys[name]
        if column not in before.columns or column not in after.columns:
            out[name] = {"comparable": False, "reason": f"{column} 없음"}
            continue
        a = before.set_index(index)[column].pipe(pd.to_numeric, errors="coerce")
        b = after.set_index(index)[column].pipe(pd.to_numeric, errors="coerce")
        shared = a.index.intersection(b.index)
        difference = (b.loc[shared] - a.loc[shared]).to_numpy(float)
        out[name] = {"comparable": True, "column": column,
                     "rows": [len(a), len(b)], "paired": int(len(shared)),
                     "median_mmag": round(float(np.nanmedian(difference)) * 1000, 3),
                     "rms_mmag": round(float(np.sqrt(np.nanmean(difference ** 2))) * 1000, 2),
                     "max_abs_mmag": round(float(np.nanmax(np.abs(difference))) * 1000, 1)}
        if name == "frame_zeropoint.csv" and "zp_scatter" in before.columns:
            out[name]["zp_scatter_median"] = [
                round(float(pd.to_numeric(t["zp_scatter"], errors="coerce").median()) * 1000, 2)
                for t in (before, after)]
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", action="append", required=True)
    ap.add_argument("--root", type=Path, default=PHASE3)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    rows = []
    for target in args.target:
        result = args.root / target / "result"
        before_dir = result / "cmd_psf_PRE_SEEDFIX"
        after_dir = result / "cmd_psf"
        if not before_dir.is_dir():
            print(f"{target}: 수정 전 사본이 없다 — 건너뜀", flush=True)
            continue
        frames = sorted(p.name.removeprefix("photometry_")
                        for p in before_dir.glob("photometry_*.tsv"))
        print(f"\n=== {target} · {len(frames)} 프레임 ===")
        print(f'{"프레임":<34}{"전":>6}{"후":>6}{"증":>5}{"감":>5}'
              f'{"중앙차":>9}{"RMS":>8}{"10mmag+":>9}{"근접n":>7}{"근접RMS":>9}')
        print("-" * 98)
        for frame in frames:
            before, after = _load(before_dir, frame), _load(after_dir, frame)
            if before is None or after is None:
                print(f"{frame:<34}(재처리본 없음)")
                continue
            row = {"target": target, "frame": frame, **_frame_delta(before, after)}
            rows.append(row)
            print(f'{frame:<34}{row["n_before"]:>6}{row["n_after"]:>6}'
                  f'{row["gained"]:>5}{row["lost"]:>5}'
                  f'{row["median_delta_mmag"]:>+8.2f}m{row["rms_delta_mmag"]:>7.1f}m'
                  f'{row["n_moved_over_10mmag"]:>9}{row["n_close"]:>7}'
                  f'{row["close_rms_delta_mmag"]:>8.1f}m')
        done = [r for r in rows if r["target"] == target]
        if done:
            frame = pd.DataFrame(done)
            print(f'  합계: 별 {frame["n_before"].sum()} -> {frame["n_after"].sum()} '
                  f'(증 {frame["gained"].sum()} · 감 {frame["lost"].sum()}) · '
                  f'10 mmag 넘게 움직인 별 {frame["n_moved_over_10mmag"].sum()} / '
                  f'{frame["n_paired"].sum()}')
        zeropoint = _zeropoint_delta(result)
        if zeropoint:
            print(f"  영점: {json.dumps(zeropoint, ensure_ascii=False)}")

    if args.output and rows:
        pd.DataFrame(rows).to_csv(args.output, index=False)
        print(f"\n표 -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
