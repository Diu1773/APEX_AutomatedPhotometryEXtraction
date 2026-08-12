"""Draw M13's colour-magnitude diagram from IRAF DAOPHOT ALLSTAR photometry.

The engine comparison measured ALLSTAR on one injected B frame; this instead
runs it on the real cluster in two bands and plots the result, to see whether
the automated DAOPHOT chain produces a cluster CMD on its own.

The magnitudes are instrumental — no zero point, no colour term, no aperture
correction. Both frames used zmag = 25, so B and V share one offset and the
colour B-V is on a consistent scale, but it is displaced from the standard
system by whatever the filter/detector zero points are. The shape of the
sequences (main sequence, red-giant branch, horizontal branch) survives that
displacement; the absolute colour and magnitude do not. This is a "does the
engine draw a cluster" check, not a calibrated CMD.

Each ALLSTAR row carries its fitted (x, y). The star's identity comes from
matching that position back to the step-7 catalogue, which carries the Gaia
source_id, so B and V are joined on the same physical star rather than on a
fragile row order.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

import matplotlib.pyplot as plt  # noqa: E402

try:
    from apex_paper_style import PALETTE, apply_paper_style, save_fig
    _STYLE = True
except Exception:  # pragma: no cover
    _STYLE = False
    PALETTE = {"black": "#111", "grey": "#777", "blue": "#0072B2"}


def with_source_id(allstar_csv: Path, step7_tsv: Path,
                   radius_px: float) -> pd.DataFrame:
    """Attach the step-7 source_id to each ALLSTAR row by position."""
    a = pd.read_csv(allstar_csv)
    a["mag"] = pd.to_numeric(a["mag"], errors="coerce")
    a = a[np.isfinite(a["mag"])].copy()

    s7 = pd.read_csv(step7_tsv, sep="\t")
    s7 = s7[["source_id", "x", "y"]].apply(pd.to_numeric, errors="coerce").dropna()

    dist, idx = cKDTree(s7[["x", "y"]].to_numpy(float)).query(
        a[["x", "y"]].to_numpy(float), k=1)
    a["source_id"] = np.where(dist <= radius_px,
                              s7["source_id"].to_numpy()[idx], -1)
    a = a[a["source_id"] > 0]
    # One ALLSTAR measurement per star; if two rows hit the same catalogue star,
    # keep the smaller error.
    return a.sort_values("merr").drop_duplicates("source_id", keep="first")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--allstar-b", required=True)
    ap.add_argument("--allstar-v", required=True)
    ap.add_argument("--step7-b", required=True)
    ap.add_argument("--step7-v", required=True)
    ap.add_argument("--match-radius-px", type=float, default=2.0)
    ap.add_argument("--outdir", default=str(REPO / "validation" / "psf_engines"))
    args = ap.parse_args()

    b = with_source_id(Path(args.allstar_b), Path(args.step7_b), args.match_radius_px)
    v = with_source_id(Path(args.allstar_v), Path(args.step7_v), args.match_radius_px)
    print(f"ALLSTAR 측광: B {len(b)}개 · V {len(v)}개")

    m = b[["source_id", "mag"]].merge(
        v[["source_id", "mag"]], on="source_id", suffixes=("_B", "_V"))
    m["BV"] = m["mag_B"] - m["mag_V"]
    print(f"B∩V 매칭(=CMD 점): {len(m)}개")

    if _STYLE:
        apply_paper_style()
    fig, ax = plt.subplots(figsize=(4.2, 5.0))
    ax.plot(m["BV"], m["mag_V"], ".", ms=2.2, mew=0, alpha=0.6,
            color=PALETTE["black"])
    ax.set_xlabel("$B-V$ (instrumental)")
    ax.set_ylabel("$V$ (instrumental)")
    ax.invert_yaxis()
    lo, hi = m["BV"].quantile([0.01, 0.99])
    ax.set_xlim(lo - 0.3, hi + 0.3)
    ax.set_title(f"M13 — IRAF DAOPHOT ALLSTAR (auto)\n{len(m)} stars", fontsize=9)
    ax.text(0.03, 0.02,
            "Moravian C3-61000 · M13 · 2026-05-15\n"
            "pp_messier13-0001-B/V · 60 s · ALLSTAR (unattended)\n"
            "instrumental mag — no zeropoint/colour term",
            transform=ax.transAxes, fontsize=6, va="bottom", color=PALETTE["grey"])
    fig.tight_layout()

    outdir = Path(args.outdir)
    if _STYLE:
        paths = save_fig(fig, "cmd_allstar_m13", outdir)
        for k, p in paths.items():
            print(f"[{k}] {p}")
    else:
        p = outdir / "cmd_allstar_m13.png"
        fig.savefig(p, dpi=150)
        print(f"[png] {p}")

    m.to_csv(outdir / "cmd_allstar_m13_points.csv", index=False)
    print(f"\nCMD 점 {len(m)}개 · B-V {m['BV'].min():.2f}..{m['BV'].max():.2f} · "
          f"V {m['mag_V'].min():.2f}..{m['mag_V'].max():.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
