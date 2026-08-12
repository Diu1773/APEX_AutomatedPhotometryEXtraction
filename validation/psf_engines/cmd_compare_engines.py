"""M13 CMD from APEX step 8 and from IRAF DAOPHOT ALLSTAR, side by side.

The same two real cluster frames (0001-B, 0001-V) go through each engine, and
each engine's stars are joined to the step-7 catalogue by position so B and V
meet on the same physical star. Both panels are instrumental — no zero point,
colour term, or aperture correction — so only the shape of the sequences is
comparable, not the absolute colour or magnitude. Both engines used a zeropoint
of 25, so the two panels sit on close (not identical) instrumental scales.

The point is not which CMD is "right" — neither is calibrated — but whether the
two engines draw the same cluster, and how the trade-off measured on artificial
stars (APEX keeps more, ALLSTAR scatters less) shows up in a real diagram.
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


def _match_source_id(x: np.ndarray, y: np.ndarray, mag: np.ndarray,
                     err: np.ndarray, step7_tsv: Path,
                     radius_px: float) -> pd.DataFrame:
    """Attach step-7 source_id to each measurement by nearest position.

    Beyond the match radius the star gets source_id -1 and is dropped. One row
    per star is kept, the one with the smaller reported error.
    """
    s7 = pd.read_csv(step7_tsv, sep="\t")
    s7 = s7[["source_id", "x", "y"]].apply(pd.to_numeric, errors="coerce").dropna()
    dist, idx = cKDTree(s7[["x", "y"]].to_numpy(float)).query(
        np.column_stack([x, y]), k=1)
    out = pd.DataFrame({
        "source_id": np.where(dist <= radius_px, s7["source_id"].to_numpy()[idx], -1),
        "mag": mag, "err": err,
    })
    out = out[out["source_id"] > 0]
    return out.sort_values("err").drop_duplicates("source_id", keep="first")


def load_allstar(csv: Path, step7: Path, radius: float) -> pd.DataFrame:
    a = pd.read_csv(csv)
    a["mag"] = pd.to_numeric(a["mag"], errors="coerce")
    a["merr"] = pd.to_numeric(a["merr"], errors="coerce")
    a = a[np.isfinite(a["mag"])]
    return _match_source_id(a["x"].to_numpy(float), a["y"].to_numpy(float),
                            a["mag"].to_numpy(float), a["merr"].to_numpy(float),
                            step7, radius)


def load_apex(tsv: Path, step7: Path, radius: float) -> pd.DataFrame:
    d = pd.read_csv(tsv, sep="\t")
    for c in ("x_fit", "y_fit", "mag_psf", "mag_psf_err", "flags_psf"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    # Clean fits only — the same gate the engine comparison used.
    d = d[(d["flags_psf"] == 0) & np.isfinite(d["mag_psf"])]
    return _match_source_id(d["x_fit"].to_numpy(float), d["y_fit"].to_numpy(float),
                            d["mag_psf"].to_numpy(float),
                            d["mag_psf_err"].to_numpy(float), step7, radius)


def cmd_points(b: pd.DataFrame, v: pd.DataFrame) -> pd.DataFrame:
    m = b[["source_id", "mag"]].merge(
        v[["source_id", "mag"]], on="source_id", suffixes=("_B", "_V"))
    m["BV"] = m["mag_B"] - m["mag_V"]
    return m


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apex-b", required=True)
    ap.add_argument("--apex-v", required=True)
    ap.add_argument("--allstar-b", required=True)
    ap.add_argument("--allstar-v", required=True)
    ap.add_argument("--step7-b", required=True)
    ap.add_argument("--step7-v", required=True)
    ap.add_argument("--match-radius-px", type=float, default=2.0)
    ap.add_argument("--outdir", default=str(REPO / "validation" / "psf_engines"))
    args = ap.parse_args()
    r = args.match_radius_px

    apex = cmd_points(load_apex(Path(args.apex_b), Path(args.step7_b), r),
                      load_apex(Path(args.apex_v), Path(args.step7_v), r))
    allstar = cmd_points(load_allstar(Path(args.allstar_b), Path(args.step7_b), r),
                         load_allstar(Path(args.allstar_v), Path(args.step7_v), r))
    print(f"CMD 점: APEX {len(apex)}개 · ALLSTAR {len(allstar)}개")

    if _STYLE:
        apply_paper_style()
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 5.2), sharey=True)
    for ax, (name, m) in zip(axes, (("APEX step 8 (empirical ePSF)", apex),
                                     ("IRAF DAOPHOT ALLSTAR (Moffat)", allstar))):
        ax.plot(m["BV"], m["mag_V"], ".", ms=2.0, mew=0, alpha=0.55,
                color=PALETTE["black"])
        ax.set_xlabel("$B-V$ (instrumental)")
        ax.set_title(f"{name}\n{len(m)} stars", fontsize=8.5)
        ax.set_xlim(-0.2, 2.2)
    axes[0].set_ylabel("$V$ (instrumental)")
    axes[0].invert_yaxis()
    fig.text(0.5, 0.005,
             "Moravian C3-61000 · M13 · 2026-05-15 · pp_messier13-0001-B/V · 60 s · "
             "instrumental mag (no zeropoint/colour term)",
             ha="center", fontsize=6, color=PALETTE["grey"])
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    outdir = Path(args.outdir)
    if _STYLE:
        paths = save_fig(fig, "cmd_compare_engines_m13", outdir)
        for k, p in paths.items():
            print(f"[{k}] {p}")
    else:
        p = outdir / "cmd_compare_engines_m13.png"
        fig.savefig(p, dpi=150)
        print(f"[png] {p}")

    for name, m in (("apex", apex), ("allstar", allstar)):
        m.to_csv(outdir / f"cmd_compare_{name}_points.csv", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
