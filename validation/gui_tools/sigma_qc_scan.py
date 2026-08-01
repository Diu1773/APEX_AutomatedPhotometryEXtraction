"""Step 4 에서 쓸 수 있는 검출 임계 QC — 음수 검출로 가짜를 센다.

배경을 뺀 영상의 부호를 뒤집으면 별은 사라지고 잡음 요동만 남는다. 거기서 같은
임계로 검출한 개수가 그 프레임의 **가짜 검출 추정치**다. 외부 카탈로그도 이론
모델도 필요 없으므로 WCS 해맞춤(Step 5) 전에 돌릴 수 있다.

    순도 추정 = (양수검출 - 음수검출) / 양수검출

Gaia 가 있으면 같이 실측해 이 추정이 맞는지 대조한다.

    python sigma_qc_scan.py <fits> [gaia.ecsv]

출력은 JSON 한 줄(stdout 마지막)로도 남겨 여러 대상 비교에 쓴다.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np
import sep
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.spatial import cKDTree

SIGMAS = (5.0, 4.0, 3.2, 2.8, 2.5, 2.2, 2.0, 1.8, 1.5, 1.2, 1.0)
MATCH_ARCSEC = 2.0
MINAREA = 3
DEBLEND_NTHRESH = 64
DEBLEND_CONT = 0.004


def _extract(img: np.ndarray, s: float, rms: float) -> int:
    try:
        return len(sep.extract(img, s, err=rms, minarea=MINAREA,
                               deblend_nthresh=DEBLEND_NTHRESH,
                               deblend_cont=DEBLEND_CONT))
    except Exception:
        return -1


def main() -> int:
    fits_path = sys.argv[1]
    gaia_path = sys.argv[2] if len(sys.argv) > 2 else None

    hdr = fits.getheader(fits_path)
    data = fits.getdata(fits_path).astype(np.float32)
    ny, nx = data.shape
    bkg = sep.Background(data, bw=61, bh=61, fw=3, fh=3)
    sub = data - bkg.back()
    neg = np.ascontiguousarray(-sub)
    rms = float(bkg.globalrms)

    tree = None
    n_truth = 0
    tol_px = 0.0
    if gaia_path and Path(gaia_path).exists():
        try:
            wcs = WCS(hdr).celestial
            cat = Table.read(gaia_path)
            gx, gy = wcs.world_to_pixel_values(np.asarray(cat["ra"], float),
                                               np.asarray(cat["dec"], float))
            inside = (gx > 0) & (gx < nx) & (gy > 0) & (gy < ny)
            n_truth = int(inside.sum())
            if n_truth:
                scale = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
                tol_px = MATCH_ARCSEC / scale
                tree = cKDTree(np.c_[gx[inside], gy[inside]])
        except Exception as exc:  # WCS 나 카탈로그가 없으면 음수 검출만 쓴다
            print(f"  (Gaia 대조 불가: {exc})")

    name = Path(fits_path).name
    print(f"\n=== {name}  {nx}x{ny}  exp={hdr.get('EXPTIME','?')}s  "
          f"filt={hdr.get('FILTER','?')}  sky rms={rms:.2f} ADU ===")
    head = f"{'sigma':>6} {'양수':>7} {'음수':>7} {'순도추정':>9}"
    if tree is not None:
        head += f" | {'Gaia가짜':>8} {'순도실측':>9}"
    print(head)

    rows = []
    for s in SIGMAS:
        n_pos = _extract(sub, s, rms)
        n_neg = _extract(neg, s, rms)
        if n_pos <= 0:
            continue
        pur_est = max(0.0, (n_pos - n_neg) / n_pos)
        line = f"{s:>6.1f} {n_pos:>7} {n_neg:>7} {pur_est:>9.3f}"
        rec = {"sigma": s, "n_pos": n_pos, "n_neg": n_neg, "purity_est": pur_est}
        if tree is not None:
            obj = sep.extract(sub, s, err=rms, minarea=MINAREA,
                              deblend_nthresh=DEBLEND_NTHRESH,
                              deblend_cont=DEBLEND_CONT)
            d, _ = tree.query(np.c_[np.asarray(obj["x"], float),
                                    np.asarray(obj["y"], float)],
                              distance_upper_bound=tol_px)
            n_match = int(np.isfinite(d).sum())
            fp_true = n_pos - n_match
            pur_true = n_match / n_pos
            line += f" | {fp_true:>8} {pur_true:>9.3f}"
            rec.update(fp_gaia=fp_true, purity_gaia=pur_true)
        rows.append(rec)
        print(line, flush=True)

    out = {"file": name, "nx": nx, "ny": ny, "rms": rms,
           "exptime": float(hdr.get("EXPTIME", 0) or 0),
           "filter": str(hdr.get("FILTER", "")), "n_gaia_fov": n_truth,
           "rows": rows}
    print("JSON " + json.dumps(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
