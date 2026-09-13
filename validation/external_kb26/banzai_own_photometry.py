"""kb26 의 밝기 치우침이 APEX 의 측광 탓인가 — 관측소 자신의 측광으로 잰다.

**세 갈래로 가른다** (`../ERROR_BUDGET.md` 7 절). kb26 은 밝은 별을 더 밝게,
어두운 별을 더 어둡게 내고 그 폭이 0.28~0.38 등급이다. Moravian 과 MuSCAT3 는
0.03 안쪽이다. 원인이 어느 층에 있는지 가르려면 층을 하나씩 바꿔 봐야 한다.

    APEX 가 원본을 보정하고 APEX 가 측광    ← 이미 잰 것. 폭 0.276 (r)
    BANZAI 가 보정하고 APEX 가 측광         ← 다시 돌린다. 보정 층만 바뀐다
    BANZAI 가 보정하고 **BANZAI 가 측광**   ← 이 스크립트. 측광 층까지 바뀐다

셋째가 이 스크립트이고 **재처리가 필요 없다.** BANZAI 프레임은 자기 천체 목록을
`CAT` 확장에 함께 담아 내놓기 때문이다 — 좌표와 밝기 값이 이미 들어 있다.

## 재는 법

프레임마다 `CAT` 을 읽어 밝기 값을 기기 등급으로 바꾸고(−2.5 log10), 별마다 여섯
장의 중앙값을 쓴다. 그것을 **APEX 가 쓴 것과 같은 기준 등급**(`ref_<밴드>`, Gaia 를
표준계로 옮긴 값)에 견주어 영점과 색항을 맞추고, 남은 잔차를 밝기 구간으로 나눈다.

    잔차 = 기준 등급 − (기기 등급 + zp + ct·색)

**적합기는 APEX 의 것을 그대로 빌린다.** 다시 구현하면 「APEX 의 측광」이 아니라
「내 재구현」을 재게 된다. 절단도 같은 방식으로 돈다.

**색은 APEX 의 보정성 표에서 가져온다.** BANZAI 목록에는 색이 없고, 색을 어디서
가져오든 그것은 두 측광에 똑같이 실리므로 **둘의 차이**에는 영향을 주지 않는다.

실행:
    python -X utf8 validation/external_kb26/banzai_own_photometry.py
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

from apex.analysis.cmd.zeropoint_runner import (  # noqa: E402
    ZeropointCalibrationRunner,
)
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

BANZAI_ROOT = Path("E:/APEX_validation_output/external_kb26/raw/banzai")
APEX_KB26 = REPO / "validation/external_kb26/results"
#: BANZAI 의 밴드 폴더 이름 → APEX 가 쓰는 필터 이름
BAND_MAP = {"B": "B", "V": "V", "rp": "r", "ip": "i"}
MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0]
MATCH_ARCSEC = 1.5
MIN_CELL = 8


class _Fitter:
    """APEX 의 적합기를 빌려 쓰기 위한 껍데기.

    `_robust_linfit` 은 러너의 메서드인데 `self` 에서 쓰는 것이 로그뿐이다.
    """

    def _log(self, *a, **k):
        pass

    fit = ZeropointCalibrationRunner._robust_linfit


FIT = _Fitter()


def mad(v) -> float:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if v.size else np.nan


def banzai_instrumental(band_dir: Path) -> pd.DataFrame:
    """그 밴드의 여섯 장에서 별마다의 기기 등급 중앙값.

    별은 하늘 좌표로 이어 붙인다. 첫 장을 기준 목록으로 삼고 나머지를 1.5 초각
    안에서 맞춘다 — BANZAI 의 측성 해가 붙어 있으므로 좌표로 이을 수 있다.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astropy.io import fits

    frames = sorted(band_dir.glob("*.fits.fz"))
    if not frames:
        return pd.DataFrame()

    base: pd.DataFrame | None = None
    cols: list[np.ndarray] = []
    for p in frames:
        with fits.open(p, memmap=False) as hdul:
            if "CAT" not in hdul:
                continue
            d = hdul["CAT"].data
            exptime = float(hdul["SCI"].header.get("EXPTIME") or 1.0)
        t = pd.DataFrame({
            "ra_deg": np.asarray(d["ra"], float),
            "dec_deg": np.asarray(d["dec"], float),
            "flux": np.asarray(d["flux"], float),
            "flag": np.asarray(d["flag"], int),
        })
        # 깃발이 선 별(포화·이웃과 겹침·잘림)은 버린다. BANZAI 자신의 판정이다.
        t = t[(t["flag"] == 0) & np.isfinite(t["flux"]) & (t["flux"] > 0)]
        if t.empty:
            continue
        # 노출로 나눠 초당 밝기로 맞춘다 — 노출이 장마다 조금씩 다르다.
        t["mag_inst"] = -2.5 * np.log10(t["flux"].to_numpy(float) / exptime)
        if base is None:
            base = t[["ra_deg", "dec_deg"]].reset_index(drop=True)
            cols.append(t["mag_inst"].to_numpy(float))
            continue
        cb = SkyCoord(base["ra_deg"].to_numpy(float),
                      base["dec_deg"].to_numpy(float), unit="deg")
        ct = SkyCoord(t["ra_deg"].to_numpy(float),
                      t["dec_deg"].to_numpy(float), unit="deg")
        k, sep, _ = cb.match_to_catalog_sky(ct)
        v = np.full(len(base), np.nan)
        ok = sep.arcsec <= MATCH_ARCSEC
        v[ok] = t["mag_inst"].to_numpy(float)[np.asarray(k)[ok]]
        cols.append(v)

    if base is None or not cols:
        return pd.DataFrame()
    arr = np.vstack(cols)
    with np.errstate(invalid="ignore"):
        base["mag_inst"] = np.nanmedian(arr, axis=0)
        base["n_frames"] = np.sum(np.isfinite(arr), axis=0)
    return base[np.isfinite(base["mag_inst"]) & (base["n_frames"] >= 3)]


def main() -> int:
    from astropy.coordinates import SkyCoord

    cal_p = APEX_KB26 / "cmd_zeropoint/gaia_sdss_calibrator_by_ID.csv"
    co_p = APEX_KB26 / "cmd_zeropoint/zp_fit_coefficients.csv"
    if not (cal_p.exists() and co_p.exists()):
        print(f"APEX 의 kb26 결과가 없다: {cal_p}")
        return 1
    cal = pd.read_csv(cal_p)
    co = pd.read_csv(co_p).set_index("filter")
    try:
        qual = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
    except Exception:  # noqa: BLE001
        qual = np.ones(len(cal), bool)

    out: list[dict] = []
    print()
    print("=== BANZAI 자신의 측광으로 잰 밝기 치우침 (kb26 M67) ===")
    print("APEX 의 측광을 빼고 잰다. 폭이 그대로면 측광 탓이 아니다.")
    print()
    hdr = "{:<6}{:>6}{:>9}".format("밴드", "별수", "잔차 MAD")
    hdr += "   " + " ".join(f"{lo:g}~{hi:g}".rjust(8)
                            for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]))
    hdr += "{:>9}{:>13}".format("폭", "APEX 의 폭")
    print(hdr)
    print("-" * len(hdr))

    #: 7 절이 APEX 측광으로 잰 폭. 이 스크립트가 견줄 상대다.
    apex_span = {"B": 0.376, "V": 0.320, "i": 0.278, "r": 0.276}

    for folder, band in BAND_MAP.items():
        bd = BANZAI_ROOT / folder
        if not bd.is_dir():
            continue
        bz = banzai_instrumental(bd)
        if bz.empty or band not in co.index:
            print(f"{band:<6} BANZAI 목록을 못 읽었다")
            continue

        colname = str(co.loc[band, "color_col"])
        cc, rc = f"color_{colname}", f"ref_{band}"
        if cc not in cal.columns or rc not in cal.columns:
            continue
        sub = cal[qual & np.isfinite(pd.to_numeric(cal[rc], errors="coerce"))
                  & np.isfinite(pd.to_numeric(cal[cc], errors="coerce"))]
        if len(sub) < 40:
            continue

        # APEX 의 보정성 별을 BANZAI 목록에 좌표로 잇는다.
        ca = SkyCoord(sub["ra_deg"].to_numpy(float),
                      sub["dec_deg"].to_numpy(float), unit="deg")
        cbz = SkyCoord(bz["ra_deg"].to_numpy(float),
                       bz["dec_deg"].to_numpy(float), unit="deg")
        k, sep, _ = ca.match_to_catalog_sky(cbz)
        ok = np.asarray(sep.arcsec) <= MATCH_ARCSEC
        if ok.sum() < 40:
            print(f"{band:<6} 이어진 별이 {int(ok.sum())} 개뿐이다")
            continue

        ref = pd.to_numeric(sub[rc], errors="coerce").to_numpy(float)[ok]
        col = pd.to_numeric(sub[cc], errors="coerce").to_numpy(float)[ok]
        inst = bz["mag_inst"].to_numpy(float)[np.asarray(k)[ok]]
        delta = ref - inst
        good = np.isfinite(delta) & np.isfinite(col) & np.isfinite(ref)
        if good.sum() < 40:
            continue

        # APEX 의 적합기로 영점과 색항을 맞춘다 — 절단까지 같은 방식이다.
        w = np.ones(int(good.sum()))
        zp, ct, n_in, _ = FIT.fit(col[good], delta[good], w=w,
                                  clip_sigma=3.0, iters=5,
                                  slope_absmax=0.8, min_n=10)
        if not np.isfinite(zp):
            print(f"{band:<6} 적합이 안 됐다")
            continue
        resid = delta[good] - (zp + ct * col[good])
        resid = resid - np.median(resid)
        mg = ref[good]

        meds, cells = [], []
        for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
            m = (mg >= lo) & (mg < hi)
            if m.sum() < MIN_CELL:
                cells.append("·".rjust(8))
                continue
            v = float(np.median(resid[m]))
            meds.append(v)
            cells.append(f"{v:+.3f}".rjust(8))
        span = float(max(meds) - min(meds)) if len(meds) >= 2 else np.nan
        print("{:<6}{:>6}{:>9.4f}   {}{:>9.4f}{:>13.3f}".format(
            band, int(good.sum()), mad(resid), " ".join(cells), span,
            apex_span.get(band, np.nan)))
        out.append(dict(band=band, n=int(good.sum()), zp=float(zp), ct=float(ct),
                        n_fit=int(n_in), resid_mad=mad(resid), span=span,
                        apex_span=apex_span.get(band)))

    print()
    print("  폭이 APEX 의 것과 비슷하면 **측광 층의 일이 아니다** — 관측소의 측광으로")
    print("  재도 같으므로 픽셀이나 검출기 쪽이다. 훨씬 작으면 APEX 의 측광이 원인이다.")

    p = REPO / "validation/external_kb26/banzai_own_photometry.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=2, default=str),
                 encoding="utf-8")
    print()
    print(f"썼다: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
