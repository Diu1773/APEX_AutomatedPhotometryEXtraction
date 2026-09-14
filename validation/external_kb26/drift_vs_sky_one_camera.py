"""카메라를 고정하고 하늘만 바꾼다 — 트랩 가설의 결정적 시험 (2026-09-14).

## 왜 이 설계인가

12 절에서 남은 물음은 하나다. 결손이 **하늘 배경을 따라 줄어드는가** — 그것이
트랩이 배경 전하로 채워진다(**fat zero**)는 가설의 핵심 예측이다. 지금까지 가진
증거는 **카메라가 다른 둘**(kb26 하늘 21.9 e- · kb27 하늘 89.0 e-)을 견준 것이라
카메라 차이와 하늘 차이가 안 갈렸고, 프레임 **안**에서 본 것은 하늘 폭이 2.56 e-
밖에 안 되는 데다 혼잡과 상관 −0.64 로 얽혀 있었다.

**그래서 카메라를 고정하고 하늘만 바꾼다.** 같은 기기가 **같은 시야**를 여러 밤에
찍은 것을 모으면 달 위상과 고도가 밤마다 달라 하늘이 몇 배씩 변한다.

이 설계의 값어치는 **기준 오차가 통째로 상쇄된다**는 데 있다. 기준 등급이 틀렸더라도
그 틀림은 프레임마다 **똑같이** 들어가므로, 어두운 밤과 밝은 밤의 **치우침 차이**에는
안 남는다. 변환식·색항·성단의 색-등급 상관 같은 것을 하나도 안 믿어도 된다.

    하늘이 밝은 밤에 치우침이 작으면   트랩이 채워진다 = fat zero 가 맞다
    하늘과 무관하면                     배경 의존이 아니라 다른 것이다

## 무엇을 쓰나

**BANZAI 가 프레임마다 붙여 놓은 `CAT` 측광을 그대로 쓴다.** APEX 를 돌릴 필요가
없고, 층을 갈아 끼워도 폭이 같다는 것은 8 절에서 이미 확인했다. 기준은 그 시야의
Gaia DR3 를 한 번 받아 캐시한다.

실행:
    python -X utf8 validation/external_kb26/drift_vs_sky_one_camera.py \\
        --frames <BANZAI 프레임 폴더> --label kb26
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits

MATCH_ARCSEC = 1.5
MIN_CELL = 8
MIN_FIT = 60
#: Gaia G 로 나눈 등급 구간. 기기 등급이 아니라 기준 등급을 축으로 쓴다.
MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0]


def gaia_for_field(ra: float, dec: float, radius_deg: float, cache: Path) -> pd.DataFrame:
    """그 시야의 Gaia DR3. 한 번 받아 두고 다음부터는 파일에서 읽는다."""
    if cache.exists():
        out = pd.read_csv(cache)
        print(f"Gaia {len(out)} 줄 (저장해 둔 것: {cache.name})")
        return out
    from astroquery.vizier import Vizier
    import astropy.units as u

    v = Vizier(columns=["RA_ICRS", "DE_ICRS", "Gmag", "BPmag", "RPmag", "BP-RP"],
               column_filters={"Gmag": "<19"}, row_limit=-1)
    res = v.query_region(SkyCoord(ra, dec, unit="deg"),
                         radius=radius_deg * u.deg, catalog="I/355/gaiadr3")
    if not res:
        raise SystemExit("Gaia 결과가 없다")
    t = res[0].to_pandas().dropna(subset=["Gmag"])
    t = t.rename(columns={"RA_ICRS": "ra_deg", "DE_ICRS": "dec_deg",
                          "Gmag": "G", "BPmag": "BP", "RPmag": "RP",
                          "BP-RP": "BP_RP"})
    cache.parent.mkdir(parents=True, exist_ok=True)
    t.to_csv(cache, index=False)
    print(f"Gaia {len(t)} 줄 (새로 받아 {cache.name} 에 저장)")
    return t


def _robust_line(x: np.ndarray, y: np.ndarray, iters: int = 5):
    """3σ 로 다듬는 직선 맞춤. 색항을 빼는 데만 쓴다."""
    keep = np.ones(len(x), bool)
    a = b = float("nan")
    for _ in range(iters):
        if int(keep.sum()) < 10:
            break
        b, a = np.polyfit(x[keep], y[keep], 1)
        r = y - (a + b * x)
        s = 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep])))
        if not np.isfinite(s) or s <= 0:
            break
        keep = np.abs(r - np.median(r[keep])) <= 3.0 * s
    return a, b


def drift_of_frame(path: Path, gaia: SkyCoord, gmag: np.ndarray,
                   colour: np.ndarray, column: str = "fluxaper3") -> dict | None:
    """프레임 한 장의 밝기 치우침과 그 프레임의 하늘·달 정보.

    **기본이 고정 조리개(`fluxaper3`)인 이유가 있다.** Kron 조리개는 별마다
    프레임마다 크기가 달라서, 시상이 나쁜 밤에는 조리개가 커지고 그만큼 하늘을 더
    담는다. 그러면 「하늘이 밝아서 덜하다」와 「조리개가 작아서 덜하다」가 안
    갈린다. 고정 반지름이면 그 변인이 사라진다.
    """
    with fits.open(path, memmap=False) as hdul:
        if "CAT" not in [h.name for h in hdul]:
            return None
        head, cat = hdul["SCI"].header, hdul["CAT"].data
        exptime = float(head.get("EXPTIME") or 1.0)
        sky = float(head.get("L1MEDIAN") or np.nan)
        fwhm = float(head.get("L1FWHM") or np.nan)
        moon_frac = float(head.get("MOONFRAC") or np.nan)
        moon_alt = float(head.get("MOONALT") or np.nan)
        filt = str(head.get("FILTER") or "?")
        night = str(head.get("DAY-OBS") or head.get("DATE-OBS") or "?")[:10]
    flag = np.asarray(cat["flag"], int)
    if column not in cat.dtype.names:
        return None
    flux = np.asarray(cat[column], float)
    keep = (flag == 0) & np.isfinite(flux) & (flux > 0)
    if int(keep.sum()) < MIN_FIT:
        return None
    cb = SkyCoord(np.asarray(cat["ra"], float)[keep],
                  np.asarray(cat["dec"], float)[keep], unit="deg")
    k, sep, _ = gaia.match_to_catalog_sky(cb)
    ok = np.asarray(sep.arcsec) <= MATCH_ARCSEC
    inst = np.full(len(gmag), np.nan)
    inst[ok] = -2.5 * np.log10(flux[keep][np.asarray(k)[ok]] / exptime)

    delta = gmag - inst
    g = np.isfinite(delta) & np.isfinite(colour) & np.isfinite(gmag)
    if int(g.sum()) < MIN_FIT:
        return None
    # **색항은 프레임마다 다시 뺀다.** 여기서 관심은 색이 아니라 밝기 축이다.
    a, b = _robust_line(colour[g], delta[g])
    resid = delta[g] - (a + b * colour[g])
    resid -= np.median(resid)
    mag = gmag[g]
    meds, cent, ns = [], [], []
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        m = (mag >= lo) & (mag < hi)
        if int(m.sum()) >= MIN_CELL:
            meds.append(float(np.median(resid[m])))
            cent.append(0.5 * (lo + hi))
            ns.append(int(m.sum()))
    if len(meds) < 3:
        return None
    slope = float(np.polyfit(cent, meds, 1)[0])
    return dict(frame=path.name, night=night, filt=filt, exptime=exptime,
                column=column,
                sky=sky, fwhm=fwhm, moon_frac=moon_frac, moon_alt=moon_alt,
                n=int(g.sum()), slope=slope,
                span=float(max(meds) - min(meds)),
                bin_medians=meds, bin_counts=ns)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="카메라를 고정하고 하늘만 바꿔 치우침을 잰다")
    ap.add_argument("--frames", type=Path, required=True,
                    help="BANZAI 프레임(.fits.fz)이 든 폴더. 하위 폴더도 훑는다")
    ap.add_argument("--label", default="kb26")
    ap.add_argument("--cache", type=Path, default=None,
                    help="Gaia 캐시 파일 (기본: frames 폴더 옆)")
    ap.add_argument("--column", default="fluxaper3",
                    help="밝기 열. 기본은 고정 조리개 — Kron(`flux`)은 시상을 탄다")
    a = ap.parse_args(argv)

    paths = sorted(a.frames.rglob("*.fits.fz"))
    if not paths:
        raise SystemExit(f"프레임이 없다: {a.frames}")
    print(f"{a.label} · 프레임 {len(paths)} 장 · 밝기 열 {a.column}")

    with fits.open(paths[0], memmap=False) as hdul:
        h = hdul["SCI"].header
        ra0, dec0 = float(h["CRVAL1"]), float(h["CRVAL2"])
    cache = a.cache or (a.frames / f"gaia_{a.label}.csv")
    gaia_tab = gaia_for_field(ra0, dec0, 0.35, cache)
    gaia = SkyCoord(gaia_tab["ra_deg"].to_numpy(float),
                    gaia_tab["dec_deg"].to_numpy(float), unit="deg")
    gmag = gaia_tab["G"].to_numpy(float)
    colour = pd.to_numeric(gaia_tab.get("BP_RP"), errors="coerce").to_numpy(float)

    rows = []
    for p in paths:
        r = drift_of_frame(p, gaia, gmag, colour, a.column)
        if r:
            rows.append(r)
    if not rows:
        raise SystemExit("잰 프레임이 하나도 없다")

    print()
    print("=== 카메라 고정 · 하늘만 바꿈 ===")
    print("기준이 틀려도 프레임마다 똑같이 틀리므로 **밤 사이의 차이**에는 안 남는다.")
    print()
    print(f"{'밤':<12}{'필터':<5}{'노출':>7}{'하늘(e-)':>10}{'달':>7}{'고도':>7}"
          f"{'별수':>7}{'기울기':>9}{'폭':>8}")
    print("-" * 74)
    for r in sorted(rows, key=lambda d: (d["filt"], d["sky"])):
        print(f"{r['night']:<12}{r['filt']:<5}{r['exptime']:>7.1f}{r['sky']:>10.1f}"
              f"{r['moon_frac']:>7.2f}{r['moon_alt']:>7.0f}{r['n']:>7}"
              f"{r['slope']:>9.4f}{r['span']:>8.4f}")

    print()
    for filt in sorted({r["filt"] for r in rows}):
        sel = [r for r in rows if r["filt"] == filt and np.isfinite(r["sky"])
               and r["sky"] > 0]
        if len(sel) < 4:
            continue
        x = np.log(np.array([r["sky"] for r in sel]))
        y = np.array([r["slope"] for r in sel])
        w = np.array([r["fwhm"] for r in sel], float)
        rho = float(np.corrcoef(x, y)[0, 1])
        b, _a = np.polyfit(x, y, 1)
        line = (f"[{filt}] 프레임 {len(sel)} 장 · 하늘 {np.exp(x.min()):.0f}~"
                f"{np.exp(x.max()):.0f} e- · 상관 {rho:+.2f} · "
                f"하늘 e 배마다 {b:+.4f}")
        # **시상을 같이 넣어야 한다** — 밤마다 시상이 달라지고, 시상이 좋으면
        # 별빛이 더 모여서 하늘 몫이 줄기 때문에 하늘과 같은 방향으로 움직인다.
        if len(sel) >= 6 and np.isfinite(w).all() and np.ptp(w) > 0:
            A = np.vstack([np.ones(len(sel)), x, w]).T
            beta, *_ = np.linalg.lstsq(A, y, rcond=None)
            # 부트스트랩으로 오차
            rng = np.random.default_rng(20260914)
            boot = np.empty((400, 3))
            for i in range(400):
                k = rng.integers(0, len(sel), len(sel))
                bb, *_ = np.linalg.lstsq(A[k], y[k], rcond=None)
                boot[i] = bb
            e = boot.std(axis=0)
            line += (f"\n     시상을 같이 넣으면  하늘 {beta[1]:+.4f} ± {e[1]:.4f}"
                     f"  ·  FWHM {beta[2]:+.4f} ± {e[2]:.4f}"
                     f"  (시상 {w.min():.2f}~{w.max():.2f}″)")
        print(line)
    print()
    print("  하늘이 밝을수록 기울기가 0 에 가까워지면(상관 +) fat zero 가 맞다.")
    print("  **시상을 같이 넣은 뒤에도 하늘 항이 남아야 한다** — 시상이 좋은 밤은")
    print("  하늘도 어두운 경향이 있어서, 안 넣으면 둘이 서로를 흉내 낸다.")

    out = a.frames / f"drift_vs_sky_{a.label}.json"
    out.write_text(json.dumps(rows, ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print(f"\n썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
