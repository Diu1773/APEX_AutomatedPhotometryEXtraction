# -*- coding: utf-8 -*-
"""LCO 교차기기 PSF 결과 그림 + 리포트.

카메라별 4패널: (a) 프레임 중앙부 (b) EPSF 모델 (c) PSF 차감 잔차 (d) PSF-구경 차이 vs 등급.
그림 안에 기기·대상·필터·노출 명세를 박는다 (데이터=점, 요약선=선).

실행: .venv-deploy\\Scripts\\python -X utf8 validation/psf_crossinstrument/report_psf_cross.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.visualization import ZScaleInterval

BASE = Path(r"E:\APEX_validation\psf_crossinstrument")

CAMS = {
    "sinistro": dict(
        frame="pp_ngc5985-0001-rp.fit",
        label="LCO 1 m elp / Sinistro fa16 (Fairchild CCD, 4-amp)",
        spec="NGC 5985 field · rp 120 s · 0.389\"/px · 2026-07-07 · archive.lco.global",
    ),
    "qhy600": dict(
        frame="pp_proximafield-0001-V.fit",
        label="LCO 0.4 m coj / sq36 = QHY600 (CMOS, single-amp)",
        spec="Proxima Cen field · V 20 s · 0.74\"/px · 2026-07-09 · archive.lco.global",
    ),
}

Z = ZScaleInterval()


def _cut(a: np.ndarray, half: int = 320) -> np.ndarray:
    cy, cx = np.array(a.shape) // 2
    return a[cy - half:cy + half, cx - half:cx + half]


def one_camera(cam: str, C: dict) -> dict | None:
    R = BASE / cam / "result"
    psf_tsv = R / "cmd_psf" / f"photometry_{C['frame']}.tsv"
    if not psf_tsv.exists():
        print(f"[{cam}] cmd_psf 없음 — 건너뜀")
        return None

    p = pd.read_csv(psf_tsv, sep="\t")
    a = pd.read_csv(R / "step7_forced_phot" / f"photometry_{C['frame']}.tsv", sep="\t")
    m = p.dropna(subset=["det_uid"]).merge(
        a.dropna(subset=["det_uid"]), on="det_uid", suffixes=("_p", "_a")
    )
    ok = m[
        np.isfinite(m.mag_psf) & np.isfinite(m.mag_inst)
        & (m.snr_psf > 20) & (m.snr > 20)
        & (m.crowding_unreliable_psf == 0) & (m.saturated_psf == 0)
    ]
    d = (ok.mag_psf - ok.mag_inst).to_numpy(float)
    med = float(np.median(d))
    mad = float(1.4826 * np.median(np.abs(d - med)))

    meta = json.loads(
        (R / "cmd_psf" / f"residual_meta_{C['frame']}.json").read_text(encoding="utf-8")
    )
    fitw = meta.get("fit_window", {})
    er = meta.get("epsf_reference", {})

    img = fits.getdata(BASE / cam / "sci" / C["frame"]).astype(float)
    resid = fits.getdata(R / "cmd_psf" / f"residual_{C['frame']}").astype(float)
    epsf_files = sorted((R / "cmd_psf").glob(f"epsf_model_*_{C['frame']}s"))
    epsf = fits.getdata(epsf_files[0]).astype(float) if epsf_files else None

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.4))
    v0, v1 = Z.get_limits(_cut(img))
    axes[0].imshow(_cut(img), origin="lower", cmap="gray", vmin=v0, vmax=v1)
    axes[0].set_title("frame (center 640 px)", fontsize=10)
    if epsf is not None:
        axes[1].imshow(np.log10(np.clip(epsf, 1e-6, None)), origin="lower", cmap="viridis")
        axes[1].set_title(
            f"EPSF model (log) | n_stars={er.get('n_selected', '?')}", fontsize=10
        )
    else:
        axes[1].text(0.5, 0.5, "no EPSF file", ha="center")
    axes[2].imshow(_cut(resid), origin="lower", cmap="gray", vmin=v0, vmax=v1)
    axes[2].set_title("PSF-subtracted residual (same stretch)", fontsize=10)

    mags = ok.mag_inst.to_numpy(float)
    axes[3].plot(mags, d, ".", ms=3, alpha=0.5, color="#1f77b4", label=f"stars SNR>20 (N={len(ok)})")
    axes[3].axhline(med, color="crimson", lw=1.2, label=f"median {med:+.3f}")
    axes[3].axhline(med + mad, color="crimson", lw=0.7, ls="--")
    axes[3].axhline(med - mad, color="crimson", lw=0.7, ls="--", label=f"±MAD {mad:.3f}")
    axes[3].set_xlabel("aperture mag_inst")
    axes[3].set_ylabel("PSF − aperture [mag]")
    axes[3].set_ylim(med - 6 * mad, med + 6 * mad)
    axes[3].legend(fontsize=8, loc="upper left")
    axes[3].set_title("PSF vs forced aperture", fontsize=10)
    for ax in axes[:3]:
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"{C['label']}  —  {C['spec']}", fontsize=11)
    fig.text(
        0.01, 0.01,
        f"APEX Step-8 headless | fit_window={fitw.get('shape_px', '?')}px "
        f"energy={fitw.get('energy_fraction', float('nan')):.3f} | "
        f"n_fit={len(p)} good={int((p.flags_psf == 0).sum())} | "
        f"Δmag(PSF−ap, SNR>20) median={med:+.3f} MAD={mad:.3f}",
        fontsize=8, color="0.35",
    )
    fig.tight_layout(rect=[0, 0.04, 1, 0.93])
    out = BASE / f"fig_psf_cross_{cam}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)
    print(f"[{cam}] fig -> {out}")

    return dict(
        cam=cam, label=C["label"], spec=C["spec"], n_fit=len(p),
        n_good=int((p.flags_psf == 0).sum()), n_cmp=len(ok),
        dmag_med=med, dmag_mad=mad,
        fit_window_px=fitw.get("shape_px"), fit_energy=fitw.get("energy_fraction", float("nan")),
        epsf_n_stars=er.get("n_selected"),
    )


def main() -> None:
    rows = [r for cam, C in CAMS.items() if (r := one_camera(cam, C))]
    if not rows:
        print("결과 없음")
        return
    lines = [
        "# LCO 교차기기 PSF 측광 결과 (APEX Step-8 헤드리스)",
        "",
        f"생성: 2026-07-29 · 데이터: archive.lco.global 공개 프레임 · "
        "보정: fig13 과 동일한 APEX 산술(BANZAI 대조 재현 확인 후 사용)",
        "",
        "| 카메라 | 프레임 | EPSF 별 | 적합 n | 양호 | Δmag(PSF−구경) med | MAD | 창(px)/에너지 |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append(
            f"| {r['label']} | {r['spec'].split(' · ')[0]} | {r['epsf_n_stars']} "
            f"| {r['n_fit']} | {r['n_good']} | {r['dmag_med']:+.3f} | {r['dmag_mad']:.3f} "
            f"| {r['fit_window_px']}/{r['fit_energy']:.2f} |"
        )
    lines += [
        "",
        "- Δmag 의 상수 오프셋은 EPSF 정규화 차이(CMD 에서는 프레임 ZP 가 흡수). 판정 기준은 MAD 와 잔차.",
        "- 한 프레임씩의 기능·일치 확인이며 다중 프레임 반복성은 범위 밖.",
        "- QHY 잔차에 남은 희미한 점원들은 `fit_init_max_sources=3000` 캡(혼잡장 탬플릿)으로 "
        "적합 대상에서 빠진 어두운 검출들이다 — 미적합이지 실패가 아니다. 밝은 별 몇 개의 "
        "도넛형 잔차는 포화·비선형 코어(예: Proxima 자체는 포화로 제외됨).",
        "- Sinistro 프레임은 형태 컷(epsf_sharp/round, Moravian 기준 튜닝)이 후보를 1개로 줄여 "
        "내장 완화(fallback) 로직이 59개로 복구했다 — 컷 기본값의 기기 의존성과 그 안전망이 실증됨.",
        "- **은하 필드 반응 (Sinistro=NGC5985, 2026-07-29 실측)**: EPSF 참조별 후보 125개 전부 "
        "은하 중심 141\" 밖(오염 0 — contamination 필터가 배제). 은하 본체(<100\") 검출은 "
        "46.2%가 품질 플래그로 표시됨(매듭·HII 는 점원이 아니라 qfit/noise 1.35 vs 필드 0.93) — "
        "조용히 오염되지 않고 표시하고 넘어간다. 은하 외곽(100–200\") 별들의 Δmag/MAD "
        "(+0.055/0.032)는 필드(+0.053/0.026)와 동일 — 은하가 주변 별 측광을 망치지 않는다. "
        "면광원 자체의 측광은 PSF 측광 범위 밖(잔차에 은하가 남는 것이 정상).",
        "- 그림: fig_psf_cross_sinistro.png / fig_psf_cross_qhy600.png",
    ]
    out = BASE / "REPORT.md"
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"report -> {out}")


if __name__ == "__main__":
    main()
