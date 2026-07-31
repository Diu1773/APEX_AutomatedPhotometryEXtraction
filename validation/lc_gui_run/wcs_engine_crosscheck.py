"""WCS 엔진 교차검증 — 같은 프레임을 internal / astnet 으로 각각 풀어 비교한다.

YZ Bootis 한 필드에서 internal 이 astnet 보다 rms 가 2.9배 좋게 나왔는데,
그것이 그 자료 특유인지 일반적인지 가리기 위한 것이다. 대상마다 Step 5 만
`--force` 로 다시 돌리고 `frame_wcs_qc.csv` 의 잔차를 모은다.

    python wcs_engine_crosscheck.py <name>=<parameters.toml> [...]

원래 step5_wcs 는 `step5_wcs_crosscheck_backup` 으로 보존했다가 마지막에 되돌린다.
"""
from __future__ import annotations

import json
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
VENV_PY = REPO / ".venv-deploy" / "Scripts" / "python.exe"


def _read_result_dir(param_file: Path) -> Path:
    for line in param_file.read_text(encoding="utf-8").splitlines():
        s = line.strip()
        if s.startswith("result_dir"):
            return Path(s.split("=", 1)[1].strip().strip('"').replace("\\\\", "\\"))
    raise SystemExit(f"result_dir 을 못 찾았다: {param_file}")


def _set_engine(param_file: Path, engine: str) -> None:
    """[wcs] 바로 아래에 engine 을 넣는다.

    [wcs] 안의 기존 engine 줄은 **먼저 지우고** 새로 넣는다. 예전 판은 삽입과
    치환을 한 번에 하려다 두 번째 엔진으로 바꿀 때 키를 중복으로 넣어 TOML 이
    깨졌고, 그래서 astnet 이 4초 만에 "실패"했다(엔진 문제가 아니었다).
    """
    lines = param_file.read_text(encoding="utf-8").splitlines()

    kept, in_wcs = [], False
    for line in lines:
        s = line.strip()
        if s.startswith("[") and s.endswith("]"):
            in_wcs = s == "[wcs]"
        elif in_wcs and re.match(r"^engine\s*=", s):
            continue          # [wcs] 안의 기존 engine 은 버린다
        kept.append(line)

    out = []
    for line in kept:
        out.append(line)
        if line.strip() == "[wcs]":
            out.append(f'engine = "{engine}"')
    param_file.write_text("\n".join(out) + "\n", encoding="utf-8")

    # 넣은 대로 읽히는지 확인 — 조용히 깨진 채로 비교하면 결론이 통째로 틀린다.
    import tomllib

    got = tomllib.loads(param_file.read_text(encoding="utf-8")).get("wcs", {}).get("engine")
    if got != engine:
        raise SystemExit(f"engine 설정 실패: 기대 {engine!r}, 실제 {got!r}")


def _run_step5(param_file: Path) -> tuple[bool, float, str]:
    t0 = time.perf_counter()
    p = subprocess.run(
        [str(VENV_PY), "-X", "utf8", "-m", "apex.cli", "run",
         "--mode", "cmd", "--config", str(param_file), "--steps", "5", "--force"],
        cwd=str(REPO), capture_output=True, text=True, encoding="utf-8",
        errors="replace", timeout=7200,
    )
    dt = time.perf_counter() - t0
    line = ""
    for ln in (p.stdout or "").splitlines():
        if ln.startswith("Step 5") and "->" in ln:
            line = ln.strip()
    if p.returncode != 0:
        tail = "\n".join(((p.stdout or "") + (p.stderr or "")).splitlines()[-6:])
        line = f"[실패] {tail}"
    return p.returncode == 0, dt, line


def _qc_stats(result_dir: Path) -> dict:
    import pandas as pd

    qc = result_dir / "step5_wcs" / "frame_wcs_qc.csv"
    if not qc.exists():
        return {"error": "frame_wcs_qc.csv 없음"}
    d = pd.read_csv(qc)
    out = {"n": len(d)}
    if "wcs_qc_pass" in d:
        out["pass"] = int(d["wcs_qc_pass"].astype(bool).sum())
    for c in ("rms_px", "resid_med_px", "resid_p99_px", "match_rate", "n_match"):
        if c in d:
            v = pd.to_numeric(d[c], errors="coerce")
            if v.notna().any():
                out[c] = round(float(v.median()), 4)
    return out


def main() -> int:
    targets = []
    for arg in sys.argv[1:]:
        name, _, path = arg.partition("=")
        targets.append((name, Path(path)))
    if not targets:
        raise SystemExit(__doc__)

    report = {}
    for name, pf in targets:
        result_dir = _read_result_dir(pf)
        step5 = result_dir / "step5_wcs"
        backup = result_dir / "step5_wcs_crosscheck_backup"
        original = pf.read_text(encoding="utf-8")
        if step5.exists() and not backup.exists():
            shutil.copytree(step5, backup)
        print(f"\n=== {name} ({result_dir}) ===", flush=True)
        report[name] = {}
        try:
            for engine in ("internal", "astnet"):
                _set_engine(pf, engine)
                ok, dt, line = _run_step5(pf)
                stats = _qc_stats(result_dir) if ok else {"error": "step5 실패"}
                stats["elapsed_s"] = round(dt, 1)
                stats["ok"] = ok
                report[name][engine] = stats
                print(f"  {engine:<9} {dt:7.1f}s  {stats}", flush=True)
                if line:
                    print(f"     {line}", flush=True)
        finally:
            pf.write_text(original, encoding="utf-8")   # 설정 원복

    out = Path(__file__).parent / "wcs_engine_crosscheck.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n[saved] {out}")

    print(f"\n{'대상':<12} {'엔진':<9} {'시간s':>8} {'rms_px':>8} {'통과':>8}")
    for name, per in report.items():
        for engine, s in per.items():
            print(f"{name:<12} {engine:<9} {s.get('elapsed_s','-'):>8} "
                  f"{s.get('rms_px','-'):>8} {str(s.get('pass','-'))+'/'+str(s.get('n','-')):>8}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
