"""The paper says what the pipeline does; the pipeline has to keep saying it.

`paper/paper.md` used to claim APEX "runs both as a PyQt5 desktop application and
as a scriptable, headless command-line pipeline" and, in the same breath, that it
branches into isochrone fitting and light curves. A referee installing the base
package would have found neither: Steps 8 and 10 imported the GUI module for
their workers and reported NOT_IMPLEMENTED without PyQt5, and the isochrone and
light-curve steps are `DeferredStep`. `docs/audit/APEX_MANUSCRIPT_CLAIM_MATRIX.md`
had already flagged both lines — "Overstated", "Mostly supported for Steps 0-7".

Steps 8 and 10 moved into `apex.analysis` on 2026-08-17, so the first half is now
true and the wording says which half is not. These tests hold the two together:
if a step changes hands, the sentence in the paper has to change with it.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from apex.pipeline.base import DeferredStep
from apex.pipeline.registry import get_steps

REPO = Path(__file__).absolute().parents[1]
PAPER = REPO / "paper/paper.md"

# What the Summary paragraph names as running without a GUI toolkit.
CLAIMED_HEADLESS = {
    "scan", "crop", "sky", "detect", "wcs", "refbuild",
    "forcedphot", "psf", "zeropoint", "isochrone",
}
# What it names as desktop-only.
CLAIMED_DESKTOP_ONLY = {"masterid", "cmdplot"}


def test_the_steps_the_paper_calls_scriptable_are_implemented():
    runnable = {step.key for step in get_steps("cmd")
                if not isinstance(step, DeferredStep)}
    missing = sorted(CLAIMED_HEADLESS - runnable)
    assert not missing, (
        f"원고는 이 스텝들이 스크립트로 돈다고 말한다: {missing} — "
        f"구현하거나 paper/paper.md 의 Summary 를 고칠 것"
    )


def test_the_steps_the_paper_calls_desktop_only_still_are():
    """The reverse direction: implementing one is good news the paper must
    carry, not a silent divergence."""
    deferred = {step.key for step in get_steps("cmd")
                if isinstance(step, DeferredStep)}
    promoted = sorted(CLAIMED_DESKTOP_ONLY - deferred)
    assert not promoted, (
        f"이제 헤드리스로 도는데 원고는 데스크톱 전용이라고 말한다: {promoted} — "
        f"paper/paper.md 의 Summary 를 고칠 것"
    )


def test_the_paper_says_which_half_needs_a_desktop():
    """The sentence a referee would test on a base install."""
    text = PAPER.read_text(encoding="utf-8")
    assert "needs no GUI toolkit installed" in text
    assert "remain desktop-only" in text
    assert "declines to run until the settings" in text


def test_a_base_install_really_gets_the_pipeline_without_qt():
    """Not "no PyQt5 import statement" — no PyQt5 module ends up loaded.

    Run in a fresh interpreter: another test may already have imported Qt, and
    in-process reloading hands stale objects to the GUI subclasses.
    """
    probe = (
        "import sys;"
        "import apex.pipeline.registry;"
        "from apex.pipeline.registry import get_steps;"
        "get_steps('cmd'); get_steps('lc');"
        "print([m for m in sys.modules if m.startswith(('PyQt5', 'apex.gui'))])"
    )
    out = subprocess.run([sys.executable, "-X", "utf8", "-c", probe],
                         cwd=REPO, capture_output=True, text=True, timeout=300)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "[]", f"Qt/GUI 가 딸려 들어온다: {out.stdout}"


@pytest.mark.parametrize("module", [
    "apex.analysis.cmd.psf_photometry_runner",
    "apex.analysis.cmd.zeropoint_runner",
])
def test_the_calculation_modules_carry_no_qt(module):
    path = REPO / (module.replace(".", "/") + ".py")
    source = path.read_text(encoding="utf-8")
    import io
    import tokenize

    code = " ".join(
        tok.string
        for tok in tokenize.generate_tokens(io.StringIO(source).readline)
        if tok.type not in (tokenize.COMMENT, tokenize.STRING)
    )
    for word in ("PyQt5", "QThread", "pyqtSignal", "apex.gui"):
        assert word not in code, f"{module} 에 {word} 가 있다"
