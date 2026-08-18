"""A result must be able to name the settings that made it.

Before 2026-08-18 a run recorded which steps ran, how long they took, and (from
that morning) the package versions. Not one parameter value. A table in
`cmd_zeropoint/` therefore rested on the `apex_config.json` beside it being
unchanged since — and this session rewrote all fifty of those three times over.

For a manuscript that is not good enough, so every run now writes the resolved
settings, the config file's SHA-256, and — measured, not declared — which step
read which setting.
"""

from __future__ import annotations

import csv
import json
import logging
from pathlib import Path
from types import SimpleNamespace

from apex.pipeline.provenance import RecordingNamespace, write_parameter_record


def test_the_recorder_notes_reads_and_stays_out_of_the_way():
    P = SimpleNamespace(alpha=1, beta=2, gamma=3)
    proxy = RecordingNamespace(P)

    assert proxy.alpha == 1
    assert proxy.gamma == 3
    assert proxy.seen == {"alpha", "gamma"}, "읽은 것만 잡아야 한다"

    proxy.delta = 9
    assert P.delta == 9, "쓰기가 진짜 객체로 가야 한다"
    assert "delta" not in proxy.seen, "쓰기는 읽기가 아니다"

    assert getattr(proxy, "_target", None) is not None
    assert "_target" not in proxy.seen, "내부 이름은 기록하지 않는다"


def test_it_survives_a_missing_attribute_the_way_the_real_one_does():
    proxy = RecordingNamespace(SimpleNamespace(a=1))
    assert getattr(proxy, "nope", "fallback") == "fallback"


def test_the_record_carries_the_settings_and_the_config_hash(tmp_path):
    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({"io": {"result_dir": "."}}), encoding="utf-8")
    params = SimpleNamespace(
        P=SimpleNamespace(apcorr_small_scale=0.8, detect_sigma=3.2, _hidden=1),
        param_file=str(config),
    )
    written = write_parameter_record(
        tmp_path, params, "cmd",
        {"forcedphot": {"apcorr_small_scale"}, "detect": {"detect_sigma"}},
        environment={"packages": {"scipy": "1.17.1"}},
    )
    assert [p.name for p in written] == ["parameters_used.json", "parameters_used.csv"]

    body = json.loads(written[0].read_text(encoding="utf-8"))
    assert body["settings"]["apcorr_small_scale"] == 0.8
    assert "_hidden" not in body["settings"], "내부 필드는 기록에서 뺀다"
    assert body["config"]["sha256"] and len(body["config"]["sha256"]) == 64
    assert body["environment"]["packages"]["scipy"] == "1.17.1"
    assert body["read_by_step"]["apcorr_small_scale"] == ["forcedphot"]
    assert body["steps_recorded"] == ["detect", "forcedphot"]


def test_the_csv_is_the_one_that_goes_in_an_appendix(tmp_path):
    params = SimpleNamespace(P=SimpleNamespace(a=1, b=2), param_file=None)
    _, csv_path = write_parameter_record(
        tmp_path, params, "cmd", {"scan": {"a"}})
    rows = list(csv.reader(csv_path.read_text(encoding="utf-8").splitlines()))
    assert rows[0] == ["setting", "value", "read_by_steps"]
    body = {r[0]: r for r in rows[1:]}
    assert body["a"][2] == "scan"
    assert body["b"][2] == "", "아무 스텝도 안 읽은 설정은 빈 칸"


def test_a_real_run_records_exactly_what_the_step_read(tmp_path):
    """End to end through the runner, not through the helper.

    Step 1 reads ten settings — the io paths, the site, the target. If the proxy
    were not installed this would record none of them; if it leaked past the
    step it would record everything.
    """
    from apex.pipeline.base import PipelineStep, StepResult, StepStatus
    from apex.pipeline.context import RunContext
    from apex.pipeline.runner import PipelineRunner

    class Probe(PipelineStep):
        index, key, name = 1, "probe", "Probe"

        def inputs(self, ctx):
            return []

        def outputs(self, ctx):
            return []

        def is_complete(self, ctx):
            return False

        def run(self, ctx):
            _ = ctx.params.P.wanted
            return StepResult(index=self.index, key=self.key,
                              status=StepStatus.OK, message="ok")

    params = SimpleNamespace(
        P=SimpleNamespace(wanted=1, ignored=2, result_dir=str(tmp_path)),
        param_file=None)
    ctx = RunContext(mode="cmd", params=params, result_dir=tmp_path,
                     data_dir=tmp_path, logger=logging.getLogger("test"),
                     force=True)
    PipelineRunner(steps=[Probe()]).run(ctx)

    body = json.loads((tmp_path / "parameters_used.json").read_text(encoding="utf-8"))
    assert body["read_by_step"].get("wanted") == ["probe"]
    assert "ignored" not in body["read_by_step"], "안 읽은 설정이 읽힌 것으로 기록됐다"
    assert params.P.wanted == 1, "실행 뒤 진짜 네임스페이스가 돌아와야 한다"
    assert not isinstance(params.P, RecordingNamespace)
