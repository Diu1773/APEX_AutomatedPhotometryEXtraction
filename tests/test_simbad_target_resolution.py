from __future__ import annotations

from types import SimpleNamespace

from astropy.table import Table

from apex.core import instrument as instrument_mod
from apex.core.instrument import InstrumentConfig


def _params(tmp_path):
    return SimpleNamespace(
        P=SimpleNamespace(
            data_dir=tmp_path,
            result_dir=tmp_path / "result",
            camera_binning=2,
            fwhm_px_min=1.0,
            fwhm_px_max=20.0,
        )
    )


def test_messier_aliases_include_catalog_names():
    assert "NGC 5272" in InstrumentConfig._target_query_aliases("M3")
    assert "NGC 2099" in InstrumentConfig._target_query_aliases("M37")
    assert "NGC 2099" in InstrumentConfig._target_query_aliases("Messier 37")


def test_known_target_fails_when_simbad_is_unavailable(tmp_path, monkeypatch):
    logs: list[str] = []
    inst = InstrumentConfig(_params(tmp_path))
    monkeypatch.setattr(inst, "_prepare_simbad", lambda: False)
    inst.resolve_targets(["M3"], log_fn=logs.append)

    joined = "\n".join(logs)
    assert inst.primary_coord is None
    assert inst.targets_resolved == []
    assert "unavailable: M3" in joined


def test_m37_attempts_ngc2099_alias_when_simbad_unavailable(tmp_path, monkeypatch):
    logs: list[str] = []
    inst = InstrumentConfig(_params(tmp_path))
    monkeypatch.setattr(inst, "_prepare_simbad", lambda: False)
    inst.resolve_targets(["M37"], log_fn=logs.append)

    joined = "\n".join(logs)
    assert inst.primary_coord is None
    assert inst.targets_resolved == []
    assert "NGC 2099" in inst.last_target_attempts


def test_m37_retries_second_simbad_server_with_ngc2099_alias(tmp_path, monkeypatch):
    inst = InstrumentConfig(_params(tmp_path))
    monkeypatch.setattr(inst, "_prepare_simbad", lambda: True)
    monkeypatch.setattr(
        inst,
        "_simbad_server_candidates",
        lambda: ["first.example", "simbad.harvard.edu"],
    )

    calls: list[tuple[str, str]] = []

    def fake_query(name: str):
        calls.append((instrument_mod.Simbad.SIMBAD_URL, name))
        if "simbad.harvard.edu" in instrument_mod.Simbad.SIMBAD_URL and name == "NGC 2099":
            return Table(
                {
                    "RA": ["05 52 17.81"],
                    "DEC": ["+32 32 42.0"],
                    "OTYPE": ["OpC"],
                }
            )
        return None

    monkeypatch.setattr(instrument_mod.Simbad, "query_object", fake_query)

    inst.resolve_targets(["M37"], log_fn=lambda _msg: None)

    assert inst.primary_target == "M37"
    assert inst.primary_coord is not None
    assert inst.targets_resolved[0]["simbad_query"] == "NGC 2099"
    assert any("first.example" in url for url, _ in calls)
    assert any("simbad.harvard.edu" in url for url, _ in calls)


def test_simbad_numeric_lowercase_ra_dec_are_degrees():
    res = Table(
        {
            "ra": [88.07416666666667],
            "dec": [32.54499999999999],
            "otype": ["OpC"],
        }
    )

    rec = InstrumentConfig._record_from_simbad_row("M37", res, "M37")

    assert abs(rec["ra_deg"] - 88.07416666666667) < 1e-9
    assert abs(rec["dec_deg"] - 32.54499999999999) < 1e-9


def test_explicit_target_resolution_does_not_append_stale_saved_targets(tmp_path, monkeypatch):
    params = _params(tmp_path)
    params.P._raw = {"target_name": "NGC457"}
    (tmp_path / "targets.txt").write_text("NGC457\n", encoding="utf-8")

    logs: list[str] = []
    inst = InstrumentConfig(params)
    monkeypatch.setattr(inst, "_prepare_simbad", lambda: False)
    inst.resolve_targets(["M37"], log_fn=logs.append)

    joined = "\n".join(logs)
    assert inst.primary_coord is None
    assert inst.targets_resolved == []
    assert "NGC457" not in joined


def test_empty_target_list_does_not_use_default_or_saved_targets(tmp_path, monkeypatch):
    params = _params(tmp_path)
    params.P._raw = {"target_name": "NGC457"}
    (tmp_path / "targets.txt").write_text("M31\n", encoding="utf-8")

    inst = InstrumentConfig(params)
    monkeypatch.setattr(inst, "_prepare_simbad", lambda: False)
    inst.resolve_targets([], log_fn=lambda _msg: None)

    assert inst.primary_coord is None
    assert inst.targets_resolved == []
    assert inst.last_target_attempts == []


def test_failed_resolution_removes_stale_simbad_sidecar(tmp_path, monkeypatch):
    params = _params(tmp_path)
    sidecar_dir = tmp_path / "result" / "step1_file_selection"
    sidecar_dir.mkdir(parents=True)
    stale = sidecar_dir / "targets_simbad.tsv"
    stale.write_text("name\tra_deg\tdec_deg\nNGC457\t298.3\t58.2\n", encoding="utf-8")

    inst = InstrumentConfig(params)
    monkeypatch.setattr(inst, "_prepare_simbad", lambda: False)
    inst.resolve_targets(["DefinitelyNotATarget"], log_fn=lambda _msg: None)

    assert inst.primary_coord is None
    assert inst.targets_resolved == []
    assert not stale.exists()
