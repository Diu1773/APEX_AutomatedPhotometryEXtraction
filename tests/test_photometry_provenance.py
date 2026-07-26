from __future__ import annotations

import json

import pandas as pd

from apex.gui.workflow.cmd.step10_zeropoint_calibration import (
    resolve_cmd_photometry_input,
    resolve_cmd_photometry_provenance,
)
from apex.core.project_state import ProjectState
from apex.utils.photometry_provenance import (
    build_photometry_provenance,
    format_photometry_provenance,
    summarize_photometry_table,
)
from apex.utils.step_paths import step7_forced_phot_dir
from apex.utils.step_paths_cmd import step8_psf_dir, step10_zp_dir


def _write_complete_psf_source(root, frame: str = "frame001.fits") -> None:
    aperture_dir = step7_forced_phot_dir(root)
    aperture_dir.mkdir(parents=True)
    aperture_index = aperture_dir / "photometry_index.csv"
    pd.DataFrame({"file": [frame], "filter": ["V"]}).to_csv(
        aperture_index, index=False
    )
    aperture_table = aperture_dir / f"photometry_{frame}.tsv"
    pd.DataFrame(
        {"det_uid": [1], "mag": [14.0], "mag_err": [0.02]}
    ).to_csv(aperture_table, sep="\t", index=False)

    psf_dir = step8_psf_dir(root)
    psf_dir.mkdir(parents=True)
    pd.DataFrame({"file": [frame], "filter": ["V"]}).to_csv(
        psf_dir / "photometry_index.csv", index=False
    )
    pd.DataFrame(
        {
            "det_uid": [1],
            "mag_psf": [13.9],
            "mag_psf_err": [0.01],
            "flags_psf": [0],
        }
    ).to_csv(psf_dir / f"photometry_{frame}.tsv", sep="\t", index=False)

    index_stat = aperture_index.stat()
    table_stat = aperture_table.stat()
    signature = {
        "frames": [frame],
        "inputs": {
            "step7_index": {
                "size": index_stat.st_size,
                "mtime_ns": index_stat.st_mtime_ns,
            },
            "frames": [
                {
                    "file": frame,
                    "step7_tsv": {
                        "size": table_stat.st_size,
                        "mtime_ns": table_stat.st_mtime_ns,
                    },
                }
            ],
        },
    }
    (psf_dir / "psf_output_signature.json").write_text(
        json.dumps(signature), encoding="utf-8"
    )


def test_psf_provenance_formats_selected_magnitude_column():
    info = build_photometry_provenance("psf", "mag_psf", "mag_psf_err")

    assert info == {
        "source": "psf",
        "mag_column": "mag_psf",
        "mag_error_column": "mag_psf_err",
    }
    assert format_photometry_provenance(info) == "Photometry: PSF | MAG: mag_psf"


def test_missing_provenance_is_not_mislabeled_as_aperture():
    assert format_photometry_provenance({}) == (
        "Photometry: Unknown | MAG: unknown"
    )


def test_table_summary_reads_filter_specific_cmd_provenance():
    table = pd.DataFrame(
        {
            "photometry_source_B": ["psf", "psf"],
            "photometry_source_V": ["psf", "psf"],
            "mag_input_column_B": ["mag_psf", "mag_psf"],
            "mag_input_column_V": ["mag_psf", "mag_psf"],
            "mag_error_input_column_B": ["mag_psf_err", "mag_psf_err"],
            "mag_error_input_column_V": ["mag_psf_err", "mag_psf_err"],
        }
    )

    assert summarize_photometry_table(table) == {
        "source": "psf",
        "mag_column": "mag_psf",
        "mag_error_column": "mag_psf_err",
    }


def test_table_summary_reports_mixed_sources_without_hiding_it():
    table = pd.DataFrame(
        {
            "photometry_source": ["aperture", "psf"],
            "mag_input_column": ["mag", "mag_psf"],
            "mag_error_input_column": ["mag_err", "mag_psf_err"],
        }
    )

    info = summarize_photometry_table(table)

    assert info["source"] == "mixed"
    assert info["mag_column"] == "mixed"
    assert format_photometry_provenance(info) == "Photometry: Mixed | MAG: mixed"


def test_cmd_resolver_prefers_provenance_saved_in_step10_output(tmp_path):
    output_dir = step10_zp_dir(tmp_path)
    output_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "ID": [1],
            "photometry_source_V": ["psf"],
            "mag_input_column_V": ["mag_psf"],
            "mag_error_input_column_V": ["mag_psf_err"],
        }
    ).to_csv(output_dir / "median_by_ID_filter_wide_cmd.csv", index=False)

    assert resolve_cmd_photometry_provenance(tmp_path) == {
        "source": "psf",
        "mag_column": "mag_psf",
        "mag_error_column": "mag_psf_err",
    }


def test_cmd_input_uses_complete_psf_output_by_default(tmp_path):
    _write_complete_psf_source(tmp_path)

    source = resolve_cmd_photometry_input(tmp_path)

    assert source["source"] == "psf"
    assert source["mag_column"] == "mag_psf"
    assert source["index_path"] == step8_psf_dir(tmp_path) / "photometry_index.csv"


def test_cmd_input_uses_aperture_when_psf_was_skipped(tmp_path):
    _write_complete_psf_source(tmp_path)
    state = ProjectState(tmp_path)
    state.store_step_data("psf_photometry", {"skip_psf": True})

    source = resolve_cmd_photometry_input(tmp_path, state)

    assert source["source"] == "aperture"
    assert source["mag_column"] == "mag"
    assert source["index_path"] == step7_forced_phot_dir(tmp_path) / "photometry_index.csv"


def test_cmd_input_falls_back_when_step7_changed_after_psf(tmp_path):
    _write_complete_psf_source(tmp_path)
    index_path = step7_forced_phot_dir(tmp_path) / "photometry_index.csv"
    index_path.write_text(
        index_path.read_text(encoding="utf-8") + "\n",
        encoding="utf-8",
    )

    source = resolve_cmd_photometry_input(tmp_path)

    assert source["source"] == "aperture"
    assert "Step 7 changed" in source["reason"]
