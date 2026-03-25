"""Tests for mzml_to_mgf_converter."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
def test_create_synthetic_mzml():
    import pyopenms as oms
    from mzml_to_mgf_converter import create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        create_synthetic_mzml(mzml_path, n_spectra=3)

        exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, exp)
        assert exp.size() == 4  # 1 MS1 + 3 MS2


@requires_pyopenms
def test_convert_mzml_to_mgf():
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=3)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path)

        assert stats["converted"] == 3
        assert stats["ms_level_spectra"] == 3

        with open(mgf_path) as fh:
            content = fh.read()
        assert content.count("BEGIN IONS") == 3
        assert content.count("END IONS") == 3
        assert "PEPMASS=" in content
        assert "CHARGE=" in content


@requires_pyopenms
def test_mgf_content_format():
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=1)
        convert_mzml_to_mgf(mzml_path, mgf_path)

        with open(mgf_path) as fh:
            lines = fh.readlines()

        # Check MGF format structure
        assert lines[0].strip() == "BEGIN IONS"
        has_title = any(line.startswith("TITLE=") for line in lines)
        has_pepmass = any(line.startswith("PEPMASS=") for line in lines)
        assert has_title
        assert has_pepmass


@requires_pyopenms
def test_min_peaks_filter():
    from mzml_to_mgf_converter import convert_mzml_to_mgf, create_synthetic_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        mgf_path = os.path.join(tmp, "test.mgf")

        create_synthetic_mzml(mzml_path, n_spectra=3)
        stats = convert_mzml_to_mgf(mzml_path, mgf_path, min_peaks=100)
        assert stats["converted"] == 0
