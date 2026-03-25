"""Tests for sirius_exporter."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_test_data(tmp_dir):
    """Create test mzML and features TSV files."""
    import pyopenms as oms

    # Create mzML with MS2 spectra
    exp = oms.MSExperiment()
    for i in range(3):
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(100.0 + i * 10)
        prec = oms.Precursor()
        prec.setMZ(500.0 + i * 50)
        prec.setCharge(1)
        ms2.setPrecursors([prec])
        ms2.set_peaks(([100.0 + j * 20 for j in range(5)], [1000.0 - j * 100 for j in range(5)]))
        exp.addSpectrum(ms2)

    mzml_path = os.path.join(tmp_dir, "test.mzML")
    oms.MzMLFile().store(mzml_path, exp)

    # Create features TSV
    features_path = os.path.join(tmp_dir, "features.tsv")
    with open(features_path, "w") as fh:
        fh.write("mz\trt\tcharge\tname\n")
        fh.write("500.0\t100.0\t1\tcompound_A\n")
        fh.write("550.0\t110.0\t1\tcompound_B\n")
        fh.write("999.0\t999.0\t1\tno_match\n")  # no matching MS2

    return mzml_path, features_path


def test_load_features_tsv():
    from sirius_exporter import load_features_tsv

    with tempfile.TemporaryDirectory() as tmp:
        features_path = os.path.join(tmp, "features.tsv")
        with open(features_path, "w") as fh:
            fh.write("mz\trt\tcharge\tname\n")
            fh.write("500.0\t100.0\t1\ttest_compound\n")

        features = load_features_tsv(features_path)
        assert len(features) == 1
        assert features[0]["mz"] == 500.0
        assert features[0]["name"] == "test_compound"


def test_find_ms2_spectra():
    from sirius_exporter import find_ms2_spectra, load_mzml

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path, _ = _create_test_data(tmp)
        exp = load_mzml(mzml_path)

        matches = find_ms2_spectra(exp, 500.0, 100.0, mz_tolerance=0.01, rt_tolerance=5.0)
        assert len(matches) == 1

        no_match = find_ms2_spectra(exp, 999.0, 999.0, mz_tolerance=0.01, rt_tolerance=5.0)
        assert len(no_match) == 0


def test_export_to_sirius():
    from sirius_exporter import export_to_sirius

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path, features_path = _create_test_data(tmp)
        output_path = os.path.join(tmp, "sirius.ms")

        stats = export_to_sirius(features_path, mzml_path, output_path)
        assert stats["features_exported"] == 3
        assert stats["features_with_ms2"] == 2  # third feature has no match

        with open(output_path) as fh:
            content = fh.read()
        assert ">compound compound_A" in content
        assert ">parentmass" in content
        assert ">ms2" in content


def test_sirius_ms_format():
    from sirius_exporter import export_to_sirius

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path, features_path = _create_test_data(tmp)
        output_path = os.path.join(tmp, "sirius.ms")

        export_to_sirius(features_path, mzml_path, output_path)

        with open(output_path) as fh:
            content = fh.read()

        # Verify .ms format structure
        assert ">compound" in content
        assert ">parentmass" in content
        assert ">rt" in content
