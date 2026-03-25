"""Tests for ms_data_to_csv_exporter."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_test_mzml(path):
    import pyopenms as oms

    exp = oms.MSExperiment()
    for i in range(3):
        s = oms.MSSpectrum()
        s.setMSLevel(1 if i == 0 else 2)
        s.setRT(10.0 + i)
        s.set_peaks(([100.0 + j * 10 for j in range(5)], [1000.0 - j * 100 for j in range(5)]))
        if i > 0:
            prec = oms.Precursor()
            prec.setMZ(500.0 + i * 50)
            prec.setCharge(2)
            s.setPrecursors([prec])
        exp.addSpectrum(s)
    oms.MzMLFile().store(path, exp)


def _create_test_featurexml(path):
    import pyopenms as oms

    fm = oms.FeatureMap()
    for i in range(3):
        f = oms.Feature()
        f.setRT(100.0 + i * 10)
        f.setMZ(500.0 + i * 50)
        f.setIntensity(10000.0 + i * 1000)
        f.setCharge(2)
        f.setOverallQuality(0.9)
        fm.push_back(f)
    oms.FeatureXMLFile().store(path, fm)


def test_export_peaks():
    from ms_data_to_csv_exporter import export_mzml_peaks

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        tsv_path = os.path.join(tmp, "peaks.tsv")
        _create_test_mzml(mzml_path)

        stats = export_mzml_peaks(mzml_path, tsv_path)
        assert stats["spectra_exported"] == 3
        assert stats["total_peaks"] == 15

        with open(tsv_path) as fh:
            lines = fh.readlines()
        assert lines[0].strip().startswith("spectrum_index")
        assert len(lines) == 16  # header + 15 peaks


def test_export_peaks_ms_level_filter():
    from ms_data_to_csv_exporter import export_mzml_peaks

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        tsv_path = os.path.join(tmp, "peaks.tsv")
        _create_test_mzml(mzml_path)

        stats = export_mzml_peaks(mzml_path, tsv_path, ms_level=2)
        assert stats["spectra_exported"] == 2
        assert stats["total_peaks"] == 10


def test_export_spectra_summary():
    from ms_data_to_csv_exporter import export_mzml_spectra_summary

    with tempfile.TemporaryDirectory() as tmp:
        mzml_path = os.path.join(tmp, "test.mzML")
        tsv_path = os.path.join(tmp, "spectra.tsv")
        _create_test_mzml(mzml_path)

        stats = export_mzml_spectra_summary(mzml_path, tsv_path)
        assert stats["spectra_exported"] == 3


def test_export_featurexml():
    from ms_data_to_csv_exporter import export_featurexml

    with tempfile.TemporaryDirectory() as tmp:
        fxml_path = os.path.join(tmp, "test.featureXML")
        tsv_path = os.path.join(tmp, "features.tsv")
        _create_test_featurexml(fxml_path)

        stats = export_featurexml(fxml_path, tsv_path)
        assert stats["features_exported"] == 3

        with open(tsv_path) as fh:
            lines = fh.readlines()
        assert lines[0].strip().startswith("feature_id")
        assert len(lines) == 4
