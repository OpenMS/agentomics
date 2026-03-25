"""Tests for gnps_fbmn_exporter."""

import csv
import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestGnpsFbmnExporter:
    def _create_test_mzml(self, path: str) -> None:
        """Create a minimal mzML file with MS1 and MS2 spectra."""
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()

        # MS1 spectrum
        ms1 = oms.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(100.0)
        ms1.set_peaks([np.array([180.0634, 200.0], dtype=np.float64),
                        np.array([1000.0, 500.0], dtype=np.float64)])
        exp.addSpectrum(ms1)

        # MS2 spectrum matching feature at mz=180.0634, rt=100.0
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(100.0)
        prec = oms.Precursor()
        prec.setMZ(180.0634)
        ms2.setPrecursors([prec])
        ms2.set_peaks([np.array([90.0, 120.0, 150.0], dtype=np.float64),
                        np.array([500.0, 200.0, 100.0], dtype=np.float64)])
        exp.addSpectrum(ms2)

        oms.MzMLFile().store(path, exp)

    def test_load_mzml(self):
        from gnps_fbmn_exporter import load_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            self._create_test_mzml(mzml_path)
            exp = load_mzml(mzml_path)
            assert exp.size() == 2

    def test_find_best_ms2(self):
        from gnps_fbmn_exporter import find_best_ms2, load_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            self._create_test_mzml(mzml_path)
            exp = load_mzml(mzml_path)

            peaks = find_best_ms2(exp, mz=180.0634, rt=100.0, mz_tol=0.01, rt_tol=30.0)
            assert peaks is not None
            assert len(peaks) == 3

    def test_find_best_ms2_no_match(self):
        from gnps_fbmn_exporter import find_best_ms2, load_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            self._create_test_mzml(mzml_path)
            exp = load_mzml(mzml_path)

            peaks = find_best_ms2(exp, mz=300.0, rt=100.0, mz_tol=0.01, rt_tol=30.0)
            assert peaks is None

    def test_write_mgf(self):
        from gnps_fbmn_exporter import write_mgf

        features = [{"feature_id": "F1", "mz": 180.0634, "rt": 100.0}]
        spectra = {"F1": [(90.0, 500.0), (120.0, 200.0)]}

        with tempfile.TemporaryDirectory() as tmpdir:
            mgf_path = os.path.join(tmpdir, "test.mgf")
            count = write_mgf(features, spectra, mgf_path)
            assert count == 1
            with open(mgf_path) as fh:
                content = fh.read()
            assert "BEGIN IONS" in content
            assert "SCANS=F1" in content
            assert "PEPMASS=180.063400" in content
            assert "RTINSECONDS=100.00" in content
            assert "END IONS" in content

    def test_write_quant_table(self):
        from gnps_fbmn_exporter import write_quant_table

        features = [{"feature_id": "F1", "mz": 180.0634, "rt": 100.0, "intensity": 1000.0}]

        with tempfile.TemporaryDirectory() as tmpdir:
            quant_path = os.path.join(tmpdir, "quant.csv")
            write_quant_table(features, quant_path)
            assert os.path.exists(quant_path)
            with open(quant_path) as fh:
                reader = csv.DictReader(fh)
                rows = list(reader)
            assert len(rows) == 1
            assert "row ID" in rows[0]

    def test_full_pipeline(self):
        from gnps_fbmn_exporter import export_fbmn

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            self._create_test_mzml(mzml_path)

            features = [{"feature_id": "F1", "mz": 180.0634, "rt": 100.0, "intensity": 1000.0}]

            mgf_path = os.path.join(tmpdir, "gnps.mgf")
            quant_path = os.path.join(tmpdir, "quant.csv")
            n_spectra, n_features = export_fbmn(mzml_path, features, mgf_path, quant_path)
            assert n_spectra == 1
            assert n_features == 1
            assert os.path.exists(mgf_path)
            assert os.path.exists(quant_path)
