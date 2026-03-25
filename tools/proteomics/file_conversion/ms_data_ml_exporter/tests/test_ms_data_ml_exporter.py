"""Tests for ms_data_ml_exporter."""

import numpy as np
import pytest

pytest.importorskip("pyopenms")


class TestMsDataMlExporter:
    def _make_experiment(self, n_spectra=5):
        """Create a synthetic MSExperiment."""
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(n_spectra):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1 if i % 2 == 0 else 2)
            spec.setRT(60.0 * i)
            mzs = np.array([100.0 + j * 10 for j in range(20)], dtype=np.float64)
            intensities = np.array([1000.0 * (j + 1) for j in range(20)], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)
        return exp

    def test_extract_features(self):
        from ms_data_ml_exporter import extract_features

        exp = self._make_experiment(n_spectra=3)
        records = extract_features(exp)
        assert len(records) == 3
        assert all("rt" in r for r in records)
        assert all("tic" in r for r in records)
        assert all("n_peaks" in r for r in records)

    def test_feature_values(self):
        from ms_data_ml_exporter import extract_features

        exp = self._make_experiment(n_spectra=1)
        records = extract_features(exp)
        r = records[0]
        assert r["n_peaks"] == 20
        assert r["mz_min"] > 0
        assert r["base_peak_intensity"] > 0
        assert r["intensity_std"] >= 0

    def test_empty_experiment(self):
        import pyopenms as oms
        from ms_data_ml_exporter import extract_features

        exp = oms.MSExperiment()
        records = extract_features(exp)
        assert records == []

    def test_write_csv(self, tmp_path):
        from ms_data_ml_exporter import extract_features, write_csv

        exp = self._make_experiment(n_spectra=3)
        records = extract_features(exp)
        out = str(tmp_path / "matrix.csv")
        write_csv(records, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 4  # header + 3 rows
        assert "spectrum_index" in lines[0]
