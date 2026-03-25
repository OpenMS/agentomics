"""Tests for mzqc_generator."""

import pytest

pytest.importorskip("pyopenms")


class TestMzqcGenerator:
    def _make_experiment(self):
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(3):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 * i)
            mzs = np.array([100.0 + j for j in range(5)], dtype=np.float64)
            ints = np.array([1000.0] * 5, dtype=np.float64)
            spec.set_peaks([mzs, ints])
            exp.addSpectrum(spec)

            ms2 = oms.MSSpectrum()
            ms2.setMSLevel(2)
            ms2.setRT(60.0 * i + 0.5)
            prec = oms.Precursor()
            prec.setMZ(500.0)
            ms2.setPrecursors([prec])
            mzs2 = np.array([200.0, 300.0], dtype=np.float64)
            ints2 = np.array([500.0, 800.0], dtype=np.float64)
            ms2.set_peaks([mzs2, ints2])
            exp.addSpectrum(ms2)

        return exp

    def test_generate_mzqc_structure(self):
        from mzqc_generator import generate_mzqc

        exp = self._make_experiment()
        mzqc = generate_mzqc(exp, "test.mzML")
        assert "mzQC" in mzqc
        assert "runQualities" in mzqc["mzQC"]
        assert len(mzqc["mzQC"]["runQualities"]) == 1

    def test_metric_count(self):
        from mzqc_generator import generate_mzqc

        exp = self._make_experiment()
        mzqc = generate_mzqc(exp)
        metrics = mzqc["mzQC"]["runQualities"][0]["qualityMetrics"]
        assert len(metrics) == 6

    def test_ms1_count_metric(self):
        from mzqc_generator import generate_mzqc

        exp = self._make_experiment()
        mzqc = generate_mzqc(exp)
        metrics = mzqc["mzQC"]["runQualities"][0]["qualityMetrics"]
        ms1_metric = next(m for m in metrics if m["accession"] == "QC:0000005")
        assert ms1_metric["value"] == 3

    def test_ms2_count_metric(self):
        from mzqc_generator import generate_mzqc

        exp = self._make_experiment()
        mzqc = generate_mzqc(exp)
        metrics = mzqc["mzQC"]["runQualities"][0]["qualityMetrics"]
        ms2_metric = next(m for m in metrics if m["accession"] == "QC:0000006")
        assert ms2_metric["value"] == 3

    def test_version_field(self):
        from mzqc_generator import generate_mzqc

        exp = self._make_experiment()
        mzqc = generate_mzqc(exp)
        assert mzqc["mzQC"]["version"] == "1.0.0"
