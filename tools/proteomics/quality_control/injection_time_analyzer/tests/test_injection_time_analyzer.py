"""Tests for injection_time_analyzer."""

import pytest

pytest.importorskip("pyopenms")


class TestInjectionTimeAnalyzer:
    def _make_experiment_with_injection_times(self):
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        for i in range(4):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1 if i % 2 == 0 else 2)
            spec.setRT(10.0 * i)
            spec.setMetaValue("MS:1000927", float(20.0 + i * 5))
            mzs = np.array([100.0], dtype=np.float64)
            ints = np.array([1000.0], dtype=np.float64)
            spec.set_peaks([mzs, ints])
            exp.addSpectrum(spec)
        return exp

    def test_extract_injection_times(self):
        from injection_time_analyzer import extract_injection_times

        exp = self._make_experiment_with_injection_times()
        records = extract_injection_times(exp)
        assert len(records) == 4
        # Check that injection times were extracted
        with_times = [r for r in records if r["injection_time_ms"] is not None]
        assert len(with_times) == 4

    def test_injection_time_values(self):
        from injection_time_analyzer import extract_injection_times

        exp = self._make_experiment_with_injection_times()
        records = extract_injection_times(exp)
        assert records[0]["injection_time_ms"] == 20.0
        assert records[1]["injection_time_ms"] == 25.0

    def test_summarize_injection_times(self):
        from injection_time_analyzer import summarize_injection_times

        records = [
            {"ms_level": 1, "injection_time_ms": 20.0},
            {"ms_level": 1, "injection_time_ms": 30.0},
            {"ms_level": 2, "injection_time_ms": 50.0},
            {"ms_level": 2, "injection_time_ms": 60.0},
        ]
        summary = summarize_injection_times(records)
        assert "MS1" in summary
        assert "MS2" in summary
        assert summary["MS1"]["mean_ms"] == 25.0
        assert summary["MS2"]["mean_ms"] == 55.0

    def test_no_injection_times(self):
        import numpy as np
        import pyopenms as oms
        from injection_time_analyzer import extract_injection_times

        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(0.0)
        spec.set_peaks([np.array([100.0], dtype=np.float64),
                       np.array([1000.0], dtype=np.float64)])
        exp.addSpectrum(spec)

        records = extract_injection_times(exp)
        assert len(records) == 1
        assert records[0]["injection_time_ms"] is None
