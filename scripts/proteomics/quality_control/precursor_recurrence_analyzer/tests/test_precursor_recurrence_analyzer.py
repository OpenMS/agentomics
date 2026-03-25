"""Tests for precursor_recurrence_analyzer."""

import numpy as np
from conftest import requires_pyopenms


@requires_pyopenms
class TestPrecursorRecurrenceAnalyzer:
    def _make_experiment(self, precursor_mzs_rts):
        """Create MSExperiment with MS2 spectra at given (mz, rt) pairs."""
        import pyopenms as oms

        exp = oms.MSExperiment()
        for mz, rt in precursor_mzs_rts:
            spec = oms.MSSpectrum()
            spec.setMSLevel(2)
            spec.setRT(rt)
            mzs = np.array([100.0, 200.0], dtype=np.float64)
            ints = np.array([500.0, 300.0], dtype=np.float64)
            spec.set_peaks([mzs, ints])
            prec = oms.Precursor()
            prec.setMZ(mz)
            prec.setCharge(2)
            spec.setPrecursors([prec])
            exp.addSpectrum(spec)
        return exp

    def test_extract_precursors(self):
        from precursor_recurrence_analyzer import extract_precursors

        exp = self._make_experiment([(500.0, 10.0), (600.0, 20.0)])
        precs = extract_precursors(exp)
        assert len(precs) == 2
        assert precs[0]["precursor_mz"] == 500.0

    def test_find_recurrent(self):
        from precursor_recurrence_analyzer import find_recurrent_precursors

        # Two precursors with same m/z within tolerance, close in RT
        precs = [
            {"spectrum_index": 0, "rt": 10.0, "precursor_mz": 500.0, "charge": 2},
            {"spectrum_index": 1, "rt": 20.0, "precursor_mz": 500.001, "charge": 2},
        ]
        groups = find_recurrent_precursors(precs, mz_tolerance_ppm=10.0, rt_tolerance_sec=30.0)
        recurrent = [g for g in groups if g["is_recurrent"]]
        assert len(recurrent) == 1
        assert recurrent[0]["count"] == 2

    def test_no_recurrence(self):
        from precursor_recurrence_analyzer import find_recurrent_precursors

        # Two precursors far apart in m/z
        precs = [
            {"spectrum_index": 0, "rt": 10.0, "precursor_mz": 500.0, "charge": 2},
            {"spectrum_index": 1, "rt": 20.0, "precursor_mz": 700.0, "charge": 2},
        ]
        groups = find_recurrent_precursors(precs, mz_tolerance_ppm=10.0, rt_tolerance_sec=30.0)
        recurrent = [g for g in groups if g["is_recurrent"]]
        assert len(recurrent) == 0

    def test_summarize(self):
        from precursor_recurrence_analyzer import find_recurrent_precursors, summarize_recurrence

        precs = [
            {"spectrum_index": 0, "rt": 10.0, "precursor_mz": 500.0, "charge": 2},
            {"spectrum_index": 1, "rt": 15.0, "precursor_mz": 500.001, "charge": 2},
            {"spectrum_index": 2, "rt": 100.0, "precursor_mz": 700.0, "charge": 2},
        ]
        groups = find_recurrent_precursors(precs, mz_tolerance_ppm=10.0, rt_tolerance_sec=30.0)
        summary = summarize_recurrence(groups)
        assert summary["total_precursor_groups"] == 2
        assert summary["recurrent_groups"] == 1

    def test_empty_input(self):
        from precursor_recurrence_analyzer import find_recurrent_precursors

        groups = find_recurrent_precursors([])
        assert groups == []

    def test_write_tsv(self, tmp_path):
        from precursor_recurrence_analyzer import write_tsv

        groups = [{"group_id": 1, "mean_mz": 500.0, "mz_range": 0.001, "rt_min": 10.0,
                    "rt_max": 15.0, "rt_span": 5.0, "count": 2, "is_recurrent": True}]
        out = str(tmp_path / "rec.tsv")
        write_tsv(groups, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
