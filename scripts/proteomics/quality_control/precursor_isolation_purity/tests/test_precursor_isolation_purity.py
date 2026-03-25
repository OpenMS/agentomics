"""Tests for precursor_isolation_purity."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestPrecursorIsolationPurity:
    def _make_experiment(self):
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()

        # MS1 with several peaks
        ms1 = oms.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(10.0)
        mzs = np.array([499.0, 500.0, 500.5, 501.0, 502.0], dtype=np.float64)
        ints = np.array([100.0, 1000.0, 200.0, 100.0, 50.0], dtype=np.float64)
        ms1.set_peaks([mzs, ints])
        exp.addSpectrum(ms1)

        # MS2 targeting 500.0
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(10.5)
        prec = oms.Precursor()
        prec.setMZ(500.0)
        ms2.setPrecursors([prec])
        ms2.set_peaks([np.array([200.0], dtype=np.float64),
                       np.array([500.0], dtype=np.float64)])
        exp.addSpectrum(ms2)

        return exp

    def test_estimate_purity(self):
        import numpy as np
        import pyopenms as oms
        from precursor_isolation_purity import estimate_purity

        ms1 = oms.MSSpectrum()
        mzs = np.array([499.0, 500.0, 501.0], dtype=np.float64)
        ints = np.array([100.0, 1000.0, 100.0], dtype=np.float64)
        ms1.set_peaks([mzs, ints])

        purity = estimate_purity(ms1, 500.0, isolation_width=2.0)
        # 500.0 has 1000, window has 499+500+501=1200
        assert purity > 0.8

    def test_compute_all_purities(self):
        from precursor_isolation_purity import compute_all_purities

        exp = self._make_experiment()
        purities = compute_all_purities(exp)
        assert len(purities) == 1
        assert purities[0]["precursor_mz"] == 500.0
        assert purities[0]["purity"] > 0.0

    def test_empty_ms1(self):
        import numpy as np
        import pyopenms as oms
        from precursor_isolation_purity import estimate_purity

        ms1 = oms.MSSpectrum()
        ms1.set_peaks([np.array([], dtype=np.float64), np.array([], dtype=np.float64)])
        purity = estimate_purity(ms1, 500.0)
        assert purity == 0.0

    def test_purity_range(self):
        from precursor_isolation_purity import compute_all_purities

        exp = self._make_experiment()
        purities = compute_all_purities(exp)
        for p in purities:
            assert 0.0 <= p["purity"] <= 1.0
