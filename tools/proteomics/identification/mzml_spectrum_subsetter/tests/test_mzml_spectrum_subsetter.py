"""Tests for mzml_spectrum_subsetter."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestMzmlSpectrumSubsetter:
    def test_subset_spectra(self):
        import pyopenms as oms
        from mzml_spectrum_subsetter import create_synthetic_mzml, subset_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "subset.mzML")
            create_synthetic_mzml(in_path, n_scans=10)
            count = subset_spectra(in_path, [0, 2, 4], out_path)
            assert count == 3

            exp = oms.MSExperiment()
            oms.MzMLFile().load(out_path, exp)
            assert exp.getNrSpectra() == 3

    def test_out_of_range(self):
        import pyopenms as oms
        from mzml_spectrum_subsetter import create_synthetic_mzml, subset_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "subset.mzML")
            create_synthetic_mzml(in_path, n_scans=5)
            count = subset_spectra(in_path, [0, 100], out_path)
            assert count == 1

            exp = oms.MSExperiment()
            oms.MzMLFile().load(out_path, exp)
            assert exp.getNrSpectra() == 1

    def test_empty_subset(self):
        from mzml_spectrum_subsetter import create_synthetic_mzml, subset_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "subset.mzML")
            create_synthetic_mzml(in_path, n_scans=5)
            count = subset_spectra(in_path, [100, 200], out_path)
            assert count == 0

    def test_preserves_rt(self):
        import pyopenms as oms
        from mzml_spectrum_subsetter import create_synthetic_mzml, subset_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "subset.mzML")
            create_synthetic_mzml(in_path, n_scans=10)
            subset_spectra(in_path, [2], out_path)

            exp = oms.MSExperiment()
            oms.MzMLFile().load(out_path, exp)
            assert exp.getNrSpectra() == 1
            # Scan index 2 has RT = 2 * 5.0 = 10.0
            assert abs(exp.getSpectrum(0).getRT() - 10.0) < 0.1
