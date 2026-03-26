"""Tests for deisotoper."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestDeisotoper:
    def test_deisotope_reduces_isotope_envelope(self):
        """Synthetic spectrum with isotope envelope should be reduced."""
        import numpy as np
        import pyopenms as oms
        from deisotoper import deisotope

        # Create a synthetic mzML with one spectrum containing an isotope envelope
        # Monoisotopic peak at m/z=500.0 with isotope peaks at +1.003, +2.006
        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(100.0)

        mono_mz = 500.0
        delta = 1.003355  # neutron mass difference
        mzs = np.array(
            [mono_mz, mono_mz + delta, mono_mz + 2 * delta],
            dtype=np.float64,
        )
        intensities = np.array([1000.0, 600.0, 200.0], dtype=np.float64)
        spec.set_peaks([mzs, intensities])
        exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.mzML")
            output_path = os.path.join(tmpdir, "output.mzML")
            oms.MzMLFile().store(input_path, exp)

            count = deisotope(input_path, output_path)
            assert count == 1  # one spectrum processed

            # Load result and check peaks were reduced
            result_exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, result_exp)
            result_spec = result_exp[0]
            result_mzs, _ = result_spec.get_peaks()

            # The deisotoper should reduce the 3 peaks to fewer peaks
            # (ideally just the monoisotopic peak)
            assert len(result_mzs) < len(mzs)

    def test_deisotope_preserves_non_isotopic_peaks(self):
        """Peaks that are not part of an isotope envelope should be preserved."""
        import numpy as np
        import pyopenms as oms
        from deisotoper import deisotope

        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(100.0)

        # Two well-separated peaks that are not isotope-related
        mzs = np.array([200.0, 500.0], dtype=np.float64)
        intensities = np.array([1000.0, 800.0], dtype=np.float64)
        spec.set_peaks([mzs, intensities])
        exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.mzML")
            output_path = os.path.join(tmpdir, "output.mzML")
            oms.MzMLFile().store(input_path, exp)

            deisotope(input_path, output_path)

            result_exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, result_exp)
            result_mzs, _ = result_exp[0].get_peaks()

            # Both non-isotopic peaks should remain
            assert len(result_mzs) == 2

    def test_deisotope_output_file_created(self):
        """Verify output mzML file is created."""
        import numpy as np
        import pyopenms as oms
        from deisotoper import deisotope

        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(50.0)
        mzs = np.array([300.0], dtype=np.float64)
        intensities = np.array([500.0], dtype=np.float64)
        spec.set_peaks([mzs, intensities])
        exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.mzML")
            output_path = os.path.join(tmpdir, "output.mzML")
            oms.MzMLFile().store(input_path, exp)

            deisotope(input_path, output_path)
            assert os.path.exists(output_path)

    def test_deisotope_multiple_spectra(self):
        """Process multiple spectra and return correct count."""
        import numpy as np
        import pyopenms as oms
        from deisotoper import deisotope

        exp = oms.MSExperiment()
        for i in range(5):
            spec = oms.MSSpectrum()
            spec.setMSLevel(2)
            spec.setRT(10.0 * i)
            mzs = np.array([300.0 + i], dtype=np.float64)
            intensities = np.array([1000.0], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "input.mzML")
            output_path = os.path.join(tmpdir, "output.mzML")
            oms.MzMLFile().store(input_path, exp)

            count = deisotope(input_path, output_path)
            assert count == 5
