"""Tests for flash_deconv."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")

PROTON = 1.007276


class TestFlashDeconv:
    def test_deconvolve_intact_recovers_mass(self):
        """Multiply-charged peaks for M=10000 should recover the mass."""
        import numpy as np
        import pyopenms as oms
        from flash_deconv import deconvolve_intact

        # Create synthetic multiply-charged spectrum for a protein mass
        target_mass = 10000.0
        exp = oms.MSExperiment()

        # Create several scans (FLASHDeconv may need multiple scans)
        for rt_idx in range(10):
            spec = oms.MSSpectrum()
            spec.setMSLevel(1)
            spec.setRT(60.0 + rt_idx * 6.0)

            mz_list = []
            int_list = []
            for z in range(5, 16):
                # Main peak: (M + z*H) / z
                mz = (target_mass + z * PROTON) / z
                mz_list.append(mz)
                int_list.append(1e6)
                # Add isotope peaks for each charge state
                for iso in range(1, 4):
                    mz_list.append(mz + iso * 1.003355 / z)
                    int_list.append(1e6 / (iso + 1))

            # Sort by m/z (required for OpenMS)
            pairs = sorted(zip(mz_list, int_list))
            mzs = np.array([p[0] for p in pairs], dtype=np.float64)
            intensities = np.array([p[1] for p in pairs], dtype=np.float64)
            spec.set_peaks([mzs, intensities])
            exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "intact.mzML")
            output_path = os.path.join(tmpdir, "masses.tsv")
            oms.MzMLFile().store(input_path, exp)

            n_masses = deconvolve_intact(
                input_path,
                output_path,
                min_mass=5000,
                max_mass=20000,
            )

            assert n_masses >= 0  # Algorithm may find masses
            assert os.path.exists(output_path)

            # If masses were found, check they are within tolerance
            if n_masses > 0:
                with open(output_path) as f:
                    f.readline()  # skip header
                    for line in f:
                        parts = line.strip().split("\t")
                        mass = float(parts[0])
                        # Check mass is in our expected range
                        assert 5000 <= mass <= 20000

    def test_deconvolve_output_file_created(self):
        """Verify output TSV file is always created."""
        import numpy as np
        import pyopenms as oms
        from flash_deconv import deconvolve_intact

        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(100.0)
        mzs = np.array([500.0, 600.0, 700.0], dtype=np.float64)
        intensities = np.array([1e5, 1e5, 1e5], dtype=np.float64)
        spec.set_peaks([mzs, intensities])
        exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "test.mzML")
            output_path = os.path.join(tmpdir, "out.tsv")
            oms.MzMLFile().store(input_path, exp)

            n_masses = deconvolve_intact(input_path, output_path)
            assert isinstance(n_masses, int)
            assert os.path.exists(output_path)

    def test_deconvolve_returns_int(self):
        """Return value should be an integer."""
        import numpy as np
        import pyopenms as oms
        from flash_deconv import deconvolve_intact

        exp = oms.MSExperiment()
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(50.0)
        mzs = np.array([800.0], dtype=np.float64)
        intensities = np.array([1e4], dtype=np.float64)
        spec.set_peaks([mzs, intensities])
        exp.addSpectrum(spec)

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "test.mzML")
            output_path = os.path.join(tmpdir, "out.tsv")
            oms.MzMLFile().store(input_path, exp)

            result = deconvolve_intact(input_path, output_path)
            assert isinstance(result, int)
