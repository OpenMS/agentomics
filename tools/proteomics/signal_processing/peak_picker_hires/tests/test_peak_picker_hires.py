"""Tests for peak_picker_hires."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestPeakPickerHiRes:
    def test_pick_peaks_returns_spectra_count(self):
        from peak_picker_hires import create_synthetic_profile_mzml, pick_peaks

        with tempfile.TemporaryDirectory() as tmp:
            profile_path = os.path.join(tmp, "profile.mzML")
            centroid_path = os.path.join(tmp, "centroid.mzML")

            create_synthetic_profile_mzml(profile_path)
            count = pick_peaks(profile_path, centroid_path, signal_to_noise=1.0)

            assert count >= 1

    def test_output_has_fewer_points(self):
        import pyopenms as oms
        from peak_picker_hires import create_synthetic_profile_mzml, pick_peaks

        with tempfile.TemporaryDirectory() as tmp:
            profile_path = os.path.join(tmp, "profile.mzML")
            centroid_path = os.path.join(tmp, "centroid.mzML")

            create_synthetic_profile_mzml(profile_path)

            # Count points in profile
            exp_in = oms.MSExperiment()
            oms.MzMLFile().load(profile_path, exp_in)
            profile_points = sum(s.size() for s in exp_in.getSpectra())

            pick_peaks(profile_path, centroid_path, signal_to_noise=1.0)

            # Count points in centroided
            exp_out = oms.MSExperiment()
            oms.MzMLFile().load(centroid_path, exp_out)
            centroid_points = sum(s.size() for s in exp_out.getSpectra())

            assert centroid_points < profile_points, (
                f"Expected fewer centroid points ({centroid_points}) "
                f"than profile points ({profile_points})"
            )

    def test_peaks_at_expected_mz(self):
        import pyopenms as oms
        from peak_picker_hires import create_synthetic_profile_mzml, pick_peaks

        with tempfile.TemporaryDirectory() as tmp:
            profile_path = os.path.join(tmp, "profile.mzML")
            centroid_path = os.path.join(tmp, "centroid.mzML")

            create_synthetic_profile_mzml(profile_path)
            pick_peaks(profile_path, centroid_path, signal_to_noise=1.0)

            exp_out = oms.MSExperiment()
            oms.MzMLFile().load(centroid_path, exp_out)

            spec = exp_out.getSpectrum(0)
            mz_values = [spec[i].getMZ() for i in range(spec.size())]

            # Expected peak centers at 400, 500, 600
            expected = [400.0, 500.0, 600.0]
            for expected_mz in expected:
                closest = min(mz_values, key=lambda x: abs(x - expected_mz))
                assert abs(closest - expected_mz) < 0.5, (
                    f"No peak found near m/z {expected_mz}, closest was {closest}"
                )
