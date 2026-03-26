"""Tests for linear_resampler."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestLinearResampler:
    def test_resample_returns_spectra_count(self):
        from linear_resampler import create_synthetic_irregular_mzml, resample_experiment

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "resampled.mzML")

            create_synthetic_irregular_mzml(input_path)
            count = resample_experiment(input_path, output_path, spacing=0.1)

            assert count == 1

    def test_output_mz_uniformly_spaced(self):
        import pyopenms as oms
        from linear_resampler import create_synthetic_irregular_mzml, resample_experiment

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "resampled.mzML")

            create_synthetic_irregular_mzml(input_path)
            spacing = 0.1
            resample_experiment(input_path, output_path, spacing=spacing)

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            spec = exp.getSpectrum(0)
            mz_values = [spec[i].getMZ() for i in range(spec.size())]

            assert len(mz_values) > 2, "Expected resampled spectrum to have points"

            # Check that consecutive m/z differences are approximately equal
            diffs = [mz_values[i + 1] - mz_values[i] for i in range(len(mz_values) - 1)]
            for diff in diffs:
                assert abs(diff - spacing) < 0.01, (
                    f"Expected uniform spacing of {spacing}, got {diff}"
                )

    def test_different_spacing_values(self):
        import pyopenms as oms
        from linear_resampler import create_synthetic_irregular_mzml, resample_experiment

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")

            create_synthetic_irregular_mzml(input_path)

            for spacing in [0.05, 0.5, 1.0]:
                output_path = os.path.join(tmp, f"resampled_{spacing}.mzML")
                resample_experiment(input_path, output_path, spacing=spacing)

                exp = oms.MSExperiment()
                oms.MzMLFile().load(output_path, exp)

                spec = exp.getSpectrum(0)
                n_points = spec.size()
                assert n_points > 0, f"No points in resampled spectrum (spacing={spacing})"

                mz_values = [spec[i].getMZ() for i in range(spec.size())]
                if len(mz_values) > 1:
                    diff = mz_values[1] - mz_values[0]
                    assert abs(diff - spacing) < 0.01
