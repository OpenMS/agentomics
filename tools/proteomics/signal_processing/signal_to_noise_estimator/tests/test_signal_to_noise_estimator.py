"""Tests for signal_to_noise_estimator."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestSignalToNoiseEstimator:
    def test_estimate_sn_returns_spectra_count(self):
        from signal_to_noise_estimator import create_synthetic_mzml_with_noise, estimate_sn

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            out_path = os.path.join(tmp, "sn_report.tsv")

            create_synthetic_mzml_with_noise(mzml_path)
            count = estimate_sn(mzml_path, out_path)

            assert count == 1

    def test_strong_peak_has_high_sn(self):
        from signal_to_noise_estimator import create_synthetic_mzml_with_noise, estimate_sn

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            out_path = os.path.join(tmp, "sn_report.tsv")

            create_synthetic_mzml_with_noise(mzml_path)
            estimate_sn(mzml_path, out_path)

            # Parse the TSV and find the peak at ~500 m/z
            max_sn = 0.0
            with open(out_path) as f:
                header = f.readline()
                assert "sn_ratio" in header
                for line in f:
                    parts = line.strip().split("\t")
                    mz = float(parts[1])
                    sn = float(parts[3])
                    if abs(mz - 500.0) < 1.0:
                        max_sn = max(max_sn, sn)

            assert max_sn > 10.0, f"Expected S/N > 10 for strong peak, got {max_sn}"

    def test_output_tsv_has_correct_columns(self):
        from signal_to_noise_estimator import create_synthetic_mzml_with_noise, estimate_sn

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            out_path = os.path.join(tmp, "sn_report.tsv")

            create_synthetic_mzml_with_noise(mzml_path)
            estimate_sn(mzml_path, out_path)

            with open(out_path) as f:
                header = f.readline().strip()
                assert header == "spectrum_index\tmz\tintensity\tsn_ratio"
                lines = f.readlines()
                assert len(lines) > 0
