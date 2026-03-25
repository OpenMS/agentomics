"""Tests for spectrum_entropy_calculator."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestSpectrumEntropyCalculator:
    def test_uniform_entropy(self):
        from spectrum_entropy_calculator import spectral_entropy

        # Uniform distribution should have entropy = 1.0
        intensities = [100.0, 100.0, 100.0, 100.0]
        entropy = spectral_entropy(intensities)
        assert abs(entropy - 1.0) < 0.001

    def test_single_peak_entropy(self):
        from spectrum_entropy_calculator import spectral_entropy

        # Single peak should have entropy = 0
        assert spectral_entropy([1000.0]) == 0.0

    def test_empty_entropy(self):
        from spectrum_entropy_calculator import spectral_entropy

        assert spectral_entropy([]) == 0.0

    def test_dominated_entropy(self):
        from spectrum_entropy_calculator import spectral_entropy

        # One dominant peak should have low entropy
        intensities = [10000.0, 1.0, 1.0, 1.0]
        entropy = spectral_entropy(intensities)
        assert 0.0 < entropy < 0.5

    def test_compute_from_mzml(self):
        from spectrum_entropy_calculator import compute_spectrum_entropies, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, n_ms2=5)
            results = compute_spectrum_entropies(mzml_path, ms_level=2)
            assert len(results) == 5
            for r in results:
                assert 0.0 <= r["entropy"] <= 1.0

    def test_result_keys(self):
        from spectrum_entropy_calculator import compute_spectrum_entropies, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, n_ms2=1)
            results = compute_spectrum_entropies(mzml_path)
            for r in results:
                assert "scan_index" in r
                assert "rt" in r
                assert "n_peaks" in r
                assert "entropy" in r
                assert "precursor_mz" in r

    def test_write_tsv(self):
        from spectrum_entropy_calculator import compute_spectrum_entropies, create_synthetic_mzml, write_tsv

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, n_ms2=3)
            results = compute_spectrum_entropies(mzml_path)
            out = os.path.join(tmpdir, "entropy.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
            with open(out) as f:
                lines = f.readlines()
            assert len(lines) == 4  # header + 3 data lines
