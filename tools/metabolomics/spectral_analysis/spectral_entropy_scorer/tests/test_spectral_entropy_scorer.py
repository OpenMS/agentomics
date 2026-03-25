"""Tests for spectral_entropy_scorer."""

import math
import os
import tempfile

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestNormalizeIntensities:
    def test_basic(self):
        from spectral_entropy_scorer import normalize_intensities

        result = normalize_intensities([100, 200, 300])
        assert sum(result) == pytest.approx(1.0)
        assert result[0] == pytest.approx(1 / 6)

    def test_zero_total(self):
        from spectral_entropy_scorer import normalize_intensities

        result = normalize_intensities([0, 0, 0])
        assert result == [0.0, 0.0, 0.0]


@requires_pyopenms
class TestSpectralEntropy:
    def test_single_peak(self):
        """Single peak spectrum has entropy 0."""
        from spectral_entropy_scorer import spectral_entropy

        assert spectral_entropy([100.0], [1000.0]) == pytest.approx(0.0)

    def test_uniform_distribution(self):
        """Uniform distribution has maximum entropy."""
        from spectral_entropy_scorer import spectral_entropy

        mzs = [100.0, 200.0, 300.0, 400.0]
        ints = [100.0, 100.0, 100.0, 100.0]
        entropy = spectral_entropy(mzs, ints)
        assert entropy == pytest.approx(math.log(4), abs=0.001)

    def test_entropy_increases_with_peaks(self):
        from spectral_entropy_scorer import spectral_entropy

        e2 = spectral_entropy([100.0, 200.0], [100.0, 100.0])
        e4 = spectral_entropy([100.0, 200.0, 300.0, 400.0], [100.0, 100.0, 100.0, 100.0])
        assert e4 > e2


@requires_pyopenms
class TestMatchPeaks:
    def test_perfect_match(self):
        from spectral_entropy_scorer import match_peaks

        mzs, a, b = match_peaks(
            [100.0, 200.0], [500, 300],
            [100.01, 200.01], [400, 200],
            tolerance=0.02,
        )
        assert len(mzs) == 2
        assert a == [500, 300]
        assert b == [400, 200]

    def test_no_match(self):
        from spectral_entropy_scorer import match_peaks

        mzs, a, b = match_peaks(
            [100.0], [500],
            [200.0], [400],
            tolerance=0.02,
        )
        assert len(mzs) == 2
        assert 0.0 in a or 0.0 in b

    def test_partial_match(self):
        from spectral_entropy_scorer import match_peaks

        mzs, a, b = match_peaks(
            [100.0, 200.0], [500, 300],
            [100.01], [400],
            tolerance=0.02,
        )
        assert len(mzs) == 2
        # 200.0 has no match in b
        assert b[1] == 0.0


@requires_pyopenms
class TestEntropySimilarity:
    def test_identical_spectra(self):
        from spectral_entropy_scorer import entropy_similarity

        score = entropy_similarity(
            [100.0, 200.0, 300.0], [500, 300, 200],
            [100.0, 200.0, 300.0], [500, 300, 200],
            tolerance=0.02,
        )
        assert score == pytest.approx(1.0, abs=0.01)

    def test_completely_different(self):
        from spectral_entropy_scorer import entropy_similarity

        score = entropy_similarity(
            [100.0], [1000],
            [500.0], [1000],
            tolerance=0.02,
        )
        assert score < 0.5

    def test_empty_spectrum(self):
        from spectral_entropy_scorer import entropy_similarity

        assert entropy_similarity([], [], [100.0], [1000], tolerance=0.02) == 0.0

    def test_score_range(self):
        from spectral_entropy_scorer import entropy_similarity

        score = entropy_similarity(
            [100.0, 200.0], [500, 300],
            [100.01, 300.0], [400, 200],
            tolerance=0.02,
        )
        assert 0.0 <= score <= 1.0


@requires_pyopenms
class TestReadPeaksFile:
    def test_roundtrip(self):
        import csv

        from spectral_entropy_scorer import read_peaks_file

        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False, newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["spectrum_id", "mz", "intensity"], delimiter="\t")
            writer.writeheader()
            writer.writerow({"spectrum_id": "S1", "mz": "100.0", "intensity": "500"})
            writer.writerow({"spectrum_id": "S1", "mz": "200.0", "intensity": "300"})
            writer.writerow({"spectrum_id": "S2", "mz": "150.0", "intensity": "1000"})
            path = f.name

        try:
            spectra = read_peaks_file(path)
            assert len(spectra) == 2
            s1 = next(s for s in spectra if s["spectrum_id"] == "S1")
            assert len(s1["mzs"]) == 2
            assert s1["mzs"][0] == 100.0
        finally:
            os.unlink(path)
