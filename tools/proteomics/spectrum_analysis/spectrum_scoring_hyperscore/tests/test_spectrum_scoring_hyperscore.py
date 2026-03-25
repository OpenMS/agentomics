"""Tests for spectrum_scoring_hyperscore."""

import json
import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestSpectrumScoringHyperscore:
    def _get_theoretical_peaks(self, sequence="PEPTIDEK", charge=1):
        """Helper to get theoretical peaks for testing."""
        import pyopenms as oms

        aa_seq = oms.AASequence.fromString(sequence)
        tsg = oms.TheoreticalSpectrumGenerator()
        param = tsg.getParameters()
        param.setValue("add_b_ions", "true")
        param.setValue("add_y_ions", "true")
        param.setValue("add_metainfo", "true")
        tsg.setParameters(param)
        spec = oms.MSSpectrum()
        tsg.getSpectrum(spec, aa_seq, charge, charge)
        mzs, _ = spec.get_peaks()
        return [float(m) for m in mzs]

    def test_matching_spectrum(self):
        from spectrum_scoring_hyperscore import compute_hyperscore

        theo_mzs = self._get_theoretical_peaks()
        intensities = [1000.0] * len(theo_mzs)
        result = compute_hyperscore(theo_mzs, intensities, "PEPTIDEK", charge=1, tolerance=0.05)
        assert result["hyperscore"] > 0
        assert result["matched_total"] > 0
        assert result["sequence"] == "PEPTIDEK"

    def test_no_match(self):
        from spectrum_scoring_hyperscore import compute_hyperscore

        result = compute_hyperscore([9999.0], [1000.0], "PEPTIDEK", charge=1, tolerance=0.02)
        assert result["hyperscore"] == 0.0
        assert result["matched_total"] == 0

    def test_result_keys(self):
        from spectrum_scoring_hyperscore import compute_hyperscore

        result = compute_hyperscore([100.0], [1000.0], "PEPTIDEK", charge=1)
        expected_keys = {"hyperscore", "matched_b", "matched_y", "matched_total", "dot_product", "sequence", "charge"}
        assert set(result.keys()) == expected_keys

    def test_output_json(self):
        from spectrum_scoring_hyperscore import compute_hyperscore

        result = compute_hyperscore([100.0], [1000.0], "PEPTIDEK", charge=1)
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "score.json")
            with open(out, "w") as f:
                json.dump(result, f)
            assert os.path.exists(out)
            with open(out) as f:
                loaded = json.load(f)
            assert "hyperscore" in loaded
