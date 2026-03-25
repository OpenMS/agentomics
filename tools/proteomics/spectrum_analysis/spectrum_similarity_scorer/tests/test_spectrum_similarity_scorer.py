"""Tests for spectrum_similarity_scorer."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")

MGF_TEMPLATE = """BEGIN IONS
TITLE={title}
PEPMASS={pepmass}
CHARGE={charge}+
{peaks}
END IONS
"""


def _write_mgf(path, spectra):
    with open(path, "w") as f:
        for s in spectra:
            peaks = "\n".join(f"{mz} {intensity}" for mz, intensity in zip(s["mzs"], s["intensities"]))
            f.write(MGF_TEMPLATE.format(title=s["title"], pepmass=s["pepmass"], charge=s["charge"], peaks=peaks))


class TestSpectrumSimilarityScorer:
    def test_parse_mgf(self):
        from spectrum_similarity_scorer import parse_mgf

        with tempfile.TemporaryDirectory() as tmpdir:
            mgf_path = os.path.join(tmpdir, "test.mgf")
            _write_mgf(mgf_path, [
                {"title": "spec1", "pepmass": 500.0, "charge": 2,
                 "mzs": [100.0, 200.0, 300.0], "intensities": [1000, 500, 200]},
            ])
            spectra = parse_mgf(mgf_path)
            assert len(spectra) == 1
            assert spectra[0]["title"] == "spec1"
            assert len(spectra[0]["mz_array"]) == 3

    def test_cosine_identical(self):
        from spectrum_similarity_scorer import cosine_similarity, mgf_to_msspectrum

        spec_dict = {"mz_array": [100.0, 200.0, 300.0], "intensity_array": [1000.0, 500.0, 200.0]}
        s1 = mgf_to_msspectrum(spec_dict)
        s2 = mgf_to_msspectrum(spec_dict)
        result = cosine_similarity(s1, s2, tolerance=0.02)
        assert result["score"] > 0.99
        assert result["matched_peaks"] == 3

    def test_cosine_no_match(self):
        from spectrum_similarity_scorer import cosine_similarity, mgf_to_msspectrum

        s1 = mgf_to_msspectrum({"mz_array": [100.0, 200.0], "intensity_array": [1000.0, 500.0]})
        s2 = mgf_to_msspectrum({"mz_array": [500.0, 600.0], "intensity_array": [1000.0, 500.0]})
        result = cosine_similarity(s1, s2, tolerance=0.02)
        assert result["score"] == 0.0
        assert result["matched_peaks"] == 0

    def test_score_spectra(self):
        from spectrum_similarity_scorer import score_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            q_path = os.path.join(tmpdir, "query.mgf")
            l_path = os.path.join(tmpdir, "library.mgf")
            _write_mgf(q_path, [
                {"title": "q1", "pepmass": 500.0, "charge": 2,
                 "mzs": [100.0, 200.0, 300.0], "intensities": [1000, 500, 200]},
            ])
            _write_mgf(l_path, [
                {"title": "l1", "pepmass": 500.0, "charge": 2,
                 "mzs": [100.0, 200.0, 300.0], "intensities": [1000, 500, 200]},
            ])
            results = score_spectra(q_path, l_path, tolerance=0.02)
            assert len(results) == 1
            assert results[0]["query_id"] == "q1"
            assert results[0]["score"] > 0.99

    def test_write_tsv(self):
        from spectrum_similarity_scorer import write_tsv

        results = [{"query_id": "q1", "library_id": "l1", "score": 0.95, "matched_peaks": 5}]
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "scores.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
            with open(out) as f:
                lines = f.readlines()
            assert len(lines) == 2
            assert "query_id" in lines[0]
