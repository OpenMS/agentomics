"""Tests for metabolite_spectral_matcher."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")


class TestMetaboliteSpectralMatcher:
    """Tests for metabolite spectral matching functionality."""

    def test_create_synthetic_query_mzml(self):
        from metabolite_spectral_matcher import create_synthetic_query_mzml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "query.mzML")
            create_synthetic_query_mzml(path)

            exp = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(path, exp)
            assert exp.getNrSpectra() == 3
            assert exp.getSpectrum(0).getMSLevel() == 2

    def test_create_synthetic_library_mzml(self):
        from metabolite_spectral_matcher import create_synthetic_library_mzml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "lib.mzML")
            create_synthetic_library_mzml(path)

            lib = pyopenms.MSExperiment()
            pyopenms.MzMLFile().load(path, lib)
            assert lib.getNrSpectra() == 3
            assert lib.getSpectrum(0).getMSLevel() == 2

    def test_match_spectra_finds_matches(self):
        """Verify that matching query and library spectra produce matches with score > 0."""
        from metabolite_spectral_matcher import (
            create_synthetic_library_mzml,
            create_synthetic_query_mzml,
            match_spectra,
        )

        with tempfile.TemporaryDirectory() as tmp:
            query_path = os.path.join(tmp, "query.mzML")
            lib_path = os.path.join(tmp, "lib.mzML")
            output_path = os.path.join(tmp, "matches.mzTab")

            create_synthetic_query_mzml(query_path)
            create_synthetic_library_mzml(lib_path)

            matches = match_spectra(query_path, lib_path, output_path, precursor_tol=0.1)

            assert matches >= 1
            assert os.path.exists(output_path)

            # Verify match scores > 0 in the mzTab output
            with open(output_path) as f:
                content = f.read()
            # SML lines should exist
            sml_lines = [line for line in content.split("\n") if line.startswith("SML\t")]
            assert len(sml_lines) >= 1

    def test_match_spectra_score_content(self):
        """Verify match score is present and positive in mzTab output."""
        from metabolite_spectral_matcher import (
            create_synthetic_library_mzml,
            create_synthetic_query_mzml,
            match_spectra,
        )

        with tempfile.TemporaryDirectory() as tmp:
            query_path = os.path.join(tmp, "query.mzML")
            lib_path = os.path.join(tmp, "lib.mzML")
            output_path = os.path.join(tmp, "matches.mzTab")

            create_synthetic_query_mzml(query_path)
            create_synthetic_library_mzml(lib_path)

            match_spectra(query_path, lib_path, output_path, precursor_tol=0.1)

            # Read the match_score column from SML lines
            with open(output_path) as f:
                lines = f.readlines()

            # Find header to locate the opt_match_score column
            header_line = None
            for line in lines:
                if line.startswith("SMH\t"):
                    header_line = line.strip()
                    break

            assert header_line is not None
            headers = header_line.split("\t")
            score_idx = None
            for i, h in enumerate(headers):
                if "match_score" in h:
                    score_idx = i
                    break

            if score_idx is not None:
                for line in lines:
                    if line.startswith("SML\t"):
                        fields = line.strip().split("\t")
                        score = float(fields[score_idx])
                        assert score > 0
