"""Tests for peptide_spectral_match_validator."""

import pytest

pytest.importorskip("pyopenms")


class TestPeptideSpectralMatchValidator:
    def test_generate_theoretical_spectrum(self):
        from peptide_spectral_match_validator import generate_theoretical_spectrum

        spec = generate_theoretical_spectrum("PEPTIDEK", charge=1)
        assert spec.size() > 0

    def test_validate_psm_with_matching_spectrum(self):
        from peptide_spectral_match_validator import generate_theoretical_spectrum, validate_psm

        # Use the theoretical spectrum itself as the experimental spectrum -> perfect match
        theo = generate_theoretical_spectrum("PEPTIDEK", charge=1)
        result = validate_psm(theo, "PEPTIDEK", tolerance=0.02, charge=1)
        assert result["coverage"] > 0
        assert result["matched_ions"] > 0
        assert result["theoretical_ions"] > 0

    def test_validate_psm_no_match(self):
        import numpy as np
        import pyopenms as oms
        from peptide_spectral_match_validator import validate_psm

        # Create a spectrum with peaks far from any theoretical ions
        spec = oms.MSSpectrum()
        mzs = np.array([10.0, 20.0, 30.0], dtype=np.float64)
        ints = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        spec.set_peaks([mzs, ints])
        result = validate_psm(spec, "PEPTIDEK", tolerance=0.02, charge=1)
        assert result["matched_ions"] == 0

    def test_validate_psms_batch(self):
        import pyopenms as oms
        from peptide_spectral_match_validator import generate_theoretical_spectrum, validate_psms

        exp = oms.MSExperiment()
        # Add a theoretical spectrum as "experimental"
        theo = generate_theoretical_spectrum("PEPTIDEK", charge=1)
        theo.setMSLevel(2)
        theo.setRT(10.0)
        exp.addSpectrum(theo)

        psms = [{"spectrum_index": 0, "sequence": "PEPTIDEK", "charge": 1}]
        results = validate_psms(exp, psms, tolerance=0.02)
        assert len(results) == 1
        assert results[0]["status"] == "valid"

    def test_invalid_index(self):
        import pyopenms as oms
        from peptide_spectral_match_validator import validate_psms

        exp = oms.MSExperiment()
        psms = [{"spectrum_index": 999, "sequence": "PEPTIDEK", "charge": 1}]
        results = validate_psms(exp, psms, tolerance=0.02)
        assert results[0]["status"] == "invalid_index"

    def test_write_tsv(self, tmp_path):
        from peptide_spectral_match_validator import write_tsv

        results = [{
            "spectrum_index": 0, "sequence": "PEPTIDEK", "theoretical_ions": 14,
            "matched_ions": 10, "coverage": 0.714, "peptide_mass": 927.4, "status": "valid",
        }]
        out = str(tmp_path / "val.tsv")
        write_tsv(results, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
