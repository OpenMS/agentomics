"""Tests for peptide_detectability_predictor."""

import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestPeptideDetectabilityPredictor:
    def test_basic_score(self):
        from peptide_detectability_predictor import calculate_detectability_score

        result = calculate_detectability_score("PEPTIDEK")
        assert 0 <= result["detectability_score"] <= 1
        assert result["length"] == 8
        assert result["monoisotopic_mass"] > 0

    def test_good_peptide_scores_high(self):
        from peptide_detectability_predictor import calculate_detectability_score

        # Typical tryptic peptide, good length, has K
        good = calculate_detectability_score("AEFGHLPQR")
        # Very short peptide
        short = calculate_detectability_score("AK")
        assert good["detectability_score"] > short["detectability_score"]

    def test_problematic_residues_lower_score(self):
        from peptide_detectability_predictor import calculate_detectability_score

        normal = calculate_detectability_score("PEPTIDEK")
        problem = calculate_detectability_score("MWMWMWMK")
        assert normal["problem_residue_score"] >= problem["problem_residue_score"]

    def test_predict_from_fasta(self):
        import pyopenms as oms
        from peptide_detectability_predictor import predict_from_fasta

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = f"{tmpdir}/test.fasta"
            entries = []
            e1 = oms.FASTAEntry()
            e1.identifier = "PROT1"
            e1.sequence = "MSPEPTIDEKAAANOTHERPEPTIDER"
            entries.append(e1)
            oms.FASTAFile().store(fasta_path, entries)

            results = predict_from_fasta(fasta_path, "Trypsin")
            assert len(results) > 0
            # Results should be sorted by detectability_score descending
            scores = [r["detectability_score"] for r in results]
            assert scores == sorted(scores, reverse=True)

    def test_score_has_all_features(self):
        from peptide_detectability_predictor import calculate_detectability_score

        result = calculate_detectability_score("PEPTIDEK")
        assert "length_score" in result
        assert "hydrophobicity_score" in result
        assert "mass_score" in result
        assert "problem_residue_score" in result
        assert "basic_residue_score" in result
