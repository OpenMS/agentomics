"""Tests for semi_tryptic_peptide_finder."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestSemiTrypticPeptideFinder:
    def test_fully_tryptic(self):
        from semi_tryptic_peptide_finder import classify_peptide

        # PEPTIDEK ends with K, preceded by protein N-term
        protein = "PEPTIDEKAVLIDR"
        result = classify_peptide("PEPTIDEK", protein, "Trypsin")
        assert result == "fully_tryptic"

    def test_semi_tryptic_n_term(self):
        from semi_tryptic_peptide_finder import classify_peptide

        # Peptide starts after K (tryptic N-term) but does not end with K/R
        protein = "PEPTIDEKAVLIDXYZ"
        result = classify_peptide("AVLID", protein, "Trypsin")
        assert result == "semi_tryptic"

    def test_semi_tryptic_c_term(self):
        from semi_tryptic_peptide_finder import classify_peptide

        # Peptide ends with K (tryptic C-term) but N-term is not after K/R
        protein = "AVLIDPEPTIDEK"
        result = classify_peptide("IDPEPTIDEK", protein, "Trypsin")
        assert result == "semi_tryptic"

    def test_non_tryptic(self):
        from semi_tryptic_peptide_finder import classify_peptide

        # Neither end matches tryptic cleavage
        protein = "AVLIDPEPTIDEGG"
        result = classify_peptide("IDPEPTIDE", protein, "Trypsin")
        assert result == "non_tryptic"

    def test_not_found(self):
        from semi_tryptic_peptide_finder import classify_peptide

        result = classify_peptide("XYZXYZ", "PEPTIDEK", "Trypsin")
        assert result == "not_found"

    def test_classify_against_fasta(self):
        from semi_tryptic_peptide_finder import classify_peptides_against_fasta

        proteins = {"P1": "PEPTIDEKAVLIDR"}
        results = classify_peptides_against_fasta(["PEPTIDEK"], proteins, "Trypsin")
        assert len(results) == 1
        assert results[0]["classification"] == "fully_tryptic"
        assert results[0]["protein"] == "P1"

    def test_digest_protein(self):
        from semi_tryptic_peptide_finder import digest_protein

        peptides = digest_protein("PEPTIDEKAVLIDR", "Trypsin", missed_cleavages=0)
        assert len(peptides) >= 2

    def test_write_tsv(self, tmp_path):
        from semi_tryptic_peptide_finder import write_tsv

        results = [{"sequence": "PEPTIDEK", "length": 8, "classification": "fully_tryptic", "protein": "P1"}]
        out = str(tmp_path / "out.tsv")
        write_tsv(results, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
        assert "fully_tryptic" in lines[1]
