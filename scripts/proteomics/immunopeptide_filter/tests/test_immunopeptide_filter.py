"""Tests for immunopeptide_filter."""

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestImmunopeptideFilter:
    def test_filter_by_length_class_i(self):
        from immunopeptide_filter import filter_peptides

        peptides = ["ACDEFGHIK", "ACDE", "ACDEFGHIKL", "ACDEFGHIKLMNPQR"]
        results = filter_peptides(peptides, min_length=8, max_length=11)
        seqs = [r["sequence"] for r in results]
        assert "ACDEFGHIK" in seqs  # 9 aa
        assert "ACDEFGHIKL" in seqs  # 10 aa
        assert "ACDE" not in seqs  # 4 aa - too short
        assert "ACDEFGHIKLMNPQR" not in seqs  # 15 aa - too long

    def test_filter_by_length_class_ii(self):
        from immunopeptide_filter import filter_peptides

        peptides = ["ACDEFGHIK", "ACDEFGHIKLMNPQR"]
        results = filter_peptides(peptides, min_length=13, max_length=25)
        seqs = [r["sequence"] for r in results]
        assert "ACDEFGHIK" not in seqs
        assert "ACDEFGHIKLMNPQR" in seqs

    def test_filter_with_motif(self):
        from immunopeptide_filter import filter_peptides

        peptides = ["ALDEFGHIK", "AXDEFGHIK"]
        results = filter_peptides(peptides, min_length=8, max_length=11, motif_pattern="^A[LIV]")
        seqs = [r["sequence"] for r in results]
        assert "ALDEFGHIK" in seqs
        assert "AXDEFGHIK" not in seqs

    def test_result_has_mass(self):
        from immunopeptide_filter import filter_peptides

        results = filter_peptides(["ACDEFGHIK"], min_length=8, max_length=11)
        assert len(results) == 1
        assert results[0]["monoisotopic_mass"] > 0
        assert results[0]["length"] == 9

    def test_empty_input(self):
        from immunopeptide_filter import filter_peptides

        results = filter_peptides([], min_length=8, max_length=11)
        assert results == []

    def test_parse_length_range(self):
        from immunopeptide_filter import parse_length_range

        assert parse_length_range("8-11") == (8, 11)
        assert parse_length_range("13-25") == (13, 25)

    def test_parse_length_range_invalid(self):
        from immunopeptide_filter import parse_length_range

        with pytest.raises(ValueError):
            parse_length_range("invalid")

    def test_read_and_write_tsv(self, tmp_path):
        from immunopeptide_filter import filter_peptides, read_peptides_from_tsv, write_tsv

        # Create input file
        inp = str(tmp_path / "input.tsv")
        with open(inp, "w") as fh:
            fh.write("sequence\n")
            fh.write("ACDEFGHIK\n")
            fh.write("ACDE\n")

        peptides = read_peptides_from_tsv(inp)
        assert len(peptides) == 2

        results = filter_peptides(peptides, min_length=8, max_length=11)
        out = str(tmp_path / "output.tsv")
        write_tsv(results, out)

        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2  # header + 1 passing peptide
