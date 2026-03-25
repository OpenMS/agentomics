"""Tests for rna_digest."""

import csv

import pytest

pytest.importorskip("pyopenms")

from rna_digest import _calculate_fragment_mass, digest_rna


class TestRnaDigest:
    def test_rnase_t1_basic(self):
        """RNase T1 cleaves after G."""
        fragments = digest_rna("AAUGCAAUGG", "RNase_T1", missed_cleavages=0)
        seqs = [f["fragment"] for f in fragments]
        # Should cleave after G at positions 4 (AAUG) and 9 (CAAUG), last G is terminal
        assert all(isinstance(f["mass"], float) for f in fragments)
        assert all(f["missed_cleavages"] == 0 for f in fragments)
        # Reconstruct sequence
        assert "".join(seqs) == "AAUGCAAUGG"

    def test_rnase_a_basic(self):
        """RNase A cleaves after C and U."""
        fragments = digest_rna("AAUGCAAUGG", "RNase_A", missed_cleavages=0)
        seqs = [f["fragment"] for f in fragments]
        assert "".join(seqs) == "AAUGCAAUGG"

    def test_missed_cleavages(self):
        fragments_0 = digest_rna("AAUGCAAUGG", "RNase_T1", missed_cleavages=0)
        fragments_1 = digest_rna("AAUGCAAUGG", "RNase_T1", missed_cleavages=1)
        assert len(fragments_1) > len(fragments_0)

    def test_no_cleavage_site(self):
        """If no cleavage site, return whole sequence."""
        fragments = digest_rna("AAAA", "RNase_T1", missed_cleavages=0)
        assert len(fragments) == 1
        assert fragments[0]["fragment"] == "AAAA"

    def test_unknown_enzyme(self):
        with pytest.raises(ValueError, match="Unknown enzyme"):
            digest_rna("AAUGC", "FakeEnzyme")

    def test_invalid_nucleotide(self):
        with pytest.raises(ValueError, match="Invalid RNA nucleotide"):
            digest_rna("AATGC", "RNase_T1")

    def test_fragment_mass_positive(self):
        mass = _calculate_fragment_mass("AAUGC")
        assert mass > 0

    def test_start_end_positions(self):
        fragments = digest_rna("AAUGCAAUGG", "RNase_T1", missed_cleavages=0)
        for f in fragments:
            assert f["start"] >= 1
            assert f["end"] <= 10
            assert f["end"] >= f["start"]

    def test_output_file(self, tmp_path):
        """Test writing to TSV file via direct function call."""
        fragments = digest_rna("AAUGCAAUGG", "RNase_T1")
        outfile = str(tmp_path / "fragments.tsv")
        with open(outfile, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["fragment", "start", "end", "missed_cleavages", "mass"],
                delimiter="\t",
            )
            writer.writeheader()
            writer.writerows(fragments)

        with open(outfile) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            rows = list(reader)
        assert len(rows) == len(fragments)
