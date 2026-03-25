"""Tests for peptide_mass_fingerprint."""

import pytest

pytest.importorskip("pyopenms")


class TestPeptideMassFingerprint:
    def _write_fasta(self, tmp_path, accession="P12345", sequence="PEPTIDEKAVLIDRACDEFGHIK"):
        """Write a simple FASTA file for testing."""
        import pyopenms as oms

        fasta_path = str(tmp_path / "test.fasta")
        entry = oms.FASTAEntry()
        entry.identifier = accession
        entry.sequence = sequence
        oms.FASTAFile().store(fasta_path, [entry])
        return fasta_path

    def test_load_fasta(self, tmp_path):
        from peptide_mass_fingerprint import load_fasta

        fasta_path = self._write_fasta(tmp_path)
        proteins = load_fasta(fasta_path)
        assert "P12345" in proteins
        assert len(proteins["P12345"]) > 0

    def test_generate_fingerprint(self):
        from peptide_mass_fingerprint import generate_fingerprint

        protein = "PEPTIDEKAVLIDRACDEFGHIK"
        fp = generate_fingerprint(protein, enzyme="Trypsin", missed_cleavages=1, min_mass=0, max_mass=10000)
        assert len(fp) > 0
        assert all("monoisotopic_mass" in p for p in fp)
        assert all("sequence" in p for p in fp)

    def test_fingerprint_sorted_by_mass(self):
        from peptide_mass_fingerprint import generate_fingerprint

        fp = generate_fingerprint("PEPTIDEKAVLIDRACDEFGHIK", min_mass=0, max_mass=10000)
        masses = [p["monoisotopic_mass"] for p in fp]
        assert masses == sorted(masses)

    def test_match_fingerprint(self):
        from peptide_mass_fingerprint import generate_fingerprint, match_fingerprint

        fp = generate_fingerprint("PEPTIDEKAVLIDRACDEFGHIK", min_mass=0, max_mass=10000)
        # Use exact theoretical mass as observed
        observed = [fp[0]["monoisotopic_mass"]]
        matches = match_fingerprint(fp, observed, tolerance_ppm=10.0)
        assert len(matches) >= 1
        assert matches[0]["ppm_error"] < 1.0

    def test_match_no_hit(self):
        from peptide_mass_fingerprint import generate_fingerprint, match_fingerprint

        fp = generate_fingerprint("PEPTIDEKAVLIDRACDEFGHIK", min_mass=0, max_mass=10000)
        matches = match_fingerprint(fp, [99999.0], tolerance_ppm=10.0)
        assert len(matches) == 0

    def test_mz_values(self):
        from peptide_mass_fingerprint import PROTON, generate_fingerprint

        fp = generate_fingerprint("PEPTIDEKAVLIDRACDEFGHIK", min_mass=0, max_mass=10000)
        for p in fp:
            expected_mz1 = (p["monoisotopic_mass"] + PROTON) / 1
            assert p["mz_1"] == pytest.approx(expected_mz1, abs=1e-4)

    def test_write_tsv(self, tmp_path):
        from peptide_mass_fingerprint import generate_fingerprint, write_tsv

        fp = generate_fingerprint("PEPTIDEKAVLIDRACDEFGHIK", min_mass=0, max_mass=10000)
        out = str(tmp_path / "fp.tsv")
        write_tsv(fp, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) > 1
        assert "sequence" in lines[0]
