"""Tests for sequence_tag_generator."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestSequenceTagGenerator:
    def test_get_residue_masses(self):
        from sequence_tag_generator import get_residue_masses

        masses = get_residue_masses()
        assert len(masses) == 20  # 20 standard amino acids
        assert "G" in masses
        assert masses["G"] > 50  # Glycine ~57 Da

    def test_match_mass_to_residue(self):
        from sequence_tag_generator import get_residue_masses, match_mass_to_residue

        masses = get_residue_masses()
        # Glycine mass ~57.02
        matches = match_mass_to_residue(57.02, masses, tolerance=0.05)
        assert "G" in matches

    def test_match_no_residue(self):
        from sequence_tag_generator import get_residue_masses, match_mass_to_residue

        masses = get_residue_masses()
        matches = match_mass_to_residue(999.0, masses, tolerance=0.02)
        assert matches == []

    def test_generate_tags_from_known_peaks(self):
        from sequence_tag_generator import generate_tags, get_residue_masses

        # Build peaks that correspond to a known tag
        # Use Glycine (~57.02) mass differences
        masses = get_residue_masses()
        g_mass = masses["G"]
        base = 200.0
        mz_values = [base, base + g_mass, base + 2 * g_mass, base + 3 * g_mass]
        tags = generate_tags(mz_values, tolerance=0.05, min_tag_length=3)
        tag_strings = [t["tag"] for t in tags]
        assert "GGG" in tag_strings

    def test_generate_tags_too_short(self):
        from sequence_tag_generator import generate_tags, get_residue_masses

        masses = get_residue_masses()
        g_mass = masses["G"]
        mz_values = [200.0, 200.0 + g_mass, 200.0 + 2 * g_mass]
        tags = generate_tags(mz_values, tolerance=0.05, min_tag_length=3)
        # Only 2 mass differences -> max tag length 2, below threshold
        assert all(t["length"] >= 3 for t in tags) or len(tags) == 0

    def test_generate_tags_empty_input(self):
        from sequence_tag_generator import generate_tags

        tags = generate_tags([], tolerance=0.02, min_tag_length=3)
        assert tags == []

    def test_write_tsv(self, tmp_path):
        from sequence_tag_generator import write_tsv

        tags = [{"tag": "GGG", "length": 3, "end_mz": 371.06}]
        out = str(tmp_path / "tags.tsv")
        write_tsv(tags, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 2
        assert "tag" in lines[0]
