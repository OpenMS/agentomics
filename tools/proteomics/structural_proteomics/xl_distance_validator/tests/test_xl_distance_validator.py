"""Tests for xl_distance_validator."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestXlDistanceValidator:
    def _create_pdb(self, tmpdir):
        """Create a minimal PDB file with CA atoms."""
        pdb_path = os.path.join(tmpdir, "structure.pdb")
        # PDB ATOM format: columns are fixed-width
        # ATOM   serial name altLoc resName chainID resSeq    x       y       z
        lines = [
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C  ",
            "ATOM      2  CA  LYS A   5      10.000   0.000   0.000  1.00  0.00           C  ",
            "ATOM      3  CA  ALA A  10      50.000   0.000   0.000  1.00  0.00           C  ",
            "ATOM      4  CA  ALA B   1       5.000   0.000   0.000  1.00  0.00           C  ",
            "END",
        ]
        with open(pdb_path, "w") as f:
            f.write("\n".join(lines) + "\n")
        return pdb_path

    def test_parse_pdb_ca_atoms(self):
        from xl_distance_validator import parse_pdb_ca_atoms
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = self._create_pdb(tmpdir)
            ca_atoms = parse_pdb_ca_atoms(pdb_path)
            assert ("A", 1) in ca_atoms
            assert ("A", 5) in ca_atoms
            assert ("B", 1) in ca_atoms
            assert ca_atoms[("A", 1)] == (0.0, 0.0, 0.0)

    def test_euclidean_distance(self):
        from xl_distance_validator import euclidean_distance
        d = euclidean_distance((0, 0, 0), (3, 4, 0))
        assert abs(d - 5.0) < 1e-6

    def test_euclidean_distance_3d(self):
        from xl_distance_validator import euclidean_distance
        d = euclidean_distance((1, 2, 3), (4, 6, 3))
        assert abs(d - 5.0) < 1e-6

    def test_validate_sequence(self):
        from xl_distance_validator import validate_sequence
        assert validate_sequence("PEPTIDEK") is True

    def test_validate_crosslink_satisfied(self):
        from xl_distance_validator import validate_crosslink
        ca_atoms = {("A", 1): (0.0, 0.0, 0.0), ("A", 5): (10.0, 0.0, 0.0)}
        result = validate_crosslink("A", 1, "A", 5, ca_atoms, max_distance=30.0)
        assert result["satisfied"] == "YES"
        assert abs(result["distance"] - 10.0) < 0.1

    def test_validate_crosslink_violated(self):
        from xl_distance_validator import validate_crosslink
        ca_atoms = {("A", 1): (0.0, 0.0, 0.0), ("A", 10): (50.0, 0.0, 0.0)}
        result = validate_crosslink("A", 1, "A", 10, ca_atoms, max_distance=30.0)
        assert result["satisfied"] == "NO"

    def test_validate_crosslink_missing_residue(self):
        from xl_distance_validator import validate_crosslink
        ca_atoms = {("A", 1): (0.0, 0.0, 0.0)}
        result = validate_crosslink("A", 1, "A", 99, ca_atoms, max_distance=30.0)
        assert result["satisfied"] == "UNKNOWN"

    def test_validate_crosslinks_batch(self):
        from xl_distance_validator import validate_crosslinks
        ca_atoms = {
            ("A", 1): (0.0, 0.0, 0.0),
            ("A", 5): (10.0, 0.0, 0.0),
            ("A", 10): (50.0, 0.0, 0.0),
        }
        crosslinks = [
            {"peptide1": "AALK", "peptide2": "KDEF", "chain1": "A", "residue1": "1",
             "chain2": "A", "residue2": "5"},
            {"peptide1": "AALK", "peptide2": "GHIK", "chain1": "A", "residue1": "1",
             "chain2": "A", "residue2": "10"},
        ]
        results = validate_crosslinks(crosslinks, ca_atoms, max_distance=30.0)
        assert len(results) == 2
        assert results[0]["satisfied"] == "YES"
        assert results[1]["satisfied"] == "NO"

    def test_write_output(self):
        from xl_distance_validator import write_output
        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = os.path.join(tmpdir, "distances.tsv")
            results = [{"chain1": "A", "residue1": 1, "distance": 10.0, "satisfied": "YES"}]
            write_output(output_path, results)
            assert os.path.exists(output_path)

    def test_full_pipeline(self):
        from xl_distance_validator import parse_pdb_ca_atoms, read_crosslinks, validate_crosslinks, write_output

        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = self._create_pdb(tmpdir)
            xl_path = os.path.join(tmpdir, "links.tsv")
            with open(xl_path, "w") as f:
                f.write("peptide1\tpeptide2\tchain1\tresidue1\tchain2\tresidue2\n")
                f.write("AALK\tKDEF\tA\t1\tA\t5\n")

            ca_atoms = parse_pdb_ca_atoms(pdb_path)
            crosslinks = read_crosslinks(xl_path)
            results = validate_crosslinks(crosslinks, ca_atoms, 30.0)
            output_path = os.path.join(tmpdir, "out.tsv")
            write_output(output_path, results)
            assert os.path.exists(output_path)
            assert results[0]["satisfied"] == "YES"
