"""Tests for metapeptide_lca_assigner."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestMetapeptideLcaAssigner:
    def _create_fasta(self, tmpdir, proteins):
        """Helper to create a FASTA file."""
        import pyopenms as oms

        fasta_path = os.path.join(tmpdir, "metadb.fasta")
        entries = []
        for acc, seq in proteins.items():
            entry = oms.FASTAEntry()
            entry.identifier = acc
            entry.sequence = seq
            entries.append(entry)
        fasta_file = oms.FASTAFile()
        fasta_file.store(fasta_path, entries)
        return fasta_path

    def _create_taxonomy(self, tmpdir, taxonomy):
        """Helper to create a taxonomy lineage file."""
        tax_path = os.path.join(tmpdir, "lineage.tsv")
        with open(tax_path, "w") as f:
            f.write("protein\tlineage\n")
            for prot, lineage in taxonomy.items():
                f.write(f"{prot}\t{lineage}\n")
        return tax_path

    def test_load_fasta(self):
        from metapeptide_lca_assigner import load_fasta
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {"P1": "ACDEFGHIK"})
            proteins = load_fasta(fasta_path)
            assert "P1" in proteins

    def test_load_taxonomy(self):
        from metapeptide_lca_assigner import load_taxonomy
        with tempfile.TemporaryDirectory() as tmpdir:
            tax_path = self._create_taxonomy(tmpdir, {
                "P1": "Bacteria;Proteobacteria;Gammaproteobacteria"
            })
            taxonomy = load_taxonomy(tax_path)
            assert "P1" in taxonomy
            assert taxonomy["P1"] == ["Bacteria", "Proteobacteria", "Gammaproteobacteria"]

    def test_strip_modifications(self):
        from metapeptide_lca_assigner import strip_modifications
        assert strip_modifications("PEPTM[147]IDEK") == "PEPTMIDEK"

    def test_find_peptide_proteins(self):
        from metapeptide_lca_assigner import find_peptide_proteins
        proteins = {"P1": "ACDEFGHIK", "P2": "XXXXACDEFYYY"}
        found = find_peptide_proteins("ACDEF", proteins)
        assert found == {"P1", "P2"}

    def test_compute_lca_same_lineage(self):
        from metapeptide_lca_assigner import compute_lca
        lineages = [
            ["Bacteria", "Proteobacteria", "Gamma"],
            ["Bacteria", "Proteobacteria", "Gamma"],
        ]
        lca = compute_lca(lineages)
        assert lca == ["Bacteria", "Proteobacteria", "Gamma"]

    def test_compute_lca_partial_overlap(self):
        from metapeptide_lca_assigner import compute_lca
        lineages = [
            ["Bacteria", "Proteobacteria", "Gammaproteobacteria"],
            ["Bacteria", "Proteobacteria", "Betaproteobacteria"],
        ]
        lca = compute_lca(lineages)
        assert lca == ["Bacteria", "Proteobacteria"]

    def test_compute_lca_no_overlap(self):
        from metapeptide_lca_assigner import compute_lca
        lineages = [
            ["Bacteria", "Proteobacteria"],
            ["Archaea", "Euryarchaeota"],
        ]
        lca = compute_lca(lineages)
        assert lca == []

    def test_compute_lca_single(self):
        from metapeptide_lca_assigner import compute_lca
        lineages = [["Bacteria", "Proteobacteria"]]
        lca = compute_lca(lineages)
        assert lca == ["Bacteria", "Proteobacteria"]

    def test_compute_lca_empty(self):
        from metapeptide_lca_assigner import compute_lca
        assert compute_lca([]) == []

    def test_assign_lca_for_peptide(self):
        from metapeptide_lca_assigner import assign_lca_for_peptide
        proteins = {"P1": "ACDEFGHIK", "P2": "XXXXACDEFYYY"}
        taxonomy = {
            "P1": ["Bacteria", "Proteobacteria", "Gamma"],
            "P2": ["Bacteria", "Proteobacteria", "Beta"],
        }
        result = assign_lca_for_peptide("ACDEF", proteins, taxonomy)
        assert result["lca"] == "Bacteria;Proteobacteria"
        assert result["lca_depth"] == 2
        assert result["num_proteins"] == 2

    def test_assign_lca_unassigned(self):
        from metapeptide_lca_assigner import assign_lca_for_peptide
        proteins = {"P1": "XXXXXXX"}
        taxonomy = {}
        result = assign_lca_for_peptide("ZZZZZ", proteins, taxonomy)
        assert result["lca"] == "unassigned"

    def test_full_pipeline(self):
        from metapeptide_lca_assigner import assign_lca_batch, load_fasta, load_taxonomy, write_output

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir, {
                "P1": "ACDEFGHIKLMNPQR",
                "P2": "XXXXACDEFYYY",
            })
            tax_path = self._create_taxonomy(tmpdir, {
                "P1": "Bacteria;Proteobacteria;Gamma",
                "P2": "Bacteria;Proteobacteria;Beta",
            })
            proteins = load_fasta(fasta_path)
            taxonomy = load_taxonomy(tax_path)
            results = assign_lca_batch(["ACDEF", "GHIKLM"], proteins, taxonomy)
            output_path = os.path.join(tmpdir, "taxonomy.tsv")
            write_output(output_path, results)
            assert os.path.exists(output_path)
            assert len(results) == 2
