"""Tests for spectral_counting_quantifier."""

import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestSpectralCountingQuantifier:
    def _create_fasta(self, tmpdir):
        import pyopenms as oms

        fasta_path = f"{tmpdir}/test.fasta"
        entries = []
        e1 = oms.FASTAEntry()
        e1.identifier = "PROT1"
        e1.sequence = "MSPEPTIDEKAAANOTHERPEPTIDER"
        entries.append(e1)
        e2 = oms.FASTAEntry()
        e2.identifier = "PROT2"
        e2.sequence = "MSGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGK"
        entries.append(e2)
        oms.FASTAFile().store(fasta_path, entries)
        return fasta_path

    def test_count_observable_peptides(self):
        from spectral_counting_quantifier import count_observable_peptides

        count = count_observable_peptides("MSPEPTIDEKAAANOTHERPEPTIDER")
        assert count >= 1

    def test_calculate_nsaf(self):
        from spectral_counting_quantifier import calculate_nsaf

        protein_data = {
            "PROT1": {"spectral_count": 10},
            "PROT2": {"spectral_count": 5},
        }
        proteins = {"PROT1": "AAAAAAAAAA", "PROT2": "AAAAA"}
        results = calculate_nsaf(protein_data, proteins)
        assert len(results) == 2
        # NSAF values should sum to ~1
        total_nsaf = sum(r["nsaf"] for r in results)
        assert abs(total_nsaf - 1.0) < 0.001

    def test_calculate_empai(self):
        from spectral_counting_quantifier import calculate_empai

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self._create_fasta(tmpdir)
            from spectral_counting_quantifier import load_fasta_proteins
            proteins = load_fasta_proteins(fasta_path)
            protein_data = {
                "PROT1": {"spectral_count": 10, "observed_peptides": 2},
            }
            results = calculate_empai(protein_data, proteins)
            assert len(results) == 1
            assert results[0]["empai"] > 0

    def test_nsaf_proportional(self):
        from spectral_counting_quantifier import calculate_nsaf

        # Same length proteins, different counts -> NSAF proportional to SpC
        proteins = {"P1": "AAAAAAAAAA", "P2": "AAAAAAAAAA"}
        data = {"P1": {"spectral_count": 20}, "P2": {"spectral_count": 10}}
        results = calculate_nsaf(data, proteins)
        nsaf_map = {r["accession"]: r["nsaf"] for r in results}
        assert nsaf_map["P1"] > nsaf_map["P2"]

    def test_load_peptide_counts(self):
        from spectral_counting_quantifier import load_peptide_counts

        with tempfile.TemporaryDirectory() as tmpdir:
            tsv_path = f"{tmpdir}/counts.tsv"
            with open(tsv_path, "w") as fh:
                fh.write("protein\tpeptide\tspectral_count\n")
                fh.write("PROT1\tPEPTIDEK\t5\n")
                fh.write("PROT1\tANOTHER\t3\n")
                fh.write("PROT2\tSEQUENCE\t2\n")
            data = load_peptide_counts(tsv_path)
            assert data["PROT1"]["spectral_count"] == 8
            assert data["PROT1"]["observed_peptides"] == 2
            assert data["PROT2"]["spectral_count"] == 2
