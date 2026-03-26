"""Tests for peptide_indexer."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")
oms = pyopenms


def _create_fasta(path):
    """Create a simple FASTA file with a known protein containing PEPTIDEK."""
    entries = []

    entry1 = oms.FASTAEntry()
    entry1.identifier = "protein1"
    entry1.description = "Test protein 1"
    entry1.sequence = "AAAPEPTIDEKLLL"
    entries.append(entry1)

    entry2 = oms.FASTAEntry()
    entry2.identifier = "protein2"
    entry2.description = "Test protein 2"
    entry2.sequence = "CCCMNPQRSTWYDDD"
    entries.append(entry2)

    oms.FASTAFile().store(path, entries)


def _create_peptide_idxml(path):
    """Create idXML with peptides that should map to proteins in the FASTA."""
    protein_id = oms.ProteinIdentification()
    protein_id.setSearchEngine("test")
    protein_id.setScoreType("pep")
    protein_id.setHigherScoreBetter(False)
    protein_id.setIdentifier("run1")
    protein_id.setHits([])

    peptide_ids = oms.PeptideIdentificationList()

    # Peptide that should map to protein1
    pep_id1 = oms.PeptideIdentification()
    pep_id1.setScoreType("pep")
    pep_id1.setHigherScoreBetter(False)
    pep_id1.setRT(100.0)
    pep_id1.setMZ(500.0)
    pep_id1.setIdentifier("run1")

    hit1 = oms.PeptideHit()
    hit1.setSequence(oms.AASequence.fromString("PEPTIDEK"))
    hit1.setScore(0.01)
    hit1.setCharge(2)
    hit1.setRank(1)
    pep_id1.setHits([hit1])
    peptide_ids.push_back(pep_id1)

    # Peptide that should map to protein2
    pep_id2 = oms.PeptideIdentification()
    pep_id2.setScoreType("pep")
    pep_id2.setHigherScoreBetter(False)
    pep_id2.setRT(200.0)
    pep_id2.setMZ(600.0)
    pep_id2.setIdentifier("run1")

    hit2 = oms.PeptideHit()
    hit2.setSequence(oms.AASequence.fromString("MNPQRSTWY"))
    hit2.setScore(0.02)
    hit2.setCharge(2)
    hit2.setRank(1)
    pep_id2.setHits([hit2])
    peptide_ids.push_back(pep_id2)

    oms.IdXMLFile().store(path, [protein_id], peptide_ids)


class TestPeptideIndexer:

    def test_index_peptides_produces_output(self):
        """Test that peptide indexing writes an output file."""
        from peptide_indexer import index_peptides

        with tempfile.TemporaryDirectory() as tmp:
            fasta_path = os.path.join(tmp, "database.fasta")
            ids_path = os.path.join(tmp, "peptides.idXML")
            output_path = os.path.join(tmp, "indexed.idXML")

            _create_fasta(fasta_path)
            _create_peptide_idxml(ids_path)

            n = index_peptides(ids_path, fasta_path, output_path)

            assert os.path.exists(output_path)
            assert n == 2

    def test_protein_reference_added(self):
        """Test that protein references are added to peptide hits."""
        from peptide_indexer import index_peptides

        with tempfile.TemporaryDirectory() as tmp:
            fasta_path = os.path.join(tmp, "database.fasta")
            ids_path = os.path.join(tmp, "peptides.idXML")
            output_path = os.path.join(tmp, "indexed.idXML")

            _create_fasta(fasta_path)
            _create_peptide_idxml(ids_path)

            index_peptides(ids_path, fasta_path, output_path)

            protein_ids = []
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(output_path, protein_ids, peptide_ids)

            # Check that at least one peptide has protein evidence
            has_evidence = False
            for pep_id in peptide_ids:
                for hit in pep_id.getHits():
                    evidences = hit.getPeptideEvidences()
                    if len(evidences) > 0:
                        has_evidence = True
                        accessions = [ev.getProteinAccession() for ev in evidences]
                        seq = hit.getSequence().toString()
                        if seq == "PEPTIDEK":
                            assert "protein1" in accessions

            assert has_evidence, "No protein evidence found after indexing"

    def test_returns_hit_count(self):
        """Test that index_peptides returns the correct hit count."""
        from peptide_indexer import index_peptides

        with tempfile.TemporaryDirectory() as tmp:
            fasta_path = os.path.join(tmp, "database.fasta")
            ids_path = os.path.join(tmp, "peptides.idXML")
            output_path = os.path.join(tmp, "indexed.idXML")

            _create_fasta(fasta_path)
            _create_peptide_idxml(ids_path)

            n = index_peptides(ids_path, fasta_path, output_path)

            assert isinstance(n, int)
            assert n == 2
