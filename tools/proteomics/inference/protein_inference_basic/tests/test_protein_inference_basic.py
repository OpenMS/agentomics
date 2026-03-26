"""Tests for protein_inference_basic."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")
oms = pyopenms


def _create_peptide_protein_idxml(path):
    """Create idXML with 3 peptides mapping to 2 proteins (one shared peptide)."""
    protein_id = oms.ProteinIdentification()
    protein_id.setSearchEngine("test")
    protein_id.setScoreType("pep")
    protein_id.setHigherScoreBetter(False)
    protein_id.setIdentifier("run1")

    prot1 = oms.ProteinHit()
    prot1.setAccession("protein1")
    prot1.setScore(0.0)

    prot2 = oms.ProteinHit()
    prot2.setAccession("protein2")
    prot2.setScore(0.0)

    protein_id.setHits([prot1, prot2])

    peptide_ids = oms.PeptideIdentificationList()

    # Peptide 1: unique to protein1
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
    ev1 = oms.PeptideEvidence()
    ev1.setProteinAccession("protein1")
    hit1.setPeptideEvidences([ev1])
    pep_id1.setHits([hit1])
    peptide_ids.push_back(pep_id1)

    # Peptide 2: shared between protein1 and protein2
    pep_id2 = oms.PeptideIdentification()
    pep_id2.setScoreType("pep")
    pep_id2.setHigherScoreBetter(False)
    pep_id2.setRT(200.0)
    pep_id2.setMZ(600.0)
    pep_id2.setIdentifier("run1")

    hit2 = oms.PeptideHit()
    hit2.setSequence(oms.AASequence.fromString("SHAREDPEPTIDE"))
    hit2.setScore(0.02)
    hit2.setCharge(2)
    hit2.setRank(1)
    ev2a = oms.PeptideEvidence()
    ev2a.setProteinAccession("protein1")
    ev2b = oms.PeptideEvidence()
    ev2b.setProteinAccession("protein2")
    hit2.setPeptideEvidences([ev2a, ev2b])
    pep_id2.setHits([hit2])
    peptide_ids.push_back(pep_id2)

    # Peptide 3: unique to protein2
    pep_id3 = oms.PeptideIdentification()
    pep_id3.setScoreType("pep")
    pep_id3.setHigherScoreBetter(False)
    pep_id3.setRT(300.0)
    pep_id3.setMZ(700.0)
    pep_id3.setIdentifier("run1")

    hit3 = oms.PeptideHit()
    hit3.setSequence(oms.AASequence.fromString("ANOTHERPEPTIDE"))
    hit3.setScore(0.03)
    hit3.setCharge(2)
    hit3.setRank(1)
    ev3 = oms.PeptideEvidence()
    ev3.setProteinAccession("protein2")
    hit3.setPeptideEvidences([ev3])
    pep_id3.setHits([hit3])
    peptide_ids.push_back(pep_id3)

    oms.IdXMLFile().store(path, [protein_id], peptide_ids)


class TestProteinInferenceBasic:

    def test_infer_proteins_produces_output(self):
        """Test that protein inference writes an output file."""
        from protein_inference_basic import infer_proteins

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "proteins.idXML")
            _create_peptide_protein_idxml(input_path)

            n = infer_proteins(input_path, output_path)

            assert os.path.exists(output_path)
            assert n >= 2

    def test_both_proteins_inferred(self):
        """Test that both proteins are present in the output."""
        from protein_inference_basic import infer_proteins

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "proteins.idXML")
            _create_peptide_protein_idxml(input_path)

            infer_proteins(input_path, output_path)

            protein_ids = []
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(output_path, protein_ids, peptide_ids)

            accessions = set()
            for prot_id in protein_ids:
                for hit in prot_id.getHits():
                    accessions.add(hit.getAccession())

            assert "protein1" in accessions
            assert "protein2" in accessions

    def test_protein_scores_assigned(self):
        """Test that protein scores are assigned after inference."""
        from protein_inference_basic import infer_proteins

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "proteins.idXML")
            _create_peptide_protein_idxml(input_path)

            infer_proteins(input_path, output_path)

            protein_ids = []
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(output_path, protein_ids, peptide_ids)

            for prot_id in protein_ids:
                for hit in prot_id.getHits():
                    score = hit.getScore()
                    assert isinstance(score, float)
