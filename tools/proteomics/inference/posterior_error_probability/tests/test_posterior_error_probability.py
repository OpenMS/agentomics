"""Tests for posterior_error_probability."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")
oms = pyopenms


def _create_bimodal_idxml(path):
    """Create idXML with bimodal score distribution (targets + decoys)."""
    protein_id = oms.ProteinIdentification()
    protein_id.setSearchEngine("test")
    protein_id.setScoreType("XCorr")
    protein_id.setHigherScoreBetter(True)
    protein_id.setIdentifier("run1")

    prot_hit = oms.ProteinHit()
    prot_hit.setAccession("protein1")
    prot_hit.setScore(100.0)

    prot_hit_decoy = oms.ProteinHit()
    prot_hit_decoy.setAccession("DECOY_protein1")
    prot_hit_decoy.setScore(0.0)

    protein_id.setHits([prot_hit, prot_hit_decoy])

    peptide_ids = oms.PeptideIdentificationList()

    # High-scoring targets (good matches)
    target_seqs = [
        "PEPTIDEK", "ACDEFTGHIK", "MNPQRSTWY",
        "DFGHIKLMN", "ACPEPTIDE", "SEQVENCELK",
        "ANOTHERPEPTIDE", "YETANOTHER", "FINALPEPTIDE", "LASTSEQK",
    ]
    for i, seq in enumerate(target_seqs):
        pep_id = oms.PeptideIdentification()
        pep_id.setScoreType("XCorr")
        pep_id.setHigherScoreBetter(True)
        pep_id.setRT(100.0 + i * 10)
        pep_id.setMZ(500.0 + i * 20)
        pep_id.setIdentifier("run1")

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq))
        hit.setScore(3.0 + i * 0.2)
        hit.setCharge(2)
        hit.setRank(1)
        hit.setMetaValue("target_decoy", "target")

        ev = oms.PeptideEvidence()
        ev.setProteinAccession("protein1")
        hit.setPeptideEvidences([ev])
        pep_id.setHits([hit])
        peptide_ids.push_back(pep_id)

    # Low-scoring decoys (random matches)
    decoy_seqs = [
        "EDITPEPK", "KIHGTFEDCA", "YWTSRQPNM",
        "NMLKIHGFD", "EDITPEPCA", "KLENECVEQS",
        "EDITPEPREHTONA", "REHTONATEYK", "EDITPEPLANIF", "KQESTSALK",
    ]
    for i, seq in enumerate(decoy_seqs):
        pep_id = oms.PeptideIdentification()
        pep_id.setScoreType("XCorr")
        pep_id.setHigherScoreBetter(True)
        pep_id.setRT(300.0 + i * 10)
        pep_id.setMZ(600.0 + i * 20)
        pep_id.setIdentifier("run1")

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq))
        hit.setScore(0.5 + i * 0.1)
        hit.setCharge(2)
        hit.setRank(1)
        hit.setMetaValue("target_decoy", "decoy")

        ev = oms.PeptideEvidence()
        ev.setProteinAccession("DECOY_protein1")
        hit.setPeptideEvidences([ev])
        pep_id.setHits([hit])
        peptide_ids.push_back(pep_id)

    oms.IdXMLFile().store(path, [protein_id], peptide_ids)


class TestPosteriorErrorProbability:

    def test_estimate_pep_produces_output(self):
        """Test that PEP estimation writes an output file."""
        from posterior_error_probability import estimate_pep

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "pep.idXML")
            _create_bimodal_idxml(input_path)

            n = estimate_pep(input_path, output_path)

            assert os.path.exists(output_path)
            assert n == 20

    def test_pep_values_between_0_and_1(self):
        """Test that PEP values are in [0, 1] range."""
        from posterior_error_probability import estimate_pep

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "pep.idXML")
            _create_bimodal_idxml(input_path)

            estimate_pep(input_path, output_path)

            protein_ids = []
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(output_path, protein_ids, peptide_ids)

            for pep_id in peptide_ids:
                for hit in pep_id.getHits():
                    score = hit.getScore()
                    assert 0.0 <= score <= 1.0, (
                        f"PEP value {score} out of [0, 1] range"
                    )

    def test_pep_returns_count(self):
        """Test that estimate_pep returns the number of scored PSMs."""
        from posterior_error_probability import estimate_pep

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "pep.idXML")
            _create_bimodal_idxml(input_path)

            n = estimate_pep(input_path, output_path)

            assert isinstance(n, int)
            assert n == 20
