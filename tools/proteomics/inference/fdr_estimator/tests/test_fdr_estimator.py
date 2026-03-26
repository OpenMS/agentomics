"""Tests for fdr_estimator."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")
oms = pyopenms


def _create_target_decoy_idxml(path):
    """Create idXML with 5 target (high score) + 5 decoy (low score) PSMs."""
    protein_id = oms.ProteinIdentification()
    protein_id.setSearchEngine("test")
    protein_id.setScoreType("q-value")
    protein_id.setHigherScoreBetter(False)
    protein_id.setIdentifier("run1")

    prot_hit_target = oms.ProteinHit()
    prot_hit_target.setAccession("protein1")
    prot_hit_target.setScore(0.0)
    prot_hit_target.setMetaValue("target_decoy", "target")

    prot_hit_decoy = oms.ProteinHit()
    prot_hit_decoy.setAccession("DECOY_protein1")
    prot_hit_decoy.setScore(1.0)
    prot_hit_decoy.setMetaValue("target_decoy", "decoy")

    protein_id.setHits([prot_hit_target, prot_hit_decoy])

    peptide_ids = oms.PeptideIdentificationList()

    # 5 target PSMs with good scores
    target_sequences = ["PEPTIDEK", "ACDEFTGHIK", "MNPQRSTWY", "DFGHIKLMN", "ACPEPTIDE"]
    for i, seq in enumerate(target_sequences):
        pep_id = oms.PeptideIdentification()
        pep_id.setScoreType("pep")
        pep_id.setHigherScoreBetter(False)
        pep_id.setRT(100.0 + i * 10)
        pep_id.setMZ(500.0 + i * 50)
        pep_id.setIdentifier("run1")

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq))
        hit.setScore(0.001 + i * 0.001)
        hit.setCharge(2)
        hit.setRank(1)
        hit.setMetaValue("target_decoy", "target")

        ev = oms.PeptideEvidence()
        ev.setProteinAccession("protein1")
        hit.setPeptideEvidences([ev])
        pep_id.setHits([hit])
        peptide_ids.push_back(pep_id)

    # 5 decoy PSMs with bad scores
    decoy_sequences = ["EDITPEPK", "KIHGTFEDCA", "YWTSRQPNM", "NMLKIHGFD", "EDITPEPCA"]
    for i, seq in enumerate(decoy_sequences):
        pep_id = oms.PeptideIdentification()
        pep_id.setScoreType("pep")
        pep_id.setHigherScoreBetter(False)
        pep_id.setRT(200.0 + i * 10)
        pep_id.setMZ(600.0 + i * 50)
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


class TestFDREstimator:

    def test_estimate_fdr_produces_output(self):
        """Test that FDR estimation writes an output file."""
        from fdr_estimator import estimate_fdr

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "fdr.idXML")
            _create_target_decoy_idxml(input_path)

            result = estimate_fdr(input_path, output_path)

            assert os.path.exists(output_path)
            # FDR apply computes q-values; targets should survive
            assert result["total_psms"] >= 5

    def test_qvalues_assigned(self):
        """Test that q-values are assigned to PSMs after FDR estimation."""
        from fdr_estimator import estimate_fdr

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "fdr.idXML")
            _create_target_decoy_idxml(input_path)

            estimate_fdr(input_path, output_path)

            protein_ids = []
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(output_path, protein_ids, peptide_ids)

            # At least target PSMs should remain with q-values
            assert len(peptide_ids) >= 5
            for pep_id in peptide_ids:
                for hit in pep_id.getHits():
                    score = hit.getScore()
                    assert 0.0 <= score <= 1.0

    def test_counts_at_thresholds(self):
        """Test that count dict is returned with expected FDR thresholds."""
        from fdr_estimator import estimate_fdr

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "fdr.idXML")
            _create_target_decoy_idxml(input_path)

            result = estimate_fdr(input_path, output_path)

            assert "counts_at_fdr" in result
            counts = result["counts_at_fdr"]
            assert 0.01 in counts
            assert 0.05 in counts
            assert 0.1 in counts
            # At looser threshold we should have at least as many PSMs
            assert counts[0.1] >= counts[0.01]

    def test_protein_level_fdr(self):
        """Test protein-level FDR estimation."""
        from fdr_estimator import estimate_fdr

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "fdr.idXML")
            _create_target_decoy_idxml(input_path)

            result = estimate_fdr(input_path, output_path, protein_level=True)

            assert os.path.exists(output_path)
            assert result["total_psms"] >= 5
