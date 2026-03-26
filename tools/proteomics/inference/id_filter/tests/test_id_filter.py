"""Tests for id_filter."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")
oms = pyopenms


def _create_mixed_idxml(path):
    """Create idXML with 10 PSMs: varied scores + decoys."""
    protein_id = oms.ProteinIdentification()
    protein_id.setSearchEngine("test")
    protein_id.setScoreType("pep")
    protein_id.setHigherScoreBetter(False)
    protein_id.setIdentifier("run1")

    prot_target = oms.ProteinHit()
    prot_target.setAccession("protein1")
    prot_target.setScore(0.0)

    prot_decoy = oms.ProteinHit()
    prot_decoy.setAccession("DECOY_protein1")
    prot_decoy.setScore(1.0)

    protein_id.setHits([prot_target, prot_decoy])

    peptide_ids = oms.PeptideIdentificationList()

    # 5 target hits with varying scores: 0.01, 0.02, 0.04, 0.08, 0.2
    target_scores = [0.01, 0.02, 0.04, 0.08, 0.2]
    target_seqs = ["PEPTIDEK", "ACDEFTGHIK", "MNPQRSTWY", "DFGHIKLMN", "ACPEPTIDE"]
    for i, (seq, score) in enumerate(zip(target_seqs, target_scores)):
        pep_id = oms.PeptideIdentification()
        pep_id.setScoreType("pep")
        pep_id.setHigherScoreBetter(False)
        pep_id.setRT(100.0 + i * 10)
        pep_id.setMZ(500.0 + i * 50)
        pep_id.setIdentifier("run1")

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq))
        hit.setScore(score)
        hit.setCharge(2)
        hit.setRank(1)
        hit.setMetaValue("target_decoy", "target")

        ev = oms.PeptideEvidence()
        ev.setProteinAccession("protein1")
        hit.setPeptideEvidences([ev])
        pep_id.setHits([hit])
        peptide_ids.push_back(pep_id)

    # 5 decoy hits with varying scores: 0.3, 0.5, 0.6, 0.7, 0.9
    decoy_scores = [0.3, 0.5, 0.6, 0.7, 0.9]
    decoy_seqs = ["EDITPEPK", "KIHGTFEDCA", "YWTSRQPNM", "NMLKIHGFD", "EDITPEPCA"]
    for i, (seq, score) in enumerate(zip(decoy_seqs, decoy_scores)):
        pep_id = oms.PeptideIdentification()
        pep_id.setScoreType("pep")
        pep_id.setHigherScoreBetter(False)
        pep_id.setRT(200.0 + i * 10)
        pep_id.setMZ(600.0 + i * 50)
        pep_id.setIdentifier("run1")

        hit = oms.PeptideHit()
        hit.setSequence(oms.AASequence.fromString(seq))
        hit.setScore(score)
        hit.setCharge(2)
        hit.setRank(1)
        hit.setMetaValue("target_decoy", "decoy")

        ev = oms.PeptideEvidence()
        ev.setProteinAccession("DECOY_protein1")
        hit.setPeptideEvidences([ev])
        pep_id.setHits([hit])
        peptide_ids.push_back(pep_id)

    oms.IdXMLFile().store(path, [protein_id], peptide_ids)


class TestIDFilter:

    def test_filter_by_score(self):
        """Test score-based filtering retains correct number of hits."""
        from id_filter import filter_ids

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "filtered.idXML")
            _create_mixed_idxml(input_path)

            # Filter at 0.05 threshold (lower is better)
            # Targets with score <= 0.05: 0.01, 0.02, 0.04 => 3 hits
            result = filter_ids(input_path, output_path, score_threshold=0.05)

            assert result["before"] == 10
            assert result["after"] == 3
            assert result["removed"] == 7
            assert os.path.exists(output_path)

    def test_remove_decoys(self):
        """Test decoy removal keeps only target hits."""
        from id_filter import filter_ids

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "filtered.idXML")
            _create_mixed_idxml(input_path)

            # Only remove decoys, no score filtering
            result = filter_ids(
                input_path, output_path,
                score_threshold=None, remove_decoys=True,
            )

            assert result["before"] == 10
            assert result["after"] == 5  # only 5 targets remain

    def test_filter_and_remove_decoys(self):
        """Test combined score filtering and decoy removal."""
        from id_filter import filter_ids

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "filtered.idXML")
            _create_mixed_idxml(input_path)

            result = filter_ids(
                input_path, output_path,
                score_threshold=0.05, remove_decoys=True,
            )

            assert result["before"] == 10
            # Should keep only targets with score <= 0.05
            assert result["after"] == 3

    def test_returns_dict(self):
        """Test that filter_ids returns a dict with expected keys."""
        from id_filter import filter_ids

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.idXML")
            output_path = os.path.join(tmp, "filtered.idXML")
            _create_mixed_idxml(input_path)

            result = filter_ids(input_path, output_path)

            assert "before" in result
            assert "after" in result
            assert "removed" in result
