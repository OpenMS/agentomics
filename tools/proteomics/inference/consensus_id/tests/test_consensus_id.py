"""Tests for consensus_id."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")
oms = pyopenms


def _create_search_idxml(path, engine_name, sequences_scores):
    """Create an idXML from a specific search engine with given peptides.

    Parameters
    ----------
    path : str
        Output path.
    engine_name : str
        Search engine name.
    sequences_scores : list of (str, float)
        List of (sequence, score) pairs.
    """
    protein_id = oms.ProteinIdentification()
    protein_id.setSearchEngine(engine_name)
    protein_id.setScoreType("pep")
    protein_id.setHigherScoreBetter(False)
    protein_id.setIdentifier("run1")

    prot_hit = oms.ProteinHit()
    prot_hit.setAccession("protein1")
    prot_hit.setScore(0.0)
    protein_id.setHits([prot_hit])

    peptide_ids = oms.PeptideIdentificationList()

    for i, (seq, score) in enumerate(sequences_scores):
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

        ev = oms.PeptideEvidence()
        ev.setProteinAccession("protein1")
        hit.setPeptideEvidences([ev])
        pep_id.setHits([hit])
        peptide_ids.push_back(pep_id)

    oms.IdXMLFile().store(path, [protein_id], peptide_ids)


class TestConsensusID:

    def test_consensus_best_produces_output(self):
        """Test that consensus ID with 'best' algorithm writes output."""
        from consensus_id import consensus_id

        with tempfile.TemporaryDirectory() as tmp:
            path1 = os.path.join(tmp, "search1.idXML")
            path2 = os.path.join(tmp, "search2.idXML")
            output = os.path.join(tmp, "consensus.idXML")

            _create_search_idxml(path1, "EngineA", [
                ("PEPTIDEK", 0.01),
                ("ACDEFTGHIK", 0.05),
            ])
            _create_search_idxml(path2, "EngineB", [
                ("PEPTIDEK", 0.02),
                ("MNPQRSTWY", 0.03),
            ])

            n = consensus_id([path1, path2], output, algorithm="best")

            assert os.path.exists(output)
            assert n > 0

    def test_consensus_merges_results(self):
        """Test that consensus merges overlapping peptide hits."""
        from consensus_id import consensus_id

        with tempfile.TemporaryDirectory() as tmp:
            path1 = os.path.join(tmp, "search1.idXML")
            path2 = os.path.join(tmp, "search2.idXML")
            output = os.path.join(tmp, "consensus.idXML")

            _create_search_idxml(path1, "EngineA", [
                ("PEPTIDEK", 0.01),
                ("ACDEFTGHIK", 0.05),
            ])
            _create_search_idxml(path2, "EngineB", [
                ("PEPTIDEK", 0.02),
                ("MNPQRSTWY", 0.03),
            ])

            consensus_id([path1, path2], output, algorithm="best")

            protein_ids = []
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(output, protein_ids, peptide_ids)

            assert len(peptide_ids) > 0

    def test_consensus_average_algorithm(self):
        """Test that 'average' algorithm works."""
        from consensus_id import consensus_id

        with tempfile.TemporaryDirectory() as tmp:
            path1 = os.path.join(tmp, "search1.idXML")
            path2 = os.path.join(tmp, "search2.idXML")
            output = os.path.join(tmp, "consensus.idXML")

            _create_search_idxml(path1, "EngineA", [("PEPTIDEK", 0.01)])
            _create_search_idxml(path2, "EngineB", [("PEPTIDEK", 0.02)])

            n = consensus_id([path1, path2], output, algorithm="average")

            assert os.path.exists(output)
            assert n > 0

    def test_consensus_ranks_algorithm(self):
        """Test that 'ranks' algorithm works."""
        from consensus_id import consensus_id

        with tempfile.TemporaryDirectory() as tmp:
            path1 = os.path.join(tmp, "search1.idXML")
            path2 = os.path.join(tmp, "search2.idXML")
            output = os.path.join(tmp, "consensus.idXML")

            _create_search_idxml(path1, "EngineA", [("PEPTIDEK", 0.01)])
            _create_search_idxml(path2, "EngineB", [("PEPTIDEK", 0.02)])

            n = consensus_id([path1, path2], output, algorithm="ranks")

            assert os.path.exists(output)
            assert n > 0

    def test_invalid_algorithm_raises(self):
        """Test that an invalid algorithm name raises ValueError."""
        from consensus_id import consensus_id

        with pytest.raises(ValueError, match="Unknown algorithm"):
            consensus_id([], "output.idXML", algorithm="invalid")
