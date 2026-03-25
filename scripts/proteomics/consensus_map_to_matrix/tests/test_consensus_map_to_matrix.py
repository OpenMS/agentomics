"""Tests for consensus_map_to_matrix."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
def test_create_synthetic_consensus():
    import pyopenms as oms
    from consensus_map_to_matrix import create_synthetic_consensus

    with tempfile.TemporaryDirectory() as tmp:
        cxml_path = os.path.join(tmp, "consensus.consensusXML")
        create_synthetic_consensus(cxml_path, n_features=5, n_maps=3)

        cmap = oms.ConsensusMap()
        oms.ConsensusXMLFile().load(cxml_path, cmap)
        assert cmap.size() == 5


@requires_pyopenms
def test_consensus_to_matrix():
    from consensus_map_to_matrix import consensus_to_matrix, create_synthetic_consensus

    with tempfile.TemporaryDirectory() as tmp:
        cxml_path = os.path.join(tmp, "consensus.consensusXML")
        tsv_path = os.path.join(tmp, "matrix.tsv")

        create_synthetic_consensus(cxml_path, n_features=5, n_maps=3)
        stats = consensus_to_matrix(cxml_path, tsv_path)

        assert stats["consensus_features"] == 5
        assert stats["n_maps"] == 3

        with open(tsv_path) as fh:
            lines = fh.readlines()
        assert len(lines) == 6  # header + 5 rows
        header = lines[0].strip().split("\t")
        assert "rt" in header
        assert "mz" in header
        assert "intensity_0" in header
        assert "intensity_2" in header


@requires_pyopenms
def test_consensus_to_matrix_single_map():
    from consensus_map_to_matrix import consensus_to_matrix, create_synthetic_consensus

    with tempfile.TemporaryDirectory() as tmp:
        cxml_path = os.path.join(tmp, "consensus.consensusXML")
        tsv_path = os.path.join(tmp, "matrix.tsv")

        create_synthetic_consensus(cxml_path, n_features=3, n_maps=1)
        stats = consensus_to_matrix(cxml_path, tsv_path)

        assert stats["n_maps"] == 1
        with open(tsv_path) as fh:
            header = fh.readline().strip().split("\t")
        assert "intensity_0" in header
