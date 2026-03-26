"""Tests for consensus_map_normalizer."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestCreateSyntheticConsensusMap:
    def test_creates_valid_file(self):
        import pyopenms as oms
        from consensus_map_normalizer import create_synthetic_consensus_map

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "consensus.consensusXML")
            create_synthetic_consensus_map(path, n_features=3, n_channels=2)

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(path, cm)
            assert cm.size() == 3

    def test_has_correct_channels(self):
        import pyopenms as oms
        from consensus_map_normalizer import create_synthetic_consensus_map

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "consensus.consensusXML")
            create_synthetic_consensus_map(path, n_features=2, n_channels=3)

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(path, cm)
            headers = cm.getColumnHeaders()
            assert len(headers) == 3

    def test_intensity_bias_applied(self):
        import pyopenms as oms
        from consensus_map_normalizer import create_synthetic_consensus_map

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "consensus.consensusXML")
            create_synthetic_consensus_map(
                path, n_features=3, n_channels=2, intensity_bias=2.0
            )

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(path, cm)

            for cf in cm:
                handles = cf.getFeatureList()
                intensities = {h.getMapIndex(): h.getIntensity() for h in handles}
                # Channel 1 should have 2x the intensity of channel 0
                ratio = intensities[1] / intensities[0]
                assert abs(ratio - 2.0) < 0.01


class TestNormalizeConsensus:
    def test_median_normalization(self):
        import pyopenms as oms
        from consensus_map_normalizer import (
            create_synthetic_consensus_map,
            normalize_consensus,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.consensusXML")
            out_path = os.path.join(tmp, "output.consensusXML")

            create_synthetic_consensus_map(
                in_path, n_features=5, n_channels=2, intensity_bias=2.0
            )

            n = normalize_consensus(in_path, out_path, method="median")
            assert n == 5
            assert os.path.exists(out_path)

            # Load normalized map and check intensities are closer
            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cm)
            for cf in cm:
                handles = cf.getFeatureList()
                intensities = [h.getIntensity() for h in handles]
                if all(i > 0 for i in intensities):
                    ratio = max(intensities) / min(intensities)
                    # After median normalization, ratio should be much closer to 1
                    assert ratio < 1.5, f"Ratio {ratio} too high after normalization"

    def test_quantile_normalization(self):
        from consensus_map_normalizer import (
            create_synthetic_consensus_map,
            normalize_consensus,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.consensusXML")
            out_path = os.path.join(tmp, "output.consensusXML")

            create_synthetic_consensus_map(
                in_path, n_features=5, n_channels=2, intensity_bias=3.0
            )

            n = normalize_consensus(in_path, out_path, method="quantile")
            assert n == 5
            assert os.path.exists(out_path)

    def test_invalid_method(self):
        with pytest.raises(ValueError, match="Unknown method"):
            from consensus_map_normalizer import normalize_consensus

            normalize_consensus("dummy.xml", "out.xml", method="invalid")

    def test_round_trip_file(self):
        import pyopenms as oms
        from consensus_map_normalizer import (
            create_synthetic_consensus_map,
            normalize_consensus,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.consensusXML")
            out_path = os.path.join(tmp, "output.consensusXML")

            create_synthetic_consensus_map(
                in_path, n_features=4, n_channels=2, intensity_bias=2.0
            )

            normalize_consensus(in_path, out_path, method="median")

            cm = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cm)
            assert cm.size() == 4
            headers = cm.getColumnHeaders()
            assert len(headers) == 2
