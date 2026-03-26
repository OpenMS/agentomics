"""Tests for feature_linker_labeled."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")

# SILAC K+8 mass shift
SILAC_K8_DELTA = 8.0142


def _create_silac_feature_map(tmp_dir):
    """Create a synthetic FeatureMap with light/heavy SILAC pairs.

    Each pair has a delta m/z of 8.0142 Da (K+8 SILAC label) at charge 1.
    Features are assigned unique IDs so the algorithm can track them.
    """
    import pyopenms as oms

    fm = oms.FeatureMap()
    fm.ensureUniqueId()
    # Light/heavy pairs with K+8 delta (8.0142 Da)
    pairs = [
        (500.0, 100.0, 1e5, 500.0 + SILAC_K8_DELTA, 100.0, 8e4),
        (600.0, 200.0, 2e5, 600.0 + SILAC_K8_DELTA, 200.0, 1.5e5),
        (700.0, 300.0, 3e5, 700.0 + SILAC_K8_DELTA, 300.0, 2.5e5),
    ]
    for light_mz, light_rt, light_int, heavy_mz, heavy_rt, heavy_int in pairs:
        f_light = oms.Feature()
        f_light.setMZ(light_mz)
        f_light.setRT(light_rt)
        f_light.setIntensity(light_int)
        f_light.setCharge(1)
        f_light.ensureUniqueId()
        fm.push_back(f_light)

        f_heavy = oms.Feature()
        f_heavy.setMZ(heavy_mz)
        f_heavy.setRT(heavy_rt)
        f_heavy.setIntensity(heavy_int)
        f_heavy.setCharge(1)
        f_heavy.ensureUniqueId()
        fm.push_back(f_heavy)

    path = os.path.join(tmp_dir, "silac_run.featureXML")
    oms.FeatureXMLFile().store(path, fm)
    return path


class TestLinkLabeledFeatures:
    def test_links_silac_pairs(self):
        """Test that SILAC light/heavy pairs are linked."""
        from feature_linker_labeled import link_labeled_features

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = _create_silac_feature_map(tmp_dir)
            output_path = os.path.join(tmp_dir, "consensus.consensusXML")

            count = link_labeled_features(
                input_path,
                output_path,
                rt_estimate=1000.0,
                mz_estimate=5.0,
                mz_pair_dists=[SILAC_K8_DELTA],
                rt_pair_dist=0.0,
            )

            # Should find 3 consensus features (one per light/heavy pair)
            assert count == 3
            assert os.path.exists(output_path)

    def test_consensus_file_loadable(self):
        """Test that the output consensusXML can be loaded back."""
        import pyopenms as oms
        from feature_linker_labeled import link_labeled_features

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = _create_silac_feature_map(tmp_dir)
            output_path = os.path.join(tmp_dir, "consensus.consensusXML")

            link_labeled_features(
                input_path,
                output_path,
                rt_estimate=1000.0,
                mz_estimate=5.0,
                mz_pair_dists=[SILAC_K8_DELTA],
                rt_pair_dist=0.0,
            )

            consensus = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(output_path, consensus)
            assert consensus.size() == 3

    def test_returns_int(self):
        """Test that the function returns an integer count."""
        from feature_linker_labeled import link_labeled_features

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = _create_silac_feature_map(tmp_dir)
            output_path = os.path.join(tmp_dir, "consensus.consensusXML")

            count = link_labeled_features(
                input_path,
                output_path,
                rt_estimate=1000.0,
                mz_estimate=5.0,
                mz_pair_dists=[SILAC_K8_DELTA],
                rt_pair_dist=0.0,
            )

            assert isinstance(count, int)
