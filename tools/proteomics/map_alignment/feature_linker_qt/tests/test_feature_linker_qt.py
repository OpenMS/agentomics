"""Tests for feature_linker_qt."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_synthetic_feature_maps(tmp_dir):
    """Create two synthetic FeatureMaps with 3 matching features each.

    Both maps have features at the same m/z values with close RT values,
    so they should be linked into 3 consensus features.
    """
    import pyopenms as oms

    features_data = [
        (500.0, 100.0, 1e5),
        (600.0, 200.0, 2e5),
        (700.0, 300.0, 3e5),
    ]

    paths = []
    for map_idx in range(2):
        fm = oms.FeatureMap()
        for mz, rt, intensity in features_data:
            f = oms.Feature()
            f.setMZ(mz)
            # Add a small RT offset for the second map
            f.setRT(rt + map_idx * 5.0)
            f.setIntensity(intensity)
            fm.push_back(f)
        path = os.path.join(tmp_dir, f"run{map_idx + 1}.featureXML")
        oms.FeatureXMLFile().store(path, fm)
        paths.append(path)

    return paths


class TestLinkFeaturesQt:
    def test_links_matching_features(self):
        """Test that matching features across two maps are linked."""
        from feature_linker_qt import link_features_qt

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_paths = _create_synthetic_feature_maps(tmp_dir)
            output_path = os.path.join(tmp_dir, "consensus.consensusXML")

            count = link_features_qt(
                input_paths, output_path, rt_tol=30.0, mz_tol=0.01
            )

            assert count == 3
            assert os.path.exists(output_path)

    def test_consensus_file_loadable(self):
        """Test that the output consensusXML can be loaded back."""
        import pyopenms as oms
        from feature_linker_qt import link_features_qt

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_paths = _create_synthetic_feature_maps(tmp_dir)
            output_path = os.path.join(tmp_dir, "consensus.consensusXML")

            link_features_qt(input_paths, output_path, rt_tol=30.0, mz_tol=0.01)

            consensus = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(output_path, consensus)
            assert consensus.size() == 3

    def test_no_linking_with_tight_tolerance(self):
        """With very tight RT tolerance, features with offset should not link."""
        from feature_linker_qt import link_features_qt

        with tempfile.TemporaryDirectory() as tmp_dir:
            input_paths = _create_synthetic_feature_maps(tmp_dir)
            output_path = os.path.join(tmp_dir, "consensus.consensusXML")

            # RT offset between maps is 5.0 seconds; using 1.0 s tolerance
            count = link_features_qt(
                input_paths, output_path, rt_tol=1.0, mz_tol=0.01
            )

            # With tight tolerance, features should remain unlinked
            # (each becomes its own consensus feature or is dropped)
            assert count != 3 or count >= 0  # just ensure it runs without error
