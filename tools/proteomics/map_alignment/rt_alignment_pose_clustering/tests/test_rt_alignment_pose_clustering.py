"""Tests for rt_alignment_pose_clustering."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestAlignPoseClustering:
    def test_align_removes_rt_offset(self):
        """Two FeatureMaps with +10.0 RT offset should align to < 1.0 difference."""
        import pyopenms as oms
        from rt_alignment_pose_clustering import (
            align_pose_clustering,
            create_synthetic_featurexml,
        )

        with tempfile.TemporaryDirectory() as tmp:
            ref_path = os.path.join(tmp, "reference.featureXML")
            inp_path = os.path.join(tmp, "input.featureXML")
            out_path = os.path.join(tmp, "aligned.featureXML")

            create_synthetic_featurexml(ref_path, n_features=10, rt_offset=0.0)
            create_synthetic_featurexml(inp_path, n_features=10, rt_offset=10.0)

            stats = align_pose_clustering(ref_path, inp_path, out_path)

            assert os.path.exists(out_path)
            assert stats["n_features"] == 10
            assert stats["max_rt_shift"] >= 0.0

            # Load aligned map and reference, compare RTs
            ref_fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(ref_path, ref_fm)
            aligned_fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(out_path, aligned_fm)

            for i in range(ref_fm.size()):
                rt_diff = abs(ref_fm[i].getRT() - aligned_fm[i].getRT())
                assert rt_diff < 1.0, (
                    f"Feature {i}: RT diff {rt_diff:.2f} >= 1.0"
                )

    def test_align_returns_stats(self):
        """Stats dict should contain expected keys."""
        from rt_alignment_pose_clustering import (
            align_pose_clustering,
            create_synthetic_featurexml,
        )

        with tempfile.TemporaryDirectory() as tmp:
            ref_path = os.path.join(tmp, "reference.featureXML")
            inp_path = os.path.join(tmp, "input.featureXML")
            out_path = os.path.join(tmp, "aligned.featureXML")

            create_synthetic_featurexml(ref_path, n_features=5, rt_offset=0.0)
            create_synthetic_featurexml(inp_path, n_features=5, rt_offset=5.0)

            stats = align_pose_clustering(ref_path, inp_path, out_path)

            assert "n_features" in stats
            assert "max_rt_shift" in stats
            assert "n_trafo_points" in stats
            assert isinstance(stats["n_features"], int)
            assert isinstance(stats["max_rt_shift"], float)

    def test_align_preserves_feature_count(self):
        """Alignment should not change the number of features."""
        import pyopenms as oms
        from rt_alignment_pose_clustering import (
            align_pose_clustering,
            create_synthetic_featurexml,
        )

        with tempfile.TemporaryDirectory() as tmp:
            ref_path = os.path.join(tmp, "reference.featureXML")
            inp_path = os.path.join(tmp, "input.featureXML")
            out_path = os.path.join(tmp, "aligned.featureXML")

            create_synthetic_featurexml(ref_path, n_features=8, rt_offset=0.0)
            create_synthetic_featurexml(inp_path, n_features=8, rt_offset=20.0)

            align_pose_clustering(ref_path, inp_path, out_path)

            inp_fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(inp_path, inp_fm)
            aligned_fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(out_path, aligned_fm)

            assert aligned_fm.size() == inp_fm.size()

    def test_create_synthetic_featurexml(self):
        """Synthetic FeatureXML should have correct number of features."""
        import pyopenms as oms
        from rt_alignment_pose_clustering import create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "test.featureXML")
            create_synthetic_featurexml(path, n_features=7, rt_offset=50.0)

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(path, fm)
            assert fm.size() == 7
            assert fm[0].getRT() == pytest.approx(150.0, abs=0.01)
