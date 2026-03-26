"""Tests for rt_alignment_kd."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestAlignKD:
    def test_align_two_maps_reduces_rt_offset(self):
        """Two FeatureMaps with systematic RT offset should have reduced difference."""
        import pyopenms as oms
        from rt_alignment_kd import align_kd, create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            p1 = os.path.join(tmp, "run1.featureXML")
            p2 = os.path.join(tmp, "run2.featureXML")
            out_dir = os.path.join(tmp, "aligned")

            create_synthetic_featurexml(p1, n_features=30, rt_offset=0.0)
            create_synthetic_featurexml(p2, n_features=30, rt_offset=15.0)

            stats = align_kd([p1, p2], out_dir)

            assert stats["n_maps"] == 2
            assert stats["n_features_per_map"] == [30, 30]

            # Check output files exist
            out1 = os.path.join(out_dir, "run1.featureXML")
            out2 = os.path.join(out_dir, "run2.featureXML")
            assert os.path.exists(out1)
            assert os.path.exists(out2)

            # Load aligned maps and check that RT offset is reduced
            fm1 = oms.FeatureMap()
            fm2 = oms.FeatureMap()
            oms.FeatureXMLFile().load(out1, fm1)
            oms.FeatureXMLFile().load(out2, fm2)

            # After alignment, average RT difference should be smaller than 15.0
            total_diff = 0.0
            for i in range(fm1.size()):
                total_diff += abs(fm1[i].getRT() - fm2[i].getRT())
            avg_diff = total_diff / fm1.size()
            assert avg_diff < 15.0, f"Average RT diff {avg_diff:.2f} not reduced"

    def test_align_returns_stats(self):
        """Stats dict should contain expected keys."""
        from rt_alignment_kd import align_kd, create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            p1 = os.path.join(tmp, "run1.featureXML")
            p2 = os.path.join(tmp, "run2.featureXML")
            out_dir = os.path.join(tmp, "aligned")

            create_synthetic_featurexml(p1, n_features=30, rt_offset=0.0)
            create_synthetic_featurexml(p2, n_features=30, rt_offset=10.0)

            stats = align_kd([p1, p2], out_dir)

            assert "n_maps" in stats
            assert "n_features_per_map" in stats
            assert "output_dir" in stats
            assert stats["n_maps"] == 2

    def test_align_preserves_feature_count(self):
        """Alignment should not change the number of features in each map."""
        import pyopenms as oms
        from rt_alignment_kd import align_kd, create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            p1 = os.path.join(tmp, "run1.featureXML")
            p2 = os.path.join(tmp, "run2.featureXML")
            out_dir = os.path.join(tmp, "aligned")

            create_synthetic_featurexml(p1, n_features=25, rt_offset=0.0)
            create_synthetic_featurexml(p2, n_features=25, rt_offset=20.0)

            align_kd([p1, p2], out_dir)

            out1 = os.path.join(out_dir, "run1.featureXML")
            out2 = os.path.join(out_dir, "run2.featureXML")

            fm1 = oms.FeatureMap()
            fm2 = oms.FeatureMap()
            oms.FeatureXMLFile().load(out1, fm1)
            oms.FeatureXMLFile().load(out2, fm2)

            assert fm1.size() == 25
            assert fm2.size() == 25

    def test_creates_output_directory(self):
        """Output directory should be created if it does not exist."""
        from rt_alignment_kd import align_kd, create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            p1 = os.path.join(tmp, "run1.featureXML")
            p2 = os.path.join(tmp, "run2.featureXML")
            out_dir = os.path.join(tmp, "new_dir", "aligned")

            create_synthetic_featurexml(p1, n_features=30, rt_offset=0.0)
            create_synthetic_featurexml(p2, n_features=30, rt_offset=10.0)

            align_kd([p1, p2], out_dir)

            assert os.path.isdir(out_dir)

    def test_create_synthetic_featurexml(self):
        """Synthetic FeatureXML should have correct number of features."""
        import pyopenms as oms
        from rt_alignment_kd import create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "test.featureXML")
            create_synthetic_featurexml(path, n_features=15, rt_offset=50.0)

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(path, fm)
            assert fm.size() == 15
            assert fm[0].getRT() == pytest.approx(150.0, abs=0.01)
