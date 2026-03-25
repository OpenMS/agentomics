"""Tests for featurexml_merger."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
def test_create_synthetic_featurexml():
    import pyopenms as oms
    from featurexml_merger import create_synthetic_featurexml

    with tempfile.TemporaryDirectory() as tmp:
        fxml_path = os.path.join(tmp, "test.featureXML")
        create_synthetic_featurexml(fxml_path, n_features=3)

        fm = oms.FeatureMap()
        oms.FeatureXMLFile().load(fxml_path, fm)
        assert fm.size() == 3


@requires_pyopenms
def test_merge_two_files():
    from featurexml_merger import create_synthetic_featurexml, merge_feature_maps

    with tempfile.TemporaryDirectory() as tmp:
        f1 = os.path.join(tmp, "f1.featureXML")
        f2 = os.path.join(tmp, "f2.featureXML")
        out = os.path.join(tmp, "merged.featureXML")

        create_synthetic_featurexml(f1, n_features=3, rt_offset=0.0)
        create_synthetic_featurexml(f2, n_features=4, rt_offset=1000.0)

        stats = merge_feature_maps([f1, f2], out)
        assert stats["total_features"] == 7
        assert stats["file_counts"][f1] == 3
        assert stats["file_counts"][f2] == 4


@requires_pyopenms
def test_merged_sorted_by_rt():
    import pyopenms as oms
    from featurexml_merger import create_synthetic_featurexml, merge_feature_maps

    with tempfile.TemporaryDirectory() as tmp:
        f1 = os.path.join(tmp, "f1.featureXML")
        f2 = os.path.join(tmp, "f2.featureXML")
        out = os.path.join(tmp, "merged.featureXML")

        create_synthetic_featurexml(f1, n_features=3, rt_offset=500.0)
        create_synthetic_featurexml(f2, n_features=3, rt_offset=0.0)

        merge_feature_maps([f1, f2], out)

        fm = oms.FeatureMap()
        oms.FeatureXMLFile().load(out, fm)
        rts = [f.getRT() for f in fm]
        assert rts == sorted(rts)
