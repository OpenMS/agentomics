"""Tests for id_mapper."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestCreateSyntheticFeatureXML:
    def test_creates_valid_file(self):
        import pyopenms as oms
        from id_mapper import create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "features.featureXML")
            create_synthetic_featurexml(path, mz=500.0, rt=100.0)

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(path, fm)
            assert fm.size() == 1
            assert abs(fm[0].getMZ() - 500.0) < 0.01
            assert abs(fm[0].getRT() - 100.0) < 0.01


class TestCreateSyntheticIdXML:
    def test_creates_valid_file(self):
        import pyopenms as oms
        from id_mapper import create_synthetic_idxml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "ids.idXML")
            create_synthetic_idxml(path, mz=500.0, rt=100.0, sequence="PEPTIDEK")

            protein_ids = []
            peptide_ids = oms.PeptideIdentificationList()
            oms.IdXMLFile().load(path, protein_ids, peptide_ids)
            assert peptide_ids.size() == 1
            assert len(peptide_ids[0].getHits()) == 1
            assert peptide_ids[0].getHits()[0].getSequence().toString() == "PEPTIDEK"


class TestMapIdsToFeatures:
    def test_annotation_within_tolerance(self):
        import pyopenms as oms
        from id_mapper import (
            create_synthetic_featurexml,
            create_synthetic_idxml,
            map_ids_to_features,
        )

        with tempfile.TemporaryDirectory() as tmp:
            feat_path = os.path.join(tmp, "features.featureXML")
            ids_path = os.path.join(tmp, "ids.idXML")
            out_path = os.path.join(tmp, "annotated.featureXML")

            # Feature at RT=100, m/z=500
            create_synthetic_featurexml(feat_path, mz=500.0, rt=100.0)
            # Peptide ID at RT=100.5, m/z=500.001 (within tolerance)
            create_synthetic_idxml(ids_path, mz=500.001, rt=100.5)

            n_annotated = map_ids_to_features(
                feat_path, ids_path, out_path, rt_tol=5.0, mz_tol=10.0
            )

            assert n_annotated == 1

            # Verify the annotation
            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(out_path, fm)
            pids = fm[0].getPeptideIdentifications()
            assert len(pids) >= 1

    def test_no_annotation_outside_tolerance(self):
        from id_mapper import (
            create_synthetic_featurexml,
            create_synthetic_idxml,
            map_ids_to_features,
        )

        with tempfile.TemporaryDirectory() as tmp:
            feat_path = os.path.join(tmp, "features.featureXML")
            ids_path = os.path.join(tmp, "ids.idXML")
            out_path = os.path.join(tmp, "annotated.featureXML")

            # Feature at RT=100, m/z=500
            create_synthetic_featurexml(feat_path, mz=500.0, rt=100.0)
            # Peptide ID far away at RT=200, m/z=600
            create_synthetic_idxml(ids_path, mz=600.0, rt=200.0)

            n_annotated = map_ids_to_features(
                feat_path, ids_path, out_path, rt_tol=5.0, mz_tol=10.0
            )

            assert n_annotated == 0

    def test_output_file_created(self):
        from id_mapper import (
            create_synthetic_featurexml,
            create_synthetic_idxml,
            map_ids_to_features,
        )

        with tempfile.TemporaryDirectory() as tmp:
            feat_path = os.path.join(tmp, "features.featureXML")
            ids_path = os.path.join(tmp, "ids.idXML")
            out_path = os.path.join(tmp, "annotated.featureXML")

            create_synthetic_featurexml(feat_path)
            create_synthetic_idxml(ids_path)

            map_ids_to_features(feat_path, ids_path, out_path)
            assert os.path.exists(out_path)
