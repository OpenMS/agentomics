"""Tests for id_conflict_resolver."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestCreateSyntheticConflictingFeatureXML:
    def test_creates_valid_file(self):
        import pyopenms as oms
        from id_conflict_resolver import create_synthetic_conflicting_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "features.featureXML")
            create_synthetic_conflicting_featurexml(path)

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(path, fm)
            assert fm.size() == 1

    def test_has_multiple_ids(self):
        import pyopenms as oms
        from id_conflict_resolver import create_synthetic_conflicting_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "features.featureXML")
            create_synthetic_conflicting_featurexml(path)

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(path, fm)
            pids = fm[0].getPeptideIdentifications()
            assert pids.size() == 2

    def test_custom_sequences_and_scores(self):
        import pyopenms as oms
        from id_conflict_resolver import create_synthetic_conflicting_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "features.featureXML")
            create_synthetic_conflicting_featurexml(
                path,
                sequences=["AAA", "BBB", "CCC"],
                scores=[0.8, 0.6, 0.3],
            )

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(path, fm)
            pids = fm[0].getPeptideIdentifications()
            assert pids.size() == 3


class TestResolveIdConflicts:
    def test_keeps_best_hit(self):
        import pyopenms as oms
        from id_conflict_resolver import (
            create_synthetic_conflicting_featurexml,
            resolve_id_conflicts,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.featureXML")
            out_path = os.path.join(tmp, "resolved.featureXML")

            # Feature with 2 hits: scores 0.9 and 0.5
            create_synthetic_conflicting_featurexml(
                in_path, scores=[0.9, 0.5]
            )

            n = resolve_id_conflicts(in_path, out_path)
            assert n == 1

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(out_path, fm)

            # After resolution, should have only 1 identification
            pids = fm[0].getPeptideIdentifications()
            assert pids.size() == 1

            # The remaining ID should be the best one
            best_hit = pids[0].getHits()[0]
            assert abs(best_hit.getScore() - 0.9) < 0.01
            assert best_hit.getSequence().toString() == "PEPTIDEK"

    def test_output_file_created(self):
        from id_conflict_resolver import (
            create_synthetic_conflicting_featurexml,
            resolve_id_conflicts,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.featureXML")
            out_path = os.path.join(tmp, "resolved.featureXML")

            create_synthetic_conflicting_featurexml(in_path)
            resolve_id_conflicts(in_path, out_path)
            assert os.path.exists(out_path)

    def test_feature_without_conflicts(self):
        import pyopenms as oms
        from id_conflict_resolver import (
            create_synthetic_conflicting_featurexml,
            resolve_id_conflicts,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.featureXML")
            out_path = os.path.join(tmp, "resolved.featureXML")

            # Single ID, no conflict
            create_synthetic_conflicting_featurexml(
                in_path,
                sequences=["PEPTIDEK"],
                scores=[0.95],
            )

            n = resolve_id_conflicts(in_path, out_path)
            assert n == 1

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(out_path, fm)
            pids = fm[0].getPeptideIdentifications()
            assert pids.size() == 1

    def test_three_way_conflict(self):
        import pyopenms as oms
        from id_conflict_resolver import (
            create_synthetic_conflicting_featurexml,
            resolve_id_conflicts,
        )

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "input.featureXML")
            out_path = os.path.join(tmp, "resolved.featureXML")

            create_synthetic_conflicting_featurexml(
                in_path,
                sequences=["AAA", "BBB", "CCC"],
                scores=[0.3, 0.8, 0.5],
            )

            resolve_id_conflicts(in_path, out_path)

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(out_path, fm)
            pids = fm[0].getPeptideIdentifications()
            assert pids.size() == 1
            best_hit = pids[0].getHits()[0]
            assert abs(best_hit.getScore() - 0.8) < 0.01
