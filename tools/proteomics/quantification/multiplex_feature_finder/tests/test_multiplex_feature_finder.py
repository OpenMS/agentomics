"""Tests for multiplex_feature_finder."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestKnownLabels:
    def test_known_labels_dict(self):
        from multiplex_feature_finder import KNOWN_LABELS

        assert "Lys8" in KNOWN_LABELS
        assert abs(KNOWN_LABELS["Lys8"] - 8.0142) < 0.01
        assert "Arg10" in KNOWN_LABELS


class TestCreateSyntheticMzML:
    def test_creates_valid_mzml(self):
        import pyopenms as oms
        from multiplex_feature_finder import create_synthetic_multiplex_mzml

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            create_synthetic_multiplex_mzml(mzml_path, n_scans=10)

            exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, exp)
            assert exp.size() == 10
            assert exp[0].getMSLevel() == 1

    def test_contains_light_and_heavy_peaks(self):
        import pyopenms as oms
        from multiplex_feature_finder import create_synthetic_multiplex_mzml

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            light_mz = 500.0
            mass_shift = 8.0142
            charge = 2
            create_synthetic_multiplex_mzml(
                mzml_path,
                light_mz=light_mz,
                mass_shift=mass_shift,
                charge=charge,
                n_scans=5,
            )

            exp = oms.MSExperiment()
            oms.MzMLFile().load(mzml_path, exp)

            # Check the center scan has peaks near light and heavy m/z
            center_spec = exp[2]
            mzs = [center_spec[i].getMZ() for i in range(center_spec.size())]

            heavy_mz = light_mz + mass_shift / charge
            has_light = any(abs(m - light_mz) < 0.1 for m in mzs)
            has_heavy = any(abs(m - heavy_mz) < 0.1 for m in mzs)
            assert has_light, f"No peak near light m/z {light_mz}"
            assert has_heavy, f"No peak near heavy m/z {heavy_mz}"


class TestFindMultiplexFeatures:
    def test_runs_without_error(self):
        from multiplex_feature_finder import (
            create_synthetic_multiplex_mzml,
            find_multiplex_features,
        )

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            out_path = os.path.join(tmp, "features.featureXML")

            create_synthetic_multiplex_mzml(
                mzml_path, n_scans=30, intensity=1e6
            )

            n_features = find_multiplex_features(
                mzml_path,
                out_path,
                labels="[][Lys8]",
                intensity_cutoff=100.0,
            )

            assert isinstance(n_features, int)
            assert n_features >= 0
            assert os.path.exists(out_path)

    def test_output_is_valid_featurexml(self):
        import pyopenms as oms
        from multiplex_feature_finder import (
            create_synthetic_multiplex_mzml,
            find_multiplex_features,
        )

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "test.mzML")
            out_path = os.path.join(tmp, "features.featureXML")

            create_synthetic_multiplex_mzml(
                mzml_path, n_scans=30, intensity=1e6
            )

            find_multiplex_features(
                mzml_path,
                out_path,
                labels="[][Lys8]",
                intensity_cutoff=100.0,
            )

            fm = oms.FeatureMap()
            oms.FeatureXMLFile().load(out_path, fm)
            # Just check we can load it; feature count depends on algorithm
            assert isinstance(fm.size(), int)
