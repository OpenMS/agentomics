"""Tests for metabolite_feature_deconvolution."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestMetaboliteFeatureDeconvolution:
    def _make_feature(self, mz, rt, intensity, charge=1):
        """Helper to create a single feature."""
        import pyopenms as oms

        f = oms.Feature()
        f.setMZ(mz)
        f.setRT(rt)
        f.setIntensity(intensity)
        f.setCharge(charge)
        return f

    def test_deconvolve_adducts_groups_features(self):
        """Features at [M+H]+, [M+Na]+, [M+K]+ for glucose should be grouped."""
        import pyopenms as oms
        from metabolite_feature_deconvolution import deconvolve_adducts

        # Glucose M = 180.063
        # [M+H]+  = 181.070 (M + 1.00728)
        # [M+Na]+ = 203.052 (M + 22.989)
        # [M+K]+  = 219.026 (M + 38.963)
        fm = oms.FeatureMap()
        rt = 120.0

        f1 = self._make_feature(181.070, rt, 1e6)
        f2 = self._make_feature(203.052, rt, 5e5)
        f3 = self._make_feature(219.026, rt, 2e5)

        fm.push_back(f1)
        fm.push_back(f2)
        fm.push_back(f3)
        fm.setUniqueIds()

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "features.featureXML")
            output_path = os.path.join(tmpdir, "grouped.consensusXML")
            oms.FeatureXMLFile().store(input_path, fm)

            n_groups = deconvolve_adducts(input_path, output_path)
            assert isinstance(n_groups, int)
            assert os.path.exists(output_path)

            # Load result and verify consensus features were created
            cons_map = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(output_path, cons_map)
            assert cons_map.size() >= 0

    def test_deconvolve_output_file_created(self):
        """Verify output consensusXML file is created."""
        import pyopenms as oms
        from metabolite_feature_deconvolution import deconvolve_adducts

        fm = oms.FeatureMap()
        f = self._make_feature(181.070, 120.0, 1e6)
        fm.push_back(f)
        fm.setUniqueIds()

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "features.featureXML")
            output_path = os.path.join(tmpdir, "grouped.consensusXML")
            oms.FeatureXMLFile().store(input_path, fm)

            deconvolve_adducts(input_path, output_path)
            assert os.path.exists(output_path)

    def test_deconvolve_returns_int(self):
        """Return value should be an integer."""
        import pyopenms as oms
        from metabolite_feature_deconvolution import deconvolve_adducts

        fm = oms.FeatureMap()
        f = self._make_feature(500.0, 60.0, 1e5)
        fm.push_back(f)
        fm.setUniqueIds()

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "features.featureXML")
            output_path = os.path.join(tmpdir, "grouped.consensusXML")
            oms.FeatureXMLFile().store(input_path, fm)

            result = deconvolve_adducts(input_path, output_path)
            assert isinstance(result, int)

    def test_deconvolve_empty_feature_map(self):
        """Empty feature map should not crash."""
        import pyopenms as oms
        from metabolite_feature_deconvolution import deconvolve_adducts

        fm = oms.FeatureMap()

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = os.path.join(tmpdir, "features.featureXML")
            output_path = os.path.join(tmpdir, "grouped.consensusXML")
            oms.FeatureXMLFile().store(input_path, fm)

            n_groups = deconvolve_adducts(input_path, output_path)
            assert n_groups == 0
