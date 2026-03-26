"""Tests for accurate_mass_searcher."""

import os
import tempfile

import pytest

pyopenms = pytest.importorskip("pyopenms")


class TestAccurateMassSearcher:
    """Tests for accurate mass searching functionality."""

    def test_create_synthetic_featurexml(self):
        from accurate_mass_searcher import create_synthetic_featurexml

        with tempfile.TemporaryDirectory() as tmp:
            fxml_path = os.path.join(tmp, "features.featureXML")
            create_synthetic_featurexml(fxml_path)

            fm = pyopenms.FeatureMap()
            pyopenms.FeatureXMLFile().load(fxml_path, fm)
            assert fm.size() == 1
            feat = fm[0]
            assert abs(feat.getMZ() - 181.070664) < 0.001
            assert feat.getCharge() == 1

    def test_create_synthetic_database(self):
        from accurate_mass_searcher import create_synthetic_database

        with tempfile.TemporaryDirectory() as tmp:
            mapping_path = os.path.join(tmp, "mapping.tsv")
            struct_path = os.path.join(tmp, "struct.tsv")
            create_synthetic_database(mapping_path, struct_path)

            with open(mapping_path) as f:
                content = f.read()
            assert "180.063388" in content
            assert "C6H12O6" in content
            assert "HMDB:HMDB0000122" in content

            with open(struct_path) as f:
                content = f.read()
            assert "Glucose" in content

    def test_search_accurate_mass_finds_glucose(self):
        """Verify that glucose [M+H]+ at m/z ~181.07 matches glucose in database."""
        from accurate_mass_searcher import (
            create_synthetic_database,
            create_synthetic_featurexml,
            search_accurate_mass,
        )

        # Glucose [M+H]+: mass 180.063388 + proton 1.007276 = 181.070664
        with tempfile.TemporaryDirectory() as tmp:
            fxml_path = os.path.join(tmp, "features.featureXML")
            mapping_path = os.path.join(tmp, "mapping.tsv")
            struct_path = os.path.join(tmp, "struct.tsv")
            output_path = os.path.join(tmp, "results.mzTab")

            create_synthetic_featurexml(fxml_path, mz=181.070664)
            create_synthetic_database(mapping_path, struct_path)

            matches = search_accurate_mass(
                fxml_path, mapping_path, output_path, mass_tol=5.0, struct_path=struct_path
            )

            assert matches >= 1
            assert os.path.exists(output_path)

            # Verify the mzTab contains glucose match
            with open(output_path) as f:
                content = f.read()
            assert "HMDB:HMDB0000122" in content
            assert "C6H12O6" in content

    def test_search_no_match_with_wrong_mass(self):
        """Verify that a feature far from glucose mass does not match."""
        from accurate_mass_searcher import (
            create_synthetic_database,
            create_synthetic_featurexml,
            search_accurate_mass,
        )

        with tempfile.TemporaryDirectory() as tmp:
            fxml_path = os.path.join(tmp, "features.featureXML")
            mapping_path = os.path.join(tmp, "mapping.tsv")
            struct_path = os.path.join(tmp, "struct.tsv")
            output_path = os.path.join(tmp, "results.mzTab")

            # Use m/z=999.0 which won't match glucose
            create_synthetic_featurexml(fxml_path, mz=999.0)
            create_synthetic_database(mapping_path, struct_path)

            matches = search_accurate_mass(
                fxml_path, mapping_path, output_path, mass_tol=5.0, struct_path=struct_path
            )

            assert matches == 0
