"""Tests for mass_defect_filter."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestMassDefectFilter:
    def test_compute_mass_defect(self):
        from mass_defect_filter import compute_mass_defect

        md = compute_mass_defect(180.0634)
        assert 0.0 < md < 1.0
        assert abs(md - 0.0634) < 0.001

    def test_compute_kendrick_mass(self):
        from mass_defect_filter import compute_kendrick_mass

        km = compute_kendrick_mass(180.0634, "CH2")
        assert km > 0

    def test_compute_kendrick_mass_defect(self):
        from mass_defect_filter import compute_kendrick_mass_defect

        kmd = compute_kendrick_mass_defect(180.0634, "CH2")
        assert -1.0 < kmd < 1.0

    def test_filter_by_mass_defect(self):
        from mass_defect_filter import filter_by_mass_defect

        features = [
            {"exact_mass": "180.0634", "name": "glucose"},
            {"exact_mass": "342.1162", "name": "sucrose"},
            {"exact_mass": "100.9000", "name": "test"},
        ]
        results = filter_by_mass_defect(features, mdf_min=0.0, mdf_max=0.2)
        assert len(results) >= 1
        for r in results:
            assert "mass_defect" in r
            assert "kendrick_mass_defect" in r

    def test_filter_excludes(self):
        from mass_defect_filter import filter_by_mass_defect

        features = [{"exact_mass": "180.0634"}]
        results = filter_by_mass_defect(features, mdf_min=0.5, mdf_max=0.9)
        assert len(results) == 0

    def test_read_write_tsv(self):
        from mass_defect_filter import read_features_tsv, write_tsv

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "features.tsv")
            with open(in_path, "w") as f:
                f.write("exact_mass\tname\n")
                f.write("180.0634\tglucose\n")
            features = read_features_tsv(in_path)
            assert len(features) == 1
            assert features[0]["exact_mass"] == "180.0634"

            out_path = os.path.join(tmpdir, "filtered.tsv")
            write_tsv([{"exact_mass": "180.0634", "mass_defect": 0.0634, "kendrick_mass_defect": 0.05}], out_path)
            assert os.path.exists(out_path)
