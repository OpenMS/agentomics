"""Tests for precursor_charge_distribution."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestPrecursorChargeDistribution:
    def test_charge_distribution(self):
        from precursor_charge_distribution import analyze_charge_distribution, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, charge_dist={2: 50, 3: 30, 4: 20})
            results = analyze_charge_distribution(mzml_path)
            assert len(results) == 3
            charge_map = {r["charge"]: r["count"] for r in results}
            assert charge_map[2] == 50
            assert charge_map[3] == 30
            assert charge_map[4] == 20

    def test_percentages(self):
        from precursor_charge_distribution import analyze_charge_distribution, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path, charge_dist={2: 50, 3: 50})
            results = analyze_charge_distribution(mzml_path)
            total_pct = sum(r["percentage"] for r in results)
            assert abs(total_pct - 100.0) < 0.1

    def test_result_keys(self):
        from precursor_charge_distribution import analyze_charge_distribution, create_synthetic_mzml

        with tempfile.TemporaryDirectory() as tmpdir:
            mzml_path = os.path.join(tmpdir, "test.mzML")
            create_synthetic_mzml(mzml_path)
            results = analyze_charge_distribution(mzml_path)
            for r in results:
                assert "charge" in r
                assert "count" in r
                assert "percentage" in r

    def test_write_tsv(self):
        from precursor_charge_distribution import write_tsv

        results = [{"charge": 2, "count": 50, "percentage": 50.0}]
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "charge_dist.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
