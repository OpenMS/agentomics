"""Tests for kovats_ri_calculator."""

import csv
import math
import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestKovatsRiCalculator:
    def test_build_alkane_table(self):
        from kovats_ri_calculator import build_alkane_table

        standards = [
            {"carbon_number": 10, "rt": 5.0},
            {"carbon_number": 12, "rt": 10.0},
            {"carbon_number": 8, "rt": 2.0},
        ]
        table = build_alkane_table(standards)
        assert len(table) == 3
        # Should be sorted by RT
        assert table[0][1] < table[1][1] < table[2][1]

    def test_calculate_kovats_ri_at_alkane(self):
        from kovats_ri_calculator import calculate_kovats_ri

        # At exactly C10 RT, RI should be ~1000
        alkane_table = [(8, 2.0), (10, 5.0), (12, 10.0)]
        ri = calculate_kovats_ri(5.0, alkane_table)
        assert ri is not None
        assert abs(ri - 1000.0) < 0.1

    def test_calculate_kovats_ri_between_alkanes(self):
        from kovats_ri_calculator import calculate_kovats_ri

        alkane_table = [(10, 5.0), (12, 10.0)]
        ri = calculate_kovats_ri(7.0, alkane_table)
        assert ri is not None
        # RI should be between 1000 and 1200
        assert 1000.0 < ri < 1200.0

    def test_calculate_kovats_ri_formula(self):
        from kovats_ri_calculator import calculate_kovats_ri

        # Verify the formula: RI = 100*n + 100*(log(RT_x) - log(RT_n))/(log(RT_n+1) - log(RT_n))
        alkane_table = [(10, 5.0), (12, 10.0)]
        rt_x = 7.0
        expected = 100 * 10 + 100 * (math.log10(7.0) - math.log10(5.0)) / (math.log10(10.0) - math.log10(5.0))
        ri = calculate_kovats_ri(rt_x, alkane_table)
        assert ri is not None
        assert abs(ri - round(expected, 2)) < 0.01

    def test_calculate_kovats_ri_out_of_range(self):
        from kovats_ri_calculator import calculate_kovats_ri

        alkane_table = [(10, 5.0), (12, 10.0)]
        # RT before first alkane
        assert calculate_kovats_ri(1.0, alkane_table) is None
        # RT after last alkane
        assert calculate_kovats_ri(15.0, alkane_table) is None

    def test_calculate_kovats_ri_negative_rt(self):
        from kovats_ri_calculator import calculate_kovats_ri

        alkane_table = [(10, 5.0), (12, 10.0)]
        assert calculate_kovats_ri(-1.0, alkane_table) is None

    def test_calculate_ri_batch(self):
        from kovats_ri_calculator import calculate_ri_batch

        alkane_table = [(8, 2.0), (10, 5.0), (12, 10.0)]
        features = [
            {"feature_id": "F1", "rt": 3.0},
            {"feature_id": "F2", "rt": 7.0},
            {"feature_id": "F3", "rt": 20.0},
        ]
        results = calculate_ri_batch(features, alkane_table)
        assert len(results) == 3
        assert results[0]["kovats_ri"] != "N/A"
        assert results[1]["kovats_ri"] != "N/A"
        assert results[2]["kovats_ri"] == "N/A"

    def test_full_pipeline(self):
        from kovats_ri_calculator import build_alkane_table, calculate_ri_batch, load_tsv, write_results

        with tempfile.TemporaryDirectory() as tmpdir:
            std_path = os.path.join(tmpdir, "alkanes.tsv")
            with open(std_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, ["carbon_number", "rt"], delimiter="\t")
                w.writeheader()
                w.writerow({"carbon_number": "8", "rt": "2.0"})
                w.writerow({"carbon_number": "10", "rt": "5.0"})
                w.writerow({"carbon_number": "12", "rt": "10.0"})

            feat_path = os.path.join(tmpdir, "features.tsv")
            with open(feat_path, "w", newline="") as fh:
                w = csv.DictWriter(fh, ["feature_id", "rt"], delimiter="\t")
                w.writeheader()
                w.writerow({"feature_id": "F1", "rt": "3.5"})

            standards = load_tsv(std_path)
            features = load_tsv(feat_path)
            alkane_table = build_alkane_table(standards)
            results = calculate_ri_batch(features, alkane_table)

            out_path = os.path.join(tmpdir, "ri.tsv")
            write_results(results, out_path)
            assert os.path.exists(out_path)
