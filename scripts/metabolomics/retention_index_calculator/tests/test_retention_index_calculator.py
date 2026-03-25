"""Tests for retention_index_calculator."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestRetentionIndexCalculator:
    def _make_standards(self):
        return [
            (8, 100.0),   # C8 at 100s
            (9, 200.0),   # C9 at 200s
            (10, 400.0),  # C10 at 400s
            (11, 800.0),  # C11 at 800s
        ]

    def test_exact_standard_ri(self):
        from retention_index_calculator import calculate_kovats_ri

        standards = self._make_standards()
        # At the C9 standard RT, RI should be 900
        ri = calculate_kovats_ri(200.0, standards)
        assert ri == 900.0

    def test_interpolation(self):
        from retention_index_calculator import calculate_kovats_ri

        standards = self._make_standards()
        ri = calculate_kovats_ri(150.0, standards)
        assert ri is not None
        assert 800 < ri < 900

    def test_out_of_range(self):
        from retention_index_calculator import calculate_kovats_ri

        standards = self._make_standards()
        ri = calculate_kovats_ri(50.0, standards)  # before first standard
        assert ri is None

    def test_calculate_all_ri(self):
        from retention_index_calculator import calculate_all_ri

        standards = self._make_standards()
        features = [
            {"mz": "100.0", "rt": "150.0", "intensity": "1000"},
            {"mz": "200.0", "rt": "300.0", "intensity": "2000"},
            {"mz": "300.0", "rt": "50.0", "intensity": "500"},  # out of range
        ]
        results = calculate_all_ri(features, standards)
        assert len(results) == 3
        assert results[0]["retention_index"] != ""
        assert results[1]["retention_index"] != ""
        assert results[2]["retention_index"] == ""

    def test_monotonic_ri(self):
        from retention_index_calculator import calculate_kovats_ri

        standards = self._make_standards()
        ri1 = calculate_kovats_ri(150.0, standards)
        ri2 = calculate_kovats_ri(300.0, standards)
        assert ri1 < ri2
