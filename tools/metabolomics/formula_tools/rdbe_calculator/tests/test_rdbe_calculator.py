"""Tests for rdbe_calculator."""


from conftest import requires_pyopenms


@requires_pyopenms
class TestCalculateRDBE:
    def test_benzene(self):
        """C6H6: RDBE = (12 + 2 - 6) / 2 = 4"""
        from rdbe_calculator import calculate_rdbe

        assert calculate_rdbe("C6H6") == 4.0

    def test_methane(self):
        """CH4: RDBE = (2 + 2 - 4) / 2 = 0"""
        from rdbe_calculator import calculate_rdbe

        assert calculate_rdbe("CH4") == 0.0

    def test_glucose(self):
        """C6H12O6: RDBE = (12 + 2 - 12) / 2 = 1"""
        from rdbe_calculator import calculate_rdbe

        assert calculate_rdbe("C6H12O6") == 1.0

    def test_naphthalene(self):
        """C10H8: RDBE = (20 + 2 - 8) / 2 = 7"""
        from rdbe_calculator import calculate_rdbe

        assert calculate_rdbe("C10H8") == 7.0

    def test_alanine(self):
        """C3H7NO2: RDBE = (6 + 2 - 7 + 1) / 2 = 1"""
        from rdbe_calculator import calculate_rdbe

        assert calculate_rdbe("C3H7NO2") == 1.0

    def test_atp_with_phosphorus(self):
        """C10H16N5O13P3: RDBE = (20 + 2 - 16 + 5 + 3) / 2 = 7"""
        from rdbe_calculator import calculate_rdbe

        assert calculate_rdbe("C10H16N5O13P3") == 7.0

    def test_ethanol(self):
        """C2H6O: RDBE = (4 + 2 - 6) / 2 = 0"""
        from rdbe_calculator import calculate_rdbe

        assert calculate_rdbe("C2H6O") == 0.0


@requires_pyopenms
class TestGetElementCounts:
    def test_glucose(self):
        from rdbe_calculator import get_element_counts

        counts = get_element_counts("C6H12O6")
        assert counts["C"] == 6
        assert counts["H"] == 12
        assert counts["O"] == 6


@requires_pyopenms
class TestCalculateRDBEBatch:
    def test_batch(self):
        from rdbe_calculator import calculate_rdbe_batch

        results = calculate_rdbe_batch(["C6H6", "CH4", "C6H12O6"])
        assert len(results) == 3
        assert results[0]["rdbe"] == 4.0
        assert results[1]["rdbe"] == 0.0
        assert results[2]["rdbe"] == 1.0

    def test_batch_has_all_fields(self):
        from rdbe_calculator import calculate_rdbe_batch

        results = calculate_rdbe_batch(["C6H6"])
        r = results[0]
        assert "formula" in r
        assert "C" in r
        assert "H" in r
        assert "rdbe" in r
