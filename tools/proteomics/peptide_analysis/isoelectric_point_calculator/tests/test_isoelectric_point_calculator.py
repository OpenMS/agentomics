"""Tests for isoelectric_point_calculator."""

import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestIsoelectricPointCalculator:
    def test_basic_pi(self):
        from isoelectric_point_calculator import calculate_pi

        pi = calculate_pi("ACDEFGHIK")
        assert 4.0 < pi < 7.0

    def test_acidic_peptide(self):
        from isoelectric_point_calculator import calculate_pi

        pi = calculate_pi("DDDDDDDK")
        assert pi < 5.0

    def test_basic_peptide(self):
        from isoelectric_point_calculator import calculate_pi

        pi = calculate_pi("KKKKKKKK")
        assert pi > 9.0

    def test_charge_at_ph(self):
        from isoelectric_point_calculator import charge_at_ph

        # At very low pH, everything is positive
        assert charge_at_ph("PEPTIDEK", 1.0) > 0
        # At very high pH, everything is negative
        assert charge_at_ph("PEPTIDEK", 14.0) < 0

    def test_different_pk_sets(self):
        from isoelectric_point_calculator import calculate_pi

        lehninger = calculate_pi("ACDEFGHIK", "lehninger")
        emboss = calculate_pi("ACDEFGHIK", "emboss")
        # Different pKa sets should give slightly different results
        assert isinstance(lehninger, float)
        assert isinstance(emboss, float)

    def test_charge_curve(self):
        from isoelectric_point_calculator import calculate_charge_curve

        curve = calculate_charge_curve("PEPTIDEK", "lehninger", 0.0, 14.0, 1.0)
        assert len(curve) == 15
        # Charge should decrease with pH
        assert curve[0]["charge"] > curve[-1]["charge"]

    def test_calculate_pi_from_sequence(self):
        from isoelectric_point_calculator import calculate_pi_from_sequence

        result = calculate_pi_from_sequence("PEPTIDEK", "lehninger")
        assert result["pk_set"] == "lehninger"
        assert result["length"] == 8
        assert abs(result["charge_at_pI"]) < 0.1  # charge near 0 at pI

    def test_fasta_processing(self):
        import pyopenms as oms
        from isoelectric_point_calculator import calculate_pi_from_sequence

        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = f"{tmpdir}/test.fasta"
            entries = []
            e1 = oms.FASTAEntry()
            e1.identifier = "PROT1"
            e1.sequence = "MSPEPTIDEKAAANOTHERPEPTIDE"
            entries.append(e1)
            oms.FASTAFile().store(fasta_path, entries)

            # Process single sequence same as would be done for FASTA entry
            result = calculate_pi_from_sequence(e1.sequence)
            assert result["pI"] > 0
