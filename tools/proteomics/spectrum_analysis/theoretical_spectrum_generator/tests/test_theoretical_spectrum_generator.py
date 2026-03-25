"""Tests for theoretical_spectrum_generator."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestTheoreticalSpectrumGenerator:
    def test_generate_by_ions(self):
        from theoretical_spectrum_generator import generate_theoretical_spectrum

        results = generate_theoretical_spectrum("PEPTIDEK", charge=1, ion_types=["b", "y"])
        assert len(results) > 0
        ion_types_found = {r["ion_type"] for r in results}
        assert "b" in ion_types_found or "y" in ion_types_found

    def test_multiple_ion_types(self):
        from theoretical_spectrum_generator import generate_theoretical_spectrum

        results = generate_theoretical_spectrum("PEPTIDEK", charge=1, ion_types=["b", "y", "a"])
        assert len(results) > 0

    def test_charge_state(self):
        from theoretical_spectrum_generator import generate_theoretical_spectrum

        r1 = generate_theoretical_spectrum("PEPTIDEK", charge=1, ion_types=["b", "y"])
        r2 = generate_theoretical_spectrum("PEPTIDEK", charge=2, ion_types=["b", "y"])
        assert len(r2) >= len(r1)

    def test_neutral_losses(self):
        from theoretical_spectrum_generator import generate_theoretical_spectrum

        r_no_loss = generate_theoretical_spectrum("PEPTIDEK", charge=1, ion_types=["b", "y"], add_losses=False)
        r_loss = generate_theoretical_spectrum("PEPTIDEK", charge=1, ion_types=["b", "y"], add_losses=True)
        assert len(r_loss) >= len(r_no_loss)

    def test_result_keys(self):
        from theoretical_spectrum_generator import generate_theoretical_spectrum

        results = generate_theoretical_spectrum("PEPTIDEK", charge=1)
        assert len(results) > 0
        for r in results:
            assert "ion_type" in r
            assert "ion_number" in r
            assert "charge" in r
            assert "mz" in r
            assert "annotation" in r

    def test_write_tsv(self):
        from theoretical_spectrum_generator import generate_theoretical_spectrum, write_tsv

        results = generate_theoretical_spectrum("PEPTIDEK", charge=1)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = os.path.join(tmpdir, "fragments.tsv")
            write_tsv(results, out_path)
            assert os.path.exists(out_path)
            with open(out_path) as f:
                lines = f.readlines()
            assert len(lines) > 1
            assert "ion_type" in lines[0]

    def test_parse_annotation(self):
        from theoretical_spectrum_generator import _parse_annotation

        result = _parse_annotation("y3+")
        assert result["ion_type"] == "y"
        assert result["ion_number"] == 3
        assert result["charge"] == 1
