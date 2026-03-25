"""Tests for rna_fragment_spectrum_generator."""

import pytest

pytest.importorskip("pyopenms")

from rna_fragment_spectrum_generator import (
    generate_a_minus_b_ions,
    generate_all_fragments,
    generate_c_ions,
    generate_w_ions,
    generate_y_ions,
)


class TestRnaFragmentSpectrumGenerator:
    def test_c_ions_count(self):
        ions = generate_c_ions("AAUGC", charge=1)
        assert len(ions) == 4  # n-1 ions for length 5
        assert all(label.startswith("c") for label, _ in ions)

    def test_y_ions_count(self):
        ions = generate_y_ions("AAUGC", charge=1)
        assert len(ions) == 4
        assert all(label.startswith("y") for label, _ in ions)

    def test_w_ions_count(self):
        ions = generate_w_ions("AAUGC", charge=1)
        assert len(ions) == 3  # from index 2 to n-1
        assert all(label.startswith("w") for label, _ in ions)

    def test_a_minus_b_ions_count(self):
        ions = generate_a_minus_b_ions("AAUGC", charge=1)
        assert len(ions) == 3
        assert all(label.startswith("a-B") for label, _ in ions)

    def test_all_fragments(self):
        fragments = generate_all_fragments("AAUGC", charge=1)
        assert len(fragments) == 14  # 4 + 4 + 3 + 3
        ion_types = {f["ion_type"] for f in fragments}
        assert ion_types == {"c", "y", "w", "a-B"}

    def test_charge_state_affects_mz(self):
        f1 = generate_all_fragments("AAUGC", charge=1)
        f2 = generate_all_fragments("AAUGC", charge=2)
        # Same ion in charge 2 should have lower m/z than charge 1
        c1_mz = next(f["mz"] for f in f1 if f["ion_label"] == "c1")
        c1_mz_z2 = next(f["mz"] for f in f2 if f["ion_label"] == "c1")
        assert c1_mz_z2 < c1_mz

    def test_all_mz_positive(self):
        fragments = generate_all_fragments("AAUGC", charge=1)
        assert all(f["mz"] > 0 for f in fragments)

    def test_invalid_nucleotide(self):
        with pytest.raises(ValueError, match="Invalid RNA nucleotide"):
            generate_all_fragments("AATGC")

    def test_short_sequence(self):
        fragments = generate_all_fragments("AU", charge=1)
        # c: 1, y: 1, w: 0, a-B: 0
        assert len(fragments) == 2
