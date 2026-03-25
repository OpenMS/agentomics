"""Tests for mid_natural_abundance_corrector."""

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestGetNumTracerAtoms:
    def test_glucose_carbons(self):
        from mid_natural_abundance_corrector import get_num_tracer_atoms

        assert get_num_tracer_atoms("C6H12O6", "13C") == 6

    def test_alanine_nitrogens(self):
        from mid_natural_abundance_corrector import get_num_tracer_atoms

        assert get_num_tracer_atoms("C3H7NO2", "15N") == 1

    def test_unsupported_tracer(self):
        from mid_natural_abundance_corrector import get_num_tracer_atoms

        with pytest.raises(ValueError, match="Unsupported tracer"):
            get_num_tracer_atoms("C6H12O6", "18O")


@requires_pyopenms
class TestBuildCorrectionMatrix:
    def test_matrix_shape(self):
        from mid_natural_abundance_corrector import build_correction_matrix

        C = build_correction_matrix(6)
        assert C.shape == (7, 7)

    def test_columns_sum_to_one(self):
        import numpy as np
        from mid_natural_abundance_corrector import build_correction_matrix

        C = build_correction_matrix(6)
        col_sums = C.sum(axis=0)
        np.testing.assert_allclose(col_sums, 1.0, atol=1e-10)

    def test_diagonal_dominant(self):
        """Natural abundance is small, so diagonal should dominate."""
        from mid_natural_abundance_corrector import build_correction_matrix

        C = build_correction_matrix(6)
        for i in range(7):
            assert C[i, i] > 0.5


@requires_pyopenms
class TestCorrectMID:
    def test_unlabeled_sample(self):
        """An unlabeled sample should have ~1.0 at M+0 after correction."""
        from mid_natural_abundance_corrector import correct_mid

        # Simulated measured MID for unlabeled glucose (with natural 13C)
        measured = [0.935, 0.061, 0.004, 0.0, 0.0, 0.0, 0.0]
        corrected = correct_mid(measured, "C6H12O6", "13C")
        assert len(corrected) == 7
        assert corrected[0] > 0.95  # M+0 should be dominant after correction
        assert sum(corrected) == pytest.approx(1.0, abs=0.01)

    def test_fully_labeled_sample(self):
        """A fully 13C-labeled glucose should have high M+6."""
        from mid_natural_abundance_corrector import correct_mid

        measured = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        corrected = correct_mid(measured, "C6H12O6", "13C")
        assert corrected[6] > 0.9

    def test_output_sums_to_one(self):
        from mid_natural_abundance_corrector import correct_mid

        measured = [0.5, 0.2, 0.15, 0.1, 0.03, 0.015, 0.005]
        corrected = correct_mid(measured, "C6H12O6", "13C")
        assert sum(corrected) == pytest.approx(1.0, abs=0.01)

    def test_non_negative(self):
        from mid_natural_abundance_corrector import correct_mid

        measured = [0.9, 0.08, 0.02, 0.0, 0.0, 0.0, 0.0]
        corrected = correct_mid(measured, "C6H12O6", "13C")
        assert all(v >= 0 for v in corrected)
