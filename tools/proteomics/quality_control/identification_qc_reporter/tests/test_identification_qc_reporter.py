"""Tests for identification_qc_reporter."""

from conftest import requires_pyopenms


@requires_pyopenms
class TestIdentificationQcReporter:
    def _make_rows(self):
        return [
            {"sequence": "PEPTIDEK", "charge": "2", "score": "0.95",
             "precursor_mz": "464.75", "modifications": ""},
            {"sequence": "ACDEFGHIK", "charge": "2", "score": "0.88",
             "precursor_mz": "502.73", "modifications": "Oxidation"},
            {"sequence": "PEPTIDEK", "charge": "3", "score": "0.72",
             "precursor_mz": "310.17", "modifications": ""},
            {"sequence": "YLAGNK", "charge": "2", "score": "0.65",
             "precursor_mz": "332.19", "modifications": "Oxidation;Phospho"},
        ]

    def test_psm_count(self):
        from identification_qc_reporter import compute_id_qc

        metrics = compute_id_qc(self._make_rows())
        assert metrics["psm_count"] == 4

    def test_unique_peptides(self):
        from identification_qc_reporter import compute_id_qc

        metrics = compute_id_qc(self._make_rows())
        assert metrics["unique_peptides"] == 3

    def test_score_range(self):
        from identification_qc_reporter import compute_id_qc

        metrics = compute_id_qc(self._make_rows())
        assert metrics["score_min"] == 0.65
        assert metrics["score_max"] == 0.95

    def test_modification_counts(self):
        from identification_qc_reporter import compute_id_qc

        metrics = compute_id_qc(self._make_rows())
        assert metrics["modification_counts"]["Oxidation"] == 2
        assert metrics["modification_counts"]["Phospho"] == 1

    def test_empty_rows(self):
        from identification_qc_reporter import compute_id_qc

        metrics = compute_id_qc([])
        assert metrics["psm_count"] == 0

    def test_mass_error_computed(self):
        from identification_qc_reporter import compute_id_qc

        metrics = compute_id_qc(self._make_rows())
        assert "mass_error_mean_ppm" in metrics
