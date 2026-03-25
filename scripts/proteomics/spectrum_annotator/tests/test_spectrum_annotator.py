"""Tests for spectrum_annotator."""

import os
import tempfile

from conftest import requires_pyopenms


@requires_pyopenms
class TestSpectrumAnnotator:
    def test_annotate_with_matching_peaks(self):
        import pyopenms as oms
        from spectrum_annotator import annotate_spectrum

        # Generate theoretical spectrum to get known m/z values
        aa_seq = oms.AASequence.fromString("PEPTIDEK")
        tsg = oms.TheoreticalSpectrumGenerator()
        param = tsg.getParameters()
        param.setValue("add_b_ions", "true")
        param.setValue("add_y_ions", "true")
        param.setValue("add_metainfo", "true")
        tsg.setParameters(param)
        spec = oms.MSSpectrum()
        tsg.getSpectrum(spec, aa_seq, 1, 1)
        theo_mzs, _ = spec.get_peaks()

        # Use first few theoretical m/z values as observed
        obs_mz = [float(theo_mzs[0]), float(theo_mzs[1])]
        obs_int = [1000.0, 500.0]

        results = annotate_spectrum(obs_mz, obs_int, "PEPTIDEK", charge=1, tolerance=0.05)
        assert len(results) == 2
        # At least one should match
        matched = [r for r in results if r["matched_ion"]]
        assert len(matched) > 0

    def test_annotate_no_match(self):
        from spectrum_annotator import annotate_spectrum

        results = annotate_spectrum([9999.0], [1000.0], "PEPTIDEK", charge=1, tolerance=0.02)
        assert len(results) == 1
        assert results[0]["matched_ion"] == ""

    def test_result_keys(self):
        from spectrum_annotator import annotate_spectrum

        results = annotate_spectrum([100.0, 200.0], [1000.0, 500.0], "PEPTIDEK", charge=1)
        for r in results:
            assert "observed_mz" in r
            assert "intensity" in r
            assert "matched_ion" in r
            assert "theoretical_mz" in r
            assert "error_da" in r

    def test_write_tsv(self):
        from spectrum_annotator import annotate_spectrum, write_tsv

        results = annotate_spectrum([100.0], [1000.0], "PEPTIDEK", charge=1)
        with tempfile.TemporaryDirectory() as tmpdir:
            out = os.path.join(tmpdir, "annotation.tsv")
            write_tsv(results, out)
            assert os.path.exists(out)
            with open(out) as f:
                lines = f.readlines()
            assert len(lines) == 2
            assert "observed_mz" in lines[0]
