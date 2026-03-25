"""Tests for spectral_library_format_converter."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestSpectralLibraryFormatConverter:
    def test_parse_msp(self):
        from spectral_library_format_converter import create_synthetic_msp, parse_msp

        with tempfile.TemporaryDirectory() as tmpdir:
            msp_path = os.path.join(tmpdir, "test.msp")
            create_synthetic_msp(msp_path, n_spectra=3)
            spectra = parse_msp(msp_path)
            assert len(spectra) == 3
            assert spectra[0]["charge"] == 2
            assert len(spectra[0]["peaks"]) == 5

    def test_msp_to_targeted_experiment(self):
        from spectral_library_format_converter import create_synthetic_msp, msp_to_targeted_experiment, parse_msp

        with tempfile.TemporaryDirectory() as tmpdir:
            msp_path = os.path.join(tmpdir, "test.msp")
            create_synthetic_msp(msp_path, n_spectra=2)
            spectra = parse_msp(msp_path)
            te = msp_to_targeted_experiment(spectra)
            assert len(te.getProteins()) == 2
            assert len(te.getPeptides()) == 2
            assert len(te.getTransitions()) == 10  # 2 spectra * 5 peaks

    def test_convert_msp_to_traml(self):
        import pyopenms as oms
        from spectral_library_format_converter import convert_msp_to_traml, create_synthetic_msp

        with tempfile.TemporaryDirectory() as tmpdir:
            msp_path = os.path.join(tmpdir, "test.msp")
            traml_path = os.path.join(tmpdir, "test.traml")
            create_synthetic_msp(msp_path, n_spectra=2)
            count = convert_msp_to_traml(msp_path, traml_path)
            assert count == 2
            assert os.path.exists(traml_path)

            # Verify TraML can be loaded back
            te = oms.TargetedExperiment()
            oms.TraMLFile().load(traml_path, te)
            assert len(te.getTransitions()) > 0

    def test_parse_msp_name_extraction(self):
        from spectral_library_format_converter import parse_msp

        with tempfile.TemporaryDirectory() as tmpdir:
            msp_path = os.path.join(tmpdir, "test.msp")
            with open(msp_path, "w") as f:
                f.write("Name: PEPTIDEK/2\n")
                f.write("MW: 500.0\n")
                f.write("Comment: Charge=2+ Parent=251.5\n")
                f.write("Num peaks: 2\n")
                f.write("100.0\t1000.0\n")
                f.write("200.0\t500.0\n")
                f.write("\n")
            spectra = parse_msp(msp_path)
            assert len(spectra) == 1
            assert spectra[0]["name"] == "PEPTIDEK/2"

    def test_create_synthetic_msp(self):
        from spectral_library_format_converter import create_synthetic_msp

        with tempfile.TemporaryDirectory() as tmpdir:
            msp_path = os.path.join(tmpdir, "test.msp")
            create_synthetic_msp(msp_path, n_spectra=5)
            assert os.path.exists(msp_path)
            with open(msp_path) as f:
                content = f.read()
            assert "Name:" in content
            assert "Num peaks:" in content
