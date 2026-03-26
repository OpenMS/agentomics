"""Tests for internal_calibrator."""

import os
import tempfile

import pytest

oms = pytest.importorskip("pyopenms")

PROTON = 1.007276


def _create_calibration_data(ppm_error=5.0, n_spectra=5):
    """Create synthetic mzML and idXML with systematic m/z error.

    Returns (mzml_path, idxml_path, tmpdir, reference_mzs) where reference_mzs
    are the true theoretical m/z values used.
    """
    peptides = ["PEPTIDEK", "EDITGPR", "ACDEK"]
    ref_data = []
    for seq_str in peptides:
        aa = oms.AASequence.fromString(seq_str)
        theo_mass = aa.getMonoWeight()
        charge = 2
        theo_mz = (theo_mass + charge * PROTON) / charge
        ref_data.append((seq_str, aa, charge, theo_mz))

    all_theo_mzs = sorted([d[3] for d in ref_data])

    # Create mzML with systematic ppm error
    exp = oms.MSExperiment()
    for i in range(n_spectra):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(float(i) * 10.0)
        obs_mzs = [mz * (1.0 + ppm_error * 1e-6) for mz in all_theo_mzs]
        spec.set_peaks((obs_mzs, [1000.0] * len(obs_mzs)))
        exp.addSpectrum(spec)

    # Create peptide IDs
    pep_ids = oms.PeptideIdentificationList()
    prot_id = oms.ProteinIdentification()
    prot_id.setIdentifier("test")

    for i in range(n_spectra):
        for seq_str, aa, charge, theo_mz in ref_data:
            pep_id = oms.PeptideIdentification()
            pep_id.setIdentifier("test")
            pep_id.setRT(float(i) * 10.0)
            obs_mz = theo_mz * (1.0 + ppm_error * 1e-6)
            pep_id.setMZ(obs_mz)
            hit = oms.PeptideHit()
            hit.setSequence(aa)
            hit.setCharge(charge)
            hit.setScore(0.99)
            pep_id.setHits([hit])
            pep_id.setScoreType("XTandem")
            pep_id.setHigherScoreBetter(True)
            pep_ids.push_back(pep_id)

    tmpdir = tempfile.mkdtemp()
    mzml_path = os.path.join(tmpdir, "run.mzML")
    idxml_path = os.path.join(tmpdir, "peptides.idXML")
    oms.MzMLFile().store(mzml_path, exp)
    oms.IdXMLFile().store(idxml_path, [prot_id], pep_ids)

    return mzml_path, idxml_path, tmpdir, all_theo_mzs


class TestInternalCalibrator:
    def test_calibrate_reduces_error(self):
        """Calibration should reduce the systematic ppm error."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, ref_mzs = _create_calibration_data(
            ppm_error=5.0
        )
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        result = calibrate_mz(mzml_path, idxml_path, out_path, model="linear")

        assert result["success"] is True
        assert abs(result["after_ppm_median"]) < abs(result["before_ppm_median"])

    def test_calibrated_output_has_low_error(self):
        """After calibration, ppm error on peaks should be near zero."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, ref_mzs = _create_calibration_data(
            ppm_error=5.0
        )
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        calibrate_mz(mzml_path, idxml_path, out_path, model="linear")

        # Load calibrated data and check errors
        cal_exp = oms.MSExperiment()
        oms.MzMLFile().load(out_path, cal_exp)

        for i in range(cal_exp.getNrSpectra()):
            mzs, _ = cal_exp.getSpectrum(i).get_peaks()
            for j, mz in enumerate(mzs):
                error_ppm = (mz - ref_mzs[j]) / ref_mzs[j] * 1e6
                assert abs(error_ppm) < 1.0, (
                    f"Spectrum {i}, peak {j}: error {error_ppm:.4f} ppm"
                )

    def test_output_file_created(self):
        """Calibrated mzML file should be created."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, _ = _create_calibration_data(ppm_error=5.0)
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        calibrate_mz(mzml_path, idxml_path, out_path, model="linear")
        assert os.path.exists(out_path)

    def test_result_contains_stats(self):
        """Result dict should contain before/after ppm stats."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, _ = _create_calibration_data(ppm_error=5.0)
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        result = calibrate_mz(mzml_path, idxml_path, out_path, model="linear")

        assert "before_ppm_median" in result
        assert "before_ppm_mad" in result
        assert "after_ppm_median" in result
        assert "after_ppm_mad" in result
        assert "n_calibrants" in result
        assert "success" in result

    def test_before_ppm_reflects_injected_error(self):
        """Before-calibration median ppm should match the injected error."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, _ = _create_calibration_data(ppm_error=5.0)
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        result = calibrate_mz(mzml_path, idxml_path, out_path, model="linear")

        # The injected error was +5 ppm; before_ppm_median should be close
        assert abs(result["before_ppm_median"] - 5.0) < 1.0

    def test_quadratic_model(self):
        """Quadratic model should also successfully calibrate."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, _ = _create_calibration_data(ppm_error=5.0)
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        result = calibrate_mz(mzml_path, idxml_path, out_path, model="quadratic")

        assert result["success"] is True
        assert abs(result["after_ppm_median"]) < abs(result["before_ppm_median"])

    def test_spectra_count_preserved(self):
        """Number of spectra should not change after calibration."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, _ = _create_calibration_data(
            ppm_error=5.0, n_spectra=7
        )
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        calibrate_mz(mzml_path, idxml_path, out_path, model="linear")

        orig_exp = oms.MSExperiment()
        oms.MzMLFile().load(mzml_path, orig_exp)
        cal_exp = oms.MSExperiment()
        oms.MzMLFile().load(out_path, cal_exp)

        assert cal_exp.getNrSpectra() == orig_exp.getNrSpectra()

    def test_calibrants_found(self):
        """Calibration should find calibrants from the peptide IDs."""
        from internal_calibrator import calibrate_mz

        mzml_path, idxml_path, tmpdir, _ = _create_calibration_data(ppm_error=5.0)
        out_path = os.path.join(tmpdir, "calibrated.mzML")

        result = calibrate_mz(mzml_path, idxml_path, out_path, model="linear")

        assert result["n_calibrants"] > 0
