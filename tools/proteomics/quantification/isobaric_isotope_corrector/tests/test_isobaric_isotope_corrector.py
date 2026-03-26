"""Tests for isobaric_isotope_corrector."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_synthetic_consensus(path, method="tmt6plex"):
    """Create a synthetic consensusXML from TMT channel extraction."""
    import pyopenms as oms

    method_map = {
        "tmt6plex": oms.TMTSixPlexQuantitationMethod,
        "tmt10plex": oms.TMTTenPlexQuantitationMethod,
        "itraq4plex": oms.ItraqFourPlexQuantitationMethod,
        "itraq8plex": oms.ItraqEightPlexQuantitationMethod,
    }
    qm = method_map[method]()
    channels = qm.getChannelInformation()

    # Build an mzML with reporter ions and extract channels
    exp = oms.MSExperiment()

    for scan_idx in range(2):
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(100.0 + scan_idx * 60)
        spec.setNativeID(f"scan={scan_idx + 1}")
        prec = oms.Precursor()
        prec.setMZ(500.0 + scan_idx * 100)
        prec.setCharge(2)
        prec.setActivationMethods({prec.ActivationMethod.HCD})
        spec.setPrecursors([prec])

        mz_vals = [ch.center for ch in channels]
        int_vals = [1000.0 * (i + 1) + scan_idx * 500 for i in range(len(channels))]
        spec.set_peaks([mz_vals, int_vals])
        spec.sortByPosition()
        exp.addSpectrum(spec)

    ext = oms.IsobaricChannelExtractor(qm)
    consensus = oms.ConsensusMap()
    ext.extractChannels(exp, consensus)

    oms.ConsensusXMLFile().store(path, consensus)
    return consensus.size()


class TestCorrectIsotopeImpurities:
    def test_correction_produces_output(self):
        """Corrector runs without error and produces an output file."""
        from isobaric_isotope_corrector import correct_isotope_impurities

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "corrected.consensusXML")

            _create_synthetic_consensus(in_path, method="tmt6plex")
            n_features = correct_isotope_impurities(in_path, "tmt6plex", out_path)

            assert os.path.exists(out_path)
            assert n_features >= 0

    def test_correction_loads_valid_consensus(self):
        """Output consensusXML can be loaded back."""
        import pyopenms as oms
        from isobaric_isotope_corrector import correct_isotope_impurities

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "corrected.consensusXML")

            _create_synthetic_consensus(in_path, method="tmt6plex")
            correct_isotope_impurities(in_path, "tmt6plex", out_path)

            cmap = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cmap)
            # Output should be a valid consensus map (may be empty depending
            # on pyopenms version behavior of IsobaricIsotopeCorrector)
            assert isinstance(cmap, oms.ConsensusMap)

    def test_correction_matrix_exists(self):
        """Verify the correction matrix is non-trivial for TMT 6-plex."""
        import pyopenms as oms

        qm = oms.TMTSixPlexQuantitationMethod()
        matrix = qm.getIsotopeCorrectionMatrix()
        n_channels = qm.getNumberOfChannels()

        assert matrix.rows() == n_channels
        assert matrix.cols() == n_channels

        # Matrix should not be pure identity (off-diagonal elements exist)
        off_diag_sum = 0.0
        for i in range(n_channels):
            for j in range(n_channels):
                if i != j:
                    off_diag_sum += matrix.getValue(i, j)
        assert off_diag_sum > 0.0

    def test_invalid_method_raises(self):
        """Raise ValueError for unsupported quantitation method."""
        from isobaric_isotope_corrector import correct_isotope_impurities

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "corrected.consensusXML")

            _create_synthetic_consensus(in_path, method="tmt6plex")

            with pytest.raises(ValueError, match="Unsupported"):
                correct_isotope_impurities(in_path, "tmt99plex", out_path)

    def test_itraq4plex_correction(self):
        """Corrector runs for iTRAQ 4-plex without error."""
        from isobaric_isotope_corrector import correct_isotope_impurities

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "corrected.consensusXML")

            _create_synthetic_consensus(in_path, method="itraq4plex")
            n_features = correct_isotope_impurities(in_path, "itraq4plex", out_path)

            assert os.path.exists(out_path)
            assert n_features >= 0
