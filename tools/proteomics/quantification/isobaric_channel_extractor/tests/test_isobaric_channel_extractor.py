"""Tests for isobaric_channel_extractor."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


def _create_synthetic_mzml(path, method="tmt6plex"):
    """Create a synthetic mzML with reporter ion peaks for the given method."""
    import pyopenms as oms

    method_map = {
        "tmt6plex": oms.TMTSixPlexQuantitationMethod,
        "tmt10plex": oms.TMTTenPlexQuantitationMethod,
        "tmt16plex": oms.TMTSixteenPlexQuantitationMethod,
        "itraq4plex": oms.ItraqFourPlexQuantitationMethod,
        "itraq8plex": oms.ItraqEightPlexQuantitationMethod,
    }
    qm = method_map[method]()
    channels = qm.getChannelInformation()

    exp = oms.MSExperiment()

    # Create MS2 spectrum with reporter ion peaks
    spec = oms.MSSpectrum()
    spec.setMSLevel(2)
    spec.setRT(120.0)
    spec.setNativeID("scan=1")

    prec = oms.Precursor()
    prec.setMZ(500.0)
    prec.setCharge(2)
    prec.setActivationMethods({prec.ActivationMethod.HCD})
    spec.setPrecursors([prec])

    mz_vals = [ch.center for ch in channels]
    int_vals = [1000.0 * (i + 1) for i in range(len(channels))]
    # Add a few extra fragment peaks
    mz_vals += [200.0, 350.0, 450.0]
    int_vals += [500.0, 700.0, 300.0]
    spec.set_peaks([mz_vals, int_vals])
    spec.sortByPosition()
    exp.addSpectrum(spec)

    oms.MzMLFile().store(path, exp)


class TestExtractChannels:
    def test_tmt6plex_extraction(self):
        """Extract TMT 6-plex reporter channels from a synthetic mzML."""
        from isobaric_channel_extractor import extract_channels

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "run.mzML")
            out_path = os.path.join(tmp, "quant.consensusXML")

            _create_synthetic_mzml(mzml_path, method="tmt6plex")
            n_features = extract_channels(mzml_path, "tmt6plex", out_path)

            assert os.path.exists(out_path)
            assert n_features >= 1

    def test_tmt6plex_channel_count(self):
        """Verify extracted consensus features contain 6 channel handles."""
        import pyopenms as oms
        from isobaric_channel_extractor import extract_channels

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "run.mzML")
            out_path = os.path.join(tmp, "quant.consensusXML")

            _create_synthetic_mzml(mzml_path, method="tmt6plex")
            extract_channels(mzml_path, "tmt6plex", out_path)

            cmap = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cmap)

            for cf in cmap:
                assert cf.size() == 6

    def test_tmt6plex_reporter_mz(self):
        """Verify reporter ion m/z values match TMT 6-plex channels."""
        import pyopenms as oms
        from isobaric_channel_extractor import extract_channels

        qm = oms.TMTSixPlexQuantitationMethod()
        expected_mz = [ch.center for ch in qm.getChannelInformation()]

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "run.mzML")
            out_path = os.path.join(tmp, "quant.consensusXML")

            _create_synthetic_mzml(mzml_path, method="tmt6plex")
            extract_channels(mzml_path, "tmt6plex", out_path)

            cmap = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cmap)

            for cf in cmap:
                handles = cf.getFeatureList()
                observed_mz = sorted(h.getMZ() for h in handles)
                for obs, exp_val in zip(observed_mz, sorted(expected_mz)):
                    assert abs(obs - exp_val) < 0.01

    def test_tmt6plex_intensities(self):
        """Verify extracted intensities correspond to input spectrum peaks."""
        import pyopenms as oms
        from isobaric_channel_extractor import extract_channels

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "run.mzML")
            out_path = os.path.join(tmp, "quant.consensusXML")

            _create_synthetic_mzml(mzml_path, method="tmt6plex")
            extract_channels(mzml_path, "tmt6plex", out_path)

            cmap = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cmap)

            for cf in cmap:
                handles = sorted(cf.getFeatureList(), key=lambda h: h.getMZ())
                # Expected: 1000, 2000, 3000, 4000, 5000, 6000
                for i, h in enumerate(handles):
                    assert h.getIntensity() == pytest.approx(
                        1000.0 * (i + 1), rel=0.01
                    )

    def test_itraq4plex_extraction(self):
        """Extract iTRAQ 4-plex reporter channels."""
        from isobaric_channel_extractor import extract_channels

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "run.mzML")
            out_path = os.path.join(tmp, "quant.consensusXML")

            _create_synthetic_mzml(mzml_path, method="itraq4plex")
            n_features = extract_channels(mzml_path, "itraq4plex", out_path)

            assert os.path.exists(out_path)
            assert n_features >= 1

    def test_multiple_spectra(self):
        """Extract channels from mzML with multiple MS2 spectra."""
        import pyopenms as oms
        from isobaric_channel_extractor import extract_channels

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "run.mzML")
            out_path = os.path.join(tmp, "quant.consensusXML")

            # Build mzML with 3 MS2 spectra
            qm = oms.TMTSixPlexQuantitationMethod()
            channels = qm.getChannelInformation()
            exp = oms.MSExperiment()

            for scan_idx in range(3):
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
                int_vals = [1000.0 * (i + 1) for i in range(len(channels))]
                spec.set_peaks([mz_vals, int_vals])
                spec.sortByPosition()
                exp.addSpectrum(spec)

            oms.MzMLFile().store(mzml_path, exp)

            n_features = extract_channels(mzml_path, "tmt6plex", out_path)
            assert n_features == 3

    def test_invalid_method_raises(self):
        """Raise ValueError for unsupported quantitation method."""
        from isobaric_channel_extractor import extract_channels

        with tempfile.TemporaryDirectory() as tmp:
            mzml_path = os.path.join(tmp, "run.mzML")
            out_path = os.path.join(tmp, "quant.consensusXML")

            _create_synthetic_mzml(mzml_path, method="tmt6plex")

            with pytest.raises(ValueError, match="Unsupported"):
                extract_channels(mzml_path, "tmt99plex", out_path)
