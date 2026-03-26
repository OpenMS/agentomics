"""Tests for isobaric_quantifier."""

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
        "tmt16plex": oms.TMTSixteenPlexQuantitationMethod,
        "itraq4plex": oms.ItraqFourPlexQuantitationMethod,
        "itraq8plex": oms.ItraqEightPlexQuantitationMethod,
    }
    qm = method_map[method]()
    channels = qm.getChannelInformation()

    # Build an mzML with reporter ions and extract channels
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
        int_vals = [1000.0 * (i + 1) + scan_idx * 500 for i in range(len(channels))]
        spec.set_peaks([mz_vals, int_vals])
        spec.sortByPosition()
        exp.addSpectrum(spec)

    ext = oms.IsobaricChannelExtractor(qm)
    consensus = oms.ConsensusMap()
    ext.extractChannels(exp, consensus)

    oms.ConsensusXMLFile().store(path, consensus)
    return consensus.size()


class TestQuantifyIsobaric:
    def test_quantification_produces_output(self):
        """Quantifier runs and produces an output file with features."""
        from isobaric_quantifier import quantify_isobaric

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "quantified.consensusXML")

            n_input = _create_synthetic_consensus(in_path, method="tmt6plex")
            n_features = quantify_isobaric(in_path, "tmt6plex", out_path)

            assert os.path.exists(out_path)
            assert n_features == n_input

    def test_quantified_intensities_differ_with_correction(self):
        """With isotope correction, quantified intensities should differ from raw."""
        import pyopenms as oms
        from isobaric_quantifier import quantify_isobaric

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "quantified.consensusXML")

            _create_synthetic_consensus(in_path, method="tmt6plex")
            quantify_isobaric(in_path, "tmt6plex", out_path)

            # Load input
            cmap_in = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(in_path, cmap_in)

            # Load output
            cmap_out = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cmap_out)

            # Both should have same number of features
            assert cmap_in.size() == cmap_out.size()

            # Compare intensities: correction should change at least some values
            any_different = False
            for cf_in, cf_out in zip(cmap_in, cmap_out):
                handles_in = sorted(cf_in.getFeatureList(), key=lambda h: h.getMZ())
                handles_out = sorted(cf_out.getFeatureList(), key=lambda h: h.getMZ())
                for h_in, h_out in zip(handles_in, handles_out):
                    if abs(h_in.getIntensity() - h_out.getIntensity()) > 0.1:
                        any_different = True
                        break
                if any_different:
                    break

            assert any_different, "Isotope correction should change intensities"

    def test_channel_count_preserved(self):
        """Each output feature should have the same number of channel handles."""
        import pyopenms as oms
        from isobaric_quantifier import quantify_isobaric

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "quantified.consensusXML")

            _create_synthetic_consensus(in_path, method="tmt6plex")
            quantify_isobaric(in_path, "tmt6plex", out_path)

            cmap = oms.ConsensusMap()
            oms.ConsensusXMLFile().load(out_path, cmap)

            for cf in cmap:
                assert cf.size() == 6

    def test_itraq4plex_quantification(self):
        """Quantifier works for iTRAQ 4-plex."""
        from isobaric_quantifier import quantify_isobaric

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "quantified.consensusXML")

            _create_synthetic_consensus(in_path, method="itraq4plex")
            n_features = quantify_isobaric(in_path, "itraq4plex", out_path)

            assert os.path.exists(out_path)
            assert n_features >= 1

    def test_multiple_features_quantified(self):
        """Quantifier processes all features from a multi-spectrum input."""
        from isobaric_quantifier import quantify_isobaric

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "quantified.consensusXML")

            n_input = _create_synthetic_consensus(in_path, method="tmt6plex")
            n_output = quantify_isobaric(in_path, "tmt6plex", out_path)

            assert n_input == 3
            assert n_output == 3

    def test_invalid_method_raises(self):
        """Raise ValueError for unsupported quantitation method."""
        from isobaric_quantifier import quantify_isobaric

        with tempfile.TemporaryDirectory() as tmp:
            in_path = os.path.join(tmp, "quant.consensusXML")
            out_path = os.path.join(tmp, "quantified.consensusXML")

            _create_synthetic_consensus(in_path, method="tmt6plex")

            with pytest.raises(ValueError, match="Unsupported"):
                quantify_isobaric(in_path, "tmt99plex", out_path)
