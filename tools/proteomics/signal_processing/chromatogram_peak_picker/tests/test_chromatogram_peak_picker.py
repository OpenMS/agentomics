"""Tests for chromatogram_peak_picker."""

import os
import tempfile

import pytest

pytest.importorskip("pyopenms")


class TestChromatogramPeakPicker:
    def test_pick_chromatogram_peaks_returns_count(self):
        from chromatogram_peak_picker import (
            create_synthetic_chromatogram_mzml,
            pick_chromatogram_peaks,
        )

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "picked.mzML")

            create_synthetic_chromatogram_mzml(input_path)
            count = pick_chromatogram_peaks(input_path, output_path)

            assert count == 2

    def test_picked_chromatograms_have_peaks(self):
        import pyopenms as oms
        from chromatogram_peak_picker import (
            create_synthetic_chromatogram_mzml,
            pick_chromatogram_peaks,
        )

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "picked.mzML")

            create_synthetic_chromatogram_mzml(input_path)
            pick_chromatogram_peaks(input_path, output_path)

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            chroms = exp.getChromatograms()
            assert len(chroms) >= 1

            # Each picked chromatogram should have at least one peak
            for chrom in chroms:
                assert chrom.size() >= 1

    def test_peak_detected_near_expected_rt(self):
        import pyopenms as oms
        from chromatogram_peak_picker import (
            create_synthetic_chromatogram_mzml,
            pick_chromatogram_peaks,
        )

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "input.mzML")
            output_path = os.path.join(tmp, "picked.mzML")

            create_synthetic_chromatogram_mzml(input_path)
            pick_chromatogram_peaks(input_path, output_path)

            exp = oms.MSExperiment()
            oms.MzMLFile().load(output_path, exp)

            chroms = exp.getChromatograms()
            # First chromatogram should have peak near RT 300
            chrom1 = chroms[0]
            rts = [chrom1[i].getRT() for i in range(chrom1.size())]
            intensities = [chrom1[i].getIntensity() for i in range(chrom1.size())]

            if len(rts) > 0 and len(intensities) > 0:
                max_idx = intensities.index(max(intensities))
                peak_rt = rts[max_idx]
                assert abs(peak_rt - 300.0) < 30.0, (
                    f"Expected peak near RT 300, found at {peak_rt}"
                )

    def test_no_chromatograms_returns_zero(self):
        import pyopenms as oms
        from chromatogram_peak_picker import pick_chromatogram_peaks

        with tempfile.TemporaryDirectory() as tmp:
            input_path = os.path.join(tmp, "empty.mzML")
            output_path = os.path.join(tmp, "picked.mzML")

            # Create mzML with no chromatograms
            exp = oms.MSExperiment()
            spectrum = oms.MSSpectrum()
            spectrum.setMSLevel(1)
            exp.addSpectrum(spectrum)
            oms.MzMLFile().store(input_path, exp)

            count = pick_chromatogram_peaks(input_path, output_path)
            assert count == 0
