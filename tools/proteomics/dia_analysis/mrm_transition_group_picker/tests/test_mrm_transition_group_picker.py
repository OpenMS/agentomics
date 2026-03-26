"""Tests for mrm_transition_group_picker."""

import os

import pytest

pyopenms = pytest.importorskip("pyopenms")
import numpy as np  # noqa: E402


def _make_group_chromatograms(path, group_id="group_0", n_transitions=3, peak_rt=50.0):
    """Create synthetic chromatograms for a transition group with Gaussian peaks."""
    exp = pyopenms.MSExperiment()
    times = np.linspace(0, 100, 100)

    for i in range(n_transitions):
        chrom = pyopenms.MSChromatogram()
        chrom.setNativeID(f"{group_id}_{i}".encode())
        # Same peak RT, slightly different intensities
        ints = np.array(
            [(1000.0 - i * 100) * np.exp(-0.5 * ((t - peak_rt) / 5.0) ** 2)
             for t in times]
        )
        chrom.set_peaks((times.tolist(), ints.tolist()))
        exp.addChromatogram(chrom)

    pyopenms.MzMLFile().store(path, exp)


class TestMRMTransitionGroupPicker:
    def test_pick_single_group(self, tmp_path):
        from mrm_transition_group_picker import pick_transition_groups

        chrom_path = str(tmp_path / "chromatograms.mzML")
        out_path = str(tmp_path / "picked.mzML")

        _make_group_chromatograms(chrom_path, "group_0", 2, 50.0)

        n = pick_transition_groups(chrom_path, out_path)
        assert n >= 1
        assert os.path.exists(out_path)

    def test_output_has_chromatograms(self, tmp_path):
        from mrm_transition_group_picker import pick_transition_groups

        chrom_path = str(tmp_path / "chromatograms.mzML")
        out_path = str(tmp_path / "picked.mzML")

        _make_group_chromatograms(chrom_path, "group_0", 3, 50.0)
        pick_transition_groups(chrom_path, out_path)

        out_exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(out_path, out_exp)
        assert out_exp.getNrChromatograms() >= 1

    def test_multiple_groups(self, tmp_path):
        from mrm_transition_group_picker import pick_transition_groups

        exp = pyopenms.MSExperiment()
        times = np.linspace(0, 100, 100)

        for gid in ["grpA", "grpB"]:
            for i in range(2):
                chrom = pyopenms.MSChromatogram()
                chrom.setNativeID(f"{gid}_{i}".encode())
                peak_rt = 40.0 if gid == "grpA" else 60.0
                ints = np.array(
                    [1000.0 * np.exp(-0.5 * ((t - peak_rt) / 5.0) ** 2)
                     for t in times]
                )
                chrom.set_peaks((times.tolist(), ints.tolist()))
                exp.addChromatogram(chrom)

        chrom_path = str(tmp_path / "chromatograms.mzML")
        out_path = str(tmp_path / "picked.mzML")
        pyopenms.MzMLFile().store(chrom_path, exp)

        n = pick_transition_groups(chrom_path, out_path)
        # Should find at least 1 feature per group
        assert n >= 2

    def test_returns_int(self, tmp_path):
        from mrm_transition_group_picker import pick_transition_groups

        chrom_path = str(tmp_path / "chromatograms.mzML")
        out_path = str(tmp_path / "picked.mzML")

        _make_group_chromatograms(chrom_path, "group_0", 2, 50.0)
        result = pick_transition_groups(chrom_path, out_path)
        assert isinstance(result, int)

    def test_parse_group_id(self):
        from mrm_transition_group_picker import _parse_group_id

        assert _parse_group_id("group_0_1") == "group_0"
        assert _parse_group_id("group_0_2") == "group_0"
        assert _parse_group_id("single") == "single"
        assert _parse_group_id("complex_name_3") == "complex_name"
