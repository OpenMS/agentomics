"""Tests for spectra_merger."""

import os
import tempfile

import pytest

oms = pytest.importorskip("pyopenms")


def _create_ms1_experiment(n_spectra=9, n_peaks=5, rt_step=10.0):
    """Create a synthetic MSExperiment with MS1 spectra."""
    exp = oms.MSExperiment()
    for i in range(n_spectra):
        spec = oms.MSSpectrum()
        spec.setMSLevel(1)
        spec.setRT(float(i) * rt_step)
        mzs = [200.0 + j * 100 for j in range(n_peaks)]
        ints = [1000.0] * n_peaks
        spec.set_peaks((mzs, ints))
        exp.addSpectrum(spec)
    return exp


def _create_ms2_experiment_shared_precursors():
    """Create MS2 spectra with shared precursor m/z values."""
    exp = oms.MSExperiment()
    # 3 MS2 spectra with precursor at 500.0
    for i in range(3):
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(float(i) * 10.0)
        prec = oms.Precursor()
        prec.setMZ(500.0)
        prec.setCharge(2)
        spec.setPrecursors([prec])
        mzs = [100.0 + j * 50 for j in range(5)]
        ints = [1000.0] * 5
        spec.set_peaks((mzs, ints))
        exp.addSpectrum(spec)
    # 2 MS2 spectra with precursor at 600.0
    for i in range(2):
        spec = oms.MSSpectrum()
        spec.setMSLevel(2)
        spec.setRT(float(i + 3) * 10.0)
        prec = oms.Precursor()
        prec.setMZ(600.0)
        prec.setCharge(2)
        spec.setPrecursors([prec])
        mzs = [200.0 + j * 50 for j in range(5)]
        ints = [500.0] * 5
        spec.set_peaks((mzs, ints))
        exp.addSpectrum(spec)
    return exp


def _store_mzml(exp, path):
    """Store an MSExperiment to mzML."""
    oms.MzMLFile().store(path, exp)


class TestSpectraMerger:
    def test_block_merge_9_spectra_block_size_3(self):
        """9 MS1 spectra merged in blocks of 3 should produce 3 spectra."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms1_experiment(n_spectra=9)
            _store_mzml(exp, in_path)

            n_after = merge_spectra(in_path, out_path, mode="block", block_size=3)
            assert n_after == 3

            # Verify output file was written and contains 3 spectra
            out_exp = oms.MSExperiment()
            oms.MzMLFile().load(out_path, out_exp)
            assert out_exp.getNrSpectra() == 3

    def test_block_merge_preserves_peaks(self):
        """Block-merged spectra should retain all m/z positions."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms1_experiment(n_spectra=6, n_peaks=5)
            _store_mzml(exp, in_path)

            merge_spectra(in_path, out_path, mode="block", block_size=3)

            out_exp = oms.MSExperiment()
            oms.MzMLFile().load(out_path, out_exp)
            for i in range(out_exp.getNrSpectra()):
                mzs, ints = out_exp.getSpectrum(i).get_peaks()
                assert len(mzs) == 5

    def test_block_merge_single_block(self):
        """When block size >= n_spectra, all spectra merge into one."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms1_experiment(n_spectra=3)
            _store_mzml(exp, in_path)

            n_after = merge_spectra(in_path, out_path, mode="block", block_size=5)
            assert n_after == 1

    def test_precursor_merge_shared_precursors(self):
        """MS2 spectra sharing a precursor m/z should merge together."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms2_experiment_shared_precursors()
            _store_mzml(exp, in_path)

            n_after = merge_spectra(in_path, out_path, mode="precursor")
            # 3 spectra at 500.0 and 2 at 600.0 -> 2 merged spectra
            assert n_after == 2

    def test_precursor_merge_output_has_precursors(self):
        """Merged MS2 spectra should retain precursor information."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms2_experiment_shared_precursors()
            _store_mzml(exp, in_path)

            merge_spectra(in_path, out_path, mode="precursor")

            out_exp = oms.MSExperiment()
            oms.MzMLFile().load(out_path, out_exp)
            precursor_mzs = set()
            for i in range(out_exp.getNrSpectra()):
                precs = out_exp.getSpectrum(i).getPrecursors()
                if precs:
                    precursor_mzs.add(round(precs[0].getMZ(), 1))
            assert 500.0 in precursor_mzs
            assert 600.0 in precursor_mzs

    def test_default_mode_is_block(self):
        """Default mode should be block merge."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms1_experiment(n_spectra=9)
            _store_mzml(exp, in_path)

            n_after = merge_spectra(in_path, out_path)
            # Default block_size=3, 9 spectra -> 3
            assert n_after == 3

    def test_output_file_created(self):
        """Output mzML file should exist after merging."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms1_experiment(n_spectra=6)
            _store_mzml(exp, in_path)

            merge_spectra(in_path, out_path, mode="block", block_size=2)
            assert os.path.exists(out_path)

    def test_ms_level_preserved_after_block_merge(self):
        """Merged spectra should keep MS level 1."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms1_experiment(n_spectra=9)
            _store_mzml(exp, in_path)

            merge_spectra(in_path, out_path, mode="block", block_size=3)

            out_exp = oms.MSExperiment()
            oms.MzMLFile().load(out_path, out_exp)
            for i in range(out_exp.getNrSpectra()):
                assert out_exp.getSpectrum(i).getMSLevel() == 1

    def test_block_merge_uneven_division(self):
        """7 spectra with block size 3 should produce 3 merged spectra."""
        from spectra_merger import merge_spectra

        with tempfile.TemporaryDirectory() as tmpdir:
            in_path = os.path.join(tmpdir, "input.mzML")
            out_path = os.path.join(tmpdir, "output.mzML")
            exp = _create_ms1_experiment(n_spectra=7)
            _store_mzml(exp, in_path)

            n_after = merge_spectra(in_path, out_path, mode="block", block_size=3)
            # 7 / 3 = 2 full blocks + 1 remainder block = 3 output spectra
            assert n_after == 3
