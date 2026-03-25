"""Tests for collision_energy_analyzer."""

import numpy as np
import pytest

pytest.importorskip("pyopenms")


class TestCollisionEnergyAnalyzer:
    def _make_experiment(self, n_ms2=5, ce_value=30.0):
        """Create a synthetic MSExperiment with MS2 spectra and CE values."""
        import pyopenms as oms

        exp = oms.MSExperiment()
        # Add an MS1 spectrum
        ms1 = oms.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(0.0)
        mzs = np.array([400.0, 500.0, 600.0], dtype=np.float64)
        ints = np.array([1000.0, 2000.0, 1500.0], dtype=np.float64)
        ms1.set_peaks([mzs, ints])
        exp.addSpectrum(ms1)

        for i in range(n_ms2):
            spec = oms.MSSpectrum()
            spec.setMSLevel(2)
            spec.setRT(10.0 * (i + 1))
            mzs = np.array([200.0, 300.0, 400.0], dtype=np.float64)
            ints = np.array([500.0, 800.0, 300.0], dtype=np.float64)
            spec.set_peaks([mzs, ints])

            prec = oms.Precursor()
            prec.setMZ(500.0 + i * 10)
            prec.setCharge(2)
            prec.setActivationEnergy(ce_value + i * 5)
            spec.setPrecursors([prec])
            exp.addSpectrum(spec)

        return exp

    def test_extract_ce(self):
        from collision_energy_analyzer import extract_collision_energies

        exp = self._make_experiment(n_ms2=3, ce_value=25.0)
        records = extract_collision_energies(exp)
        assert len(records) == 3
        assert all(r["collision_energy"] != "N/A" for r in records)

    def test_extract_ce_values(self):
        from collision_energy_analyzer import extract_collision_energies

        exp = self._make_experiment(n_ms2=3, ce_value=25.0)
        records = extract_collision_energies(exp)
        assert records[0]["collision_energy"] == 25.0
        assert records[1]["collision_energy"] == 30.0

    def test_summarize(self):
        from collision_energy_analyzer import extract_collision_energies, summarize_ce

        exp = self._make_experiment(n_ms2=4, ce_value=20.0)
        records = extract_collision_energies(exp)
        summary = summarize_ce(records)
        assert summary["total_ms2"] == 4
        assert summary["with_ce"] == 4
        assert summary["min_ce"] == 20.0

    def test_empty_experiment(self):
        import pyopenms as oms
        from collision_energy_analyzer import extract_collision_energies, summarize_ce

        exp = oms.MSExperiment()
        records = extract_collision_energies(exp)
        assert records == []
        summary = summarize_ce(records)
        assert summary["total_ms2"] == 0

    def test_write_tsv(self, tmp_path):
        from collision_energy_analyzer import extract_collision_energies, write_tsv

        exp = self._make_experiment(n_ms2=2)
        records = extract_collision_energies(exp)
        out = str(tmp_path / "ce.tsv")
        write_tsv(records, out)
        with open(out) as fh:
            lines = fh.readlines()
        assert len(lines) == 3  # header + 2 rows
