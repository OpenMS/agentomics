"""Tests for volcano_plot_data_generator."""

import math

from conftest import requires_pyopenms
from volcano_plot_data_generator import generate_volcano_data, read_de_results, summarize_volcano


@requires_pyopenms
class TestVolcanoPlotDataGenerator:
    def _make_de_results(self):
        return [
            {"feature": "prot1", "log2fc": 2.0, "pvalue": 0.001},    # up
            {"feature": "prot2", "log2fc": -1.5, "pvalue": 0.01},    # down
            {"feature": "prot3", "log2fc": 0.5, "pvalue": 0.001},    # ns (fc too low)
            {"feature": "prot4", "log2fc": 2.0, "pvalue": 0.1},      # ns (pval too high)
            {"feature": "prot5", "log2fc": float("nan"), "pvalue": float("nan")},  # ns
        ]

    def test_classification(self):
        results = self._make_de_results()
        volcano = generate_volcano_data(results, fc_threshold=1.0, pvalue_threshold=0.05)
        regs = {v["feature"]: v["regulation"] for v in volcano}
        assert regs["prot1"] == "up"
        assert regs["prot2"] == "down"
        assert regs["prot3"] == "ns"
        assert regs["prot4"] == "ns"
        assert regs["prot5"] == "ns"

    def test_neg_log10_pvalue(self):
        results = [{"feature": "p1", "log2fc": 1.0, "pvalue": 0.01}]
        volcano = generate_volcano_data(results)
        assert abs(volcano[0]["neg_log10_pvalue"] - 2.0) < 0.01

    def test_summarize(self):
        results = self._make_de_results()
        volcano = generate_volcano_data(results, fc_threshold=1.0, pvalue_threshold=0.05)
        counts = summarize_volcano(volcano)
        assert counts["up"] == 1
        assert counts["down"] == 1
        assert counts["ns"] == 3

    def test_custom_thresholds(self):
        results = self._make_de_results()
        volcano = generate_volcano_data(results, fc_threshold=0.3, pvalue_threshold=0.5)
        regs = {v["feature"]: v["regulation"] for v in volcano}
        assert regs["prot3"] == "up"  # fc=0.5 > 0.3 threshold
        assert regs["prot4"] == "up"  # pval=0.1 < 0.5 threshold

    def test_read_de_results(self, tmp_path):
        infile = str(tmp_path / "de.tsv")
        with open(infile, "w") as fh:
            fh.write("feature\tlog2fc\tadj_pvalue\n")
            fh.write("p1\t1.5\t0.01\n")
            fh.write("p2\tNA\tNA\n")
        results = read_de_results(infile)
        assert len(results) == 2
        assert results[0]["log2fc"] == 1.5
        assert math.isnan(results[1]["log2fc"])

    def test_empty_input(self):
        volcano = generate_volcano_data([])
        assert volcano == []
