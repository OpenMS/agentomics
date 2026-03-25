"""Tests for search_result_merger."""

import pytest
from conftest import requires_pyopenms
from search_result_merger import _make_key, merge_results, read_identification_tsv


@requires_pyopenms
class TestSearchResultMerger:
    def _write_tsv(self, tmp_path, name, rows):
        filepath = str(tmp_path / name)
        with open(filepath, "w") as fh:
            if rows:
                keys = list(rows[0].keys())
                fh.write("\t".join(keys) + "\n")
                for row in rows:
                    fh.write("\t".join(str(row[k]) for k in keys) + "\n")
        return filepath

    def test_union_merge(self, tmp_path):
        f1 = self._write_tsv(tmp_path, "e1.tsv", [
            {"peptide": "PEPTIDEK", "charge": "2", "score": "0.99"},
            {"peptide": "TESTPEP", "charge": "2", "score": "0.95"},
        ])
        f2 = self._write_tsv(tmp_path, "e2.tsv", [
            {"peptide": "PEPTIDEK", "charge": "2", "score": "0.98"},
            {"peptide": "ANOTHERPEP", "charge": "3", "score": "0.90"},
        ])
        _, merged = merge_results([f1, f2], method="union")
        peptides = [r["peptide"] for r in merged]
        assert "PEPTIDEK" in peptides
        assert "TESTPEP" in peptides
        assert "ANOTHERPEP" in peptides
        assert len(merged) == 3

    def test_intersection_merge(self, tmp_path):
        f1 = self._write_tsv(tmp_path, "e1.tsv", [
            {"peptide": "PEPTIDEK", "charge": "2"},
            {"peptide": "TESTPEP", "charge": "2"},
        ])
        f2 = self._write_tsv(tmp_path, "e2.tsv", [
            {"peptide": "PEPTIDEK", "charge": "2"},
            {"peptide": "ANOTHERPEP", "charge": "3"},
        ])
        _, merged = merge_results([f1, f2], method="intersection")
        assert len(merged) == 1
        assert merged[0]["peptide"] == "PEPTIDEK"

    def test_n_engines_count(self, tmp_path):
        f1 = self._write_tsv(tmp_path, "e1.tsv", [{"peptide": "PEP", "charge": "2"}])
        f2 = self._write_tsv(tmp_path, "e2.tsv", [{"peptide": "PEP", "charge": "2"}])
        _, merged = merge_results([f1, f2], method="union")
        assert merged[0]["n_engines"] == "2"

    def test_unknown_method(self, tmp_path):
        f1 = self._write_tsv(tmp_path, "e1.tsv", [{"peptide": "PEP"}])
        with pytest.raises(ValueError, match="Unknown method"):
            merge_results([f1], method="invalid")

    def test_make_key(self):
        row = {"peptide": "PEPTIDEK", "charge": "2", "spectrum": "scan1"}
        key = _make_key(row)
        assert "PEPTIDEK" in key
        assert "2" in key
        assert "scan1" in key

    def test_empty_input(self, tmp_path):
        f1 = self._write_tsv(tmp_path, "e1.tsv", [])
        _, merged = merge_results([f1], method="union")
        assert len(merged) == 0

    def test_read_tsv(self, tmp_path):
        filepath = self._write_tsv(tmp_path, "test.tsv", [
            {"peptide": "PEPTIDEK", "score": "0.99"},
        ])
        rows = read_identification_tsv(filepath)
        assert len(rows) == 1
        assert rows[0]["peptide"] == "PEPTIDEK"
