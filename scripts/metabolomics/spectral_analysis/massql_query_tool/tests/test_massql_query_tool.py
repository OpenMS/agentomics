"""Tests for massql_query_tool."""

import pytest
from conftest import requires_pyopenms


@requires_pyopenms
class TestMassqlQueryTool:
    def _make_experiment(self):
        import numpy as np
        import pyopenms as oms

        exp = oms.MSExperiment()
        # MS1 with peaks at 180.0, 181.0
        ms1 = oms.MSSpectrum()
        ms1.setMSLevel(1)
        ms1.setRT(60.0)
        ms1.set_peaks([
            np.array([180.0, 181.0, 250.0], dtype=np.float64),
            np.array([1000.0, 200.0, 500.0], dtype=np.float64),
        ])
        exp.addSpectrum(ms1)

        # MS2 with precursor at 500.0 and product at 226.18
        ms2 = oms.MSSpectrum()
        ms2.setMSLevel(2)
        ms2.setRT(61.0)
        prec = oms.Precursor()
        prec.setMZ(500.0)
        ms2.setPrecursors([prec])
        ms2.set_peaks([
            np.array([226.18, 300.0], dtype=np.float64),
            np.array([800.0, 400.0], dtype=np.float64),
        ])
        exp.addSpectrum(ms2)

        return exp

    def test_parse_ms2prod_query(self):
        from massql_query_tool import parse_query

        q = parse_query("MS2PROD=226.18")
        assert q["query_type"] == "MS2PROD"
        assert q["target_mz"] == 226.18

    def test_parse_ms1mz_query(self):
        from massql_query_tool import parse_query

        q = parse_query("MS1MZ=180.06")
        assert q["query_type"] == "MS1MZ"

    def test_parse_precmz_query(self):
        from massql_query_tool import parse_query

        q = parse_query("PRECMZ=500.0")
        assert q["query_type"] == "PRECMZ"

    def test_invalid_query(self):
        from massql_query_tool import parse_query

        with pytest.raises(ValueError):
            parse_query("INVALID=123")

    def test_ms2prod_search(self):
        from massql_query_tool import execute_query, parse_query

        exp = self._make_experiment()
        query = parse_query("MS2PROD=226.18")
        results = execute_query(exp, query, tolerance_da=0.5)
        assert len(results) == 1
        assert abs(results[0]["matched_mz"] - 226.18) < 0.01

    def test_ms1mz_search(self):
        from massql_query_tool import execute_query, parse_query

        exp = self._make_experiment()
        query = parse_query("MS1MZ=180.0")
        results = execute_query(exp, query, tolerance_da=0.5)
        assert len(results) == 1

    def test_precmz_search(self):
        from massql_query_tool import execute_query, parse_query

        exp = self._make_experiment()
        query = parse_query("PRECMZ=500.0")
        results = execute_query(exp, query, tolerance_da=0.5)
        assert len(results) == 1

    def test_no_match(self):
        from massql_query_tool import execute_query, parse_query

        exp = self._make_experiment()
        query = parse_query("MS2PROD=999.0")
        results = execute_query(exp, query, tolerance_da=0.1)
        assert len(results) == 0
