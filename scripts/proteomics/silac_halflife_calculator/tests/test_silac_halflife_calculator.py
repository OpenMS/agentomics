"""Tests for silac_halflife_calculator."""

import csv
import math
import sys

from conftest import requires_pyopenms


@requires_pyopenms
def test_exponential_decay():
    import numpy as np
    from silac_halflife_calculator import exponential_decay

    t = np.array([0, 1, 2, 3])
    result = exponential_decay(t, r0=10.0, k=0.5)
    assert abs(result[0] - 10.0) < 0.01
    assert result[1] < result[0]  # decaying


@requires_pyopenms
def test_fit_halflife_perfect():
    from silac_halflife_calculator import fit_halflife

    # Generate perfect exponential data
    k_true = 0.05
    r0_true = 10.0
    t = [0, 6, 12, 24, 48]
    ratios = [r0_true * math.exp(-k_true * ti) for ti in t]

    result = fit_halflife(t, ratios)
    assert result is not None
    assert abs(result["k"] - k_true) < 0.01
    expected_halflife = math.log(2) / k_true
    assert abs(result["halflife"] - expected_halflife) < 1.0
    assert result["r_squared"] > 0.99


@requires_pyopenms
def test_fit_halflife_insufficient_data():
    from silac_halflife_calculator import fit_halflife

    result = fit_halflife([0], [10.0])
    assert result is None


@requires_pyopenms
def test_compute_halflives():
    import math

    from silac_halflife_calculator import compute_halflives

    k = 0.05
    timepoints = [0, 6, 12, 24, 48]
    proteins = {
        "P1": [10.0 * math.exp(-k * t) for t in timepoints],
        "P2": [5.0 * math.exp(-0.1 * t) for t in timepoints],
    }

    results = compute_halflives(proteins, timepoints)
    assert len(results) == 2
    assert results[0]["status"] == "ok"
    assert results[1]["status"] == "ok"
    # P2 has higher k -> shorter halflife
    assert results[1]["halflife"] < results[0]["halflife"]


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    import math

    from silac_halflife_calculator import main

    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "output.tsv"
    timepoints = [0, 6, 12, 24, 48]
    k = 0.05

    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        cols = ["protein_id"] + [f"t{t}" for t in timepoints]
        writer.writerow(cols)
        ratios = [10.0 * math.exp(-k * t) for t in timepoints]
        writer.writerow(["P1"] + [f"{r:.4f}" for r in ratios])

    sys.argv = [
        "silac_halflife_calculator.py",
        "--input", str(input_file),
        "--timepoints", ",".join(str(t) for t in timepoints),
        "--output", str(output_file),
    ]
    main()

    assert output_file.exists()
    with open(output_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert len(rows) == 1
    assert rows[0]["status"] == "ok"
