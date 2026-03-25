"""Tests for isobaric_purity_corrector."""

import csv
import sys

import numpy as np
from conftest import requires_pyopenms


@requires_pyopenms
def test_get_channels():
    from isobaric_purity_corrector import get_channels

    channels = get_channels("TMT6plex")
    assert len(channels) == 6
    assert "126" in channels

    channels_16 = get_channels("TMT16plex")
    assert len(channels_16) == 16


@requires_pyopenms
def test_get_channels_unknown():
    import pytest
    from isobaric_purity_corrector import get_channels

    with pytest.raises(ValueError):
        get_channels("UnknownLabel")


@requires_pyopenms
def test_correct_intensities_identity():
    from isobaric_purity_corrector import correct_intensities

    # Identity matrix = no correction
    purity = np.eye(3)
    observed = np.array([[100.0, 200.0, 300.0]])
    corrected = correct_intensities(observed, purity)
    np.testing.assert_allclose(corrected, observed, atol=1e-6)


@requires_pyopenms
def test_correct_intensities_with_crosstalk():
    from isobaric_purity_corrector import correct_intensities

    # Purity matrix with some crosstalk
    purity = np.array([
        [0.95, 0.03, 0.02],
        [0.02, 0.94, 0.04],
        [0.01, 0.03, 0.96],
    ])
    # True intensities
    true = np.array([[100.0, 200.0, 300.0]])
    # Observed = purity @ true (simulate measurement)
    observed = true @ purity.T
    # Correct should recover true
    corrected = correct_intensities(observed, purity)
    np.testing.assert_allclose(corrected, true, atol=1.0)


@requires_pyopenms
def test_load_purity_matrix(tmp_path):
    from isobaric_purity_corrector import load_purity_matrix

    purity_file = tmp_path / "purity.csv"
    with open(purity_file, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([0.95, 0.03, 0.02])
        writer.writerow([0.02, 0.94, 0.04])
        writer.writerow([0.01, 0.03, 0.96])

    matrix = load_purity_matrix(str(purity_file))
    assert matrix.shape == (3, 3)
    assert abs(matrix[0, 0] - 0.95) < 0.001


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    from isobaric_purity_corrector import main

    channels = ["126", "127", "128"]
    purity = np.eye(3)

    purity_file = tmp_path / "purity.csv"
    with open(purity_file, "w", newline="") as fh:
        writer = csv.writer(fh)
        for row in purity:
            writer.writerow(row.tolist())

    input_file = tmp_path / "quant.tsv"
    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["spectrum_id"] + channels)
        writer.writerow(["s1", "100.0", "200.0", "300.0"])

    output_file = tmp_path / "corrected.tsv"
    sys.argv = [
        "isobaric_purity_corrector.py",
        "--input", str(input_file),
        "--label", "TMT6plex",
        "--purity-matrix", str(purity_file),
        "--output", str(output_file),
    ]
    # TMT6plex has 6 channels but our test data only has 3 columns + wrong purity size
    # Use a proper 6-channel test instead
    channels_6 = ["126", "127", "128", "129", "130", "131"]
    purity_6 = np.eye(6)
    with open(purity_file, "w", newline="") as fh:
        writer = csv.writer(fh)
        for row in purity_6:
            writer.writerow(row.tolist())
    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["spectrum_id"] + channels_6)
        writer.writerow(["s1"] + ["100.0"] * 6)

    sys.argv = [
        "isobaric_purity_corrector.py",
        "--input", str(input_file),
        "--label", "TMT6plex",
        "--purity-matrix", str(purity_file),
        "--output", str(output_file),
    ]
    main()
    assert output_file.exists()
