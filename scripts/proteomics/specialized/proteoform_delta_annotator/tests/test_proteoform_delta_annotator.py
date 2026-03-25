"""Tests for proteoform_delta_annotator."""

import csv
import sys

from conftest import requires_pyopenms


@requires_pyopenms
def test_build_ptm_mass_table():
    from proteoform_delta_annotator import build_ptm_mass_table

    table = build_ptm_mass_table()
    assert len(table) > 0
    assert all("name" in entry and "mass_shift" in entry for entry in table)


@requires_pyopenms
def test_annotate_delta_phospho():
    from proteoform_delta_annotator import annotate_delta, build_ptm_mass_table

    table = build_ptm_mass_table()
    # Phosphorylation mass shift is ~79.966 Da
    matches = annotate_delta(79.966, table, 0.5)
    assert len(matches) > 0
    # At least one match should contain "Phospho"
    names = [m["name"] for m in matches]
    assert any("Phospho" in n for n in names), f"Expected Phospho in {names}"


@requires_pyopenms
def test_annotate_delta_no_match():
    from proteoform_delta_annotator import annotate_delta, build_ptm_mass_table

    table = build_ptm_mass_table()
    matches = annotate_delta(9999.0, table, 0.01)
    assert len(matches) == 0


@requires_pyopenms
def test_annotate_proteoform_deltas():
    from proteoform_delta_annotator import annotate_proteoform_deltas

    masses = [
        ("P1_unmod", 10000.0),
        ("P1_phospho", 10079.966),
    ]
    results = annotate_proteoform_deltas(masses, tolerance=0.5)
    assert len(results) == 2
    assert results[0]["delta"] == 0.0
    assert abs(results[1]["delta"] - 79.966) < 0.001
    assert "Phospho" in results[1]["annotations"]


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    from proteoform_delta_annotator import main

    input_file = tmp_path / "input.tsv"
    output_file = tmp_path / "output.tsv"

    with open(input_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["proteoform_id", "mass"])
        writer.writerow(["P1_unmod", "10000.0"])
        writer.writerow(["P1_phospho", "10079.966"])

    sys.argv = [
        "proteoform_delta_annotator.py",
        "--input", str(input_file),
        "--tolerance", "0.5",
        "--output", str(output_file),
    ]
    main()

    assert output_file.exists()
    with open(output_file) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    assert len(rows) == 2
