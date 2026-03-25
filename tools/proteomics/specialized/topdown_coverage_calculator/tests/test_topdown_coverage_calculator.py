"""Tests for topdown_coverage_calculator."""

import csv

from conftest import requires_pyopenms


@requires_pyopenms
def test_theoretical_fragments():
    from topdown_coverage_calculator import theoretical_fragments

    frags = theoretical_fragments("ACDEFGHIK")
    assert "b" in frags and "y" in frags
    # For a 9-residue peptide, there are 8 possible b and y ions
    assert len(frags["b"]) == 8
    assert len(frags["y"]) == 8
    # b-ion numbers should be 1..8
    assert [f[0] for f in frags["b"]] == list(range(1, 9))


@requires_pyopenms
def test_match_fragments_exact():
    from topdown_coverage_calculator import match_fragments, theoretical_fragments

    seq = "PEPTIDE"
    frags = theoretical_fragments(seq)
    # Use exact theoretical masses as observed
    observed = [f[1] for f in frags["b"][:3]]  # first 3 b-ions
    matches = match_fragments(frags, observed, tolerance_ppm=10.0)
    assert len(matches["b"]) >= 3


@requires_pyopenms
def test_bond_coverage():
    from topdown_coverage_calculator import (
        bond_coverage,
        match_fragments,
        theoretical_fragments,
    )

    seq = "PEPTIDE"
    frags = theoretical_fragments(seq)
    # Match all b-ions
    observed = [f[1] for f in frags["b"]]
    matches = match_fragments(frags, observed, tolerance_ppm=10.0)
    cov = bond_coverage(seq, matches)
    assert len(cov) == len(seq) - 1
    # All bonds should be covered via b-ions
    assert all(c["covered"] for c in cov)


@requires_pyopenms
def test_coverage_summary():
    from topdown_coverage_calculator import coverage_summary

    bond_cov = [
        {"covered": True}, {"covered": True}, {"covered": False},
        {"covered": True}, {"covered": False},
    ]
    summary = coverage_summary(bond_cov)
    assert summary["total_bonds"] == 5
    assert summary["covered_bonds"] == 3
    assert abs(summary["coverage_fraction"] - 0.6) < 0.01


@requires_pyopenms
def test_cli_roundtrip(tmp_path):
    from click.testing import CliRunner
    from topdown_coverage_calculator import main, theoretical_fragments

    seq = "PEPTIDE"
    frags = theoretical_fragments(seq)

    frag_file = tmp_path / "fragments.tsv"
    output_file = tmp_path / "coverage.tsv"

    with open(frag_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["mass"])
        for _, mass in frags["b"][:3]:
            writer.writerow([f"{mass:.6f}"])

    runner = CliRunner()
    result = runner.invoke(main, [
        "--sequence", seq,
        "--fragments", str(frag_file),
        "--tolerance", "10",
        "--output", str(output_file),
    ])
    assert result.exit_code == 0, f"CLI failed: {result.output}\n{result.exception}"

    assert output_file.exists()
