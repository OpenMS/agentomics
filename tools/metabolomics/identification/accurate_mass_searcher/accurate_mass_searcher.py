"""
Accurate Mass Searcher
======================
Search metabolomics features against a database using accurate mass matching.

Wraps pyopenms.AccurateMassSearchEngine to identify metabolites from
featureXML input by matching observed m/z values against a mass database.

Usage
-----
    python accurate_mass_searcher.py --input features.featureXML \
        --database mapping.tsv --struct struct.tsv --output results.mzTab \
        --mass-tol 5.0
"""

import sys

import click

try:
    import pyopenms as oms
except ImportError:
    sys.exit("pyopenms is required. Install it with:  pip install pyopenms")


def search_accurate_mass(
    input_path: str,
    database_path: str,
    output_path: str,
    mass_tol: float = 5.0,
    struct_path: str = "",
    ionization_mode: str = "positive",
) -> int:
    """Search features against a mass database using AccurateMassSearchEngine.

    Parameters
    ----------
    input_path : str
        Path to input featureXML file.
    database_path : str
        Path to database mapping TSV file (mass, formula, identifier columns).
    output_path : str
        Path to output mzTab file.
    mass_tol : float
        Mass tolerance in ppm (default 5.0).
    struct_path : str
        Path to structure mapping TSV file (identifier, name, INCHI, formula).
        If empty, a minimal temporary file is created.
    ionization_mode : str
        Ionization mode: 'positive' or 'negative'.

    Returns
    -------
    int
        Number of matches found.
    """
    import os
    import tempfile

    # Load features
    fm = oms.FeatureMap()
    oms.FeatureXMLFile().load(input_path, fm)

    # Create a minimal struct file if not provided
    tmp_dir = None
    if not struct_path:
        tmp_dir = tempfile.mkdtemp()
        struct_path = os.path.join(tmp_dir, "struct.tsv")
        with open(struct_path, "w") as f:
            f.write("")  # Empty struct file

    try:
        engine = oms.AccurateMassSearchEngine()
        params = engine.getDefaults()
        params.setValue("mass_error_value", mass_tol)
        params.setValue("mass_error_unit", "ppm")
        params.setValue("ionization_mode", ionization_mode)
        params.setValue("db:mapping", [database_path])
        params.setValue("db:struct", [struct_path])
        params.setValue("keep_unidentified_masses", "false")
        engine.setParameters(params)
        engine.init()

        mztab = oms.MzTab()
        engine.run(fm, mztab)

        # Store results
        oms.MzTabFile().store(output_path, mztab)

        # Count matches by reading the mzTab file for SML lines
        match_count = 0
        with open(output_path) as f:
            for line in f:
                if line.startswith("SML\t"):
                    match_count += 1

        return match_count
    finally:
        if tmp_dir is not None:
            import shutil

            shutil.rmtree(tmp_dir, ignore_errors=True)


def create_synthetic_featurexml(output_path: str, mz: float = 181.070664, rt: float = 100.0) -> None:
    """Create a synthetic featureXML with a single feature for testing."""
    fm = oms.FeatureMap()
    feat = oms.Feature()
    feat.setRT(rt)
    feat.setMZ(mz)
    feat.setIntensity(10000.0)
    feat.setCharge(1)
    feat.setOverallQuality(1.0)
    fm.push_back(feat)
    oms.FeatureXMLFile().store(output_path, fm)


def create_synthetic_database(mapping_path: str, struct_path: str) -> None:
    """Create synthetic database files with a glucose entry for testing.

    mapping.tsv format: header lines + mass<tab>formula<tab>identifier
    struct.tsv format: identifier<tab>name<tab>inchi<tab>formula
    """
    mapping_content = (
        "database_name\tHMDB\n"
        "database_version\ttest\n"
        "180.063388\tC6H12O6\tHMDB:HMDB0000122\n"
    )
    with open(mapping_path, "w") as f:
        f.write(mapping_content)

    struct_content = "HMDB:HMDB0000122\tGlucose\tInChI=1S/C6H12O6\tC6H12O6\n"
    with open(struct_path, "w") as f:
        f.write(struct_content)


@click.command(help="Search features against a mass database for metabolite identification.")
@click.option("--input", "input_path", required=True, help="Input featureXML file")
@click.option("--database", "database_path", required=True, help="Database mapping TSV file")
@click.option("--struct", "struct_path", default="", help="Structure mapping TSV file")
@click.option("--output", "output_path", required=True, help="Output mzTab file")
@click.option("--mass-tol", default=5.0, type=float, help="Mass tolerance in ppm")
@click.option(
    "--ionization-mode",
    default="positive",
    type=click.Choice(["positive", "negative"]),
    help="Ionization mode",
)
def main(input_path, database_path, struct_path, output_path, mass_tol, ionization_mode) -> None:
    matches = search_accurate_mass(
        input_path, database_path, output_path, mass_tol, struct_path, ionization_mode
    )
    click.echo(f"Found {matches} matches, saved to {output_path}")


if __name__ == "__main__":
    main()
