# Proteomics Scripts

A collection of agentic-created Python scripts for common proteomics tasks
using [pyopenms](https://pyopenms.readthedocs.io/).

## Scripts

### `peptide_mass_calculator.py`
Calculate monoisotopic and average masses for peptide sequences, and compute
b-ion / y-ion fragment series.

```bash
# Plain sequence
python peptide_mass_calculator.py --sequence PEPTIDEK

# Modified sequence (oxidised Met) at charge 2
python peptide_mass_calculator.py --sequence PEPTM[147]IDEK --charge 2

# Also print fragment ions
python peptide_mass_calculator.py --sequence ACDEFGHIK --fragments
```

### `protein_digest.py`
Perform in-silico enzymatic digestion of a protein sequence and report
the resulting peptides with their masses.

```bash
# Tryptic digest
python protein_digest.py --sequence MKVLWAALLVTFLAGCQAK... --enzyme Trypsin

# Lys-C with up to 2 missed cleavages
python protein_digest.py --sequence MKVLWAALLVTFLAGCQAK... --enzyme Lys-C --missed-cleavages 2

# List all available enzymes
python protein_digest.py --list-enzymes
```

### `spectrum_file_info.py`
Summarise the contents of an mzML file: number of spectra by MS level,
retention time range, m/z range, and TIC statistics.

```bash
python spectrum_file_info.py --input sample.mzML
python spectrum_file_info.py --input sample.mzML --tic   # also print per-spectrum TIC
```

### `feature_detection_proteomics.py`
Detect peptide isotope features in centroided LC-MS/MS data using the
`FeatureFinderCentroided` algorithm.  Output is written as a featureXML file.

```bash
python feature_detection_proteomics.py --input sample.mzML
python feature_detection_proteomics.py --input sample.mzML --output features.featureXML
```
