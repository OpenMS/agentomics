# Immunopeptidome QC

Quality control for immunopeptidomics data: peptide length distribution, anchor residue frequencies, and per-position information content.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python immunopeptidome_qc.py --input hla_peptides.tsv --hla-class I \
    --output length_dist.tsv --motifs anchor_freq.tsv
```

### Input format

Tab-separated file with a `sequence` column containing peptide sequences:

```
sequence
AAFGIILPK
AAGIGILTV
GILGFVFTL
```

### Parameters

| Flag | Description |
|------|-------------|
| `--input` | Input TSV with `sequence` column |
| `--hla-class` | HLA class: `I` (8-12aa) or `II` (12-25aa) |
| `--output` | Output TSV for length distribution and QC metrics |
| `--motifs` | Output TSV for anchor residue frequencies |
