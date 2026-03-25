# Library Coverage Estimator

Given a spectral library and a FASTA proteome, compute proteome coverage at both peptide and protein level.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python library_coverage_estimator.py --library lib.tsv --fasta proteome.fasta \
    --enzyme Trypsin --output coverage.tsv
```

### Input format

The library file is a TSV with a `PeptideSequence` or `sequence` column. The FASTA file is a standard protein FASTA.

### Parameters

| Flag | Description |
|------|-------------|
| `--library` | Spectral library TSV |
| `--fasta` | Proteome FASTA file |
| `--enzyme` | Enzyme name (default: Trypsin) |
| `--missed-cleavages` | Missed cleavages (default: 1) |
| `--output` | Output coverage TSV |
