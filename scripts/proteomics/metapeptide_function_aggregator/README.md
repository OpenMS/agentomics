# Metapeptide Function Aggregator

Aggregate GO/KEGG functional annotations from peptide-to-protein mappings for metaproteomics.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python metapeptide_function_aggregator.py --peptides identified.tsv \
    --annotations go_terms.tsv --output function.tsv
```

### Input format

Peptides TSV with `peptide` and `protein` columns:

```
peptide	protein
AGIILTK	P12345
PEPTIDEK	P12345;P67890
```

Annotations TSV with `protein`, `term_id`, `term_name` columns:

```
protein	term_id	term_name
P12345	GO:0006412	translation
```

### Parameters

| Flag | Description |
|------|-------------|
| `--peptides` | TSV with peptide-protein mappings |
| `--annotations` | TSV with protein-function annotations |
| `--output` | Output functional aggregation TSV |
