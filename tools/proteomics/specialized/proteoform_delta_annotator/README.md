# Proteoform Delta Annotator

Annotate mass differences between proteoforms with known PTMs from the pyopenms ModificationsDB.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python proteoform_delta_annotator.py --input proteoform_masses.tsv --tolerance 0.5 --output annotated.tsv
```

### Input format

Tab-separated file with `proteoform_id` and `mass` columns:

```
proteoform_id	mass
P1_unmod	12345.678
P1_phospho	12425.644
```

### Parameters

| Flag | Description |
|------|-------------|
| `--input` | Input TSV with proteoform IDs and masses |
| `--tolerance` | Mass tolerance in Da (default: 0.5) |
| `--output` | Output annotated TSV |
