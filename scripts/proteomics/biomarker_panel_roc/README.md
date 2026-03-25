# Biomarker Panel ROC

Compute ROC curves and AUC values for individual protein biomarkers and multi-marker panels.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python biomarker_panel_roc.py --input protein_quant.tsv --groups case,control --output roc.tsv
```

### Input format

Tab-separated protein quantification matrix (rows=proteins, columns=samples):

```
protein_id	case_1	case_2	control_1	control_2
P12345	100.5	120.3	50.2	45.8
```

### Parameters

| Flag | Description |
|------|-------------|
| `--input` | Protein quantification TSV |
| `--groups` | Comma-separated group names: positive,negative |
| `--group-file` | Optional TSV mapping sample_id to group |
| `--output` | Output ROC/AUC TSV |
