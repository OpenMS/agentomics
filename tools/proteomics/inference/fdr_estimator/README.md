# FDR Estimator

Estimate false discovery rates for peptide and protein identifications using the target-decoy approach.

## Usage

```bash
# Peptide-level FDR
python fdr_estimator.py --input peptides.idXML --output fdr.idXML

# Include protein-level FDR
python fdr_estimator.py --input peptides.idXML --output fdr.idXML --protein
```

## Requirements

```
pip install pyopenms click
```

## How it works

Wraps `pyopenms.FalseDiscoveryRate` to compute q-values from target-decoy search results. Input PSMs must have `target_decoy` meta values set (target or decoy). The tool assigns q-values and reports counts at common FDR thresholds (1%, 5%, 10%).
