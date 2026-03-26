# Basic Protein Inference

Infer proteins from peptide identifications using the basic aggregation algorithm.

## Usage

```bash
python protein_inference_basic.py --input peptides.idXML --output proteins.idXML
```

## Requirements

```
pip install pyopenms click
```

## How it works

Wraps `pyopenms.BasicProteinInferenceAlgorithm` to aggregate peptide-level scores into protein-level scores. Handles shared peptides by distributing evidence across all matching proteins.
