# Bayesian Protein Inference

Infer proteins from peptide identifications using a Bayesian probabilistic model.

## Usage

```bash
python protein_inference_bayesian.py --input peptides.idXML --output proteins.idXML
```

## Requirements

```
pip install pyopenms click
```

## How it works

Wraps `pyopenms.BayesianProteinInferenceAlgorithm` to compute posterior probabilities for each protein given the observed peptide evidence. Handles shared peptides probabilistically, providing calibrated protein-level confidence scores.
