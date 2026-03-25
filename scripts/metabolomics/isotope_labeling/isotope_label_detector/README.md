# Isotope Label Detector

Detect 13C/15N-labeled metabolites by pairing unlabeled and labeled features based on RT proximity and expected mass shift.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python isotope_label_detector.py --unlabeled features_ctrl.tsv --labeled features_13c.tsv \
    --tracer 13C --ppm 5 --output pairs.tsv
```

## Input format

Tab-separated files with columns: id, mz, rt.

## Output format

Tab-separated file with columns: unlabeled_id, unlabeled_mz, labeled_id, labeled_mz, mass_diff, n_labels, rt_diff, ppm_error.
