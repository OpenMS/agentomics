# Lipid ECN-RT Predictor

Predict lipid retention times from Equivalent Carbon Number (ECN = total_carbons - 2 * double_bonds) using linear regression per lipid class.

## Usage

```bash
python lipid_ecn_rt_predictor.py --input lipids.tsv --calibration standards.tsv --output predictions.tsv
```

### Input formats

**standards.tsv** (calibration, tab-separated):
```
lipid_class	total_carbons	double_bonds	rt
PC	32	0	10.0
PC	34	0	12.0
PC	36	0	14.0
```

**lipids.tsv** (tab-separated):
```
lipid_class	total_carbons	double_bonds
PC	34	1
```
