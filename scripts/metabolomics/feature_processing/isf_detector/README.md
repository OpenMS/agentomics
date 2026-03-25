# In-Source Fragmentation Detector

Detect in-source fragmentation (ISF) artifacts by identifying coeluting features whose mass difference matches common neutral losses (H2O, CO2, NH3, CO, HCOOH, CH3OH).

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python isf_detector.py --input features.tsv --rt-tolerance 3 --output isf_annotated.tsv
```

## Input format

Tab-separated file with columns: id, mz, rt, intensity.

```
id	mz	rt	intensity
F1	180.0634	120.5	50000
F2	162.0528	121.0	15000
```

## Output format

Original columns plus: isf_flag, isf_role (precursor/fragment), isf_partner, isf_loss.
