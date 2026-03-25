# SCP Reporter QC

Single-cell proteomics QC: compute sample-to-carrier ratio per spectrum for carrier-based SCP experiments.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
python scp_reporter_qc.py --input reporter_ions.tsv --carrier-channel 131C --output qc.tsv
```

### Input format

Tab-separated file with `spectrum_id` and one column per reporter ion channel:

```
spectrum_id	126	127N	127C	128N	131C
spec1	100.5	95.2	110.3	88.7	50000.0
```

### Parameters

| Flag | Description |
|------|-------------|
| `--input` | Input TSV with reporter ion intensities |
| `--carrier-channel` | Carrier channel name (e.g. `131C`) |
| `--output` | Output QC TSV |
