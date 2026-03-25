# GNPS FBMN Exporter

Export MS2 spectra and a quantification table in GNPS Feature-Based Molecular Networking (FBMN) format.

## Usage

```bash
python gnps_fbmn_exporter.py --mzml data.mzML --features features.tsv --output-mgf gnps.mgf --output-quant quant.csv
```

### Input formats

**features.tsv** (tab-separated):
```
feature_id	mz	rt	intensity
F1	180.0634	100.0	5000
```

### Output

- **MGF file**: MS2 spectra with SCANS=feature_id, PEPMASS, RTINSECONDS headers
- **Quantification CSV**: Feature table with row ID, row m/z, row retention time columns for GNPS
