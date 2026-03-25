# MS Data to CSV Exporter

Export mzML or featureXML data to flat CSV/TSV files.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Export peaks from mzML
python ms_data_to_csv_exporter.py --input data.mzML --type peaks --output peaks.tsv

# Export MS2 peaks only
python ms_data_to_csv_exporter.py --input data.mzML --type peaks --ms-level 2 --output ms2_peaks.tsv

# Export spectrum summaries
python ms_data_to_csv_exporter.py --input data.mzML --type spectra --output spectra.tsv

# Export features from featureXML
python ms_data_to_csv_exporter.py --input features.featureXML --type features --output features.tsv
```
