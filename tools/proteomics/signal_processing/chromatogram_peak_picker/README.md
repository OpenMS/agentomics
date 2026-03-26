# Chromatogram Peak Picker

Pick peaks from chromatographic data using pyopenms PeakPickerChromatogram.

## Usage

```bash
python chromatogram_peak_picker.py --input run.mzML --output picked.mzML
```

## API

```python
from chromatogram_peak_picker import pick_chromatogram_peaks

count = pick_chromatogram_peaks("run.mzML", "picked.mzML")
```
