# Peak Picker HiRes

Pick peaks from profile (raw) mass spectra using the high-resolution peak picker.

## Usage

```bash
python peak_picker_hires.py --input profile.mzML --output centroid.mzML --signal-to-noise 1.0
```

## API

```python
from peak_picker_hires import pick_peaks

count = pick_peaks("profile.mzML", "centroid.mzML", signal_to_noise=1.0)
```
