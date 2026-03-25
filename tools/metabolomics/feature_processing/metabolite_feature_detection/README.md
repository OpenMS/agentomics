# Metabolite Feature Detection

Detect small-molecule features in centroided LC-MS data using the
`FeatureFinderMetabo` pipeline. Output is written as a featureXML file.

## Usage

```bash
python metabolite_feature_detection.py --input sample.mzML
python metabolite_feature_detection.py --input sample.mzML --output features.featureXML --noise 1e5
```
