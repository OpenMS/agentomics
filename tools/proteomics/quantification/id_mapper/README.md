# ID Mapper

Map peptide identifications onto features by matching RT and m/z values within configurable tolerances.

## Usage

```bash
# Default tolerances (RT=5s, m/z=10 Da)
python id_mapper.py --features features.featureXML --ids peptides.idXML --output annotated.featureXML

# Custom tolerances
python id_mapper.py --features features.featureXML --ids peptides.idXML --output annotated.featureXML \
    --rt-tol 10.0 --mz-tol 20.0
```

## Dependencies

```
pyopenms
click
```
