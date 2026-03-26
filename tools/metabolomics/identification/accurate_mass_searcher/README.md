# Accurate Mass Searcher

Search metabolomics features against a database using accurate mass matching.

Wraps `pyopenms.AccurateMassSearchEngine` to identify metabolites from featureXML input by matching observed m/z values against a mass database.

## Usage

```bash
python accurate_mass_searcher.py --input features.featureXML \
    --database hmdb_mapping.tsv --struct hmdb_struct.tsv \
    --output results.mzTab --mass-tol 5.0
```

## Options

| Option | Default | Description |
|---|---|---|
| `--input` | required | Input featureXML file |
| `--database` | required | Database mapping TSV (mass, formula, identifier) |
| `--struct` | optional | Structure mapping TSV (identifier, name, INCHI, formula) |
| `--output` | required | Output mzTab file |
| `--mass-tol` | 5.0 | Mass tolerance in ppm |
| `--ionization-mode` | positive | Ionization mode (positive/negative) |

## Database Format

**Mapping file** (tab-separated):
```
database_name	HMDB
database_version	4.0
180.063388	C6H12O6	HMDB:HMDB0000122
```

**Structure file** (tab-separated):
```
HMDB:HMDB0000122	Glucose	InChI=1S/C6H12O6	C6H12O6
```

## Python API

```python
from accurate_mass_searcher import search_accurate_mass

matches = search_accurate_mass(
    "features.featureXML",
    "mapping.tsv",
    "results.mzTab",
    mass_tol=5.0,
    struct_path="struct.tsv",
)
print(f"Found {matches} matches")
```
