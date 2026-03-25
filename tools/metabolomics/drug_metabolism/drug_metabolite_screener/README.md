# Drug Metabolite Screener

Predict drug metabolites from Phase I/II biotransformation reactions and optionally screen mzML files for matching ions.

## Installation

```bash
pip install -r requirements.txt
```

## Usage

```bash
# Predict metabolites only
python drug_metabolite_screener.py --parent-formula C17H14ClN3O --reactions phase1,phase2 --output metabolites.tsv

# Predict and screen mzML
python drug_metabolite_screener.py --parent-formula C17H14ClN3O --reactions phase1,phase2 \
    --input run.mzML --ppm 5 --output metabolites.tsv
```

## Built-in reactions

**Phase I:** oxidation (+O), demethylation (-CH2), hydroxylation (+O), dehydrogenation (-H2), reduction (+H2)

**Phase II:** glucuronidation (+C6H8O6), sulfation (+SO3), glutathione (+C10H15N3O6S), acetylation (+C2H2O), methylation (+CH2)
