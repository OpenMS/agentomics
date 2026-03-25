# Missing Value Imputation

Impute missing values in quantification matrices using MinDet, MinProb, or KNN methods.

## Usage

```bash
python missing_value_imputation.py --input matrix.tsv --method mindet --output imputed.tsv
python missing_value_imputation.py --input matrix.tsv --method knn --k 5 --output imputed.tsv
python missing_value_imputation.py --input matrix.tsv --method minprob --output imputed.tsv
```

## Methods

- **mindet** - Replace missing values with the minimum detected value per column
- **minprob** - Random draws from a low-intensity Gaussian distribution
- **knn** - K-nearest-neighbor imputation using observed features
