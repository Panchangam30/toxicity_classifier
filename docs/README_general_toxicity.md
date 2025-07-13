# General Toxicity Module

## Purpose

Predicts general toxicity endpoints using machine learning models trained on public datasets.

## Endpoints

- Carcinogenicity (Carcinogens_Lagunin): real model (Random Forest)
- ToxCast: real model (XGBoost)
- LD50: TODO (not yet implemented)

## Model Files

- `carcinogen_classifier_rf.joblib` (required)
- `toxcast_classifier_xgb.joblib` (required)

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.general_toxicity import predict_carcinogenicity, predict_toxcast
print(predict_carcinogenicity('CCO'))
print(predict_toxcast('CCO'))
```

## Output

- Dict with 'prediction' (class) and 'probability' (confidence)

## Notes

- Add new endpoints as models become available
- See main README for setup and usage
