# Morphological Cytotoxicity Module

## Purpose

Predicts cell morphological cytotoxicity risk (stubbed for now).

## Implementation

- Returns a random value between 0 and 1
- TODO: Integrate real IMPA model

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.morphological_cytotoxicity import predict_morphological_cytotoxicity
print(predict_morphological_cytotoxicity({'ECFP': [...] }))
```

## Output

- Dict with 'morphological_cytotoxicity' key

## Notes

- Replace stub with real model as available
- See main README for setup and usage
