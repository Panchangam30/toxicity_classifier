# Immunotoxicity Module

## Purpose

Predicts immunotoxicity risk (stubbed for now).

## Implementation

- Returns a random value between 0 and 1
- TODO: Integrate real immunotoxicity model

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.immunotoxicity import predict_immunotoxicity
print(predict_immunotoxicity({'ECFP': [...] }))
```

## Output

- Dict with 'immunotoxicity' key

## Notes

- Replace stub with real model as available
- See main README for setup and usage
