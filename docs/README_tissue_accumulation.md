# Tissue Accumulation Module

## Purpose

Predicts tissue accumulation levels for key organs (stubbed for now).

## Implementation

- Returns random values ('low', 'moderate', 'high') for liver and brain
- TODO: Integrate real pkCSM model

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.tissue_accumulation import predict_tissue_accumulation
print(predict_tissue_accumulation({'ECFP': [...] }))
```

## Output

- Dict with 'liver' and 'brain' keys

## Notes

- Replace stub with real model as available
- See main README for setup and usage
