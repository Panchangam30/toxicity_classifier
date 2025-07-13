# Neurotoxicity Module

## Purpose

Predicts CNS toxicity risk (stubbed for now).

## Implementation

- Returns a random value between 0 and 1
- TODO: Integrate real CONVERGE model

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.neurotoxicity import predict_cns_toxicity_converge
print(predict_cns_toxicity_converge({'ECFP': [...] }))
```

## Output

- Dict with 'CNS_toxicity_CONVERGE' key

## Notes

- Replace stub with real model as available
- See main README for setup and usage
