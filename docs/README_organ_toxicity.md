# Organ-Specific Toxicity Module

## Purpose

Predicts organ-specific toxicity risks (cardiotoxicity, hepatotoxicity, nephrotoxicity) and extracts features from histology images (stubs for now).

## Implementation

- Cardiotoxicity, hepatotoxicity, nephrotoxicity: stubbed (random values)
- H-optimus-0 and UNI: feature extraction code present, but not integrated with real models (TODO)
- Merlin: stubbed (TODO)

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.organ_toxicity import predict_organ_toxicity
print(predict_organ_toxicity({'histology_image_path': 'path/to/image.png'}))
```

## Output

- Dict with per-organ risks and extracted features (if image provided)

## Notes

- Replace stubs with real models as they become available
- See main README for setup and usage
