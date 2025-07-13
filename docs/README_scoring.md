# Scoring & Aggregation Module

## Purpose

Aggregates module outputs into a composite toxicity score and generates risk flags.

## Implementation

- Composite score formula:
  ```
  score = (
      0.15 * general_tox +
      0.2 * organ_tox_avg +
      0.15 * neurotox +
      0.1 * mito_tox +
      0.1 * morpho_tox +
      0.1 * accumulation_penalty +
      0.1 * immunotox +
      0.1 * structural_alert_penalty
  )
  ```
- Normalizes score to [0, 1]
- Adds flags for high-risk endpoints (e.g., if any risk > 0.8)
- Weights and thresholds can be tuned as needed

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.scoring import aggregate_scores
scores = aggregate_scores({
    'general_tox': {'carcinogenicity': {'prediction': 1}, 'toxcast': {'prediction': 0}},
    'organ_tox': {'cardiotoxicity': 0.9, 'hepatotoxicity': 0.2, 'nephrotoxicity': 0.1},
    'alerts': ['BRENK: BRENK_A']
})
print(scores)
```

## Output

- Dict with 'composite_score' and 'flags' keys

## Notes

- Adjust weights/thresholds as needed for your application
- See main README for setup and usage
