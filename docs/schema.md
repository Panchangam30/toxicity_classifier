# Toxicity Report JSON Schema

This document defines the standard schema for the final toxicity report output (`toxicity_report.json`).

## Fields

- **toxicity_flag** (`bool`)

  - **Description:** Indicates if the compound is considered toxic based on model aggregation logic.
  - **Example:** `true`

- **organ_systems_flagged** (`list[str]`)

  - **Description:** List of organ systems flagged as at risk by any model (e.g., `["heart", "general"]`).
  - **Example:** `["heart"]`

- **num_models_flagged** (`int`)

  - **Description:** Number of models that flagged the compound as toxic.
  - **Example:** `3`

- **model_scores** (`dict[str, float|bool]`)
  - **Description:** Dictionary mapping model names to their toxicity score or flag.
  - **Example:**
    ```json
    {
      "DeepTox": 0.82,
      "hERG_Karim-CNN": 0.51,
      "Mutagenicity": false
    }
    ```

## Example: High Toxicity

```json
{
  "toxicity_flag": true,
  "organ_systems_flagged": ["heart"],
  "num_models_flagged": 3,
  "model_scores": {
    "DeepTox": 0.82,
    "hERG_Karim-CNN": 0.51,
    "Mutagenicity": 0.67
  }
}
```

## Example: Low Toxicity

```json
{
  "toxicity_flag": false,
  "organ_systems_flagged": [],
  "num_models_flagged": 1,
  "model_scores": {
    "DeepTox": 0.12,
    "hERG_Karim-CNN": 0.21,
    "Mutagenicity": 0.18
  }
}
```
