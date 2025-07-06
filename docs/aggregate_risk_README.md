# Risk Aggregation Engine (`aggregate_risk.py`)

## Description

`aggregate_risk.py` aggregates the results from multiple toxicity prediction models and synthesizes them into a single toxicity risk summary. It applies decision logic to count flagged models, map risks to organ systems, and trigger an overall toxicity flag.

## Prerequisites

- Python 3.7+
- All model output JSON files (e.g., `output_DeepTox.json`, `output_Mutagenicity.json`, `output_hERG.json`) should be present in the same folder (e.g., `examples/`).
- Dependencies: Standard Python libraries (no special install required for aggregation).

## How to Run

From your project root, run:

```bash
python utils/aggregate_risk.py examples
```

## Input

- A folder containing individual model output JSON files, each with the format:
  ```json
  {
    "model_name": "DeepTox",
    "score": 0.82,
    "flag": true,
    "organ": "general"
  }
  ```

## Output

- Prints the aggregated summary to the console.
- Writes the final summary to `toxicity_summary.json` in the same folder.

## Example Usage

```bash
python utils/aggregate_risk.py examples
```

Example output file: `examples/toxicity_summary.json`

## Notes

- Run each model script first to generate the required output files.
- The script will skip any malformed or unrelated JSON files in the folder.
