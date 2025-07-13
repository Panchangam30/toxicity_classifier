# Structural Alerts Module

## Purpose

Detects structural alerts (PAINS, BRENK) in small molecules using SMARTS patterns.

## Implementation

- Uses a pure-Python SMARTS-based approach for PAINS and BRENK alerts
- Does NOT use RDKit FilterCatalog (avoids C++/Python binding issues)
- Alerts are triggered if the molecule matches any SMARTS pattern
- Easy to add new patterns to `PAINS_SMARTS` or `BRENK_SMARTS` lists

## Usage

Called automatically by the main pipeline. Example:

```python
from modules.structural_alerts import check_structural_alerts
alerts = check_structural_alerts('c1ccc(cc1)[N+](=O)[O-]')
print(alerts)  # ['BRENK: BRENK_A']
```

## Output

- List of triggered alerts (e.g., ['PAINS: PAINS_A', 'BRENK: BRENK_A'])
- Empty list if no alerts

## Notes

- Add new SMARTS patterns as needed for broader coverage
- See main README for setup and usage
