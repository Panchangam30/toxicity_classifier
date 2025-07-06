# Model Thresholds

This document describes the decision thresholds for each toxicity prediction model used in the ensemble.

---

## DeepTox

- **Score meaning:** Probability of toxicity (0 to 1)
- **Threshold:** 0.5
- **Flag logic:** `flag = score >= 0.5`
- **Notes:** Higher score indicates higher predicted toxicity. Model: mikemayuare/SMILY-BPE-tox21

---

## Mutagenicity (Ames)

- **Score meaning:** Probability of mutagenicity (0 to 1)
- **Threshold:** 0.5
- **Flag logic:** `flag = score >= 0.5`
- **Notes:** LABEL_1 = toxic, LABEL_0 = non-toxic. Model: ML4chemistry/Toxicity_Prediction_of_Ames_test

---

## hERG Blocker

- **Score meaning:** Probability of hERG channel blocking (0 to 1)
- **Threshold:** 0.5
- **Flag logic:** `flag = score >= 0.5`
- **Notes:** Used as a proxy for cardiac toxicity. Model: hERG_Karim-CNN

---
