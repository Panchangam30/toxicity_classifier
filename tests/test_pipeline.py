import json
import os
import sys
import pytest
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils.aggregate_risk import get_toxicity_summary
from utils.format_output import format_output, validate_schema
from utils.error_handler import handle_invalid_smiles, handle_missing_model_output

def test_full_pipeline():
    """
    full pipeline: aggregate -> format -> validate -> save.
    """
    summary = get_toxicity_summary("examples")
    report = format_output(summary)
    validate_schema(report)
    with open("examples/toxicity_report.json", "w") as f:
        json.dump(report, f, indent=2)
    print("toxicity_report.json created and validated!")

def test_handle_invalid_smiles():
    assert handle_invalid_smiles("") is False
    assert handle_invalid_smiles(None) is False
    assert handle_invalid_smiles("CCO") is True

def test_handle_missing_model_output():
    missing = handle_missing_model_output("FakeModel")
    assert missing["model_name"] == "FakeModel"
    assert missing["score"] is None
    assert missing["flag"] is None
    assert "error" in missing

def test_conflicting_predictions_sets_low_confidence():
    
    fake_summary = {
        "toxicity_flag": False,
        "organ_risks": {},
        "flagged_models": 1,
        "model_results": [
            {"model_name": "A", "score": 0.9, "flag": True},
            {"model_name": "B", "score": 0.1, "flag": False},
            {"model_name": "C", "score": 0.5, "flag": False}
        ]
    }
    report = format_output(fake_summary)
    assert report["confidence_flag"] == "low"

def test_schema_validation_catches_errors():
    bad_report = {"toxicity_flag": "yes", "organ_systems_flagged": [], "num_models_flagged": 1, "model_scores": {}}
    try:
        validate_schema(bad_report)
        assert False, "Should have raised ValueError"
    except ValueError:
        pass

def test_edge_case_examples():
    with open("examples/edge_case_examples.json") as f:
        cases = json.load(f)
    for case in cases:
        validate_schema(case)
        assert case.get("confidence_flag") == "low"

if __name__ == "__main__":
    test_full_pipeline()
    test_handle_invalid_smiles()
    test_handle_missing_model_output()
    test_conflicting_predictions_sets_low_confidence()
    test_schema_validation_catches_errors()
    test_edge_case_examples()
    print("All tests passed!")
