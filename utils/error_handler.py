def handle_invalid_smiles(smiles: str) -> bool:
    """
    Checks if a SMILES string is valid. (Placeholder: returns False if empty or None)
    """
    return bool(smiles and isinstance(smiles, str) and smiles.strip())

def handle_missing_model_output(model_name: str) -> dict:
    """
    Returns a default output dict for a missing model output.
    """
    return {
        "model_name": model_name,
        "score": None,
        "flag": None,
        "error": "Model output missing"
    }

def log_error(message: str):
    """
    Logs or prints an error message.
    """
    print(f"[ERROR] {message}")
