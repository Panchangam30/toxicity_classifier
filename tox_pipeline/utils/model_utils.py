import os
import joblib

def check_model_file(path):
    """
    Check if a model file exists. Raise FileNotFoundError if not.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Model file not found: {path}")
    return True

def load_joblib_model(path):
    """
    Load a joblib model from the given path, with error handling.
    """
    check_model_file(path)
    try:
        model = joblib.load(path)
        return model
    except Exception as e:
        raise RuntimeError(f"Failed to load model from {path}: {e}") 