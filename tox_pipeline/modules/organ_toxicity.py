# Organ-Specific Toxicity Module
import numpy as np
import torch
import timm
from torchvision import transforms
from huggingface_hub import login
import os
from PIL import Image
from timm.data import resolve_data_config
from timm.data.transforms_factory import create_transform

# Login to Hugging Face if token is in environment
if "HF_TOKEN" in os.environ:
    login(token=os.environ["HF_TOKEN"])

# Load H-optimus-0 model only once
def get_hoptimus_model():
    if not hasattr(get_hoptimus_model, "model"):
        model = timm.create_model(
            "hf-hub:bioptimus/H-optimus-0", pretrained=True, init_values=1e-5, dynamic_img_size=False
        )
        device = "cuda" if torch.cuda.is_available() else "cpu"
        model.to(device)
        model.eval()
        get_hoptimus_model.model = model
        get_hoptimus_model.device = device
    return get_hoptimus_model.model, get_hoptimus_model.device

# Define the transform for H-optimus-0
hoptimus_transform = transforms.Compose([
    transforms.ToTensor(),
    transforms.Normalize(
        mean=(0.707223, 0.578729, 0.703617),
        std=(0.211883, 0.230117, 0.177517)
    ),
])

def extract_hoptimus_features(image_path):
    """
    Extract features from a histology image using H-optimus-0.
    Args:
        image_path (str): Path to a 224x224 histology image.
    Returns:
        features (np.ndarray): 1536-dim feature vector.
    """
    model, device = get_hoptimus_model()
    img = Image.open(image_path).convert("RGB").resize((224, 224))
    tensor = hoptimus_transform(img).unsqueeze(0).to(device)
    with torch.autocast(device_type=device, dtype=torch.float16 if device == "cuda" else torch.float32):
        with torch.inference_mode():
            features = model(tensor)
    return features.cpu().numpy().flatten()

# Load UNI model only once
def get_uni_model():
    if not hasattr(get_uni_model, "model"):
        model = timm.create_model(
            "hf-hub:MahmoodLab/UNI", pretrained=True, init_values=1e-5, dynamic_img_size=True
        )
        config = resolve_data_config(model.pretrained_cfg, model=model)
        transform = create_transform(**config)
        device = "cuda" if torch.cuda.is_available() else "cpu"
        model.to(device)
        model.eval()
        get_uni_model.model = model
        get_uni_model.transform = transform
        get_uni_model.device = device
    return get_uni_model.model, get_uni_model.transform, get_uni_model.device

def extract_uni_features(image_path):
    """
    Extract features from a histology image using MahmoodLab/UNI.
    Args:
        image_path (str): Path to a histology image.
    Returns:
        features (np.ndarray): 1024-dim feature vector.
    """
    model, transform, device = get_uni_model()
    img = Image.open(image_path).convert("RGB")
    tensor = transform(img).unsqueeze(0).to(device)
    with torch.inference_mode():
        features = model(tensor)
    return features.cpu().numpy().flatten()

# TODO: Replace with real H-optimus-0 model integration
# TODO: Replace with real UNI model integration
# TODO: Replace with real Merlin model integration
def predict_merlin_stub(descriptors):
    """
    Stub for Merlin model (CT-based preclinical toxicity).
    TODO: Replace with real model prediction.
    """
    return float(np.random.uniform(0, 1))


def predict_organ_toxicity(descriptors):
    """
    Organ-specific toxicity predictions using H-optimus-0 and UNI (if image provided), and stub for Merlin.
    Args:
        descriptors (dict): Should include 'histology_image_path' (str) for H-optimus-0 and UNI, or will skip if not present.
    Returns:
        dict: Organ toxicity predictions and model-specific outputs.
    """
    # TODO: Replace stubs with real model predictions for cardiotoxicity, hepatotoxicity, nephrotoxicity
    cardiotoxicity = float(np.random.uniform(0, 1))
    hepatotoxicity = float(np.random.uniform(0, 1))
    nephrotoxicity = float(np.random.uniform(0, 1))

    # H-optimus-0 feature extraction (if image path provided)
    h_optimus_0_features = None
    if 'histology_image_path' in descriptors:
        try:
            h_optimus_0_features = extract_hoptimus_features(descriptors['histology_image_path'])
        except Exception as e:
            h_optimus_0_features = f"Error: {e}"
    else:
        h_optimus_0_features = None

    # UNI feature extraction (if image path provided)
    uni_features = None
    if 'histology_image_path' in descriptors:
        try:
            uni_features = extract_uni_features(descriptors['histology_image_path'])
        except Exception as e:
            uni_features = f"Error: {e}"
    else:
        uni_features = None

    # Merlin stub
    merlin_score = predict_merlin_stub(descriptors)

    return {
        "cardiotoxicity": cardiotoxicity,
        "hepatotoxicity": hepatotoxicity,
        "nephrotoxicity": nephrotoxicity,
        "H-optimus-0_tissue_features": h_optimus_0_features,
        "UNI_rare_damage_features": uni_features,
        "Merlin_CT_toxicity_stub": merlin_score
    } 