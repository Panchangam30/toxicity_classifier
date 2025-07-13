# general_toxicity_train.py

from tdc.single_pred import Tox
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
import joblib
import requests
import pubchempy as pcp
import os
import pandas as pd


def smiles_to_ecfp(smiles, radius=2, n_bits=2048):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return np.zeros(n_bits)
    return np.array(rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits))

def name_to_smiles_opsin_web(name):
    try:
        url = f'https://opsin.ch.cam.ac.uk/opsin/{requests.utils.quote(str(name))}.smi'
        response = requests.get(url, timeout=10)
        if response.status_code == 200 and response.text.strip():
            return response.text.strip()
    except Exception:
        pass
    return None

def load_carc_cache():
    if os.path.exists(CARC_CACHE_FILE):
        cache_df = pd.read_csv(CARC_CACHE_FILE)
        return dict(zip(cache_df['Drug'], cache_df['smiles']))
    else:
        return {}

def save_carc_cache(cache):
    pd.DataFrame(list(cache.items()), columns=['Drug', 'smiles']).to_csv(CARC_CACHE_FILE, index=False)

def name_to_smiles(name, cache):
    if name in cache:
        return cache[name]
    print(f"Querying Drug: {name}")
    # Try PubChem name
    try:
        compounds = pcp.get_compounds(name, 'name')
        if compounds and compounds[0].isomeric_smiles:
            smiles = compounds[0].isomeric_smiles
            cache[name] = smiles
            save_carc_cache(cache)
            return smiles
    except Exception:
        pass
    # Try OPSIN web
    smiles = name_to_smiles_opsin_web(name)
    if smiles:
        cache[name] = smiles
        save_carc_cache(cache)
        return smiles
    cache[name] = None
    save_carc_cache(cache)
    return None

# 1. Carcinogens_Lagunin
print("Loading Carcinogens_Lagunin...")
carc = Tox(name='Carcinogens_Lagunin')
carc_data = carc.get_data()
print(f"Rows: {len(carc_data)}, Columns: {carc_data.columns}")

# For Carcinogens_Lagunin, the 'Drug' column is already SMILES
carc_data['smiles'] = carc_data['Drug']
print(f"Number of drugs with valid SMILES: {carc_data['smiles'].notnull().sum()}")

print("Featurizing Carcinogens_Lagunin...")
X_carc = np.stack(carc_data['smiles'].apply(smiles_to_ecfp))
y_carc = carc_data['Y'].values

X_train, X_test, y_train, y_test = train_test_split(X_carc, y_carc, test_size=0.2, random_state=42)

from sklearn.model_selection import GridSearchCV
# from sklearn.ensemble import RandomForestClassifier
from lightgbm import LGBMClassifier
from xgboost import XGBClassifier

# param_grid_rf = {
#     'n_estimators': [100, 200, 500],
#     'max_depth': [None, 10, 20],
#     'max_features': ['sqrt', 'log2', None]
# }
# clf_carc = RandomForestClassifier(random_state=42)
# grid_search_carc = GridSearchCV(clf_carc, param_grid_rf, cv=5, n_jobs=-1, verbose=2)
# grid_search_carc.fit(X_train, y_train)
# print("Best parameters (Carcinogen RF):", grid_search_carc.best_params_)
# print("Best CV score (Carcinogen RF):", grid_search_carc.best_score_)
# best_clf_carc = grid_search_carc.best_estimator_
# print(f"Carcinogen classifier test accuracy (RF): {best_clf_carc.score(X_test, y_test):.3f}")
# joblib.dump(best_clf_carc, "carcinogen_classifier_rf.joblib")

# LightGBM tuning
param_grid_lgbm = {
    'num_leaves': [15, 31, 63],
    'n_estimators': [100, 200, 500],
    'max_depth': [-1, 10, 20]
}
# Set verbose=-1 to suppress LightGBM warnings
clf_carc_lgbm = LGBMClassifier(random_state=42, verbose=-1)
grid_search_carc_lgbm = GridSearchCV(clf_carc_lgbm, param_grid_lgbm, cv=5, n_jobs=-1, verbose=2)
grid_search_carc_lgbm.fit(X_train, y_train)
print("Best parameters (Carcinogen LGBM):", grid_search_carc_lgbm.best_params_)
print("Best CV score (Carcinogen LGBM):", grid_search_carc_lgbm.best_score_)
best_clf_carc_lgbm = grid_search_carc_lgbm.best_estimator_
print(f"Carcinogen classifier test accuracy (LGBM): {best_clf_carc_lgbm.score(X_test, y_test):.3f}")
joblib.dump(best_clf_carc_lgbm, "carcinogen_classifier_lgbm.joblib")

# XGBoost tuning
param_grid_xgb = {
    'n_estimators': [100, 200, 500],
    'max_depth': [3, 6, 10],
    'learning_rate': [0.01, 0.1, 0.3]
}
clf_carc_xgb = XGBClassifier(random_state=42, use_label_encoder=False, eval_metric='logloss', verbosity=0)
grid_search_carc_xgb = GridSearchCV(clf_carc_xgb, param_grid_xgb, cv=5, n_jobs=-1, verbose=2)
grid_search_carc_xgb.fit(X_train, y_train)
print("Best parameters (Carcinogen XGB):", grid_search_carc_xgb.best_params_)
print("Best CV score (Carcinogen XGB):", grid_search_carc_xgb.best_score_)
best_clf_carc_xgb = grid_search_carc_xgb.best_estimator_
print(f"Carcinogen classifier test accuracy (XGB): {best_clf_carc_xgb.score(X_test, y_test):.3f}")
joblib.dump(best_clf_carc_xgb, "carcinogen_classifier_xgb.joblib")

# 2. ToxCast (pick a target, e.g., 'ACEA_T47D_80hr_Negative')
print("Loading ToxCast...")
toxcast = Tox(name='ToxCast', label_name='ACEA_T47D_80hr_Negative')
toxcast_data = toxcast.get_data()
print(f"Rows: {len(toxcast_data)}, Columns: {toxcast_data.columns}")

# The 'Drug' column is SMILES
toxcast_data['smiles'] = toxcast_data['Drug']

print("Featurizing ToxCast...")
X_toxcast = np.stack(toxcast_data['smiles'].apply(smiles_to_ecfp))
y_toxcast = toxcast_data['Y'].values

X_train, X_test, y_train, y_test = train_test_split(X_toxcast, y_toxcast, test_size=0.2, random_state=42)

# clf_toxcast = RandomForestClassifier(random_state=42)
# grid_search_toxcast = GridSearchCV(clf_toxcast, param_grid_rf, cv=5, n_jobs=-1, verbose=2)
# grid_search_toxcast.fit(X_train, y_train)
# print("Best parameters (ToxCast RF):", grid_search_toxcast.best_params_)
# print("Best CV score (ToxCast RF):", grid_search_toxcast.best_score_)
# best_clf_toxcast = grid_search_toxcast.best_estimator_
# print(f"ToxCast classifier test accuracy (RF): {best_clf_toxcast.score(X_test, y_test):.3f}")
# joblib.dump(best_clf_toxcast, "toxcast_classifier_rf.joblib")

clf_toxcast_lgbm = LGBMClassifier(random_state=42, verbose=-1)
grid_search_toxcast_lgbm = GridSearchCV(clf_toxcast_lgbm, param_grid_lgbm, cv=5, n_jobs=-1, verbose=2)
grid_search_toxcast_lgbm.fit(X_train, y_train)
print("Best parameters (ToxCast LGBM):", grid_search_toxcast_lgbm.best_params_)
print("Best CV score (ToxCast LGBM):", grid_search_toxcast_lgbm.best_score_)
best_clf_toxcast_lgbm = grid_search_toxcast_lgbm.best_estimator_
print(f"ToxCast classifier test accuracy (LGBM): {best_clf_toxcast_lgbm.score(X_test, y_test):.3f}")
joblib.dump(best_clf_toxcast_lgbm, "toxcast_classifier_lgbm.joblib")

clf_toxcast_xgb = XGBClassifier(random_state=42, use_label_encoder=False, eval_metric='logloss', verbosity=0)
grid_search_toxcast_xgb = GridSearchCV(clf_toxcast_xgb, param_grid_xgb, cv=5, n_jobs=-1, verbose=2)
grid_search_toxcast_xgb.fit(X_train, y_train)
print("Best parameters (ToxCast XGB):", grid_search_toxcast_xgb.best_params_)
print("Best CV score (ToxCast XGB):", grid_search_toxcast_xgb.best_score_)
best_clf_toxcast_xgb = grid_search_toxcast_xgb.best_estimator_
print(f"ToxCast classifier test accuracy (XGB): {best_clf_toxcast_xgb.score(X_test, y_test):.3f}")
joblib.dump(best_clf_toxcast_xgb, "toxcast_classifier_xgb.joblib") 