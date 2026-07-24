"""
retrain_app_models.py
======================
Re-entrena y exporta los modelos servidos por app/app.py, usando
UNICAMENTE los descriptores que la propia app puede calcular en
tiempo real a partir de un SMILES (calc_descriptors_2d + Morgan +
MACCS, definidos exactamente igual que en app.py, mas MACCS keys que
se anaden aqui porque son baratas de calcular con RDKit y no
requieren DFT).

IMPORTANTE: estos modelos son distintos (y en general algo menos
precisos) que los modelos "de investigacion" usados para las cifras
de la memoria (que usan el conjunto completo de 1097 descriptores,
incluyendo sigma-profiles COSMO-RS y descriptores 3D de Mordred que
requieren un calculo DFT previo y por tanto no se pueden calcular en
vivo para una molecula introducida por el usuario). Este es un modelo
"ligero" pensado para la demo interactiva.
"""

import os
import numpy as np
import pandas as pd
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator, MACCSkeys
from sklearn.ensemble import RandomForestRegressor, GradientBoostingClassifier
from sklearn.preprocessing import RobustScaler
from sklearn.feature_selection import SelectKBest, f_regression, f_classif
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.model_selection import cross_val_score, KFold, StratifiedKFold
import warnings
warnings.filterwarnings("ignore")

SEED = 42
N_FOLDS = 5

BASE = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")
OUT_DIR = os.path.join(BASE, "app_models_ligeros")
os.makedirs(OUT_DIR, exist_ok=True)

DOCKING_COLS = ["dG_P-gp_MDR1", "dG_CYP3A4", "dG_TfR1", "dG_FRalpha", "dG_Lisozima", "dG_HSA"]
# Nombres de fichero que espera app.py (con guion bajo en vez de guion medio)
APP_KEYS = {
    "dG_P-gp_MDR1": "dG_P_gp_MDR1",
    "dG_CYP3A4":    "dG_CYP3A4",
    "dG_TfR1":      "dG_TfR1",
    "dG_FRalpha":   "dG_FRalpha",
    "dG_Lisozima":  "dG_Lisozima",
    "dG_HSA":       "dG_HSA",
}

# ── Descriptores EXACTAMENTE como los calcula app.py ───────────
def calc_descriptors_2d(mol):
    desc = {}
    desc["MW"]             = Descriptors.MolWt(mol)
    desc["ExactMW"]        = Descriptors.ExactMolWt(mol)
    desc["LogP"]           = Descriptors.MolLogP(mol)
    desc["TPSA"]           = Descriptors.TPSA(mol)
    desc["HBD"]            = rdMolDescriptors.CalcNumHBD(mol)
    desc["HBA"]            = rdMolDescriptors.CalcNumHBA(mol)
    desc["RotBonds"]       = rdMolDescriptors.CalcNumRotatableBonds(mol)
    desc["AromaticRings"]  = rdMolDescriptors.CalcNumAromaticRings(mol)
    desc["HeavyAtomCount"] = mol.GetNumHeavyAtoms()
    desc["NumAtoms"]       = mol.GetNumAtoms()
    desc["NumBonds"]       = mol.GetNumBonds()
    desc["BertzCT"]        = Descriptors.BertzCT(mol)
    desc["Chi0"]           = Descriptors.Chi0(mol)
    desc["Chi1"]           = Descriptors.Chi1(mol)
    desc["Kappa1"]         = Descriptors.Kappa1(mol)
    desc["Kappa2"]         = Descriptors.Kappa2(mol)
    desc["FractionCSP3"]   = rdMolDescriptors.CalcFractionCSP3(mol)
    desc["MolMR"]          = Descriptors.MolMR(mol)
    return desc

def calc_fingerprints(mol):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp  = gen.GetFingerprintAsNumPy(mol)
    return {f"Morgan_{i}": int(fp[i]) for i in range(2048)}

def calc_maccs(mol):
    fp = MACCSkeys.GenMACCSKeys(mol)
    return {f"MACCS_{i}": int(fp[i]) for i in range(167)}

def build_feature_vector_dict(mol):
    d2  = calc_descriptors_2d(mol)
    fps = calc_fingerprints(mol)
    mac = calc_maccs(mol)
    return {**d2, **fps, **mac}

# ── Cargar dataset corregido y recalcular features "app-compatibles" ──
df = pd.read_csv(os.path.join(BASE, "dataset_ML.csv"))
print(f"Dataset: {df.shape[0]} moleculas")

rows = []
valid_idx = []
for i, smi in enumerate(df["smiles"]):
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        print(f"  AVISO: SMILES invalido para {df['nombre'].iloc[i]}, se omite")
        continue
    rows.append(build_feature_vector_dict(mol))
    valid_idx.append(i)

df_feat = pd.DataFrame(rows)
df = df.iloc[valid_idx].reset_index(drop=True)
df_feat = df_feat.reset_index(drop=True)
feature_cols = list(df_feat.columns)
X_raw = df_feat.values
print(f"Features app-compatibles (2D + Morgan + MACCS): {len(feature_cols)}")

def make_pipeline_reg(model, k=200):
    return Pipeline([("imputer", SimpleImputer(strategy="median")),
                      ("scaler", RobustScaler()),
                      ("selector", SelectKBest(f_regression, k=min(k, len(feature_cols)))),
                      ("model", model)])

def make_pipeline_clf(model, k=200):
    return Pipeline([("imputer", SimpleImputer(strategy="median")),
                      ("scaler", RobustScaler()),
                      ("selector", SelectKBest(f_classif, k=min(k, len(feature_cols)))),
                      ("model", model)])

# ── Modelo 1: regresion RF por proteina (igual que espera app.py: rf_*) ──
kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
print("\n--- Regresion (RF, features app-compatibles) ---")
for col in DOCKING_COLS:
    y = df[col].values
    model = RandomForestRegressor(n_estimators=200, random_state=SEED, n_jobs=-1)
    pipe = make_pipeline_reg(model)
    r2 = cross_val_score(pipe, X_raw, y, cv=kf, scoring="r2")
    print(f"  {col:18s} R2={r2.mean():.3f}+-{r2.std():.3f}")
    pipe.fit(X_raw, y)  # fit final sobre todos los datos para exportar
    app_key = APP_KEYS[col]
    joblib.dump(pipe, os.path.join(OUT_DIR, f"rf_{app_key}.joblib"))

# ── Modelo 3: clasificador GB (igual que espera app.py: gb_classifier) ──
df_iti = pd.read_csv(os.path.join(BASE, "resultados_ML", "modelo2_ITI_corregido.csv"))
df_iti_idx = df_iti.set_index("nombre")
y_clf = (df_iti_idx.loc[df["nombre"], "ITI_score"] >= 50).astype(int).values

print("\n--- Clasificador (GB, features app-compatibles) ---")
model_clf = GradientBoostingClassifier(n_estimators=200, random_state=SEED)
skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
pipe_clf = make_pipeline_clf(model_clf)
auc = cross_val_score(pipe_clf, X_raw, y_clf, cv=skf, scoring="roc_auc")
acc = cross_val_score(pipe_clf, X_raw, y_clf, cv=skf, scoring="accuracy")
print(f"  GradientBoosting  AUC={auc.mean():.3f}+-{auc.std():.3f}  Acc={acc.mean():.3f}")
pipe_clf.fit(X_raw, y_clf)
joblib.dump(pipe_clf, os.path.join(OUT_DIR, "gb_classifier.joblib"))

# ── Exportar feature_cols.joblib ──
joblib.dump(feature_cols, os.path.join(OUT_DIR, "feature_cols.joblib"))

print(f"\nOK. Modelos exportados en: {OUT_DIR}")
print("Archivos:", os.listdir(OUT_DIR))
