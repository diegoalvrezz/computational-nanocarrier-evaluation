"""
modelos_ML.py  (version sobre dataset_ML.csv CORREGIDO)
=========================================================
Identico en metodologia al script original (mismos modelos,
mismos hiperparametros, mismo N_FOLDS=5, mismo SEED=42,
mismo SelectKBest k=200) pero ejecutado sobre el dataset_ML.csv
regenerado por combinar_descriptores.py (nombres normalizados,
sin duplicados de fingerprints, con sigma-profiles y Mordred 3D
realmente incluidos).

Salida: resultados_ML/
    modelo1_regresion_resultados.csv
    modelo2_ITI.csv        (usa calcular_ITI.py aparte para el ITI con mw_factor)
    modelo3_clasificador_resultados.csv
    shap_importancia_Pgp.csv
"""

import os
import warnings
import numpy as np
import pandas as pd
from sklearn.ensemble import (RandomForestRegressor, RandomForestClassifier,
                               GradientBoostingRegressor, GradientBoostingClassifier)
from sklearn.svm import SVR, SVC
from sklearn.neural_network import MLPRegressor, MLPClassifier
from sklearn.preprocessing import RobustScaler
from sklearn.impute import SimpleImputer
from sklearn.model_selection import cross_val_score, KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_regression, f_classif

warnings.filterwarnings('ignore')

BASE      = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")
INPUT_CSV = os.path.join(BASE, "dataset_ML.csv")
OUT_DIR   = os.path.join(BASE, "resultados_ML")
os.makedirs(OUT_DIR, exist_ok=True)

DOCKING_COLS = ["dG_P-gp_MDR1", "dG_CYP3A4", "dG_TfR1", "dG_FRalpha", "dG_Lisozima", "dG_HSA"]
META_COLS    = ["nombre", "grupo", "smiles"]
SEED         = 42
N_FOLDS      = 5

print("=" * 60)
print("MODELOS ML - TFM Diego Vallina (dataset corregido)")
print("=" * 60)

df = pd.read_csv(INPUT_CSV)
print(f"\nDataset: {df.shape[0]} moleculas x {df.shape[1]} columnas")

feature_cols = [c for c in df.columns
                if c not in META_COLS + DOCKING_COLS
                and df[c].dtype in [np.float64, np.int64, float, int]]
X_raw = df[feature_cols].values
print(f"Features disponibles: {len(feature_cols)}")
print(f"  (incluye sigma-profiles y descriptores 3D de Mordred, "
      f"perdidos en la version anterior del pipeline)")

# ── Pipeline de preprocesado (se anade Imputer: el dataset corregido
#    conserva 1 molecula con NaN residuales en columnas Mordred) ──
def make_pipeline(model, k_features=200):
    return Pipeline([
        ("imputer",  SimpleImputer(strategy="median")),
        ("scaler",   RobustScaler()),
        ("selector", SelectKBest(f_regression, k=min(k_features, len(feature_cols)))),
        ("model",    model),
    ])

def make_pipeline_clf(model, k_features=200):
    return Pipeline([
        ("imputer",  SimpleImputer(strategy="median")),
        ("scaler",   RobustScaler()),
        ("selector", SelectKBest(f_classif, k=min(k_features, len(feature_cols)))),
        ("model",    model),
    ])

# ── MODELO 1: Regresion de afinidad dG ────────────────────────
print("\n" + "=" * 60)
print("MODELO 1: Regresion de afinidad dG por proteina")
print("=" * 60)

modelos_reg = {
    "Random Forest":    RandomForestRegressor(n_estimators=200, random_state=SEED, n_jobs=-1),
    "GradientBoosting": GradientBoostingRegressor(n_estimators=200, random_state=SEED),
    "SVR":              SVR(kernel="rbf", C=10, epsilon=0.1),
    "MLP":              MLPRegressor(hidden_layer_sizes=(256, 128, 64), max_iter=500,
                                      random_state=SEED, early_stopping=True),
}

kf = KFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
resultados_reg = []

for prot in DOCKING_COLS:
    y = df[prot].values
    print(f"\n  Proteina: {prot}")
    for nombre_modelo, modelo in modelos_reg.items():
        pipe = make_pipeline(modelo)
        r2_scores   = cross_val_score(pipe, X_raw, y, cv=kf, scoring="r2")
        mae_scores  = -cross_val_score(pipe, X_raw, y, cv=kf, scoring="neg_mean_absolute_error")
        rmse_scores = np.sqrt(-cross_val_score(pipe, X_raw, y, cv=kf, scoring="neg_mean_squared_error"))
        resultados_reg.append({
            "proteina": prot, "modelo": nombre_modelo,
            "R2_mean": round(r2_scores.mean(), 4), "R2_std": round(r2_scores.std(), 4),
            "MAE_mean": round(mae_scores.mean(), 4), "RMSE_mean": round(rmse_scores.mean(), 4),
        })
        print(f"    {nombre_modelo:20s} R2={r2_scores.mean():.3f}+-{r2_scores.std():.3f}  "
              f"MAE={mae_scores.mean():.3f}  RMSE={rmse_scores.mean():.3f}")

df_reg = pd.DataFrame(resultados_reg)
df_reg.to_csv(os.path.join(OUT_DIR, "modelo1_regresion_resultados.csv"), index=False)
print(f"\n  OK modelo1_regresion_resultados.csv guardado")

# ── MODELO 2: ITI -- delegado a calcular_ITI.py (usa columna MW) ──
print("\n" + "=" * 60)
print("MODELO 2: Indice de Transporte (ITI) -- ver calcular_ITI.py")
print("=" * 60)
from calcular_ITI import calcular_iti
df_iti = calcular_iti(INPUT_CSV, os.path.join(OUT_DIR, "modelo2_ITI_corregido.csv"))

# ── MODELO 3: Clasificador favorable/desfavorable ─────────────
print("\n" + "=" * 60)
print("MODELO 3: Clasificador favorable/desfavorable")
print("=" * 60)

y_clf = (df_iti.set_index("nombre").loc[df["nombre"], "ITI_score"] >= 50).astype(int).values
print(f"\n  Distribucion clases: Favorable={y_clf.sum()} | Desfavorable={len(y_clf)-y_clf.sum()}")

modelos_clf = {
    "Random Forest":    RandomForestClassifier(n_estimators=200, random_state=SEED, n_jobs=-1),
    "GradientBoosting": GradientBoostingClassifier(n_estimators=200, random_state=SEED),
    "SVC":              SVC(kernel="rbf", C=10, probability=True, random_state=SEED),
    "MLP":              MLPClassifier(hidden_layer_sizes=(256, 128), max_iter=500,
                                       random_state=SEED, early_stopping=True),
}

skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=SEED)
resultados_clf = []

for nombre_modelo, modelo in modelos_clf.items():
    pipe = make_pipeline_clf(modelo)
    acc_scores = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="accuracy")
    auc_scores = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="roc_auc")
    f1_scores  = cross_val_score(pipe, X_raw, y_clf, cv=skf, scoring="f1")
    resultados_clf.append({
        "modelo": nombre_modelo,
        "Accuracy": round(acc_scores.mean(), 4), "Acc_std": round(acc_scores.std(), 4),
        "ROC_AUC": round(auc_scores.mean(), 4), "AUC_std": round(auc_scores.std(), 4),
        "F1": round(f1_scores.mean(), 4), "F1_std": round(f1_scores.std(), 4),
    })
    print(f"  {nombre_modelo:20s} Acc={acc_scores.mean():.3f}+-{acc_scores.std():.3f}  "
          f"AUC={auc_scores.mean():.3f}  F1={f1_scores.mean():.3f}")

df_clf = pd.DataFrame(resultados_clf)
df_clf.to_csv(os.path.join(OUT_DIR, "modelo3_clasificador_resultados.csv"), index=False)
print(f"\n  OK modelo3_clasificador_resultados.csv guardado")

# ── SHAP para P-gp (Random Forest) ─────────────────────────────
print("\n" + "=" * 60)
print("SHAP: importancia de descriptores para P-gp MDR1")
print("=" * 60)
try:
    import shap
    y_pgp = df["dG_P-gp_MDR1"].values
    pipe = make_pipeline(RandomForestRegressor(n_estimators=200, random_state=SEED, n_jobs=-1))
    pipe.fit(X_raw, y_pgp)
    selected_mask = pipe.named_steps["selector"].get_support()
    selected_features = np.array(feature_cols)[selected_mask]
    X_selected = pipe.named_steps["scaler"].transform(
        pipe.named_steps["imputer"].transform(X_raw))[:, selected_mask]
    explainer = shap.TreeExplainer(pipe.named_steps["model"])
    shap_values = explainer.shap_values(X_selected)
    shap_mean_abs = np.abs(shap_values).mean(axis=0)
    df_shap = pd.DataFrame({"feature": selected_features, "shap_mean_abs": shap_mean_abs})
    df_shap = df_shap.sort_values("shap_mean_abs", ascending=False)
    df_shap.to_csv(os.path.join(OUT_DIR, "shap_importancia_Pgp.csv"), index=False)
    print(df_shap.head(15).to_string(index=False))
except ImportError:
    print("  (paquete 'shap' no disponible en este entorno; instalar con pip install shap)")

print("\n" + "=" * 60)
print("HECHO. Comparar resultados_ML/*.csv con las cifras de la memoria.")
print("=" * 60)
