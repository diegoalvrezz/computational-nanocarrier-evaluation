"""
combinar_descriptores.py  (CORREGIDO)
========================
Combina todos los datasets de descriptores en un único DataFrame
listo para los modelos de Machine Learning.

CORRECCIÓN respecto a la versión original:
    docking_energias.csv usa guiones bajos en los nombres de molécula
    ("Triethylene_glycol"), mientras que el resto de archivos usa
    espacios ("Triethylene glycol"). Al hacer merge por "nombre" sin
    normalizar, solo emparejaban 32 de 57 moléculas y el resto de la
    energía de docking se perdía silenciosamente (NaN), lo que a su
    vez hacía que el filtro de limpieza posterior descartara por
    error las familias de sigma-profiles y descriptores 3D de Mordred.
    Aquí se normaliza el nombre (guion bajo -> espacio, strip, espacios
    dobles) ANTES de cualquier merge.

Datasets:
    - descriptores_2D.csv          (60 mol x 38 desc)
    - sigma_profiles.csv           (58 mol x 100 desc)
    - fingerprints_maccs.csv       (60 mol x 167 desc)
    - fingerprints_morgan.csv      (60 mol x 2048 desc)
    - descriptores_3D_mordred.csv  (56 mol x 1826 desc)
    - docking_energias.csv         (57 mol x 6 desc)
    - base_molecular_pubchem.csv   (metadata: grupo)

Salida:
    dataset_completo.csv     -> todos los descriptores
    dataset_ML.csv           -> solo moleculas con todos los datos + limpieza
"""

import os
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')

BASE = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")
OUT_DIR = BASE

# ── Normalización de nombres ──────────────────────────────────
def norm_name(s):
    """Normaliza nombres de molécula: guion bajo -> espacio, strip,
    colapsa espacios múltiples. Debe aplicarse a TODOS los archivos
    antes de cualquier merge por 'nombre'."""
    s = str(s).strip().replace("_", " ")
    while "  " in s:
        s = s.replace("  ", " ")
    return s


def cargar_y_normalizar(path, col):
    df = pd.read_csv(path)
    df[col] = df[col].apply(norm_name)
    return df


print("Cargando datasets...")

df_2d     = cargar_y_normalizar(os.path.join(BASE, "descriptores_2D.csv"), "nombre")
df_sigma  = cargar_y_normalizar(os.path.join(BASE, "sigma_profiles.csv"), "nombre")
df_maccs  = cargar_y_normalizar(os.path.join(BASE, "fingerprints_maccs.csv"), "nombre")
df_morgan = cargar_y_normalizar(os.path.join(BASE, "fingerprints_morgan.csv"), "nombre")
df_3d     = cargar_y_normalizar(os.path.join(BASE, "descriptores_3D_mordred.csv"), "nombre")
df_dock   = cargar_y_normalizar(os.path.join(BASE, "docking_energias.csv"), "molecula")
df_dock   = df_dock.rename(columns={"molecula": "nombre"})
df_meta   = pd.read_csv(os.path.join(BASE, "base_molecular_pubchem.csv"))
df_meta   = df_meta[["nombre_entrada", "grupo"]].rename(columns={"nombre_entrada": "nombre"})
df_meta["nombre"] = df_meta["nombre"].apply(norm_name)

print(f"  2D RDKit:       {df_2d.shape}")
print(f"  Sigma profiles: {df_sigma.shape}")
print(f"  MACCS:          {df_maccs.shape}")
print(f"  Morgan:         {df_morgan.shape}")
print(f"  3D Mordred:     {df_3d.shape}")
print(f"  Docking dG:     {df_dock.shape}")

# ── Comprobación de cobertura del join (para detectar futuras
#    discrepancias de nombres antes de que se propaguen en silencio) ──
nombres_2d = set(df_2d["nombre"])
nombres_dock = set(df_dock["nombre"])
sin_match = nombres_dock - nombres_2d
if sin_match:
    print(f"\n  AVISO: {len(sin_match)} nombres de docking_energias.csv no "
          f"encuentran pareja en descriptores_2D.csv tras normalizar:")
    for n in sorted(sin_match):
        print(f"    - {n}")
else:
    print(f"\n  OK: todos los nombres de docking_energias.csv casan con "
          f"descriptores_2D.csv tras normalizar ({len(nombres_dock)} moleculas).")

# ── Merge progresivo ──────────────────────────────────────────
print("\nCombinando datasets...")

df = df_2d.copy()

df = df.merge(df_meta[["nombre", "grupo"]], on="nombre", how="left", suffixes=("", "_meta"))
if "grupo_meta" in df.columns:
    df["grupo"] = df["grupo"].fillna(df["grupo_meta"])
    df.drop(columns=["grupo_meta"], inplace=True)

df = df.merge(df_sigma, on="nombre", how="left", suffixes=("", "_sigma"))

df_maccs_clean = df_maccs.drop(columns=["grupo"], errors="ignore")
df = df.merge(df_maccs_clean, on="nombre", how="left", suffixes=("", "_maccs"))

df_morgan_clean = df_morgan.drop(columns=["grupo"], errors="ignore")
df = df.merge(df_morgan_clean, on="nombre", how="left", suffixes=("", "_morgan"))

# Nota: 4 columnas de Mordred (TPSA, BertzCT, LabuteASA, MW) coinciden en
# nombre con columnas ya presentes en descriptores_2D.csv. pandas las
# renombra automaticamente con sufijo "_3d"; se dejan así (deliberado)
# para conservar ambas versiones (RDKit-2D vs Mordred) sin colisión.
df = df.merge(df_3d, on="nombre", how="left", suffixes=("", "_3d"))

df = df.merge(df_dock, on="nombre", how="left", suffixes=("", "_dock"))

print(f"  Dataset completo: {df.shape}")

df.to_csv(os.path.join(OUT_DIR, "dataset_completo.csv"), index=False)
print(f"  OK dataset_completo.csv guardado")

# ── Dataset ML (solo moleculas con docking completo) ──────────
docking_cols = [c for c in df.columns if c.startswith("dG_")]
df_ml = df.dropna(subset=docking_cols).copy()

# Convertir a numerico todo lo que no sea metadato (ANTES de filtrar
# por NaN/varianza, para no arrastrar columnas "object" con strings
# de error como las que a veces genera Mordred)
meta_cols = ["nombre", "grupo", "smiles"]
for col in df_ml.columns:
    if col not in meta_cols:
        df_ml[col] = pd.to_numeric(df_ml[col], errors="coerce")

# Eliminar columnas con >50% NaN
threshold = len(df_ml) * 0.5
df_ml = df_ml.dropna(axis=1, thresh=threshold)

# Eliminar columnas de varianza cero
numeric_cols = df_ml.select_dtypes(include=[np.number]).columns
stds = df_ml[numeric_cols].std()
cols_a_mantener = list(stds[stds != 0].index) + [c for c in meta_cols if c in df_ml.columns] + docking_cols
cols_a_mantener = list(dict.fromkeys(cols_a_mantener))  # sin duplicados, preserva orden
df_ml = df_ml[cols_a_mantener]

df_ml.to_csv(os.path.join(OUT_DIR, "dataset_ML.csv"), index=False)

n_sigma  = sum(c.startswith("sigma_") for c in df_ml.columns)
n_morgan = sum(c.startswith("Morgan_") for c in df_ml.columns)
n_maccs  = sum(c.startswith("MACCS_") for c in df_ml.columns)
n_dock   = len(docking_cols)
n_meta   = sum(c in meta_cols for c in df_ml.columns)
n_resto  = df_ml.shape[1] - n_sigma - n_morgan - n_maccs - n_dock - n_meta

print(f"\n{'='*60}")
print(f"  Moleculas en dataset ML:  {len(df_ml)}")
print(f"  Features totales (ML):    {df_ml.shape[1] - n_dock - n_meta}")
print(f"    - sigma profiles:       {n_sigma}")
print(f"    - Morgan/ECFP4:         {n_morgan}")
print(f"    - MACCS:                {n_maccs}")
print(f"    - 2D + 3D Mordred:      {n_resto}")
print(f"\n  Columnas de docking:")
for col in docking_cols:
    vals = df_ml[col].dropna()
    print(f"    {col}: {len(vals)} valores, media={vals.mean():.2f}")
print(f"\n  OK Archivos guardados:")
print(f"    - dataset_completo.csv")
print(f"    - dataset_ML.csv")
print(f"{'='*60}")
