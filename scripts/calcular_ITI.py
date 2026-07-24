"""
calcular_ITI.py
===============
Indice de Eficiencia de Transporte (ITI), incluyendo el factor
corrector por peso molecular mencionado en la memoria (seccion 2.7,
Modelo 2) pero que no existia como script versionado en el
repositorio -- solo se conservaba el CSV de salida
(modelo2_ITI_corregido.csv). Reconstruido aqui a partir de la
relacion MW -> mw_factor observada en ese CSV.

Formula:
    score_targeting  = media(dG_TfR1_norm, dG_FRalpha_norm)
    score_eflujo     = media(dG_P-gp_MDR1_norm, dG_CYP3A4_norm)
    score_transporte = dG_HSA_norm
    score_biocompat  = dG_Lisozima_norm

    ITI_raw = 0.40*score_targeting + 0.25*score_transporte
            + 0.20*score_biocompat - 0.15*score_eflujo

    mw_factor:
        MW <  60           -> 0.1
        60  <= MW < 100     -> 0.4
        100 <= MW < 150     -> 0.7
        150 <= MW <= 900    -> 1.0
        MW >  900           -> 0.8

    ITI_final = ITI_raw * mw_factor  (renormalizado 0-100)

NOTA: el factor mw_factor no esta documentado en la memoria mas alla
de "incluye un factor corrector por peso molecular" -- se recomienda
anadir esta tabla de umbrales al texto (seccion 2.7) para que el
Modelo 2 sea reproducible sin depender de este script de recuperacion.
"""

import os
import pandas as pd
import numpy as np

BASE = os.path.expanduser("~/MASTER_SALUD_2026/DIEGO/Archivos_Trabajo")

DOCKING_COLS = ["dG_P-gp_MDR1", "dG_CYP3A4", "dG_TfR1", "dG_FRalpha", "dG_Lisozima", "dG_HSA"]


def mw_factor(mw):
    if mw < 60:
        return 0.1
    elif mw < 100:
        return 0.4
    elif mw < 150:
        return 0.7
    elif mw <= 900:
        return 1.0
    else:
        return 0.8


def calcular_iti(dataset_ml_path, out_path):
    df = pd.read_csv(dataset_ml_path)
    cols_needed = ["nombre", "grupo", "MW"] + DOCKING_COLS
    df_iti = df[cols_needed].copy()

    for col in DOCKING_COLS:
        mn, mx = df_iti[col].min(), df_iti[col].max()
        df_iti[col + "_norm"] = (df_iti[col] - mn) / (mx - mn)

    df_iti["score_targeting"] = df_iti[["dG_TfR1_norm", "dG_FRalpha_norm"]].mean(axis=1)
    df_iti["score_eflujo"] = df_iti[["dG_P-gp_MDR1_norm", "dG_CYP3A4_norm"]].mean(axis=1)
    df_iti["score_transporte"] = df_iti["dG_HSA_norm"]
    df_iti["score_biocompat"] = df_iti["dG_Lisozima_norm"]

    df_iti["ITI_raw"] = (
        0.40 * df_iti["score_targeting"]
        + 0.25 * df_iti["score_transporte"]
        + 0.20 * df_iti["score_biocompat"]
        - 0.15 * df_iti["score_eflujo"]
    )

    df_iti["mw_factor"] = df_iti["MW"].apply(mw_factor)
    df_iti["ITI_corr"] = df_iti["ITI_raw"] * df_iti["mw_factor"]

    mn, mx = df_iti["ITI_corr"].min(), df_iti["ITI_corr"].max()
    df_iti["ITI_score"] = ((df_iti["ITI_corr"] - mn) / (mx - mn) * 100).round(2)

    df_iti["perfil"] = pd.cut(
        df_iti["ITI_score"], bins=[0, 33, 66, 100],
        labels=["Desfavorable", "Moderado", "Favorable"], include_lowest=True,
    )

    df_iti = df_iti.sort_values("ITI_score", ascending=False)
    df_iti.to_csv(out_path, index=False)

    print(f"Top 10 ITI:")
    print(df_iti[["nombre", "ITI_score", "perfil"]].head(10).to_string(index=False))
    return df_iti


if __name__ == "__main__":
    calcular_iti(
        os.path.join(BASE, "dataset_ML.csv"),
        os.path.join(BASE, "modelo2_ITI_corregido.csv"),
    )
