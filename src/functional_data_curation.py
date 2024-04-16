import os

import pandas as pd

## curated functional data not only
# define a path to save out the curated data

def wrangle_brca1_functional(df_path):
    df = pd.read_excel(df_path, sheet_name=0, skiprows=2)
    df = df[
        [
            "gene",
            "chromosome",
            "transcript_ID",
            "transcript_variant",
            "protein_variant",
            "consequence",
            "function.score.mean",
            "func.class",
        ]
    ]
    df["protein_variant"] = df["protein_variant"].fillna("NA")
    return df


def wrangle_msh2_functional(df_path):
    df = pd.read_excel(df_path, sheet_name=4)
    df = df[["Variant", "Position", "LOF score"]]
    df = df.rename(
        columns={
            "Variant": "protein_variant",
            "LOF score": "lof_score",
            "Position": "chromosome",
        }
    )
    df["protein_variant"] = df["protein_variant"].fillna("NA")
    df["lof_score"] = df["lof_score"].fillna(0)
    df["func.class"] = ""
    df.loc[df["lof_score"] > 0, "func.class"] = "LOF"
    df.loc[df["lof_score"] < 0, "func.class"] = "FUNC"
    df.loc[df["lof_score"] == 0, "func.class"] = "INT"
    return df


def wrangle_tp53_data(df_path, gene_files):
    for f in gene_files:
        f_path = os.path.join(df_path, f)
        df = pd.read_excel(f_path, sheet_name=0, skiprows=1)
        print(df)
