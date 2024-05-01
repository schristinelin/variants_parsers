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


def wrangle_tp53_data(df_path):
    files = os.listdir(df_path)
    df_giaco = pd.read_excel(
        os.path.join(
            df_path,
            [f for f in files if "giacomelli" in f.lower() and not f.startswith("~")][
                0
            ],
        ),
        sheet_name=0,
        skiprows=1,
    )
    score_cols = [col for col in df_giaco.columns if "score" in col]
    score_cols = ["Allele"] + score_cols
    df_giaco = df_giaco[score_cols].rename(columns={"Allele": "protein_variant"})

    # df_bot = pd.read_excel(
    #     os.path.join(
    #         df_path, [f for f in files if "bottcher" in f and not f.startswith("~")][0]
    #     ),
    #     sheet_name=0,
    #     index_col=0,
    # )

    # df_bot["protein_variant"] = (
    #     df_bot["Wt_aa"] + df_bot["POS"].astype(str) + df_bot["Vt_aa"]
    # )
    # df_bot = df_bot[["protein_variant", "R1_Nutlin_Ratio_Lo_Hi"]]

    # df_all = df_bot.merge(df_giaco, on="protein_variant")

    df_fayer = pd.read_excel(
        os.path.join(
            df_path, [f for f in files if "fayer" in f and not f.startswith("~")][0]
        ),
        sheet_name=3,
        skiprows=1,
    )
    df_fayer["Variant"] = df_fayer["Variant"].str.replace("p.", "")
    df_fayer = df_fayer[["Variant", "DN_reporter_score"]].rename(
        columns={"Variant": "protein_variant"}
    )

    df_all = df_fayer.merge(df_giaco, on="protein_variant")

    return df_all


def wrangle_pten_data(df_path):
    files = os.listdir(df_path)
    ## not seeing a way to automate this as different papers genearte their data in different ways
    df_matreyek = pd.read_csv(
        os.path.join(df_path, [f for f in files if "matreyek" in f][0]),
        sep="\t",
        engine="python",
    )[["variant", "score"]].rename(
        columns={"variant": "protein_variant", "score": "matreyek_score"}
    )
    df_matreyek = df_matreyek[~df_matreyek["matreyek_score"].isna()]

    df_mighell = pd.read_excel(
        os.path.join(
            df_path, [f for f in files if "mighell" in f and not f.startswith("~")][0]
        ),
        sheet_name=2,
        skiprows=1,
    )
    df_mighell = df_mighell[~df_mighell["Cum_score"].isna()]

    df_mighell = df_mighell[["Variant (one letter)", "Cum_score"]].rename(
        columns={
            "Variant (one letter)": "protein_variant",
            "Cum_score": "mighell_score",
        }
    )

    df_all = df_mighell.merge(df_matreyek, on=["protein_variant"])

    return df_all
