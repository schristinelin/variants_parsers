import os

import pandas as pd

from constants import classification_map, mutations_dict


def wrangle_clinvar_txt(df):
    # brand new data container
    df = df[
        ~df["Germline review status"].isin(
            ["no assertion criteria provided", "no classification provided"]
        )
    ]
    # df = df[df['Name'].str.contains('\(p.')] ## do not filter this out, put NA in occurences that are null
    df_new = pd.DataFrame()

    # built under the assumption that data retrieved from clinvar all have the same column names
    df_new["transcript_variant"] = df["Name"].str.split(":").str[1]
    df_new[["transcript_variant", "protein_variant"]] = df_new[
        "transcript_variant"
    ].str.split(" ", expand=True)
    df_new["protein_variant"] = df_new["protein_variant"].fillna("NA")
    df_new["protein_variant"] = df_new["protein_variant"].str.replace(
        r"\(|\)", "", regex=True
    )
    df_new["classification"] = df["Germline classification"]
    df_new["classification"] = (
        df_new["classification"]
        .apply(lambda i: [k for k, v in classification_map.items() if i in v])
        .str[0]
    )
    df_new = df_new[df_new["classification"].isin(["Pathogenic", "Benign"])]
    df_new[["suffix", "mid_string", "str1", "str2"]] = ""
    df_new.loc[(df_new["protein_variant"] != "NA"), "suffix"] = (
        df_new.loc[df_new["protein_variant"] != "NA", "protein_variant"]
        .str.split(".")
        .str[0]
        + "."
    )
    df_new["mid_string"] = df_new.loc[
        df_new["protein_variant"] != "NA", "protein_variant"
    ].str.extract("(\d+)")
    df_new[["str1", "str2"]] = (
        df_new.loc[df_new["protein_variant"] != "NA", "protein_variant"]
        .str.split("\d+", expand=True)
        .iloc[:, 0:2]
    )
    df_new["str1"] = df_new["str1"].str.replace("p.", "")
    df_new = df_new[df_new["protein_variant"] != "p.?"]  ## what's with these?
    df_new.loc[(df_new["protein_variant"] != "NA"), "protein_variant"] = (
        df_new.loc[(df_new["protein_variant"] != "NA"), "suffix"]
        + df_new.loc[(df_new["protein_variant"] != "NA"), "str1"]
        .apply(lambda i: [k for k, v in mutations_dict.items() if i in v])
        .str[0]
        + df_new.loc[(df_new["protein_variant"] != "NA"), "mid_string"]
        + df_new.loc[(df_new["protein_variant"] != "NA"), "str2"]
        .apply(lambda i: [k for k, v in mutations_dict.items() if i in v])
        .str[0]
    )

    df_new = df_new.drop(["suffix", "mid_string", "str1", "str2"], axis=1)
    return df_new




def map_genome_codes(str_val):
    import re

    suffix = str_val.split(".")[0] + "."
    mid_string = re.findall("\d+", str_val)[0]
    print(str_val.split("\d+"))
    str1, str2 = str_val.split("\d+")[0], str_val.split("\d+")[1]
    str1 = str1.replace("p.", "")
    new_str_val = (
        suffix
        + str1.apply(lambda i: [k for k, v in mutations_dict.items() if i in v]).str[0]
        + mid_string
        + str2.apply(lambda i: [k for k, v in mutations_dict.items() if i in v]).str[0]
    )

    return new_str_val


def calc_odds_path(df):
    ## p1
    sample_size = len(df)
    num_pathogenic = len(
        df[df["classification"].str.lower().str.contains("pathogenic")]
    )
    num_benign = len(df[df["classification"].str.lower().str.contains("benign")])

    pathogenic_prop_p1 = num_pathogenic / sample_size
    benign_prop_p1 = num_benign / sample_size

    ## p2
    ## FUNC
    func_df = df[df["func.class"] == "FUNC"]
    sample_size_1 = len(func_df)
    num_pathogenic_func = len(
        func_df[func_df["classification"].str.lower().str.contains("pathogenic")]
    )
    num_benign_func = len(
        func_df[func_df["classification"].str.lower().str.contains("benign")]
    )

    pathogenic_prop_func_p2 = num_pathogenic_func / sample_size_1
    benign_prop_func_p2 = num_benign_func / sample_size_1

    ## LOF
    lof_df = df[df["func.class"] == "LOF"]
    sample_size_2 = len(lof_df)
    num_pathogenic_lof = len(
        lof_df[lof_df["classification"].str.lower().str.contains("pathogenic")]
    )
    num_benign_lof = len(
        lof_df[lof_df["classification"].str.lower().str.contains("benign")]
    )

    pathogenic_prop_lof_p2 = num_pathogenic_lof / sample_size_2
    benign_prop_lof_p2 = num_pathogenic_lof / sample_size_2

    ## func oddspath
    func_oddspath = (pathogenic_prop_func_p2 * (1 - pathogenic_prop_p1)) / (
        (1 - pathogenic_prop_func_p2) * pathogenic_prop_p1
    )
    func_oddspath_evidence = oddspath_strength_evidence(func_oddspath)

    ## lof oddspath
    lof_oddspath = (pathogenic_prop_lof_p2 * (1 - pathogenic_prop_p1)) / (
        (1 - pathogenic_prop_lof_p2) * pathogenic_prop_p1
    )
    lof_oddspath_evidence = oddspath_strength_evidence(lof_oddspath)

    row_list = [
        func_oddspath,
        func_oddspath_evidence,
        lof_oddspath,
        lof_oddspath_evidence,
    ]

    return row_list


def oddspath_strength_evidence(oddspath_val):
    if oddspath_val < 0.053:
        evidence_strength = "BS3"
    elif oddspath_val >= 0.053 and oddspath_val < 0.23:
        evidence_strength = "BS3_moderate"
    elif oddspath_val >= 0.23 and oddspath_val < 0.48:
        evidence_strength = "BS3_supporting"
    elif oddspath_val >= 0.48 and oddspath_val <= 2.1:
        evidence_strength = "Indetermine"
    elif oddspath_val > 2.1 and oddspath_val <= 4.3:
        evidence_strength = "PS3_supporting"
    elif oddspath_val > 4.3 and oddspath_val <= 18.7:
        evidence_strength = "PS3_moderate"
    elif oddspath_val > 18.7 and oddspath_val <= 350:
        evidence_strength = "PS3"
    elif oddspath_val > 350:
        evidence_strength = "PS3_very_strong"

    return evidence_strength



def plot_scatter(x, y, color, x_label, fpath, title, fname):
    import matplotlib.patches
    import matplotlib.pyplot as plt

    levels, categories = pd.factorize(color)
    colors = [plt.cm.tab10(i) for i in levels]  # using the "tab10" colormap
    handles = [
        matplotlib.patches.Patch(color=plt.cm.tab10(i), label=c)
        for i, c in enumerate(categories)
    ]

    plt.scatter(x, y, c=colors)
    plt.legend(handles=handles)
    plt.gca().set(title=title, xlabel=x_label, ylabel="functional data")
    plt.savefig(os.path.join(fpath, str(fname)))
    plt.clf()


def plot_hist_pathogenic(score1, score2, model_name, fpath, fname):
    import matplotlib.pyplot as plt

    plt.hist(score1, alpha=0.5, label="benign")
    plt.hist(score2, alpha=0.5, label="pathogenic")
    plt.legend(loc="upper right")
    plt.gca().set(title="Pathogenecity Histogram of " + model_name)
    plt.savefig(os.path.join(fpath, str(fname)))
    plt.clf()