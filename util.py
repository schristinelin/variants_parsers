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


def alphamissense_data_pull(filepath, output_dir, genome_coord, gene_name):
    ## Create output directory if not exist
    if not os.path.exists(output_dir):
        # create if not exist
        os.makedirs(output_dir)
    else:
        pass

    code_dict = {
        "BRCA1": "P38398",
        "MSH2": "P43246",
        "TP53": "P04637",
    }  ## get from uniprot database: https://www.uniprot.org/uniprotkb?dir=ascend&query=%28reviewed%3Atrue%29+AND+%28gene%3ATP53%29&sort=organism_name
    import csv
    import gzip

    with gzip.open(filepath, "rt") as f:
        tsv_reader = csv.reader(f, delimiter="\t")

        # total number of lines = 69716659
        number_of_lines = 69716659
        df_header = []
        df_list = []
        # print(next(tsv_reader))

        for _ in range(number_of_lines):
            row = next(tsv_reader)
            if "#CHROM" in row[0]:
                df_header.append(row)
            elif len(row) == 10 and code_dict[gene_name] in row[5]:
                df_list.append(row)

        df_out = pd.DataFrame(df_list, columns=df_header[0])
        df_out.to_csv(
            os.path.join(
                output_dir,
                gene_name,
                str("alphamissense_data_" + genome_coord + ".csv"),
            ),
            index=False,
        )

        return df_out


def popeve_data_pull(filepath):
    df_fp = os.path.join(filepath, [f for f in os.listdir(filepath) if "csv" in f][0])
    df = pd.read_csv(df_fp)[["mutant", "popEVE"]]
    df = df.rename(columns={"mutant": "protein_variant"})
    return df


def eve_data_pull(filepath):
    df_fp = os.path.join(filepath, [f for f in os.listdir(filepath) if "csv" in f][0])
    df = pd.read_csv(df_fp)[["wt_aa", "position", "mt_aa", "EVE_scores_ASM"]]
    df["protein_variant"] = df["wt_aa"] + df["position"].astype(str) + df["mt_aa"]
    df = df[["protein_variant", "EVE_scores_ASM"]]

    return df


def varity_data_pull(filepath):
    df_fp = os.path.join(filepath, [f for f in os.listdir(filepath) if "csv" in f][0])
    df = pd.read_csv(df_fp)[["aa_pos", "aa_ref", "aa_alt", "VARITY_ER"]]
    df["protein_variant"] = df["aa_ref"] + df["aa_pos"].astype(str) + df["aa_alt"]
    df = df[["protein_variant", "VARITY_ER"]]

    return df


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
