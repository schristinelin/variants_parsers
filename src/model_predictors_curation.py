import os

import pandas as pd


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
        "PTEN": "P60484",
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
