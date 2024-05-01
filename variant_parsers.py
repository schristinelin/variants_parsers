import os

import click
import pandas as pd

from src.functional_data_curation import (
    wrangle_brca1_functional,
    wrangle_msh2_functional,
    wrangle_pten_data,
    wrangle_tp53_data,
)
from src.model_predictors_curation import (
    alphamissense_data_pull,
    eve_data_pull,
    popeve_data_pull,
    varity_data_pull,
)
from src.util import (
    calc_odds_path,
    gaussian_mixture_model,
    plot_hist_pathogenic,
    plot_scatter,
    train_kmeans_model,
    train_naive_bayes_model,
    wrangle_clinvar_txt,
)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("gene_name", type=click.STRING)  # Example: BRCA1
@click.option(
    "--regen_output_files",
    "-of",
    is_flag=True,
    default=False,
    help="Set to true to regenerate output files",
)
@click.option(
    "--regen_alphamissense_data",
    "-am",
    is_flag=True,
    default=False,
    help="Set to true to regenerate alphamissense data",
)
def variant_parser(gene_name, regen_output_files, regen_alphamissense_data):
    print("###############################################")
    print("Running pipeline for " + gene_name)
    ## envs and paths
    wd = os.getcwd()
    clinvar_dir = os.path.join(wd, "input_data", "clinvar_inputs", gene_name)
    functional_data_dir = os.path.join(wd, "input_data", "functional_inputs", gene_name)
    am_data_dir = os.path.join(wd, "input_data", "alphamissence_inputs")
    eve_data_dir = os.path.join(wd, "input_data", "eve_inputs", gene_name)
    popeve_data_dir = os.path.join(wd, "input_data", "popeve_inputs", gene_name)
    varity_data_dir = os.path.join(wd, "input_data", "varity_inputs", gene_name)

    # output directory
    output_data_dir = os.path.join(wd, "output_data")
    output_gene_data_dir = os.path.join(output_data_dir, gene_name)
    if not os.path.exists(output_data_dir):
        # create if not exist
        os.makedirs(output_data_dir)
    elif not os.path.exists(output_gene_data_dir):
        os.makedirs(output_gene_data_dir)
    else:
        pass

    ## get gene-related datasets
    # should only have 1 file for now, code needs to be modified if it ends up being multiple
    clinvar_df = pd.read_csv(
        os.path.join(
            clinvar_dir,
            [f for f in os.listdir(clinvar_dir) if not f.startswith(".")][0],
        ),
        sep="\t",
        engine="python",
    )

    alphamissense_data_files = [
        f for f in os.listdir(am_data_dir) if not f.startswith(".") and "tsv" in f
    ]

    if gene_name == "BRCA1":
        functional_df_path = os.path.join(
            functional_data_dir,
            [f for f in os.listdir(functional_data_dir) if not f.startswith(".")][0],
        )
        functional_df = wrangle_brca1_functional(functional_df_path)
    elif gene_name == "MSH2":
        functional_df_path = os.path.join(
            functional_data_dir,
            [f for f in os.listdir(functional_data_dir) if not f.startswith(".")][0],
        )
        functional_df = wrangle_msh2_functional(functional_df_path)
    elif gene_name == "TP53":
        functional_df = wrangle_tp53_data(functional_data_dir)
    elif gene_name == "PTEN":
        functional_df = wrangle_pten_data(functional_data_dir)

    clinvar_df = wrangle_clinvar_txt(clinvar_df)

    ## train the data for PTEN and TP53
    if gene_name in ["TP53", "PTEN"]:
        train_naive_bayes_model(
            functional_df, clinvar_df, output_gene_data_dir, "gaussian_nb_plot.png"
        )
        train_kmeans_model(
            functional_df, output_gene_data_dir, "kmeans_clustering_plot.png"
        )
        functional_df = gaussian_mixture_model(
            functional_df, output_gene_data_dir, "gaussian_mixture_plot.png"
        )
    else:
        pass

    # pull alphamissense data if it doesn't exist
    if (
        not any(
            fname.endswith(".csv")
            for fname in os.listdir(os.path.join(am_data_dir, gene_name))
        )
        or regen_alphamissense_data
    ):
        hg19_fp = os.path.join(
            am_data_dir, [f for f in alphamissense_data_files if "hg19.tsv" in f][0]
        )
        hg38_fp = os.path.join(
            am_data_dir, [f for f in alphamissense_data_files if "hg38.tsv" in f][0]
        )  # perhaps not great to hard-code these if there are other conditions
        am_df_19 = alphamissense_data_pull(hg19_fp, am_data_dir, "hg19", gene_name)
        am_df_38 = alphamissense_data_pull(hg38_fp, am_data_dir, "hg38", gene_name)
        am_df = pd.concat([am_df_19, am_df_38])
    else:
        am_df_19 = pd.read_csv(
            os.path.join(
                am_data_dir,
                gene_name,
                [f for f in os.listdir(am_data_dir) if "hg19.csv" in f][0],
            )
        )
        am_df_38 = pd.read_csv(
            os.path.join(
                am_data_dir,
                gene_name,
                [f for f in os.listdir(am_data_dir) if "hg38.csv" in f][0],
            )
        )
        am_df = pd.concat([am_df_19, am_df_38])

    # do calculations for clinvar data
    if gene_name == "BRCA1":
        classified_df = functional_df.merge(
            clinvar_df, on=["transcript_variant", "protein_variant"]
        ).drop_duplicates()
    elif gene_name in ["MSH2", "PTEN", "TP53"]:
        clinvar_df.loc[clinvar_df["protein_variant"] != "NA", "protein_variant"] = (
            clinvar_df["protein_variant"].str.replace("p.", "")
        )
        classified_df = functional_df.merge(
            clinvar_df, on=["protein_variant"]
        ).drop_duplicates()

    clinvar_output_dir = os.path.join(
        output_gene_data_dir, "functional_clinvar_merged_data.csv"
    )
    if not os.path.exists(clinvar_output_dir) or regen_output_files:
        classified_df.to_csv(clinvar_output_dir, index=False)

    df_calc_results = pd.DataFrame(
        columns=[
            "gene_name",
            "func_oddspath",
            "func_evidence",
            "lof_oddspath",
            "lof_evidence",
            "predictor",
        ]
    )  # container to store results

    # calculate oddspath for clinvar
    df_calc_results.loc[len(df_calc_results)] = sum(
        [[gene_name], calc_odds_path(classified_df), ["clinvar"]], []
    )

    # do calculations for alphamissense data
    if gene_name == "BRCA1":
        functional_df_am = functional_df.copy()
        functional_df_am["protein_variant"] = functional_df_am[
            "protein_variant"
        ].str.replace("p.", "")
        am_calc = functional_df_am.merge(
            am_df, on=["protein_variant"]
        ).drop_duplicates()
        y_val = "function.score.mean"
    elif gene_name == "MSH2":
        am_calc = functional_df.merge(am_df, on=["protein_variant"]).drop_duplicates()
        y_val = "lof_score"
    elif gene_name == "PTEN":
        am_calc = functional_df.merge(am_df, on=["protein_variant"]).drop_duplicates()

    am_output_dir = os.path.join(
        output_gene_data_dir, "functional_alphamissense_merged_data.csv"
    )
    if not os.path.exists(am_output_dir) or regen_output_files:
        am_calc.to_csv(am_output_dir, index=False)

    print("Histogram for Alphamissense pathogenecity saved to " + output_gene_data_dir)

    ## merge with clinvar to get threshold for pathogenecity
    clinvar_df_calc = clinvar_df.copy()
    clinvar_df_calc["protein_variant"] = clinvar_df_calc["protein_variant"].str.replace(
        "p.", ""
    )
    am_clinvar_all = am_df.merge(clinvar_df_calc, on=["protein_variant"])

    plot_hist_pathogenic(
        am_clinvar_all[
            am_clinvar_all["classification"].str.lower().str.contains("benign")
        ]["am_pathogenicity"],
        am_clinvar_all[
            am_clinvar_all["classification"].str.lower().str.contains("pathogenic")
        ]["am_pathogenicity"],
        "alphamissense",
        output_gene_data_dir,
        "alphamissense_histograms_pathogenecity.png",
    )

    # calculate oddspath for alphamissense
    ## prompt user for input
    pathogenic_threshold = input(
        "Please review the histogram for AlphaMissense and enter the threshold for pathogenic:"
    )
    benign_threshold = input(
        "Please review the histogram for AlphaMissense and enter the threshold for benign:"
    )
    ## set category based on threshold
    am_calc["classification"] = ""
    am_calc.loc[
        am_calc["am_pathogenicity"] > float(pathogenic_threshold), "classification"
    ] = "pathogenic"
    am_calc.loc[
        am_calc["am_pathogenicity"] <= float(benign_threshold), "classification"
    ] = "benign"
    am_calc["classification"] = am_calc["classification"].replace(
        r"^\s*$", "ambiguous", regex=True
    )
    # calculate
    df_calc_results.loc[len(df_calc_results)] = sum(
        [[gene_name], calc_odds_path(am_calc), ["AlphaMissense"]], []
    )
    # print(df_calc_results)

    if gene_name == "PTEN":
        print("generating scatter plot for Mighell data")
        plot_scatter(
            am_calc["am_pathogenicity"],
            am_calc["mighell_score"],
            am_calc["classification"],
            "alphamissense_pathogenicity_score",
            "functional data",
            output_gene_data_dir,
            "Alphamissense vs. functional data",
            "alphamissense_plot_mighell.png",
        )
        print("generating scatter plot for Metreyek data")
        plot_scatter(
            am_calc["am_pathogenicity"],
            am_calc["matreyek_score"],
            am_calc["classification"],
            "alphamissense_pathogenicity_score",
            "functional data",
            output_gene_data_dir,
            "Alphamissense vs. functional data",
            "alphamissense_plot_matreyek.png",
        )
    else:
        # generate scatter plot
        print("generating scatter plot")
        plot_scatter(
            am_calc["am_pathogenicity"],
            am_calc[y_val],
            am_calc["classification"],
            "alphamissense_pathogenicity_score",
            "functional data",
            output_gene_data_dir,
            "Alphamissense vs. functional data",
            "alphamissense_plot.png",
        )

    ## pull regular eve data
    eve_df = eve_data_pull(eve_data_dir)
    # calculations for regular eve data
    if gene_name == "BRCA1":
        functional_df_pe = functional_df.copy()
        functional_df_pe["protein_variant"] = functional_df_pe[
            "protein_variant"
        ].str.replace("p.", "")
        eve_calc = functional_df_pe.merge(eve_df, on=["protein_variant"])
        y_val = "function.score.mean"
    elif gene_name == "MSH2":
        eve_calc = functional_df.merge(eve_df, on=["protein_variant"])
        y_val = "lof_score"
    elif gene_name == "PTEN":
        eve_calc = functional_df.merge(eve_df, on=["protein_variant"])

    eve_output_dir = os.path.join(
        output_gene_data_dir, "functional_eve_merged_data.csv"
    )
    if not os.path.exists(eve_output_dir) or regen_output_files:
        eve_calc.to_csv(eve_output_dir, index=False)

    ## merge with clinvar to get threshold for pathogenecity
    eve_clinvar_all = eve_df.merge(clinvar_df_calc, on=["protein_variant"])

    plot_hist_pathogenic(
        eve_clinvar_all[
            eve_clinvar_all["classification"].str.lower().str.contains("benign")
        ]["EVE_scores_ASM"],
        eve_clinvar_all[
            eve_clinvar_all["classification"].str.lower().str.contains("pathogenic")
        ]["EVE_scores_ASM"],
        "EVE_scores_ASM",
        output_gene_data_dir,
        "reg_eve_histograms_pathogenecity.png",
    )
    print("Histogram for EVE pathogenecity saved to " + output_gene_data_dir)

    # calculate oddspath for regular eve
    ## prompt user for input
    pathogenic_threshold = input(
        "Please review the histogram for EVE and enter the threshold for pathogenic:"
    )
    benign_threshold = input(
        "Please review the histogram for EVE and enter the threshold for benign:"
    )
    ## set category based on threshold
    eve_calc["classification"] = ""
    eve_calc.loc[
        eve_calc["EVE_scores_ASM"] >= float(benign_threshold), "classification"
    ] = "benign"
    eve_calc.loc[
        eve_calc["EVE_scores_ASM"] <= float(pathogenic_threshold), "classification"
    ] = "pathogenic"
    eve_calc["classification"] = eve_calc["classification"].replace(
        r"^\s*$", "ambiguous", regex=True
    )
    # calculate
    df_calc_results.loc[len(df_calc_results)] = sum(
        [[gene_name], calc_odds_path(eve_calc), ["Regular EVE"]], []
    )

    if gene_name == "PTEN":
        print("generating scatter plot for Mighell data")
        plot_scatter(
            eve_calc["EVE_scores_ASM"],
            eve_calc["mighell_score"],
            eve_calc["classification"],
            "eve_score",
            "functional data",
            output_gene_data_dir,
            "Regular EVE vs. functional data",
            "reg_eve_plot_mighell.png",
        )
        print("generating scatter plot for Metreyek data")
        plot_scatter(
            eve_calc["EVE_scores_ASM"],
            eve_calc["matreyek_score"],
            eve_calc["classification"],
            "eve_score",
            "functional data",
            output_gene_data_dir,
            "Regular EVE vs. functional data",
            "reg_eve_plot_matreyek.png",
        )
    else:
        # generate scatter plot
        print("generating scatter plot")
        plot_scatter(
            eve_calc["EVE_scores_ASM"],
            eve_calc[y_val],
            eve_calc["classification"],
            "eve_score",
            "functional data",
            output_gene_data_dir,
            "Regular EVE vs. functional data",
            "reg_eve_plot.png",
        )

    ## pull popeve data
    if gene_name == "PTEN":
        print("popEVE does not have PTEN data.")
    else:
        popeve_df = popeve_data_pull(popeve_data_dir)
        # do calculations for popeve data
        if gene_name == "BRCA1":
            functional_df_pe = functional_df.copy()
            functional_df_pe["protein_variant"] = functional_df_pe[
                "protein_variant"
            ].str.replace("p.", "")
            popeve_calc = functional_df_pe.merge(popeve_df, on=["protein_variant"])
            y_val = "function.score.mean"
        elif gene_name == "MSH2":
            popeve_calc = functional_df.merge(popeve_df, on=["protein_variant"])
            y_val = "lof_score"

        popeve_output_dir = os.path.join(
            output_gene_data_dir, "functional_popeve_merged_data.csv"
        )
        if not os.path.exists(popeve_output_dir) or regen_output_files:
            popeve_calc.to_csv(popeve_output_dir, index=False)

        ## merge with clinvar to get threshold for pathogenecity
        popeve_clinvar_all = popeve_df.merge(clinvar_df_calc, on=["protein_variant"])

        plot_hist_pathogenic(
            popeve_clinvar_all[
                popeve_clinvar_all["classification"].str.lower().str.contains("benign")
            ]["popEVE"],
            popeve_clinvar_all[
                popeve_clinvar_all["classification"]
                .str.lower()
                .str.contains("pathogenic")
            ]["popEVE"],
            "popEVE",
            output_gene_data_dir,
            "popeve_histograms_pathogenecity.png",
        )
        print("Histogram for popEVE pathogenecity saved to " + output_gene_data_dir)

        # calculate oddspath for popeve
        ## prompt user for input
        pathogenic_threshold = input(
            "Please review the histogram for popEVE and enter the threshold for pathogenic:"
        )
        benign_threshold = input(
            "Please review the histogram for popEVE and enter the threshold for benign:"
        )
        ## set category based on threshold
        popeve_calc["classification"] = ""
        popeve_calc.loc[
            popeve_calc["popEVE"] >= float(benign_threshold), "classification"
        ] = "benign"
        popeve_calc.loc[
            popeve_calc["popEVE"] <= float(pathogenic_threshold), "classification"
        ] = "pathogenic"
        popeve_calc["classification"] = popeve_calc["classification"].replace(
            r"^\s*$", "ambiguous", regex=True
        )
        # calculate
        df_calc_results.loc[len(df_calc_results)] = sum(
            [[gene_name], calc_odds_path(popeve_calc), ["popEVE"]], []
        )

        # generate scatter plot
        print("generating scatter plot")
        plot_scatter(
            popeve_calc["popEVE"],
            popeve_calc[y_val],
            popeve_calc["classification"],
            "popeve_score",
            "functional data",
            output_gene_data_dir,
            "popEVE vs. functional data",
            "popeve_plot.png",
        )

    ## pull varity data
    varity_df = varity_data_pull(varity_data_dir)

    # do calculations for popeve data
    if gene_name == "BRCA1":
        functional_df_pe = functional_df.copy()
        functional_df_pe["protein_variant"] = functional_df_pe[
            "protein_variant"
        ].str.replace("p.", "")
        varity_calc = functional_df_pe.merge(varity_df, on=["protein_variant"])
        y_val = "function.score.mean"
    elif gene_name == "MSH2":
        varity_calc = functional_df.merge(varity_df, on=["protein_variant"])
        y_val = "lof_score"

    varity_output_dir = os.path.join(
        output_gene_data_dir, "functional_varity_merged_data.csv"
    )
    if not os.path.exists(varity_output_dir) or regen_output_files:
        varity_calc.to_csv(varity_output_dir, index=False)

    ## merge with clinvar to get threshold for pathogenecity
    varity_clinvar_all = varity_df.merge(clinvar_df_calc, on=["protein_variant"])

    plot_hist_pathogenic(
        varity_clinvar_all[
            varity_clinvar_all["classification"].str.lower().str.contains("benign")
        ]["VARITY_ER"],
        varity_clinvar_all[
            varity_clinvar_all["classification"].str.lower().str.contains("pathogenic")
        ]["VARITY_ER"],
        "Varity",
        output_gene_data_dir,
        "varity_histograms_pathogenecity.png",
    )
    print("Histogram for Varity pathogenecity saved to " + output_gene_data_dir)

    # calculate oddspath for popeve
    ## prompt user for input
    pathogenic_threshold = input(
        "Please review the histogram for Varity and enter the threshold for pathogenic:"
    )
    benign_threshold = input(
        "Please review the histogram for Varity and enter the threshold for benign:"
    )
    ## set category based on threshold
    varity_calc["classification"] = ""
    varity_calc.loc[
        varity_calc["VARITY_ER"] >= float(benign_threshold), "classification"
    ] = "benign"
    varity_calc.loc[
        varity_calc["VARITY_ER"] <= float(pathogenic_threshold), "classification"
    ] = "pathogenic"
    varity_calc["classification"] = varity_calc["classification"].replace(
        r"^\s*$", "ambiguous", regex=True
    )
    # calculate
    df_calc_results.loc[len(df_calc_results)] = sum(
        [[gene_name], calc_odds_path(varity_calc), ["Varity"]], []
    )

    # generate scatter plot
    print("generating scatter plot")
    plot_scatter(
        varity_calc["VARITY_ER"],
        varity_calc[y_val],
        varity_calc["classification"],
        "varity_score",
        "functional data",
        output_gene_data_dir,
        "Varity vs. functional data",
        "varity_plot.png",
    )

    print(df_calc_results)
    df_calc_results.to_csv(
        os.path.join(output_gene_data_dir, "calculated_calibrations_all.csv"),
        index=False,
    )


if __name__ == "__main__":
    variant_parser()
