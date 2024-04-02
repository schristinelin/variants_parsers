import os

import click
import pandas as pd

from util import (alphamissense_data_pull, calc_odds_path,
                  wrangle_brca1_functional, wrangle_clinvar_txt,
                  wrangle_msh2_functional)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("gene_name", type=click.STRING) # Example: BRCA1
@click.argument("output_file_name", type=click.STRING) # Example: BRCA1
@click.option("--regen_alphamissense_data", "-am", is_flag=True, default=False, help="Set to true to regenerate alphamissense data")
def variant_parser(gene_name, output_file_name, regen_alphamissense_data):
    ## env
    wd = os.getcwd()
    clinvar_dir = os.path.join(wd, 'input_data', 'clinvar_inputs', gene_name)
    functional_data_dir = os.path.join(wd, 'input_data', 'functional_inputs', gene_name)
    am_data_dir = os.path.join(wd, 'input_data', 'alphamissence_inputs')

    # output directory
    output_data_dir = os.path.join(wd, 'output_data')
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
    clinvar_df = pd.read_csv(os.path.join(clinvar_dir, [f for f in os.listdir(clinvar_dir) if not f.startswith('.')][0]), sep = '\t', engine='python')
    functional_df_path = os.path.join(functional_data_dir, [f for f in os.listdir(functional_data_dir) if not f.startswith('.')][0])
    alphamissense_data_files = [f for f in os.listdir(am_data_dir) if not f.startswith('.') and 'tsv' in f]

    if gene_name == 'BRCA1':
        functional_df = wrangle_brca1_functional(functional_df_path)
    elif gene_name == 'MSH2':
        functional_df = wrangle_msh2_functional(functional_df_path)

    clinvar_df = wrangle_clinvar_txt(clinvar_df)

    ## if user sets flag to true, repull data regardless of existence
    # if regen_alphamissense_data:
    #     hg19_fp = os.path.join(am_data_dir, [f for f in alphamissense_data_files if 'hg19' in f][0])
    #     hg38_fp = os.path.join(am_data_dir, [f for f in alphamissense_data_files if 'hg38' in f][0]) # perhaps not great to hard-code these if there are other conditions
    #     am_df_19 = alphamissense_data_pull(hg19_fp, am_data_dir, 'hg19')
    #     am_df_38 = alphamissense_data_pull(hg38_fp, am_data_dir, 'hg38')
    #     am_df = pd.concat([am_df_19, am_df_38])
    # else:
    #     pass

    # pull alphamissense data if it doesn't exist
    if not any(fname.endswith('.csv') for fname in os.listdir(am_data_dir)) or regen_alphamissense_data:
        hg19_fp = os.path.join(am_data_dir, [f for f in alphamissense_data_files if 'hg19.tsv' in f][0])
        hg38_fp = os.path.join(am_data_dir, [f for f in alphamissense_data_files if 'hg38.tsv' in f][0]) # perhaps not great to hard-code these if there are other conditions
        am_df_19 = alphamissense_data_pull(hg19_fp, am_data_dir, 'hg19')
        am_df_38 = alphamissense_data_pull(hg38_fp, am_data_dir, 'hg38')
        am_df = pd.concat([am_df_19, am_df_38])
    else:
        am_df_19 = pd.read_csv(os.path.join(am_data_dir, [f for f in os.listdir(am_data_dir)if 'hg19.csv' in f][0]))
        am_df_38 = pd.read_csv(os.path.join(am_data_dir, [f for f in os.listdir(am_data_dir)if 'hg38.csv' in f][0]))
        am_df = pd.concat([am_df_19, am_df_38])

    
    # do calculations for clinvar data
    if gene_name == 'BRCA1': 
        classified_df = functional_df.merge(clinvar_df, on = ['transcript_variant', 'protein_variant']).drop_duplicates()
    elif gene_name == 'MSH2':
        clinvar_df.loc[clinvar_df['protein_variant'] != 'NA', 'protein_variant'] = clinvar_df['protein_variant'].str.replace('p.', '')
        classified_df = functional_df.merge(clinvar_df, on = ['protein_variant']).drop_duplicates()
    
    classified_df.to_csv(os.path.join(output_gene_data_dir, output_file_name), index= False)

    df_calc_results = pd.DataFrame(columns=['gene_name', 'func_oddspath', 'func_evidence', 'lof_oddspath', 'lof_evidence', 'predictor']) # container to store results

    # calculate oddspath for clinvar
    df_calc_results.loc[len(df_calc_results)] = sum([[gene_name], calc_odds_path(classified_df), ['clinvar']], [])

    # do calculations for alphamissense data
    if gene_name == 'BRCA1':
        functional_df_am = functional_df.copy()
        functional_df_am['protein_variant'] = functional_df_am['protein_variant'].str.replace('p.', '')
        am_calc = functional_df_am.merge(am_df, on = ['protein_variant']).drop_duplicates()
        am_calc = am_calc.rename(columns={'am_class':'classification'})
    elif gene_name == 'MSH2':
        am_calc = functional_df.merge(am_df, on = ['protein_variant']).drop_duplicates()
        am_calc = am_calc.rename(columns={'am_class':'classification'})

    # calculate oddspath for alphamissense
    df_calc_results.loc[len(df_calc_results)] = sum([[gene_name], calc_odds_path(am_calc), ['alphamissense']], [])

    print(df_calc_results)




if __name__ == "__main__":
    variant_parser()