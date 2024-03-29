import os

import click
import pandas as pd

from util import (calc_odds_path, wrangle_brca1_functional,
                  wrangle_clinvar_txt, wrangle_msh2_functional)


@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.argument("gene_name", type=click.STRING) # Example: BRCA1
@click.argument("output_file_name", type=click.STRING) # Example: BRCA1
def variant_parser(gene_name, output_file_name):
    ## env
    wd = os.getcwd()
    clinvar_dir = os.path.join(wd, 'input_data', 'clinvar_inputs', gene_name)
    functional_data_dir = os.path.join(wd, 'input_data', 'functional_inputs', gene_name)

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

    if gene_name == 'BRCA1':
        functional_df = wrangle_brca1_functional(functional_df_path)
    elif gene_name == 'MSH2':
        functional_df = wrangle_msh2_functional(functional_df_path)

    clinvar_df = wrangle_clinvar_txt(clinvar_df)

    if gene_name == 'BRCA1': 
        classified_df = functional_df.merge(clinvar_df, on = ['transcript_variant', 'protein_variant']).drop_duplicates()
    elif gene_name == 'MSH2':
        clinvar_df.loc[clinvar_df['protein_variant'] != 'NA', 'protein_variant'] = clinvar_df['protein_variant'].str.replace('p.', '')
        classified_df = functional_df.merge(clinvar_df, on = ['protein_variant']).drop_duplicates()
    
    classified_df.to_csv(os.path.join(output_gene_data_dir, output_file_name), index= False)

    # calculate oddspath
    calc_odds_path(classified_df)




if __name__ == "__main__":
    variant_parser()