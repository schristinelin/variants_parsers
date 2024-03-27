
def wrangle_clinvar_txt(df):
    import pandas as pd

    from constants import classification_map, mutations_dict

    # brand new data container
    df = df[~df['Germline review status'].isin(['no assertion criteria provided', 'no classification provided'])]
    df = df[df['Name'].str.contains('\(p.')]

    df_new = pd.DataFrame()

    # built under the assumption that data retrieved from clinvar all have the same column names
    df_new['transcript_variant'] = df['Name'].str.split(':').str[1]
    df_new[['transcript_variant','protein_variant']] = df_new['transcript_variant'].str.split(' ',expand=True)
    df_new['protein_variant'] = df_new['protein_variant'].str.replace(r"\(|\)", "", regex=True)
    df_new['classification'] = df['Germline classification']
    df_new['classification'] = df_new['classification'].apply(lambda i:[k for k, v in classification_map.items() if i in v]).str[0]
    df_new = df_new[df_new['classification'].isin(['Pathogenic', 'Benign'])]

    df_new['suffix'] = df_new['protein_variant'].str.split('.').str[0]+'.'
    df_new['mid_string'] = df_new['protein_variant'].str.extract('(\d+)')
    df_new[['str1', 'str2']] = df_new['protein_variant'].str.split('\d+', expand = True)
    df_new['str1'] = df_new['str1'].str.replace('p.', '')
    df_new['protein_variant'] = df_new['suffix']+df_new['str1'].apply(lambda i:[k for k, v in mutations_dict.items() if i in v]).str[0]+df_new['mid_string']+df_new['str2'].apply(lambda i:[k for k, v in mutations_dict.items() if i in v]).str[0]


    df_new = df_new.drop(['suffix', 'mid_string', 'str1', 'str2'], axis=1)
    return(df_new)

def wrangle_brca1_functional(df):
    df = df[['gene', 'chromosome',
       'transcript_ID', 'transcript_variant', 
       'protein_variant', 'consequence', 'function.score.mean', 'func.class']]
    return(df)

def calc_odds_path(df):
    ## p1
    sample_size = len(df)
    num_pathogenic = len(df[df['classification'] == 'Pathogenic'])
    num_benign = len(df[df['classification'] == 'Benign'])

    pathogenic_prop_p1 = num_pathogenic/sample_size
    benign_prop_p1 = num_benign/sample_size

    ## p2
    df_func_norm_and_abnorm = df[df['func.class'].isin(['FUNC', 'LOF'])]
    sample_size_p2 = len(df_func_norm_and_abnorm)
    num_pathogenic_p2 = len(df_func_norm_and_abnorm[df_func_norm_and_abnorm['classification'] == 'Pathogenic'])
    num_benign_p2 = len(df_func_norm_and_abnorm[df_func_norm_and_abnorm['classification'] == 'Benign'])

    pathogenic_prop_p2 = num_pathogenic_p2/sample_size_p2
    benign_prop_p2 = num_benign_p2/sample_size_p2

    oddspath = (pathogenic_prop_p2*(1-pathogenic_prop_p1))/((1-pathogenic_prop_p2)*pathogenic_prop_p1)
    oddspath_evidence = oddspath_strength_evidence(oddspath)
    print(oddspath_evidence)
    print(oddspath)


def oddspath_strength_evidence(oddspath_val):
    if oddspath_val < 0.053:
        evidence_strength = 'BS3'
    elif oddspath_val >= 0.053 and oddspath_val < 0.23:
        evidence_strength = 'BS3_moderate'
    elif oddspath_val >= 0.23 and oddspath_val < 0.48:
        evidence_strength = 'BS3_moderate'
    elif oddspath_val >= 0.48 and oddspath_val <= 2.1:
        evidence_strength = 'Indetermine'
    elif oddspath_val > 2.1 and oddspath_val <= 4.3:
        evidence_strength = 'PS3_supporting'
    elif oddspath_val > 4.3 and oddspath_val <= 18.7:
        evidence_strength = 'PS3_moderate'
    elif oddspath_val > 18.7 and oddspath_val <= 350:
        evidence_strength = 'PS3'
    elif oddspath_val > 350:
        evidence_strength = 'PS3_very_strong'

    return evidence_strength

