
def wrangle_clinvar_txt(df):
    import pandas as pd

    from constants import classification_map, mutations_dict

    # brand new data container
    print('question 1')
    print(df['Germline review status'].unique()) ## ask about this
    df = df[df['Germline review status'].isin(['reviewed by expert panel', 'criteria provided, single submitter', 'criteria provided, multiple submitters, no conflicts'])]
    df = df[df['Name'].str.contains('\(p.')]

    df_new = pd.DataFrame()

    # built under the assumption that data retrieved from clinvar all have the same column names
    df_new['transcript_variant'] = df['Name'].str.split(':').str[1]
    df_new[['transcript_variant','protein_variant']] = df_new['transcript_variant'].str.split(' ',expand=True)
    df_new['protein_variant'] = df_new['protein_variant'].str.replace(r"\(|\)", "", regex=True)
    df_new['classification'] = df['Germline classification']
    print('question 2')
    print(df_new['classification'].unique()) # ask about the 'unknown
    df_new['classification'] = df_new['classification'].apply(lambda i:[k for k, v in classification_map.items() if i in v]).str[0]


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


