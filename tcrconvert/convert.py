import pandas as pd
from importlib.resources import files

# Standard column names for different sources of TCR data
col_ref = {'adaptive': ('v_resolved', 'd_resolved', 'j_resolved'),
           'adaptive_v2': ('vMaxResolved', 'dMaxResolved', 'jMaxResolved'),
           'imgt': ('v_gene', 'd_gene', 'j_gene'),
           'tenx': ('v_gene', 'd_gene', 'j_gene')}


# TODO: gracefully handle missing/NA values (shouldn't be printed as not in reference)
def convert_vdj(df, frm, to, frm_cols=None, species='human'):
    '''Convert TCR V, D, and J gene names from one naming convention to another.

    :param df: Dataframe of TCRs with gene names
    :type df: dataframe
    :param frm: Starting format of TCR data ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
    :type frm: str
    :param to: Format to convert TCR data to ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
    :type to: str
    :param frm_cols: List of current V, D and J gene column names
    :type frm_cols: list, optional
    :param species: Name of data folder with desired lookup tables
    :type species: str, optional
    :return: Converted TCR data
    :rtype: dataframe
    '''

    # Determine if we're using the 10X lookup table
    if frm == 'tenx':
        lookup_f = files('tcrconvert') / 'data' / species / 'lookup_from_tenx.csv'
        print("Warning: Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.")
    else:
        lookup_f = files('tcrconvert') / 'data' / species / 'lookup.csv'

    # Load lookup table
    try:
        lookup = pd.read_csv(lookup_f)
    except FileNotFoundError:
        print('Lookup table not found, please download IMGT reference FASTAs and run build_lookup_from_fastas()')

    # Figure out columns to use
    if frm == 'imgt' and not frm_cols:
        print('No column names provided for IMGT data, will use 10X column names:\n' 
                      + str(col_ref['tenx']))
        cols_from = 'tenx'
    else:
        cols_from = frm

    if frm_cols:
        v_from = frm_cols[0]
        d_from = frm_cols[1]
        j_from = frm_cols[2]
    else:
        v_from = col_ref[cols_from][0]
        d_from = col_ref[cols_from][1]
        j_from = col_ref[cols_from][2]

    # Convert V genes
    v_genes = df[[v_from]].merge(lookup[[frm, to]], how='left',
                                 left_on=v_from, right_on=frm).\
        drop(columns=frm)
    v_genes_bad = v_genes[v_genes[to].isna()][v_from].tolist()

    # Convert D genes
    d_genes = df[[d_from]].merge(lookup[[frm, to]], how='left',
                                 left_on=d_from, right_on=frm).\
        drop(columns=frm)
    d_genes_bad = d_genes[d_genes[to].isna()][d_from].tolist()

    # Convert J genes
    j_genes = df[[j_from]].merge(lookup[[frm, to]], how='left',
                                 left_on=j_from, right_on=frm).\
        drop(columns=frm)
    j_genes_bad = j_genes[j_genes[to].isna()][j_from].tolist()

    # Display genes we couldn't convert
    # TODO: have option to keep mangled gene names
    if len(v_genes_bad + d_genes_bad + j_genes_bad) > 0:
        print('Warning: These genes are not in the IMGT reference and have replaced with NA. '
          'Please fix the gene names and re-run:\n',
          v_genes_bad + d_genes_bad + j_genes_bad)

    # Swap out data in original dataframe
    out_df = df.copy()
    out_df[v_from] = v_genes[to].values
    out_df[d_from] = d_genes[to].values
    out_df[j_from] = j_genes[to].values

    return out_df
