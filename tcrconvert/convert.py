import pandas as pd
from importlib.resources import files

# Standard column names for different sources of TCR data
col_ref = {'tenx': ('v_gene', 'j_gene', 'cdr3'),
           'adaptive': ('v_resolved', 'j_resolved', 'cdr3_amino_acid'),
           'adaptive_v2': ('vMaxResolved', 'jMaxResolved', 'aminoAcid'),
           'imgt': ('v_gene', 'j_gene', 'cdr3')}


def convert_tcr(df,
                fmt_from,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                fmt_to,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                cols_use=None, extract_tcr=False, convert_cols=False):
    '''Convert TCR gene names from one naming convention to another.

    :param df: Dataframe of TCRs
    :type df: dataframe
    :param fmt_from: Starting format of TCR data
    :type fmt_from: str
    :param fmt_to: Format to convert TCR data to
    :type fmt_to: str
    :param cols_use: List or tuple of current V, J, and CDR3 column names
    :type cols_use: list or tuple, optional
    :param extract_tcr: Whether to return only the V, J, and CDR3 columns
    :type extract_tcr: bool, optional
    :param convert_cols: Whether to convert column names to match new format
    :type convert_cols: bool, optional
    :return: Converted TCR data
    :rtype: dataframe
    '''

    # Determine if we're using the 10X lookup table
    if fmt_from == 'tenx':
        lookup_f = files('tcrconvert.data').joinpath('lookup_from_tenx.csv')
        print("CONVERTING FROM 10X: CHOOSING *01 AS ALLELE FOR ALL GENES")
    else:
        lookup_f = files('tcrconvert.data').joinpath('lookup.csv')

    # Load lookup table
    try:
        lookup = pd.read_csv(lookup_f)
    except FileNotFoundError:
        print('Lookup table not found, please download IMGT reference FASTAs and run build_lookup.py')

    # Figure out columns to use
    if fmt_from == 'adaptive' and 'vMaxResolved' in df.columns:
        cols_from = 'adaptive_v2'
    elif fmt_from == 'imgt' and not cols_use:
        print('No column names provided for IMGT data, will use 10X column names:\n' 
                      + str(col_ref['tenx']))
        cols_from = 'tenx'
    else:
        cols_from = fmt_from

    # Column names to use for input and output
    # TODO: What if converting to imgt and want to update columns?
    if cols_use:
        v_from = cols_use[0]
        j_from = cols_use[1]
        cdr3_from = cols_use[2]
    else:
        v_from = col_ref[cols_from][0]
        j_from = col_ref[cols_from][1]
        cdr3_from = col_ref[cols_from][2]

    if convert_cols:
        v_to = col_ref[fmt_to][0]
        j_to = col_ref[fmt_to][1]
        cdr3_to = col_ref[fmt_to][2]
    else:
        v_to = v_from
        j_to = j_from
        cdr3_to = cdr3_from

    # Convert V and J genes
    v_genes = df[[v_from]].merge(lookup[[fmt_from, fmt_to]], how='left', 
                                 left_on=v_from, right_on=fmt_from).\
        drop(columns=fmt_from)
    v_genes_bad = v_genes[v_genes[fmt_to].isna()][v_from].tolist()

    j_genes = df[[j_from]].merge(lookup[[fmt_from, fmt_to]], how='left', 
                                 left_on=j_from, right_on=fmt_from).\
        drop(columns=fmt_from)
    j_genes_bad = j_genes[j_genes[fmt_to].isna()][j_from].tolist()

    # Display genes we couldn't convert
    # TODO: have option to keep mangled gene names
    if len(v_genes_bad + j_genes_bad) > 0:
        print('These genes are not in the IMGT reference and have replaced with NA. '
          'Please repair gene names manually and re-run:\n',
          v_genes_bad + j_genes_bad)

    # Consolidate converted data
    out = pd.DataFrame()
    out[v_to] = v_genes[fmt_to]
    out[j_to] = j_genes[fmt_to]
    out[cdr3_to] = df[cdr3_from]

    if extract_tcr:
        return out
    else:
        out_df = df.copy()
        out_df[v_from] = out[v_to].values
        out_df[j_from] = out[j_to].values
        out_df[cdr3_from] = out[cdr3_to].values
        out_df.rename(columns={v_from: v_to, 
                               j_from: j_to,
                               cdr3_from: cdr3_to}, inplace = True)
        return out_df
