import pandas as pd
import sys
from importlib.resources import files

# Standard column names for different sources of TCR data
# TODO: Generate 'resolved' column from 'gene' and 'allele' if missing in Adaptive data
col_ref = {'adaptive': ('v_resolved', 'd_resolved', 'j_resolved'),
           'adaptivev2': ('vMaxResolved', 'dMaxResolved', 'jMaxResolved'),
           'imgt': ('v_gene', 'd_gene', 'j_gene', 'c_gene'),
           'tenx': ('v_gene', 'd_gene', 'j_gene', 'c_gene')}


def convert_gene(df, frm, to, frm_cols=[], species='human'):
    '''Convert TCR V, D, J and C gene names from one naming convention to another.

    For AIRR use: frm_cols=['v_call', 'd_call', 'j_call', 'c_call']

    :param df: Dataframe of TCRs with gene names
    :type df: Pandas dataframe
    :param frm: Starting format of TCR data ['tenx', 'adaptive', 'adaptivev2', 'imgt']
    :type frm: str
    :param to: Format to convert TCR data to ['tenx', 'adaptive', 'adaptivev2', 'imgt']
    :type to: str
    :param frm_cols: List of custom gene column names. For AIRR use: ['v_call', 'd_call', 'j_call', 'c_call']
    :type frm_cols: list of str, optional
    :param species: Name of data folder with desired lookup tables
    :type species: str, optional
    :return: Converted TCR data
    :rtype: Pandas dataframe
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
        sys.exit()  # Quit

    # Warn about no Adaptive C genes
    if to == 'adaptive' or to == 'adaptivev2':
        print('Warning: Adaptive only captures VDJ genes, any C genes will become NA.')

    # Figure out columns to use
    if frm == 'imgt' and not frm_cols:
        print('No column names provided for IMGT data, will assume 10X column names:\n' 
                      + str(col_ref['tenx']))
        cols_from = 'tenx'
    if frm_cols:
        # Ensure at least one non-empty string to use if using custom columns
        cols_from = [s for s in frm_cols if s]
        if not cols_from:
            print('Please include at least one colunn name if using custom columns.')
        # Ensure not too many columns
        elif len(cols_from) > 4:
            print('Please only include colunn names for V, D, J, and/or C genes.')
        else:
            print('Warning: Using these custom column names:\n',
                  frm_cols)
    else:
        cols_from = col_ref[frm]

    # TODO: Ensure all columns are in input dataframe

    # Loop through gene name columns, getting converted genes
    new_genes = {}
    bad_genes = []

    for col in cols_from:
        if col in df.columns:
            new_genes[col] = df[[col]].\
                merge(lookup[[frm, to]], how='left', left_on=col, right_on=frm).\
                drop(columns=frm)
            bad_genes = bad_genes + new_genes[col][new_genes[col][to].isna()][col].dropna().tolist()
        else:
            continue

    # Display genes we couldn't convert
    # TODO: still warn but keep original gene names instead of replace with NA
    if bad_genes:
        print('Warning: These genes are not in IMGT and will be replaced with NA.\n',
                set(bad_genes)
        )

    # Swap out data in original dataframe
    out_df = df.copy()
    for col in new_genes:
        out_df[col] = new_genes[col][to].values
    
    # Replace NoData and np.nan with pd.NA
    out_df = out_df.replace('NoData', pd.NA).fillna(pd.NA)
    out_df = out_df

    return out_df
