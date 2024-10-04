import pandas as pd
import sys
from importlib.resources import files
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Standard column names for different sources of TCR data
col_ref = {'adaptive': ('v_resolved', 'd_resolved', 'j_resolved'),
           'adaptivev2': ('vMaxResolved', 'dMaxResolved', 'jMaxResolved'),
           'imgt': ('v_gene', 'd_gene', 'j_gene', 'c_gene'),
           'tenx': ('v_gene', 'd_gene', 'j_gene', 'c_gene')}


def convert_gene(df, frm, to, species='human', frm_cols=[]):
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

    # Check that required input is ok
    if frm == to:
        logger.error('"frm" and "to" formats should be different.')
        sys.exit(1)

    if df.empty:
        logger.error('Input dataframe is empty.')
        sys.exit(1)

    # Determine which lookup table to use
    if frm == 'tenx':
        lookup_f = files('tcrconvert') / 'data' / species / 'lookup_from_tenx.csv'
        logger.warning('Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.')
    elif frm == 'adaptive' or frm == 'adaptivev2':
        lookup_f = files('tcrconvert') / 'data' / species / 'lookup_from_adaptive.csv'
        if to == 'imgt':
            logger.warning('Converting from Adaptive to IMGT. If a gene lacks allele, will choose *01 as allele.')
    else:
        lookup_f = files('tcrconvert') / 'data' / species / 'lookup.csv'

    # Load lookup table
    try:
        lookup = pd.read_csv(lookup_f)
    except FileNotFoundError:
        logger.error('Lookup table not found, please download IMGT reference FASTAs and run build_lookup_from_fastas()')
        sys.exit(1)

    # Warn about no Adaptive C genes
    if to == 'adaptive' or to == 'adaptivev2':
        logger.warning('Adaptive only captures VDJ genes, any C genes will become NA.')

    # Figure out columns to use
    if frm == 'imgt' and not frm_cols:
        logger.info('No column names provided for IMGT data, will assume 10X column names: %s',
                    str(col_ref['tenx']))
        cols_from = 'tenx'
    if frm_cols:
        missing_cols = set(frm_cols) - set(df.columns)
        if missing_cols:
            logger.error('These columns are not in the input dataframe: %s', str(missing_cols))
            sys.exit(1)
        else:
            logger.info('Using these custom column names: %s',
                  str(frm_cols))
    else:
        cols_from = col_ref[frm]

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
    # TODO: still warn but perhaps keep original gene names instead of replace with NA?
    if bad_genes:
        logger.warning('These genes are not in IMGT and will be replaced with NA:\n %s',
                str(set(bad_genes))
        )

    # Swap out data in original dataframe
    out_df = df.copy()
    for col in new_genes:
        out_df[col] = new_genes[col][to].values
    
    # Replace NoData and np.nan with pd.NA
    out_df = out_df.replace('NoData', pd.NA).fillna(pd.NA)
    out_df = out_df

    return out_df
