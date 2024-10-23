import pandas as pd
from importlib.resources import files
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Standard column names for different sources of TCR data
col_ref = {'adaptive': ['v_resolved', 'd_resolved', 'j_resolved'],
           'adaptivev2': ['vMaxResolved', 'dMaxResolved', 'jMaxResolved'],
           'imgt': ['v_gene', 'd_gene', 'j_gene', 'c_gene'],
           'tenx': ['v_gene', 'd_gene', 'j_gene', 'c_gene']}


def choose_lookup(frm, to, species='human'):
    '''Determine which lookup table to use.

    :param frm: Input format of TCR data 
    :type frm: str
    :param to: Output format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type to: str
    :param species: Species folder name under ``tcrconvert/data/``
    :type species: str, optional
    :return: Path to correct lookup table
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.choose_lookup('tenx', 'imgt')
    '/path/to/data/human/lookup_from_tenx.csv'
    '''

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

    # Check path
    if lookup_f.is_file():
        return lookup_f
    else:
        logger.error('Lookup table not found, please download IMGT reference FASTAs and run build_lookup_from_fastas()')
        raise(FileNotFoundError)


def which_frm_cols(df, frm, frm_cols=[]):
    '''Determine input columns to use for converting gene names.

    :param df: Dataframe containing TCR gene names
    :type df: DataFrame
    :param frm: Input format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type frm: str
    :param frm_cols: Custom column names to use
    :type frm_cols: list of str, optional
    :return: Column names to use
    :rtype: list of str

    :Example:

    >>> import pandas as pd
    >>> import tcrconvert
    >>> df = pd.read_csv('path/to/10x_tcrs.csv')
    >>> tcrconvert.which_frm_cols(df, 'tenx')
    ['v_gene', 'd_gene', 'j_gene', 'c_gene']
    '''

    if frm == 'imgt' and not frm_cols:
        cols_from = col_ref['tenx']
        logger.info('No column names provided for IMGT data, will assume 10X column names: %s',
                    str(cols_from))
    if frm_cols:
        missing_cols = set(frm_cols) - set(df.columns)
        if missing_cols:
            logger.error('These columns are not in the input dataframe: %s', str(missing_cols))
            raise(ValueError)
        else:
            cols_from = frm_cols
            logger.info('Using these custom column names: %s', str(cols_from))
    else:
        cols_from = col_ref[frm]
    
    return cols_from


def convert_gene(df, frm, to, species='human', frm_cols=[], quiet=False):
    '''Convert T-cell receptor V, D, J, and/or C gene names from one naming convention to another.

    :param df: Dataframe containing TCR gene names
    :type df: DataFrame
    :param frm: Input format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type frm: str
    :param to: Output format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type to: str
    :param frm_cols: Custom V/D/J/C gene column names.
    :type frm_cols: list of str, optional
    :param species: Species folder name under ``tcrconvert/data/``.
    :type species: str, optional
    :param quiet: Whether to suppress warning messages.
    :type quiet: bool, optional
    :return: Converted TCR data
    :rtype: DataFrame

    :Example:

    >>> import pandas as pd
    >>> import tcrconvert
    >>> tenx_df = pd.read_csv('path/to/10x_tcrs.csv')
    >>> tenx_df
            v_gene  d_gene   j_gene  c_gene    cdr3
    0     TRAV12-1    <NA>   TRAJ16    TRAC  CAVLIF
    1       TRBV15   TRBD1  TRBJ2-5   TRBC2  CASSGF
    >>> tcrconvert.convert_gene(tenx_df, 'tenx', 'adaptive')
              v_gene         d_gene         j_gene  c_gene    cdr3
    0  TCRAV12-01*01           <NA>  TCRAJ16-01*01    <NA>  CAVLIF
    1  TCRBV15-01*01  TCRBD01-01*01  TCRBJ02-05*01    <NA>  CASSGF
    '''

    # Set logging level
    if quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.WARNING)

    # Check that required input is ok
    if frm == to:
        logger.error('"frm" and "to" formats should be different.')
        raise(ValueError)
    if df.empty:
        logger.error('Input dataframe is empty.')
        raise(ValueError)

    # Warn about no Adaptive C genes if needed
    if to == 'adaptive' or to == 'adaptivev2':
        logger.warning('Adaptive only captures VDJ genes, any C genes will become NA.')

    # Load lookup table
    lookup_f = choose_lookup(frm, to, species)
    lookup = pd.read_csv(lookup_f)

    # Figure out columns to use
    cols_from = which_frm_cols(df, frm, frm_cols)

    # Loop over gene columns, doing a pandas merge to get converted gene names
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
    if bad_genes:
        logger.warning('These genes are not in IMGT and will be replaced with NA:\n %s',
                str(set(bad_genes)))

    # Swap out data in original dataframe
    out_df = df.copy()
    for col in new_genes:
        out_df[col] = new_genes[col][to].values

    # Replace NoData and np.nan with pd.NA
    out_df = out_df.replace('NoData', pd.NA).fillna(pd.NA)

    return out_df
