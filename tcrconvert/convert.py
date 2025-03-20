import pandas as pd
from importlib.resources import files
import logging
import click
import os
import platformdirs

# Set up logging
logging.basicConfig(level=logging.WARNING, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Standard column names for different sources of TCR data
col_ref = {
    'adaptive': ['v_resolved', 'd_resolved', 'j_resolved'],
    'adaptivev2': ['vMaxResolved', 'dMaxResolved', 'jMaxResolved'],
    'imgt': ['v_gene', 'd_gene', 'j_gene', 'c_gene'],
    'tenx': ['v_gene', 'd_gene', 'j_gene', 'c_gene'],
}


def choose_lookup(frm, to, species='human', verbose=True):
    """Choose lookup table

    Determine which CSV lookup table to use based on the the input format
    (``frm``) and returns the path to that file.

    :param frm: Input format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type frm: str
    :param to: Output format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type to: str
    :param species: Species, defaults to ``'human'``
    :type species: str, optional
    :param verbose: Whether to show all messages, defaults to ``True``
    :type verbose: bool, optional
    :return: Path to correct lookup table
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.convert.choose_lookup('imgt', 'adaptive', verbose=False) # doctest: +ELLIPSIS
    '.../tcrconvert/data/human/lookup.csv'
    """

    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    # Determine where to find lookup tables
    if species in ['human', 'mouse', 'rhesus']:
        lookup_dir = os.path.join(files('tcrconvert'), 'data')
    else:
        lookup_dir = platformdirs.user_data_dir('tcrconvert', 'Emmma Bishop')

    data_path = os.path.join(lookup_dir, species)

    if frm == 'tenx':
        lookup_f = os.path.join(data_path, 'lookup_from_tenx.csv')
        logger.info('Converting from 10X. Using *01 as allele for all genes.')
    elif frm == 'adaptive' or frm == 'adaptivev2':
        lookup_f = os.path.join(data_path, 'lookup_from_adaptive.csv')
        if to == 'imgt':
            logger.info(
                'Converting from Adaptive to IMGT. Using *01 for genes lacking alleles.'
            )
    else:
        lookup_f = os.path.join(data_path, 'lookup.csv')

    if os.path.exists(lookup_f):
        return lookup_f
    else:
        logger.error('Lookup table not found, please run build_lookup_from_fastas().')
        raise (FileNotFoundError)


def which_frm_cols(df, frm, frm_cols=[], verbose=True):
    """Determine input columns to use

    Determine the columns that are expected to hold gene name information in
    the input file based on the input format (``frm``). Returns a list of
    those column names.

    :param df: Dataframe containing TCR gene names
    :type df: DataFrame
    :param frm: Input format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type frm: str
    :param frm_cols: Custom column names to use
    :type frm_cols: list of str, optional
    :param verbose: Whether to show all messages, defaults to ``True``
    :type verbose: bool, optional
    :return: Column names to use
    :rtype: list of str

    :Example:

    >>> import pandas as pd
    >>> import tcrconvert
    >>> tcr_file = tcrconvert.get_example_path('tenx.csv')
    >>> df = pd.read_csv(tcr_file)
    >>> tcrconvert.convert.which_frm_cols(df, 'tenx')
    ['v_gene', 'd_gene', 'j_gene', 'c_gene']
    """

    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    if frm == 'imgt' and not frm_cols:
        cols_from = col_ref['tenx']
        logger.warning(
            f'No column names for IMGT data. Using 10X columns: {str(cols_from)}'
        )
    if frm_cols:
        missing_cols = set(frm_cols) - set(df.columns)
        if missing_cols:
            logger.error(
                f'These columns are not in the input dataframe: {str(missing_cols)}'
            )
            raise (ValueError)
        else:
            cols_from = frm_cols
            logger.info(f'Using custom column names: {str(cols_from)}')
    else:
        cols_from = col_ref[frm]

    return cols_from


def convert_gene(df, frm, to, species='human', frm_cols=[], verbose=True):
    """Convert gene names

    Convert T-cell receptor (TCR) gene names between the IMGT, 10X, and Adaptive
    formats. Determines the columns to convert based on the input format
    (``frm``) unless specified by the user (``frm_cols``). Returns a modified
    version of the input data frame with converted gene names while preserving
    row order.

    Behavioral Notes:

    - If a gene name cannot be mapped, it is replaced with ``NaN``, and a warning is issued.
    - If ``frm`` is ``'imgt'`` and ``frm_cols`` is not provided, 10X column names are assumed.
    - Constant (C) genes are set to ``NaN`` when converting to Adaptive formats, as Adaptive does not capture constant regions.
    - The input does not need to include all gene types; partial inputs (e.g., only V genes) are supported.
    - If no values in a custom column can be mapped (e.g., a CDR3 column) it is skipped and a warning is raised.

    Standard Column Names:

    If ``frm_cols`` is not provided, these column names will be used if present:

    - **IMGT**: ``'v_gene'``, ``'d_gene'``, ``'j_gene'``, ``'c_gene'``
    - **10X**: ``'v_gene'``, ``'d_gene'``, ``'j_gene'``, ``'c_gene'``
    - **Adaptive**: ``'v_resolved'``, ``'d_resolved'``, ``'j_resolved'``
    - **Adaptive v2**: ``'vMaxResolved'``, ``'dMaxResolved'``, ``'jMaxResolved'``

    :param df: Dataframe containing TCR gene names
    :type df: DataFrame
    :param frm: Input format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type frm: str
    :param to: Output format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type to: str
    :param species: Species name. Defaults to ``'human'``.
    :type species: str, optional
    :param frm_cols: Custom gene column names.
    :type frm_cols: list of str, optional
    :param verbose: Whether to show all messages. Defaults to ``True``.
    :type verbose: bool, optional
    :return: Converted TCR data
    :rtype: DataFrame

    :Example:

    >>> import pandas as pd
    >>> import tcrconvert
    >>> tcr_file = tcrconvert.get_example_path('tenx.csv')
    >>> df = pd.read_csv(tcr_file)[['barcode', 'v_gene', 'j_gene', 'cdr3']]
    >>> df
                  barcode        v_gene   j_gene             cdr3
    0  AAACCTGAGACCACGA-1    TRAV29/DV5   TRAJ12     CAVMDSSYKLIF
    1  AAACCTGAGACCACGA-1  TRBV20/OR9-2  TRBJ2-1  CASSGLAGGYNEQFF
    2  AAACCTGAGGCTCTTA-1         TRDV2    TRDJ3  CASSGVAGGTDTQYF
    3  AAACCTGAGGCTCTTA-1         TRGV9    TRGJ1     CAVKDSNYQLIW
    >>> tcrconvert.convert_gene(df, 'tenx', 'adaptive', verbose=False)
                  barcode              v_gene         j_gene             cdr3
    0  AAACCTGAGACCACGA-1       TCRAV29-01*01  TCRAJ12-01*01     CAVMDSSYKLIF
    1  AAACCTGAGACCACGA-1  TCRBV20-or09_02*01  TCRBJ02-01*01  CASSGLAGGYNEQFF
    2  AAACCTGAGGCTCTTA-1       TCRDV02-01*01  TCRDJ03-01*01  CASSGVAGGTDTQYF
    3  AAACCTGAGGCTCTTA-1       TCRGV09-01*01  TCRGJ01-01*01     CAVKDSNYQLIW
    """

    if verbose:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)

    if frm == to:
        logger.error('"frm" and "to" formats should be different.')
        raise (ValueError)
    if not isinstance(df, pd.DataFrame):
        logger.error('Input is not a pandas DataFrame.')
        raise (TypeError)
    if df.empty:
        logger.error('Input data is empty.')
        raise (ValueError)
    if to == 'adaptive' or to == 'adaptivev2':
        logger.warning('Adaptive only captures VDJ genes; C genes will be NA.')

    # Load lookup table and determine input columns
    lookup_f = choose_lookup(frm, to, species, verbose)
    lookup = pd.read_csv(lookup_f)
    cols_from = which_frm_cols(df, frm, frm_cols, verbose)

    # Loop over gene columns, doing a merge to get converted gene names
    new_genes = {}
    bad_genes = []

    for col in cols_from:
        if col in df.columns:
            good_genes = (
                df[[col]]
                .merge(lookup[[frm, to]], how='left', left_on=col, right_on=frm)
                .drop(columns=frm)
            )
            # Note genes where the merge produced an NA on the 'to' format side
            new_bad_genes = good_genes[good_genes[to].isna()][col].dropna().tolist()
            # We don't expect the entire column of genes to be empty.
            if len(new_bad_genes) < len(good_genes):
                new_genes[col] = good_genes
                bad_genes += new_bad_genes
            else:
                logger.warning(
                    f"The input column '{col}' doesn't contain any valid genes and was skipped."
                )
                continue

    # Display genes we couldn't convert
    if bad_genes:
        sorted_list = sorted(list(set(bad_genes)))
        logger.warning(
            f'These genes are not in IMGT for this species and will be replaced with NA:\n {str(sorted_list)}'
        )

    # Swap out data in original dataframe
    out_df = df.copy()
    for col in new_genes:
        out_df[col] = new_genes[col][to].values
        # Replace NoData and np.nan with pd.NA
        out_df[col] = out_df[col].replace('NoData', pd.NA)

    return out_df


# Command-line version of convert_gene()
@click.command(name='convert', no_args_is_help=True)
@click.option(
    '-i',
    '--input',
    help='Input file (CSV or TSV)',
    required=True,
    type=click.Path(exists=True),
)
@click.option('-o', '--output', help='Output file (CSV or TSV)', required=True)
@click.option(
    '-f',
    '--frm',
    help='Input TCR gene format',
    required=True,
    type=click.Choice(['tenx', 'adaptive', 'adaptivev2', 'imgt'], case_sensitive=False),
)
@click.option(
    '-t',
    '--to',
    help='Output TCR gene format',
    required=True,
    type=click.Choice(['tenx', 'adaptive', 'adaptivev2', 'imgt'], case_sensitive=False),
)
@click.option(
    '-s', '--species', default='human', help='Species name', show_default=True
)
@click.option(
    '-c',
    '--column',
    default=[],
    help='Custom gene column name',
    show_default=True,
    multiple=True,
)
@click.option(
    '-v',
    '--verbose',
    default=True,
    help='Show INFO-level messages',
    show_default=True,
)
def convert_gene_cli(input, output, frm, to, species, column, verbose):
    """Convert T-cell receptor V/D/J/C gene names.

    :Example:

    Using custom input columns 'myVgene', 'myDgene', 'myJgene'.

    .. code-block:: bash

       \b
       $ tcrconvert convert \\
           --input tcrconvert/examples/customcols.csv \\
           --output tcrconvert/examples/custom2adapt.tsv \\
           --frm tenx \\
           --to adaptive \\
           -c myVgene \\
           -c myDgene \\
           -c myJgene
    """

    # Check that input and output paths are CSV/TSV
    if not input.endswith(('csv', 'tsv')):
        raise click.BadParameter('"input" must be a .csv or .tsv file')

    if not output.endswith(('csv', 'tsv')):
        raise click.BadParameter('"output" must be a .csv or .tsv file')

    # Load data
    # For our purposes, read in every column as string so that boolean values
    # don't get converted from uppercase to capitalized, etc.
    if verbose:
        click.echo(f'Reading input TCR data from: {os.path.abspath(input)}')
    sep_in = ',' if input.endswith('csv') else '\t'
    df = pd.read_csv(input, sep=sep_in, dtype=str)

    # Convert gene names
    # Cast frm_cols as list because will be read in from command line as tuple
    if verbose:
        click.echo(f'Converting gene nomenclature from "{frm}" to "{to}"')
    out_df = convert_gene(df, frm, to, species, list(column), verbose)

    # Save output
    if verbose:
        click.echo(
            f'Writing TCR data with converted gene names to: {os.path.abspath(output)}'
        )
    sep_out = ',' if output.endswith('csv') else '\t'
    out_df.to_csv(output, sep=sep_out, index=False)
