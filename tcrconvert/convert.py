import pandas as pd
from importlib.resources import files
import logging
import click
import os
import platformdirs

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# Standard column names for different sources of TCR data
col_ref = {
    'adaptive': ['v_resolved', 'd_resolved', 'j_resolved'],
    'adaptivev2': ['vMaxResolved', 'dMaxResolved', 'jMaxResolved'],
    'imgt': ['v_gene', 'd_gene', 'j_gene', 'c_gene'],
    'tenx': ['v_gene', 'd_gene', 'j_gene', 'c_gene'],
}


def choose_lookup(frm, to, species='human'):
    """Determine which lookup table to use and get filepath.

    :param frm: Input format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type frm: str
    :param to: Output format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type to: str
    :param species: Species
    :type species: str, optional
    :return: Path to correct lookup table
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.convert.choose_lookup('imgt', 'adaptive')
    '.../tcrconvert/data/human/lookup.csv'
    """

    # Determine where to find lookup tables
    if species in ['human', 'mouse', 'rhesus']:  # Built-in
        lookup_dir = os.path.join(files('tcrconvert'), 'data')
    else:
        lookup_dir = platformdirs.user_data_dir('tcrconvert', 'Emmma Bishop')

    species_dir = os.path.join(lookup_dir, species)

    # Determine which lookup table to use
    if frm == 'tenx':
        lookup_f = os.path.join(species_dir, 'lookup_from_tenx.csv')
        logger.info(
            'Converting from 10X which lacks allele info. Choosing *01 as allele for all genes.'
        )
    elif frm == 'adaptive' or frm == 'adaptivev2':
        lookup_f = os.path.join(species_dir, 'lookup_from_adaptive.csv')
        if to == 'imgt':
            logger.info(
                'Converting from Adaptive to IMGT. If a gene lacks allele, will choose *01 as allele.'
            )
    else:
        lookup_f = os.path.join(species_dir, 'lookup.csv')

    # Check path
    if os.path.exists(lookup_f):
        return lookup_f
    else:
        logger.error(
            'Lookup table not found, please download IMGT reference FASTAs and run build_lookup_from_fastas()'
        )
        raise (FileNotFoundError)


def which_frm_cols(df, frm, frm_cols=[]):
    """Determine input columns to use for converting gene names.

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
    >>> tcr_file = tcrconvert.get_example_path('tenx.csv')
    >>> df = pd.read_csv(tcr_file)
    >>> tcrconvert.convert.which_frm_cols(df, 'tenx')
    ['v_gene', 'd_gene', 'j_gene', 'c_gene']
    """

    if frm == 'imgt' and not frm_cols:
        cols_from = col_ref['tenx']
        logger.info(
            f'No column names provided for IMGT data, will assume 10X column names: {str(cols_from)}'
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
            logger.info(f'Using these custom column names: {str(cols_from)}')
    else:
        cols_from = col_ref[frm]

    return cols_from


def convert_gene(df, frm, to, species='human', frm_cols=[], quiet=False):
    """Convert T-cell receptor V, D, J, and/or C gene names from one naming convention to another.

    :param df: Dataframe containing TCR gene names
    :type df: DataFrame
    :param frm: Input format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type frm: str
    :param to: Output format of TCR data ``['tenx', 'adaptive', 'adaptivev2', 'imgt']``
    :type to: str
    :param species: Species
    :type species: str, optional
    :param frm_cols: Custom V/D/J/C gene column names.
    :type frm_cols: list of str, optional
    :param quiet: Whether to suppress warning messages.
    :type quiet: bool, optional
    :return: Converted TCR data
    :rtype: DataFrame

    :Example:

    >>> import pandas as pd
    >>> import tcrconvert
    >>> original = pd.read_csv(tcrconvert.get_example_path('tenx.csv'))
    >>> original[['v_gene', 'd_gene', 'j_gene', 'c_gene', 'cdr3']]
        v_gene d_gene   j_gene c_gene             cdr3
    0  TRAV1-2  TRBD1   TRAJ12   TRAC     CAVMDSSYKLIF
    1  TRBV6-1  TRBD2  TRBJ2-1  TRBC2  CASSGLAGGYNEQFF
    2  TRBV6-4  TRBD2  TRBJ2-3  TRBC2  CASSGVAGGTDTQYF
    3  TRAV1-2  TRBD1   TRAJ33   TRAC     CAVKDSNYQLIW
    4    TRBV2  TRBD1  TRBJ1-2  TRBC1    CASNQGLNYGYTF
    >>> converted = tcrconvert.convert_gene(original, 'tenx', 'adaptive', quiet=True)
    >>> converted[['v_gene', 'd_gene', 'j_gene', 'c_gene', 'cdr3']]
              v_gene         d_gene         j_gene c_gene             cdr3
    0  TCRAV01-02*01  TCRBD01-01*01  TCRAJ12-01*01   <NA>     CAVMDSSYKLIF
    1  TCRBV06-01*01  TCRBD02-01*01  TCRBJ02-01*01   <NA>  CASSGLAGGYNEQFF
    2  TCRBV06-04*01  TCRBD02-01*01  TCRBJ02-03*01   <NA>  CASSGVAGGTDTQYF
    3  TCRAV01-02*01  TCRBD01-01*01  TCRAJ33-01*01   <NA>     CAVKDSNYQLIW
    4  TCRBV02-01*01  TCRBD01-01*01  TCRBJ01-02*01   <NA>    CASNQGLNYGYTF
    """

    # Set logging level
    if quiet:
        logger.setLevel(logging.ERROR)
    else:
        logger.setLevel(logging.INFO)

    # Check that input is ok
    if frm == to:
        logger.error('"frm" and "to" formats should be different.')
        raise (ValueError)
    if df.empty:
        logger.error('Input data is empty.')
        raise (ValueError)

    # Warn about no Adaptive C genes if needed
    if to == 'adaptive' or to == 'adaptivev2':
        logger.info(
            'Adaptive only captures VDJ genes. Converted C genes will become NA.'
        )

    # Load lookup table
    lookup_f = choose_lookup(frm, to, species)
    lookup = pd.read_csv(lookup_f)

    # Determine columns to use
    cols_from = which_frm_cols(df, frm, frm_cols)

    # Loop over gene columns, doing a pandas merge to get converted gene names
    new_genes = {}
    bad_genes = []

    for col in cols_from:
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
    out_df = out_df.replace('NoData', pd.NA).fillna(pd.NA)

    return out_df


# Command-line version of convert_gene()
@click.command(name='convert', no_args_is_help=True)
@click.option(
    '-i',
    '--infile',
    help='Path to input CSV or TSV',
    required=True,
    type=click.Path(exists=True),
)
@click.option('-o', '--outfile', help='Path to output CSV or TSV', required=True)
@click.option(
    '-f',
    '--frm',
    help='Input format of TCR data',
    required=True,
    type=click.Choice(['tenx', 'adaptive', 'adaptivev2', 'imgt'], case_sensitive=False),
)
@click.option(
    '-t',
    '--to',
    help='Output format of TCR data',
    required=True,
    type=click.Choice(['tenx', 'adaptive', 'adaptivev2', 'imgt'], case_sensitive=False),
)
@click.option(
    '-s', '--species', default='human', help='Species name.', show_default=True
)
@click.option(
    '-c',
    '--frm_cols',
    default=[],
    help='List of custom V/D/J/C gene column names.',
    show_default=True,
    multiple=True,
)
@click.option(
    '-q',
    '--quiet',
    is_flag=True,
    default=False,
    help='Whether to suppress warning messages.',
    show_default=True,
)
def convert_gene_cli(infile, outfile, frm, to, species, frm_cols, quiet):
    """Convert T-cell receptor V/D/J/C gene names.

    :Example:

    Using custom input columns 'myVgene', 'myDgene', 'myJgene'.

    .. code-block:: bash

       \b
       $ tcrconvert convert --infile tcrconvert/examples/customcols.csv \\
           --outfile tcrconvert/examples/custom2adapt.tsv \\
           --frm tenx \\
           --to adaptive \\
           -c myVgene \\
           -c myDgene \\
           -c myJgene
    """

    # Check that input and output paths are CSV/TSV
    if not infile.endswith(('csv', 'tsv')):
        raise click.BadParameter('"infile" must be a .csv or .tsv file')

    if not outfile.endswith(('csv', 'tsv')):
        raise click.BadParameter('"outfile" must be a .csv or .tsv file')

    # Load data
    # For our purposes, read in every column as string so that boolean values
    # don't get converted from uppercase to capitalized, etc.
    if not quiet:
        click.echo(f'Reading input file {infile}')
    sep_in = ',' if infile.endswith('csv') else '\t'
    df = pd.read_csv(infile, sep=sep_in, dtype=str)

    # Convert gene names
    # Cast frm_cols as list because will be read in from command line as tuple
    if not quiet:
        click.echo(f'Converting gene names from {frm} to {to}')
    out_df = convert_gene(df, frm, to, species, list(frm_cols), quiet)

    # Save output
    if not quiet:
        click.echo(f'Writing output to {outfile}')
    sep_out = ',' if outfile.endswith('csv') else '\t'
    out_df.to_csv(outfile, sep=sep_out, index=False)
