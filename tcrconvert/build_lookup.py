import os
import re
import pandas as pd
import click
import platformdirs
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


def parse_imgt_fasta(infile):
    """Extract gene names from a reference FASTA

    Extracts the second element from a "|"-delimited FASTA header, which will
    be the gene name for IMGT reference FASTAs.

    :param infile: Path to FASTA file
    :type infile: str
    :return: Gene names
    :rtype: list of str

    :Example:

    Given a FASTA file containing this header:

    .. code-block:: text

       \b
       >SomeText|TRBV29-1*01|MoreText|
       >SomeText|TRBV29-1*02|MoreText|
       >SomeText|TRBV29/OR9-2*01|MoreText|

    >>> import tcrconvert
    >>> fasta = tcrconvert.get_example_path('fasta_dir/test_trbv.fa')
    >>> tcrconvert.build_lookup.parse_imgt_fasta(fasta)
    ['TRBV29-1*01', 'TRBV29-1*02', 'TRBV29/OR9-2*01']
    """

    with open(infile, 'r') as f:
        lines = f.readlines()

    # Extract the second element from lines starting with ">"
    imgt_list = []
    for line in lines:
        if line.startswith('>'):
            gene = line.split('|')[1]
            imgt_list.append(gene)

    return imgt_list


def extract_imgt_genes(data_dir):
    """Extract all gene names from a folder of FASTAs

    First run ``parse_imgt_fasta()`` on all FASTA files in a given folder to
    pull out the gene names. Then return those names in an alphabetically
    sorted dataframe.

    :param data_dir: Path to directory containing FASTA files
    :type data_dir: str
    :return: Gene names
    :rtype: DataFrame

    :Example:

    Given a folder with FASTA files containing these headers:

    .. code-block:: text

       \b
        >SomeText|TRAC*01|MoreText|
        >SomeText|TRAV1-1*01|MoreText|
        >SomeText|TRAV1-1*02|MoreText|
        >SomeText|TRAV1-2*01|MoreText|
        >SomeText|TRAV14/DV4*01|MoreText|
        >SomeText|TRAV38-1*01|MoreText|
        >SomeText|TRAV38-2/DV8*01|MoreText|
        >SomeText|TRBV29-1*01|MoreText|
        >SomeText|TRBV29-1*02|MoreText|
        >SomeText|TRBV29/OR9-2*01|MoreText|

    >>> import tcrconvert
    >>> fastadir = tcrconvert.get_example_path('fasta_dir')
    >>> tcrconvert.build_lookup.extract_imgt_genes(fastadir)
                  imgt
    0          TRAC*01
    1       TRAV1-1*01
    2       TRAV1-1*02
    3       TRAV1-2*01
    4    TRAV14/DV4*01
    5      TRAV38-1*01
    6  TRAV38-2/DV8*01
    7      TRBV29-1*01
    8      TRBV29-1*02
    9  TRBV29/OR9-2*01
    """

    fastas = []
    for file in os.listdir(data_dir):
        if file.endswith('.fa') | file.endswith('.fasta'):
            fastas.append(os.path.join(data_dir, file))
    imgt = []
    for fa in fastas:
        imgt = imgt + parse_imgt_fasta(fa)

    # Create and sort output data frame
    lookup = pd.DataFrame({'imgt': imgt})
    lookup_sorted = lookup.sort_values('imgt').reset_index(drop=True)

    return lookup_sorted


def add_dash_one(gene_str):
    """Add ``-01`` to gene names lacking gene-level info

    Some genes just have the IMGT subgroup (e.g. TRBV2) and allele (e.g. *01)
    designation. The Adaptive format always includes an IMGT gene (e.g. -01)
    designation, with "-01" as the apparent default. `add_dash_one()` adds a
    default gene-level designation if it's missing.

    :param gene_str: Gene name
    :type gene_str: str
    :return: Updated gene name
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.build_lookup.add_dash_one('TRBV2*01')
    'TRBV2-01*01'
    """

    if '-' not in gene_str:
        return gene_str.replace('*', '-01*')
    return gene_str


def pad_single_digit(gene_str):
    """Add a ``0`` to single-digit gene-level designation

    Take a gene name and ensure that any single-digit number following a
    sequence of letters is padded with a leading zero. This is to match the
    Adaptive format.

    :param gene_str: Gene name
    :type gene_str: str
    :return: Updated gene name
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.build_lookup.pad_single_digit('TCRBV1-2')
    'TCRBV01-2'
    """

    updated_string = re.sub(r'([A-Za-z]+)(\d)([-\*])', r'\g<1>0\g<2>\g<3>', gene_str)
    return updated_string


def save_lookup(df, savedir, name):
    """Save a lookup table to a CSV file

    Save a dataframe as a CSV file (without row names) in the specified directory.

    :param df: Dataframe containing the lookup table data
    :type df: pandas DataFrame
    :param savedir: Path to the save directory
    :type savedir: str
    :param name: File name (should end in `.csv`)
    :type name: str
    :return: None

    :Example:

    >>> import tcrconvert
    >>> import tempfile
    >>> dat = pd.read_csv(tcrconvert.get_example_path("fasta_dir/lookup.csv"))
    >>> save_dir = os.path.join(tempfile.gettempdir(), 'tcrconvert_tmp')
    >>> tcrconvert.build_lookup.save_lookup(dat, save_dir, "newlookup.csv")
    """

    # Ensure valid inputs
    if not os.path.exists(savedir):
        os.makedirs(savedir, exist_ok=True)

    if not isinstance(df, pd.DataFrame):
        raise ValueError("'df' must be a pandas DataFrame")

    file_path = os.path.join(savedir, name)
    df.to_csv(file_path, index=False)


def build_lookup_from_fastas(data_dir, species):
    """Create lookup tables

    Process IMGT reference FASTA files in a given folder to generate lookup
    tables used for making gene name conversions. It extracts all gene names
    and transforms them into 10X and Adaptive formats following predefined
    conversion rules. The resulting files are created:

    - ``lookup.csv``: IMGT gene names and their 10X and Adaptive equivalents.
    - ``lookup_from_tenx.csv``: Gene names aggregated by their 10X identifiers, with one representative allele (``*01``) for each.
    - ``lookup_from_adaptive.csv``: Adaptive gene names, with or without alleles, and their IMGT and 10X equivalents.

    The files are saved in a given subfolder (``species``) within the appropriate
    application folder via ``platformdirs``. For example:

    - MacOS: ``~/Library/Application Support/<AppName>``
    - Windows: ``C:\\Documents and Settings\\<User>\\Application Data\\Local Settings\\<AppAuthor>\\<AppName>``
    - Linux: ``~/.local/share/<AppName>``

    If a folder named ``species`` already exists in that location, it will be replaced.

    Key transformations from IMGT:

    - **10X:**
        - Remove allele information (e.g., ``*01``) and modify ``/DV`` occurrences.
    - **Adaptive:**
        - Apply renaming rules such as adding gene-level designations and zero-padding single-digit numbers.
        - Convert constant genes to ``'NoData'`` (Adaptive only captures VDJ) which will become ``NA`` after the merge in ``convert_gene()``.

    :param data_dir: Directory containing FASTA files
    :type data_dir: str
    :param species: Name of species that will be used when running TCRconvert with these lookup tables.
    :type species: str
    :return: Path to the new lookup directory
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> fastadir = tcrconvert.get_example_path('fasta_dir')
    >>> tcrconvert.build_lookup.build_lookup_from_fastas(fastadir, 'rabbit') # doctest: +ELLIPSIS
    '...tcrconvert/rabbit'
    """

    # Check that species can be a valid folder name
    forbidden_char = r'[/\\:*?\"<>|~`\n\t]'
    if re.search(forbidden_char, species):
        sanitized = re.sub(forbidden_char, '_', species)
        raise ValueError(
            f"Proposed folder name '{species}' contains invalid characters.\n"
            f'Suggestion: {sanitized}'
        )

    # Get the user data directory for saving lookup tables
    user_dir = platformdirs.user_data_dir('tcrconvert', 'Emmma Bishop')
    save_dir = os.path.join(user_dir, species)
    os.makedirs(save_dir, exist_ok=True)

    lookup = extract_imgt_genes(data_dir)

    # Create 10X column
    lookup['tenx'] = lookup['imgt'].apply(lambda x: x[:-3])
    lookup['tenx'] = lookup['tenx'].apply(
        lambda x: x.replace('TRAV13-4/DV7', 'TRAV13-4-DV7')
        .replace('TRAV14D-3/DV8', 'TRAV14D-3-DV8')
        .replace('TRAV15-1/DV6-1', 'TRAV15-1-DV6-1')
        .replace('TRAV15-2/DV6-2', 'TRAV15-2-DV6-2')
        .replace('TRAV15D-1/DV6D-1', 'TRAV15D-1-DV6D-1')
        .replace('TRAV15D-2/DV6D-2', 'TRAV15D-2-DV6D-2')
        .replace('TRAV16D/DV11', 'TRAV16D-DV11')
        .replace('TRAV21/DV12', 'TRAV21-DV12')
        .replace('TRAV4-4/DV10', 'TRAV4-4-DV10')
        .replace('TRAV6-7/DV9', 'TRAV6-7-DV9')
    )

    # Create Adaptive columns
    lookup['adaptive'] = lookup['imgt'].apply(
        lambda x: x.replace('TRAV14/DV4', 'TRAV14-1')
        .replace('TRAV23/DV6', 'TRAV23-1')
        .replace('TRAV29/DV5', 'TRAV29-1')
        .replace('TRAV36/DV7', 'TRAV36-1')
        .replace('TRAV38-2/DV8', 'TRAV38-2')
        .replace('TRAV4-4/DV10', 'TRAV4-4/')
        .replace('TRAV6-7/DV9', 'TRAV6-7')
        .replace('TRAV13-4/DV7', 'TRAV13-4')
        .replace('TRAV14D-3/DV8', 'TRAV14D-3')
        .replace('TRAV15D-1/DV6D-1', 'TRAV15D-1')
        .replace('TRAV15-1/DV6-1', 'TRAV15-1')
        .replace('TRAV16D/DV11', 'TRAV16D-1')
        .replace('TRAV21/DV12', 'TRAV21-1')
        .replace('TRAV15-2/DV6-2', 'TRAV15-2')
        .replace('TRAV15D-2/DV6D-2', 'TRAV15D-2')
        .replace('TR', 'TCR')
        .replace('-', '-0')
        .replace('/OR9-02', '-or09_02')
    )
    lookup['adaptive'] = lookup['adaptive'].apply(lambda x: add_dash_one(x))
    lookup['adaptive'] = lookup['adaptive'].apply(lambda x: pad_single_digit(x))
    lookup['adaptivev2'] = lookup['adaptive']
    lookup.loc[lookup['imgt'].str.contains('C'), ['adaptive', 'adaptivev2']] = 'NoData'

    # If converting from 10X will just need the *01 allele
    from_tenx = lookup.groupby('tenx').first().reset_index()

    # Start Adaptive tables
    lookup2 = lookup[~lookup.adaptive.str.contains('NoData')]

    # Adaptive: Gene-level info but not allele-level (e.g. TCRAJ03-01)
    adapt_no_allele = lookup2[['adaptivev2', 'imgt', 'tenx']]
    adapt_no_allele['adaptive'] = adapt_no_allele['adaptivev2'].apply(lambda x: x[:-3])
    adapt_no_allele = (
        adapt_no_allele.drop(columns='adaptivev2')
        .groupby('adaptive')
        .first()
        .reset_index()
    )

    # Start Adaptive: No gene-level info where unneeded, with and without allele-level (e.g. TCRAV14*01 and TCRAV14)
    subgroup_only = lookup2[['adaptivev2', 'imgt', 'tenx']]
    subgroup_only['tenx_prefix'] = subgroup_only['tenx'].str.split('-', expand=True)[0]

    # Group by 'tenx_prefix', keeping groups with only one unique 'tenx' value
    agg_data = subgroup_only.groupby('tenx_prefix')['tenx'].nunique().reset_index()
    agg_data['tenx'] = agg_data['tenx'] == 1
    agg_filtered = agg_data[agg_data['tenx']]

    # Make a DataFrame with subgroup-level Adaptive gene names
    merged_data = pd.merge(subgroup_only, agg_filtered, on='tenx_prefix')

    subgroup_only_rows_start = pd.DataFrame(
        {
            'adaptive': merged_data['adaptivev2'].str.replace(
                r'-\d+.*', '', regex=True
            ),
            'imgt': merged_data['imgt'],
            'tenx': merged_data['tenx_x'],
            'tenx_prefix': merged_data['tenx_prefix'],
        }
    )
    subgroup_only_rows = (
        subgroup_only_rows_start.groupby('adaptive').first().reset_index()
    )

    subgroup_only_with_allele_rows = pd.DataFrame(
        {
            'adaptive': merged_data['adaptivev2'].str.replace(
                r'-0[0-9]', '', regex=True
            ),
            'imgt': merged_data['imgt'],
            'tenx': merged_data['tenx_x'],
            'tenx_prefix': merged_data['tenx_prefix'],
        }
    )

    # Combine the new rows with the original data
    from_adaptive_updated = pd.concat(
        [adapt_no_allele, subgroup_only_rows, subgroup_only_with_allele_rows],
        ignore_index=True,
    )
    from_adaptive_updated = from_adaptive_updated.drop(columns=['tenx_prefix'])

    # Final polishing
    from_adaptive_updated['adaptivev2'] = from_adaptive_updated['adaptive']
    from_adaptive = pd.concat([lookup2, from_adaptive_updated], ignore_index=True)
    from_adaptive = from_adaptive[['adaptive', 'adaptivev2', 'imgt', 'tenx']]
    from_adaptive = from_adaptive.sort_values(by='adaptive').reset_index(drop=True)

    # Remove duplicate rows
    lookup.drop_duplicates(inplace=True)
    from_tenx.drop_duplicates(inplace=True)
    from_adaptive.drop_duplicates(inplace=True)

    # Save
    logger.info(f'Writing lookup tables to: {save_dir}')
    save_lookup(lookup, save_dir, 'lookup.csv')
    save_lookup(from_tenx, save_dir, 'lookup_from_tenx.csv')
    save_lookup(from_adaptive, save_dir, 'lookup_from_adaptive.csv')

    return save_dir


# Command-line version of build_lookup_from_fastas()
@click.command(name='build', no_args_is_help=True)
@click.option(
    '-i',
    '--input',
    help='Path to folder of FASTA files',
    required=True,
    type=click.Path(exists=True),
)
@click.option('-s', '--species', help='Species name.', required=True)
def build_lookup_from_fastas_cli(input, species):
    """Create lookup tables
    :Example:

    .. code-block:: bash

       $ tcrconvert build -i tcrconvert/examples/fasta_dir/ -s rabbit
    """

    file_path = build_lookup_from_fastas(input, species)
    click.echo(f'Lookup table written to: {file_path}')
