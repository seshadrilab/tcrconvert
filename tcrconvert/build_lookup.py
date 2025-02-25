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
    '''Extract gene names from a reference FASTA.

    :param infile: Path to FASTA file
    :type infile: str
    :return: Gene names
    :rtype: list of str

    :Example:

    Given a FASTA file containing this header:

    .. code-block:: text

       \b
       >SomeText|TRAV1*01|MoreText|
       >SomeText|TRAV14/DV4*01|MoreText|
       >SomeText|TRAV38-2/DV8*01|MoreText|

    >>> import tcrconvert
    >>> fasta = tcrconvert.get_example_path('fasta_dir/test_trav.fa')
    >>> tcrconvert.build_lookup.parse_imgt_fasta(fasta)
    ['TRAV1*01', 'TRAV14/DV4*01', 'TRAV38-2/DV8*01']
    '''
    
    # Read the file and extract lines starting with ">"
    with open(infile, "r") as f:
        lines = f.readlines()

    # Extract the second element from lines starting with ">"
    imgt_list = []
    for line in lines:
        if line.startswith(">"):
            gene = line.split("|")[1]
            imgt_list.append(gene)

    return imgt_list


def extract_imgt_genes(data_dir):
    '''Extract gene names from all reference FASTA files in a folder.

    :param data_dir: Path to directory containing FASTA files
    :type data_dir: str
    :return: Gene names
    :rtype: DataFrame

    :Example:

    Given a folder with FASTA files containing these headers:

    .. code-block:: text

       \b
       >SomeText|TRAV1*01|MoreText|
       >SomeText|TRAV14/DV4*01|MoreText|
       >SomeText|TRAV38-2/DV8*01|MoreText|
       >SomeText|TRBV29/OR9-2*01|MoreText|
       >SomeText|TRBVA/OR9-2*01|MoreText|

    >>> import tcrconvert
    >>> fastadir = tcrconvert.get_example_path('fasta_dir') + '/'
    >>> tcrconvert.build_lookup.extract_imgt_genes(fastadir)
                  imgt
    0         TRAV1*01
    1    TRAV14/DV4*01
    2  TRAV38-2/DV8*01
    3  TRBV29/OR9-2*01
    4   TRBVA/OR9-2*01
    '''

    fastas = []
    for file in os.listdir(data_dir):
        if file.endswith('.fa') | file.endswith('.fasta'):
            fastas.append(data_dir + file)

    # Extract gene names and put into a dataframe
    imgt = []
    for fa in fastas:
        imgt = imgt + parse_imgt_fasta(fa)
    lookup = pd.DataFrame({'imgt': imgt})
    lookup_sorted = lookup.sort_values('imgt').reset_index(drop=True)

    return lookup_sorted


def add_dash_one(gene_str):
    '''Add a ``-01`` to genes without IMGT gene-level designation.

    :param gene_str: Gene name
    :type gene_str: str
    :return: Updated gene name
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.build_lookup.add_dash_one('TRBV2*01')
    'TRBV2-01*01'
    '''

    if '-' not in gene_str:
        # Add -1 before allele
        return gene_str.replace('*', '-01*')
    return gene_str


def pad_single_digit(gene_str):
    '''Add a zero to single-digit gene-level designatinon in gene names.

    :param gene_str: Gene name
    :type gene_str: str
    :return: Gene name
    :rtype: str

    :Example:

    >>> import tcrconvert
    >>> tcrconvert.build_lookup.pad_single_digit('TCRBV1-2')
    'TCRBV01-2'
    '''

    # Use regex to find a single digit preceded by letters and followed by a hyphen or asterisk
    updated_string = re.sub(r'([A-Za-z]+)(\d)([-\*])', r'\g<1>0\g<2>\g<3>', gene_str)
    return updated_string


def build_lookup_from_fastas(data_dir, species):
    '''Create these lookup tables from a directory of FASTA files:

    - lookup.csv
    - lookup_from_tenx.csv
    - lookup_from_adaptive.csv

    The lookup tables are stored in an application data folder via platformdirs. For example:

    - MacOS: ``~/Library/Application Support/tcrconvert/<species>``
    - Windows: ``C:/Documents and Settings/<User>/Application Data/Local Settings/Emma Bishop/tcrconvert/<species>``
    - Linux: ``~/.local/share/tcrconvert/<species>``

    :param data_dir: Directory containing FASTA files
    :type data_dir: str
    :param species: Name of the species, used to create a species-specific folder for storing lookup tables.
    :type species: str
    :return: None

    :Example:

    >>> import tcrconvert
    >>> fastadir = tcrconvert.get_example_path('fasta_dir') + '/'
    >>> tcrconvert.build_lookup.build_lookup_from_fastas(fastadir, 'rabbit')
    '''

    # Get the user data directory for saving lookup tables
    user_dir = platformdirs.user_data_dir("tcrconvert", "Emmma Bishop")
    save_dir = os.path.join(user_dir, species)
    os.makedirs(save_dir, exist_ok=True)

    # Extract IMGT gene names and put into a dataframe
    lookup = extract_imgt_genes(data_dir)

    # Create 10X column by removing allele info (e.g. *01) and slash from "/DV"
    lookup['tenx'] = lookup['imgt'].apply(lambda x: x[:-3].replace('/DV', 'DV'))

    # Create Adaptive columns by adding letters, 0's, removing /DV and renaming /OR
    lookup['adaptive'] = lookup['imgt'].apply(lambda x: x.\
                                            replace('TRAV14/DV4', 'TRAV14-1').\
                                            replace('TRAV23/DV6', 'TRAV23-1').\
                                            replace('TRAV29/DV5', 'TRAV29-1').\
                                            replace('TRAV36/DV7', 'TRAV36-1').\
                                            replace('TRAV38-2/DV8', 'TRAV38-2').\
                                            replace('TRAV4-4/DV10', 'TRAV4-4/').\
                                            replace('TRAV6-7/DV9', 'TRAV6-7').\
                                            replace('TRAV13-4/DV7', 'TRAV13-4').\
                                            replace('TRAV14D-3/DV8', 'TRAV14D-3').\
                                            replace('TRAV15D-1/DV6D-1', 'TRAV15D-1').\
                                            replace('TRAV15-1/DV6-1', 'TRAV15-1').\
                                            replace('TRAV16D/DV11', 'TRAV16D-1').\
                                            replace('TRAV21/DV12', 'TRAV21-1').\
                                            replace('TRAV15-2/DV6-2', 'TRAV15-2').\
                                            replace('TRAV15D-2/DV6D-2', 'TRAV15D-2').\
                                            replace('TR', 'TCR').\
                                            replace('-', '-0').\
                                            replace('/OR9-02', '-or09_02'))
    lookup['adaptive'] = lookup['adaptive'].apply(lambda x: add_dash_one(x))
    lookup['adaptive'] = lookup['adaptive'].apply(lambda x: pad_single_digit(x))
    lookup['adaptivev2'] = lookup['adaptive']

    # Set Adaptive columns to 'NoData' for constant genes (Adaptive only captures VDJ)
    # 'NoData' gets converted to NA by convert.convert_gene()
    lookup.loc[lookup['imgt'].str.contains('C'), ['adaptive', 'adaptivev2']] = 'NoData'

    # If converting from 10X will just need the first *01 allele
    from_tenx = lookup.groupby('tenx').first().reset_index()

    # Make table for Adaptive genes with or without allele
    lookup2 = lookup[~lookup.adaptivev2.str.contains('NoData')]
    from_adapt = lookup2[['adaptivev2', 'imgt', 'tenx']]
    from_adapt['adaptive'] = from_adapt['adaptivev2'].apply(lambda x: x[:-3])
    from_adapt = from_adapt.drop(columns='adaptivev2').groupby('adaptive').first().reset_index()
    from_adapt['adaptivev2'] = from_adapt['adaptive']
    from_adaptive = pd.concat([lookup2, from_adapt])[['adaptive', 'adaptivev2', 'imgt', 'tenx']]

    # Remove any duplicate rows
    lookup_path = os.path.join(save_dir, 'lookup.csv')
    from_tenx_path = os.path.join(save_dir, 'lookup_from_tenx.csv')
    from_adaptive_path = os.path.join(save_dir, 'lookup_from_adaptive.csv')

    # Save
    logger.info(f'Writing lookup tables to: {save_dir}')
    lookup.drop_duplicates().to_csv(lookup_path, index=False)
    from_tenx.drop_duplicates().to_csv(from_tenx_path, index=False)
    from_adaptive.drop_duplicates().to_csv(from_adaptive_path, index=False)


# Command-line version of build_lookup_from_fastas()
@click.command(name='build', no_args_is_help=True)
@click.option('-i', '--input', help='Path to folder of FASTA files', required=True, type=click.Path(exists=True))
@click.option('-s', '--species', help='Species name.', required=True)
def build_lookup_from_fastas_cli(input, species):
    '''Create lookup tables from within a folder of FASTA files.

    :Example:

    .. code-block:: bash

       $ tcrconvert build -i tcrconvert/examples/fasta_dir/ -s rabbit
    '''

    build_lookup_from_fastas(input, species)