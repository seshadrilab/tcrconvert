
import os
import re
import pandas as pd

here = os.path.dirname(__file__)
data_dir = os.path.join(here, 'data/')

def main():
    # Extract IMGT gene names and put into a dataframe
    lookup = extract_imgt_genes(data_dir)

    # Create 10X column by removing allele info (e.g. *01)
    lookup['tenx'] = lookup['imgt'].apply(lambda x: x[:-3])

    # Create Adaptive columns by adding letters and 0's
    lookup['adaptive'] = lookup['imgt'].apply(lambda x: x.replace('TR', 'TCR').replace('-', '-0'))
    lookup['adaptive'] = lookup['adaptive'].apply(lambda x: pad_single_digit(x))
    lookup['adaptive_v2'] = lookup['adaptive']

    # If converting from 10X will just need the first *01 allele
    from_tenx = lookup.groupby('tenx').first()

    # Save
    lookup.to_csv(data_dir + '/lookup.csv')
    from_tenx.to_csv(data_dir + '/lookup_from_tenx.csv')


def parse_imgt_fasta(infile):
    '''Extract gene names from an IMGT reference fasta.

    :param infile: IMGT reference fasta for TRAV, TRAJ, TRBV, or TRBJ
    :type infile: str
    :return: Gene names
    :rtype: list
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
    '''Extract gene names from all IMGT reference fasta files in a folder.

    :param data_dir: Directory containing IMGT reference fasta files
    :type data_dir: str
    :return: Gene names
    :rtype: dataframe
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


def pad_single_digit(s):
    '''Add a '0' to single-digit gene names to match double-digit adaptive format.
    '''
    # Match strings that are six characters before a '-' etc. (e.g. TCRBV1-)
    pattern = r'^\w{6}(?=\W)'
    match = re.match(pattern, s)
    if match:
        # Add a '0' after the letters (e.g.g TCRBV becomes TCRBV0)
        return s.replace(match.group(0)[4], match.group(0)[4] + '0')
    else:
        return s


if __name__ == "__main__":
    main()