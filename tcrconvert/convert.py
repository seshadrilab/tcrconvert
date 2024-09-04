import pandas as pd
from importlib.resources import files

# Standard column names for different sources of TCR data
col_ref = {'adaptive': ('v_resolved', 'd_resolved', 'j_resolved', 'cdr3_amino_acid'),
           'adaptive_v2': ('vMaxResolved', 'dMaxResolved', 'jMaxResolved', 'aminoAcid'),
           'imgt': ('v_gene', 'd_gene', 'j_gene', 'cdr3'),
           'tenx': ('v_gene', 'd_gene', 'j_gene', 'cdr3')}







def convert_df(df,
                frm,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                to,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                rename_cols=True,
                frm_cols=None,
                species='human'):
    '''Convert TCR gene names from one naming convention to another.

    :param df: Dataframe of TCRs
    :type df: dataframe
    :param frm: Starting format of TCR data
    :type frm: str
    :param to: Format to convert TCR data to
    :type to: str
    :param rename_cols: Whether to convert column names to match new format
    :type rename_cols: bool, optional
    :param frm_cols: List or tuple of current V, J, and CDR3 column names
    :type frm_cols: list or tuple, optional
    :param species: Species data folder to use
    :type species: str, optional
    :return: Converted TCR data
    :rtype: dataframe
    '''

    # Determine if we're using the 10X lookup table
    if frm == 'tenx':
        lookup_f = files('tcrconvert') / 'data' / species / 'lookup_from_tenx.csv'
        print("CONVERTING FROM 10X: CHOOSING *01 AS ALLELE FOR ALL GENES")
    else:
        lookup_f = files('tcrconvert') / 'data' / species / 'lookup.csv'

    # Load lookup table
    try:
        lookup = pd.read_csv(lookup_f)
    except FileNotFoundError:
        print('Lookup table not found, please download IMGT reference FASTAs and run build_lookup_from_fastas()')

    # Figure out columns to use
    if frm == 'adaptive' and 'vMaxResolved' in df.columns:
        cols_from = 'adaptive_v2'
    elif frm == 'imgt' and not frm_cols:
        print('No column names provided for IMGT data, will use 10X column names:\n' 
                      + str(col_ref['tenx']))
        cols_from = 'tenx'
    else:
        cols_from = frm


# 'v_gene', 'd_gene', 'j_gene', 'cdr3'

    # Column names to use for input and output
    # TODO: What if converting to imgt and want to update columns?
    if frm_cols:
        v_from = frm_cols[0]
        d_from = frm_cols[1]
        j_from = frm_cols[2]
        cdr3_from = frm_cols[3]
    else:
        v_from = col_ref[cols_from][0]
        d_from = col_ref[cols_from][1]
        j_from = col_ref[cols_from][2]
        cdr3_from = col_ref[cols_from][3]

    if rename_cols:
        v_to = col_ref[to][0]
        d_to = col_ref[to][1]
        j_to = col_ref[to][2]
        cdr3_to = col_ref[to][3]
    else:
        v_to = v_from
        d_to = d_from
        j_to = j_from
        cdr3_to = cdr3_from

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
        print('These genes are not in the IMGT reference and have replaced with NA. '
          'Please repair gene names manually and re-run:\n',
          v_genes_bad + d_genes_bad + j_genes_bad)

    # Consolidate converted data
    out = pd.DataFrame()
    out[v_to] = v_genes[to]
    out[d_to] = d_genes[to]
    out[j_to] = j_genes[to]
    out[cdr3_to] = df[cdr3_from]

    # Swap out data in original dataframe
    out_df = df.copy()
    out_df[v_from] = out[v_to].values
    out_df[d_from] = out[d_to].values
    out_df[j_from] = out[j_to].values
    out_df[cdr3_from] = out[cdr3_to].values
    out_df.rename(columns={v_from: v_to,
                           j_from: j_to,
                           d_from: d_to,
                           cdr3_from: cdr3_to}, inplace = True)
    return out_df


def convert_file(infile,
                 outfile,
                 frm,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                 to,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                 rename_cols=True,
                 frm_cols=None,
                 species='human'):
    '''Convert file with TCR data from one naming convention to another.

    :param infile: Path to input file (e.g. '~/Downloads/indata.csv')
    :type infile: str
    :param outfile: Path for output file (e.g. '~/Downloads/outdata.tsv')
    :type outfile: str
    :param df: Dataframe of TCRs
    :type df: dataframe
    :param frm: Starting format of TCR data
    :type frm: str
    :param to: Format to convert TCR data to
    :type to: str
    :param rename_cols: Whether to convert column names to match new format
    :type rename_cols: bool, optional
    :param frm_cols: List or tuple of current V, J, and CDR3 column names
    :type frm_cols: list or tuple, optional
    :param species: Species data folder to use
    :type species: str, optional
    :return: Converted TCR data
    :rtype: dataframe
    '''
    # Load input
    if infile.endswith('.csv'):
        dat = pd.read_csv(infile)
    elif infile.endswith('.tsv'):
        dat = pd.read_csv(infile, sep='\t')
    elif infile.endswith('.xls') or infile.endswith('.xlsx') or infile.endswith('.ods'):
        dat = pd.read_excel(infile)
    else:
        print('Please input a CSV, TSV, Excel, or ODF Spreadsheet file')

    # Convert
    out_dat = convert_df(dat, frm=frm, to=to, rename_cols=rename_cols,
                         frm_cols=frm_cols, species=species)

    # Save output
    if outfile.endswith('.csv'):
        out_dat.to_csv(outfile, index=False)
    elif outfile.endswith('.tsv'):
        out_dat.to_csv(outfile, index=False, sep='\t')
    elif outfile.endswith('.xls') or outfile.endswith('.xlsx') or outfile.endswith('.ods'):
        out_dat.to_excel(outfile, index=False)
    else:
        print('Please use a CSV, TSV, Excel, or ODF Spreadsheet file as output.')


def load_convert(infile,
                 frm,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                 to,  # ['tenx', 'adaptive', 'adaptive_v2', 'imgt']
                 rename_cols=True,
                 frm_cols=None,
                 species='human'):
    '''Load a file with TCR data and convert from one naming convention to another.

    :param infile: Path to input file (e.g. '~/Downloads/indata.csv')
    :type infile: str
    :param df: Dataframe of TCRs
    :type df: dataframe
    :param frm: Starting format of TCR data
    :type frm: str
    :param to: Format to convert TCR data to
    :type to: str
    :param rename_cols: Whether to convert column names to match new format
    :type rename_cols: bool, optional
    :param frm_cols: List or tuple of current V, J, and CDR3 column names
    :type frm_cols: list or tuple, optional
    :param species: Species data folder to use
    :type species: str, optional
    :return: Converted TCR data
    :rtype: dataframe
    '''

    # Load input
    if infile.endswith('.csv'):
        dat = pd.read_csv(infile)
    elif infile.endswith('.tsv'):
        dat = pd.read_csv(infile, sep='\t')
    elif infile.endswith('.xls') or infile.endswith('.xlsx') or infile.endswith('.ods'):
        dat = pd.read_excel(infile)
    else:
        print('Please input a CSV, TSV, Excel, or ODF Spreadsheet file')

    # Convert
    out_dat = convert_df(dat, frm=frm, to=to, rename_cols=rename_cols,
                         frm_cols=frm_cols, species=species)

    return(out_dat)