import pytest
import pandas as pd
from importlib.resources import files
import tcrconvert

imgt_df = pd.DataFrame({'v_gene': ['TRAV1-2*01', 'TRBV6-1*01'],
                        'j_gene': ['TRAJ12*01', 'TRBJ2-1*01'],
                        'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

tenx_df = pd.DataFrame({'v_gene': ['TRAV1-2', 'TRBV6-1'],
                        'j_gene': ['TRAJ12', 'TRBJ2-1'],
                        'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_df = pd.DataFrame({'v_resolved': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                         'j_resolved': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                         'cdr3_amino_acid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_v2_df = pd.DataFrame({'vMaxResolved': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                            'jMaxResolved': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                            'aminoAcid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

custom_df = pd.DataFrame({'myV': ['TRAV1-2*01', 'TRBV6-1*01'],
                          'myJ': ['TRAJ12*01', 'TRBJ2-1*01'],
                          'myCDR3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

extract_df = pd.DataFrame({'v_gene': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                           'j_gene': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                           'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

example_adaptive = pd.DataFrame([['AAACCTGAGACCACGA-1', True, 'AAACCTGAGACCACGA-1_contig_1', True, 521, 'TRA', 'TCRAV01-02*01', 'TRBD1', 'TCRAJ12*01', 'TRAC', True, True, 'CAVMDSSYKLIF', 'TGTGCTGTGATGGATAGCAGCTATAAATTGATCTTC', 1569, 2, 'clonotype16', 'clonotype16_consensus_1'],
                                 ['AAACCTGAGACCACGA-1', True, 'AAACCTGAGACCACGA-1_contig_3', True, 584, 'TRB', 'TCRBV06-01*01', 'TRBD2', 'TCRBJ02-01*01', 'TRBC2', True, True, 'CASSGLAGGYNEQFF', 'TGTGCCAGCAGTGGACTAGCGGGGGGCTACAATGAGCAGTTCTTC', 5238, 7, 'clonotype16', 'clonotype16_consensus_2'],
                                 ['AAACCTGAGGCTCTTA-1', True, 'AAACCTGAGGCTCTTA-1_contig_1', True, 551, 'TRB', 'TCRBV06-04*01', 'TRBD2', 'TCRBJ02-03*01', 'TRBC2', True, True, 'CASSGVAGGTDTQYF', 'TGTGCCAGCAGTGGGGTAGCGGGAGGCACAGATACGCAGTATTTT', 3846, 4, 'clonotype26', 'clonotype26_consensus_2'],
                                 ['AAACCTGAGGCTCTTA-1', True, 'AAACCTGAGGCTCTTA-1_contig_2', True, 518, 'TRA', 'TCRAV01-02*01', 'TRBD1', 'TCRAJ33*01', 'TRAC', True, True, 'CAVKDSNYQLIW', 'TGTGCTGTGAAGGATAGCAACTATCAGTTAATCTGG', 2019, 2, 'clonotype26', 'clonotype26_consensus_1'],
                                 ['AAACCTGAGTGAACGC-1', True, 'AAACCTGAGTGAACGC-1_contig_1', True, 674, 'TRB', 'TCRBV02*01', 'TRBD1', 'TCRBJ01-02*01', 'TRBC1', True, True, 'CASNQGLNYGYTF', 'TGTGCCAGCAATCAGGGCCTTAACTATGGCTACACCTTC', 3002, 6, 'clonotype81', 'clonotype81_consensus_2']],
                                columns=['barcode', 'is_cell', 'contig_id',
                                         'high_confidence', 'length', 'chain',
                                         'v_resolved', 'd_gene', 'j_resolved', 'c_gene',
                                         'full_length', 'productive', 'cdr3_amino_acid',
                                         'cdr3_nt', 'reads', 'umis',
                                         'raw_clonotype_id', 'raw_consensus_id'])


@pytest.mark.parametrize('df, frm, to, rename_cols, frm_cols, out', [
    # 10X <-> Adaptive
    (tenx_df, 'tenx', 'adaptive', True, None, adapt_df),
    (tenx_df, 'tenx', 'adaptive_v2', True, None, adapt_v2_df),
    (adapt_df, 'adaptive', 'tenx', True, None, tenx_df),
    (adapt_v2_df, 'adaptive_v2', 'tenx', True, None, tenx_df),
    # Adaptive <-> IMGT
    (adapt_df, 'adaptive', 'imgt', True, None, imgt_df),
    (adapt_v2_df, 'adaptive_v2', 'imgt', True, None, imgt_df),
    (imgt_df, 'imgt', 'adaptive', True, None, adapt_df),
    (imgt_df, 'imgt', 'adaptive_v2', True, None, adapt_v2_df),
    # IMGT <-> 10X
    (tenx_df, 'tenx', 'imgt', True, None, imgt_df),
    (imgt_df, 'imgt', 'tenx', True, None, tenx_df),
    # Custom column names
    (custom_df, 'imgt', 'tenx', True, ['myV', 'myJ', 'myCDR3'], tenx_df),
    # Adaptive V2 format but user inputs "adaptive"
    (adapt_v2_df, 'adaptive', 'tenx', True, None, tenx_df)])
def test_convert_df(df, frm, to, rename_cols, frm_cols, out):
    result = tcrconvert.convert_df(df, frm, to, rename_cols, frm_cols)
    pd.testing.assert_frame_equal(result, out)


def test_load_convert():
    f = str(files('tcrconvert') / 'examples' / 'example_human_10x.csv')
    result = tcrconvert.load_convert(f, frm='tenx', to='adaptive')
    pd.testing.assert_frame_equal(result, example_adaptive)