import os
import pandas as pd
from tcrconverter import build_lookup

def test_parse_imgt_fasta():
    f = os.path.dirname(__file__) + '/test_traj.fa'
    assert build_lookup.parse_imgt_fasta(f) == ['TRAJ1*01',
                                                'TRAJ10*01',
                                                'TRAJ11*01',
                                                'TRAJ12*01',
                                                'TRAJ13*01',
                                                'TRAJ13*02',
                                                'TRAJ14*01',
                                                'TRAJ15*01']


def test_extract_imgt_genes():
    testdir = os.path.dirname(__file__) + '/'
    outdf = pd.DataFrame({'imgt': ['TRAJ1*01', 'TRAJ10*01', 'TRAJ11*01', 
                                   'TRAJ12*01', 'TRAJ13*01', 'TRAJ13*02', 
                                   'TRAJ14*01', 'TRAJ15*01', 'TRAV1-1*01', 
                                   'TRAV1-1*02', 'TRAV1-2*01']})
    testdf = build_lookup.extract_imgt_genes(testdir)
    pd.testing.assert_frame_equal(outdf, testdf)


def test_pad_single_digit():
    s1 = 'TCRBV1-2'
    s2 = 'TCRBV11-2'
    assert build_lookup.pad_single_digit(s1) == 'TCRBV01-2'
    assert build_lookup.pad_single_digit(s2) == 'TCRBV11-2'