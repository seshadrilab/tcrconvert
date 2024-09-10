import os
import pandas as pd
from tcrconvert import build_lookup

def test_parse_imgt_fasta():
    f = os.path.dirname(__file__) + '/data/test_traj.fa'
    assert build_lookup.parse_imgt_fasta(f) == ['TRAJ1*01',
                                                'TRAJ10*01',
                                                'TRAJ11*01',
                                                'TRAJ12*01',
                                                'TRAJ13*01',
                                                'TRAJ13*02',
                                                'TRAJ14*01',
                                                'TRAJ15*01']


def test_extract_imgt_genes():
    testdir = os.path.dirname(__file__) + '/data/'
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


def test_build_lookup_from_fastas():
    testdir = os.path.dirname(__file__) + '/data/'
    build_lookup.build_lookup_from_fastas(testdir)

    with open(testdir + 'lookup_from_tenx.csv') as lookup10x: 
        assert lookup10x.read() == 'tenx,imgt,adaptive,adaptivev2\n' \
                                    'TRAJ1,TRAJ1*01,TCRAJ01*01,TCRAJ01*01\n' \
                                    'TRAJ10,TRAJ10*01,TCRAJ10*01,TCRAJ10*01\n' \
                                    'TRAJ11,TRAJ11*01,TCRAJ11*01,TCRAJ11*01\n' \
                                    'TRAJ12,TRAJ12*01,TCRAJ12*01,TCRAJ12*01\n' \
                                    'TRAJ13,TRAJ13*01,TCRAJ13*01,TCRAJ13*01\n' \
                                    'TRAJ14,TRAJ14*01,TCRAJ14*01,TCRAJ14*01\n' \
                                    'TRAJ15,TRAJ15*01,TCRAJ15*01,TCRAJ15*01\n' \
                                    'TRAV1-1,TRAV1-1*01,TCRAV01-01*01,TCRAV01-01*01\n' \
                                    'TRAV1-2,TRAV1-2*01,TCRAV01-02*01,TCRAV01-02*01\n'

    with open(testdir + 'lookup.csv') as lookup: 
        assert lookup.read() == 'imgt,tenx,adaptive,adaptivev2\n' \
                                'TRAJ1*01,TRAJ1,TCRAJ01*01,TCRAJ01*01\n' \
                                'TRAJ10*01,TRAJ10,TCRAJ10*01,TCRAJ10*01\n' \
                                'TRAJ11*01,TRAJ11,TCRAJ11*01,TCRAJ11*01\n' \
                                'TRAJ12*01,TRAJ12,TCRAJ12*01,TCRAJ12*01\n' \
                                'TRAJ13*01,TRAJ13,TCRAJ13*01,TCRAJ13*01\n' \
                                'TRAJ13*02,TRAJ13,TCRAJ13*02,TCRAJ13*02\n' \
                                'TRAJ14*01,TRAJ14,TCRAJ14*01,TCRAJ14*01\n' \
                                'TRAJ15*01,TRAJ15,TCRAJ15*01,TCRAJ15*01\n' \
                                'TRAV1-1*01,TRAV1-1,TCRAV01-01*01,TCRAV01-01*01\n' \
                                'TRAV1-1*02,TRAV1-1,TCRAV01-01*02,TCRAV01-01*02\n' \
                                'TRAV1-2*01,TRAV1-2,TCRAV01-02*01,TCRAV01-02*01\n'
