import os
import pandas as pd
from tcrconvert import build_lookup

def test_parse_imgt_fasta():
    f = os.path.dirname(__file__) + '/data/test.fa'
    assert build_lookup.parse_imgt_fasta(f) == ['TRAV1*01',
                                                'TRAV14/DV4*01',
                                                'TRAV38-2/DV8*01',
                                                'TRBV29/OR9-2*01',
                                                'TRBVA/OR9-2*01']


def test_extract_imgt_genes():
    testdir = os.path.dirname(__file__) + '/data/'
    testdf = build_lookup.extract_imgt_genes(testdir)
    outdf = pd.DataFrame({'imgt': ['TRAV1*01',
                                    'TRAV14/DV4*01',
                                    'TRAV38-2/DV8*01',
                                    'TRBV29/OR9-2*01',
                                    'TRBVA/OR9-2*01']})
    pd.testing.assert_frame_equal(testdf, outdf)


def test_add_dash_one():
    s1 = 'TRBV2*01'
    s2 = 'TRBV1-01*01'
    assert build_lookup.add_dash_one(s1) == 'TRBV2-01*01'
    assert build_lookup.add_dash_one(s2) == s2


def test_pad_single_digit():
    s1 = 'TCRBV1-2'
    s2 = 'TCRBV11-2'
    assert build_lookup.pad_single_digit(s1) == 'TCRBV01-2'
    assert build_lookup.pad_single_digit(s2) == s2


def test_build_lookup_from_fastas():
    testdir = os.path.dirname(__file__) + '/data/'
    build_lookup.build_lookup_from_fastas(testdir)

    with open(testdir + 'lookup_from_tenx.csv') as lookup10x: 
        assert lookup10x.read() == 'tenx,imgt,adaptive,adaptivev2\n' \
                                   'TRAV1,TRAV1*01,TCRAV01-01*01,TCRAV01-01*01\n' \
                                   'TRAV14DV4,TRAV14/DV4*01,TCRAV14-01*01,TCRAV14-01*01\n' \
                                   'TRAV38-2DV8,TRAV38-2/DV8*01,TCRAV38-02*01,TCRAV38-02*01\n' \
                                   'TRBV29/OR9-2,TRBV29/OR9-2*01,TCRBV29-or09_02*01,TCRBV29-or09_02*01\n' \
                                   'TRBVA/OR9-2,TRBVA/OR9-2*01,TCRBVA-or09_02*01,TCRBVA-or09_02*01\n'

    with open(testdir + 'lookup.csv') as lookup: 
        assert lookup.read() == 'imgt,tenx,adaptive,adaptivev2\n' \
                                'TRAV1*01,TRAV1,TCRAV01-01*01,TCRAV01-01*01\n' \
                                'TRAV14/DV4*01,TRAV14DV4,TCRAV14-01*01,TCRAV14-01*01\n' \
                                'TRAV38-2/DV8*01,TRAV38-2DV8,TCRAV38-02*01,TCRAV38-02*01\n' \
                                'TRBV29/OR9-2*01,TRBV29/OR9-2,TCRBV29-or09_02*01,TCRBV29-or09_02*01\n' \
                                'TRBVA/OR9-2*01,TRBVA/OR9-2,TCRBVA-or09_02*01,TCRBVA-or09_02*01\n'
