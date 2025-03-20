import pandas as pd
import os
import tempfile
import shutil
import pytest
from unittest.mock import patch
from tcrconvert import build_lookup, utils


def test_parse_imgt_fasta():
    fasta = utils.get_example_path('fasta_dir/test_trav.fa')
    assert build_lookup.parse_imgt_fasta(fasta) == [
        'TRAV1-1*01',
        'TRAV1-1*02',
        'TRAV1-2*01',
        'TRAV14/DV4*01',
        'TRAV38-1*01',
        'TRAV38-2/DV8*01',
        'TRAC*01',
    ]


def test_extract_imgt_genes():
    fastadir = utils.get_example_path('fasta_dir')
    df = build_lookup.extract_imgt_genes(fastadir)
    outdf = pd.DataFrame(
        {
            'imgt': [
                'TRAC*01',
                'TRAV1-1*01',
                'TRAV1-1*02',
                'TRAV1-2*01',
                'TRAV14/DV4*01',
                'TRAV38-1*01',
                'TRAV38-2/DV8*01',
                'TRBV29-1*01',
                'TRBV29-1*02',
                'TRBV29/OR9-2*01',
            ]
        }
    )
    pd.testing.assert_frame_equal(df, outdf)


def test_add_dash_one():
    gene_str1 = 'TRBV2*01'
    gene_str2 = 'TRBV1-01*01'
    assert build_lookup.add_dash_one(gene_str1) == 'TRBV2-01*01'
    assert build_lookup.add_dash_one(gene_str2) == gene_str2


def test_pad_single_digit():
    gene_str1 = 'TCRBV1-2'
    gene_str2 = 'TCRBV11-2'
    assert build_lookup.pad_single_digit(gene_str1) == 'TCRBV01-2'
    assert build_lookup.pad_single_digit(gene_str2) == gene_str2


def test_build_lookup_from_fastas():
    fastadir = utils.get_example_path('fasta_dir')

    # Create mock folder in temporary directory to write to
    mock_path = os.path.join(tempfile.gettempdir(), 'tcrconvert_tmp')
    os.makedirs(mock_path, exist_ok=True)

    with patch('platformdirs.user_data_dir', return_value=str(mock_path)):
        build_lookup.build_lookup_from_fastas(fastadir, species='rabbit')

    with open(mock_path + '/rabbit/lookup_from_adaptive.csv') as lookupadapt:
        assert (
            lookupadapt.read() == 'adaptive,adaptivev2,imgt,tenx\n'
            'TCRAV01-01,TCRAV01-01,TRAV1-1*01,TRAV1-1\n'
            'TCRAV01-01*01,TCRAV01-01*01,TRAV1-1*01,TRAV1-1\n'
            'TCRAV01-01*02,TCRAV01-01*02,TRAV1-1*02,TRAV1-1\n'
            'TCRAV01-02,TCRAV01-02,TRAV1-2*01,TRAV1-2\n'
            'TCRAV01-02*01,TCRAV01-02*01,TRAV1-2*01,TRAV1-2\n'
            'TCRAV14,TCRAV14,TRAV14/DV4*01,TRAV14/DV4\n'
            'TCRAV14*01,TCRAV14*01,TRAV14/DV4*01,TRAV14/DV4\n'
            'TCRAV14-01,TCRAV14-01,TRAV14/DV4*01,TRAV14/DV4\n'
            'TCRAV14-01*01,TCRAV14-01*01,TRAV14/DV4*01,TRAV14/DV4\n'
            'TCRAV38-01,TCRAV38-01,TRAV38-1*01,TRAV38-1\n'
            'TCRAV38-01*01,TCRAV38-01*01,TRAV38-1*01,TRAV38-1\n'
            'TCRAV38-02,TCRAV38-02,TRAV38-2/DV8*01,TRAV38-2/DV8\n'
            'TCRAV38-02*01,TCRAV38-02*01,TRAV38-2/DV8*01,TRAV38-2/DV8\n'
            'TCRBV29,TCRBV29,TRBV29-1*01,TRBV29-1\n'
            'TCRBV29*01,TCRBV29*01,TRBV29-1*01,TRBV29-1\n'
            'TCRBV29*02,TCRBV29*02,TRBV29-1*02,TRBV29-1\n'
            'TCRBV29-01,TCRBV29-01,TRBV29-1*01,TRBV29-1\n'
            'TCRBV29-01*01,TCRBV29-01*01,TRBV29-1*01,TRBV29-1\n'
            'TCRBV29-01*02,TCRBV29-01*02,TRBV29-1*02,TRBV29-1\n'
            'TCRBV29-or09_02,TCRBV29-or09_02,TRBV29/OR9-2*01,TRBV29/OR9-2\n'
            'TCRBV29-or09_02*01,TCRBV29-or09_02*01,TRBV29/OR9-2*01,TRBV29/OR9-2\n'
        )

    with open(mock_path + '/rabbit/lookup_from_tenx.csv') as lookup10x:
        assert (
            lookup10x.read() == 'tenx,imgt,adaptive,adaptivev2\n'
            'TRAC,TRAC*01,NoData,NoData\n'
            'TRAV1-1,TRAV1-1*01,TCRAV01-01*01,TCRAV01-01*01\n'
            'TRAV1-2,TRAV1-2*01,TCRAV01-02*01,TCRAV01-02*01\n'
            'TRAV14/DV4,TRAV14/DV4*01,TCRAV14-01*01,TCRAV14-01*01\n'
            'TRAV38-1,TRAV38-1*01,TCRAV38-01*01,TCRAV38-01*01\n'
            'TRAV38-2/DV8,TRAV38-2/DV8*01,TCRAV38-02*01,TCRAV38-02*01\n'
            'TRBV29-1,TRBV29-1*01,TCRBV29-01*01,TCRBV29-01*01\n'
            'TRBV29/OR9-2,TRBV29/OR9-2*01,TCRBV29-or09_02*01,TCRBV29-or09_02*01\n'
        )

    with open(mock_path + '/rabbit/lookup.csv') as lookup:
        assert (
            lookup.read() == 'imgt,tenx,adaptive,adaptivev2\n'
            'TRAC*01,TRAC,NoData,NoData\n'
            'TRAV1-1*01,TRAV1-1,TCRAV01-01*01,TCRAV01-01*01\n'
            'TRAV1-1*02,TRAV1-1,TCRAV01-01*02,TCRAV01-01*02\n'
            'TRAV1-2*01,TRAV1-2,TCRAV01-02*01,TCRAV01-02*01\n'
            'TRAV14/DV4*01,TRAV14/DV4,TCRAV14-01*01,TCRAV14-01*01\n'
            'TRAV38-1*01,TRAV38-1,TCRAV38-01*01,TCRAV38-01*01\n'
            'TRAV38-2/DV8*01,TRAV38-2/DV8,TCRAV38-02*01,TCRAV38-02*01\n'
            'TRBV29-1*01,TRBV29-1,TCRBV29-01*01,TCRBV29-01*01\n'
            'TRBV29-1*02,TRBV29-1,TCRBV29-01*02,TCRBV29-01*02\n'
            'TRBV29/OR9-2*01,TRBV29/OR9-2,TCRBV29-or09_02*01,TCRBV29-or09_02*01\n'
        )


def test_reject_invalid_species():
    # Create mock folder in temporary directory to write to
    mock_path = os.path.join(tempfile.gettempdir(), 'tcrconvert_tmp')
    os.makedirs(mock_path, exist_ok=True)

    # Generate a temporary directory to put input files in
    fastadir = os.path.join(mock_path, 'fastas')
    os.makedirs(fastadir, exist_ok=True)
    # Copy test files into it
    example_dir = utils.get_example_path(fastadir)
    example_fastas = [f for f in os.listdir(example_dir) if f.endswith('.fa')]
    for fasta in example_fastas:
        shutil.copy(os.path.join(example_dir, fasta), fastadir)

    # Mock `user_data_dir` inside `build_lookup_from_fastas`
    with patch('platformdirs.user_data_dir', return_value=str(mock_path)):
        invalid_species = [
            'rabbit\\',
            'rabbit/',
            'rabbit:',
            'rabbit*',
            'rabbit?',
            'rabbit<',
            'rabbit>',
            'rabbit|',
            'rabbit~',
            'rabbit`',
            'rabbit\n',
            'rabbi\t',
        ]

        # Attempt to create lookup tables with invalid species folder names
        for species in invalid_species:
            with pytest.raises(ValueError):
                build_lookup.build_lookup_from_fastas(fastadir, species)

    # Delete temp directories
    shutil.rmtree(mock_path, ignore_errors=True)
    shutil.rmtree(fastadir, ignore_errors=True)


def test_save_lookup_tables():
    mock_path = os.path.join(tempfile.gettempdir(), 'tcrconvert_tmp')
    os.makedirs(mock_path, exist_ok=True)
    dat = pd.read_csv(utils.get_example_path('fasta_dir/lookup.csv'))

    build_lookup.save_lookup(dat, mock_path, 'newlookup.csv')

    # Check that the file exists
    assert os.path.exists(os.path.join(mock_path, 'newlookup.csv'))

    # Delete temp directory
    shutil.rmtree(mock_path, ignore_errors=True)
