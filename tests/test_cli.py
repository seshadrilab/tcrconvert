import os
import tempfile
from click.testing import CliRunner
from unittest.mock import patch
from tcrconvert import cli, utils

in_csv = utils.get_example_path('customcols.csv')
out_tsv = tempfile.gettempdir() + '/custom2adapt.tsv'


def test_build_lookup_from_fastas_cli():
    fastadir = utils.get_example_path('fasta_dir')

    # Create mock folder in temporary directory to write to
    mock_path = os.path.join(tempfile.gettempdir(), 'mock_data')
    os.makedirs(mock_path, exist_ok=True)

    with patch('platformdirs.user_data_dir', return_value=str(mock_path)):
        result = CliRunner().invoke(
            cli.entry_point,
            [
                'build',
                '-i',
                fastadir,
                '-s',
                'rabbit',
            ],
            catch_exceptions=False,
        )

    assert result.exit_code == 0

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


def test_convert_gene_cli(caplog):
    result = CliRunner().invoke(
        cli.entry_point,
        [
            'convert',
            '--input',
            in_csv,
            '--output',
            out_tsv,
            '--frm',
            'tenx',
            '--to',
            'adaptive',
            '--species',
            'mouse',
            '-c',
            'myVgene',
            '-c',
            'myDgene',
            '-c',
            'myJgene',
            '-c',
            'myCgene',
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0

    # Expected warning messages
    assert 'Adaptive only captures VDJ genes; C genes will be NA.' in caplog.text
    assert 'Converting from 10X. Using *01 as allele for all genes.' in caplog.text
    assert (
        'These genes are not in IMGT for this species and will be replaced with NA:'
        in caplog.text
    )
    assert " ['TRAV1-2', 'TRBV6-1', 'TRBV6-4']" in caplog.text

    # Expected output file contents
    with open(out_tsv) as output:
        assert (
            output.read()
            == 'myVgene	myDgene	myJgene	myCgene	myCDR3	antigen\n'
            '	TCRBD01-01*01	TCRAJ12-01*01		CAVMDSSYKLIF	Flu\n'
            '	TCRBD02-01*01	TCRBJ02-01*01		CASSGLAGGYNEQFF	Flu\n'
            '	TCRBD02-01*01	TCRBJ02-03*01		CASSGVAGGTDTQYF	CMV\n'
            '	TCRBD01-01*01	TCRAJ33-01*01		CAVKDSNYQLIW	CMV\n'
            'TCRBV02-01*01	TCRBD01-01*01	TCRBJ01-02*01		CASNQGLNYGYTF	CMV\n'
        )


def test_convert_gene_cli_errors():
    badinfile = os.path.dirname(__file__) + '/data/badinput.txt'
    badoutfile = os.path.dirname(__file__) + '/data/badoutput.txt'

    # Input not CSV/TSV
    result_in = CliRunner().invoke(
        cli.entry_point,
        [
            'convert',
            '--input',
            badinfile,
            '--output',
            out_tsv,
            '--frm',
            'tenx',
            '--to',
            'adaptive',
            '--species',
            'mouse',
            '-c',
            'myVgene',
            '-c',
            'myDgene',
            '-c',
            'myJgene',
            '-c',
            'myCgene',
        ],
        catch_exceptions=False,
    )

    assert result_in.exit_code != 0
    assert '"input" must be a .csv or .tsv file' in result_in.output

    # Output not CSV/TSV
    result_out = CliRunner().invoke(
        cli.entry_point,
        [
            'convert',
            '--input',
            in_csv,
            '--output',
            badoutfile,
            '--frm',
            'tenx',
            '--to',
            'adaptive',
            '--species',
            'mouse',
            '-c',
            'myVgene',
            '-c',
            'myDgene',
            '-c',
            'myJgene',
            '-c',
            'myCgene',
        ],
        catch_exceptions=False,
    )

    assert result_out.exit_code != 0
    assert '"output" must be a .csv or .tsv file' in result_out.output
