import os
from click.testing import CliRunner
from tcrconvert import cli, utils

in_csv = utils.get_example_path('customcols.csv')
out_tsv = utils.get_example_path('converted_custom_adapt.tsv')


def test_convert_gene_cli(caplog):
    result = CliRunner().invoke(cli.entry_point, [
        'convert',
        '--infile', in_csv,
        '--outfile', out_tsv,
        '--frm', 'tenx',
        '--to', 'adaptive',
        '--species', 'mouse',
        '-c', 'myVgene',
        '-c', 'myDgene',
        '-c', 'myJgene',
        '-c', 'myCgene'
    ], catch_exceptions=False)

    assert result.exit_code == 0

    # Expected warning messages
    assert "Adaptive only captures VDJ genes, any C genes will become NA." in caplog.text
    assert "Converting from 10X which lacks allele info. Choosing *01 as allele for all genes." in caplog.text
    assert "These genes are not in IMGT for this species and will be replaced with NA:" in caplog.text
    assert " ['TRAV1-2', 'TRBV6-1', 'TRBV6-4']" in caplog.text

    # Expected output file contents
    with open(out_tsv) as output:
        assert output.read() == 'myVgene	myDgene	myJgene	myCgene	myCDR3	antigen\n' \
                                '	TCRBD01-01*01	TCRAJ12-01*01		CAVMDSSYKLIF	Flu\n' \
                                '	TCRBD02-01*01	TCRBJ02-01*01		CASSGLAGGYNEQFF	Flu\n' \
                                '	TCRBD02-01*01	TCRBJ02-03*01		CASSGVAGGTDTQYF	CMV\n' \
                                '	TCRBD01-01*01	TCRAJ33-01*01		CAVKDSNYQLIW	CMV\n' \
                                'TCRBV02-01*01	TCRBD01-01*01	TCRBJ01-02*01		CASNQGLNYGYTF	CMV\n'


def test_convert_gene_cli_errors():
    badinfile = os.path.dirname(__file__) + '/data/badinput.txt'
    badoutfile = os.path.dirname(__file__) + '/data/badoutput.txt'

    # Input not CSV/TSV
    result_in = CliRunner().invoke(cli.entry_point, [
            'convert',
            '--infile', badinfile,
            '--outfile', out_tsv,
            '--frm', 'tenx',
            '--to', 'adaptive',
            '--species', 'mouse',
            '-c', 'myVgene',
            '-c', 'myDgene',
            '-c', 'myJgene',
            '-c', 'myCgene'
        ], catch_exceptions=False)
    
    assert result_in.exit_code != 0
    assert '"infile" must be a .csv or .tsv file' in result_in.output

    # Output not CSV/TSV
    result_out = CliRunner().invoke(cli.entry_point, [
            'convert',
            '--infile', in_csv,
            '--outfile', badoutfile,
            '--frm', 'tenx',
            '--to', 'adaptive',
            '--species', 'mouse',
            '-c', 'myVgene',
            '-c', 'myDgene',
            '-c', 'myJgene',
            '-c', 'myCgene'
        ], catch_exceptions=False)

    assert result_out.exit_code != 0
    assert '"outfile" must be a .csv or .tsv file' in result_out.output
