import pytest
import pandas as pd
import os
from importlib.resources import files
import logging
from tcrconvert import convert

imgt_df = pd.DataFrame(
    {
        'v_gene': ['TRAV12-1*01', 'TRBV15*01'],
        'd_gene': [pd.NA, 'TRBD1*01'],
        'j_gene': ['TRAJ16*01', 'TRBJ2-5*01'],
        'c_gene': ['TRAC*01', 'TRBC2*01'],
        'cdr3': ['CAVLIF', 'CASSGF'],
    }
)

tenx_df = pd.DataFrame(
    {
        'v_gene': ['TRAV12-1', 'TRBV15'],
        'd_gene': [pd.NA, 'TRBD1'],
        'j_gene': ['TRAJ16', 'TRBJ2-5'],
        'c_gene': ['TRAC', 'TRBC2'],
        'cdr3': ['CAVLIF', 'CASSGF'],
    }
)

adapt_df = pd.DataFrame(
    {
        'v_resolved': ['TCRAV12-01*01', 'TCRBV15-01*01'],
        'd_resolved': [pd.NA, 'TCRBD01-01*01'],
        'j_resolved': ['TCRAJ16-01*01', 'TCRBJ02-05*01'],
        'cdr3_amino_acid': ['CAVLIF', 'CASSGF'],
    }
)

adapt_v2_df = adapt_df.rename(
    columns={
        'v_resolved': 'vMaxResolved',
        'd_resolved': 'dMaxResolved',
        'j_resolved': 'jMaxResolved',
        'cdr3_amino_acid': 'aminoAcid',
    }
)

custom_df = imgt_df.rename(
    columns={
        'v_gene': 'myV',
        'd_gene': 'myD',
        'j_gene': 'myJ',
        'c_gene': 'myC',
        'cdr3': 'myCDR3',
    }
)

custom_vj_tenx_df = pd.DataFrame(
    {
        'myV': ['TRAV12-1', 'TRBV15'],
        'myD': [pd.NA, 'TRBD1*01'],
        'myJ': ['TRAJ16', 'TRBJ2-5'],
        'myC': ['TRAC*01', 'TRBC2*01'],
        'myCDR3': ['CAVLIF', 'CASSGF'],
    }
)

tenx_to_adapt_df = adapt_df.rename(
    columns={
        'v_resolved': 'v_gene',
        'd_resolved': 'd_gene',
        'j_resolved': 'j_gene',
        'cdr3_amino_acid': 'cdr3',
    }
)
tenx_to_adapt_df.insert(3, 'c_gene', [pd.NA, pd.NA])

adapt_to_tenx_df = tenx_df.rename(
    columns={
        'v_gene': 'v_resolved',
        'd_gene': 'd_resolved',
        'j_gene': 'j_resolved',
        'cdr3': 'cdr3_amino_acid',
    }
).drop(columns='c_gene')

adapt_to_imgt_df = imgt_df.rename(
    columns={
        'v_gene': 'v_resolved',
        'd_gene': 'd_resolved',
        'j_gene': 'j_resolved',
        'cdr3': 'cdr3_amino_acid',
    }
).drop(columns='c_gene')

adaptv2_to_tenx_df = tenx_df.rename(
    columns={
        'v_gene': 'vMaxResolved',
        'd_gene': 'dMaxResolved',
        'j_gene': 'jMaxResolved',
        'cdr3': 'aminoAcid',
    }
).drop(columns='c_gene')

adaptv2_to_imgt_df = imgt_df.rename(
    columns={
        'v_gene': 'vMaxResolved',
        'd_gene': 'dMaxResolved',
        'j_gene': 'jMaxResolved',
        'cdr3': 'aminoAcid',
    }
).drop(columns='c_gene')

custom_to_tenx_df = tenx_df.rename(
    columns={
        'v_gene': 'myV',
        'd_gene': 'myD',
        'j_gene': 'myJ',
        'c_gene': 'myC',
        'cdr3': 'myCDR3',
    }
)

adapt_no_allele_df = pd.DataFrame(
    {
        'v_resolved': ['TCRAV12-01', 'TCRBV15-01*01'],
        'd_resolved': [pd.NA, 'TCRBD01-01'],
        'j_resolved': ['TCRAJ16-01*01', 'TCRBJ02-05'],
        'cdr3_amino_acid': ['CAVLIF', 'CASSGF'],
    }
)


@pytest.mark.parametrize(
    'df, frm, to, species, frm_cols, out',
    [
        # 10X <-> Adaptive
        (tenx_df, 'tenx', 'adaptive', 'human', None, tenx_to_adapt_df),
        (tenx_df, 'tenx', 'adaptivev2', 'human', None, tenx_to_adapt_df),
        (adapt_df, 'adaptive', 'tenx', 'human', None, adapt_to_tenx_df),
        (adapt_v2_df, 'adaptivev2', 'tenx', 'human', None, adaptv2_to_tenx_df),
        # 10X <-> IMGT
        (tenx_df, 'tenx', 'imgt', 'human', None, imgt_df),
        (imgt_df, 'imgt', 'tenx', 'human', None, tenx_df),
        # IMGT <-> Adaptive
        (imgt_df, 'imgt', 'adaptive', 'human', None, tenx_to_adapt_df),
        (imgt_df, 'imgt', 'adaptivev2', 'human', None, tenx_to_adapt_df),
        (adapt_df, 'adaptive', 'imgt', 'human', None, adapt_to_imgt_df),
        (adapt_v2_df, 'adaptivev2', 'imgt', 'human', None, adaptv2_to_imgt_df),
        # Custom column names
        (
            custom_df,
            'imgt',
            'tenx',
            'human',
            ['myV', 'myD', 'myJ', 'myC'],
            custom_to_tenx_df,
        ),
        # MOUSE
        (tenx_df, 'tenx', 'adaptive', 'mouse', None, tenx_to_adapt_df),
        (tenx_df, 'tenx', 'adaptivev2', 'mouse', None, tenx_to_adapt_df),
        (adapt_df, 'adaptive', 'tenx', 'mouse', None, adapt_to_tenx_df),
        (adapt_v2_df, 'adaptivev2', 'tenx', 'mouse', None, adaptv2_to_tenx_df),
        (tenx_df, 'tenx', 'imgt', 'mouse', None, imgt_df),
        (imgt_df, 'imgt', 'tenx', 'mouse', None, tenx_df),
        (imgt_df, 'imgt', 'adaptive', 'mouse', None, tenx_to_adapt_df),
        (imgt_df, 'imgt', 'adaptivev2', 'mouse', None, tenx_to_adapt_df),
        (adapt_df, 'adaptive', 'imgt', 'mouse', None, adapt_to_imgt_df),
        (adapt_v2_df, 'adaptivev2', 'imgt', 'mouse', None, adaptv2_to_imgt_df),
        (
            custom_df,
            'imgt',
            'tenx',
            'mouse',
            ['myV', 'myD', 'myJ', 'myC'],
            custom_to_tenx_df,
        ),
        # RHESUS MACAQUE
        (tenx_df, 'tenx', 'adaptive', 'rhesus', None, tenx_to_adapt_df),
        (tenx_df, 'tenx', 'adaptivev2', 'rhesus', None, tenx_to_adapt_df),
        (adapt_df, 'adaptive', 'tenx', 'rhesus', None, adapt_to_tenx_df),
        (adapt_v2_df, 'adaptivev2', 'tenx', 'rhesus', None, adaptv2_to_tenx_df),
        (tenx_df, 'tenx', 'imgt', 'rhesus', None, imgt_df),
        (imgt_df, 'imgt', 'tenx', 'rhesus', None, tenx_df),
        (imgt_df, 'imgt', 'adaptive', 'rhesus', None, tenx_to_adapt_df),
        (imgt_df, 'imgt', 'adaptivev2', 'rhesus', None, tenx_to_adapt_df),
        (adapt_df, 'adaptive', 'imgt', 'rhesus', None, adapt_to_imgt_df),
        (adapt_v2_df, 'adaptivev2', 'imgt', 'rhesus', None, adaptv2_to_imgt_df),
        (
            custom_df,
            'imgt',
            'tenx',
            'rhesus',
            ['myV', 'myD', 'myJ', 'myC'],
            custom_to_tenx_df,
        ),
        # Some Adaptive genes without allele
        (adapt_no_allele_df, 'adaptive', 'imgt', 'human', None, adapt_to_imgt_df),
        # Confirm won't convert non-VDJC gene column to NAs
        (
            custom_df,
            'imgt',
            'tenx',
            'human',
            ['myV', 'myJ', 'myCDR3'],
            custom_vj_tenx_df,
        ),
    ],
)
def test_convert_gene(df, frm, to, species, frm_cols, out):
    result = convert.convert_gene(df, frm, to, species, frm_cols)
    # Standardize the NA values so we can check for equality
    test_result = result.fillna('blank')
    test_out = out.fillna('blank')
    pd.testing.assert_frame_equal(test_result, test_out)


def test_convert_gene_input():
    # Same input and output format
    with pytest.raises(ValueError):
        convert.convert_gene(tenx_df, 'tenx', 'tenx')

    # Empty input dataframe
    with pytest.raises(TypeError):
        convert.convert_gene(pd.DataFrame, 'tenx', 'imgt')


def test_choose_lookup():
    lookup_dir = os.path.join(files('tcrconvert'), 'data', 'human')
    lookup_tenx = os.path.join(lookup_dir, 'lookup_from_tenx.csv')
    lookup_adapt = os.path.join(lookup_dir, 'lookup_from_adaptive.csv')
    lookup_imgt = os.path.join(lookup_dir, 'lookup.csv')

    # Species we don't have lookups for
    with pytest.raises(FileNotFoundError):
        convert.choose_lookup('tenx', 'imgt', 'non-existent-species')

    # Different 'frm' formats
    assert convert.choose_lookup('tenx', 'imgt') == lookup_tenx
    assert convert.choose_lookup('adaptivev2', 'imgt') == lookup_adapt
    assert convert.choose_lookup('imgt', 'tenx') == lookup_imgt


def test_which_frm_cols():
    col_ref = {
        'adaptive': ['v_resolved', 'd_resolved', 'j_resolved'],
        'adaptivev2': ['vMaxResolved', 'dMaxResolved', 'jMaxResolved'],
        'imgt': ['v_gene', 'd_gene', 'j_gene', 'c_gene'],
        'tenx': ['v_gene', 'd_gene', 'j_gene', 'c_gene'],
    }

    assert convert.which_frm_cols(adapt_df, 'adaptive') == col_ref['adaptive']
    assert convert.which_frm_cols(adapt_v2_df, 'adaptivev2') == col_ref['adaptivev2']
    assert convert.which_frm_cols(imgt_df, 'imgt') == col_ref['imgt']
    assert convert.which_frm_cols(tenx_df, 'tenx') == col_ref['tenx']

    # Custom columns
    custom_col = ['myV', 'myD', 'myJ', 'myC']
    assert convert.which_frm_cols(custom_df, 'tenx', frm_cols=custom_col) == custom_col

    # Non-existent column
    with pytest.raises(ValueError):
        convert.which_frm_cols(tenx_df, 'tenx', frm_cols=['v_gene', 'j_gene', 'x_gene'])


def test_choose_lookup_verbose(caplog):
    # Capture messages when verbose=True
    with caplog.at_level(logging.INFO):
        convert.choose_lookup('tenx', 'adaptive')
        assert 'Converting from 10X. Using *01 as allele for all genes.' in caplog.text
        caplog.clear()

        convert.choose_lookup('adaptive', 'imgt')
        assert (
            'Converting from Adaptive to IMGT. Using *01 for genes lacking alleles.'
            in caplog.text
        )

    # Ensure no messages when verbose=False
    caplog.clear()
    with caplog.at_level(logging.INFO):
        convert.choose_lookup('tenx', 'adaptive', verbose=False)
        assert not caplog.text

        convert.choose_lookup('adaptive', 'imgt', verbose=False)
        assert not caplog.text


def test_which_frm_cols_verbose(caplog):
    custom_df = pd.DataFrame(
        {
            'myV': ['TRAV12-1*01', 'TRBV15*01'],
            'myD': [None, 'TRBD1*01'],
            'myJ': ['TRAJ16*0', 'TRBJ2-5*01'],
            'myC': ['TRAC*01', 'TRBC2*01'],
            'myCDR3': ['CAVLIF', 'CASSGF'],
        }
    )

    # Custom columns
    with caplog.at_level(logging.INFO):
        convert.which_frm_cols(custom_df, 'imgt', frm_cols=['myV', 'myJ'])
        assert "Using custom column names: ['myV', 'myJ']" in caplog.text

    # Custom columns with verbose=False
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        convert.which_frm_cols(
            custom_df, 'imgt', frm_cols=['myV', 'myJ'], verbose=False
        )
        assert not caplog.text

    # IMGT format without column names, still expect warnings with verbose=False
    caplog.clear()
    with caplog.at_level(logging.WARNING):
        convert.which_frm_cols(custom_df, 'imgt', verbose=False)
        assert (
            "No column names for IMGT data. Using 10X columns: ['v_gene', 'd_gene', 'j_gene', 'c_gene']"
            in caplog.text
        )


def test_convert_gene_verbose(caplog):
    tenx_df_bad = pd.DataFrame(
        {
            'v_gene': ['TRAV12-1', 'TRBV15'],
            'd_gene': [None, 'TRBD1'],
            'j_gene': ['TRAJ16', 'BAD_J_GENE'],
            'c_gene': ['TRAC', 'TRBC2'],
            'cdr3': ['CAVLIF', 'CASSGF'],
        }
    )

    with caplog.at_level(logging.WARNING):
        convert.convert_gene(tenx_df_bad, 'tenx', 'adaptive', verbose=False)

        # Verify expected warnings
        assert 'Adaptive only captures VDJ genes; C genes will be NA.' in caplog.text
        assert (
            'These genes are not in IMGT for this species and will be replaced with NA'
            in caplog.text
        )
