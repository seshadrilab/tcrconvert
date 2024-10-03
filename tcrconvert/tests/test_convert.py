import pytest
import pandas as pd
import tcrconvert

imgt_df = pd.DataFrame({'v_gene': ['TRAV12-1*01', 'TRBV15*01'],
                        'd_gene': [pd.NA, 'TRBD1*01'],
                        'j_gene': ['TRAJ16*01', 'TRBJ2-5*01'],
                        'c_gene': ['TRAC*01', 'TRBC2*01'],
                        'cdr3': ['CAVLIF', 'CASSGF']})

tenx_df = pd.DataFrame({'v_gene': ['TRAV12-1', 'TRBV15'],
                        'd_gene': [pd.NA, 'TRBD1'],
                        'j_gene': ['TRAJ16', 'TRBJ2-5'],
                        'c_gene': ['TRAC', 'TRBC2'],
                        'cdr3': ['CAVLIF', 'CASSGF']})

adapt_df = pd.DataFrame({'v_resolved': ['TCRAV12-01*01', 'TCRBV15-01*01'],
                         'd_resolved': [pd.NA, 'TCRBD01-01*01'],
                         'j_resolved': ['TCRAJ16-01*01', 'TCRBJ02-05*01'],
                         'cdr3_amino_acid': ['CAVLIF', 'CASSGF']})

adapt_v2_df = adapt_df.rename(columns={'v_resolved': 'vMaxResolved', 
                                       'd_resolved': 'dMaxResolved',
                                       'j_resolved': 'jMaxResolved',
                                       'cdr3_amino_acid': 'aminoAcid'})

custom_df = imgt_df.rename(columns={'v_gene': 'myV',
                                    'd_gene': 'myD',
                                    'j_gene': 'myJ',
                                    'c_gene': 'myC',
                                    'cdr3': 'myCDR3'})

tenx_to_adapt_df = adapt_df.rename(columns={'v_resolved': 'v_gene', 
                                            'd_resolved': 'd_gene',
                                            'j_resolved': 'j_gene',
                                            'cdr3_amino_acid': 'cdr3'})
tenx_to_adapt_df.insert(3, 'c_gene', [pd.NA, pd.NA])

adapt_to_tenx_df = tenx_df.rename(columns={'v_gene': 'v_resolved',
                                           'd_gene': 'd_resolved',
                                           'j_gene': 'j_resolved',
                                           'cdr3': 'cdr3_amino_acid'}).\
                            drop(columns='c_gene')

adapt_to_imgt_df = imgt_df.rename(columns={'v_gene': 'v_resolved',
                                           'd_gene': 'd_resolved',
                                           'j_gene': 'j_resolved',
                                           'cdr3': 'cdr3_amino_acid'}).\
                            drop(columns='c_gene')

adaptv2_to_tenx_df = tenx_df.rename(columns={'v_gene': 'vMaxResolved',
                                             'd_gene': 'dMaxResolved',
                                             'j_gene': 'jMaxResolved',
                                             'cdr3': 'aminoAcid'}).\
                            drop(columns='c_gene')

adaptv2_to_imgt_df = imgt_df.rename(columns={'v_gene': 'vMaxResolved',
                                             'd_gene': 'dMaxResolved',
                                             'j_gene': 'jMaxResolved',
                                             'cdr3': 'aminoAcid'}).\
                            drop(columns='c_gene')

custom_to_tenx_df = tenx_df.rename(columns={'v_gene': 'myV',
                                            'd_gene': 'myD',
                                            'j_gene': 'myJ',
                                            'c_gene': 'myC',
                                            'cdr3': 'myCDR3'})

adapt_no_allele_df = pd.DataFrame({'v_resolved': ['TCRAV12-01', 'TCRBV15-01*01'],
                                    'd_resolved': [pd.NA, 'TCRBD01-01'],
                                    'j_resolved': ['TCRAJ16-01*01', 'TCRBJ02-05'],
                                    'cdr3_amino_acid': ['CAVLIF', 'CASSGF']})

@pytest.mark.parametrize('df, frm, to, species, frm_cols, out', [
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
    (custom_df, 'imgt', 'tenx', 'human', ['myV', 'myD', 'myJ', 'myC'], custom_to_tenx_df),
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
    (custom_df, 'imgt', 'tenx', 'mouse', ['myV', 'myD', 'myJ', 'myC'], custom_to_tenx_df),
    # RHESUS MACAQUE
    (tenx_df, 'tenx', 'adaptive', 'macaque', None, tenx_to_adapt_df),
    (tenx_df, 'tenx', 'adaptivev2', 'macaque', None, tenx_to_adapt_df),
    (adapt_df, 'adaptive', 'tenx', 'macaque', None, adapt_to_tenx_df),
    (adapt_v2_df, 'adaptivev2', 'tenx', 'macaque', None, adaptv2_to_tenx_df),
    (tenx_df, 'tenx', 'imgt', 'macaque', None, imgt_df),
    (imgt_df, 'imgt', 'tenx', 'macaque', None, tenx_df),
    (imgt_df, 'imgt', 'adaptive', 'macaque', None, tenx_to_adapt_df),
    (imgt_df, 'imgt', 'adaptivev2', 'macaque', None, tenx_to_adapt_df),
    (adapt_df, 'adaptive', 'imgt', 'macaque', None, adapt_to_imgt_df),
    (adapt_v2_df, 'adaptivev2', 'imgt', 'macaque', None, adaptv2_to_imgt_df),
    (custom_df, 'imgt', 'tenx', 'macaque', ['myV', 'myD', 'myJ', 'myC'], custom_to_tenx_df),
    # Some Adaptive genes without allele
    (adapt_no_allele_df, 'adaptive', 'imgt', 'human', None, adapt_to_imgt_df)])
def test_convert_gene(df, frm, to, species, frm_cols, out):
    result = tcrconvert.convert_gene(df, frm, to, species, frm_cols)
    # Standardize the NA values so we can check for equality
    test_result = result.fillna('blank')
    test_out = out.fillna('blank')
    pd.testing.assert_frame_equal(test_result, test_out)