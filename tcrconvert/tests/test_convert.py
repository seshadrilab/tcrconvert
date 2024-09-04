import pytest
import pandas as pd
import tcrconvert

tenx_df = pd.DataFrame({'v_gene': ['TRAV1-2', 'TRBV6-1'],
                        'd_gene': [pd.NA, 'TRBD2'],
                        'j_gene': ['TRAJ12', 'TRBJ2-1'],
                        'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

imgt_df = pd.DataFrame({'v_gene': ['TRAV1-2*01', 'TRBV6-1*01'],
                        'd_gene': [pd.NA, 'TRBD2*02'],
                        'j_gene': ['TRAJ12*01', 'TRBJ2-1*01'],
                        'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_df = pd.DataFrame({'v_resolved': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                         'd_resolved': [pd.NA, 'TCRBD02*02'],
                         'j_resolved': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                         'cdr3_amino_acid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_v2_df = pd.DataFrame({'vMaxResolved': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                            'dMaxResolved': [pd.NA, 'TCRBD02*02'],
                            'jMaxResolved': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                            'aminoAcid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

tenx_to_imgt_df = pd.DataFrame({'v_gene': ['TRAV1-2*01', 'TRBV6-1*01'],
                                'd_gene': [pd.NA, 'TRBD2*01'],
                                'j_gene': ['TRAJ12*01', 'TRBJ2-1*01'],
                                'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

tenx_to_adapt_df = pd.DataFrame({'v_gene': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                                 'd_gene': [pd.NA, 'TCRBD02*01'],
                                 'j_gene': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                                 'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

imgt_to_adapt_df = pd.DataFrame({'v_gene': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                                 'd_gene': [pd.NA, 'TCRBD02*02'],
                                 'j_gene': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                                 'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_to_tenx_df = pd.DataFrame({'v_resolved': ['TRAV1-2', 'TRBV6-1'],
                                 'd_resolved': [pd.NA, 'TRBD2'],
                                 'j_resolved': ['TRAJ12', 'TRBJ2-1'],
                                 'cdr3_amino_acid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_to_imgt_df = pd.DataFrame({'v_resolved': ['TRAV1-2*01', 'TRBV6-1*01'],
                                 'd_resolved': [pd.NA, 'TRBD2*02'],
                                 'j_resolved': ['TRAJ12*01', 'TRBJ2-1*01'],
                                 'cdr3_amino_acid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adaptv2_to_tenx_df = pd.DataFrame({'vMaxResolved': ['TRAV1-2', 'TRBV6-1'],
                                   'dMaxResolved': [pd.NA, 'TRBD2'],
                                   'jMaxResolved': ['TRAJ12', 'TRBJ2-1'],
                                   'aminoAcid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adaptv2_to_imgt_df = pd.DataFrame({'vMaxResolved': ['TRAV1-2*01', 'TRBV6-1*01'],
                                   'dMaxResolved': [pd.NA, 'TRBD2*02'],
                                   'jMaxResolved': ['TRAJ12*01', 'TRBJ2-1*01'],
                                   'aminoAcid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})
 
custom_df = pd.DataFrame({'myV': ['TRAV1-2*01', 'TRBV6-1*01'],
                          'myD': [pd.NA, 'TRBD2*02'],
                          'myJ': ['TRAJ12*01', 'TRBJ2-1*01'],
                          'myCDR3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

custom_to_tenx_df = pd.DataFrame({'myV': ['TRAV1-2', 'TRBV6-1'],
                                  'myD': [pd.NA, 'TRBD2'],
                                  'myJ': ['TRAJ12', 'TRBJ2-1'],
                                  'myCDR3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})


@pytest.mark.parametrize('df, frm, to, frm_cols, out', [
    # 10X <-> Adaptive
    (tenx_df, 'tenx', 'adaptive', None, tenx_to_adapt_df),
    (tenx_df, 'tenx', 'adaptive_v2', None, tenx_to_adapt_df),
    (adapt_df, 'adaptive', 'tenx', None, adapt_to_tenx_df),
    (adapt_v2_df, 'adaptive_v2', 'tenx', None, adaptv2_to_tenx_df),
    # 10X <-> IMGT
    (tenx_df, 'tenx', 'imgt', None, tenx_to_imgt_df),
    (imgt_df, 'imgt', 'tenx', None, tenx_df),
    # IMGT <-> Adaptive
    (imgt_df, 'imgt', 'adaptive', None, imgt_to_adapt_df),
    (imgt_df, 'imgt', 'adaptive_v2', None, imgt_to_adapt_df),
    (adapt_df, 'adaptive', 'imgt', None, adapt_to_imgt_df),
    (adapt_v2_df, 'adaptive_v2', 'imgt', None, adaptv2_to_imgt_df),
    # Custom column names
    (custom_df, 'imgt', 'tenx', ['myV', 'myD', 'myJ', 'myCDR3'], custom_to_tenx_df)])
    # TODO: test rhesus and mouse
def test_convert_vdj(df, frm, to, frm_cols, out):
    result = tcrconvert.convert_vdj(df, frm, to, frm_cols)
    # Standardize the NA values so we can check for equality
    result.fillna('blank', inplace=True)
    out.fillna('blank', inplace=True)
    pd.testing.assert_frame_equal(result, out)


# TODO: Test bad gene names (mangled)