import pytest
import pandas as pd
import tcrconvert

# TODO: Properly test C genes, gamma J genes (these are species-specific)

tenx_df = pd.DataFrame({'v_gene': ['TRAV12-1', 'TRBV16', 'TRGV3', 'TRDV3'],
                        'd_gene': [pd.NA, 'TRBD2', pd.NA, 'TRDD2'],
                        'j_gene': ['TRAJ12', 'TRBJ2-1', 'TRGJ1', 'TRDJ2'],
                        'cdr3': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

imgt_df = pd.DataFrame({'v_gene': ['TRAV12-1*01', 'TRBV16*02', 'TRGV3*01', 'TRDV3*01'],
                        'd_gene': [pd.NA, 'TRBD2*01', pd.NA, 'TRDD2*01'],
                        'j_gene': ['TRAJ12*01', 'TRBJ2-1*01', 'TRGJ1*01', 'TRDJ2*01'],
                        'cdr3': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

adapt_df = pd.DataFrame({'v_resolved': ['TCRAV12-01*01', 'TCRBV16*02', 'TCRGV03*01', 'TCRDV03*01'],
                         'd_resolved': [pd.NA, 'TCRBD02*01', pd.NA,'TCRD0D02*01'],
                         'j_resolved': ['TCRAJ12*01', 'TCRBJ02-01*01', 'TCRGJ01*01', 'TCRDJ02*01'],
                         'cdr3_amino_acid': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

adapt_v2_df = pd.DataFrame({'vMaxResolved': ['TCRAV12-01*01', 'TCRBV16*02', 'TCRGV03*01', 'TCRDV03*01'],
                            'dMaxResolved': [pd.NA, 'TCRBD02*01', pd.NA, 'TCRD0D02*01'],
                            'jMaxResolved': ['TCRAJ12*01', 'TCRBJ02-01*01', 'TCRGJ01*01', 'TCRDJ02*01'],
                            'aminoAcid': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

tenx_to_imgt_df = pd.DataFrame({'v_gene': ['TRAV12-1*01', 'TRBV16*01', 'TRGV3*01','TRDV3*01'],
                                'd_gene': [pd.NA, 'TRBD2*01', pd.NA, 'TRDD2*01'],
                                'j_gene': ['TRAJ12*01', 'TRBJ2-1*01', 'TRGJ1*01', 'TRDJ2*01'],
                                'cdr3': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

tenx_to_adapt_df = pd.DataFrame({'v_gene': ['TCRAV12-01*01', 'TCRBV16*01', 'TCRGV03*01', 'TCRDV03*01'],
                                 'd_gene': [pd.NA, 'TCRBD02*01', pd.NA, 'TCRD0D02*01'],
                                 'j_gene': ['TCRAJ12*01', 'TCRBJ02-01*01', 'TCRGJ01*01', 'TCRDJ02*01'],
                                 'cdr3': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

imgt_to_adapt_df = pd.DataFrame({'v_gene': ['TCRAV12-01*01', 'TCRBV16*02', 'TCRGV03*01', 'TCRDV03*01'],
                                 'd_gene': [pd.NA, 'TCRBD02*01', pd.NA, 'TCRD0D02*01'],
                                 'j_gene': ['TCRAJ12*01', 'TCRBJ02-01*01', 'TCRGJ01*01', 'TCRDJ02*01'],
                                 'cdr3': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

adapt_to_tenx_df = pd.DataFrame({'v_resolved': ['TRAV12-1', 'TRBV16', 'TRGV3', 'TRDV3'],
                                 'd_resolved': [pd.NA, 'TRBD2', pd.NA, 'TRDD2'],
                                 'j_resolved': ['TRAJ12', 'TRBJ2-1', 'TRGJ1', 'TRDJ2'],
                                 'cdr3_amino_acid': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

adapt_to_imgt_df = pd.DataFrame({'v_resolved': ['TRAV12-1*01', 'TRBV16*02', 'TRGV3*01', 'TRDV3*01'],
                                 'd_resolved': [pd.NA, 'TRBD2*01', pd.NA, 'TRDD2*01'],
                                 'j_resolved': ['TRAJ12*01', 'TRBJ2-1*01', 'TRGJ1*01', 'TRDJ2*01'],
                                 'cdr3_amino_acid': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

adaptv2_to_tenx_df = pd.DataFrame({'vMaxResolved': ['TRAV12-1', 'TRBV16', 'TRGV3', 'TRDV3'],
                                   'dMaxResolved': [pd.NA, 'TRBD2', pd.NA, 'TRDD2'],
                                   'jMaxResolved': ['TRAJ12', 'TRBJ2-1', 'TRGJ1', 'TRDJ2'],
                                   'aminoAcid': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

adaptv2_to_imgt_df = pd.DataFrame({'vMaxResolved': ['TRAV12-1*01', 'TRBV16*02', 'TRGV3*01', 'TRDV3*01'],
                                   'dMaxResolved': [pd.NA, 'TRBD2*01', pd.NA, 'TRDD2*01'],
                                   'jMaxResolved': ['TRAJ12*01', 'TRBJ2-1*01', 'TRGJ1*01', 'TRDJ2*01'],
                                   'aminoAcid': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})
 
custom_df = pd.DataFrame({'myV': ['TRAV12-1*01', 'TRBV16*02', 'TRGV3*01', 'TRDV3*01'],
                          'myD': [pd.NA, 'TRBD2*01', pd.NA, 'TRDD2*01'],
                          'myJ': ['TRAJ12*01', 'TRBJ2-1*01', 'TRGJ1*01', 'TRDJ2*01'],
                          'myCDR3': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})

custom_to_tenx_df = pd.DataFrame({'myV': ['TRAV12-1', 'TRBV16', 'TRGV3', 'TRDV3'],
                                  'myD': [pd.NA, 'TRBD2', pd.NA, 'TRDD2'],
                                  'myJ': ['TRAJ12', 'TRBJ2-1', 'TRGJ1', 'TRDJ2'],
                                  'myCDR3': ['CAVLIF', 'CASSGF', 'CASAAF', 'CASLLF']})


@pytest.mark.parametrize('df, frm, to, frm_cols, species, out', [
    # 10X <-> Adaptive
    (tenx_df, 'tenx', 'adaptive', None, 'human', tenx_to_adapt_df),
    (tenx_df, 'tenx', 'adaptivev2', None, 'human', tenx_to_adapt_df),
    (adapt_df, 'adaptive', 'tenx', None, 'human', adapt_to_tenx_df),
    (adapt_v2_df, 'adaptivev2', 'tenx', None, 'human', adaptv2_to_tenx_df),
    # 10X <-> IMGT
    (tenx_df, 'tenx', 'imgt', None, 'human', tenx_to_imgt_df),
    (imgt_df, 'imgt', 'tenx', None, 'human', tenx_df),
    # IMGT <-> Adaptive
    (imgt_df, 'imgt', 'adaptive', None, 'human', imgt_to_adapt_df),
    (imgt_df, 'imgt', 'adaptivev2', None, 'human', imgt_to_adapt_df),
    (adapt_df, 'adaptive', 'imgt', None, 'human', adapt_to_imgt_df),
    (adapt_v2_df, 'adaptivev2', 'imgt', None, 'human', adaptv2_to_imgt_df),
    # Custom column names
    (custom_df, 'imgt', 'tenx', ['myV', 'myD', 'myJ', 'myCDR3'], 'human', custom_to_tenx_df),
    # MOUSE
    (tenx_df, 'tenx', 'adaptive', None, 'mouse', tenx_to_adapt_df),
    (tenx_df, 'tenx', 'adaptivev2', None, 'mouse', tenx_to_adapt_df),
    (adapt_df, 'adaptive', 'tenx', None, 'mouse', adapt_to_tenx_df),
    (adapt_v2_df, 'adaptivev2', 'tenx', None, 'mouse', adaptv2_to_tenx_df),
    (tenx_df, 'tenx', 'imgt', None, 'mouse', tenx_to_imgt_df),
    (imgt_df, 'imgt', 'tenx', None, 'mouse', tenx_df),
    (imgt_df, 'imgt', 'adaptive', None, 'mouse', imgt_to_adapt_df),
    (imgt_df, 'imgt', 'adaptivev2', None, 'mouse', imgt_to_adapt_df),
    (adapt_df, 'adaptive', 'imgt', None, 'mouse', adapt_to_imgt_df),
    (adapt_v2_df, 'adaptivev2', 'imgt', None, 'mouse', adaptv2_to_imgt_df),
    (custom_df, 'imgt', 'tenx', ['myV', 'myD', 'myJ', 'myCDR3'], 'mouse', custom_to_tenx_df),
    # RHESUS MACAQUE
    (tenx_df, 'tenx', 'adaptive', None, 'macaque', tenx_to_adapt_df),
    (tenx_df, 'tenx', 'adaptivev2', None, 'macaque', tenx_to_adapt_df),
    (adapt_df, 'adaptive', 'tenx', None, 'macaque', adapt_to_tenx_df),
    (adapt_v2_df, 'adaptivev2', 'tenx', None, 'macaque', adaptv2_to_tenx_df),
    (tenx_df, 'tenx', 'imgt', None, 'macaque', tenx_to_imgt_df),
    (imgt_df, 'imgt', 'tenx', None, 'macaque', tenx_df),
    (imgt_df, 'imgt', 'adaptive', None, 'macaque', imgt_to_adapt_df),
    (imgt_df, 'imgt', 'adaptivev2', None, 'macaque', imgt_to_adapt_df),
    (adapt_df, 'adaptive', 'imgt', None, 'macaque', adapt_to_imgt_df),
    (adapt_v2_df, 'adaptivev2', 'imgt', None, 'macaque', adaptv2_to_imgt_df),
    (custom_df, 'imgt', 'tenx', ['myV', 'myD', 'myJ', 'myCDR3'], 'macaque', custom_to_tenx_df)])
def test_convert_vdj(df, frm, to, frm_cols, species, out):
    result = tcrconvert.convert_vdj(df, frm, to, frm_cols, species)
    # Standardize the NA values so we can check for equality
    result.fillna('blank', inplace=True)
    out.fillna('blank', inplace=True)
    pd.testing.assert_frame_equal(result, out)


# TODO: Test bad gene names (mangled)