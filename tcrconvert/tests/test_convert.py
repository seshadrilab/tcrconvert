import pytest
import pandas as pd
import tcrconvert

imgt_df = pd.DataFrame({'v_gene': ['TRAV1-2*01', 'TRBV6-1*01'],
                        'j_gene': ['TRAJ12*01', 'TRBJ2-1*01'],
                        'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

tenx_df = pd.DataFrame({'v_gene': ['TRAV1-2', 'TRBV6-1'],
                        'j_gene': ['TRAJ12', 'TRBJ2-1'],
                        'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_df = pd.DataFrame({'v_resolved': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                         'j_resolved': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                         'cdr3_amino_acid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

adapt_v2_df = pd.DataFrame({'vMaxResolved': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                            'jMaxResolved': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                            'aminoAcid': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

custom_df = pd.DataFrame({'myV': ['TRAV1-2*01', 'TRBV6-1*01'],
                          'myJ': ['TRAJ12*01', 'TRBJ2-1*01'],
                          'myCDR3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})

extract_df = pd.DataFrame({'v_gene': ['TCRAV01-02*01', 'TCRBV06-01*01'],
                           'j_gene': ['TCRAJ12*01', 'TCRBJ02-01*01'],
                           'cdr3': ['CAVMDSSYKLIF', 'CASSGLAGGYNEQFF']})


@pytest.mark.parametrize('df,fmt_from,fmt_to,cols_use,extract_tcr,convert_cols,out', [
    # 10X <-> Adaptive
    (tenx_df,'tenx','adaptive',None,False,True,adapt_df),
    (tenx_df,'tenx','adaptive_v2',None,False,True,adapt_v2_df),
    (adapt_df,'adaptive','tenx',None,False,True,tenx_df),
    (adapt_v2_df,'adaptive_v2','tenx',None,False,True,tenx_df),
    # Adaptive <-> IMGT
    # (adapt_df,'adaptive','imgt',None,False,True,imgt_df),  # KeyError: 'imgt'
    # (adapt_v2_df,'adaptive_v2','imgt',None,False,True,imgt_df),  # KeyError: 'imgt'
    (imgt_df,'imgt','adaptive',None,False,True,adapt_df),
    (imgt_df,'imgt','adaptive_v2',None,False,True,adapt_v2_df),
    # IMGT <-> 10X
    # (tenx_df,'tenx','imgt',None,False,True,imgt_df),  # KeyError: 'imgt'
    (imgt_df,'imgt','tenx',None,False,True,tenx_df),
    # Custom column names
    (custom_df,'imgt','tenx',['myV', 'myJ', 'myCDR3'],False,True,tenx_df),
    # Adaptive V2 format but user inputs "adaptive"
    (adapt_v2_df,'adaptive','tenx',None,False,True,tenx_df),
    # Extract TCRs only
    (tenx_df,'tenx','adaptive',None,True,False,extract_df)])
def test_convert_tcr(df,fmt_from,fmt_to,cols_use,extract_tcr,convert_cols,out):
    result = tcrconvert.convert_tcr(df,fmt_from,fmt_to,cols_use,extract_tcr,convert_cols)
    pd.testing.assert_frame_equal(result, out)
