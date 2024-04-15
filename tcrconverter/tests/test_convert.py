import pytest
import pandas as pd
from tcrconverter import convert

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


@pytest.mark.parametrize('df,fmt_from,fmt_to,cols_use,extract_tcr,convert_cols,out', [
    (tenx_df,'tenx','adaptive',None,False,True,adapt_df),
    (tenx_df,'tenx','adaptive_v2',None,False,True,adapt_v2_df),
    (adapt_df,'adaptive','tenx',None,False,True,tenx_df),
    (adapt_v2_df,'adaptive_v2','tenx',None,False,True,tenx_df),
    (custom_df,'imgt','tenx',['myV', 'myJ', 'myCDR3'],False,True,tenx_df),
    (imgt_df,'imgt','adaptive',None,False,True,adapt_df)])
def test_convert_tcr(df,fmt_from,fmt_to,cols_use,extract_tcr,convert_cols,out):
    result = convert.convert_tcr(df,fmt_from,fmt_to,cols_use,extract_tcr,convert_cols)
    pd.testing.assert_frame_equal(result, out)

