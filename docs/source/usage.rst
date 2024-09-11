Usage
=====

The input and output for TCRconvert is a ``pandas`` dataframe. 
TCRconvert assumes that standard 10X, Adaptive, or Adaptive V2 column names are being used in the input data (and that 10X column names are used for IMGT data).


Convert TCR annotations between 10X, Adaptive, and IMGT formats
---------------------------------------------------------------

Load your 10X, Adaptive, or IMGT-formatted TCR data into a ``pandas`` dataframe. It may have come from files such as:

* 10X: ``filtered_contig_annotations.csv``
* Adaptive: ``SAMPLE_TCRB.tsv``
* IMGT: A custom CSV file

Then, convert to your desired format:

.. code-block:: python

    import pandas as pd
    from tcrconvert import convert

    df = pd.read_csv("filtered_contig_annotations.csv")
    converted = convert.convert_tcr(df, frm='tenx', to='adaptive')


Using the Adaptive Immune Receptor Repertoire (AIRR) format
-----------------------------------------------------------

To use Adaptive Immune Receptor Repertoire (AIRR) formatted files, simply specify the column names:

.. code-block:: python
    import pandas as pd
    import tcrconvert

    df = pd.read_csv("airr_rearrangement.tsv", sep='\t')
    converted = tcrconvert.convert_gene(df, frm='tenx', to='adaptive',
                                        frm_cols=['v_call', 'd_call', 'j_call', 'c_call'])


Gamma-delta TCRs
----------------

TCRconvert supports alpha, beta, gamma, and delta TCR genes by defualt. All are contained in the lookup tables under ``data/``


Mouse and rhesus macaque TCRs
-----------------------------

TCRconvert supports mouse and rhesus macaque TCRs, simply specify the species (since ``human`` is the default):

.. code-block:: python

    import pandas as pd
    from tcrconvert import convert

    df = pd.read_csv("filtered_contig_annotations.csv")
    converted = convert.convert_tcr(df, frm='tenx', to='adaptive', species='macaque')  # or species='mouse'
