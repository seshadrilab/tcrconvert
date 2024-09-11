Data Handling
=============

Here are some extra notes on how TCRconvert works.

Missing V, D, J, and/or C columns
---------------------------------

TCRconvert converts data for the gene columns that are in your data, not all gene columns are required.


Missing genes
-------------

Genes that are ``NA`` in the input dataframe will be ``NA`` in the output dataframe. 
Genes that are not in the lookup table (which is based on the IMGT reference) will be converted to ``NA``.


Converting from 10X
-------------------

Becuase 10X does not provide allele-level information, all genes are assigned allele ``*01``.


Converting from IMGT
--------------------

Becuase IMGT does not have standard column names, TCRconvert assumes 10X column names are used: 'v_gene', 'd_gene', 'j_gene', 'c_gene'. 
If not, please specify using ``frm_cols=``.


Converting to Adaptive
----------------------

Adaptive does not capture C genes, so all C genes become NA in the output data.


AIRR data
---------

The AIRR community specifies the columns to use in file of TCR data, but not the gene naming conventions. Those are still based on the TCR platform. 
To use AIRR-formatted files, specify the AIRR column names using ``frm_cols=``.


Custom column names
-------------------

If you're not using standard 10X, Adaptive, or Adaptive V2 column names, please specify your column names using e.g. ``frm_cols=['myV', 'myD', 'myJ']``.


Rhesus macaque TCRs
-------------------

Genes names for rhesus macaque were loaded from IMGT and `Gerdemann et al. (Front Immunol, 2022)<https://doi.org/10.3389/fimmu.2021.804932>`_.

The rehsus macaque lookup table does not currently contain delta chain constant genes, so these are converted to NA.

